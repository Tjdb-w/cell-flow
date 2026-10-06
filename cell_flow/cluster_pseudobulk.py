"""最终簇内按生物学重复（样本）汇总的 pseudobulk 差异表达。

仅在提供 ``--replicate-metadata`` 并启用 ``--cluster-pseudobulk-de``
无值开关时执行。对每个最终簇和样本，只汇总该簇内通过质控的保留细胞在
保留基因上的原始计数，把簇内样本文库（该样本在该簇内、在保留基因上的
计数和）归一到 10000 后取 log1p，再以样本为观测值对 ``group`` 做
one-versus-rest 与按 group 升序的两两比较（双侧 Welch t、均值差 log_fc、
每个簇的每项比较内跨保留基因 BH 校正）。

某簇仅在至少有两个 group 且每个 group 至少有两个含保留细胞的有效样本时
进入分析；不满足的簇跳过（无细胞样本不纳入也不补零）。全部簇都不可检验时
抛 :class:`cell_flow.errors.CellFlowDataError`（退出码 4），且不触碰结果
目录。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .errors import CellFlowDataError
from .io import ExpressionMatrix
from .markers import GroupComparison, run_group_comparisons
from .qc import QCResult
from .replicate import NORMALIZATION_TARGET, PseudobulkData


@dataclass(frozen=True)
class ClusterPseudobulkResult:
    """单个最终簇的簇内 pseudobulk 与逐比较差异表达。"""

    cluster: int
    sample_ids: List[str]          # 簇内有保留细胞的样本（首次出现顺序）
    sample_groups: List[str]       # 与 sample_ids 对齐的 group
    values: List[List[float]]      # values[gene][sample] = ln(count/total*1e4 + 1)
    comparisons: List[GroupComparison]


@dataclass(frozen=True)
class ClusterPseudobulkDE:
    """全部最终簇的簇内 pseudobulk 差异表达（仅可检验簇）与跳过簇信息。

    ``results`` 按簇标签升序；``clusters_tested`` 为可检验簇（同序），
    ``clusters_skipped`` 为不满足可检验条件的簇（同序）。``comparisons_flat``
    为 (cluster, GroupComparison) 扁平列表，按 cluster 升序、簇内先全部
    one-vs-rest（group 升序）再全部 pairwise（group 升序）。
    """

    gene_ids: List[str]
    results: List[ClusterPseudobulkResult]
    clusters_tested: List[int]
    clusters_skipped: List[int]

    @property
    def comparisons_flat(self) -> List[Tuple[int, GroupComparison]]:
        """(cluster, 比较四元组) 扁平列表，簇升序、簇内 ovr 后 pairwise。"""
        flat: List[Tuple[int, GroupComparison]] = []
        for result in self.results:
            for comp in result.comparisons:
                flat.append((result.cluster, comp))
        return flat

    @property
    def comparison_count(self) -> int:
        return sum(len(r.comparisons) for r in self.results)

    @property
    def test_count(self) -> int:
        return self.comparison_count * len(self.gene_ids)

    @property
    def min_p_value_adj(self) -> float:
        return min(
            record.p_value_adj
            for _, (_, _, _, records) in self.comparisons_flat
            for record in records
        )


def compute_cluster_pseudobulk_de(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: List[int],
    cell_samples: List[str],
    bulk: PseudobulkData,
) -> ClusterPseudobulkDE:
    """按最终簇分别汇总簇内样本 pseudobulk 并做样本级分组差异表达。

    ``labels`` 为最终簇标签、``cell_samples`` 为每个保留细胞的样本 ID，
    两者均与质控后保留细胞（``qc.kept_cells``）列序对齐；``bulk`` 给出有
    保留细胞的样本（首次出现顺序）及其分组。原始计数只在 ``qc.kept_genes``
    上、且只统计该簇内的保留细胞求和；样本是否有效按簇内是否有保留细胞
    决定，无细胞样本不纳入也不补零。
    """
    gene_ids = [matrix.gene_ids[g] for g in qc.kept_genes]
    group_of = dict(zip(bulk.sample_ids, bulk.sample_groups))
    # bulk.sample_ids 已按输入细胞首次出现排序；簇内样本次序沿用它
    sample_order = list(bulk.sample_ids)

    # 每个最终簇 -> 该簇内保留细胞的原矩阵列下标（保留细胞原顺序）
    cells_of_cluster: Dict[int, List[int]] = {}
    for c, label in enumerate(labels):
        cells_of_cluster.setdefault(label, []).append(qc.kept_cells[c])

    results: List[ClusterPseudobulkResult] = []
    clusters_tested: List[int] = []
    clusters_skipped: List[int] = []

    for cluster in sorted(cells_of_cluster):
        # 该簇内每个有效样本的保留细胞原矩阵列下标（按首次出现收集）
        cells_of_sample: Dict[str, List[int]] = {}
        seen = set()
        for c in cells_of_cluster[cluster]:
            sample_id = cell_samples[c]
            if sample_id not in seen:
                seen.add(sample_id)
                cells_of_sample[sample_id] = []
            cells_of_sample[sample_id].append(c)
        present = [s for s in sample_order if s in seen]

        # 可检验性：至少两个 group，且每个 group 至少两个簇内有效样本
        group_to_samples: Dict[str, List[str]] = {}
        for sample_id in present:
            group_to_samples.setdefault(group_of[sample_id], []).append(sample_id)
        present_groups = sorted(group_to_samples)
        if len(present_groups) < 2 or any(
            len(group_to_samples[group]) < 2 for group in present_groups
        ):
            clusters_skipped.append(cluster)
            continue
        clusters_tested.append(cluster)

        # counts[g][s]：保留基因 g 在该簇、样本 s 内的原始计数和；
        # 只加总该簇内的保留细胞
        counts: List[List[int]] = []
        for g in qc.kept_genes:
            row = matrix.counts[g]
            bulk_row: List[int] = []
            for sample_id in present:
                total = 0
                for c in cells_of_sample[sample_id]:
                    total += row[c]
                bulk_row.append(total)
            counts.append(bulk_row)

        # 簇内样本文库：该样本在该簇内、保留基因上的计数和；归一到 10000 log1p
        n_samples = len(present)
        sample_totals = [0] * n_samples
        for row in counts:
            for s, value in enumerate(row):
                sample_totals[s] += value
        values: List[List[float]] = []
        for row in counts:
            normalized_row: List[float] = []
            for s, count in enumerate(row):
                total = sample_totals[s]
                scaled = count / total * NORMALIZATION_TARGET if total > 0 else 0.0
                normalized_row.append(math.log1p(scaled))
            values.append(normalized_row)

        sample_groups = [group_of[sample_id] for sample_id in present]
        comparisons = run_group_comparisons(gene_ids, values, sample_groups)
        results.append(
            ClusterPseudobulkResult(
                cluster=cluster,
                sample_ids=present,
                sample_groups=sample_groups,
                values=values,
                comparisons=comparisons,
            )
        )

    if not results:
        raise CellFlowDataError(
            "没有任何最终簇满足簇内 pseudobulk 差异表达的可检验条件"
            "（至少两个 group 且每个 group 至少两个含保留细胞的有效样本）"
        )

    return ClusterPseudobulkDE(
        gene_ids=gene_ids,
        results=results,
        clusters_tested=clusters_tested,
        clusters_skipped=clusters_skipped,
    )
