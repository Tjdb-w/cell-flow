"""最终簇内按生物学重复（样本）汇总的 pseudobulk 差异表达。

仅在提供 ``--replicate-metadata`` 并启用 ``--cluster-pseudobulk-de``
无值开关时执行。对每个最终簇与每个样本，只汇总该簇内通过质控的保留细胞在
**质控保留基因**上的原始计数；样本文库（该簇内该样本在保留基因上的计数和）
归一到 10000 后取 log1p。差异表达以样本为观测单位，在**每个簇内部**按
``group`` 做 one-vs-rest（group 升序）与按 group 升序的两两比较，采用双侧
Welch t 检验、均值差 ``log_fc_a_vs_b``，并在每个簇的每项比较内跨保留基因做
Benjamini-Hochberg 校正（口径与 :mod:`cell_flow.markers` 的分组差异一致）。

某个簇只有在至少有两个 group、且每个 group 至少有两个“在该簇内含保留细胞”
的有效样本时才进入分析；不满足的簇记为跳过。无细胞的样本不纳入该簇，也不
补零。所有最终簇都不可检验时由调用方抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .io import ExpressionMatrix
from .markers import run_group_comparisons
from .qc import QCResult
from .replicate import ReplicateMetadata, NORMALIZATION_TARGET


@dataclass(frozen=True)
class ClusterPseudobulkRecord:
    """一个（簇, 比较, 基因）的簇内 pseudobulk 差异检验记录。"""

    cluster: int
    comparison_type: str
    group_a: str
    group_b: str                 # one-vs-rest 时为空
    gene_id: str
    mean_in_a: float
    mean_in_b: float
    log_fc_a_vs_b: float
    t_stat: float
    p_value: float
    p_value_adj: float


# (cluster, comparison_type, group_a, group_b, records)；
# one-vs-rest 的 group_b 为 ""
ClusterPseudobulkComparison = Tuple[
    int, str, str, str, List[ClusterPseudobulkRecord]
]


@dataclass(frozen=True)
class ClusterPseudobulkDE:
    """全部最终簇的簇内 pseudobulk 差异结果。

    ``comparisons`` 按 cluster 升序，每簇内先全部 one-vs-rest（group 升序，
    ``group_b`` 为空）再全部 pairwise（group 升序）；每个比较内记录按校正
    P 值升序、log_fc_a_vs_b 降序、gene_id 升序排列。``tested_clusters`` 与
    ``skipped_clusters`` 互补且各自升序，合起来即全部最终簇。
    """

    gene_ids: List[str]                       # 质控保留基因（保留行序）
    tested_clusters: List[int]
    skipped_clusters: List[int]
    comparisons: List[ClusterPseudobulkComparison]

    @property
    def comparison_count(self) -> int:
        return len(self.comparisons)

    @property
    def test_count(self) -> int:
        """全部（簇, 比较, 基因）检验数。"""
        return sum(len(records) for _, _, _, _, records in self.comparisons)

    @property
    def min_p_value_adj(self) -> float:
        return min(
            record.p_value_adj
            for _, _, _, _, records in self.comparisons
            for record in records
        )


def _cluster_pseudobulk_values(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: Sequence[int],
    replicate: ReplicateMetadata,
    cluster: int,
    sample_order: List[str],
) -> Tuple[List[str], List[str], List[List[float]]]:
    """汇总单个簇内：返回（有效样本、对齐 group、基因 × 样本 log1p 值）。

    只计入属于该簇的保留细胞；无该簇细胞的样本不在返回列中，也不补零。
    ``sample_order`` 为全样本（有保留细胞）首次出现次序，本簇取其子序列。
    """
    # 该簇每个样本中的保留细胞（原矩阵列下标）
    cells_of_sample: Dict[str, List[int]] = {}
    for c, kept_index in enumerate(qc.kept_cells):
        if labels[c] != cluster:
            continue
        sample_id = replicate.samples[matrix.cell_ids[kept_index]]
        cells_of_sample.setdefault(sample_id, []).append(kept_index)

    active_samples = [s for s in sample_order if s in cells_of_sample]
    active_groups = [
        replicate.sample_group[sample_id] for sample_id in active_samples
    ]

    # counts[g][s]：保留基因 g 在该簇有效样本 s 中的原始计数和
    counts: List[List[int]] = []
    for g in qc.kept_genes:
        row = matrix.counts[g]
        bulk_row: List[int] = []
        for sample_id in active_samples:
            total = 0
            for cell_index in cells_of_sample[sample_id]:
                total += row[cell_index]
            bulk_row.append(total)
        counts.append(bulk_row)

    # 该簇内样本文库总计数（仅在保留基因上）
    sample_totals = [0] * len(active_samples)
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

    return active_samples, active_groups, values


def _eligible(active_groups: Sequence[str]) -> bool:
    """簇准入：至少两个 group，且每个出现的 group 至少两个含保留细胞的样本。"""
    group_counts: Dict[str, int] = {}
    for group in active_groups:
        group_counts[group] = group_counts.get(group, 0) + 1
    return len(group_counts) >= 2 and all(n >= 2 for n in group_counts.values())


def compute_cluster_pseudobulk_de(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: Sequence[int],
    replicate: ReplicateMetadata,
) -> ClusterPseudobulkDE:
    """按最终簇逐簇汇总样本 pseudobulk 并做簇内分组差异表达。

    ``labels`` 与 ``qc.kept_cells`` 对齐（每个保留细胞一个最终簇标签）。
    样本的全局次序沿用 pseudobulk 基线：按其输入细胞在矩阵列序中的首次出现；
    每个簇只取其中在该簇内有保留细胞的样本子序列。
    """
    # 全样本（质控后有保留细胞）首次出现顺序；与 build_pseudobulk 一致
    sample_order: List[str] = []
    seen = set()
    for kept_index in qc.kept_cells:
        sample_id = replicate.samples[matrix.cell_ids[kept_index]]
        if sample_id not in seen:
            seen.add(sample_id)
            sample_order.append(sample_id)

    gene_ids = [matrix.gene_ids[g] for g in qc.kept_genes]

    tested_clusters: List[int] = []
    skipped_clusters: List[int] = []
    comparisons: List[ClusterPseudobulkComparison] = []

    for cluster in sorted(set(labels)):
        _, active_groups, values = _cluster_pseudobulk_values(
            matrix, qc, labels, replicate, cluster, sample_order
        )
        if not _eligible(active_groups):
            skipped_clusters.append(cluster)
            continue
        tested_clusters.append(cluster)

        # 以样本为观测复用分组比较口径（one-vs-rest + 升序两两、Welch t、
        # 比较内 BH）；记录字段与 group marker 完全一致
        group_comparisons = run_group_comparisons(
            gene_ids, values, active_groups
        )
        for comparison_type, group_a, group_b, records in group_comparisons:
            cluster_records = [
                ClusterPseudobulkRecord(
                    cluster=cluster,
                    comparison_type=record.comparison_type,
                    group_a=record.group_a,
                    group_b=record.group_b,
                    gene_id=record.gene_id,
                    mean_in_a=record.mean_in_a,
                    mean_in_b=record.mean_in_b,
                    log_fc_a_vs_b=record.log_fc_a_vs_b,
                    t_stat=record.t_stat,
                    p_value=record.p_value,
                    p_value_adj=record.p_value_adj,
                )
                for record in records
            ]
            comparisons.append(
                (cluster, comparison_type, group_a, group_b, cluster_records)
            )

    return ClusterPseudobulkDE(
        gene_ids=gene_ids,
        tested_clusters=tested_clusters,
        skipped_clusters=skipped_clusters,
        comparisons=comparisons,
    )
