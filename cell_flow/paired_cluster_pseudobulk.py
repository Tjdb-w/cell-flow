"""最终簇内按配对生物学重复（样本）的 pseudobulk 差异表达。

仅在同时提供 ``--replicate-metadata`` 与 ``--paired-replicate-metadata``、
并启用 ``--paired-cluster-pseudobulk-de`` 无值开关时执行。计数保留范围、
汇总方式与文库归一化完全沿用簇内 pseudobulk 基线
（:mod:`cell_flow.cluster_pseudobulk`）：对每个最终簇与样本，只汇总该簇内
通过质控的保留细胞在质控保留基因上的原始计数；样本文库（该簇内该样本在
保留基因上的计数和）归一到 10000 后取 log1p。

质控后恰有两个 group 时按字典序确定 a、b；每个 pair 必须各含一个 a、b
样本且样本不得重复（由配对元数据逐样本唯一行保证），否则按 pair 归属错误
抛 :class:`cell_flow.errors.CellFlowInputError`（退出码 2）。每个簇只用两个
成员在该簇内都有保留细胞的完整 pair；完整 pair 少于两个的簇跳过，无任何
可检验簇时抛 :class:`cell_flow.errors.CellFlowDataError`（退出码 4）。

每个基因对 a 减 b 的 pair 差值做双侧配对 t 检验：``t_stat`` 为差值均值除以
差值样本标准差再乘 pair 数平方根，自由度为 pair 数减一，簇内跨基因做
Benjamini-Hochberg 校正。零方差且均值为零时 ``t_stat=0``、``p_value=1``；
零方差且均值非零时 ``t_stat`` 为 ``inf``/``-inf``、``p_value=0``。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .cluster_pseudobulk import _cluster_pseudobulk_values
from .errors import CellFlowDataError, CellFlowInputError
from .io import ExpressionMatrix
from .linalg import benjamini_hochberg, t_sf_two_sided, unbiased_variance
from .paired_replicate import PairedReplicateMetadata
from .qc import QCResult
from .replicate import ReplicateMetadata


@dataclass(frozen=True)
class PairedClusterPseudobulkRecord:
    """一个（簇, 基因）的簇内配对 pseudobulk 差异检验记录。"""

    cluster: int
    group_a: str
    group_b: str
    pair_count: int
    gene_id: str
    mean_difference: float
    t_stat: float
    p_value: float
    p_value_adj: float


# (cluster, pair_count, records)
PairedClusterResult = Tuple[int, int, List[PairedClusterPseudobulkRecord]]


@dataclass(frozen=True)
class PairedClusterPseudobulkDE:
    """全部最终簇的簇内配对 pseudobulk 差异结果。

    ``results`` 按 cluster 升序；每个簇内记录按校正 P 值升序、
    mean_difference 降序、gene_id 升序排列。``tested_clusters`` 与
    ``skipped_clusters`` 互补且各自升序，合起来即全部最终簇。
    """

    gene_ids: List[str]                       # 质控保留基因（保留行序）
    group_a: str
    group_b: str
    tested_clusters: List[int]
    skipped_clusters: List[int]
    results: List[PairedClusterResult]

    @property
    def cluster_count(self) -> int:
        return len(self.results)

    @property
    def test_count(self) -> int:
        """全部（簇, 基因）检验数。"""
        return sum(len(records) for _, _, records in self.results)

    @property
    def min_p_value_adj(self) -> float:
        return min(
            record.p_value_adj
            for _, _, records in self.results
            for record in records
        )


def resolve_pair_groups(
    sample_ids: Sequence[str],
    sample_group: Dict[str, str],
    paired: PairedReplicateMetadata,
) -> Tuple[str, str, Dict[str, List[str]]]:
    """校验质控后恰两个 group，并把样本按 pair 归集。

    返回 ``(group_a, group_b, members_of_pair)``：a、b 按 group 字典序确定；
    ``members_of_pair`` 按 pair 标识字典序排列。任何 pair 不是恰含一个 a、
    一个 b 样本（单成员、同组两成员、多于两个成员）都属 pair 归属错误，
    抛 :class:`CellFlowInputError`；质控后 group 数不为 2 同样无法定义
    a、b，按 pair 归属错误处理。
    """
    groups = sorted({sample_group[sample_id] for sample_id in sample_ids})
    if len(groups) != 2:
        raise CellFlowInputError(
            f"配对差异表达要求质控后恰有两个 group，实际为 {len(groups)} 个："
            f"{groups}"
        )
    group_a, group_b = groups

    members_of_pair: Dict[str, List[str]] = {}
    for pair_id in paired.pair_order:
        members_of_pair.setdefault(pair_id, [])
    for sample_id in sample_ids:
        pair_id = paired.sample_pair[sample_id]
        members_of_pair.setdefault(pair_id, []).append(sample_id)

    for pair_id in sorted(members_of_pair):
        members = members_of_pair[pair_id]
        in_a = [s for s in members if sample_group[s] == group_a]
        in_b = [s for s in members if sample_group[s] == group_b]
        if len(members) != 2 or len(in_a) != 1 or len(in_b) != 1:
            detail = "、".join(
                f"{s}({sample_group[s]})" for s in sorted(members)
            )
            raise CellFlowInputError(
                f"配对 {pair_id!r} 必须各含一个 {group_a!r} 与一个 "
                f"{group_b!r} 样本，实际成员为：{detail or '空'}"
            )

    return group_a, group_b, members_of_pair


def _paired_t_test(differences: Sequence[float]) -> Tuple[float, float]:
    """双侧配对 t 检验，返回 (t 统计量, 双侧 P 值)。

    t = mean(d) / sd(d) * sqrt(n)，自由度 n-1。零方差退化：均值为 0 时
    t=0、p=1；均值非零时 t 为与均值同号的 inf、p=0。
    """
    n = len(differences)
    mean_d = sum(differences) / n
    variance = unbiased_variance(differences, mean_d)
    if variance <= 0.0:
        if mean_d == 0.0:
            return 0.0, 1.0
        return (math.inf if mean_d > 0.0 else -math.inf), 0.0
    t_stat = mean_d / math.sqrt(variance) * math.sqrt(n)
    return t_stat, t_sf_two_sided(t_stat, n - 1)


def compute_paired_cluster_pseudobulk_de(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: Sequence[int],
    replicate: ReplicateMetadata,
    paired: PairedReplicateMetadata,
    sample_order: Sequence[str],
) -> PairedClusterPseudobulkDE:
    """按最终簇逐簇汇总样本 pseudobulk 并做簇内配对差异表达。

    ``labels`` 与 ``qc.kept_cells`` 对齐（每个保留细胞一个最终簇标签）；
    ``sample_order`` 为质控后有保留细胞的全部样本（pseudobulk 基线样本序，
    即样本输入细胞在矩阵列序中的首次出现），配对表已与之覆盖一致。
    """
    group_a, group_b, members_of_pair = resolve_pair_groups(
        sample_order, replicate.sample_group, paired
    )

    gene_ids = [matrix.gene_ids[g] for g in qc.kept_genes]

    tested_clusters: List[int] = []
    skipped_clusters: List[int] = []
    results: List[PairedClusterResult] = []

    for cluster in sorted(set(labels)):
        active_samples, _, values = _cluster_pseudobulk_values(
            matrix, qc, labels, replicate, cluster, sample_order
        )
        active = set(active_samples)

        # 只用两个成员在该簇内都有保留细胞的完整 pair；pair 按标识升序，
        # 使差值序列（均值/样本方差与其顺序无关）的构造完全确定
        complete_pairs: List[Tuple[str, str, str]] = []
        for pair_id in sorted(members_of_pair):
            members = members_of_pair[pair_id]
            sample_a = next(s for s in members if replicate.sample_group[s] == group_a)
            sample_b = next(s for s in members if replicate.sample_group[s] == group_b)
            if sample_a in active and sample_b in active:
                complete_pairs.append((pair_id, sample_a, sample_b))

        if len(complete_pairs) < 2:
            skipped_clusters.append(cluster)
            continue
        tested_clusters.append(cluster)

        index_of = {sample_id: i for i, sample_id in enumerate(active_samples)}
        n_pairs = len(complete_pairs)

        records: List[PairedClusterPseudobulkRecord] = []
        pvalues: List[float] = []
        for g, gene_id in enumerate(gene_ids):
            row = values[g]
            differences = [
                row[index_of[sample_a]] - row[index_of[sample_b]]
                for _, sample_a, sample_b in complete_pairs
            ]
            mean_difference = sum(differences) / n_pairs
            t_stat, p_value = _paired_t_test(differences)
            pvalues.append(p_value)
            records.append(
                PairedClusterPseudobulkRecord(
                    cluster=cluster,
                    group_a=group_a,
                    group_b=group_b,
                    pair_count=n_pairs,
                    gene_id=gene_id,
                    mean_difference=mean_difference,
                    t_stat=t_stat,
                    p_value=p_value,
                    p_value_adj=0.0,
                )
            )

        adjusted = benjamini_hochberg(pvalues)
        records = [
            PairedClusterPseudobulkRecord(
                cluster=r.cluster,
                group_a=r.group_a,
                group_b=r.group_b,
                pair_count=r.pair_count,
                gene_id=r.gene_id,
                mean_difference=r.mean_difference,
                t_stat=r.t_stat,
                p_value=r.p_value,
                p_value_adj=adjusted[g],
            )
            for g, r in enumerate(records)
        ]
        records.sort(
            key=lambda r: (r.p_value_adj, -r.mean_difference, r.gene_id)
        )
        results.append((cluster, n_pairs, records))

    if not tested_clusters:
        raise CellFlowDataError(
            "没有任何最终簇具备至少两个成员均含保留细胞的完整 pair，"
            "无法进行簇内配对 pseudobulk 差异表达"
        )

    return PairedClusterPseudobulkDE(
        gene_ids=gene_ids,
        group_a=group_a,
        group_b=group_b,
        tested_clusters=tested_clusters,
        skipped_clusters=skipped_clusters,
        results=results,
    )
