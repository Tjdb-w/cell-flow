"""簇级样本差异丰度分析（``--differential-abundance``）。

仅在提供 ``--replicate-metadata`` 并启用 ``--differential-abundance``
无值开关时执行。只用质控后保留细胞、最终簇标签与样本分组：对每个
（样本, 簇）组合统计保留细胞数 ``n_cells``、该样本保留细胞总数
``total_cells`` 与比例 ``proportion = n_cells / total_cells``（样本中
不属于该簇时 ``n_cells`` 计 0），随后以样本为观测单位，对每个簇的
样本比例做分组差异检验。

比较口径沿用分组差异基线（:func:`cell_flow.markers.run_group_comparisons`
的次序）：每个 group 先 one-vs-rest（group 升序，``group_b`` 为空），
再按 group 升序两两比较；检验为样本比例的双侧 Welch t，差异值为
a 组均值减 b 组均值，BH 校正以单个比较内的全部簇为一个校正家族，
每个比较内按校正 P 值升序、差异值降序、簇升序排列。与逐基因检验
唯一的差别：两组样本比例的方差均为 0 时，无论均值是否相等，
``t_stat = 0``、``p_value = 1``。

质控后样本无保留细胞、非空 group 不足两个或任一 group 有效重复不足
两个，由 pseudobulk 基线（``build_pseudobulk``）在本分析之前先行抛出
数据错误（退出码 4），此处不再重复校验。
"""

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .linalg import benjamini_hochberg, unbiased_variance, welch_ttest
from .markers import ONE_VS_REST, PAIRWISE
from .replicate import PseudobulkData


@dataclass(frozen=True)
class ClusterAbundanceRow:
    """一个（样本, 簇）组合的丰度记录。"""

    sample_id: str
    group: str
    cluster: int
    n_cells: int
    total_cells: int
    proportion: float


@dataclass(frozen=True)
class AbundanceComparisonRecord:
    """一个（比较, 簇）的差异丰度检验记录。"""

    comparison_type: str
    group_a: str
    group_b: str                 # one-vs-rest 时为空
    cluster: int
    mean_in_a: float
    mean_in_b: float
    difference: float            # mean_in_a - mean_in_b
    t_stat: float
    p_value: float
    p_value_adj: float


# (comparison_type, group_a, group_b, records)；one-vs-rest 的 group_b 为 ""
AbundanceComparison = Tuple[str, str, str, List[AbundanceComparisonRecord]]


@dataclass(frozen=True)
class DifferentialAbundance:
    """簇级样本丰度表与逐比较差异丰度结果。"""

    rows: List[ClusterAbundanceRow]        # 样本 × 簇全组合，sample_id、cluster 升序
    clusters: List[int]                    # 最终簇升序
    sample_ids: List[str]                  # pseudobulk 样本次序（首次出现顺序）
    sample_groups: List[str]               # 与 sample_ids 对齐的 group
    comparisons: List[AbundanceComparison]


def _proportion_ttest(
    values_a: Sequence[float], values_b: Sequence[float]
) -> Tuple[float, float]:
    """样本比例的 Welch t 检验；两组方差均为 0 时恒为 t=0、p=1。"""
    mu_a = sum(values_a) / len(values_a)
    mu_b = sum(values_b) / len(values_b)
    if (
        unbiased_variance(values_a, mu_a) == 0.0
        and unbiased_variance(values_b, mu_b) == 0.0
    ):
        return 0.0, 1.0
    return welch_ttest(values_a, values_b)


def _cluster_comparison(
    clusters: List[int],
    proportions: List[List[float]],
    obs_a: List[int],
    obs_b: List[int],
    comparison_type: str,
    group_a: str,
    group_b: str,
) -> List[AbundanceComparisonRecord]:
    """两组样本之间逐簇的比例检验：差异值即均值差，比较内跨簇 BH 校正。

    ``proportions[k][s]`` 为簇 ``clusters[k]`` 在第 ``s`` 个样本上的比例；
    ``obs_a``/``obs_b`` 为两组样本的列下标。"""
    records: List[AbundanceComparisonRecord] = []
    pvalues: List[float] = []
    for k, cluster in enumerate(clusters):
        row = proportions[k]
        values_a = [row[s] for s in obs_a]
        values_b = [row[s] for s in obs_b]
        mean_a = sum(values_a) / len(values_a)
        mean_b = sum(values_b) / len(values_b)
        difference = mean_a - mean_b
        t_stat, p_value = _proportion_ttest(values_a, values_b)
        pvalues.append(p_value)
        records.append(
            AbundanceComparisonRecord(
                comparison_type=comparison_type,
                group_a=group_a,
                group_b=group_b,
                cluster=cluster,
                mean_in_a=mean_a,
                mean_in_b=mean_b,
                difference=difference,
                t_stat=t_stat,
                p_value=p_value,
                p_value_adj=0.0,
            )
        )

    adjusted = benjamini_hochberg(pvalues)
    records = [
        AbundanceComparisonRecord(
            comparison_type=r.comparison_type,
            group_a=r.group_a,
            group_b=r.group_b,
            cluster=r.cluster,
            mean_in_a=r.mean_in_a,
            mean_in_b=r.mean_in_b,
            difference=r.difference,
            t_stat=r.t_stat,
            p_value=r.p_value,
            p_value_adj=adjusted[k],
        )
        for k, r in enumerate(records)
    ]
    records.sort(key=lambda r: (r.p_value_adj, -r.difference, r.cluster))
    return records


def compute_differential_abundance(
    labels: List[int],
    cell_samples: List[str],
    bulk: PseudobulkData,
) -> DifferentialAbundance:
    """按最终簇标签与样本分组计算簇级丰度与样本级分组差异。

    ``labels`` 为最终簇标签、``cell_samples`` 为每个保留细胞的样本 ID，
    两者均与质控后保留细胞列序对齐；``bulk.sample_ids``/``sample_groups``
    给出有保留细胞的样本（首次出现顺序）及其分组。样本与分组规模约束
    （无空样本、至少两个 group、每组至少两个重复）由 ``bulk`` 的构建
    阶段保证。
    """
    clusters = sorted(set(labels))
    cluster_position = {cluster: k for k, cluster in enumerate(clusters)}
    n_clusters = len(clusters)
    n_samples = len(bulk.sample_ids)
    sample_position = {
        sample_id: s for s, sample_id in enumerate(bulk.sample_ids)
    }

    # counts[k][s]：簇 k 在样本 s 中的保留细胞数；total_cells[s]：样本保留细胞总数
    counts = [[0] * n_samples for _ in range(n_clusters)]
    total_cells = [0] * n_samples
    for label, sample_id in zip(labels, cell_samples):
        counts[cluster_position[label]][sample_position[sample_id]] += 1
        total_cells[sample_position[sample_id]] += 1

    proportions = [
        [counts[k][s] / total_cells[s] for s in range(n_samples)]
        for k in range(n_clusters)
    ]

    group_of = dict(zip(bulk.sample_ids, bulk.sample_groups))
    rows: List[ClusterAbundanceRow] = []
    for sample_id in sorted(bulk.sample_ids):
        s = sample_position[sample_id]
        for k, cluster in enumerate(clusters):
            rows.append(
                ClusterAbundanceRow(
                    sample_id=sample_id,
                    group=group_of[sample_id],
                    cluster=cluster,
                    n_cells=counts[k][s],
                    total_cells=total_cells[s],
                    proportion=proportions[k][s],
                )
            )

    group_to_obs: Dict[str, List[int]] = {}
    for s, group in enumerate(bulk.sample_groups):
        group_to_obs.setdefault(group, []).append(s)

    present = sorted(group_to_obs)
    comparisons: List[AbundanceComparison] = []
    for group in present:
        obs_in = group_to_obs[group]
        obs_out = [s for s in range(n_samples) if bulk.sample_groups[s] != group]
        records = _cluster_comparison(
            clusters, proportions, obs_in, obs_out, ONE_VS_REST, group, ""
        )
        comparisons.append((ONE_VS_REST, group, "", records))
    for i, group_a in enumerate(present):
        for group_b in present[i + 1:]:
            records = _cluster_comparison(
                clusters,
                proportions,
                group_to_obs[group_a],
                group_to_obs[group_b],
                PAIRWISE,
                group_a,
                group_b,
            )
            comparisons.append((PAIRWISE, group_a, group_b, records))

    return DifferentialAbundance(
        rows=rows,
        clusters=clusters,
        sample_ids=list(bulk.sample_ids),
        sample_groups=list(bulk.sample_groups),
        comparisons=comparisons,
    )
