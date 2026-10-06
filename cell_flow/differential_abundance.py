"""簇级样本差异丰度：质控后细胞在最终簇上的样本构成与分组差异检验。

仅在启用 ``--differential-abundance``（且必带 ``--replicate-metadata``）时
执行。只用质控后保留细胞（启用双细胞识别时为其最终细胞）、最终簇标签与
重复元数据给出的样本分组：

- 每个样本在每个最终簇上的细胞数 ``n_cells``（缺簇计 0）、该样本保留细胞
  总数 ``total_cells`` 与比例 ``proportion = n_cells / total_cells``，
  覆盖样本 × 簇的全部组合；
- 以样本为观测单位，对每个 group 先做 one-vs-rest、再按 group 升序两两
  比较，对每个簇执行双侧 Welch t 检验，差异值为 a 组均值减 b 组均值，
  并在单个比较内跨簇做 Benjamini-Hochberg 校正。

质控后样本无保留细胞、非空 group 不足两个或任一 group 有效重复不足两个，
一律抛 :class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple

from .errors import CellFlowDataError
from .io import ExpressionMatrix
from .linalg import benjamini_hochberg, welch_ttest
from .qc import QCResult
from .replicate import ReplicateMetadata

ONE_VS_REST = "one_vs_rest"
PAIRWISE = "pairwise"

TOP_N_ABUNDANCE = 20


@dataclass(frozen=True)
class ClusterAbundanceRow:
    """一个（样本, 簇）组合上的细胞丰度；缺簇同样给出一行（n_cells=0）。"""

    sample_id: str
    group: str
    cluster: int
    n_cells: int
    total_cells: int
    proportion: float


@dataclass(frozen=True)
class ClusterDifferentialRow:
    """一个比较中一个簇的样本比例差异检验结果。"""

    comparison_type: str
    group_a: str
    group_b: str  # one-vs-rest 时为空
    cluster: int
    mean_in_a: float
    mean_in_b: float
    proportion_difference_a_vs_b: float
    t_stat: float
    p_value: float
    p_value_adj: float


# (comparison_type, group_a, group_b, records)；one-vs-rest 的 group_b 为 ""
DifferentialComparison = Tuple[str, str, str, List[ClusterDifferentialRow]]


@dataclass(frozen=True)
class ClusterAbundanceResult:
    """样本 × 簇丰度表与逐比较的簇差异检验结果。

    ``rows`` 按 sample_id、cluster 升序覆盖全部组合；``comparisons`` 先全部
    one-vs-rest（group 升序）再全部 pairwise（group 升序），比较内记录按
    p_value_adj 升序、差异值降序、cluster 升序排列。``sample_ids`` 按
    sample_id 升序，``cluster_ids`` 按簇编号升序。
    """

    rows: List[ClusterAbundanceRow]
    comparisons: List[DifferentialComparison]
    sample_ids: List[str]
    cluster_ids: List[int]


def _data_fail(message: str) -> None:
    raise CellFlowDataError(message)


def build_cluster_abundance(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: List[int],
    replicate: ReplicateMetadata,
) -> ClusterAbundanceResult:
    """汇总质控后细胞在最终簇上的样本丰度并做分组 Welch t/BH 差异检验。

    ``labels`` 与 ``qc.kept_cells``（保留细胞原列序）对齐；样本分组取自
    重复元数据。观测顺序统一按 sample_id 升序，保证浮点求和与输出确定。
    """
    # 保留细胞 -> 样本；labels 与保留细胞列序一一对应
    cells_of_sample: Dict[str, List[int]] = {}
    for pos, c in enumerate(qc.kept_cells):
        sample_id = replicate.samples[matrix.cell_ids[c]]
        cells_of_sample.setdefault(sample_id, []).append(pos)

    empty_samples = [
        sample_id
        for sample_id in replicate.sample_order
        if sample_id not in cells_of_sample
    ]
    if empty_samples:
        _data_fail(
            f"质控后样本 {empty_samples[0]!r} 等 {len(empty_samples)} 个"
            f"无保留细胞，无法进行簇级差异丰度分析"
        )

    sample_ids = sorted(cells_of_sample)
    group_of_sample = {
        sample_id: replicate.sample_group[sample_id] for sample_id in sample_ids
    }
    group_to_samples: Dict[str, List[str]] = {}
    for sample_id in sample_ids:
        group_to_samples.setdefault(group_of_sample[sample_id], []).append(sample_id)
    present_groups = sorted(group_to_samples)
    if len(present_groups) < 2:
        _data_fail(
            f"质控后非空分组仅 {len(present_groups)} 个，不足两个，"
            f"无法进行簇级差异丰度分析"
        )
    single_replicate_groups = [
        group for group in present_groups if len(group_to_samples[group]) < 2
    ]
    if single_replicate_groups:
        group = single_replicate_groups[0]
        _data_fail(
            f"分组 {group!r} 质控后有效重复仅 "
            f"{len(group_to_samples[group])} 个，少于两个，"
            f"无法进行簇级差异丰度分析"
        )

    cluster_ids = sorted(set(labels))

    # 每样本保留细胞总数；n[sample][cluster] 为该样本在该簇的细胞数
    totals = {sample_id: len(cells_of_sample[sample_id]) for sample_id in sample_ids}
    counts: Dict[Tuple[str, int], int] = {}
    for sample_id in sample_ids:
        for pos in cells_of_sample[sample_id]:
            key = (sample_id, labels[pos])
            counts[key] = counts.get(key, 0) + 1

    # 观测矩阵 proportions[sample_index][cluster_index]，行列均按升序
    proportions: List[List[float]] = []
    rows: List[ClusterAbundanceRow] = []
    for sample_id in sample_ids:
        group = group_of_sample[sample_id]
        total = totals[sample_id]
        sample_proportions: List[float] = []
        for cluster in cluster_ids:
            n_cells = counts.get((sample_id, cluster), 0)
            proportion = n_cells / total
            sample_proportions.append(proportion)
            rows.append(
                ClusterAbundanceRow(
                    sample_id=sample_id,
                    group=group,
                    cluster=cluster,
                    n_cells=n_cells,
                    total_cells=total,
                    proportion=proportion,
                )
            )
        proportions.append(sample_proportions)

    comparisons = _run_comparisons(
        cluster_ids, proportions, [group_of_sample[s] for s in sample_ids]
    )

    return ClusterAbundanceResult(
        rows=rows,
        comparisons=comparisons,
        sample_ids=sample_ids,
        cluster_ids=cluster_ids,
    )


def _run_comparisons(
    cluster_ids: List[int],
    proportions: List[List[float]],
    obs_groups: List[str],
) -> List[DifferentialComparison]:
    """以样本为观测的簇比例分组比较：先 one-vs-rest 再两两，组内升序。

    双侧 Welch t，差异值为 a 组簇比例均值减 b 组均值；BH 校正以单个比较内
    的全部簇为一个校正家族。记录按 p_value_adj 升序、差异值降序、
    cluster 升序排列。
    """
    group_to_obs: Dict[str, List[int]] = {}
    for i, group in enumerate(obs_groups):
        group_to_obs.setdefault(group, []).append(i)

    present = sorted(group_to_obs)
    comparisons: List[DifferentialComparison] = []

    def compare(
        obs_a: List[int],
        obs_b: List[int],
        comparison_type: str,
        group_a: str,
        group_b: str,
    ) -> List[ClusterDifferentialRow]:
        records: List[ClusterDifferentialRow] = []
        pvalues: List[float] = []
        for k, cluster in enumerate(cluster_ids):
            values_a = [proportions[i][k] for i in obs_a]
            values_b = [proportions[i][k] for i in obs_b]
            mean_a = sum(values_a) / len(values_a)
            mean_b = sum(values_b) / len(values_b)
            difference = mean_a - mean_b
            # 两组方差均为 0 且均值相同时 welch_ttest 给出 t=0、p=1
            t_stat, p_value = welch_ttest(values_a, values_b)
            pvalues.append(p_value)
            records.append(
                ClusterDifferentialRow(
                    comparison_type=comparison_type,
                    group_a=group_a,
                    group_b=group_b,
                    cluster=cluster,
                    mean_in_a=mean_a,
                    mean_in_b=mean_b,
                    proportion_difference_a_vs_b=difference,
                    t_stat=t_stat,
                    p_value=p_value,
                    p_value_adj=0.0,
                )
            )
        adjusted = benjamini_hochberg(pvalues)
        records = [
            ClusterDifferentialRow(
                comparison_type=r.comparison_type,
                group_a=r.group_a,
                group_b=r.group_b,
                cluster=r.cluster,
                mean_in_a=r.mean_in_a,
                mean_in_b=r.mean_in_b,
                proportion_difference_a_vs_b=r.proportion_difference_a_vs_b,
                t_stat=r.t_stat,
                p_value=r.p_value,
                p_value_adj=adjusted[k],
            )
            for k, r in enumerate(records)
        ]
        records.sort(
            key=lambda r: (
                r.p_value_adj,
                -r.proportion_difference_a_vs_b,
                r.cluster,
            )
        )
        return records

    for group in present:
        obs_in = list(group_to_obs[group])
        obs_out = [
            i for i in range(len(obs_groups)) if obs_groups[i] != group
        ]
        records = compare(obs_in, obs_out, ONE_VS_REST, group, "")
        comparisons.append((ONE_VS_REST, group, "", records))
    for i, group_a in enumerate(present):
        for group_b in present[i + 1:]:
            records = compare(
                group_to_obs[group_a],
                group_to_obs[group_b],
                PAIRWISE,
                group_a,
                group_b,
            )
            comparisons.append((PAIRWISE, group_a, group_b, records))
    return comparisons
