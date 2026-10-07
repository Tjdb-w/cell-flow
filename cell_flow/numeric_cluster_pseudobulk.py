"""最终簇内数值型样本协变量的 pseudobulk 独立检验。

仅在同时提供 ``--replicate-metadata``、``--cluster-pseudobulk-de``、
``--cluster-pseudobulk-adjusted-de``（必带 ``--pseudobulk-covariates``）
并给出 ``--pseudobulk-numeric-covariates`` 时执行。计数保留范围、汇总方式、
文库归一化与簇/比较枚举次序完全沿用簇内 pseudobulk 基线
（:mod:`cell_flow.cluster_pseudobulk`）：对每个最终簇与样本，只汇总该簇内
通过质控的保留细胞在质控保留基因上的原始计数；样本文库（该簇内该样本在
保留基因上的计数和）归一到 10000 后取 log1p。本模块产出独立于
``cluster_pseudobulk_adjusted_de.tsv`` 的检验结果，不改写任何既有文件。

每个簇与每项比较（one-vs-rest 按 group 升序、再按 group 升序两两比较，
次序与既有簇内 pseudobulk 结果一致）只用该簇内的有效样本（该簇内含保留
细胞的样本；两两比较只用两个目标 group 的样本），以样本为观测拟合线性
模型：截距 + 组别项 + 每个分类协变量以其字典序最小水平为参照的哑变量
（与协变量校正基线一致）+ 每个数值协变量一个连续列（表头顺序）。对
**组别项与每个数值协变量列**分别检验：系数即 ``effect``，同时输出标准误、
残差自由度、t 统计量与双侧 P 值，并在每个簇比较内跨全部（保留基因 ×
被检验项）记录做 Benjamini-Hochberg 校正。记录按校正 P 值升序、
``effect`` 降序、``gene_id`` 升序排列。

设计矩阵秩不足或残差自由度不大于零时抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .adjusted_cluster_pseudobulk import _invert_matrix
from .cluster_pseudobulk import _cluster_pseudobulk_values, _eligible
from .errors import CellFlowDataError
from .io import ExpressionMatrix
from .linalg import benjamini_hochberg, t_sf_two_sided
from .markers import ONE_VS_REST, PAIRWISE
from .pseudobulk_covariates import PseudobulkCovariates
from .pseudobulk_numeric_covariates import PseudobulkNumericCovariates
from .qc import QCResult
from .replicate import ReplicateMetadata

# 组别被检验项在 covariate 列中的固定标签
GROUP_TERM = "group"


@dataclass(frozen=True)
class NumericCovariateRecord:
    """一个（簇, 比较, 被检验项, 基因）的数值协变量 pseudobulk 检验记录。"""

    cluster: int
    comparison_type: str
    group_a: str
    group_b: str                 # one-vs-rest 时为空
    covariate: str               # 被检验项：组别为 "group"，否则为数值协变量列名
    gene_id: str
    effect: float
    std_error: float
    df: int
    t_stat: float
    p_value: float
    p_value_adj: float


# (cluster, comparison_type, group_a, group_b, records)；
# one-vs-rest 的 group_b 为 ""
NumericCovariateComparison = Tuple[
    int, str, str, str, List[NumericCovariateRecord]
]


@dataclass(frozen=True)
class NumericCovariateClusterPseudobulkDE:
    """全部最终簇的数值协变量簇内 pseudobulk 检验结果。

    ``comparisons`` 按 cluster 升序，每簇内先全部 one-vs-rest（group 升序，
    ``group_b`` 为空）再全部 pairwise（group 升序），与既有簇内 pseudobulk
    结果的比较次序一致；每个比较内记录按校正 P 值升序、effect 降序、
    gene_id 升序排列。``tested_clusters`` 与 ``skipped_clusters`` 互补且
    各自升序，合起来即全部最终簇。
    """

    gene_ids: List[str]                       # 质控保留基因（保留行序）
    covariate_names: List[str]                # 参与检验的数值协变量（表头顺序）
    tested_clusters: List[int]
    skipped_clusters: List[int]
    comparisons: List[NumericCovariateComparison]

    @property
    def comparison_count(self) -> int:
        return len(self.comparisons)

    @property
    def test_count(self) -> int:
        """全部（簇, 比较, 被检验项, 基因）检验数。"""
        return sum(len(records) for _, _, _, _, records in self.comparisons)

    @property
    def min_p_value_adj(self) -> float:
        return min(
            record.p_value_adj
            for _, _, _, _, records in self.comparisons
            for record in records
        )


def _design_matrix(
    obs: Sequence[int],
    obs_group_indicator: Sequence[int],
    active_samples: Sequence[str],
    covariates: PseudobulkCovariates,
    numeric: PseudobulkNumericCovariates,
) -> List[List[float]]:
    """构造设计矩阵：截距、组别项、分类协变量哑变量，再按表头顺序接连续列。

    分类协变量的水平取本比较实际使用样本中的取值，字典序最小水平为参照，
    其余每个水平一列（按字典序）；只出现一个水平的协变量不产生哑变量列
    （与协变量校正基线一致）。数值协变量每个一列，取样本的原始有限取值。
    """
    # 每个分类协变量的非参照水平（哑变量列），按表头顺序、水平字典序展开
    dummies: List[Tuple[str, str]] = []
    for name in covariates.covariate_names:
        present = sorted(
            {covariates.sample_levels[active_samples[j]][name] for j in obs}
        )
        for level in present[1:]:
            dummies.append((name, level))

    rows: List[List[float]] = []
    for position, i in enumerate(obs):
        sample_id = active_samples[i]
        levels = covariates.sample_levels[sample_id]
        values = numeric.sample_values[sample_id]
        row = [1.0, float(obs_group_indicator[position])]
        for name, level in dummies:
            row.append(1.0 if levels[name] == level else 0.0)
        for name in numeric.covariate_names:
            row.append(values[name])
        rows.append(row)
    return rows


def _ols_fit(
    design: List[List[float]],
    inverse: List[List[float]],
    y: Sequence[float],
    df: int,
) -> Tuple[List[float], float]:
    """对单个基因拟合 OLS，返回（全部系数, 残差方差 sigma2）。"""
    n_params = len(design[0])
    xty = [0.0] * n_params
    for i, row in enumerate(design):
        value = y[i]
        for k in range(n_params):
            xty[k] += row[k] * value
    beta = [
        sum(inverse[k][j] * xty[j] for j in range(n_params))
        for k in range(n_params)
    ]
    rss = 0.0
    for i, row in enumerate(design):
        fitted = sum(row[k] * beta[k] for k in range(n_params))
        residual = y[i] - fitted
        rss += residual * residual
    return beta, rss / df


def _term_test(
    effect: float, variance: float, df: int
) -> Tuple[float, float, float]:
    """由系数与其方差求（标准误, t, 双侧 P）。

    标准误为 0 的退化情形与 Welch 口径一致：系数为 0 时 t=0、P=1，
    否则 t 为同号 inf、P=0。
    """
    std_error = math.sqrt(variance) if variance > 0.0 else 0.0
    if std_error == 0.0:
        if effect == 0.0:
            return 0.0, 0.0, 1.0
        return 0.0, math.inf if effect > 0.0 else -math.inf, 0.0
    t_stat = effect / std_error
    return std_error, t_stat, t_sf_two_sided(t_stat, df)


def compute_numeric_covariate_cluster_pseudobulk_de(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: Sequence[int],
    replicate: ReplicateMetadata,
    covariates: PseudobulkCovariates,
    numeric: PseudobulkNumericCovariates,
) -> NumericCovariateClusterPseudobulkDE:
    """按最终簇逐簇汇总样本 pseudobulk 并对组别与数值协变量做独立检验。

    ``labels`` 与 ``qc.kept_cells`` 对齐（每个保留细胞一个最终簇标签）。
    样本次序、有效样本界定与簇准入（至少两个 group、每组至少两个有效样本）
    与簇内 pseudobulk 基线完全一致；簇与比较的枚举次序亦相同。
    """
    # 全样本（质控后有保留细胞）首次出现顺序；与簇内 pseudobulk 基线一致
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
    comparisons: List[NumericCovariateComparison] = []

    for cluster in sorted(set(labels)):
        active_samples, active_groups, values = _cluster_pseudobulk_values(
            matrix, qc, labels, replicate, cluster, sample_order
        )
        if not _eligible(active_groups):
            skipped_clusters.append(cluster)
            continue
        tested_clusters.append(cluster)

        group_to_obs: Dict[str, List[int]] = {}
        for i, group in enumerate(active_groups):
            group_to_obs.setdefault(group, []).append(i)
        present = sorted(group_to_obs)

        # 与基线相同的比较枚举：先全部 one-vs-rest（group 升序），
        # 再全部 (a, b)（a < b，group 升序）两两比较
        plan: List[Tuple[str, str, str, List[int], List[int]]] = []
        for group in present:
            obs_a = group_to_obs[group]
            obs_b = [i for i in range(len(active_groups)) if active_groups[i] != group]
            plan.append((ONE_VS_REST, group, "", obs_a, obs_b))
        for x, group_a in enumerate(present):
            for group_b in present[x + 1:]:
                plan.append(
                    (
                        PAIRWISE,
                        group_a,
                        group_b,
                        group_to_obs[group_a],
                        group_to_obs[group_b],
                    )
                )

        for comparison_type, group_a, group_b, obs_a, obs_b in plan:
            obs = list(obs_a) + list(obs_b)
            indicator = [1] * len(obs_a) + [0] * len(obs_b)
            design = _design_matrix(
                obs, indicator, active_samples, covariates, numeric
            )
            n_obs = len(design)
            n_params = len(design[0])
            df = n_obs - n_params

            gram = [
                [
                    sum(row[p] * row[q] for row in design)
                    for q in range(n_params)
                ]
                for p in range(n_params)
            ]
            inverse = _invert_matrix(gram)
            if inverse is None:
                raise CellFlowDataError(
                    f"簇 {cluster} 比较 {comparison_type} "
                    f"{group_a!r}/{group_b!r} 的设计矩阵秩不足"
                    f"（{n_obs} 个观测、{n_params} 个参数），"
                    f"无法进行数值协变量检验"
                )
            if df <= 0:
                raise CellFlowDataError(
                    f"簇 {cluster} 比较 {comparison_type} "
                    f"{group_a!r}/{group_b!r} 的残差自由度为 {df}"
                    f"（{n_obs} 个观测、{n_params} 个参数），"
                    f"无法进行数值协变量检验"
                )

            # 被检验项：组别项（第 2 列）与每个数值协变量连续列（表头顺序）
            n_dummies = n_params - 2 - len(numeric.covariate_names)
            terms: List[Tuple[str, int]] = [(GROUP_TERM, 1)]
            for offset, name in enumerate(numeric.covariate_names):
                terms.append((name, 2 + n_dummies + offset))

            records: List[NumericCovariateRecord] = []
            pvalues: List[float] = []
            for g, gene_id in enumerate(gene_ids):
                row = values[g]
                y = [row[i] for i in obs]
                beta, sigma2 = _ols_fit(design, inverse, y, df)
                for term_name, term_index in terms:
                    effect = beta[term_index]
                    std_error, t_stat, p_value = _term_test(
                        effect, sigma2 * inverse[term_index][term_index], df
                    )
                    pvalues.append(p_value)
                    records.append(
                        NumericCovariateRecord(
                            cluster=cluster,
                            comparison_type=comparison_type,
                            group_a=group_a,
                            group_b=group_b,
                            covariate=term_name,
                            gene_id=gene_id,
                            effect=effect,
                            std_error=std_error,
                            df=df,
                            t_stat=t_stat,
                            p_value=p_value,
                            p_value_adj=0.0,
                        )
                    )

            # 每个簇比较内跨全部（保留基因 × 被检验项）记录做 BH 校正
            adjusted = benjamini_hochberg(pvalues)
            records = [
                NumericCovariateRecord(
                    cluster=r.cluster,
                    comparison_type=r.comparison_type,
                    group_a=r.group_a,
                    group_b=r.group_b,
                    covariate=r.covariate,
                    gene_id=r.gene_id,
                    effect=r.effect,
                    std_error=r.std_error,
                    df=r.df,
                    t_stat=r.t_stat,
                    p_value=r.p_value,
                    p_value_adj=adjusted[g],
                )
                for g, r in enumerate(records)
            ]
            # 稳定排序：键相同（同 p_value_adj、同 effect、同 gene_id）时
            # 保持生成次序（基因外序、被检验项内序：组别在前、连续列按表头）
            records.sort(
                key=lambda r: (r.p_value_adj, -r.effect, r.gene_id)
            )
            comparisons.append(
                (cluster, comparison_type, group_a, group_b, records)
            )

    return NumericCovariateClusterPseudobulkDE(
        gene_ids=gene_ids,
        covariate_names=list(numeric.covariate_names),
        tested_clusters=tested_clusters,
        skipped_clusters=skipped_clusters,
        comparisons=comparisons,
    )
