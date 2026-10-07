"""最终簇内按样本协变量校正的 pseudobulk 差异表达。

仅在同时提供 ``--replicate-metadata`` 与 ``--pseudobulk-covariates``、
并启用 ``--cluster-pseudobulk-adjusted-de`` 无值开关（且已启用
``--cluster-pseudobulk-de``）时执行。计数保留范围、汇总方式与文库归一化
完全沿用簇内 pseudobulk 基线（:mod:`cell_flow.cluster_pseudobulk`）：
对每个最终簇与样本，只汇总该簇内通过质控的保留细胞在质控保留基因上的
原始计数；样本文库（该簇内该样本在保留基因上的计数和）归一到 10000 后
取 log1p。本模块是对 ``cluster_pseudobulk_de.tsv`` 的补充，不改写其
Welch 结果。

每个簇与每项比较（one-vs-rest 按 group 升序、再按 group 升序两两比较，
次序与既有簇内 pseudobulk 结果一致）只用该簇内的有效样本（该簇内含保留
细胞的样本；两两比较只用两个目标 group 的样本），以样本为观测拟合线性
模型：截距 + 组别项 + 每个协变量以其字典序最小水平为参照的哑变量。
组别系数即 ``log_fc_adjusted``，同时输出标准误、自由度、t 统计量与双侧
P 值，并在每个簇比较内跨保留基因做 Benjamini-Hochberg 校正。

设计矩阵秩不足或残差自由度不大于零时抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
"""

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

from .cluster_pseudobulk import _cluster_pseudobulk_values, _eligible
from .errors import CellFlowDataError
from .io import ExpressionMatrix
from .linalg import benjamini_hochberg, t_sf_two_sided
from .markers import ONE_VS_REST, PAIRWISE
from .pseudobulk_covariates import PseudobulkCovariates
from .qc import QCResult
from .replicate import ReplicateMetadata


@dataclass(frozen=True)
class AdjustedClusterPseudobulkRecord:
    """一个（簇, 比较, 基因）的协变量校正 pseudobulk 差异检验记录。"""

    cluster: int
    comparison_type: str
    group_a: str
    group_b: str                 # one-vs-rest 时为空
    n_samples_a: int
    n_samples_b: int
    gene_id: str
    mean_in_a: float
    mean_in_b: float
    log_fc_adjusted: float
    std_error: float
    df: int
    t_stat: float
    p_value: float
    p_value_adj: float


# (cluster, comparison_type, group_a, group_b, records)；
# one-vs-rest 的 group_b 为 ""
AdjustedClusterPseudobulkComparison = Tuple[
    int, str, str, str, List[AdjustedClusterPseudobulkRecord]
]


@dataclass(frozen=True)
class AdjustedClusterPseudobulkDE:
    """全部最终簇的协变量校正簇内 pseudobulk 差异结果。

    ``comparisons`` 按 cluster 升序，每簇内先全部 one-vs-rest（group 升序，
    ``group_b`` 为空）再全部 pairwise（group 升序），与既有簇内 pseudobulk
    结果的比较次序一致；每个比较内记录按校正 P 值升序、log_fc_adjusted
    降序、gene_id 升序排列。``tested_clusters`` 与 ``skipped_clusters``
    互补且各自升序，合起来即全部最终簇。
    """

    gene_ids: List[str]                       # 质控保留基因（保留行序）
    covariate_names: List[str]                # 参与校正的协变量（表头顺序）
    tested_clusters: List[int]
    skipped_clusters: List[int]
    comparisons: List[AdjustedClusterPseudobulkComparison]

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


def _is_singular(matrix: Sequence[Sequence[float]]) -> bool:
    """精确判定方阵奇异：Bareiss 无分式高斯消元（Fraction 精确算术）。

    设计矩阵只含 0/1 哑变量，Gram 矩阵元素为整数值的 float，
    ``Fraction`` 转换无损；浮点主元近乎为零但不恰为零的情形由此避免。
    """
    n = len(matrix)
    a = [[Fraction(value) for value in row] for row in matrix]
    previous = Fraction(1)
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = None
            for r in range(k + 1, n):
                if a[r][k] != 0:
                    pivot = r
                    break
            if pivot is None:
                return True
            a[k], a[pivot] = a[pivot], a[k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) / previous
        previous = a[k][k]
    return a[n - 1][n - 1] == 0


def _invert_matrix(matrix: Sequence[Sequence[float]]) -> Optional[List[List[float]]]:
    """Gauss-Jordan 部分主元求逆；奇异（经精确判定）时返回 None。"""
    if _is_singular(matrix):
        return None
    n = len(matrix)
    a = [
        list(matrix[i]) + [1.0 if i == j else 0.0 for j in range(n)]
        for i in range(n)
    ]
    for col in range(n):
        pivot = col
        for r in range(col + 1, n):
            if abs(a[r][col]) > abs(a[pivot][col]):
                pivot = r
        if a[pivot][col] == 0.0:
            return None
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
        pivot_value = a[col][col]
        row = a[col]
        for j in range(2 * n):
            row[j] /= pivot_value
        for r in range(n):
            if r == col:
                continue
            factor = a[r][col]
            if factor == 0.0:
                continue
            other = a[r]
            for j in range(2 * n):
                other[j] -= factor * row[j]
    return [row[n:] for row in a]


def _design_matrix(
    obs: Sequence[int],
    obs_group_indicator: Sequence[int],
    active_samples: Sequence[str],
    covariates: PseudobulkCovariates,
) -> List[List[float]]:
    """构造设计矩阵：截距、组别项，再按表头顺序展开各协变量哑变量。

    每个协变量的水平取本比较实际使用样本中的取值，字典序最小水平为参照，
    其余每个水平一列（按字典序）；只出现一个水平的协变量不产生哑变量列。
    """
    # 每个协变量的非参照水平（哑变量列），按表头顺序、水平字典序展开
    dummies: List[Tuple[str, str]] = []
    for name in covariates.covariate_names:
        present = sorted(
            {covariates.sample_levels[active_samples[j]][name] for j in obs}
        )
        for level in present[1:]:
            dummies.append((name, level))

    rows: List[List[float]] = []
    for position, i in enumerate(obs):
        levels = covariates.sample_levels[active_samples[i]]
        row = [1.0, float(obs_group_indicator[position])]
        for name, level in dummies:
            row.append(1.0 if levels[name] == level else 0.0)
        rows.append(row)
    return rows


def _ols_group_effect(
    design: List[List[float]],
    inverse: List[List[float]],
    y: Sequence[float],
    df: int,
) -> Tuple[float, float, float, float]:
    """对单个基因拟合 OLS，返回（组别系数, 标准误, t, 双侧 P）。

    ``inverse`` 为 (XᵀX)⁻¹；组别项固定为第 2 列（下标 1）。标准误为 0 的
    退化情形与 Welch 口径一致：系数为 0 时 t=0、P=1，否则 t 为同号 inf、
    P=0。
    """
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
    sigma2 = rss / df
    effect = beta[1]
    variance = sigma2 * inverse[1][1]
    std_error = math.sqrt(variance) if variance > 0.0 else 0.0
    if std_error == 0.0:
        if effect == 0.0:
            return effect, 0.0, 0.0, 1.0
        return (
            effect,
            0.0,
            math.inf if effect > 0.0 else -math.inf,
            0.0,
        )
    t_stat = effect / std_error
    return effect, std_error, t_stat, t_sf_two_sided(t_stat, df)


def compute_adjusted_cluster_pseudobulk_de(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: Sequence[int],
    replicate: ReplicateMetadata,
    covariates: PseudobulkCovariates,
) -> AdjustedClusterPseudobulkDE:
    """按最终簇逐簇汇总样本 pseudobulk 并做协变量校正的簇内差异表达。

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
    comparisons: List[AdjustedClusterPseudobulkComparison] = []

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
            design = _design_matrix(obs, indicator, active_samples, covariates)
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
                    f"无法进行协变量校正的差异表达"
                )
            if df <= 0:
                raise CellFlowDataError(
                    f"簇 {cluster} 比较 {comparison_type} "
                    f"{group_a!r}/{group_b!r} 的残差自由度为 {df}"
                    f"（{n_obs} 个观测、{n_params} 个参数），"
                    f"无法进行协变量校正的差异表达"
                )

            n_a = len(obs_a)
            n_b = len(obs_b)
            records: List[AdjustedClusterPseudobulkRecord] = []
            pvalues: List[float] = []
            for g, gene_id in enumerate(gene_ids):
                row = values[g]
                y = [row[i] for i in obs]
                mean_a = sum(y[:n_a]) / n_a
                mean_b = sum(y[n_a:]) / n_b
                effect, std_error, t_stat, p_value = _ols_group_effect(
                    design, inverse, y, df
                )
                pvalues.append(p_value)
                records.append(
                    AdjustedClusterPseudobulkRecord(
                        cluster=cluster,
                        comparison_type=comparison_type,
                        group_a=group_a,
                        group_b=group_b,
                        n_samples_a=n_a,
                        n_samples_b=n_b,
                        gene_id=gene_id,
                        mean_in_a=mean_a,
                        mean_in_b=mean_b,
                        log_fc_adjusted=effect,
                        std_error=std_error,
                        df=df,
                        t_stat=t_stat,
                        p_value=p_value,
                        p_value_adj=0.0,
                    )
                )

            adjusted = benjamini_hochberg(pvalues)
            records = [
                AdjustedClusterPseudobulkRecord(
                    cluster=r.cluster,
                    comparison_type=r.comparison_type,
                    group_a=r.group_a,
                    group_b=r.group_b,
                    n_samples_a=r.n_samples_a,
                    n_samples_b=r.n_samples_b,
                    gene_id=r.gene_id,
                    mean_in_a=r.mean_in_a,
                    mean_in_b=r.mean_in_b,
                    log_fc_adjusted=r.log_fc_adjusted,
                    std_error=r.std_error,
                    df=r.df,
                    t_stat=r.t_stat,
                    p_value=r.p_value,
                    p_value_adj=adjusted[g],
                )
                for g, r in enumerate(records)
            ]
            records.sort(
                key=lambda r: (r.p_value_adj, -r.log_fc_adjusted, r.gene_id)
            )
            comparisons.append(
                (cluster, comparison_type, group_a, group_b, records)
            )

    return AdjustedClusterPseudobulkDE(
        gene_ids=gene_ids,
        covariate_names=list(covariates.covariate_names),
        tested_clusters=tested_clusters,
        skipped_clusters=skipped_clusters,
        comparisons=comparisons,
    )
