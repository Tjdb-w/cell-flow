"""最终簇内校正样本协变量的 pseudobulk 差异表达。

仅在提供 ``--replicate-metadata``、启用 ``--cluster-pseudobulk-de``、
提供 ``--pseudobulk-covariates`` 并启用 ``--cluster-pseudobulk-adjusted-de``
无值开关时执行，作为 ``cluster_pseudobulk_de.tsv`` 的补充（不改写其 Welch
结果）。计数汇总、保留基因范围、有效样本界定与样本文库归一到 10000 后取
log1p 的值完全沿用簇内 pseudobulk 基线（:mod:`cell_flow.cluster_pseudobulk`）。

每个簇与每项比较只用该簇内的有效样本（one-vs-rest 用全部有效样本，
两两比较只用两个 group 的有效样本）。对每个基因拟合普通最小二乘模型：
截距 + 组别项 + 各协变量的哑变量（每个协变量以参与该模型的样本中字典序
最小水平为参照，其余每个水平一列，协变量按离散标签处理）。组别系数即
``log_fc_adjusted``；同时输出标准误、残差自由度（有效样本数减设计矩阵
列数）、t 统计量与双侧 P 值，并在每个簇的每项比较内跨保留基因做
Benjamini-Hochberg 校正。设计矩阵秩不足或残差自由度不大于零时抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
"""

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from .cluster_pseudobulk import _cluster_pseudobulk_values, _eligible
from .errors import CellFlowDataError
from .io import ExpressionMatrix
from .linalg import benjamini_hochberg, t_sf_two_sided
from .markers import ONE_VS_REST, PAIRWISE
from .pseudobulk_covariates import PseudobulkCovariates
from .qc import QCResult
from .replicate import ReplicateMetadata

# 主元小于该相对阈值即视为设计矩阵秩不足（X'X 为整数矩阵，合法主元不会
# 落到该量级；阈值同时吸收浮点消元残差）
_PIVOT_EPS = 1e-12


@dataclass(frozen=True)
class AdjustedClusterPseudobulkRecord:
    """一个（簇, 比较, 基因）的协变量校正 pseudobulk 差异检验记录。"""

    cluster: int
    comparison_type: str
    group_a: str
    group_b: str                 # one-vs-rest 时为空
    gene_id: str
    n_samples: int               # 参与该模型的有效样本数
    mean_in_a: float
    mean_in_b: float
    log_fc_adjusted: float       # 组别项的最小二乘系数
    std_error: float             # 组别系数的标准误
    df: int                      # 残差自由度
    t_stat: float
    p_value: float
    p_value_adj: float


# (cluster, comparison_type, group_a, group_b, records)；
# one-vs-rest 的 group_b 为 ""
AdjustedClusterPseudobulkComparison = Tuple[
    int, str, str, str, List[AdjustedClusterPseudobulkRecord]
]


@dataclass(frozen=True)
class ClusterPseudobulkAdjustedDE:
    """全部最终簇的协变量校正簇内 pseudobulk 差异结果。

    ``comparisons`` 的簇与比较次序与簇内 pseudobulk 基线一致（cluster 升序，
    每簇内先全部 one-vs-rest 再全部 pairwise）；每个比较内记录按校正 P 值
    升序、log_fc_adjusted 降序、gene_id 升序排列。``tested_clusters`` 与
    ``skipped_clusters`` 互补且各自升序，合起来即全部最终簇。
    """

    gene_ids: List[str]                       # 质控保留基因（保留行序）
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


def _invert(matrix: List[List[float]]) -> Optional[List[List[float]]]:
    """Gauss-Jordan 求逆（部分主元）；秩不足时返回 ``None``。

    纯 float64 且主元选取确定（列内首个最大绝对值），结果确定。
    """
    n = len(matrix)
    a = [
        list(row) + [1.0 if i == j else 0.0 for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    scale = max(abs(value) for row in matrix for value in row)
    if scale <= 0.0:
        return None
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) <= _PIVOT_EPS * scale:
            return None
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
        pivot_row = [value / a[col][col] for value in a[col]]
        a[col] = pivot_row
        for r in range(n):
            if r == col:
                continue
            factor = a[r][col]
            if factor == 0.0:
                continue
            a[r] = [
                value - factor * pivot_value
                for value, pivot_value in zip(a[r], pivot_row)
            ]
    return [row[n:] for row in a]


def _build_design(
    samples: Sequence[str],
    group_flags: Sequence[int],
    covariates: PseudobulkCovariates,
) -> List[List[float]]:
    """构造设计矩阵：截距、组别项、各协变量哑变量（字典序最小水平为参照）。

    每个协变量的水平取自参与本模型的样本；只有一个水平时不产生哑变量列。
    """
    n = len(samples)
    columns: List[List[float]] = [
        [1.0] * n,
        [float(flag) for flag in group_flags],
    ]
    for name in covariates.covariate_columns:
        levels = sorted({covariates.values[s][name] for s in samples})
        for level in levels[1:]:
            columns.append(
                [1.0 if covariates.values[s][name] == level else 0.0
                 for s in samples]
            )
    # 转置为逐样本行
    return [[column[i] for column in columns] for i in range(n)]


def _fit_gene(
    design: List[List[float]],
    inverse_xtx: List[List[float]],
    y: Sequence[float],
    df: int,
) -> Tuple[float, float, float, float]:
    """单基因 OLS：返回（组别系数, 标准误, t 统计量, 双侧 P 值）。

    组别项固定为设计矩阵第 2 列（下标 1）。零残差方差退化：系数为 0 时
    t=0、p=1；系数非零时 t 为与系数同号的 inf、p=0。
    """
    p = len(design[0])
    xty = [
        sum(row[j] * y[i] for i, row in enumerate(design)) for j in range(p)
    ]
    beta = [
        sum(inverse_xtx[j][k] * xty[k] for k in range(p)) for j in range(p)
    ]
    rss = 0.0
    for i, row in enumerate(design):
        fitted = sum(row[j] * beta[j] for j in range(p))
        residual = y[i] - fitted
        rss += residual * residual
    sigma2 = rss / df
    var_beta = sigma2 * inverse_xtx[1][1]
    coefficient = beta[1]
    if var_beta <= 0.0:
        if coefficient == 0.0:
            return coefficient, 0.0, 0.0, 1.0
        return (
            coefficient,
            0.0,
            math.inf if coefficient > 0.0 else -math.inf,
            0.0,
        )
    std_error = math.sqrt(var_beta)
    t_stat = coefficient / std_error
    return coefficient, std_error, t_stat, t_sf_two_sided(t_stat, df)


def compute_cluster_pseudobulk_adjusted_de(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: Sequence[int],
    replicate: ReplicateMetadata,
    covariates: PseudobulkCovariates,
    sample_order: Sequence[str],
) -> ClusterPseudobulkAdjustedDE:
    """按最终簇逐簇汇总样本 pseudobulk 并做协变量校正的分组差异表达。

    ``labels`` 与 ``qc.kept_cells`` 对齐；``sample_order`` 为质控后有保留
    细胞的全部样本（pseudobulk 基线样本序），协变量表已与之覆盖一致。
    簇准入口径与簇内 pseudobulk 基线一致；进入分析的簇若任一比较的设计
    矩阵秩不足或残差自由度不大于零，抛 :class:`CellFlowDataError`。
    """
    gene_ids = [matrix.gene_ids[g] for g in qc.kept_genes]

    tested_clusters: List[int] = []
    skipped_clusters: List[int] = []
    comparisons: List[AdjustedClusterPseudobulkComparison] = []

    for cluster in sorted(set(labels)):
        active_samples, active_groups, values = _cluster_pseudobulk_values(
            matrix, qc, labels, replicate, cluster, list(sample_order)
        )
        if not _eligible(active_groups):
            skipped_clusters.append(cluster)
            continue
        tested_clusters.append(cluster)

        present = sorted(set(active_groups))
        # 与基线相同的比较次序：先全部 one-vs-rest（group 升序），
        # 再全部 (a, b)（a < b，按 group 升序）两两比较
        plan: List[Tuple[str, str, str, List[int]]] = []
        for group in present:
            plan.append(
                (ONE_VS_REST, group, "", list(range(len(active_samples))))
            )
        for i, group_a in enumerate(present):
            for group_b in present[i + 1:]:
                obs = [
                    s for s, group in enumerate(active_groups)
                    if group == group_a or group == group_b
                ]
                plan.append((PAIRWISE, group_a, group_b, obs))

        for comparison_type, group_a, group_b, obs in plan:
            samples = [active_samples[s] for s in obs]
            group_flags = [
                1 if active_groups[s] == group_a else 0 for s in obs
            ]
            design = _build_design(samples, group_flags, covariates)
            n_samples = len(samples)
            n_columns = len(design[0])
            df = n_samples - n_columns
            if df <= 0:
                raise CellFlowDataError(
                    f"簇 {cluster} 比较 {group_a!r} vs {group_b!r} 的校正模型"
                    f"残差自由度为 {df}（有效样本 {n_samples} 个、设计矩阵 "
                    f"{n_columns} 列），不大于零，无法进行协变量校正差异表达"
                )
            xtx = [
                [
                    sum(row[j] * row[k] for row in design)
                    for k in range(n_columns)
                ]
                for j in range(n_columns)
            ]
            inverse_xtx = _invert(xtx)
            if inverse_xtx is None:
                raise CellFlowDataError(
                    f"簇 {cluster} 比较 {group_a!r} vs {group_b!r} 的校正模型"
                    f"设计矩阵秩不足（有效样本 {n_samples} 个、设计矩阵 "
                    f"{n_columns} 列），无法进行协变量校正差异表达"
                )

            index_a = [i for i, flag in enumerate(group_flags) if flag == 1]
            index_b = [i for i, flag in enumerate(group_flags) if flag == 0]

            records: List[AdjustedClusterPseudobulkRecord] = []
            pvalues: List[float] = []
            for g, gene_id in enumerate(gene_ids):
                row = values[g]
                y = [row[s] for s in obs]
                mean_a = sum(y[i] for i in index_a) / len(index_a)
                mean_b = sum(y[i] for i in index_b) / len(index_b)
                coefficient, std_error, t_stat, p_value = _fit_gene(
                    design, inverse_xtx, y, df
                )
                pvalues.append(p_value)
                records.append(
                    AdjustedClusterPseudobulkRecord(
                        cluster=cluster,
                        comparison_type=comparison_type,
                        group_a=group_a,
                        group_b=group_b,
                        gene_id=gene_id,
                        n_samples=n_samples,
                        mean_in_a=mean_a,
                        mean_in_b=mean_b,
                        log_fc_adjusted=coefficient,
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
                    gene_id=r.gene_id,
                    n_samples=r.n_samples,
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

    return ClusterPseudobulkAdjustedDE(
        gene_ids=gene_ids,
        tested_clusters=tested_clusters,
        skipped_clusters=skipped_clusters,
        comparisons=comparisons,
    )
