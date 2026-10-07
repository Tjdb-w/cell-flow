"""最终簇内含数值型样本协变量的 pseudobulk 差异检验。

仅在同时提供 ``--replicate-metadata``、启用 ``--cluster-pseudobulk-de``
与 ``--cluster-pseudobulk-adjusted-de`` 并提供
``--pseudobulk-numeric-covariates`` 时执行（可再并用
``--pseudobulk-covariates`` 分类协变量）。计数保留范围、汇总方式与文库
归一化完全沿用簇内 pseudobulk 基线（:mod:`cell_flow.cluster_pseudobulk`），
是对 ``cluster_pseudobulk_de.tsv`` 与
``cluster_pseudobulk_adjusted_de.tsv`` 的独立补充，不改写既有结果。

每个簇与每项比较（one-vs-rest 按 group 升序、再按 group 升序两两比较，
次序与既有簇内 pseudobulk 结果一致）只用该簇内的有效样本（两两比较只用
两个目标 group 的样本），以样本为观测拟合一个线性模型：截距 + 组别项 +
每个（可选的）分类协变量以其字典序最小水平为参照的哑变量 + 每个数值协变量
原值连续列。对簇内 pseudobulk log1p 的每个保留基因，输出组别项与每个
数值协变量连续列的系数（``effect``）、标准误、自由度、t 统计量与双侧
P 值；每个簇比较内把“组别 + 全部连续列”跨保留基因的全部 P 值作为一个
校正家族做 Benjamini-Hochberg 校正。

设计矩阵秩不足或残差自由度不大于零时抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .adjusted_cluster_pseudobulk import _invert_matrix
from .cluster_pseudobulk import _cluster_pseudobulk_values, _eligible
from .errors import CellFlowDataError
from .io import ExpressionMatrix
from .linalg import benjamini_hochberg, t_sf_two_sided
from .markers import ONE_VS_REST, PAIRWISE
from .numeric_pseudobulk_covariates import PseudobulkNumericCovariates
from .pseudobulk_covariates import PseudobulkCovariates
from .qc import QCResult
from .replicate import ReplicateMetadata

#: 组别项在输出 ``covariate`` 列中的标识
GROUP_COVARIATE = "group"


@dataclass(frozen=True)
class NumericClusterPseudobulkRecord:
    """一个（簇, 比较, 检验项, 基因）的数值协变量 pseudobulk 检验记录。"""

    cluster: int
    comparison_type: str
    group_a: str
    group_b: str                 # one-vs-rest 时为空
    covariate: str               # 组别项为 "group"，否则为数值协变量列名
    gene_id: str
    effect: float
    std_error: float
    df: int
    t_stat: float
    p_value: float
    p_value_adj: float


# (cluster, comparison_type, group_a, group_b, records)；
# one-vs-rest 的 group_b 为 ""；records 内先组别项再各连续列（表头顺序），
# 每个检验项内按校正 P 值升序、effect 降序、gene_id 升序排列
NumericClusterPseudobulkComparison = Tuple[
    int, str, str, str, List[NumericClusterPseudobulkRecord]
]


@dataclass(frozen=True)
class NumericClusterPseudobulkDE:
    """全部最终簇的含数值协变量簇内 pseudobulk 检验结果。

    ``comparisons`` 按 cluster 升序，每簇内先全部 one-vs-rest（group 升序，
    ``group_b`` 为空）再全部 pairwise（group 升序），与既有簇内 pseudobulk
    结果的比较次序一致。``tested_clusters`` 与 ``skipped_clusters`` 互补且
    各自升序，合起来即全部最终簇。
    """

    gene_ids: List[str]                       # 质控保留基因（保留行序）
    covariate_names: List[str]                # 数值协变量列名（表头顺序）
    tested_clusters: List[int]
    skipped_clusters: List[int]
    comparisons: List[NumericClusterPseudobulkComparison]

    @property
    def comparison_count(self) -> int:
        return len(self.comparisons)

    @property
    def test_count(self) -> int:
        """全部（簇, 比较, 检验项, 基因）检验数。"""
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
    categorical: Optional[PseudobulkCovariates],
    numeric: PseudobulkNumericCovariates,
) -> Tuple[List[List[float]], List[Tuple[str, int]]]:
    """构造设计矩阵并返回（矩阵, 受检连续项列表）。

    列依次为：截距、组别项、各分类协变量哑变量（仅当并给分类协变量表时，
    口径同 :mod:`cell_flow.adjusted_cluster_pseudobulk`）、各数值协变量原值
    连续列（数值协变量表头顺序）。受检项固定先组别项（下标 1）再各连续列。
    """
    # 每个分类协变量的非参照水平（哑变量列），按表头顺序、水平字典序展开
    dummies: List[Tuple[str, str]] = []
    if categorical is not None:
        for name in categorical.covariate_names:
            present = sorted(
                {categorical.sample_levels[active_samples[j]][name] for j in obs}
            )
            for level in present[1:]:
                dummies.append((name, level))

    tested_columns: List[Tuple[str, int]] = [(GROUP_COVARIATE, 1)]
    rows: List[List[float]] = []
    for position, i in enumerate(obs):
        sample_id = active_samples[i]
        row = [1.0, float(obs_group_indicator[position])]
        if categorical is not None:
            levels = categorical.sample_levels[sample_id]
            for name, level in dummies:
                row.append(1.0 if levels[name] == level else 0.0)
        rows.append(row)

    # 连续列追加在哑变量之后，记录其列下标供逐基因检验
    next_index = 2 + len(dummies)
    for position, name in enumerate(numeric.covariate_names):
        column_index = next_index + position
        tested_columns.append((name, column_index))
        for row_position, i in enumerate(obs):
            sample_id = active_samples[i]
            rows[row_position].append(numeric.sample_values[sample_id][name])

    return rows, tested_columns


def _ols_coefficients(
    design: List[List[float]],
    inverse: List[List[float]],
    y: Sequence[float],
) -> Tuple[List[float], float]:
    """拟合一次 OLS，返回（全部系数 beta, 残差平方和 rss）。"""
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
    return beta, rss


def _coefficient_test(
    effect: float,
    inverse_diag: float,
    sigma2: float,
    df: int,
) -> Tuple[float, float, float]:
    """由系数、(XᵀX)⁻¹ 对角元与残差方差返回（标准误, t, 双侧 P）。

    标准误为 0 的退化情形与既有 OLS/Welch 口径一致：系数为 0 时 t=0、P=1，
    否则 t 为同号 inf、P=0。
    """
    variance = sigma2 * inverse_diag
    std_error = math.sqrt(variance) if variance > 0.0 else 0.0
    if std_error == 0.0:
        if effect == 0.0:
            return 0.0, 0.0, 1.0
        return 0.0, math.inf if effect > 0.0 else -math.inf, 0.0
    t_stat = effect / std_error
    return std_error, t_stat, t_sf_two_sided(t_stat, df)


def compute_numeric_cluster_pseudobulk_de(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: Sequence[int],
    replicate: ReplicateMetadata,
    numeric: PseudobulkNumericCovariates,
    categorical: Optional[PseudobulkCovariates] = None,
) -> NumericClusterPseudobulkDE:
    """按最终簇逐簇汇总样本 pseudobulk 并做含数值协变量的簇内检验。

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
    comparisons: List[NumericClusterPseudobulkComparison] = []

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
            obs_b = [
                i for i in range(len(active_groups)) if active_groups[i] != group
            ]
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
            design, tested_columns = _design_matrix(
                obs, indicator, active_samples, categorical, numeric
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
                    f"无法进行含数值协变量的 pseudobulk 差异检验"
                )
            if df <= 0:
                raise CellFlowDataError(
                    f"簇 {cluster} 比较 {comparison_type} "
                    f"{group_a!r}/{group_b!r} 的残差自由度为 {df}"
                    f"（{n_obs} 个观测、{n_params} 个参数），"
                    f"无法进行含数值协变量的 pseudobulk 差异检验"
                )

            # 每个检验项（组别 + 各连续列）一列、每列按保留基因行序放 P 值；
            # 展平下标 = 检验项序号 * 基因数 + 基因序号，BH 家族为整项比较
            n_genes = len(gene_ids)
            n_terms = len(tested_columns)
            pvalue_table: List[List[float]] = [
                [0.0] * n_genes for _ in range(n_terms)
            ]
            # 暂存每检验项每基因的 (effect, se, t, p) 供建记录
            stats_table: List[List[Tuple[float, float, float, float]]] = [
                [] for _ in range(n_terms)
            ]
            for g in range(n_genes):
                y = [values[g][i] for i in obs]
                beta, rss = _ols_coefficients(design, inverse, y)
                sigma2 = rss / df
                for term_index, (_, column_index) in enumerate(tested_columns):
                    effect = beta[column_index]
                    std_error, t_stat, p_value = _coefficient_test(
                        effect, inverse[column_index][column_index], sigma2, df
                    )
                    pvalue_table[term_index][g] = p_value
                    stats_table[term_index].append(
                        (effect, std_error, t_stat, p_value)
                    )

            flat_pvalues = [
                pvalue_table[term_index][g]
                for term_index in range(n_terms)
                for g in range(n_genes)
            ]
            flat_adjusted = benjamini_hochberg(flat_pvalues)

            comparison_records: List[NumericClusterPseudobulkRecord] = []
            for term_index, (covariate_name, _) in enumerate(tested_columns):
                block: List[NumericClusterPseudobulkRecord] = []
                for g, gene_id in enumerate(gene_ids):
                    effect, std_error, t_stat, p_value = stats_table[
                        term_index
                    ][g]
                    p_value_adj = flat_adjusted[term_index * n_genes + g]
                    block.append(
                        NumericClusterPseudobulkRecord(
                            cluster=cluster,
                            comparison_type=comparison_type,
                            group_a=group_a,
                            group_b=group_b,
                            covariate=covariate_name,
                            gene_id=gene_id,
                            effect=effect,
                            std_error=std_error,
                            df=df,
                            t_stat=t_stat,
                            p_value=p_value,
                            p_value_adj=p_value_adj,
                        )
                    )
                block.sort(key=lambda r: (r.p_value_adj, -r.effect, r.gene_id))
                comparison_records.extend(block)
            comparisons.append(
                (cluster, comparison_type, group_a, group_b, comparison_records)
            )

    return NumericClusterPseudobulkDE(
        gene_ids=gene_ids,
        covariate_names=list(numeric.covariate_names),
        tested_clusters=tested_clusters,
        skipped_clusters=skipped_clusters,
        comparisons=comparisons,
    )
