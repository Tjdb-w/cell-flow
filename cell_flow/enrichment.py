"""Marker 基因集富集：在逐簇 one-versus-rest 差异表达命中上做超几何检验。

每簇命中基因（marker hit）为 ``markers.tsv`` 同一口径的 Welch t 检验结果中
``p_value_adj <= alpha`` 且 ``log_fc >= min_log_fc`` 的保留基因；背景为
全部质控后保留基因。对每个 ``set_id``：

- ``n_set_total``：集合在基因集文件中的成员总数（含矩阵外、QC 剔除基因）；
- ``n_set_used``：集合与背景（保留基因）的交集大小；
- ``n_markers``：该簇命中数；
- ``n_overlap``：命中且属于集合的基因数。

单侧 P 值为超几何分布的生存函数
``P(X >= n_overlap; N=背景, K=n_set_used, n=n_markers)``；
每个簇内跨集合做 Benjamini-Hochberg 校正。无重叠（``n_overlap == 0``）时
``p_value = 1``、``fold_enrichment = 0``。

``expected_overlap = n_markers * n_set_used / N``；
``fold_enrichment = n_overlap / expected_overlap``，除零取 0；
``odds_ratio`` 按“命中 × 是否属于集合”的 2×2 表四格各加 0.5 计算。
"""

import math
from dataclasses import dataclass
from typing import Dict, List

from .gene_sets import GeneSets
from .linalg import benjamini_hochberg
from .markers import MarkerRecord

TOP_N_ENRICHMENT = 20


@dataclass(frozen=True)
class EnrichmentRecord:
    cluster: int
    set_id: str
    n_set_total: int
    n_set_used: int
    n_markers: int
    n_overlap: int
    expected_overlap: float
    fold_enrichment: float
    odds_ratio: float
    p_value: float
    p_value_adj: float


@dataclass(frozen=True)
class MarkerEnrichment:
    """全部簇 × 全部集合的富集结果。

    ``rows`` 按 cluster 升序、p_value_adj 升序、set_id 升序，
    odds_ratio 降序（末位仅在三者并列时生效）排列；
    ``set_order`` 为升序 set_id 列表；``cluster_marker_counts`` 为
    每簇命中基因数（n_markers）。
    """

    rows: List[EnrichmentRecord]
    set_order: List[str]
    cluster_marker_counts: Dict[int, int]


def _hypergeom_sf(k: int, n_population: int, n_success: int, n_draws: int) -> float:
    """超几何单侧生存函数 P(X >= k)。

    总体 ``n_population`` 个，其中成功 ``n_success`` 个，不放回抽
    ``n_draws`` 个，观测成功数不少于 ``k`` 的概率。以对数组合数与
    max-log-sum-exp 计算，纯 float64，结果确定。
    """
    n = n_population
    k_pop = n_success
    m = n_draws
    if k <= 0:
        return 1.0
    upper = min(k_pop, m)
    if k > upper:
        return 0.0
    lower = max(0, m - (n - k_pop))
    if k < lower:
        # 观测低于理论下限时，全部支撑均 >= k
        return 1.0

    def log_pmf(j: int) -> float:
        # log C(K,j) + log C(N-K, M-j) - log C(N, M)
        return (
            math.lgamma(k_pop + 1)
            - math.lgamma(j + 1)
            - math.lgamma(k_pop - j + 1)
            + math.lgamma(n - k_pop + 1)
            - math.lgamma(m - j + 1)
            - math.lgamma(n - k_pop - m + j + 1)
            - math.lgamma(n + 1)
            + math.lgamma(m + 1)
            + math.lgamma(n - m + 1)
        )

    log_terms = [log_pmf(j) for j in range(k, upper + 1)]
    max_log = max(log_terms)
    total = sum(math.exp(term - max_log) for term in log_terms)
    p = math.exp(max_log) * total
    if p < 0.0:
        return 0.0
    if p > 1.0:
        return 1.0
    return p


def enrich_marker_gene_sets(
    gene_sets: GeneSets,
    markers: Dict[int, List[MarkerRecord]],
    *,
    kept_gene_ids: List[str],
    alpha: float,
    min_log_fc: float,
) -> MarkerEnrichment:
    """计算逐簇 marker 命中对各基因集的超几何富集。

    命中口径与 ``markers.tsv`` 的 one-versus-rest Welch t 检验完全一致：
    ``p_value_adj <= alpha`` 且 ``log_fc >= min_log_fc``。背景为
    ``kept_gene_ids``（质控后保留基因）。集合成员与背景精确匹配，
    不做别名转换；背景外交员计入 ``n_set_total``。
    """
    background = set(kept_gene_ids)
    n_background = len(kept_gene_ids)

    set_order = sorted(gene_sets.members)
    set_used: Dict[str, set] = {}
    set_total: Dict[str, int] = {}
    for set_id in set_order:
        members = gene_sets.members[set_id]
        set_total[set_id] = len(members)
        set_used[set_id] = {gene_id for gene_id in members if gene_id in background}

    rows: List[EnrichmentRecord] = []
    cluster_marker_counts: Dict[int, int] = {}
    for cluster in sorted(markers):
        hit_genes = {
            r.gene_id
            for r in markers[cluster]
            if r.p_value_adj <= alpha and r.log_fc >= min_log_fc
        }
        n_markers = len(hit_genes)
        cluster_marker_counts[cluster] = n_markers

        cluster_rows: List[EnrichmentRecord] = []
        pvalues: List[float] = []
        for set_id in set_order:
            used = set_used[set_id]
            n_set_used = len(used)
            overlap = len(hit_genes & used)
            expected = n_markers * n_set_used / n_background
            fold = overlap / expected if expected > 0.0 else 0.0
            if overlap == 0:
                # 无命中：生存函数恒为 1，fold 恒为 0
                p_value = 1.0
            else:
                p_value = _hypergeom_sf(
                    overlap, n_background, n_set_used, n_markers
                )
            pvalues.append(p_value)

            # “命中 × 是否属于集合”的 2×2 表，四格统一加 0.5
            a = overlap
            b = n_markers - overlap
            c = n_set_used - overlap
            d = n_background - n_markers - n_set_used + overlap
            odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

            cluster_rows.append(
                EnrichmentRecord(
                    cluster=cluster,
                    set_id=set_id,
                    n_set_total=set_total[set_id],
                    n_set_used=n_set_used,
                    n_markers=n_markers,
                    n_overlap=overlap,
                    expected_overlap=expected,
                    fold_enrichment=fold,
                    odds_ratio=odds_ratio,
                    p_value=p_value,
                    p_value_adj=0.0,
                )
            )

        adjusted = benjamini_hochberg(pvalues)
        cluster_rows = [
            EnrichmentRecord(
                cluster=r.cluster,
                set_id=r.set_id,
                n_set_total=r.n_set_total,
                n_set_used=r.n_set_used,
                n_markers=r.n_markers,
                n_overlap=r.n_overlap,
                expected_overlap=r.expected_overlap,
                fold_enrichment=r.fold_enrichment,
                odds_ratio=r.odds_ratio,
                p_value=r.p_value,
                p_value_adj=adjusted[i],
            )
            for i, r in enumerate(cluster_rows)
        ]
        rows.extend(cluster_rows)

    # 排序：cluster、p_value_adj、set_id 均升序；odds_ratio 降序作为末位
    # （set_id 在同簇内唯一，该键仅在三者完全并列时参与，保持确定）。
    # 两次稳定排序实现混合升降序：先 odds_ratio 降序，再主键升序。
    rows.sort(key=lambda r: -r.odds_ratio)
    rows.sort(key=lambda r: (r.cluster, r.p_value_adj, r.set_id))
    return MarkerEnrichment(
        rows=rows,
        set_order=set_order,
        cluster_marker_counts=cluster_marker_counts,
    )
