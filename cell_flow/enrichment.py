"""marker 基因集富集：逐簇显著 marker 对基因集的超几何富集检验。

仅当显式启用 ``--enrich-markers``（须同时提供 ``--gene-sets``）时计算。
沿用最终细胞、保留基因、最终簇与 ``markers.tsv`` 的 one-versus-rest
Welch t 检验结果：每簇命中基因为 ``p_value_adj <= alpha`` 且
``log_fc >= min_log_fc`` 的基因，背景为全部保留基因。每个集合取集合与
背景的交集（``n_set_used``），以命中数不小于观测值的超几何单侧 P 值
度量富集，并在每簇内跨集合做 BH 校正；矩阵外与被 QC 剔除的集合成员
只计入 ``n_set_total``。
"""

import math
from dataclasses import dataclass
from typing import Dict, List

from .gene_sets import GeneSets
from .linalg import benjamini_hochberg
from .markers import MarkerRecord

DEFAULT_ENRICHMENT_ALPHA = 0.05
DEFAULT_ENRICHMENT_MIN_LOG_FC = 0.0


@dataclass(frozen=True)
class EnrichmentRecord:
    """一个（簇, 基因集）的富集检验结果。"""

    cluster: int
    set_id: str
    n_set_total: int        # 集合成员总数（含矩阵外与 QC 剔除基因）
    n_set_used: int         # 集合与背景（保留基因）的交集大小
    n_markers: int          # 该簇命中基因数
    n_overlap: int          # 命中基因中属于该集合的数量
    expected_overlap: float
    fold_enrichment: float
    odds_ratio: float
    p_value: float
    p_value_adj: float


@dataclass(frozen=True)
class EnrichmentResult:
    """全部（簇, 基因集）富集结果。

    ``rows`` 按 cluster、p_value_adj、set_id 升序、odds_ratio 降序排列；
    ``n_marker_hits`` 为各簇命中基因数之和；``n_significant_sets`` 为
    校正后仍显著（p_value_adj <= alpha）的（簇, 集合）行数。
    """

    rows: List[EnrichmentRecord]
    n_marker_hits: int
    n_significant_sets: int


def _log_choose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeometric_sf(
    observed: int, n_background: int, n_set: int, n_drawn: int
) -> float:
    """超几何分布单侧生存函数 P(X >= observed)。

    X ~ Hypergeometric(总体 n_background，成功态 n_set，抽取 n_drawn)。
    对数空间逐点求和（math.fsum 保证求和确定），结果截断到 [0, 1]。
    """
    lo = max(0, n_drawn - (n_background - n_set))
    hi = min(n_set, n_drawn)
    if observed <= lo:
        return 1.0
    if observed > hi:
        return 0.0
    log_denominator = _log_choose(n_background, n_drawn)
    log_terms = [
        _log_choose(n_set, i)
        + _log_choose(n_background - n_set, n_drawn - i)
        - log_denominator
        for i in range(observed, hi + 1)
    ]
    p_value = math.fsum(math.exp(t) for t in log_terms)
    if p_value < 0.0:
        return 0.0
    if p_value > 1.0:
        return 1.0
    return p_value


def enrich_marker_gene_sets(
    gene_sets: GeneSets,
    markers: Dict[int, List[MarkerRecord]],
    *,
    kept_gene_ids: List[str],
    alpha: float,
    min_log_fc: float,
) -> EnrichmentResult:
    """逐簇 marker 对基因集的超几何富集检验。

    ``markers`` 为最终簇的 one-versus-rest Welch t 检验结果（与
    ``markers.tsv`` 同源）；背景为 ``kept_gene_ids``（保留基因）。
    集合与背景交集为空的基因集在评分阶段已报数据错误（退出码 4），
    这里不会遇到。
    """
    background = set(kept_gene_ids)
    n_background = len(kept_gene_ids)

    rows: List[EnrichmentRecord] = []
    n_marker_hits = 0
    for cluster in sorted(markers):
        hits = {
            r.gene_id
            for r in markers[cluster]
            if r.p_value_adj <= alpha and r.log_fc >= min_log_fc
        }
        n_markers = len(hits)
        n_marker_hits += n_markers

        cluster_rows: List[EnrichmentRecord] = []
        pvalues: List[float] = []
        for set_id in sorted(gene_sets.members):
            members = gene_sets.members[set_id]
            n_set_total = len(members)
            used = [gene_id for gene_id in members if gene_id in background]
            n_set_used = len(used)
            n_overlap = sum(1 for gene_id in used if gene_id in hits)
            if n_background > 0:
                expected_overlap = n_markers * n_set_used / n_background
            else:
                expected_overlap = 0.0
            # 除零取 0：期望为 0 时倍数富集定义为 0
            fold_enrichment = (
                n_overlap / expected_overlap if expected_overlap > 0.0 else 0.0
            )
            # 命中与否 × 属于集合与否的 2x2 表，四项均加 0.5
            a = n_overlap
            b = n_markers - n_overlap
            c = n_set_used - n_overlap
            d = n_background - n_markers - n_set_used + n_overlap
            odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
            p_value = hypergeometric_sf(
                n_overlap, n_background, n_set_used, n_markers
            )
            pvalues.append(p_value)
            cluster_rows.append(
                EnrichmentRecord(
                    cluster=cluster,
                    set_id=set_id,
                    n_set_total=n_set_total,
                    n_set_used=n_set_used,
                    n_markers=n_markers,
                    n_overlap=n_overlap,
                    expected_overlap=expected_overlap,
                    fold_enrichment=fold_enrichment,
                    odds_ratio=odds_ratio,
                    p_value=p_value,
                    p_value_adj=0.0,
                )
            )

        # 每簇内跨集合 BH 校正
        adjusted = benjamini_hochberg(pvalues)
        rows.extend(
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
        )

    rows.sort(
        key=lambda r: (r.cluster, r.p_value_adj, r.set_id, -r.odds_ratio)
    )
    n_significant_sets = sum(1 for r in rows if r.p_value_adj <= alpha)
    return EnrichmentResult(
        rows=rows,
        n_marker_hits=n_marker_hits,
        n_significant_sets=n_significant_sets,
    )
