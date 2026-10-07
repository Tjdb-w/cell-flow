"""结果文件写出：全部内容先在内存构建，再经同级隐藏暂存目录事务性发布。

约束：
- 预检先于任何创建：目标路径不是目录、结果目录已存在且非空、父目录不存在，
  一律抛 ``OutputPathError`` 且不创建或改动任何目录；
- 所有结果先写入父目录下的一个隐藏暂存目录，再以一次 ``os.replace`` 原子改名
  发布；``run.json`` 与其余文件同时出现，其存在即标志结果完整，不可能单独可见；
- 暂存或发布任一步失败都整体删除暂存目录并抛 ``OutputPathError``，目标路径
  恢复调用前状态（原不存在则不残留、原空目录仍为空、原非空目录原样保留）；
- 不覆盖任何已有文件；
- 数值用 ``repr`` 最短往返表示，文本与行序在同输入同版本下逐字节一致。
"""

import json
import math
import os
import shutil
import tempfile
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .errors import OutputPathError
from .abundance import DifferentialAbundance
from .adjusted_cluster_pseudobulk import AdjustedClusterPseudobulkDE
from .batch import BatchSummaryRow
from .cell_metadata import CellBatchReport
from .cell_types import CellTypeAnnotations
from .cluster_pseudobulk import ClusterPseudobulkDE
from .doublets import DoubletResult
from .enrichment import MarkerEnrichment
from .gene_sets import GeneSetScores
from .kmeans import KMeansResult
from .markers import GroupComparison, MarkerRecord, PairwiseMarkerRecord
from .normalize import NormalizedData
from .numeric_cluster_pseudobulk import NumericClusterPseudobulkDE
from .paired_cluster_pseudobulk import PairedClusterPseudobulkDE
from .pca import PCALoadings, PCAResult
from .pseudobulk_gene_sets import PseudobulkGeneSetScores
from .qc import QCResult
from .replicate import PseudobulkData
from .selection import ClusterSelectionResult
from .stability import (
    STABILITY_FILE_DISTRIBUTION,
    STABILITY_FILE_SCORES,
    STABILITY_FILE_SUMMARY,
    StabilityResult,
)

TOP_N_MARKERS = 20
TOP_N_LOADINGS = 20
HIST_BINS = 20


def fmt_float(value: float) -> str:
    """float64 的最短往返文本；非有限值用小写 inf/-inf/nan。"""
    return repr(float(value))


def fmt_bool(value: bool) -> str:
    return "true" if value else "false"


def tsv_row(values: Sequence[Any]) -> str:
    return "\t".join(str(v) for v in values)


@dataclass(frozen=True)
class Artifacts:
    qc: QCResult
    data: NormalizedData
    pca: PCAResult
    clustering: KMeansResult
    markers: Dict[int, List[MarkerRecord]]
    pairwise_markers: List[Tuple[int, int, List[PairwiseMarkerRecord]]]
    run_info: Dict[str, Any]
    pca_loadings: Optional[PCALoadings] = None
    cluster_selection: Optional[ClusterSelectionResult] = None
    group_markers: Optional[List[GroupComparison]] = None
    batch_summary: Optional[List[BatchSummaryRow]] = None
    gene_set_scores: Optional[GeneSetScores] = None
    cell_batch_report: Optional[CellBatchReport] = None
    pseudobulk: Optional[PseudobulkData] = None
    doublets: Optional[DoubletResult] = None
    stability: Optional[StabilityResult] = None
    cell_type_annotations: Optional[CellTypeAnnotations] = None
    marker_enrichment: Optional[MarkerEnrichment] = None
    pseudobulk_gene_set_scores: Optional[PseudobulkGeneSetScores] = None
    differential_abundance: Optional[DifferentialAbundance] = None
    cluster_pseudobulk: Optional[ClusterPseudobulkDE] = None
    paired_cluster_pseudobulk: Optional[PairedClusterPseudobulkDE] = None
    adjusted_cluster_pseudobulk: Optional[AdjustedClusterPseudobulkDE] = None
    numeric_cluster_pseudobulk: Optional[NumericClusterPseudobulkDE] = None


def _cells_tsv(qc: QCResult) -> str:
    lines = [
        tsv_row(
            [
                "cell_id",
                "total_counts",
                "detected_genes",
                "mitochondrial_counts",
                "mitochondrial_fraction",
                "retained",
            ]
        )
    ]
    for c in qc.cell_qc:
        lines.append(
            tsv_row(
                [
                    c.cell_id,
                    c.total_counts,
                    c.detected_genes,
                    c.mitochondrial_counts,
                    fmt_float(c.mitochondrial_fraction),
                    fmt_bool(c.retained),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _genes_tsv(qc: QCResult) -> str:
    lines = [
        tsv_row(["gene_id", "total_counts", "detected_cells", "retained"])
    ]
    for g in qc.gene_qc:
        lines.append(
            tsv_row(
                [
                    g.gene_id,
                    g.total_counts,
                    g.detected_cells,
                    fmt_bool(g.retained),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _pca_tsv(pca: PCAResult) -> str:
    header = ["cell_id"] + [f"PC{k + 1}" for k in range(pca.n_pcs)]
    lines = [tsv_row(header)]
    for c, cell_id in enumerate(pca.cell_ids):
        row = [cell_id] + [fmt_float(pca.scores[c][k]) for k in range(pca.n_pcs)]
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _clusters_tsv(pca: PCAResult, clustering: KMeansResult) -> str:
    lines = [tsv_row(["cell_id", "cluster"])]
    for c, cell_id in enumerate(pca.cell_ids):
        lines.append(tsv_row([cell_id, clustering.labels[c]]))
    return "\n".join(lines) + "\n"


def _markers_tsv(markers: Dict[int, List[MarkerRecord]]) -> str:
    lines = [
        tsv_row(
            [
                "cluster",
                "gene_id",
                "mean_in_cluster",
                "mean_out_cluster",
                "log_fc",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for cluster in sorted(markers):
        for r in markers[cluster]:
            lines.append(
                tsv_row(
                    [
                        cluster,
                        r.gene_id,
                        fmt_float(r.mean_in_cluster),
                        fmt_float(r.mean_out_cluster),
                        fmt_float(r.log_fc),
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _histogram(
    metric: str,
    retained: bool,
    values: List[float],
    lo: float,
    hi: float,
) -> List[List[Any]]:
    rows: List[List[Any]] = []
    width = (hi - lo) / HIST_BINS
    counts = [0] * HIST_BINS
    for value in values:
        if width <= 0.0:
            idx = 0
        else:
            idx = int((value - lo) / width)
            if idx < 0:
                idx = 0
            elif idx >= HIST_BINS:
                idx = HIST_BINS - 1
        counts[idx] += 1
    for b in range(HIST_BINS):
        start = lo + b * width
        end = lo + (b + 1) * width if width > 0.0 else hi
        rows.append(
            [
                metric,
                fmt_bool(retained),
                b,
                fmt_float(start),
                fmt_float(end),
                counts[b],
            ]
        )
    return rows


def _qc_chart_tsv(qc: QCResult) -> str:
    lines = [
        tsv_row(
            ["metric", "retained", "bin_index", "bin_start", "bin_end", "count"]
        )
    ]
    groups = ((True, [c for c in qc.cell_qc if c.retained]),
              (False, [c for c in qc.cell_qc if not c.retained]))
    for retained, cells in groups:
        detected = [float(c.detected_genes) for c in cells]
        totals = [float(c.total_counts) for c in cells]
        fractions = [c.mitochondrial_fraction for c in cells]

        max_detected = max(detected, default=0.0)
        max_total = max(totals, default=0.0)
        for rows in (
            _histogram("detected_genes", retained, detected, 0.0, max_detected),
            _histogram("total_counts", retained, totals, 0.0, max_total),
            _histogram("mitochondrial_fraction", retained, fractions, 0.0, 1.0),
        ):
            for row in rows:
                lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _pca_scatter_tsv(pca: PCAResult, clustering: KMeansResult) -> str:
    lines = [tsv_row(["cell_id", "PC1", "PC2", "cluster"])]
    for c, cell_id in enumerate(pca.cell_ids):
        pc1 = pca.scores[c][0] if pca.n_pcs >= 1 else 0.0
        pc2 = pca.scores[c][1] if pca.n_pcs >= 2 else 0.0
        lines.append(
            tsv_row(
                [
                    cell_id,
                    fmt_float(pc1),
                    fmt_float(pc2),
                    clustering.labels[c],
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _cluster_annotations_tsv(annotations: CellTypeAnnotations) -> str:
    """按 cluster 升序的逐簇注释结果。

    最高分不大于 0 时 annotation_status 为 unassigned，但 cell_type/score
    仍取最高分候选；n_markers_total 为该类型参考表标记总数，
    n_markers_used 为命中保留基因、实际参与平均的标记数。
    """
    lines = [
        tsv_row(
            [
                "cluster",
                "n_cells",
                "cell_type",
                "annotation_status",
                "score",
                "n_markers_total",
                "n_markers_used",
            ]
        )
    ]
    for a in annotations.clusters:
        lines.append(
            tsv_row(
                [
                    a.cluster,
                    a.n_cells,
                    a.cell_type,
                    a.annotation_status,
                    fmt_float(a.score),
                    a.n_markers_total,
                    a.n_markers_used,
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _cluster_annotation_chart_tsv(
    annotations: CellTypeAnnotations,
    pca: PCAResult,
    clustering: KMeansResult,
) -> str:
    """在 pca_scatter.tsv 的细胞顺序与 cell_id、PC1、PC2、cluster 四列上
    追加 cell_type、annotation_status、score；仅含最终细胞。"""
    by_cluster = {a.cluster: a for a in annotations.clusters}
    lines = [
        tsv_row(
            [
                "cell_id",
                "PC1",
                "PC2",
                "cluster",
                "cell_type",
                "annotation_status",
                "score",
            ]
        )
    ]
    for c, cell_id in enumerate(pca.cell_ids):
        pc1 = pca.scores[c][0] if pca.n_pcs >= 1 else 0.0
        pc2 = pca.scores[c][1] if pca.n_pcs >= 2 else 0.0
        a = by_cluster[clustering.labels[c]]
        lines.append(
            tsv_row(
                [
                    cell_id,
                    fmt_float(pc1),
                    fmt_float(pc2),
                    clustering.labels[c],
                    a.cell_type,
                    a.annotation_status,
                    fmt_float(a.score),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _pca_variance_tsv(pca: PCAResult) -> str:
    # 比例以本次实际输出的全部主成分解释方差之和为分母
    total = sum(pca.explained_variance)
    lines = [
        tsv_row(
            [
                "component",
                "explained_variance",
                "explained_variance_ratio",
                "cumulative_explained_variance_ratio",
            ]
        )
    ]
    cumulative = 0.0
    for k in range(pca.n_pcs):
        variance = pca.explained_variance[k]
        ratio = variance / total if total > 0.0 else 0.0
        # 分母为 0 时累计比例同样保持 0
        cumulative = cumulative + ratio if total > 0.0 else 0.0
        lines.append(
            tsv_row(
                [
                    f"PC{k + 1}",
                    fmt_float(variance),
                    fmt_float(ratio),
                    fmt_float(cumulative),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _pca_loadings_tsv(loadings: PCALoadings) -> str:
    """高变基因（既有选择顺序）× PC1..PCn 的带符号载荷宽表。"""
    header = ["gene_id"] + [f"PC{k + 1}" for k in range(loadings.n_pcs)]
    lines = [tsv_row(header)]
    for j, gene_id in enumerate(loadings.gene_ids):
        row = [gene_id] + [fmt_float(value) for value in loadings.loadings[j]]
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _pca_loading_chart_tsv(loadings: PCALoadings) -> str:
    """每个 PC 取绝对载荷最大的前 20 个高变基因。

    rank 自 1 起在每个 PC 内单独编号；绝对载荷并列时按 gene_id 升序。
    contribution 为该基因载荷平方除以该 PC 全部载荷平方和，分母为 0 时取 0。
    """
    lines = [tsv_row(["component", "rank", "gene_id", "loading", "contribution"])]
    n_genes = len(loadings.gene_ids)
    for k in range(loadings.n_pcs):
        column = [loadings.loadings[j][k] for j in range(n_genes)]
        total_sq = sum(value * value for value in column)
        order = sorted(
            range(n_genes),
            key=lambda j: (-abs(column[j]), loadings.gene_ids[j]),
        )
        for rank, j in enumerate(order[:TOP_N_LOADINGS], start=1):
            value = column[j]
            contribution = value * value / total_sq if total_sq > 0.0 else 0.0
            lines.append(
                tsv_row(
                    [
                        f"PC{k + 1}",
                        rank,
                        loadings.gene_ids[j],
                        fmt_float(value),
                        fmt_float(contribution),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _top_markers_tsv(markers: Dict[int, List[MarkerRecord]]) -> str:
    lines = [
        tsv_row(
            [
                "cluster",
                "rank",
                "gene_id",
                "mean_in_cluster",
                "mean_out_cluster",
                "log_fc",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for cluster in sorted(markers):
        for rank, r in enumerate(markers[cluster][:TOP_N_MARKERS], start=1):
            lines.append(
                tsv_row(
                    [
                        cluster,
                        rank,
                        r.gene_id,
                        fmt_float(r.mean_in_cluster),
                        fmt_float(r.mean_out_cluster),
                        fmt_float(r.log_fc),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _top_marker_expression_tsv(
    qc: QCResult,
    data: NormalizedData,
    pca: PCAResult,
    clustering: KMeansResult,
    markers: Dict[int, List[MarkerRecord]],
) -> str:
    lines = [
        tsv_row(["cluster", "gene_id", "cell_id", "in_cluster", "expression"])
    ]
    labels = clustering.labels
    gene_position = {gene_id: g for g, gene_id in enumerate(data.gene_ids)}
    # PCA 细胞顺序 = 保留细胞顺序
    for cluster in sorted(markers):
        top = markers[cluster][:TOP_N_MARKERS]
        for r in top:
            g = gene_position[r.gene_id]
            for c, cell_id in enumerate(pca.cell_ids):
                lines.append(
                    tsv_row(
                        [
                            cluster,
                            r.gene_id,
                            cell_id,
                            fmt_bool(labels[c] == cluster),
                            fmt_float(data.analysis_values[g][c]),
                        ]
                    )
                )
    return "\n".join(lines) + "\n"


def _neg_log10_p_adj(p_value_adj: float) -> str:
    if p_value_adj == 0.0:
        return "inf"
    if p_value_adj == 1.0:
        return "0"
    return fmt_float(-math.log10(p_value_adj))


def _pairwise_markers_tsv(
    comparisons: List[Tuple[int, int, List[PairwiseMarkerRecord]]],
) -> str:
    lines = [
        tsv_row(
            [
                "cluster_a",
                "cluster_b",
                "gene_id",
                "mean_in_a",
                "mean_in_b",
                "log_fc_a_vs_b",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for cluster_a, cluster_b, records in comparisons:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        cluster_a,
                        cluster_b,
                        r.gene_id,
                        fmt_float(r.mean_in_a),
                        fmt_float(r.mean_in_b),
                        fmt_float(r.log_fc_a_vs_b),
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _pairwise_marker_chart_tsv(
    comparisons: List[Tuple[int, int, List[PairwiseMarkerRecord]]],
) -> str:
    lines = [
        tsv_row(
            [
                "cluster_a",
                "cluster_b",
                "gene_id",
                "log_fc_a_vs_b",
                "p_value",
                "p_value_adj",
                "neg_log10_p_adj",
            ]
        )
    ]
    for cluster_a, cluster_b, records in comparisons:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        cluster_a,
                        cluster_b,
                        r.gene_id,
                        fmt_float(r.log_fc_a_vs_b),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                        _neg_log10_p_adj(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _group_markers_tsv(comparisons: List[GroupComparison]) -> str:
    lines = [
        tsv_row(
            [
                "comparison_type",
                "group_a",
                "group_b",
                "gene_id",
                "mean_in_a",
                "mean_in_b",
                "log_fc_a_vs_b",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for _, _, _, records in comparisons:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        r.gene_id,
                        fmt_float(r.mean_in_a),
                        fmt_float(r.mean_in_b),
                        fmt_float(r.log_fc_a_vs_b),
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _group_marker_chart_tsv(comparisons: List[GroupComparison]) -> str:
    lines = [
        tsv_row(
            [
                "comparison_type",
                "group_a",
                "group_b",
                "rank",
                "gene_id",
                "log_fc_a_vs_b",
                "p_value",
                "p_value_adj",
                "neg_log10_p_adj",
            ]
        )
    ]
    for _, _, _, records in comparisons:
        for rank, r in enumerate(records[:TOP_N_MARKERS], start=1):
            lines.append(
                tsv_row(
                    [
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        rank,
                        r.gene_id,
                        fmt_float(r.log_fc_a_vs_b),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                        _neg_log10_p_adj(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _pseudobulk_expression_tsv(bulk: PseudobulkData) -> str:
    """质控后保留基因（原行序）× 样本（输入细胞首次出现顺序）的
    pseudobulk log 归一化表达。

    数值为按样本对保留基因原始计数求和后，除以该样本文库总计数乘 10000、
    再取 ln(x + 1)；总计数为零的样本按 0.0 写出。
    """
    lines = [tsv_row(["gene_id"] + list(bulk.sample_ids))]
    for g, gene_id in enumerate(bulk.gene_ids):
        row = [gene_id] + [fmt_float(value) for value in bulk.values[g]]
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _pseudobulk_gene_set_scores_tsv(scores: PseudobulkGeneSetScores) -> str:
    """set_id 升序、集合内按 pseudobulk 样本顺序的逐样本评分。

    沿用 gene_set_scores.tsv 的评分布局，细胞列改为 sample_id。
    """
    lines = [
        tsv_row(
            ["set_id", "n_genes_total", "n_genes_used", "sample_id", "score"]
        )
    ]
    for r in scores.rows:
        lines.append(
            tsv_row(
                [
                    r.set_id,
                    r.n_total,
                    r.n_used,
                    r.sample_id,
                    fmt_float(r.score),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _pseudobulk_gene_set_score_chart_tsv(
    scores: PseudobulkGeneSetScores,
    bulk: PseudobulkData,
) -> str:
    """宽表：以 sample_id、group 开头，集合列按 set_id 升序；
    样本顺序沿用 pseudobulk 样本次序（输入细胞在矩阵列序中的首次出现）。"""
    score_by_sample: Dict[Tuple[str, str], float] = {
        (r.set_id, r.sample_id): r.score for r in scores.rows
    }
    set_order = scores.set_order
    lines = [tsv_row(["sample_id", "group"] + set_order)]
    for s, sample_id in enumerate(bulk.sample_ids):
        row: List[Any] = [sample_id, bulk.sample_groups[s]]
        row.extend(
            fmt_float(score_by_sample[(set_id, sample_id)])
            for set_id in set_order
        )
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _pseudobulk_gene_set_de_tsv(scores: PseudobulkGeneSetScores) -> str:
    """pseudobulk 基因集分组差异，沿用 group marker 布局：
    set_id 取代 gene_id、score_difference_a_vs_b 取代 log_fc_a_vs_b。

    比较先全部 one-vs-rest（group 升序）再全部 pairwise；每个比较内
    记录按校正 P 值升序、评分均值差降序、set_id 升序排列（统计模块已定序）。
    """
    lines = [
        tsv_row(
            [
                "comparison_type",
                "group_a",
                "group_b",
                "set_id",
                "mean_in_a",
                "mean_in_b",
                "score_difference_a_vs_b",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for _, _, _, records in scores.comparisons:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        r.gene_id,
                        fmt_float(r.mean_in_a),
                        fmt_float(r.mean_in_b),
                        fmt_float(r.log_fc_a_vs_b),
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _pseudobulk_gene_set_chart_tsv(scores: PseudobulkGeneSetScores) -> str:
    """每个比较取前 20 个集合：前三列同差异表，group_b 后加 rank，
    末列为 neg_log10_p_adj；rank 自 1 起在每个比较内单独编号。"""
    lines = [
        tsv_row(
            [
                "comparison_type",
                "group_a",
                "group_b",
                "rank",
                "set_id",
                "score_difference_a_vs_b",
                "p_value",
                "p_value_adj",
                "neg_log10_p_adj",
            ]
        )
    ]
    for _, _, _, records in scores.comparisons:
        for rank, r in enumerate(records[:TOP_N_MARKERS], start=1):
            lines.append(
                tsv_row(
                    [
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        rank,
                        r.gene_id,
                        fmt_float(r.log_fc_a_vs_b),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                        _neg_log10_p_adj(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _cluster_abundance_tsv(abundance: DifferentialAbundance) -> str:
    """样本 × 簇全组合的丰度表：按 sample_id、cluster 升序（统计模块已定序）。

    proportion 为该样本保留细胞中属于该簇的比例（n_cells / total_cells），
    样本中不属于该簇时 n_cells 计 0。
    """
    lines = [
        tsv_row(
            [
                "sample_id",
                "group",
                "cluster",
                "n_cells",
                "total_cells",
                "proportion",
            ]
        )
    ]
    for r in abundance.rows:
        lines.append(
            tsv_row(
                [
                    r.sample_id,
                    r.group,
                    r.cluster,
                    r.n_cells,
                    r.total_cells,
                    fmt_float(r.proportion),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _cluster_differential_abundance_tsv(abundance: DifferentialAbundance) -> str:
    """簇级样本差异丰度：沿用分组差异布局，cluster 取代 gene_id、
    proportion_difference_a_vs_b 取代 log_fc_a_vs_b。

    比较先全部 one-vs-rest（group 升序，group_b 为空）再全部 pairwise；
    每个比较内记录按校正 P 值升序、差异值降序、cluster 升序排列
    （统计模块已定序）。
    """
    lines = [
        tsv_row(
            [
                "comparison_type",
                "group_a",
                "group_b",
                "cluster",
                "mean_in_a",
                "mean_in_b",
                "proportion_difference_a_vs_b",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for _, _, _, records in abundance.comparisons:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        r.cluster,
                        fmt_float(r.mean_in_a),
                        fmt_float(r.mean_in_b),
                        fmt_float(r.difference),
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _cluster_abundance_chart_tsv(abundance: DifferentialAbundance) -> str:
    """每簇按比例降序、样本升序取前 20 个样本，rank 自 1 起逐簇编号；
    其余列沿用 cluster_abundance.tsv。"""
    lines = [
        tsv_row(
            [
                "cluster",
                "rank",
                "sample_id",
                "group",
                "n_cells",
                "total_cells",
                "proportion",
            ]
        )
    ]
    rows_by_cluster: Dict[int, List[Any]] = {}
    for r in abundance.rows:
        rows_by_cluster.setdefault(r.cluster, []).append(r)
    for cluster in sorted(rows_by_cluster):
        top = sorted(
            rows_by_cluster[cluster],
            key=lambda r: (-r.proportion, r.sample_id),
        )[:TOP_N_MARKERS]
        for rank, r in enumerate(top, start=1):
            lines.append(
                tsv_row(
                    [
                        cluster,
                        rank,
                        r.sample_id,
                        r.group,
                        r.n_cells,
                        r.total_cells,
                        fmt_float(r.proportion),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _cluster_pseudobulk_de_tsv(result: ClusterPseudobulkDE) -> str:
    """最终簇内 pseudobulk 分组差异表达：前缀为 cluster、comparison_type、
    group_a、group_b、gene_id，其余统计列沿用 group marker 口径。

    行序按 cluster 升序，每簇先全部 one-vs-rest（group 升序，group_b 为空）
    再全部 pairwise；每个比较内记录按校正 P 值升序、log_fc_a_vs_b 降序、
    gene_id 升序排列（统计模块已定序）。
    """
    lines = [
        tsv_row(
            [
                "cluster",
                "comparison_type",
                "group_a",
                "group_b",
                "gene_id",
                "mean_in_a",
                "mean_in_b",
                "log_fc_a_vs_b",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for _, _, _, _, records in result.comparisons:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        r.cluster,
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        r.gene_id,
                        fmt_float(r.mean_in_a),
                        fmt_float(r.mean_in_b),
                        fmt_float(r.log_fc_a_vs_b),
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _cluster_pseudobulk_de_chart_tsv(result: ClusterPseudobulkDE) -> str:
    """每个簇每项比较取前 20 个基因：沿用差异表列序，在 group_b 后加入
    自 1 起、每比较内单独编号的 rank，末列为 neg_log10_p_adj。"""
    lines = [
        tsv_row(
            [
                "cluster",
                "comparison_type",
                "group_a",
                "group_b",
                "rank",
                "gene_id",
                "log_fc_a_vs_b",
                "p_value",
                "p_value_adj",
                "neg_log10_p_adj",
            ]
        )
    ]
    for cluster, _, _, _, records in result.comparisons:
        for rank, r in enumerate(records[:TOP_N_MARKERS], start=1):
            lines.append(
                tsv_row(
                    [
                        cluster,
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        rank,
                        r.gene_id,
                        fmt_float(r.log_fc_a_vs_b),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                        _neg_log10_p_adj(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _paired_cluster_pseudobulk_de_tsv(
    result: PairedClusterPseudobulkDE,
) -> str:
    """最终簇内配对 pseudobulk 差异表达：列为 cluster、group_a、group_b、
    pair_count、gene_id、mean_difference、t_stat、p_value、p_value_adj。

    行序按 cluster 升序，每簇内按校正 P 值升序、mean_difference 降序、
    gene_id 升序排列（统计模块已定序）。
    """
    lines = [
        tsv_row(
            [
                "cluster",
                "group_a",
                "group_b",
                "pair_count",
                "gene_id",
                "mean_difference",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for cluster, pair_count, records in result.results:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        cluster,
                        r.group_a,
                        r.group_b,
                        pair_count,
                        r.gene_id,
                        fmt_float(r.mean_difference),
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _paired_cluster_pseudobulk_de_chart_tsv(
    result: PairedClusterPseudobulkDE,
) -> str:
    """每个簇取差异表前 20 行，在 pair_count 后加入自 1 起、每簇内单独
    编号的 rank；其余列与差异表一致。"""
    lines = [
        tsv_row(
            [
                "cluster",
                "group_a",
                "group_b",
                "pair_count",
                "rank",
                "gene_id",
                "mean_difference",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for cluster, pair_count, records in result.results:
        for rank, r in enumerate(records[:TOP_N_MARKERS], start=1):
            lines.append(
                tsv_row(
                    [
                        cluster,
                        r.group_a,
                        r.group_b,
                        pair_count,
                        rank,
                        r.gene_id,
                        fmt_float(r.mean_difference),
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _cluster_pseudobulk_adjusted_de_tsv(
    result: AdjustedClusterPseudobulkDE,
) -> str:
    """最终簇内协变量校正的 pseudobulk 差异表达：前缀为 cluster、
    comparison_type、group_a、group_b 与两组样本规模，效应与统计列为
    校正模型的组别系数及其标准误、自由度、t、P 与校正 P。

    行序按 cluster 升序，每簇先全部 one-vs-rest（group 升序，group_b 为空）
    再全部 pairwise，与 cluster_pseudobulk_de.tsv 的比较次序一致；每个比较内
    记录按校正 P 值升序、log_fc_adjusted 降序、gene_id 升序排列
    （统计模块已定序）。
    """
    lines = [
        tsv_row(
            [
                "cluster",
                "comparison_type",
                "group_a",
                "group_b",
                "n_samples_a",
                "n_samples_b",
                "gene_id",
                "mean_in_a",
                "mean_in_b",
                "log_fc_adjusted",
                "std_error",
                "df",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for _, _, _, _, records in result.comparisons:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        r.cluster,
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        r.n_samples_a,
                        r.n_samples_b,
                        r.gene_id,
                        fmt_float(r.mean_in_a),
                        fmt_float(r.mean_in_b),
                        fmt_float(r.log_fc_adjusted),
                        fmt_float(r.std_error),
                        r.df,
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _cluster_pseudobulk_adjusted_de_chart_tsv(
    result: AdjustedClusterPseudobulkDE,
) -> str:
    """每个簇每项比较取前 20 个基因：沿用差异表列序，在 group_b 后加入
    自 1 起、每比较内单独编号的 rank，末列为 neg_log10_p_adj。"""
    lines = [
        tsv_row(
            [
                "cluster",
                "comparison_type",
                "group_a",
                "group_b",
                "rank",
                "gene_id",
                "log_fc_adjusted",
                "p_value",
                "p_value_adj",
                "neg_log10_p_adj",
            ]
        )
    ]
    for cluster, _, _, _, records in result.comparisons:
        for rank, r in enumerate(records[:TOP_N_MARKERS], start=1):
            lines.append(
                tsv_row(
                    [
                        cluster,
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        rank,
                        r.gene_id,
                        fmt_float(r.log_fc_adjusted),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                        _neg_log10_p_adj(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _cluster_pseudobulk_numeric_covariate_de_tsv(
    result: NumericClusterPseudobulkDE,
) -> str:
    """最终簇内含数值型协变量的 pseudobulk 检验：前缀为 cluster、
    comparison_type、group_a、group_b、covariate、gene_id，统计列为对应
    模型项（组别项或数值协变量连续列）的系数 effect 及其标准误、自由度、
    t、P 与校正 P。

    行序按 cluster 升序，每簇先全部 one-vs-rest（group 升序，group_b 为空）
    再全部 pairwise，与 cluster_pseudobulk_de.tsv 的比较次序一致；每比较内
    先组别项（covariate 为 "group"）再各数值协变量（表头顺序），每个检验项
    内按校正 P 值升序、effect 降序、gene_id 升序排列（统计模块已定序）。
    """
    lines = [
        tsv_row(
            [
                "cluster",
                "comparison_type",
                "group_a",
                "group_b",
                "covariate",
                "gene_id",
                "effect",
                "std_error",
                "df",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for _, _, _, _, records in result.comparisons:
        for r in records:
            lines.append(
                tsv_row(
                    [
                        r.cluster,
                        r.comparison_type,
                        r.group_a,
                        r.group_b,
                        r.covariate,
                        r.gene_id,
                        fmt_float(r.effect),
                        fmt_float(r.std_error),
                        r.df,
                        fmt_float(r.t_stat),
                        fmt_float(r.p_value),
                        fmt_float(r.p_value_adj),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _cluster_pseudobulk_numeric_covariate_chart_tsv(
    result: NumericClusterPseudobulkDE,
) -> str:
    """每个（簇, 比较, 检验项）取差异表前 20 行：沿用差异表列，在 cluster
    之后加入自 1 起、每检验项内单独编号的 rank。

    检验项即组别项（covariate 为 "group"）与各数值协变量连续列；记录在统计
    模块中已按检验项分块、每块按校正 P 值升序、effect 降序、gene_id 升序
    定序，这里按 covariate 切换分块取各自前 20 条。
    """
    lines = [
        tsv_row(
            [
                "cluster",
                "rank",
                "comparison_type",
                "group_a",
                "group_b",
                "covariate",
                "gene_id",
                "effect",
                "std_error",
                "df",
                "t_stat",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for _, _, _, _, records in result.comparisons:
        block_start = 0
        while block_start < len(records):
            covariate = records[block_start].covariate
            block_end = block_start + 1
            while (
                block_end < len(records)
                and records[block_end].covariate == covariate
            ):
                block_end += 1
            for rank, r in enumerate(
                records[block_start:block_end][:TOP_N_MARKERS], start=1
            ):
                lines.append(
                    tsv_row(
                        [
                            r.cluster,
                            rank,
                            r.comparison_type,
                            r.group_a,
                            r.group_b,
                            r.covariate,
                            r.gene_id,
                            fmt_float(r.effect),
                            fmt_float(r.std_error),
                            r.df,
                            fmt_float(r.t_stat),
                            fmt_float(r.p_value),
                            fmt_float(r.p_value_adj),
                        ]
                    )
                )
            block_start = block_end
    return "\n".join(lines) + "\n"



def _normalized_expression_tsv(data: NormalizedData) -> str:
    """质控后保留基因（原行序）× 保留细胞（原列序）的 log 归一化表达。

    数值即管线中间矩阵 ``NormalizedData.values``：原始计数除以该细胞在
    全部输入基因上的总计数后乘 10000，再取 ln(x + 1)；总计数为零的保留
    细胞在归一化阶段已得到 0.0，这里按最短往返表示原样写出。
    """
    lines = [tsv_row(["gene_id"] + list(data.cell_ids))]
    for g, gene_id in enumerate(data.gene_ids):
        row = [gene_id] + [fmt_float(value) for value in data.values[g]]
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _batch_corrected_expression_tsv(data: NormalizedData) -> str:
    """批次均值中心化后的表达；行列顺序与 normalized_expression.tsv 一致。"""
    corrected = data.corrected_values
    if corrected is None:
        raise ValueError("无批次校正值可写出")
    lines = [tsv_row(["gene_id"] + list(data.cell_ids))]
    for g, gene_id in enumerate(data.gene_ids):
        row = [gene_id] + [fmt_float(value) for value in corrected[g]]
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _batch_summary_tsv(rows: List[BatchSummaryRow]) -> str:
    """按 batch 升序的批次汇总；后三列为批次内保留细胞均值。"""
    lines = [
        tsv_row(
            [
                "batch_id",
                "n_cells",
                "total_counts",
                "detected_genes",
                "mitochondrial_fraction",
            ]
        )
    ]
    for r in rows:
        lines.append(
            tsv_row(
                [
                    r.batch_id,
                    r.n_cells,
                    fmt_float(r.mean_total_counts),
                    fmt_float(r.mean_detected_genes),
                    fmt_float(r.mean_mitochondrial_fraction),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _gene_set_scores_tsv(scores: GeneSetScores) -> str:
    """set_id 升序、集合内按保留细胞原顺序的逐细胞评分。"""
    lines = [
        tsv_row(
            ["set_id", "n_genes_total", "n_genes_used", "cell_id", "score"]
        )
    ]
    for r in scores.rows:
        lines.append(
            tsv_row(
                [
                    r.set_id,
                    r.n_total,
                    r.n_used,
                    r.cell_id,
                    fmt_float(r.score),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _gene_set_score_chart_tsv(
    scores: GeneSetScores,
    pca: PCAResult,
    clustering: KMeansResult,
) -> str:
    """宽表：以 cell_id、cluster 开头，集合列按 set_id 升序；
    细胞顺序沿用 clusters.tsv（即保留细胞原顺序）。"""
    score_by_cell: Dict[Tuple[str, str], float] = {
        (r.set_id, r.cell_id): r.score for r in scores.rows
    }
    set_order = scores.set_order
    lines = [tsv_row(["cell_id", "cluster"] + set_order)]
    for c, cell_id in enumerate(pca.cell_ids):
        row: List[Any] = [cell_id, clustering.labels[c]]
        row.extend(
            fmt_float(score_by_cell[(set_id, cell_id)])
            for set_id in set_order
        )
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _marker_gene_set_enrichment_tsv(enrichment: MarkerEnrichment) -> str:
    """逐簇 × 集合的完整超几何富集结果。

    行序沿用统计模块：cluster、p_value_adj、set_id 升序，odds_ratio 降序。
    """
    lines = [
        tsv_row(
            [
                "cluster",
                "set_id",
                "n_set_total",
                "n_set_used",
                "n_markers",
                "n_overlap",
                "expected_overlap",
                "fold_enrichment",
                "odds_ratio",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    for r in enrichment.rows:
        lines.append(
            tsv_row(
                [
                    r.cluster,
                    r.set_id,
                    r.n_set_total,
                    r.n_set_used,
                    r.n_markers,
                    r.n_overlap,
                    fmt_float(r.expected_overlap),
                    fmt_float(r.fold_enrichment),
                    fmt_float(r.odds_ratio),
                    fmt_float(r.p_value),
                    fmt_float(r.p_value_adj),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _marker_gene_set_enrichment_chart_tsv(
    enrichment: MarkerEnrichment,
) -> str:
    """每簇取富集排序前 20 行，追加在每簇内自 1 起的 rank 列。"""
    lines = [
        tsv_row(
            [
                "cluster",
                "rank",
                "set_id",
                "n_set_total",
                "n_set_used",
                "n_markers",
                "n_overlap",
                "expected_overlap",
                "fold_enrichment",
                "odds_ratio",
                "p_value",
                "p_value_adj",
            ]
        )
    ]
    current_cluster = None
    rank = 0
    for r in enrichment.rows:
        if r.cluster != current_cluster:
            current_cluster = r.cluster
            rank = 1
        else:
            rank += 1
        if rank > TOP_N_MARKERS:
            continue
        lines.append(
            tsv_row(
                [
                    r.cluster,
                    rank,
                    r.set_id,
                    r.n_set_total,
                    r.n_set_used,
                    r.n_markers,
                    r.n_overlap,
                    fmt_float(r.expected_overlap),
                    fmt_float(r.fold_enrichment),
                    fmt_float(r.odds_ratio),
                    fmt_float(r.p_value),
                    fmt_float(r.p_value_adj),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _cell_metadata_tsv(report: CellBatchReport) -> str:
    """保留细胞的元数据透传：表头与字段保持输入原样，行序同 clusters.tsv。"""
    lines = [tsv_row(report.columns)]
    for row in report.rows:
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _batch_pca_scatter_tsv(report: CellBatchReport, pca: PCAResult) -> str:
    """按批次着色的降维坐标：实际用于聚类的 PCA 空间（校正后为校正空间）。"""
    lines = [tsv_row(["cell_id", "batch", "PC1", "PC2"])]
    for c, cell_id in enumerate(pca.cell_ids):
        pc1 = pca.scores[c][0] if pca.n_pcs >= 1 else 0.0
        pc2 = pca.scores[c][1] if pca.n_pcs >= 2 else 0.0
        lines.append(
            tsv_row(
                [
                    cell_id,
                    report.batches[c],
                    fmt_float(pc1),
                    fmt_float(pc2),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _batch_mixing_tsv(report: CellBatchReport) -> str:
    """校正前后各一个批次混合分数（15 最近邻中其他批次占比的均值）。"""
    lines = [tsv_row(["stage", "mixing_score"])]
    lines.append(tsv_row(["before", fmt_float(report.mixing_before)]))
    lines.append(tsv_row(["after", fmt_float(report.mixing_after)]))
    return "\n".join(lines) + "\n"


def _doublet_scores_tsv(result: DoubletResult) -> str:
    """QC 候选细胞（原列序）的逐细胞双细胞评分；rank 自 1 起。"""
    lines = [
        tsv_row(
            [
                "cell_id",
                "doublet_score",
                "doublet_rank",
                "doublet_flag",
                "retained_after_doublet_filter",
            ]
        )
    ]
    for i, cell_id in enumerate(result.cell_ids):
        lines.append(
            tsv_row(
                [
                    cell_id,
                    fmt_float(result.scores[i]),
                    result.ranks[i],
                    fmt_bool(result.flags[i]),
                    fmt_bool(not result.flags[i]),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _doublet_score_chart_tsv(result: DoubletResult) -> str:
    """doublet_score 的 20 个等宽箱；空箱保留，边界由分数最小/最大值确定。"""
    lines = [tsv_row(["bin_start", "bin_end", "cell_count"])]
    scores = result.scores
    lo = min(scores)
    hi = max(scores)
    width = (hi - lo) / HIST_BINS
    counts = [0] * HIST_BINS
    for value in scores:
        if width <= 0.0:
            idx = 0
        else:
            idx = int((value - lo) / width)
            if idx < 0:
                idx = 0
            elif idx >= HIST_BINS:
                idx = HIST_BINS - 1
        counts[idx] += 1
    for b in range(HIST_BINS):
        start = lo + b * width
        end = lo + (b + 1) * width if width > 0.0 else hi
        lines.append(
            tsv_row([fmt_float(start), fmt_float(end), counts[b]])
        )
    return "\n".join(lines) + "\n"


def _run_json(run_info: Dict[str, Any]) -> str:
    return json.dumps(run_info, indent=2, ensure_ascii=False) + "\n"


def _stability_scores_csv(result: StabilityResult) -> str:
    """逐次抽样的稳定性分数 CSV：一行一次抽样，按抽样次序排列。"""
    lines = ["sample,n_sampled_cells,n_intersection_cells,adjusted_rand_index"]
    for i in range(result.n_samples):
        lines.append(
            ",".join(
                [
                    str(i + 1),
                    str(result.sample_sizes[i]),
                    str(result.intersection_sizes[i]),
                    fmt_float(result.scores[i]),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _stability_summary_json(result: StabilityResult) -> str:
    """summary.json：抽样次数、抽样比例、随机种子与均值/中位数/最小/最大值。"""
    from .stability import summarize_scores

    stats = summarize_scores(result.scores)
    payload = {
        "n_samples": result.config.n_samples,
        "sample_fraction": result.config.sample_fraction,
        "seed": result.config.seed,
        "n_available_cells": result.n_available_cells,
        "sample_size": result.sample_sizes[0],
        "mean": stats["mean"],
        "median": stats["median"],
        "min": stats["min"],
        "max": stats["max"],
    }
    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"


def _stability_distribution_json(result: StabilityResult) -> str:
    """供绘图使用的分数分布数据：逐次分数序列与 20 个等宽箱。

    箱边界由分数最小/最大值确定（与既有等宽直方图口径一致），空箱保留。
    """
    scores = result.scores
    lo = min(scores)
    hi = max(scores)
    width = (hi - lo) / HIST_BINS
    counts = [0] * HIST_BINS
    for value in scores:
        if width <= 0.0:
            idx = 0
        else:
            idx = int((value - lo) / width)
            if idx < 0:
                idx = 0
            elif idx >= HIST_BINS:
                idx = HIST_BINS - 1
        counts[idx] += 1
    bins = []
    for b in range(HIST_BINS):
        start = lo + b * width
        end = lo + (b + 1) * width if width > 0.0 else hi
        bins.append(
            {
                "bin_index": b,
                "bin_start": float(start),
                "bin_end": float(end),
                "count": counts[b],
            }
        )
    payload = {
        "metric": "adjusted_rand_index",
        "n_samples": result.config.n_samples,
        "sample_fraction": result.config.sample_fraction,
        "seed": result.config.seed,
        "n_available_cells": result.n_available_cells,
        "sample_size": result.sample_sizes[0],
        "scores": [float(value) for value in scores],
        "bins": bins,
    }
    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"


def _cluster_selection_tsv(selection: ClusterSelectionResult) -> str:
    # 候选按 k 升序；布尔 true/false；浮点最短往返；无效行两项指标 nan
    lines = [
        tsv_row(
            [
                "k",
                "formed_clusters",
                "valid",
                "iterations",
                "within_cluster_sse",
                "mean_silhouette",
                "selected",
            ]
        )
    ]
    for row in selection.rows:
        lines.append(
            tsv_row(
                [
                    row.k,
                    row.formed_clusters,
                    fmt_bool(row.valid),
                    row.iterations,
                    fmt_float(row.within_cluster_sse),
                    fmt_float(row.mean_silhouette),
                    fmt_bool(row.selected),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _validate_target(output_dir: str) -> None:
    """预检目标路径，不创建、不改动任何目录。

    - 路径为空：拒绝；
    - 已存在且不是目录：拒绝；
    - 已存在且是目录：必须为空，复用该空目录（非空则拒绝）；
    - 不存在：其父目录必须已存在（沿用“父目录不存在”的既有判定）。
    """
    if output_dir is None or output_dir == "":
        raise OutputPathError("结果目录路径为空")
    if os.path.lexists(output_dir):
        if not os.path.isdir(output_dir) or os.path.islink(output_dir):
            # isdir 会跟随符号链接：指向目录的链接仍属“不是（普通）目录”，
            # 不能把结果发布到链接背后，故同样拒绝
            raise OutputPathError(f"输出路径已存在且不是目录：{output_dir}")
        contents = os.listdir(output_dir)
        if contents:
            raise OutputPathError(
                f"结果目录已存在且非空：{output_dir}（含 {len(contents)} 个条目）"
            )
    else:
        parent = os.path.dirname(os.path.abspath(output_dir))
        if not os.path.isdir(parent):
            raise OutputPathError(f"输出目录的父目录不存在：{parent}")


def _remove_stage(stage_dir: str) -> None:
    """删除暂存目录；删除本身失败不再向上传播（调用方即将报原始写出故障）。"""
    shutil.rmtree(stage_dir, ignore_errors=True)


def _write_staging_payload(
    stage_dir: str, payload: List[Tuple[str, str]]
) -> None:
    """把全部文件内容写入暂存目录；失败统一抛 :class:`OutputPathError`。"""
    for name, content in payload:
        stage_path = os.path.join(stage_dir, name)
        fd = None
        try:
            fd = os.open(stage_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
            with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
                fd = None
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
        except OSError as exc:
            if fd is not None:
                try:
                    os.close(fd)
                except OSError:
                    pass
            raise OutputPathError(
                f"结果文件写出失败：{name}（{exc}）"
            ) from exc


def _is_existing_result_dir(output_dir: str) -> bool:
    """目标是否为一个已发布的分析目录（普通目录且含 run.json 完整性标志）。"""
    if not os.path.lexists(output_dir):
        return False
    if not os.path.isdir(output_dir) or os.path.islink(output_dir):
        return False
    return os.path.isfile(os.path.join(output_dir, "run.json"))


def _publish_stability_overlay(
    output_dir: str,
    stability_payload: List[Tuple[str, str]],
    run_info: Dict[str, Any],
) -> List[str]:
    """向已有分析目录增量发布稳定性结果。

    仅覆盖三个稳定性文件与 run.json，目录中的其他产物原样保留、字节不动。
    每个文件先在目标目录内写隐藏临时文件并刷盘，再以一次 ``os.replace``
    原子改名覆盖同名文件；run.json 最后替换。任一文件失败都清理本函数产生
    的临时文件并抛 :class:`OutputPathError`，已存在的其他文件不受影响。
    """
    overlay = list(stability_payload) + [("run.json", _run_json(run_info))]
    written: List[Tuple[str, str]] = []
    temp_paths: List[str] = []
    try:
        for name, content in overlay:
            fd = None
            tmp_path = None
            try:
                fd, tmp_path = tempfile.mkstemp(prefix=f"..{name}.", dir=output_dir)
                with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
                    fd = None
                    handle.write(content)
                    handle.flush()
                    os.fsync(handle.fileno())
            except OSError as exc:
                if fd is not None:
                    try:
                        os.close(fd)
                    except OSError:
                        pass
                raise OutputPathError(
                    f"稳定性结果文件写出失败：{name}（{exc}）"
                ) from exc
            temp_paths.append(tmp_path)
            written.append((name, tmp_path))

        for name, tmp_path in written:
            try:
                os.replace(tmp_path, os.path.join(output_dir, name))
            except OSError as exc:
                raise OutputPathError(
                    f"稳定性结果文件发布失败：{name}（{exc}）"
                ) from exc
            # 已发布：该临时路径即目标路径，不得再清理
            temp_paths[temp_paths.index(tmp_path)] = ""
    finally:
        for tmp_path in temp_paths:
            if tmp_path:
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass
    return [name for name, _ in overlay]


def publish_results(output_dir: str, artifacts: Artifacts) -> List[str]:
    """事务性发布全部结果，返回文件名列表（run.json 在最后）。

    全新运行：全部文件先落到父目录下的一个同级隐藏暂存目录，随后一次
    ``os.replace`` 原子改名为目标目录。任何创建、写入、刷盘、改名或最终
    发布失败都清理暂存目录并抛 :class:`OutputPathError`，目标路径回到
    调用前状态。

    稳定性结果重跑（目标目录已存在且含 run.json，且本次启用稳定性分析）：
    不再重建主分析，而是增量覆盖三个稳定性文件与 run.json，目录中的其他
    分析产物原样保留。
    """
    stability_enabled = artifacts.stability is not None
    if stability_enabled and _is_existing_result_dir(output_dir):
        stability_payload = [
            (STABILITY_FILE_SCORES, _stability_scores_csv(artifacts.stability)),
            (
                STABILITY_FILE_SUMMARY,
                _stability_summary_json(artifacts.stability),
            ),
            (
                STABILITY_FILE_DISTRIBUTION,
                _stability_distribution_json(artifacts.stability),
            ),
        ]
        return _publish_stability_overlay(
            output_dir, stability_payload, artifacts.run_info
        )

    _validate_target(output_dir)

    parent = os.path.dirname(os.path.abspath(output_dir))
    final_names = [
        "cells.tsv",
        "genes.tsv",
        "pca.tsv",
        "clusters.tsv",
        "markers.tsv",
        "qc.tsv",
        "pca_scatter.tsv",
        "pca_variance.tsv",
        "top_markers.tsv",
        "top_marker_expression.tsv",
        "pairwise_markers.tsv",
        "pairwise_marker_chart.tsv",
        "normalized_expression.tsv",
    ]
    if artifacts.pca_loadings is not None:
        # 仅 --pca-loadings 运行产出载荷产物；未启用时文件集与基线一致
        final_names.append("pca_loadings.tsv")
        final_names.append("pca_loading_chart.tsv")
    if artifacts.cluster_selection is not None:
        # 仅 --n-clusters auto 产出候选评估表；显式整数模式文件集与基线一致
        final_names.append("cluster_selection.tsv")
    if artifacts.group_markers is not None:
        # 仅 --metadata 运行产出分组差异表达；无元数据运行文件集与基线一致
        final_names.append("group_markers.tsv")
        final_names.append("group_marker_chart.tsv")
    if artifacts.data.corrected_values is not None:
        # 仅 --batch-metadata 运行产出校正表达；无批次运行文件集与基线一致
        final_names.append("batch_corrected_expression.tsv")
    if artifacts.batch_summary is not None:
        final_names.append("batch_summary.tsv")
    if artifacts.gene_set_scores is not None:
        # 仅 --gene-sets 运行产出基因集评分；未提供时文件集与基线一致
        final_names.append("gene_set_scores.tsv")
        final_names.append("gene_set_score_chart.tsv")
    if artifacts.marker_enrichment is not None:
        # 仅 --enrich-markers（必带 --gene-sets）运行产出富集结果
        final_names.append("marker_gene_set_enrichment.tsv")
        final_names.append("marker_gene_set_enrichment_chart.tsv")
    if artifacts.cell_batch_report is not None:
        # 仅 --cell-metadata 运行产出跨样本批次校正产物；未提供时文件集不变
        final_names.append("cell_metadata.tsv")
        final_names.append("batch_pca_scatter.tsv")
        final_names.append("batch_mixing.tsv")
    if artifacts.pseudobulk is not None:
        # 仅 --replicate-metadata 运行产出 pseudobulk 产物；未提供时文件集不变
        final_names.append("pseudobulk_expression.tsv")
        final_names.append("pseudobulk_group_markers.tsv")
        final_names.append("pseudobulk_marker_chart.tsv")
    if artifacts.pseudobulk_gene_set_scores is not None:
        # 仅 --pseudobulk-gene-set-de（必带 --gene-sets 与
        # --replicate-metadata）运行产出 pseudobulk 基因集产物
        final_names.append("pseudobulk_gene_set_scores.tsv")
        final_names.append("pseudobulk_gene_set_score_chart.tsv")
        final_names.append("pseudobulk_gene_set_de.tsv")
        final_names.append("pseudobulk_gene_set_chart.tsv")
    if artifacts.differential_abundance is not None:
        # 仅 --differential-abundance（必带 --replicate-metadata）运行产出
        # 簇级样本差异丰度产物；未启用时文件集不变
        final_names.append("cluster_abundance.tsv")
        final_names.append("cluster_differential_abundance.tsv")
        final_names.append("cluster_abundance_chart.tsv")
    if artifacts.cluster_pseudobulk is not None:
        # 仅 --cluster-pseudobulk-de（必带 --replicate-metadata）运行产出
        # 簇内 pseudobulk 差异表达产物；未启用时文件集不变
        final_names.append("cluster_pseudobulk_de.tsv")
        final_names.append("cluster_pseudobulk_de_chart.tsv")
    if artifacts.paired_cluster_pseudobulk is not None:
        # 仅 --paired-cluster-pseudobulk-de（必带 --replicate-metadata、
        # --cluster-pseudobulk-de 与 --paired-replicate-metadata）运行产出
        # 簇内配对 pseudobulk 差异表达产物；未启用时文件集不变
        final_names.append("paired_cluster_pseudobulk_de.tsv")
        final_names.append("paired_cluster_pseudobulk_de_chart.tsv")
    if artifacts.adjusted_cluster_pseudobulk is not None:
        # 仅 --cluster-pseudobulk-adjusted-de（必带 --replicate-metadata、
        # --cluster-pseudobulk-de 与 --pseudobulk-covariates）运行产出
        # 协变量校正的簇内 pseudobulk 差异表达产物；未启用时文件集不变
        final_names.append("cluster_pseudobulk_adjusted_de.tsv")
        final_names.append("cluster_pseudobulk_adjusted_de_chart.tsv")
    if artifacts.numeric_cluster_pseudobulk is not None:
        # 仅提供 --pseudobulk-numeric-covariates（必带 --replicate-metadata、
        # --cluster-pseudobulk-de 与 --cluster-pseudobulk-adjusted-de）运行
        # 产出含数值协变量的簇内 pseudobulk 检验产物；未提供时文件集不变
        final_names.append("cluster_pseudobulk_numeric_covariate_de.tsv")
        final_names.append(
            "cluster_pseudobulk_numeric_covariate_chart.tsv"
        )
    if artifacts.doublets is not None:
        # 仅 --detect-doublets 运行产出双细胞评分产物；未启用时文件集不变
        final_names.append("doublet_scores.tsv")
        final_names.append("doublet_score_chart.tsv")
    if artifacts.cell_type_annotations is not None:
        # 仅 --cell-type-reference 运行产出注释产物；未提供时文件集不变
        final_names.append("cluster_annotations.tsv")
        final_names.append("cluster_annotation_chart.tsv")
    stability_payload: List[Tuple[str, str]] = []
    if artifacts.stability is not None:
        # 仅 --stability-analysis 运行产出三个稳定性文件；未启用时文件集不变
        stability_payload = [
            (
                STABILITY_FILE_SCORES,
                _stability_scores_csv(artifacts.stability),
            ),
            (
                STABILITY_FILE_SUMMARY,
                _stability_summary_json(artifacts.stability),
            ),
            (
                STABILITY_FILE_DISTRIBUTION,
                _stability_distribution_json(artifacts.stability),
            ),
        ]
        final_names.extend(name for name, _ in stability_payload)
    final_names.append("run.json")

    payload: List[Tuple[str, str]] = [
        ("cells.tsv", _cells_tsv(artifacts.qc)),
        ("genes.tsv", _genes_tsv(artifacts.qc)),
        ("pca.tsv", _pca_tsv(artifacts.pca)),
        ("clusters.tsv", _clusters_tsv(artifacts.pca, artifacts.clustering)),
        ("markers.tsv", _markers_tsv(artifacts.markers)),
        ("qc.tsv", _qc_chart_tsv(artifacts.qc)),
        ("pca_scatter.tsv", _pca_scatter_tsv(artifacts.pca, artifacts.clustering)),
        ("pca_variance.tsv", _pca_variance_tsv(artifacts.pca)),
        ("top_markers.tsv", _top_markers_tsv(artifacts.markers)),
        (
            "top_marker_expression.tsv",
            _top_marker_expression_tsv(
                artifacts.qc,
                artifacts.data,
                artifacts.pca,
                artifacts.clustering,
                artifacts.markers,
            ),
        ),
        (
            "pairwise_markers.tsv",
            _pairwise_markers_tsv(artifacts.pairwise_markers),
        ),
        (
            "pairwise_marker_chart.tsv",
            _pairwise_marker_chart_tsv(artifacts.pairwise_markers),
        ),
        (
            "normalized_expression.tsv",
            _normalized_expression_tsv(artifacts.data),
        ),
    ]
    if artifacts.pca_loadings is not None:
        payload.append(
            ("pca_loadings.tsv", _pca_loadings_tsv(artifacts.pca_loadings))
        )
        payload.append(
            (
                "pca_loading_chart.tsv",
                _pca_loading_chart_tsv(artifacts.pca_loadings),
            )
        )
    if artifacts.cluster_selection is not None:
        payload.append(
            (
                "cluster_selection.tsv",
                _cluster_selection_tsv(artifacts.cluster_selection),
            )
        )
    if artifacts.group_markers is not None:
        payload.append(
            ("group_markers.tsv", _group_markers_tsv(artifacts.group_markers))
        )
        payload.append(
            (
                "group_marker_chart.tsv",
                _group_marker_chart_tsv(artifacts.group_markers),
            )
        )
    if artifacts.data.corrected_values is not None:
        payload.append(
            (
                "batch_corrected_expression.tsv",
                _batch_corrected_expression_tsv(artifacts.data),
            )
        )
    if artifacts.batch_summary is not None:
        payload.append(
            ("batch_summary.tsv", _batch_summary_tsv(artifacts.batch_summary))
        )
    if artifacts.gene_set_scores is not None:
        payload.append(
            (
                "gene_set_scores.tsv",
                _gene_set_scores_tsv(artifacts.gene_set_scores),
            )
        )
        payload.append(
            (
                "gene_set_score_chart.tsv",
                _gene_set_score_chart_tsv(
                    artifacts.gene_set_scores,
                    artifacts.pca,
                    artifacts.clustering,
                ),
            )
        )
    if artifacts.marker_enrichment is not None:
        payload.append(
            (
                "marker_gene_set_enrichment.tsv",
                _marker_gene_set_enrichment_tsv(artifacts.marker_enrichment),
            )
        )
        payload.append(
            (
                "marker_gene_set_enrichment_chart.tsv",
                _marker_gene_set_enrichment_chart_tsv(
                    artifacts.marker_enrichment
                ),
            )
        )
    if artifacts.cell_batch_report is not None:
        payload.append(
            ("cell_metadata.tsv", _cell_metadata_tsv(artifacts.cell_batch_report))
        )
        payload.append(
            (
                "batch_pca_scatter.tsv",
                _batch_pca_scatter_tsv(artifacts.cell_batch_report, artifacts.pca),
            )
        )
        payload.append(
            ("batch_mixing.tsv", _batch_mixing_tsv(artifacts.cell_batch_report))
        )
    if artifacts.pseudobulk is not None:
        payload.append(
            (
                "pseudobulk_expression.tsv",
                _pseudobulk_expression_tsv(artifacts.pseudobulk),
            )
        )
        payload.append(
            (
                "pseudobulk_group_markers.tsv",
                _group_markers_tsv(artifacts.pseudobulk.comparisons),
            )
        )
        payload.append(
            (
                "pseudobulk_marker_chart.tsv",
                _group_marker_chart_tsv(artifacts.pseudobulk.comparisons),
            )
        )
    if artifacts.pseudobulk_gene_set_scores is not None:
        payload.append(
            (
                "pseudobulk_gene_set_scores.tsv",
                _pseudobulk_gene_set_scores_tsv(artifacts.pseudobulk_gene_set_scores),
            )
        )
        payload.append(
            (
                "pseudobulk_gene_set_score_chart.tsv",
                _pseudobulk_gene_set_score_chart_tsv(
                    artifacts.pseudobulk_gene_set_scores,
                    artifacts.pseudobulk,
                ),
            )
        )
        payload.append(
            (
                "pseudobulk_gene_set_de.tsv",
                _pseudobulk_gene_set_de_tsv(artifacts.pseudobulk_gene_set_scores),
            )
        )
        payload.append(
            (
                "pseudobulk_gene_set_chart.tsv",
                _pseudobulk_gene_set_chart_tsv(
                    artifacts.pseudobulk_gene_set_scores
                ),
            )
        )
    if artifacts.differential_abundance is not None:
        payload.append(
            (
                "cluster_abundance.tsv",
                _cluster_abundance_tsv(artifacts.differential_abundance),
            )
        )
        payload.append(
            (
                "cluster_differential_abundance.tsv",
                _cluster_differential_abundance_tsv(
                    artifacts.differential_abundance
                ),
            )
        )
        payload.append(
            (
                "cluster_abundance_chart.tsv",
                _cluster_abundance_chart_tsv(artifacts.differential_abundance),
            )
        )
    if artifacts.cluster_pseudobulk is not None:
        payload.append(
            (
                "cluster_pseudobulk_de.tsv",
                _cluster_pseudobulk_de_tsv(artifacts.cluster_pseudobulk),
            )
        )
        payload.append(
            (
                "cluster_pseudobulk_de_chart.tsv",
                _cluster_pseudobulk_de_chart_tsv(artifacts.cluster_pseudobulk),
            )
        )
    if artifacts.paired_cluster_pseudobulk is not None:
        payload.append(
            (
                "paired_cluster_pseudobulk_de.tsv",
                _paired_cluster_pseudobulk_de_tsv(
                    artifacts.paired_cluster_pseudobulk
                ),
            )
        )
        payload.append(
            (
                "paired_cluster_pseudobulk_de_chart.tsv",
                _paired_cluster_pseudobulk_de_chart_tsv(
                    artifacts.paired_cluster_pseudobulk
                ),
            )
        )
    if artifacts.adjusted_cluster_pseudobulk is not None:
        payload.append(
            (
                "cluster_pseudobulk_adjusted_de.tsv",
                _cluster_pseudobulk_adjusted_de_tsv(
                    artifacts.adjusted_cluster_pseudobulk
                ),
            )
        )
        payload.append(
            (
                "cluster_pseudobulk_adjusted_de_chart.tsv",
                _cluster_pseudobulk_adjusted_de_chart_tsv(
                    artifacts.adjusted_cluster_pseudobulk
                ),
            )
        )
    if artifacts.numeric_cluster_pseudobulk is not None:
        payload.append(
            (
                "cluster_pseudobulk_numeric_covariate_de.tsv",
                _cluster_pseudobulk_numeric_covariate_de_tsv(
                    artifacts.numeric_cluster_pseudobulk
                ),
            )
        )
        payload.append(
            (
                "cluster_pseudobulk_numeric_covariate_chart.tsv",
                _cluster_pseudobulk_numeric_covariate_chart_tsv(
                    artifacts.numeric_cluster_pseudobulk
                ),
            )
        )
    if artifacts.doublets is not None:
        payload.append(
            ("doublet_scores.tsv", _doublet_scores_tsv(artifacts.doublets))
        )
        payload.append(
            (
                "doublet_score_chart.tsv",
                _doublet_score_chart_tsv(artifacts.doublets),
            )
        )
    if artifacts.cell_type_annotations is not None:
        payload.append(
            (
                "cluster_annotations.tsv",
                _cluster_annotations_tsv(artifacts.cell_type_annotations),
            )
        )
        payload.append(
            (
                "cluster_annotation_chart.tsv",
                _cluster_annotation_chart_tsv(
                    artifacts.cell_type_annotations,
                    artifacts.pca,
                    artifacts.clustering,
                ),
            )
        )
    payload.extend(stability_payload)
    payload.append(("run.json", _run_json(artifacts.run_info)))

    stage_dir: Optional[str] = None
    try:
        try:
            stage_dir = tempfile.mkdtemp(
                prefix="..cell-flow-staging.", dir=parent
            )
            # mkdtemp 默认 0700；发布后的结果目录沿用 makedirs 的默认权限（受 umask 约束）
            umask = os.umask(0)
            os.umask(umask)
            os.chmod(stage_dir, 0o777 & ~umask)
        except OSError as exc:
            raise OutputPathError(f"结果暂存目录创建失败：{exc}") from exc

        _write_staging_payload(stage_dir, payload)

        try:
            os.replace(stage_dir, os.path.abspath(output_dir))
        except OSError as exc:
            raise OutputPathError(f"结果目录发布失败：{output_dir}（{exc}）") from exc
        stage_dir = None  # 已发布：暂存目录即目标目录，不得再清理
    finally:
        if stage_dir is not None:
            # 发布前失败：仅删除我们自己的暂存目录，绝不触碰目标目录。
            _remove_stage(stage_dir)

    return final_names
