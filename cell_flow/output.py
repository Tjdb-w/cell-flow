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
from .kmeans import KMeansResult
from .markers import GroupComparison, MarkerRecord, PairwiseMarkerRecord
from .normalize import NormalizedData
from .pca import PCAResult
from .qc import QCResult
from .selection import ClusterSelectionResult

TOP_N_MARKERS = 20
HIST_BINS = 20


def fmt_float(value: float) -> str:
    """float64 的最短往返文本；非有限值用小写 inf/-inf/nan。"""
    return repr(float(value))


def fmt_bool(value: bool) -> str:
    return "true" if value else "false"


def tsv_row(values: Sequence[Any]) -> str:
    return "\t".join(str(v) for v in values)


@dataclass(frozen=True)
class BatchSummaryRow:
    """批次汇总中的一行：批次内保留细胞数与其三项指标的批次均值。"""

    batch_id: str
    n_cells: int
    total_counts: float              # 批次内保留细胞平均总计数
    detected_genes: float            # 批次内保留细胞平均检出基因数
    mitochondrial_fraction: float    # 批次内保留细胞平均线粒体比例


@dataclass(frozen=True)
class Artifacts:
    qc: QCResult
    data: NormalizedData
    pca: PCAResult
    clustering: KMeansResult
    markers: Dict[int, List[MarkerRecord]]
    pairwise_markers: List[Tuple[int, int, List[PairwiseMarkerRecord]]]
    run_info: Dict[str, Any]
    cluster_selection: Optional[ClusterSelectionResult] = None
    group_markers: Optional[List[GroupComparison]] = None
    batch_summary: Optional[List[BatchSummaryRow]] = None


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
    expression = data.analysis_values
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
                            fmt_float(expression[g][c]),
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


def _normalized_expression_tsv(data: NormalizedData) -> str:
    """质控后保留基因（原行序）× 保留细胞（原列序）的 log 归一化表达。

    数值即管线中间矩阵 ``NormalizedData.values``：原始计数除以该细胞在
    全部输入基因上的总计数后乘 10000，再取 ln(x + 1)；总计数为零的保留
    细胞在归一化阶段已得到 0.0，这里按最短往返表示原样写出。
    即使提供批次元数据，本文件仍写未校正原值；校正值见
    ``batch_corrected_expression.tsv``。
    """
    lines = [tsv_row(["gene_id"] + list(data.cell_ids))]
    for g, gene_id in enumerate(data.gene_ids):
        row = [gene_id] + [fmt_float(value) for value in data.values[g]]
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _batch_corrected_expression_tsv(data: NormalizedData) -> str:
    """与 ``normalized_expression.tsv`` 同形的批次均值中心化表达。

    每个值为对应 log 归一化值减去所属批次的该基因保留细胞均值、再加回
    全部保留细胞的该基因总均值；行序（保留基因原行序）与列序（保留细胞
    原列序）均不改变。
    """
    corrected = data.batch_corrected_values
    lines = [tsv_row(["gene_id"] + list(data.cell_ids))]
    for g, gene_id in enumerate(data.gene_ids):
        row = [gene_id] + [fmt_float(value) for value in corrected[g]]
        lines.append(tsv_row(row))
    return "\n".join(lines) + "\n"


def _batch_summary_tsv(rows: List[BatchSummaryRow]) -> str:
    """按 batch_id 升序的批次汇总；后三列为批次内保留细胞均值。"""
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
    for row in rows:
        lines.append(
            tsv_row(
                [
                    row.batch_id,
                    row.n_cells,
                    fmt_float(row.total_counts),
                    fmt_float(row.detected_genes),
                    fmt_float(row.mitochondrial_fraction),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _run_json(run_info: Dict[str, Any]) -> str:
    return json.dumps(run_info, indent=2, ensure_ascii=False) + "\n"


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


def publish_results(output_dir: str, artifacts: Artifacts) -> List[str]:
    """事务性发布全部结果，返回文件名列表（run.json 在最后）。

    全部文件先落到父目录下的一个同级隐藏暂存目录，随后一次 ``os.replace``
    原子改名为目标目录。任何创建、写入、刷盘、改名或最终发布失败都清理暂存
    目录并抛 :class:`OutputPathError`，目标路径回到调用前状态。
    """
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
    if artifacts.batch_summary is not None:
        # 仅 --batch-metadata 运行产出批次校正值与批次汇总；
        # 无批次元数据运行文件集与基线一致
        final_names.append("batch_corrected_expression.tsv")
        final_names.append("batch_summary.tsv")
    if artifacts.cluster_selection is not None:
        # 仅 --n-clusters auto 产出候选评估表；显式整数模式文件集与基线一致
        final_names.append("cluster_selection.tsv")
    if artifacts.group_markers is not None:
        # 仅 --metadata 运行产出分组差异表达；无元数据运行文件集与基线一致
        final_names.append("group_markers.tsv")
        final_names.append("group_marker_chart.tsv")
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
    if artifacts.batch_summary is not None:
        payload.append(
            (
                "batch_corrected_expression.tsv",
                _batch_corrected_expression_tsv(artifacts.data),
            )
        )
        payload.append(
            ("batch_summary.tsv", _batch_summary_tsv(artifacts.batch_summary))
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
