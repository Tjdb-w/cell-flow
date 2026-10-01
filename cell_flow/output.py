"""结果文件写出：全部内容先在内存构建，再以临时文件 + 原子替换落盘。

约束：
- 不覆盖任何已有文件；目标目录已存在且非空时由上游抛出 ``OutputPathError``；
- ``run.json`` 最后写出，其存在即标志结果完整；
- 数值用 ``repr`` 最短往返表示，文本与行序在同输入同版本下逐字节一致。
"""

import json
import os
import tempfile
from dataclasses import dataclass
from typing import Any, Dict, List, Sequence

from .errors import OutputPathError
from .kmeans import KMeansResult
from .markers import MarkerRecord
from .normalize import NormalizedData
from .pca import PCAResult
from .qc import QCResult

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
class Artifacts:
    qc: QCResult
    data: NormalizedData
    pca: PCAResult
    clustering: KMeansResult
    markers: Dict[int, List[MarkerRecord]]
    run_info: Dict[str, Any]


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
                            fmt_float(data.values[g][c]),
                        ]
                    )
                )
    return "\n".join(lines) + "\n"


def _run_json(run_info: Dict[str, Any]) -> str:
    return json.dumps(run_info, indent=2, ensure_ascii=False) + "\n"


def ensure_output_dir(output_dir: str) -> None:
    """目录不存在则创建；存在且非空则抛错；存在且为空则复用。"""
    if output_dir is None or output_dir == "":
        raise OutputPathError("结果目录路径为空")
    if os.path.exists(output_dir):
        if not os.path.isdir(output_dir):
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
        os.makedirs(output_dir)


def write_results(output_dir: str, artifacts: Artifacts) -> List[str]:
    """原子写出全部结果，返回写出的文件名（run.json 在最后）。"""
    payload: List[tuple[str, str]] = [
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
        ("run.json", _run_json(artifacts.run_info)),
    ]

    temp_paths: List[str] = []
    final_names: List[str] = []
    try:
        for name, content in payload:
            final_path = os.path.join(output_dir, name)
            if os.path.exists(final_path):
                # 绝不覆盖已有文件
                raise OutputPathError(f"结果文件已存在，拒绝覆盖：{final_path}")
            fd, tmp_path = tempfile.mkstemp(
                prefix=f".{name}.", suffix=".tmp", dir=output_dir
            )
            temp_paths.append(tmp_path)
            os.chmod(tmp_path, 0o644)
            try:
                with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
                    handle.write(content)
                    handle.flush()
                    os.fsync(handle.fileno())
            except OSError as exc:
                raise OutputPathError(f"结果文件写出失败：{name}（{exc}）")
        for (name, _), tmp_path in list(zip(payload, temp_paths)):
            os.replace(tmp_path, os.path.join(output_dir, name))
            final_names.append(name)
        temp_paths = []
    finally:
        for tmp_path in temp_paths:
            try:
                os.remove(tmp_path)
            except OSError:
                pass
    return final_names
