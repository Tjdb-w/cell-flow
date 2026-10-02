"""结果文件的事务性发布：全部内容先在内存构建，再暂存、最后原子发布。

约束：
- 不覆盖任何已有文件；目标已存在且非空、目标不是目录、父目录不存在时抛 ``OutputPathError``；
- 所有结果先写入目标父目录下的隐藏暂存目录（与目标同一文件系统），逐文件
  临时文件 + 改名 + fsync，随后一次性发布：目标不存在则整目录 rename（原子），
  目标为空目录则逐条目 rename（run.json 最后），失败即撤回；
- ``run.json`` 最后发布，其存在即标志结果完整；任何中途失败都不会让它单独可见；
- 任一文件系统故障都把目标路径恢复到调用前状态（不存在/仍为空/不动非空目录），
  不残留临时文件，统一以 ``OutputPathError`` 上报，绝不把底层 OSError 透出；
- 数值用 ``repr`` 最短往返表示，文本与行序在同输入同版本下逐字节一致。
"""

import contextlib
import json
import math
import os
import shutil
import tempfile
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .errors import OutputPathError
from .kmeans import KMeansResult
from .markers import MarkerRecord, PairwiseMarkerRecord
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
class Artifacts:
    qc: QCResult
    data: NormalizedData
    pca: PCAResult
    clustering: KMeansResult
    markers: Dict[int, List[MarkerRecord]]
    pairwise_markers: List[Tuple[int, int, List[PairwiseMarkerRecord]]]
    run_info: Dict[str, Any]
    cluster_selection: Optional[ClusterSelectionResult] = None


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


def _validate_target(output_dir: str) -> str:
    """校验目标路径并返回其绝对路径；不创建任何东西。

    目标已存在但不是目录、已存在且非空、或父目录不存在时抛 ``OutputPathError``。
    目标不存在或存在且为空都允许（空目录稍后原地发布、复用）。
    """
    if output_dir is None or output_dir == "":
        raise OutputPathError("结果目录路径为空")
    abs_output = os.path.abspath(output_dir)
    if os.path.lexists(abs_output):
        if not os.path.isdir(abs_output):
            raise OutputPathError(f"输出路径已存在且不是目录：{output_dir}")
        try:
            contents = os.listdir(abs_output)
        except OSError as exc:
            raise OutputPathError(f"结果目录无法访问：{output_dir}（{exc}）")
        if contents:
            raise OutputPathError(
                f"结果目录已存在且非空：{output_dir}（含 {len(contents)} 个条目）"
            )
    else:
        parent = os.path.dirname(abs_output)
        if not os.path.isdir(parent):
            raise OutputPathError(f"输出目录的父目录不存在：{parent}")
    return abs_output


def _build_payload(artifacts: "Artifacts") -> List[Tuple[str, str]]:
    """按既有顺序构建全部结果文本（纯内存计算，不触碰文件系统）。"""
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
    ]
    # 仅 --n-clusters auto 产出候选评估表；显式整数模式文件集与基线一致
    if artifacts.cluster_selection is not None:
        payload.append(
            (
                "cluster_selection.tsv",
                _cluster_selection_tsv(artifacts.cluster_selection),
            )
        )
    # run.json 永远最后：它的存在即标志整组结果完整
    payload.append(("run.json", _run_json(artifacts.run_info)))
    return payload


def _stage_file(staging_dir: str, name: str, content: str) -> None:
    """把单个结果以临时文件 + 改名 + fsync 的方式落到暂存目录。"""
    final_path = os.path.join(staging_dir, name)
    fd, tmp_path = tempfile.mkstemp(
        prefix=f".{name}.", suffix=".tmp", dir=staging_dir
    )
    try:
        os.chmod(tmp_path, 0o644)
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, final_path)
    except OSError as exc:
        with contextlib.suppress(OSError):
            os.remove(tmp_path)
        raise OutputPathError(f"结果文件写出失败：{name}（{exc}）")


def _fsync_dir(path: str) -> None:
    """刷目录项，使其中的改名在崩溃后仍持久；不支持时静默跳过。"""
    try:
        fd = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(fd)
    except OSError:
        pass
    finally:
        os.close(fd)


def _remove_temp_tree(path: str) -> None:
    """尽力删除暂存目录，忽略清理自身的故障（原始失败优先上报）。"""
    shutil.rmtree(path, ignore_errors=True)


def _publish_new(abs_output: str, staging_dir: str) -> None:
    """目标原本不存在：整目录原子改名到目标，失败则删除暂存、目标不出现。"""
    try:
        os.rename(staging_dir, abs_output)
    except OSError as exc:
        _remove_temp_tree(staging_dir)
        raise OutputPathError(f"结果目录发布失败：{abs_output}（{exc}）")
    _fsync_dir(os.path.dirname(abs_output))


def _publish_into_empty(
    abs_output: str, staging_dir: str, names: List[str]
) -> None:
    """目标原本为空目录：按序逐条目改名进入，任一失败立即撤回已发布条目。

    run.json 在 names 末尾，故总是最后可见；撤回时它也最先被移除。
    """
    moved: List[str] = []
    current = ""
    try:
        for name in names:
            current = name
            os.rename(os.path.join(staging_dir, name), os.path.join(abs_output, name))
            moved.append(name)
        _fsync_dir(abs_output)
    except OSError as exc:
        # 逆序撤回：run.json（若已发布）最先消失，最终目标恢复为空
        for moved_name in reversed(moved):
            with contextlib.suppress(OSError):
                os.replace(
                    os.path.join(abs_output, moved_name),
                    os.path.join(staging_dir, moved_name),
                )
        _fsync_dir(abs_output)
        _remove_temp_tree(staging_dir)
        raise OutputPathError(f"结果文件发布失败：{current}（{exc}）")
    # 全部条目已就位，删除现已清空的暂存目录；
    # 此刻 run.json 已随最后一次改名持久落盘，提交已经完成，
    # 清理仅剩的空暂存目录不应让一次完整成功的运行转为失败
    try:
        os.rmdir(staging_dir)
    except OSError:
        _remove_temp_tree(staging_dir)


def publish_results(output_dir: str, artifacts: "Artifacts") -> List[str]:
    """事务性发布全部结果，返回按写出顺序排列的文件名（run.json 在最后）。

    先校验目标、再在目标父目录下的同文件系统暂存目录中写好全部文件，
    最后一次性发布。任一步骤失败都把目标恢复到调用前状态并抛
    ``OutputPathError``，绝不透出底层 OSError。
    """
    abs_output = _validate_target(output_dir)
    parent = os.path.dirname(abs_output)
    payload = _build_payload(artifacts)
    names = [name for name, _ in payload]

    # 暂存目录与目标同父目录 => 同一文件系统，后续 rename 为同卷原子改名
    try:
        staging_dir = tempfile.mkdtemp(prefix=".cell-flow-staging-", dir=parent)
    except OSError as exc:
        raise OutputPathError(f"暂存目录创建失败：{parent}（{exc}）")
    try:
        os.chmod(staging_dir, 0o755)
    except OSError as exc:
        _remove_temp_tree(staging_dir)
        raise OutputPathError(f"暂存目录创建失败：{parent}（{exc}）")

    try:
        for name, content in payload:
            _stage_file(staging_dir, name, content)
        _fsync_dir(staging_dir)

        if os.path.isdir(abs_output):
            # 空目录：校验阶段已确认其为空，发布前再确认一次以避免竞态写入
            try:
                if os.listdir(abs_output):
                    raise OutputPathError(f"结果目录已存在且非空：{output_dir}")
            except OSError as exc:
                raise OutputPathError(f"结果目录无法访问：{output_dir}（{exc}）")
            _publish_into_empty(abs_output, staging_dir, names)
        else:
            _publish_new(abs_output, staging_dir)
    except OutputPathError:
        # 发布/撤回例程已自行回收暂存目录；更早阶段的失败在此统一回收
        _remove_temp_tree(staging_dir)
        raise
    except OSError as exc:
        _remove_temp_tree(staging_dir)
        raise OutputPathError(f"结果写出失败：{output_dir}（{exc}）")

    return names
