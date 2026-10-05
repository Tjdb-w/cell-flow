"""自动细胞类型注释：读入与计算。

``--cell-type-reference`` 指向 UTF-8 制表符文本，表头恰为 ``cell_type``、
``gene_id`` 两列，每行一个标记关系；字段非空、``(cell_type, gene_id)``
组合唯一，至少一条数据行。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，
与文件名无关）。基因 ID 与质控后保留基因精确匹配，不做大小写、别名、前缀
转换；矩阵外或被剔除（质控剔除/非高变）标记只计总数，不参与评分。任何
不合法（空路径、文件不存在/不可读、非普通文件、UTF-8、表头、空值、重复、
数据行不足，或 gzip 多成员、尾随数据、截断、CRC/长度错误）都抛
:class:`cell_flow.errors.CellFlowInputError`（退出码 2），且在触碰输出
目录之前失败。

评分只用质控后保留的细胞与基因：有批次校正取实际用于聚类的批次均值中心化
值，否则取 log 归一化值。对每个（簇, 类型），仅使用同时列于该类型标记表
且属于保留基因的标记基因，逐基因求该簇均值减其余簇均值，再对这些标记
取算术平均。无匹配保留基因的类型在该簇不参与竞争；允许同一类型注释多个簇。
每簇取分数最高的类型，并列按 cell_type 的 Unicode 码点升序取第一；最高
分不大于 0 时状态为 ``unassigned``，但仍保留候选类型与分数。全部标记与
保留基因无交集时报数据错误（退出码 4）。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowDataError, CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single

STATUS_ASSIGNED = "assigned"
STATUS_UNASSIGNED = "unassigned"


@dataclass(frozen=True)
class CellTypeReference:
    """细胞类型标记参考，按文件行序保留类型出现顺序与每个类型的标记。"""

    path: str
    name: str                    # 文件名（路径末段）
    sha256: str                  # 原始字节（gzip 即压缩字节）的 SHA-256
    markers: Dict[str, List[str]]  # cell_type -> 标记基因 ID（文件行序）


@dataclass(frozen=True)
class TypeScore:
    """一个类型在一个簇上的候选分数。"""

    cell_type: str
    n_total: int                 # 该类型参考表中的标记总数
    n_used: int                  # 命中保留基因、实际参与平均的标记数
    score: float


@dataclass(frozen=True)
class ClusterAnnotation:
    """一个最终簇的注释结果（含全部候选）。"""

    cluster: int
    n_cells: int
    cell_type: str               # 获胜类型；unassigned 时仍给出最高候选
    annotation_status: str       # assigned / unassigned
    score: float                 # 获胜候选分数；无候选时为 0.0
    n_markers_total: int         # 获胜类型标记总数；无候选时为 0
    n_markers_used: int          # 获胜类型使用标记数；无候选时为 0
    candidates: List[TypeScore]  # 全部参与竞争类型，按 cell_type 升序


@dataclass(frozen=True)
class CellTypeAnnotations:
    """全部最终簇的注释结果，按 cluster 升序。"""

    clusters: List[ClusterAnnotation]
    type_order: List[str]        # 参考表中全部类型，按 cell_type 升序
    totals: Dict[str, int]       # cell_type -> 参考表标记总数


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_cell_type_reference(path: str) -> CellTypeReference:
    """读取并严格校验细胞类型参考文件（不依赖表达矩阵，可在触碰输出目录前失败）。"""
    if path is None or path == "":
        _fail("细胞类型参考文件路径为空")
    if not os.path.exists(path):
        _fail(f"细胞类型参考文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"细胞类型参考路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"细胞类型参考文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"细胞类型参考文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"细胞类型参考文件不是合法的 UTF-8 文本：{path}")

    # splitlines 兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("细胞类型参考表头缺失：文件为空")
    header = lines[0].split("\t")
    if header != ["cell_type", "gene_id"]:
        _fail("细胞类型参考表头必须恰为 cell_type 和 gene_id 两列")

    markers: Dict[str, List[str]] = {}
    seen_pairs = set()
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"细胞类型参考第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != 2:
            _fail(
                f"细胞类型参考第 {offset} 行列数为 {len(fields)}，"
                f"与表头两列不一致"
            )
        cell_type, gene_id = fields
        if cell_type == "":
            _fail(f"细胞类型参考第 {offset} 行 cell_type 为空")
        if gene_id == "":
            _fail(f"细胞类型参考第 {offset} 行 gene_id 为空")
        pair = (cell_type, gene_id)
        if pair in seen_pairs:
            _fail(
                "细胞类型标记关系重复："
                f"cell_type={cell_type!r}、gene_id={gene_id!r}"
            )
        seen_pairs.add(pair)
        markers.setdefault(cell_type, []).append(gene_id)

    if not markers:
        _fail("细胞类型参考没有任何数据行")

    return CellTypeReference(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        markers=markers,
    )


def annotate_clusters(
    reference: CellTypeReference,
    *,
    kept_gene_ids: List[str],
    labels: List[int],
    analysis_values: List[List[float]],
) -> CellTypeAnnotations:
    """计算每个最终簇对各细胞类型的标记平均分并选出获胜类型。

    ``analysis_values[g][c]`` 与 ``kept_gene_ids``、``labels[c]`` 对齐，
    为批次中心化值（有批次校正）或 log 归一化值。仅使用命中保留基因的
    标记，逐基因求簇内均值减其余簇均值后取算术平均；矩阵外或被剔除的
    标记只计总数。全部类型与保留基因都无交集时无法注释，抛数据错误
    （退出码 4）。
    """
    gene_position: Dict[str, int] = {
        gene_id: g for g, gene_id in enumerate(kept_gene_ids)
    }

    n_cells = len(labels)
    cells_by_cluster: Dict[int, List[int]] = {}
    for c, label in enumerate(labels):
        cells_by_cluster.setdefault(label, []).append(c)

    # 每个类型的标记总数与命中保留基因的行下标（参考表已保证
    # (cell_type, gene_id) 唯一，直接建下标映射）
    totals: Dict[str, int] = {}
    used_positions: Dict[str, List[int]] = {}
    for cell_type in sorted(reference.markers):
        entries = reference.markers[cell_type]
        positions = [
            gene_position[gene_id]
            for gene_id in entries
            if gene_id in gene_position
        ]
        totals[cell_type] = len(entries)
        used_positions[cell_type] = positions

    if not any(used_positions.values()):
        raise CellFlowDataError(
            "细胞类型参考标记与质控后保留基因的交集为空，无法注释"
        )

    # 每个 (簇, 基因) 的簇内均值与其余簇均值只计算一次，
    # 供共享同一标记基因的各类型复用
    present = sorted(cells_by_cluster)
    mean_in: Dict[int, List[float]] = {}
    mean_out: Dict[int, List[float]] = {}
    for cluster in present:
        in_cells = cells_by_cluster[cluster]
        out_cells = [c for c in range(n_cells) if labels[c] != cluster]
        n_in = len(in_cells)
        n_out = len(out_cells)
        in_means: List[float] = [0.0] * len(kept_gene_ids)
        out_means: List[float] = [0.0] * len(kept_gene_ids)
        for g in range(len(kept_gene_ids)):
            row = analysis_values[g]
            in_means[g] = sum(row[c] for c in in_cells) / n_in
            out_means[g] = sum(row[c] for c in out_cells) / n_out
        mean_in[cluster] = in_means
        mean_out[cluster] = out_means

    annotations: List[ClusterAnnotation] = []
    for cluster in present:
        in_means = mean_in[cluster]
        out_means = mean_out[cluster]
        candidates: List[TypeScore] = []
        for cell_type in sorted(reference.markers):
            positions = used_positions[cell_type]
            if not positions:
                # 无可用标记的类型在该簇不参与竞争
                continue
            # 沿用管线既有口径（内置 sum）累加，保证与归一化/批次校正同精度
            total = sum(in_means[g] - out_means[g] for g in positions)
            n_used = len(positions)
            candidates.append(
                TypeScore(
                    cell_type=cell_type,
                    n_total=totals[cell_type],
                    n_used=n_used,
                    score=total / n_used,
                )
            )

        if candidates:
            # 分数降序、cell_type 的 Unicode 码点升序取第一；
            # Python 字符串默认按 Unicode 码点比较
            winner = min(
                candidates, key=lambda t: (-t.score, t.cell_type)
            )
            cell_type = winner.cell_type
            score = winner.score
            n_total = winner.n_total
            n_used = winner.n_used
            status = STATUS_ASSIGNED if score > 0.0 else STATUS_UNASSIGNED
        else:
            # 理论上不会发生：已在上方拒绝全部类型无交集的情形
            cell_type = ""
            score = 0.0
            n_total = 0
            n_used = 0
            status = STATUS_UNASSIGNED

        annotations.append(
            ClusterAnnotation(
                cluster=cluster,
                n_cells=len(cells_by_cluster[cluster]),
                cell_type=cell_type,
                annotation_status=status,
                score=score,
                n_markers_total=n_total,
                n_markers_used=n_used,
                candidates=candidates,
            )
        )

    return CellTypeAnnotations(
        clusters=annotations,
        type_order=sorted(reference.markers),
        totals=totals,
    )
