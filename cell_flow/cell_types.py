"""细胞类型注释：参考标记读入与簇评分。

``--cell-type-reference`` 指向 UTF-8 制表符文本，表头恰为 ``cell_type``、
``gene_id`` 两列，每行一个“类型-标记”成员关系；字段非空、
(cell_type, gene_id) 组合唯一，至少一条数据行。文件可为纯文本或单成员
gzip（按 gzip 魔数识别，与文件名无关）。基因 ID 与表达矩阵首列精确匹配，
不做大小写、别名、前缀转换；矩阵外或被质控剔除的标记不计分但计入总数。
任何不合法（空路径、UTF-8、表头、空值、重复、数据行不足，或 gzip 多成员、
尾随数据、截断、CRC/长度错误）都抛
:class:`cell_flow.errors.CellFlowInputError`（退出码 2），且不改动结果。

评分只用质控后保留细胞与基因的表达值（有批次校正时取实际用于聚类的批次
均值中心化值，否则取 log 归一化值）与最终簇标签：对某簇某类型，逐标记
基因求“簇内均值减簇外（其余全部簇）均值”，再对可用标记取算术平均。
每簇取最高分类型，并列按 cell_type 的 Unicode 码点升序取第一；无可用
标记的类型不参与竞争，允许多个簇注释为同一类型。最高分不大于 0 时状态为
``unassigned``（仍保留候选类型与分数），否则为 ``assigned``。全部类型的
标记与保留基因都没有交集时抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
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
    """细胞类型标记参考，按文件行序保留每个类型的标记基因。"""

    path: str
    name: str                    # 文件名（路径末段）
    sha256: str                  # 原始字节（gzip 即压缩字节）的 SHA-256
    markers: Dict[str, List[str]]  # cell_type -> 标记基因 ID（文件行序）


@dataclass(frozen=True)
class ClusterAnnotation:
    """一个簇的注释结果；unassigned 时同样保留候选类型与分数。"""

    cluster: int
    n_cells: int
    cell_type: str
    status: str               # STATUS_ASSIGNED 或 STATUS_UNASSIGNED
    score: float
    n_markers_total: int      # 候选类型在参考中的标记总数
    n_markers_used: int       # 其中与保留基因匹配、实际参与评分的标记数


@dataclass(frozen=True)
class CellTypeAnnotation:
    """全部簇的注释结果；``rows`` 按 cluster 升序排列。"""

    reference: CellTypeReference
    rows: List[ClusterAnnotation]
    by_cluster: Dict[int, ClusterAnnotation]


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
                f"细胞类型参考第 {offset} 行列数为 {len(fields)}，与表头两列不一致"
            )
        cell_type, gene_id = fields
        if cell_type == "":
            _fail(f"细胞类型参考第 {offset} 行 cell_type 为空")
        if gene_id == "":
            _fail(f"细胞类型参考第 {offset} 行 gene_id 为空")
        pair = (cell_type, gene_id)
        if pair in seen_pairs:
            _fail(
                f"细胞类型参考成员关系重复：cell_type={cell_type!r}、gene_id={gene_id!r}"
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


def annotate_cell_types(
    reference: CellTypeReference,
    *,
    kept_gene_ids: List[str],
    analysis_values: List[List[float]],
    labels: List[int],
) -> CellTypeAnnotation:
    """按参考标记为每个最终簇评分并选出候选细胞类型。

    ``analysis_values[g][c]`` 与 ``kept_gene_ids``、``labels``（最终簇标签，
    保留细胞列序）对齐，为批次中心化值（实际施加批次校正时）或 log 归一化
    值。只用与保留基因匹配的标记；矩阵外或被剔除的标记只计入
    ``n_markers_total``。全部类型的标记与保留基因交集都为空时抛数据错误
    （退出码 4），不产出任何结果。
    """
    gene_position: Dict[str, int] = {
        gene_id: g for g, gene_id in enumerate(kept_gene_ids)
    }

    # 每个类型的可用标记（保留基因行下标，沿用参考文件行序）；
    # 无可用标记的类型不参与竞争
    used_positions: Dict[str, List[int]] = {}
    for cell_type in sorted(reference.markers):
        positions = [
            gene_position[gene_id]
            for gene_id in reference.markers[cell_type]
            if gene_id in gene_position
        ]
        if positions:
            used_positions[cell_type] = positions
    if not used_positions:
        raise CellFlowDataError(
            "细胞类型参考的标记与质控后保留基因的交集为空，无法注释"
        )

    n_cells = len(labels)
    cluster_members: Dict[int, List[int]] = {}
    for c, label in enumerate(labels):
        cluster_members.setdefault(label, []).append(c)

    rows: List[ClusterAnnotation] = []
    by_cluster: Dict[int, ClusterAnnotation] = {}
    for cluster in sorted(cluster_members):
        members = cluster_members[cluster]
        others = [c for c in range(n_cells) if labels[c] != cluster]
        n_in = len(members)
        n_out = len(others)

        # 每个相关基因先求“簇内均值减其余簇均值”，再按类型对可用标记取平均；
        # 沿用管线既有口径（内置 sum）累加，保证与归一化/批次校正同精度
        needed = sorted({p for ps in used_positions.values() for p in ps})
        diffs: Dict[int, float] = {}
        for p in needed:
            row = analysis_values[p]
            mean_in = sum(row[c] for c in members) / n_in
            mean_out = sum(row[c] for c in others) / n_out
            diffs[p] = mean_in - mean_out

        # 按 cell_type 的 Unicode 码点升序竞争，严格更大才替换：
        # 分数并列时排最前的类型胜出
        best_type = None
        best_score = None
        for cell_type in sorted(used_positions):
            positions = used_positions[cell_type]
            score = sum(diffs[p] for p in positions) / len(positions)
            if best_score is None or score > best_score:
                best_type = cell_type
                best_score = score

        assert best_type is not None and best_score is not None
        status = (
            STATUS_ASSIGNED if best_score > 0.0 else STATUS_UNASSIGNED
        )
        row = ClusterAnnotation(
            cluster=cluster,
            n_cells=n_in,
            cell_type=best_type,
            status=status,
            score=best_score,
            n_markers_total=len(reference.markers[best_type]),
            n_markers_used=len(used_positions[best_type]),
        )
        rows.append(row)
        by_cluster[cluster] = row

    return CellTypeAnnotation(
        reference=reference,
        rows=rows,
        by_cluster=by_cluster,
    )
