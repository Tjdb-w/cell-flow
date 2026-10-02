"""批次元数据读取、按批次的基因均值中心化校正，以及批次汇总。

``--batch-metadata`` 指向 UTF-8 制表符文本，表头恰为 ``cell_id`` 与
``batch`` 两列；``cell_id`` 唯一且与表达矩阵的全部细胞一一对应，
``batch`` 非空。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与
文件名无关）。任何不合法（表头不符、ID 重复、覆盖不全、多出矩阵之外
的细胞、空批次，或 gzip 多成员、尾随数据、截断、校验失败）都抛
:class:`cell_flow.errors.CellFlowInputError`，且不改动任何已有结果。

校正按基因进行：每个 log 归一化值减去对应批次保留细胞的基因均值，
再加回全部保留细胞的基因总均值（batch mean centering）。校正值进入
高变基因选择、PCA、聚类、markers、成对 markers、分组差异表达与图表
数据；``normalized_expression.tsv`` 仍写未校正值。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single


@dataclass(frozen=True)
class BatchMetadata:
    """细胞批次：每个矩阵细胞恰好归属一个非空批次。"""

    path: str
    name: str                  # 文件名（路径末段）
    sha256: str                # 原始字节（gzip 即压缩字节）的 SHA-256
    batches: Dict[str, str]    # cell_id -> batch
    batch_sizes: Dict[str, int]  # batch -> 细胞数（按 batch 升序）


@dataclass(frozen=True)
class BatchSummaryRow:
    """一个批次在质控后保留细胞上的汇总（均值按批次内保留细胞计算）。"""

    batch_id: str
    n_cells: int
    mean_total_counts: float
    mean_detected_genes: float
    mean_mitochondrial_fraction: float


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_batch_metadata(path: str, cell_ids: List[str]) -> BatchMetadata:
    """读取并校验批次元数据；``cell_ids`` 为表达矩阵的全部细胞 ID（列序）。"""
    if path is None or path == "":
        _fail("批次元数据路径为空")
    if not os.path.exists(path):
        _fail(f"批次元数据文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"批次元数据路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"批次元数据文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"批次元数据文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"批次元数据文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("批次元数据表头缺失：文件为空")
    header = lines[0].split("\t")
    if header != ["cell_id", "batch"]:
        _fail("批次元数据表头必须恰为 cell_id 和 batch 两列")

    batches: Dict[str, str] = {}
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"批次元数据第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != 2:
            _fail(
                f"批次元数据第 {offset} 行列数为 {len(fields)}，与表头两列不一致"
            )
        cell_id, batch = fields
        if cell_id == "":
            _fail(f"批次元数据第 {offset} 行细胞 ID 为空")
        if batch == "":
            _fail(f"批次元数据第 {offset} 行批次为空（细胞 {cell_id!r}）")
        if cell_id in batches:
            _fail(f"批次元数据细胞 ID 重复：{cell_id!r}")
        batches[cell_id] = batch

    missing = [cell_id for cell_id in cell_ids if cell_id not in batches]
    if missing:
        _fail(
            f"批次元数据未覆盖矩阵全部细胞：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个"
        )
    matrix_cells = set(cell_ids)
    extra = [cell_id for cell_id in batches if cell_id not in matrix_cells]
    if extra:
        _fail(
            f"批次元数据包含矩阵之外的细胞：{extra[0]!r} 等 {len(extra)} 个"
        )

    sizes: Dict[str, int] = {}
    for batch in batches.values():
        sizes[batch] = sizes.get(batch, 0) + 1
    batch_sizes = {batch: sizes[batch] for batch in sorted(sizes)}

    return BatchMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        batches=batches,
        batch_sizes=batch_sizes,
    )


def center_by_batch(
    values: List[List[float]], batch_labels: List[str]
) -> List[List[float]]:
    """按基因做批次均值中心化。

    ``values[g][c]`` 为保留基因 × 保留细胞的 log 归一化表达，
    ``batch_labels[c]`` 为对应细胞的批次。校正值 =
    原值 - 该细胞所在批次的基因均值 + 全部保留细胞的基因总均值。
    均值按列序累加，结果确定。
    """
    n_cells = len(batch_labels)
    corrected: List[List[float]] = []
    for row in values:
        grand_mean = sum(row) / n_cells
        batch_sums: Dict[str, float] = {}
        batch_counts: Dict[str, int] = {}
        for c, batch in enumerate(batch_labels):
            batch_sums[batch] = batch_sums.get(batch, 0.0) + row[c]
            batch_counts[batch] = batch_counts.get(batch, 0) + 1
        batch_means = {
            batch: batch_sums[batch] / batch_counts[batch]
            for batch in batch_sums
        }
        corrected.append(
            [
                row[c] - batch_means[batch_labels[c]] + grand_mean
                for c in range(n_cells)
            ]
        )
    return corrected
