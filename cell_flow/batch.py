"""批次元数据读取与严格校验。

``--batch-metadata`` 指向 UTF-8 制表符文本，表头恰为 ``cell_id`` 与
``batch`` 两列；``cell_id`` 唯一且与表达矩阵的全部细胞一一对应，
``batch`` 非空。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与
文件名无关）。任何不合法（空文件、不可读、编码错误、表头不符、
行宽不一致、空 batch、空/重复 cell_id、覆盖不全或超出矩阵细胞，
或 gzip 多成员、尾随数据、截断、校验失败）都抛
:class:`cell_flow.errors.CellFlowInputError`，且不改动任何已有结果。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single


@dataclass(frozen=True)
class BatchMetadata:
    """批次归属：每个矩阵细胞恰好归属一个非空批次。"""

    path: str
    name: str                  # 文件名（路径末段）
    sha256: str                # 原始字节（gzip 即压缩字节）的 SHA-256
    assignments: Dict[str, str]  # cell_id -> batch


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

    assignments: Dict[str, str] = {}
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"批次元数据第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != 2:
            _fail(
                f"批次元数据第 {offset} 行列数为 {len(fields)}，"
                f"与表头两列不一致"
            )
        cell_id, batch = fields
        if cell_id == "":
            _fail(f"批次元数据第 {offset} 行细胞 ID 为空")
        if batch == "":
            _fail(f"批次元数据第 {offset} 行批次为空（细胞 {cell_id!r}）")
        if cell_id in assignments:
            _fail(f"批次元数据细胞 ID 重复：{cell_id!r}")
        assignments[cell_id] = batch

    missing = [cell_id for cell_id in cell_ids if cell_id not in assignments]
    if missing:
        _fail(
            f"批次元数据未覆盖矩阵全部细胞：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个"
        )
    extra = [cell_id for cell_id in assignments if cell_id not in set(cell_ids)]
    if extra:
        _fail(
            f"批次元数据包含矩阵之外的细胞：{extra[0]!r} 等 {len(extra)} 个"
        )

    return BatchMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        assignments=assignments,
    )
