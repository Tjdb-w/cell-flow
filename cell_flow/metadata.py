"""细胞分组元数据读取与严格校验。

``--metadata`` 指向 UTF-8 制表符文本，表头恰为 ``cell_id`` 与 ``group``
两列；``cell_id`` 唯一且与表达矩阵的全部细胞一一对应，``group`` 非空。
文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名无关）。
任何不合法（表头不符、ID 重复、覆盖不全、多出矩阵之外的细胞、空分组，
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
class CellMetadata:
    """细胞分组：每个矩阵细胞恰好归属一个非空分组。"""

    path: str
    name: str                  # 文件名（路径末段）
    sha256: str                # 原始字节（gzip 即压缩字节）的 SHA-256
    groups: Dict[str, str]     # cell_id -> group
    group_sizes: Dict[str, int]  # group -> 细胞数（按 group 升序）


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_metadata(path: str, cell_ids: List[str]) -> CellMetadata:
    """读取并校验分组元数据；``cell_ids`` 为表达矩阵的全部细胞 ID（列序）。"""
    if path is None or path == "":
        _fail("元数据路径为空")
    if not os.path.exists(path):
        _fail(f"元数据文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"元数据路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"元数据文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"元数据文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"元数据文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("元数据表头缺失：文件为空")
    header = lines[0].split("\t")
    if header != ["cell_id", "group"]:
        _fail("元数据表头必须恰为 cell_id 和 group 两列")

    groups: Dict[str, str] = {}
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"元数据第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != 2:
            _fail(
                f"元数据第 {offset} 行列数为 {len(fields)}，与表头两列不一致"
            )
        cell_id, group = fields
        if cell_id == "":
            _fail(f"元数据第 {offset} 行细胞 ID 为空")
        if group == "":
            _fail(f"元数据第 {offset} 行分组为空（细胞 {cell_id!r}）")
        if cell_id in groups:
            _fail(f"元数据细胞 ID 重复：{cell_id!r}")
        groups[cell_id] = group

    missing = [cell_id for cell_id in cell_ids if cell_id not in groups]
    if missing:
        _fail(
            f"元数据未覆盖矩阵全部细胞：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个"
        )
    matrix_cells = set(cell_ids)
    extra = [cell_id for cell_id in groups if cell_id not in matrix_cells]
    if extra:
        _fail(
            f"元数据包含矩阵之外的细胞：{extra[0]!r} 等 {len(extra)} 个"
        )

    sizes: Dict[str, int] = {}
    for group in groups.values():
        sizes[group] = sizes.get(group, 0) + 1
    group_sizes = {group: sizes[group] for group in sorted(sizes)}

    return CellMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        groups=groups,
        group_sizes=group_sizes,
    )
