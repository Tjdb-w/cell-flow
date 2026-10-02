"""细胞分组元数据读取与严格校验。

输入为制表符分隔文本：表头恰为 ``cell_id`` 与 ``group`` 两列，其后每行
给出一个细胞 ID 与其分组。文件承载方式与表达矩阵一致：UTF-8 纯文本或
单成员 gzip（按 gzip 魔数自动识别，与文件名无关）。

合法条件：
- 表头恰为 ``cell_id\\tgroup``，列数不多不少；
- 每个数据行恰为两列，``cell_id`` 非空且唯一，``group`` 非空；
- ``cell_id`` 集合与表达矩阵的全部细胞一一对应（不缺、不多）。

任何不合法情形（含不可读、非 UTF-8、gzip 校验失败、多成员、尾随数据）
都抛 :class:`cell_flow.errors.CellFlowInputError`。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single

HEADER_CELL_ID = "cell_id"
HEADER_GROUP = "group"


@dataclass(frozen=True)
class CellGroupMetadata:
    """细胞分组：``groups`` 把每个细胞 ID 映射到其（非空）分组名。

    ``sha256`` 按磁盘上的原始字节计算（gzip 即压缩字节）；``name`` 为
    调用时给定路径的文件名部分，用于在 run.json 中记录来源。
    """

    path: str
    name: str
    sha256: str
    groups: Dict[str, str]


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_metadata(path: str, expected_cell_ids: List[str]) -> CellGroupMetadata:
    if path is None or path == "":
        _fail("元数据路径为空")
    if not os.path.exists(path):
        _fail(f"元数据文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"元数据路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"元数据文件不可读：{path}")

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

    lines = text.splitlines()
    if not lines:
        _fail("元数据表头缺失：文件为空")

    header = lines[0].split("\t")
    if header != [HEADER_CELL_ID, HEADER_GROUP]:
        _fail(
            f"元数据表头必须恰为 {HEADER_CELL_ID}\\t{HEADER_GROUP}，"
            f"得到 {header!r}"
        )

    groups: Dict[str, str] = {}
    for offset, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != 2:
            _fail(
                f"元数据第 {offset} 行需要恰为 2 列（cell_id、group），"
                f"得到 {len(fields)} 列"
            )
        cell_id, group = fields
        if cell_id == "":
            _fail(f"元数据第 {offset} 行 cell_id 为空")
        if group == "":
            _fail(f"元数据第 {offset} 行（细胞 {cell_id!r}）group 为空")
        if cell_id in groups:
            _fail(f"元数据 cell_id 重复：{cell_id!r}")
        groups[cell_id] = group

    if not groups:
        _fail("元数据没有任何细胞分组数据行")

    expected = set(expected_cell_ids)
    actual = set(groups)
    if actual != expected:
        missing = [c for c in expected_cell_ids if c not in actual]
        extra = [c for c in groups if c not in expected]
        parts = []
        if missing:
            parts.append(f"缺少 {len(missing)} 个细胞（如 {missing[0]!r}）")
        if extra:
            parts.append(f"出现 {len(extra)} 个矩阵外细胞（如 {extra[0]!r}）")
        _fail("元数据 cell_id 未覆盖矩阵全部细胞：" + "；".join(parts))

    return CellGroupMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        groups=groups,
    )
