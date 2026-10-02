"""标准 10x MatrixMarket 稀疏输入读取与严格校验。

输入目录须含三个文件，三者必须同为未压缩或同为 gzip 单成员压缩：
- ``matrix.mtx``（或 ``matrix.mtx.gz``）：MatrixMarket 坐标格式，仅接受
  ``coordinate integer general`` 或 ``coordinate real general``；
  行是基因、列是细胞、索引从 1 开始；零值可省略，显式值必须能
  无损解析为有限非负整数；
- ``barcodes.tsv``（或 ``barcodes.tsv.gz``）：每个非空行是一个细胞 ID，
  按文件行序进入分析；
- ``features.tsv``（或 ``features.tsv.gz``）：每个非空数据行取第一列
  作为基因 ID，按文件行序进入分析。

压缩只改变承载方式，解压后文本口径与未压缩输入一致。声明维度必须与
ID 数量匹配；缺文件、压缩与未压缩混用、gzip 无法无损还原、格式头不支持、
越界索引、重复坐标、空或重复 ID、非法数值一律抛
:class:`cell_flow.errors.CellFlowInputError`。
"""

import hashlib
import os
from decimal import Decimal, InvalidOperation
from typing import Dict, List, Optional, Tuple

from .errors import CellFlowInputError
from .io import ExpressionMatrix, InputFile, _decompress_gzip

MATRIX_NAME = "matrix.mtx"
BARCODES_NAME = "barcodes.tsv"
FEATURES_NAME = "features.tsv"
MEMBER_NAMES = (MATRIX_NAME, BARCODES_NAME, FEATURES_NAME)
GZIP_SUFFIX = ".gz"

_BANNER = "%%matrixmarket"


def _fail(message: str):
    raise CellFlowInputError(message)


def _resolve_members(directory: str) -> Tuple[List[str], bool]:
    """确定三个成员的实际文件名与是否 gzip 承载。

    三者必须同为未压缩原名或同为 ``.gz`` 名；混用（含同一角色原名与
    ``.gz`` 并存）或缺文件一律抛
    :class:`cell_flow.errors.CellFlowInputError`。
    """
    plain_present = [
        os.path.exists(os.path.join(directory, name)) for name in MEMBER_NAMES
    ]
    gz_present = [
        os.path.exists(os.path.join(directory, name + GZIP_SUFFIX))
        for name in MEMBER_NAMES
    ]
    if all(plain_present) and not any(gz_present):
        return list(MEMBER_NAMES), False
    if all(gz_present) and not any(plain_present):
        return [name + GZIP_SUFFIX for name in MEMBER_NAMES], True
    if any(plain_present) and any(gz_present):
        _fail(
            "MTX 三个输入文件必须同为未压缩或同为 gzip："
            f"{directory} 中混用了两种形式"
        )
    # 无混用但存在缺失：沿用“缺少文件”口径；目录中已出现 .gz 成员时
    # 按 gzip 形式报告缺失名，否则按原名报告
    names = MEMBER_NAMES
    present = plain_present
    if any(gz_present):
        names = tuple(name + GZIP_SUFFIX for name in MEMBER_NAMES)
        present = gz_present
    for name, ok in zip(names, present):
        if not ok:
            _fail(f"MTX 输入缺少文件：{os.path.join(directory, name)}")
    # 不可达：上面的分支已覆盖全部存在性组合
    raise AssertionError("成员存在性组合未被覆盖")


def _read_file(directory: str, name: str, compressed: bool) -> Tuple[bytes, str]:
    path = os.path.join(directory, name)
    if not os.path.exists(path):
        _fail(f"MTX 输入缺少文件：{path}")
    if not os.path.isfile(path):
        _fail(f"MTX 输入路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"MTX 输入文件不可读：{path}")
    try:
        with open(path, "rb") as handle:
            raw = handle.read()
    except OSError as exc:
        _fail(f"MTX 输入文件不可读：{path}（{exc}）")
    payload = _decompress_gzip(raw, path) if compressed else raw
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"MTX 输入文件不是合法的 UTF-8 文本：{path}")
    return raw, text


def _read_id_lines(
    text: str,
    path: str,
    *,
    first_column_only: bool,
    kind: str,
) -> List[str]:
    """读取 ID 列表。

    barcodes：每个非空行整体是一个细胞 ID；
    features：每个非空数据行取第一列作为基因 ID。
    空行跳过，空 ID 与重复 ID 均非法。ID 按文件行序排列。
    """
    ids: List[str] = []
    seen = set()
    for line_no, line in enumerate(text.splitlines(), start=1):
        if line == "":
            continue
        token = line.split("\t", 1)[0] if first_column_only else line
        if token == "":
            _fail(f"{path} 第 {line_no} 行{kind} ID 为空")
        if token in seen:
            _fail(f"{kind} ID 重复：{token!r}")
        seen.add(token)
        ids.append(token)
    return ids


def _parse_nonnegative_integer(
    token: str, field_type: str, entry_no: int, matrix_name: str
) -> Optional[int]:
    """把一个显式矩阵元素解析为有限非负整数，要求无损。

    integer 字段只接受十进制非负整数字面量；real 字段接受浮点写法，
    但数值必须有限、非负且与其取整结果严格相等（无损）。
    返回 None 表示该元素是可省略的零。
    """
    if field_type == "integer":
        if not _is_decimal_uint(token):
            _fail(
                f"{matrix_name} 第 {entry_no} 个数据元素 {token!r} "
                f"不是非负整数"
            )
        value = int(token)
    else:  # real
        # 用 Decimal 精确解析，避免经 float64 在 2^53 以上发生舍入而
        # 破坏“无损”；显式值必须有限、非负且恰为整数
        try:
            number = Decimal(token)
        except InvalidOperation:
            number = Decimal("nan")
        if not number.is_finite():
            _fail(
                f"{matrix_name} 第 {entry_no} 个数据元素 {token!r} "
                f"不是有限数值"
            )
        if number < 0:
            _fail(
                f"{matrix_name} 第 {entry_no} 个数据元素 {token!r} 为负数"
            )
        if number != number.to_integral_value():
            _fail(
                f"{matrix_name} 第 {entry_no} 个数据元素 {token!r} "
                f"不是整数值，无法无损解析为计数"
            )
        value = int(number)
    if value == 0:
        return None
    return value


def _is_decimal_uint(token: str) -> bool:
    if not token:
        return False
    # 拒绝前导符号、小数点、指数、空白等；仅接受 0 或非零开头的十进制数
    if token == "0":
        return True
    if token[0] == "0" or not token[0].isdigit():
        return False
    return all(ch.isdigit() for ch in token)


def _parse_matrix(
    text: str, path: str, matrix_name: str
) -> Tuple[int, int, Dict[Tuple[int, int], int]]:
    lines = text.splitlines()

    # 1) 校验 banner 头
    if not lines:
        _fail(f"{path} 为空，缺少 MatrixMarket 头部")
    header_tokens = lines[0].split()
    if not header_tokens or header_tokens[0].lower() != _BANNER:
        _fail(
            f"{path} 格式头不支持：首行必须以 {_BANNER!r} 开头，"
            f"得到 {lines[0]!r}"
        )
    header = [t.lower() for t in header_tokens]
    # 标准 object/format/field/symmetry 五个关键字
    if len(header) != 5:
        _fail(f"{path} 格式头不支持：关键字数量不是 5：{lines[0]!r}")
    _, object_, format_, field_type, symmetry = header
    if object_ != "matrix":
        _fail(f"{path} 格式头不支持：object 必须是 matrix，得到 {object_!r}")
    if format_ != "coordinate":
        _fail(
            f"{path} 格式头不支持：仅接受 coordinate 格式，得到 {format_!r}"
        )
    if field_type not in ("integer", "real"):
        _fail(
            f"{path} 格式头不支持：仅接受 integer 或 real 字段，"
            f"得到 {field_type!r}"
        )
    if symmetry != "general":
        _fail(
            f"{path} 格式头不支持：仅接受 general 对称性，得到 {symmetry!r}"
        )

    # 2) 跳过注释、找到尺寸行
    body_index = 1
    dimensions: Optional[List[str]] = None
    while body_index < len(lines):
        stripped = lines[body_index].strip()
        body_index += 1
        if stripped == "":
            continue
        if stripped.startswith("%"):
            continue
        dimensions = stripped.split()
        break
    if dimensions is None:
        _fail(f"{path} 缺少维度声明行（rows columns entries）")
    if len(dimensions) != 3:
        _fail(
            f"{path} 维度声明行需要 3 个整数（rows columns entries），"
            f"得到 {dimensions!r}"
        )
    dims: List[int] = []
    for label, token in zip(("rows", "columns", "entries"), dimensions):
        if not _is_decimal_uint(token) or int(token) < 0:
            _fail(f"{path} 维度声明 {label} 不是非负整数：{token!r}")
        dims.append(int(token))
    n_rows, n_cols, n_entries = dims

    # 3) 逐行解析坐标元素
    # entries 只存非零值；seen_coords 记录全部已出现坐标（含显式零），
    # 以便对零值重复坐标同样报“重复坐标”
    entries: Dict[Tuple[int, int], int] = {}
    seen_coords = set()
    entry_no = 0
    data_seen = 0
    # lines[body_index] 即维度行之后的第一行，其 1 基行号为 body_index + 1
    for line_no, raw_line in enumerate(lines[body_index:], start=body_index + 1):
        stripped = raw_line.strip()
        if stripped == "" or stripped.startswith("%"):
            continue
        entry_no += 1
        parts = stripped.split()
        if len(parts) != 3:
            _fail(
                f"{path} 第 {line_no} 行数据需要 3 列"
                f"（row column value），得到 {len(parts)} 列"
            )
        row_tok, col_tok, val_tok = parts
        if not _is_decimal_uint(row_tok) or not _is_decimal_uint(col_tok):
            _fail(
                f"{path} 第 {line_no} 行坐标不是正整数："
                f"{row_tok!r}, {col_tok!r}"
            )
        row_idx = int(row_tok)
        col_idx = int(col_tok)
        if row_idx < 1 or row_idx > n_rows:
            _fail(
                f"{path} 第 {entry_no} 个数据元素行索引 {row_idx} 越界"
                f"（声明行数 {n_rows}）"
            )
        if col_idx < 1 or col_idx > n_cols:
            _fail(
                f"{path} 第 {entry_no} 个数据元素列索引 {col_idx} 越界"
                f"（声明列数 {n_cols}）"
            )
        value = _parse_nonnegative_integer(val_tok, field_type, entry_no, matrix_name)
        key = (row_idx, col_idx)
        if key in seen_coords:
            _fail(
                f"{path} 出现重复坐标（行 {row_idx}、列 {col_idx}）"
            )
        seen_coords.add(key)
        data_seen += 1
        if value is not None:
            entries[key] = value

    if data_seen != n_entries:
        _fail(
            f"{path} 实际数据元素 {data_seen} 个与声明 {n_entries} 个不一致"
        )

    return n_rows, n_cols, entries


def read_mtx_directory(directory: str) -> ExpressionMatrix:
    """读取并严格校验一个 10x MatrixMarket 目录。"""
    if directory is None or directory == "":
        _fail("输入路径为空")
    if not os.path.exists(directory):
        _fail(f"MTX 输入目录不存在：{directory}")
    if not os.path.isdir(directory):
        _fail(f"MTX 输入路径不是目录：{directory}")

    member_names, compressed = _resolve_members(directory)
    matrix_name, barcodes_name, features_name = member_names

    matrix_raw, matrix_text = _read_file(directory, matrix_name, compressed)
    barcodes_raw, barcodes_text = _read_file(directory, barcodes_name, compressed)
    features_raw, features_text = _read_file(directory, features_name, compressed)

    cell_ids = _read_id_lines(
        barcodes_text,
        os.path.join(directory, barcodes_name),
        first_column_only=False,
        kind="细胞",
    )
    gene_ids = _read_id_lines(
        features_text,
        os.path.join(directory, features_name),
        first_column_only=True,
        kind="基因",
    )
    if not cell_ids:
        _fail(f"{barcodes_name} 没有任何非空细胞 ID 行")
    if not gene_ids:
        _fail(f"{features_name} 没有任何非空基因 ID 数据行")

    n_rows, n_cols, entries = _parse_matrix(
        matrix_text, os.path.join(directory, matrix_name), matrix_name
    )
    if n_rows != len(gene_ids):
        _fail(
            f"{matrix_name} 声明行数 {n_rows} 与 {features_name} 基因数 "
            f"{len(gene_ids)} 不一致"
        )
    if n_cols != len(cell_ids):
        _fail(
            f"{matrix_name} 声明列数 {n_cols} 与 {barcodes_name} 细胞数 "
            f"{len(cell_ids)} 不一致"
        )

    # 以基因/ID 行序展开为稠密计数矩阵，缺失坐标即零
    counts: List[List[int]] = [[0] * n_cols for _ in range(n_rows)]
    total = 0
    for (row_idx, col_idx), value in entries.items():
        counts[row_idx - 1][col_idx - 1] = value
        total += value

    files = [
        InputFile(name=matrix_name, sha256=hashlib.sha256(matrix_raw).hexdigest()),
        InputFile(name=barcodes_name, sha256=hashlib.sha256(barcodes_raw).hexdigest()),
        InputFile(name=features_name, sha256=hashlib.sha256(features_raw).hexdigest()),
    ]
    return ExpressionMatrix(
        gene_ids=gene_ids,
        cell_ids=cell_ids,
        counts=counts,
        sha256=files[0].sha256,
        path=directory,
        total_counts=total,
        input_format="mtx",
        files=files,
    )
