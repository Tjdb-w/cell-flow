"""标准 10x MatrixMarket 稀疏输入读取与严格校验。

输入目录须含三个文件，每个文件可保持原名，也可整体改用 gzip 压缩名
（三个文件必须同为未压缩或同为 gzip，混用报输入错误）：
- ``matrix.mtx`` / ``matrix.mtx.gz``：MatrixMarket 坐标格式，仅接受
  ``coordinate integer general`` 或 ``coordinate real general``；
  行是基因、列是细胞、索引从 1 开始；零值可省略，显式值必须能
  无损解析为有限非负整数；
- ``barcodes.tsv`` / ``barcodes.tsv.gz``：每个非空行是一个细胞 ID，
  按文件行序进入分析；
- ``features.tsv`` / ``features.tsv.gz``：每个非空数据行取第一列作为
  基因 ID，按文件行序进入分析。

gzip 载体仅允许单成员且无尾随数据；gzip 头、CRC、长度校验失败或流被
截断一律报输入错误。声明维度必须与 ID 数量匹配；缺文件、两种文件名
并存、格式头不支持、越界索引、重复坐标、空或重复 ID、非法数值一律
抛 :class:`cell_flow.errors.CellFlowInputError`。
"""

import hashlib
import os
from decimal import Decimal, InvalidOperation
from typing import Dict, List, Optional, Tuple

from .errors import CellFlowInputError
from .io import GZIP_MAGIC, ExpressionMatrix, InputFile, gunzip_single

MATRIX_NAME = "matrix.mtx"
BARCODES_NAME = "barcodes.tsv"
FEATURES_NAME = "features.tsv"
GZIP_SUFFIX = ".gz"

_BANNER = "%%matrixmarket"


def _fail(message: str):
    raise CellFlowInputError(message)


def _resolve_name(directory: str, base_name: str) -> Tuple[str, bool]:
    """在目录中解析 ``<name>`` 或 ``<name>.gz``，二者必须恰好存在一个。

    返回（实际文件名，是否 gzip）。此步只做布局校验，不读取文件内容，
    以保证缺文件与“压缩/未压缩混用”的报错先于任何内容损坏报错。
    """
    gz_name = base_name + GZIP_SUFFIX
    plain_path = os.path.join(directory, base_name)
    gz_path = os.path.join(directory, gz_name)
    plain_exists = os.path.exists(plain_path)
    gz_exists = os.path.exists(gz_path)
    if plain_exists and gz_exists:
        _fail(
            f"MTX 输入同时存在 {base_name} 与 {gz_name}，"
            f"二者只能保留其一：{directory}"
        )
    if gz_exists:
        return gz_name, True
    if plain_exists:
        return base_name, False
    _fail(f"MTX 输入缺少文件：{plain_path}（或 {gz_name}）")


def _read_named_file(
    directory: str, name: str, *, compressed: bool
) -> Tuple[bytes, str]:
    path = os.path.join(directory, name)
    if not os.path.isfile(path):
        _fail(f"MTX 输入路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"MTX 输入文件不可读：{path}")
    try:
        with open(path, "rb") as handle:
            raw = handle.read()
    except OSError as exc:
        _fail(f"MTX 输入文件不可读：{path}（{exc}）")

    if compressed:
        # .gz 是压缩状态的显式声明：内容必须确实是 gzip，由
        # gunzip_single 严格校验魔数、CRC、长度、截断与单成员
        payload = gunzip_single(raw, path)
    elif raw[:2] == GZIP_MAGIC:
        # 未压缩名下出现 gzip 字节属于布局不合法，避免静默按内容嗅探
        _fail(
            f"MTX 输入文件 {name} 的内容是 gzip，但文件名未带 "
            f"{GZIP_SUFFIX} 后缀：{path}"
        )
    else:
        payload = raw

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
    allow_duplicates: bool = False,
) -> List[str]:
    """读取 ID 列表。

    barcodes：每个非空行整体是一个细胞 ID；
    features：每个非空数据行取第一列作为基因 ID。
    空行跳过，空 ID 非法；默认重复 ID 也非法，``allow_duplicates``
    为真时保留重复（供批次功能读取后统一按 :class:`ValueError` 报告）。
    ID 按文件行序排列。
    """
    ids: List[str] = []
    seen = set()
    for line_no, line in enumerate(text.splitlines(), start=1):
        if line == "":
            continue
        token = line.split("\t", 1)[0] if first_column_only else line
        if token == "":
            _fail(f"{path} 第 {line_no} 行{kind} ID 为空")
        if token in seen and not allow_duplicates:
            _fail(f"{kind} ID 重复：{token!r}")
        seen.add(token)
        ids.append(token)
    return ids


def _parse_nonnegative_integer(
    token: str, field_type: str, entry_no: int
) -> Optional[int]:
    """把一个显式矩阵元素解析为有限非负整数，要求无损。

    integer 字段只接受十进制非负整数字面量；real 字段接受浮点写法，
    但数值必须有限、非负且与其取整结果严格相等（无损）。
    返回 None 表示该元素是可省略的零。
    """
    if field_type == "integer":
        if not _is_decimal_uint(token):
            _fail(
                f"{MATRIX_NAME} 第 {entry_no} 个数据元素 {token!r} "
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
                f"{MATRIX_NAME} 第 {entry_no} 个数据元素 {token!r} "
                f"不是有限数值"
            )
        if number < 0:
            _fail(
                f"{MATRIX_NAME} 第 {entry_no} 个数据元素 {token!r} 为负数"
            )
        if number != number.to_integral_value():
            _fail(
                f"{MATRIX_NAME} 第 {entry_no} 个数据元素 {token!r} "
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
    text: str, path: str
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
        value = _parse_nonnegative_integer(val_tok, field_type, entry_no)
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


def read_mtx_directory(
    directory: str, *, allow_duplicate_cells: bool = False
) -> ExpressionMatrix:
    """读取并严格校验一个 10x MatrixMarket 目录。

    ``allow_duplicate_cells`` 为真时保留 barcodes 中的重复细胞条码
    （供批次功能读取后统一按 :class:`ValueError` 报告条码冲突）。
    """
    if directory is None or directory == "":
        _fail("输入路径为空")
    if not os.path.exists(directory):
        _fail(f"MTX 输入目录不存在：{directory}")
    if not os.path.isdir(directory):
        _fail(f"MTX 输入路径不是目录：{directory}")

    # 先只解析三个文件名并校验“同为压缩或同为未压缩”，
    # 使布局错误先于任何文件内容错误暴露
    matrix_actual, matrix_gz = _resolve_name(directory, MATRIX_NAME)
    barcodes_actual, barcodes_gz = _resolve_name(directory, BARCODES_NAME)
    features_actual, features_gz = _resolve_name(directory, FEATURES_NAME)
    if not (matrix_gz == barcodes_gz == features_gz):
        _fail(
            f"MTX 三个输入文件必须同为未压缩或同为 gzip，当前混用："
            f"{matrix_actual}、{barcodes_actual}、{features_actual}"
        )

    matrix_raw, matrix_text = _read_named_file(
        directory, matrix_actual, compressed=matrix_gz
    )
    barcodes_raw, barcodes_text = _read_named_file(
        directory, barcodes_actual, compressed=barcodes_gz
    )
    features_raw, features_text = _read_named_file(
        directory, features_actual, compressed=features_gz
    )

    cell_ids = _read_id_lines(
        barcodes_text,
        os.path.join(directory, barcodes_actual),
        first_column_only=False,
        kind="细胞",
        allow_duplicates=allow_duplicate_cells,
    )
    gene_ids = _read_id_lines(
        features_text,
        os.path.join(directory, features_actual),
        first_column_only=True,
        kind="基因",
    )
    if not cell_ids:
        _fail(f"{barcodes_actual} 没有任何非空细胞 ID 行")
    if not gene_ids:
        _fail(f"{features_actual} 没有任何非空基因 ID 数据行")

    n_rows, n_cols, entries = _parse_matrix(
        matrix_text, os.path.join(directory, matrix_actual)
    )
    if n_rows != len(gene_ids):
        _fail(
            f"{matrix_actual} 声明行数 {n_rows} 与 {features_actual} 基因数 "
            f"{len(gene_ids)} 不一致"
        )
    if n_cols != len(cell_ids):
        _fail(
            f"{matrix_actual} 声明列数 {n_cols} 与 {barcodes_actual} 细胞数 "
            f"{len(cell_ids)} 不一致"
        )

    # 以基因/ID 行序展开为稠密计数矩阵，缺失坐标即零
    counts: List[List[int]] = [[0] * n_cols for _ in range(n_rows)]
    total = 0
    for (row_idx, col_idx), value in entries.items():
        counts[row_idx - 1][col_idx - 1] = value
        total += value

    files = [
        InputFile(name=matrix_actual, sha256=hashlib.sha256(matrix_raw).hexdigest()),
        InputFile(name=barcodes_actual, sha256=hashlib.sha256(barcodes_raw).hexdigest()),
        InputFile(name=features_actual, sha256=hashlib.sha256(features_raw).hexdigest()),
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
