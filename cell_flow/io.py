"""输入矩阵读取与严格校验。

支持两种输入格式：

- ``tsv``（默认）：单个制表符分隔矩阵，首列是唯一基因 ID，表头其余列名
  是唯一细胞 ID，取值为非负整数 UMI 计数；
- ``mtx``：标准 10x MatrixMarket 目录，须含 ``matrix.mtx``、
  ``barcodes.tsv``、``features.tsv`` 三个文件。``matrix.mtx`` 仅接受
  ``coordinate integer general`` 或 ``coordinate real general``，
  行＝基因、列＝细胞、索引从 1 开始。

两种读取器产出同构的 :class:`ExpressionMatrix`（基因行序、细胞列序、
逐格整数计数完全一致），任何不合法情形都抛出
:class:`cell_flow.errors.CellFlowInputError`。
"""

import hashlib
import os
import re
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from typing import List, Tuple

from .errors import CellFlowConfigError, CellFlowInputError

_INTEGER_RE = re.compile(r"^(0|[1-9][0-9]*)$")

MTX_FORMAT = "mtx"
TSV_FORMAT = "tsv"
INPUT_FORMATS = (TSV_FORMAT, MTX_FORMAT)

MTX_MATRIX_FILE = "matrix.mtx"
MTX_BARCODES_FILE = "barcodes.tsv"
MTX_FEATURES_FILE = "features.tsv"
MTX_FILE_NAMES = (MTX_MATRIX_FILE, MTX_BARCODES_FILE, MTX_FEATURES_FILE)


@dataclass(frozen=True)
class SourceFile:
    """一个输入文件的来源记录：文件名与其原始字节的 SHA-256。"""

    name: str
    sha256: str


@dataclass(frozen=True)
class ExpressionMatrix:
    """原始计数矩阵（行＝基因，列＝细胞）。"""

    gene_ids: List[str]
    cell_ids: List[str]
    counts: List[List[int]]  # counts[gene_index][cell_index]
    sha256: str
    path: str
    total_counts: int
    input_format: str = TSV_FORMAT
    source_files: Tuple[SourceFile, ...] = field(default_factory=tuple)

    @property
    def n_genes(self) -> int:
        return len(self.gene_ids)

    @property
    def n_cells(self) -> int:
        return len(self.cell_ids)


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def _parse_count(raw: str, gene_id: str, cell_id: str, row_no: int) -> int:
    if not _INTEGER_RE.match(raw):
        _fail(
            f"计数非法（第 {row_no} 行基因 {gene_id!r}、细胞 {cell_id!r}）："
            f"{raw!r} 不是非负整数"
        )
    return int(raw)


def read_matrix(path: str) -> ExpressionMatrix:
    """读取制表符分隔的稠密计数矩阵（基线格式）。"""
    if path is None or path == "":
        _fail("输入路径为空")
    if not os.path.exists(path):
        _fail(f"输入文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"输入路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"输入文件不可读：{path}")

    raw_bytes, text = _read_text(path)

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("表头缺失：文件为空")

    header = lines[0].split("\t")
    if len(header) < 2:
        _fail("表头缺失：至少需要基因 ID 列与一个细胞 ID 列")
    cell_ids = header[1:]
    if any(cell == "" for cell in cell_ids):
        _fail("表头存在空细胞 ID")
    if len(set(cell_ids)) != len(cell_ids):
        dup = _first_duplicate(cell_ids)
        _fail(f"细胞 ID 重复：{dup!r}")

    gene_ids: List[str] = []
    counts: List[List[int]] = []
    seen_genes = set()
    n_cells = len(cell_ids)

    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != n_cells + 1:
            _fail(
                f"第 {offset} 行列数为 {len(fields) - 1}，"
                f"与表头细胞数 {n_cells} 不一致"
            )
        gene_id = fields[0]
        if gene_id == "":
            _fail(f"第 {offset} 行基因 ID 为空")
        if gene_id in seen_genes:
            _fail(f"基因 ID 重复：{gene_id!r}")
        seen_genes.add(gene_id)
        row = [
            _parse_count(fields[c + 1], gene_id, cell_ids[c], offset)
            for c in range(n_cells)
        ]
        gene_ids.append(gene_id)
        counts.append(row)

    if not gene_ids:
        _fail("空矩阵：没有任何基因数据行")

    total = sum(sum(row) for row in counts)

    return ExpressionMatrix(
        gene_ids=gene_ids,
        cell_ids=cell_ids,
        counts=counts,
        sha256=hashlib.sha256(raw_bytes).hexdigest(),
        path=path,
        total_counts=total,
        input_format=TSV_FORMAT,
    )


def _read_text(path: str) -> Tuple[bytes, str]:
    """读取普通文件的原始字节并按 UTF-8 解码，失败统一为输入错误。"""
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"输入文件不可读：{path}（{exc}）")
    try:
        text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"输入文件不是合法的 UTF-8 文本：{path}")
    return raw_bytes, text


def _require_input_file(dir_path: str, name: str) -> Tuple[bytes, str]:
    file_path = os.path.join(dir_path, name)
    if not os.path.exists(file_path):
        _fail(f"mtx 输入目录缺少文件：{name}（目录：{dir_path}）")
    if not os.path.isfile(file_path):
        _fail(f"mtx 输入条目不是普通文件：{name}（目录：{dir_path}）")
    if not os.access(file_path, os.R_OK):
        _fail(f"mtx 输入文件不可读：{file_path}")
    return _read_text(file_path)


def _read_id_lines(
    text: str,
    file_name: str,
    *,
    first_column_only: bool,
) -> List[str]:
    """按文件行序读取 ID。

    - 空行跳过（“每个非空行”进入分析）；
    - ``first_column_only`` 为真时（features.tsv）取每行第一列；
    - 空 ID、重复 ID 均为输入错误。
    """
    ids: List[str] = []
    seen = set()
    for offset, line in enumerate(text.splitlines(), start=1):
        if line == "":
            continue
        identifier = line.split("\t", 1)[0] if first_column_only else line
        if identifier == "":
            _fail(f"{file_name} 第 {offset} 行存在空 ID")
        if identifier in seen:
            _fail(f"{file_name} 中 ID 重复：{identifier!r}")
        seen.add(identifier)
        ids.append(identifier)
    return ids


def _parse_nonnegative_int(raw: str, what: str) -> int:
    if not _INTEGER_RE.match(raw):
        _fail(f"{what} 必须是非负整数，得到 {raw!r}")
    return int(raw)


def _parse_mtx_value(raw: str, kind: str, line_no: int) -> int:
    """解析数据行中的计数值。

    integer 头仅接受严格非负整数字面；real 头接受任何“无损解析为有限
    非负整数”的写法（如 1、1.0、5e2），拒绝小数、负数、nan/inf。

    real 分支用 :class:`decimal.Decimal` 精确判定，避免大整数经 float
    中转时丢精度（如 99999999999999999999.0）。
    """
    if _INTEGER_RE.match(raw):
        return int(raw)
    if kind == "integer":
        _fail(
            f"matrix.mtx 第 {line_no} 行计数值 {raw!r} 不是非负整数"
            "（integer 格式不接受小数或指数写法）"
        )
    try:
        value = Decimal(raw)
    except InvalidOperation:
        value = Decimal("nan")
    # 显式负号（含 -0.0）与 TSV 读取器口径一致，一律拒绝：非负计数没有
    # 合法的负号写法
    if (
        raw[:1] == "-"
        or not value.is_finite()
        or value < 0
        or value != value.to_integral_value()
    ):
        _fail(
            f"matrix.mtx 第 {line_no} 行计数值 {raw!r} 非法："
            "显式值须无损解析为有限非负整数"
        )
    return int(value)


def _parse_mtx_matrix(
    text: str,
    n_genes: int,
    n_cells: int,
) -> List[List[int]]:
    lines = text.splitlines()
    if not lines:
        _fail("matrix.mtx 为空：缺少 MatrixMarket 格式头")

    banner = lines[0].strip()
    if not banner.startswith("%%MatrixMarket"):
        _fail(
            "matrix.mtx 格式头不支持：首行须以 %%MatrixMarket 开头，"
            f"得到 {lines[0]!r}"
        )
    parts = banner.split()
    if (
        len(parts) != 5
        or parts[1] != "matrix"
        or parts[2] != "coordinate"
        or parts[3] not in ("integer", "real")
        or parts[4] != "general"
    ):
        _fail(
            "matrix.mtx 格式头不支持：仅接受 coordinate integer general 或 "
            f"coordinate real general，得到 {banner!r}"
        )
    kind = parts[3]

    # 尺寸行之前只允许注释行（% 开头）与空行
    index = 1
    while index < len(lines) and (
        lines[index].startswith("%") or lines[index].strip() == ""
    ):
        index += 1
    if index >= len(lines):
        _fail("matrix.mtx 缺少尺寸行（rows columns entries）")

    size_fields = lines[index].split()
    if len(size_fields) != 3:
        _fail(
            f"matrix.mtx 尺寸行须含三个整数（rows columns entries），"
            f"得到 {lines[index]!r}"
        )
    n_rows = _parse_nonnegative_int(size_fields[0], "matrix.mtx 尺寸行行数")
    n_cols = _parse_nonnegative_int(size_fields[1], "matrix.mtx 尺寸行列数")
    n_entries = _parse_nonnegative_int(size_fields[2], "matrix.mtx 尺寸行条目数")

    if n_rows != n_genes:
        _fail(
            f"matrix.mtx 声明行数 {n_rows} 与 features.tsv 基因数 "
            f"{n_genes} 不一致"
        )
    if n_cols != n_cells:
        _fail(
            f"matrix.mtx 声明列数 {n_cols} 与 barcodes.tsv 细胞数 "
            f"{n_cells} 不一致"
        )

    counts = [[0] * n_cells for _ in range(n_rows)]
    seen_coords = set()
    data_lines = lines[index + 1:]
    n_seen = 0
    for offset, line in enumerate(data_lines, start=index + 2):
        fields = line.split()
        if len(fields) != 3:
            _fail(
                f"matrix.mtx 第 {offset} 行须含 “行 列 值” 三个字段，"
                f"得到 {line!r}"
            )
        row_idx = _parse_nonnegative_int(fields[0], f"matrix.mtx 第 {offset} 行行索引")
        col_idx = _parse_nonnegative_int(fields[1], f"matrix.mtx 第 {offset} 行列索引")
        if not 1 <= row_idx <= n_rows:
            _fail(
                f"matrix.mtx 第 {offset} 行行索引 {row_idx} 越界"
                f"（合法范围 1..{n_rows}）"
            )
        if not 1 <= col_idx <= n_cols:
            _fail(
                f"matrix.mtx 第 {offset} 行列索引 {col_idx} 越界"
                f"（合法范围 1..{n_cols}）"
            )
        coord = (row_idx, col_idx)
        if coord in seen_coords:
            _fail(
                f"matrix.mtx 第 {offset} 行坐标 ({row_idx}, {col_idx}) 重复"
            )
        seen_coords.add(coord)
        counts[row_idx - 1][col_idx - 1] = _parse_mtx_value(
            fields[2], kind, offset
        )
        n_seen += 1

    if n_seen != n_entries:
        _fail(
            f"matrix.mtx 尺寸行声明 {n_entries} 个条目，实际有 {n_seen} 个"
        )

    return counts


def read_mtx_directory(path: str) -> ExpressionMatrix:
    """读取 10x MatrixMarket 目录并展开为稠密 :class:`ExpressionMatrix`。"""
    if path is None or path == "":
        _fail("输入路径为空")
    if not os.path.exists(path):
        _fail(f"mtx 输入目录不存在：{path}")
    if not os.path.isdir(path):
        _fail(f"mtx 输入路径不是目录：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"mtx 输入目录不可读：{path}")

    # 先确认三个文件齐备，再逐个读取并计算 SHA-256
    matrix_raw, matrix_text = _require_input_file(path, MTX_MATRIX_FILE)
    barcodes_raw, barcodes_text = _require_input_file(path, MTX_BARCODES_FILE)
    features_raw, features_text = _require_input_file(path, MTX_FEATURES_FILE)

    cell_ids = _read_id_lines(
        barcodes_text, MTX_BARCODES_FILE, first_column_only=False
    )
    gene_ids = _read_id_lines(
        features_text, MTX_FEATURES_FILE, first_column_only=True
    )
    if not cell_ids:
        _fail(f"{MTX_BARCODES_FILE} 没有任何非空细胞 ID 行")
    if not gene_ids:
        _fail(f"{MTX_FEATURES_FILE} 没有任何非空基因 ID 行")

    counts = _parse_mtx_matrix(
        matrix_text, n_genes=len(gene_ids), n_cells=len(cell_ids)
    )
    total = sum(sum(row) for row in counts)

    source_files = (
        SourceFile(MTX_MATRIX_FILE, hashlib.sha256(matrix_raw).hexdigest()),
        SourceFile(MTX_BARCODES_FILE, hashlib.sha256(barcodes_raw).hexdigest()),
        SourceFile(MTX_FEATURES_FILE, hashlib.sha256(features_raw).hexdigest()),
    )
    return ExpressionMatrix(
        gene_ids=gene_ids,
        cell_ids=cell_ids,
        counts=counts,
        sha256=source_files[0].sha256,
        path=path,
        total_counts=total,
        input_format=MTX_FORMAT,
        source_files=source_files,
    )


def load_matrix(path: str, input_format: str) -> ExpressionMatrix:
    """按指定格式分派读取（``tsv`` 或 ``mtx``）。"""
    if input_format == TSV_FORMAT:
        return read_matrix(path)
    if input_format == MTX_FORMAT:
        return read_mtx_directory(path)
    # 正常情况下命令行与 validate_config 已拦截；程序化调用绕过配置层时
    # 仍按配置错误（而非输入错误）报告
    raise CellFlowConfigError(
        f"未知输入格式：{input_format!r}（仅支持 {TSV_FORMAT!r} 或 {MTX_FORMAT!r}）"
    )


def _first_duplicate(values: List[str]) -> str:
    seen = set()
    for value in values:
        if value in seen:
            return value
        seen.add(value)
    return values[0]
