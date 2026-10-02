"""输入矩阵读取与严格校验。

输入为制表符分隔矩阵：首列是唯一基因 ID，其余列名是唯一细胞 ID，
取值为非负整数 UMI 计数。任何不合法情形都抛出
:class:`cell_flow.errors.CellFlowInputError`。
"""

import hashlib
import os
import re
from dataclasses import dataclass, field
from typing import List

from .errors import CellFlowInputError

_INTEGER_RE = re.compile(r"^(0|[1-9][0-9]*)$")


@dataclass(frozen=True)
class InputFile:
    """一个输入文件的来源记录：文件名与原始字节的 SHA-256。"""

    name: str
    sha256: str


@dataclass(frozen=True)
class ExpressionMatrix:
    """原始计数矩阵（行＝基因，列＝细胞）。

    ``input_format`` 为 ``"tsv"`` 时沿用 ``path``/``sha256`` 单一文件记录；
    为 ``"mtx"`` 时 ``files`` 依次记录 matrix.mtx、barcodes.tsv、features.tsv。
    """

    gene_ids: List[str]
    cell_ids: List[str]
    counts: List[List[int]]  # counts[gene_index][cell_index]
    sha256: str
    path: str
    total_counts: int
    input_format: str = "tsv"
    files: List[InputFile] = field(default_factory=list)

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
    if path is None or path == "":
        _fail("输入路径为空")
    if not os.path.exists(path):
        _fail(f"输入文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"输入路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"输入文件不可读：{path}")

    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"输入文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    try:
        text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"输入文件不是合法的 UTF-8 文本：{path}")

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

    total = 0
    for row in counts:
        total += sum(row)

    return ExpressionMatrix(
        gene_ids=gene_ids,
        cell_ids=cell_ids,
        counts=counts,
        sha256=digest,
        path=path,
        total_counts=total,
    )


def _first_duplicate(values: List[str]) -> str:
    seen = set()
    for value in values:
        if value in seen:
            return value
        seen.add(value)
    return values[0]
