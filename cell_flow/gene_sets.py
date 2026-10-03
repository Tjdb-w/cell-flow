"""基因集读取、严格校验与逐细胞集合评分。

``--gene-sets`` 指向 UTF-8 制表符文本，表头恰为 ``set_id`` 与 ``gene_id``
两列；每行一个成员关系，字段非空、``(set_id, gene_id)`` 组合唯一，且至少
有一条数据行。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名
无关）。任何不合法（UTF-8 解码失败、表头不符、空字段、组合重复、无数据行，
或 gzip 多成员、尾随数据、截断、CRC/长度校验失败）都抛
:class:`cell_flow.errors.CellFlowInputError`，且不改动任何已有结果。

评分只用质控后保留的细胞与基因：有 ``--batch-metadata`` 时取批次均值
中心化值，否则取 log 归一化值（即 ``NormalizedData.analysis_values``）。
基因 ID 与表达矩阵首列精确匹配，不做大小写、别名或前缀转换；矩阵之外
（含被质控滤除）的基因不计分但计入 ``n_genes_total``。score 为集合与
保留基因交集内表达值的算术平均；任一集合交集为空抛
:class:`cell_flow.errors.CellFlowDataError`。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .errors import CellFlowDataError, CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single
from .normalize import NormalizedData


@dataclass(frozen=True)
class GeneSet:
    """一个基因集：ID 与按文件行序排列的成员基因（集合内唯一）。"""

    set_id: str
    gene_ids: Tuple[str, ...]


@dataclass(frozen=True)
class GeneSets:
    """基因集来源记录与全部集合（按 set_id 升序）。"""

    path: str
    name: str                  # 文件名（路径末段）
    sha256: str                # 原始字节（gzip 即压缩字节）的 SHA-256
    sets: Tuple[GeneSet, ...]


@dataclass(frozen=True)
class GeneSetScore:
    """一个集合在全部保留细胞上的评分（与保留细胞列序对齐）。"""

    set_id: str
    n_genes_total: int         # 集合列出的全部基因数（含矩阵外基因）
    n_genes_used: int          # 与保留基因的交集大小
    scores: Tuple[float, ...]


@dataclass(frozen=True)
class GeneSetScores:
    """全部集合评分：集合按 set_id 升序，细胞按保留细胞原顺序。"""

    cell_ids: List[str]
    sets: List[GeneSetScore]


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_gene_sets(path: str) -> GeneSets:
    """读取并校验基因集文件。"""
    if path is None or path == "":
        _fail("基因集路径为空")
    if not os.path.exists(path):
        _fail(f"基因集文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"基因集路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"基因集文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"基因集文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"基因集文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("基因集表头缺失：文件为空")
    header = lines[0].split("\t")
    if header != ["set_id", "gene_id"]:
        _fail("基因集表头必须恰为 set_id 和 gene_id 两列")

    members: Dict[str, List[str]] = {}
    seen_pairs = set()
    n_rows = 0
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"基因集第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != 2:
            _fail(
                f"基因集第 {offset} 行列数为 {len(fields)}，与表头两列不一致"
            )
        set_id, gene_id = fields
        if set_id == "":
            _fail(f"基因集第 {offset} 行集合 ID 为空")
        if gene_id == "":
            _fail(f"基因集第 {offset} 行基因 ID 为空（集合 {set_id!r}）")
        pair = (set_id, gene_id)
        if pair in seen_pairs:
            _fail(f"基因集成员关系重复：{set_id!r} × {gene_id!r}")
        seen_pairs.add(pair)
        members.setdefault(set_id, []).append(gene_id)
        n_rows += 1

    if n_rows == 0:
        _fail("基因集没有任何成员数据行")

    sets = tuple(
        GeneSet(set_id=set_id, gene_ids=tuple(members[set_id]))
        for set_id in sorted(members)
    )
    return GeneSets(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        sets=sets,
    )


def score_gene_sets(gene_sets: GeneSets, data: NormalizedData) -> GeneSetScores:
    """在质控后保留细胞与基因上为每个集合打分。

    表达矩阵取 ``data.analysis_values``（有批次校正时为批次均值中心化值，
    否则为 log 归一化值）；score 为集合与保留基因交集内表达值的算术平均。
    任一集合与保留基因交集为空时抛 :class:`CellFlowDataError`。
    """
    values = data.analysis_values
    gene_row = {gene_id: g for g, gene_id in enumerate(data.gene_ids)}
    n_cells = len(data.cell_ids)

    scored: List[GeneSetScore] = []
    for gene_set in gene_sets.sets:
        rows = [
            gene_row[gene_id]
            for gene_id in gene_set.gene_ids
            if gene_id in gene_row
        ]
        if not rows:
            raise CellFlowDataError(
                f"基因集 {gene_set.set_id!r} 与质控后保留基因交集为空，"
                f"无法评分"
            )
        n_used = len(rows)
        scores: List[float] = []
        for c in range(n_cells):
            total = 0.0
            for g in rows:
                total += values[g][c]
            scores.append(total / n_used)
        scored.append(
            GeneSetScore(
                set_id=gene_set.set_id,
                n_genes_total=len(gene_set.gene_ids),
                n_genes_used=n_used,
                scores=tuple(scores),
            )
        )
    return GeneSetScores(cell_ids=list(data.cell_ids), sets=scored)
