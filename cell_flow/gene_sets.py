"""基因集评分：读入与计算。

``--gene-sets`` 指向 UTF-8 制表符文本，表头恰为 ``set_id``、``gene_id``
两列，每行一个成员关系；字段非空、(set_id, gene_id) 组合唯一，至少一条
数据行。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名无关）。
基因 ID 与表达矩阵首列精确匹配，不做大小写、别名、前缀转换；矩阵外基因
不计分但计入总数。任何不合法（空路径、UTF-8、表头、空值、重复、数据行
不足，或 gzip 多成员、尾随数据、截断、CRC/长度错误）都抛
:class:`cell_flow.errors.CellFlowInputError`（退出码 2），且不改动结果。

评分只用质控后保留细胞与基因：有 ``--batch-metadata`` 时取批次均值中心化
值，否则取 log 归一化值。score 为集合与保留基因交集内表达值的算术平均；
任一集合与保留基因交集为空都抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowDataError, CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single


@dataclass(frozen=True)
class GeneSets:
    """基因集成员关系，按文件行序保留集合出现顺序与每个集合的成员。"""

    path: str
    name: str                 # 文件名（路径末段）
    sha256: str               # 原始字节（gzip 即压缩字节）的 SHA-256
    members: Dict[str, List[str]]  # set_id -> 成员基因 ID（文件行序）


@dataclass(frozen=True)
class GeneSetScore:
    """一个基因集对一个保留细胞的评分结果。"""

    set_id: str
    n_total: int
    n_used: int
    cell_id: str
    score: float


@dataclass(frozen=True)
class GeneSetScores:
    """全部基因集在全部保留细胞上的评分。

    ``rows`` 按 set_id 升序、集合内按保留细胞原顺序排列；
    ``set_order`` 为升序 set_id 列表，供图表列序使用。
    """

    rows: List[GeneSetScore]
    set_order: List[str]
    totals: Dict[str, int]    # set_id -> n_genes_total
    used: Dict[str, int]      # set_id -> n_genes_used


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_gene_sets(path: str) -> GeneSets:
    """读取并严格校验基因集文件（不依赖表达矩阵，可在触碰输出目录前失败）。"""
    if path is None or path == "":
        _fail("基因集文件路径为空")
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

    # splitlines 兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("基因集表头缺失：文件为空")
    header = lines[0].split("\t")
    if header != ["set_id", "gene_id"]:
        _fail("基因集表头必须恰为 set_id 和 gene_id 两列")

    members: Dict[str, List[str]] = {}
    seen_pairs = set()
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
            _fail(f"基因集第 {offset} 行 set_id 为空")
        if gene_id == "":
            _fail(f"基因集第 {offset} 行 gene_id 为空")
        pair = (set_id, gene_id)
        if pair in seen_pairs:
            _fail(
                f"基因集成员关系重复：set_id={set_id!r}、gene_id={gene_id!r}"
            )
        seen_pairs.add(pair)
        members.setdefault(set_id, []).append(gene_id)

    if not members:
        _fail("基因集没有任何数据行")

    return GeneSets(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        members=members,
    )


def score_gene_sets(
    gene_sets: GeneSets,
    *,
    kept_gene_ids: List[str],
    cell_ids: List[str],
    analysis_values: List[List[float]],
) -> GeneSetScores:
    """计算各基因集对每个保留细胞的算术平均评分。

    ``analysis_values[g][c]`` 与 ``kept_gene_ids``/``cell_ids`` 对齐，
    为批次中心化值（有批次元数据）或 log 归一化值。矩阵外基因不计分但
    计入 ``n_genes_total``；任一集合与保留基因交集为空即无法评分，
    抛数据错误（退出码 4），不产出任何结果。
    """
    gene_position: Dict[str, int] = {
        gene_id: g for g, gene_id in enumerate(kept_gene_ids)
    }

    totals: Dict[str, int] = {}
    used: Dict[str, int] = {}
    set_positions: Dict[str, List[int]] = {}
    for set_id in sorted(gene_sets.members):
        entries = gene_sets.members[set_id]
        # 集合内成员在文件中已保证 (set_id, gene_id) 唯一，直接建下标映射
        positions = [
            gene_position[gene_id]
            for gene_id in entries
            if gene_id in gene_position
        ]
        totals[set_id] = len(entries)
        used[set_id] = len(positions)
        if not positions:
            raise CellFlowDataError(
                f"基因集 {set_id!r} 与质控后保留基因的交集为空，无法评分"
            )
        set_positions[set_id] = positions

    rows: List[GeneSetScore] = []
    for set_id in sorted(gene_sets.members):
        positions = set_positions[set_id]
        n_used = len(positions)
        n_total = len(gene_sets.members[set_id])
        for c, cell_id in enumerate(cell_ids):
            # 沿用管线既有口径（内置 sum）累加，保证与归一化/批次校正同精度
            total = sum(analysis_values[g][c] for g in positions)
            rows.append(
                GeneSetScore(
                    set_id=set_id,
                    n_total=n_total,
                    n_used=n_used,
                    cell_id=cell_id,
                    score=total / n_used,
                )
            )

    return GeneSetScores(
        rows=rows,
        set_order=sorted(gene_sets.members),
        totals=totals,
        used=used,
    )
