"""跨样本批次校正的细胞元数据读取、严格校验与批次混合分数。

``--cell-metadata`` 指向 UTF-8 制表符文本：首列必须是细胞条码列
``cell_id``，并须包含样本标识列（默认 ``sample_id``，可用
``--sample-column`` 指定）与批次标签列（默认 ``batch``，可用
``--batch-column`` 指定）；其余列任意，原样保留并透传到公开细胞结果，
不覆盖、不重命名任何既有字段。``cell_id`` 唯一且与表达矩阵的全部细胞
一一对应（不多不少），样本标识与批次标签均非空（空批次不静默归入
未知批次）。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名
无关）。任何不合法（表头缺列、列名重复、条码重复、细胞集合不一致、
空样本、空批次，或 gzip 多成员、尾随数据、截断、校验失败）都抛
:class:`cell_flow.errors.CellFlowInputError`（即 ``ValueError``），
且在启动降维聚类之前失败，不生成任何部分成功结果。

校正发生在质量控制之后、降维聚类之前：质控与 log 归一化沿用既有规则，
随后按基因做批次均值中心化（与 :mod:`cell_flow.batch` 同一口径），
校正值进入高变基因选择、PCA、聚类、markers 与差异表达。质控后仅一个
批次时不施加任何扰动，下游结果与无批次运行逐字节一致。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single

BARCODE_COLUMN = "cell_id"
DEFAULT_BATCH_COLUMN = "batch"
DEFAULT_SAMPLE_COLUMN = "sample_id"
MIXING_NEIGHBORS = 15


@dataclass(frozen=True)
class CellMetadata:
    """细胞元数据：每个矩阵细胞恰好一行，样本标识与批次标签均非空。"""

    path: str
    name: str                    # 文件名（路径末段）
    sha256: str                  # 原始字节（gzip 即压缩字节）的 SHA-256
    columns: List[str]           # 表头全部列名（原顺序，含 cell_id）
    rows: Dict[str, List[str]]   # cell_id -> 与 columns 对齐的整行字段
    batch_column: str
    sample_column: str
    batches: Dict[str, str]      # cell_id -> batch
    samples: Dict[str, str]      # cell_id -> sample
    batch_sizes: Dict[str, int]  # batch -> 细胞数（按 batch 升序）


@dataclass(frozen=True)
class CellBatchReport:
    """一次批次校正的公开产物数据（均与保留细胞顺序对齐）。"""

    columns: List[str]           # 元数据表头（原顺序）
    rows: List[List[str]]        # 保留细胞的元数据整行（保留细胞顺序）
    batches: List[str]           # 保留细胞的原始批次标签
    mixing_before: float         # 校正前批次混合分数
    mixing_after: float          # 校正后批次混合分数


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_cell_metadata(
    path: str,
    cell_ids: List[str],
    *,
    batch_column: str = DEFAULT_BATCH_COLUMN,
    sample_column: str = DEFAULT_SAMPLE_COLUMN,
) -> CellMetadata:
    """读取并校验细胞元数据；``cell_ids`` 为表达矩阵的全部细胞 ID（列序）。"""
    if path is None or path == "":
        _fail("细胞元数据路径为空")
    if not os.path.exists(path):
        _fail(f"细胞元数据文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"细胞元数据路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"细胞元数据文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"细胞元数据文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"细胞元数据文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("细胞元数据表头缺失：文件为空")
    columns = lines[0].split("\t")
    if any(column == "" for column in columns):
        _fail("细胞元数据表头存在空列名")
    if len(set(columns)) != len(columns):
        seen = set()
        dup = next(c for c in columns if c in seen or seen.add(c))
        _fail(f"细胞元数据表头列名重复：{dup!r}")
    if BARCODE_COLUMN not in columns:
        _fail(f"细胞元数据表头缺少细胞条码列 {BARCODE_COLUMN}")
    if sample_column not in columns:
        _fail(f"细胞元数据表头缺少样本标识列：{sample_column!r}")
    if batch_column not in columns:
        _fail(f"细胞元数据表头缺少批次标签列：{batch_column!r}")

    barcode_idx = columns.index(BARCODE_COLUMN)
    sample_idx = columns.index(sample_column)
    batch_idx = columns.index(batch_column)
    n_columns = len(columns)

    rows: Dict[str, List[str]] = {}
    batches: Dict[str, str] = {}
    samples: Dict[str, str] = {}
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"细胞元数据第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != n_columns:
            _fail(
                f"细胞元数据第 {offset} 行列数为 {len(fields)}，"
                f"与表头 {n_columns} 列不一致"
            )
        cell_id = fields[barcode_idx]
        if cell_id == "":
            _fail(f"细胞元数据第 {offset} 行细胞条码为空")
        if cell_id in rows:
            _fail(f"细胞元数据细胞条码重复：{cell_id!r}")
        sample = fields[sample_idx]
        if sample == "":
            _fail(f"细胞元数据第 {offset} 行样本标识为空（细胞 {cell_id!r}）")
        batch = fields[batch_idx]
        if batch == "":
            # 空批次不静默归入未知批次，直接报输入冲突
            _fail(f"细胞元数据第 {offset} 行批次标签为空（细胞 {cell_id!r}）")
        rows[cell_id] = fields
        samples[cell_id] = sample
        batches[cell_id] = batch

    missing = [cell_id for cell_id in cell_ids if cell_id not in rows]
    if missing:
        _fail(
            f"细胞元数据与矩阵细胞集合不一致：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个细胞"
        )
    matrix_cells = set(cell_ids)
    extra = [cell_id for cell_id in rows if cell_id not in matrix_cells]
    if extra:
        _fail(
            f"细胞元数据与矩阵细胞集合不一致：包含矩阵之外的细胞 "
            f"{extra[0]!r} 等 {len(extra)} 个"
        )

    sizes: Dict[str, int] = {}
    for batch in batches.values():
        sizes[batch] = sizes.get(batch, 0) + 1
    batch_sizes = {batch: sizes[batch] for batch in sorted(sizes)}

    return CellMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        columns=columns,
        rows=rows,
        batch_column=batch_column,
        sample_column=sample_column,
        batches=batches,
        samples=samples,
        batch_sizes=batch_sizes,
    )


def batch_mixing_score(
    scores: List[List[float]],
    batch_labels: List[str],
    *,
    n_neighbors: int = MIXING_NEIGHBORS,
) -> float:
    """批次混合分数：每个细胞最近邻中来自其他批次的占比，对全部细胞求平均。

    ``scores[c]`` 为细胞 ``c`` 的低维坐标（PCA 得分），``batch_labels[c]``
    为其批次。最近邻按欧氏距离（比较时用平方距离，次序等价）选取，
    距离并列时按细胞下标升序，结果确定；邻居数取 ``n_neighbors`` 与
    细胞数减一的较小值。全部细胞同属一个批次时分数为 0。
    """
    n_cells = len(scores)
    if n_cells < 2:
        return 0.0
    k = min(n_neighbors, n_cells - 1)
    total = 0.0
    for i in range(n_cells):
        own = scores[i]
        distances = []
        for j in range(n_cells):
            if j == i:
                continue
            other = scores[j]
            dist_sq = 0.0
            for a, b in zip(own, other):
                diff = a - b
                dist_sq += diff * diff
            distances.append((dist_sq, j))
        # 距离并列按下标升序，保证同输入同结果
        distances.sort(key=lambda item: (item[0], item[1]))
        own_batch = batch_labels[i]
        foreign = sum(
            1 for _, j in distances[:k] if batch_labels[j] != own_batch
        )
        total += foreign / k
    return total / n_cells
