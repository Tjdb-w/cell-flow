"""批次元数据读取、按批次的基因均值中心化校正、批次混合分数与批次汇总。

``--batch-metadata`` 指向 UTF-8 制表符文本，必须包含三列：``cell_id``
（细胞条码）、``sample_id``（样本标识）与 ``batch``（批次标签）；列名可
通过 ``--sample-column`` / ``--batch-column`` 改用其他名字，``cell_id``
列名固定。表头中其余列原样保留并透传到逐细胞公开结果。文件可为纯文本或
单成员 gzip（按 gzip 魔数识别，与文件名无关）。

校验规则（任一不满足都抛 :class:`ValueError`，消息指出具体冲突类型，
且不改动任何已有结果）：

- 表头必须同时含条码列、样本列与批次列；
- ``cell_id`` 与样本/批次字段均非空，条码不得重复；
- 细胞集合必须与表达矩阵全部细胞一一对应：既不缺少也不多余。

单批次数据不报错、不校正，按原流程继续分析；多批次数据在质控、log
归一化之后做按基因的批次均值中心化：每个归一化值减去对应批次保留细胞
的基因均值，再加回全部保留细胞的基因总均值（batch mean centering）。
校正值进入高变基因选择、PCA、聚类、markers、成对 markers、分组差异表达
与图表数据；``normalized_expression.tsv`` 始终写未校正值。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .errors import CellFlowBatchError
from .io import GZIP_MAGIC, gunzip_single

# 批次混合分数所用的最近邻个数
MIXING_N_NEIGHBORS = 15


@dataclass(frozen=True)
class BatchMetadata:
    """细胞批次元数据：每个矩阵细胞恰好归属一个非空样本与非空批次。"""

    path: str
    name: str                    # 文件名（路径末段）
    sha256: str                  # 原始字节（gzip 即压缩字节）的 SHA-256
    sample_column: str           # 样本列实际列名
    batch_column: str            # 批次列实际列名
    samples: Dict[str, str]      # cell_id -> sample_id
    batches: Dict[str, str]      # cell_id -> batch
    batch_sizes: Dict[str, int]  # batch -> 全部输入细胞数（按 batch 升序）
    extra_columns: Tuple[str, ...]              # 透传列列名（表头原顺序）
    extra_values: Dict[str, Dict[str, str]]     # 列名 -> (cell_id -> 值)


@dataclass(frozen=True)
class BatchSummaryRow:
    """一个批次在质控后保留细胞上的汇总（均值按批次内保留细胞计算）。"""

    batch_id: str
    n_cells: int
    mean_total_counts: float
    mean_detected_genes: float
    mean_mitochondrial_fraction: float


@dataclass(frozen=True)
class BatchMixing:
    """校正前后的批次混合分数（越大表示批次在最近邻中混合越充分）。"""

    before_score: float   # 未校正低维表示上的平均混合比例
    after_score: float    # 校正后实际用于聚类的低维表示上的平均混合比例
    applied: bool         # 本次是否实际执行了批次中心化
    n_neighbors: int      # 每个细胞实际使用的最近邻个数（细胞不足时为 n-1）


def _fail(message: str) -> None:
    raise CellFlowBatchError(message)


def read_batch_metadata(
    path: str,
    cell_ids: Sequence[str],
    *,
    sample_column: str = "sample_id",
    batch_column: str = "batch",
) -> BatchMetadata:
    """读取并校验批次元数据；``cell_ids`` 为表达矩阵的全部细胞 ID（列序）。

    列名冲突（条码列缺失、样本/批次列缺失或重名等）、条码重复、字段为空、
    细胞集合与矩阵不一致，均抛 :class:`ValueError`。
    """
    if path is None or path == "":
        _fail("批次元数据路径为空")
    if not isinstance(sample_column, str) or sample_column == "":
        _fail("样本列名必须是非空字符串")
    if not isinstance(batch_column, str) or batch_column == "":
        _fail("批次列名必须是非空字符串")
    if sample_column == "cell_id":
        _fail("样本列名不能与固定条码列名 cell_id 相同")
    if batch_column == "cell_id":
        _fail("批次列名不能与固定条码列名 cell_id 相同")
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
    duplicate_columns = {name for name in header if header.count(name) > 1}
    if duplicate_columns:
        _fail(
            f"批次元数据表头列名重复：{sorted(duplicate_columns)[0]!r}，"
            f"无法无歧义地读取列"
        )

    def _require_column(role: str, column: str) -> int:
        matches = [i for i, name in enumerate(header) if name == column]
        if not matches:
            _fail(
                f"批次元数据缺少{role}列 {column!r}：表头为 {header}，"
                f"无法关联细胞批次信息"
            )
        if len(matches) > 1:
            _fail(f"批次元数据{role}列 {column!r} 在表头中重复出现")
        return matches[0]

    id_pos = _require_column("细胞条码", "cell_id")
    sample_pos = _require_column("样本标识", sample_column)
    batch_pos = _require_column("批次标签", batch_column)
    if sample_column == batch_column:
        _fail(
            f"样本列与批次列不能同名：{sample_column!r} 同时被指定为"
            f" --sample-column 与 --batch-column"
        )

    reserved = {id_pos, sample_pos, batch_pos}
    extra_columns = tuple(
        name for i, name in enumerate(header) if i not in reserved
    )
    extra_values: Dict[str, Dict[str, str]] = {
        name: {} for name in extra_columns
    }

    samples: Dict[str, str] = {}
    batches: Dict[str, str] = {}
    n_columns = len(header)
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"批次元数据第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != n_columns:
            _fail(
                f"批次元数据第 {offset} 行列数为 {len(fields)}，"
                f"与表头列数 {n_columns} 不一致"
            )
        cell_id = fields[id_pos]
        sample = fields[sample_pos]
        batch = fields[batch_pos]
        if cell_id == "":
            _fail(f"批次元数据第 {offset} 行细胞条码为空")
        if sample == "":
            _fail(
                f"批次元数据第 {offset} 行样本标识（列 {sample_column!r}）"
                f"为空（细胞 {cell_id!r}）"
            )
        if batch == "":
            _fail(
                f"批次元数据第 {offset} 行批次标签缺失或为空"
                f"（列 {batch_column!r}，细胞 {cell_id!r}），"
                f"不能归入未知批次"
            )
        if cell_id in batches:
            _fail(f"批次元数据细胞条码重复：{cell_id!r}")
        samples[cell_id] = sample
        batches[cell_id] = batch
        # 其余列原样保留：按表头下标取值，值随细胞透传，不做任何转换
        for pos, name in enumerate(header):
            if pos in reserved:
                continue
            extra_values[name][cell_id] = fields[pos]

    missing = [cell_id for cell_id in cell_ids if cell_id not in batches]
    if missing:
        _fail(
            f"批次元数据与表达矩阵的细胞集合不一致：元数据缺少矩阵细胞 "
            f"{missing[0]!r} 等 {len(missing)} 个"
        )
    matrix_cells = set(cell_ids)
    extra = [cell_id for cell_id in batches if cell_id not in matrix_cells]
    if extra:
        _fail(
            f"批次元数据与表达矩阵的细胞集合不一致：元数据包含矩阵之外的"
            f"细胞 {extra[0]!r} 等 {len(extra)} 个"
        )

    sizes: Dict[str, int] = {}
    for batch in batches.values():
        sizes[batch] = sizes.get(batch, 0) + 1
    batch_sizes = {batch: sizes[batch] for batch in sorted(sizes)}

    return BatchMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        sample_column=sample_column,
        batch_column=batch_column,
        samples=samples,
        batches=batches,
        batch_sizes=batch_sizes,
        extra_columns=extra_columns,
        extra_values=extra_values,
    )


def center_by_batch(
    values: List[List[float]], batch_labels: List[str]
) -> List[List[float]]:
    """按基因做批次均值中心化。

    ``values[g][c]`` 为保留基因 × 保留细胞的 log 归一化表达，
    ``batch_labels[c]`` 为对应细胞的批次。校正值 =
    原值 - 该细胞所在批次的基因均值 + 全部保留细胞的基因总均值。
    均值按列序累加，结果确定。单批次输入不应调用本函数（调用即原样
    不校正的快速路径在管线中处理）。
    """
    n_cells = len(batch_labels)
    corrected: List[List[float]] = []
    for row in values:
        grand_mean = sum(row) / n_cells
        batch_sums: Dict[str, float] = {}
        batch_counts: Dict[str, int] = {}
        for c, batch in enumerate(batch_labels):
            batch_sums[batch] = batch_sums.get(batch, 0.0) + row[c]
            batch_counts[batch] = batch_counts.get(batch, 0) + 1
        batch_means = {
            batch: batch_sums[batch] / batch_counts[batch]
            for batch in batch_sums
        }
        corrected.append(
            [
                row[c] - batch_means[batch_labels[c]] + grand_mean
                for c in range(n_cells)
            ]
        )
    return corrected


def _squared_distance(a: Sequence[float], b: Sequence[float]) -> float:
    total = 0.0
    for x, y in zip(a, b):
        d = x - y
        total += d * d
    return total


def batch_mixing_score(
    points: Sequence[Sequence[float]],
    batch_labels: Sequence[str],
    n_neighbors: int = MIXING_N_NEIGHBORS,
) -> float:
    """批次混合分数：每个细胞取最近的 ``n_neighbors`` 个其他细胞，
    统计其中批次标签与该细胞不同者的占比，再对全部细胞求平均。

    最近邻按欧氏平方距离升序选取，距离并列时以细胞下标升序打破，
    因此结果与输入顺序一一对应、完全确定，不使用随机数。点数不足
    ``n_neighbors + 1`` 时取全部其他细胞。
    """
    n = len(points)
    if n < 2:
        return 0.0
    k = min(n_neighbors, n - 1)
    total_fraction = 0.0
    for c in range(n):
        own = points[c]
        # (平方距离, 下标) 升序：距离并列取较小下标，跳过自身
        ranked = sorted(
            (j for j in range(n) if j != c),
            key=lambda j: (_squared_distance(own, points[j]), j),
        )
        neighbors = ranked[:k]
        different = sum(1 for j in neighbors if batch_labels[j] != batch_labels[c])
        total_fraction += different / k
    return total_fraction / n
