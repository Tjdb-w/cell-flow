"""生物学重复（样本）元数据读取，以及按样本汇总的 pseudobulk 差异表达。

``--replicate-metadata`` 指向 UTF-8 制表符文本，表头恰为 ``cell_id``、
``sample_id``、``group`` 三列；每个输入细胞恰好一行，``cell_id`` 唯一且与
表达矩阵的全部细胞一一对应（不多不少），``sample_id`` 与 ``group`` 非空，
同一 ``sample_id`` 只能归属一个 ``group``。文件可为纯文本或单成员 gzip
（按 gzip 魔数识别，与文件名无关）。任何不合法（表头不符、ID 重复、
覆盖不全、多出矩阵之外的细胞、空样本或空分组、同一样本归属多个分组，
或 gzip 多成员、尾随数据、截断、校验失败）都抛
:class:`cell_flow.errors.CellFlowInputError`，且不改动任何已有结果。

分析只用既有质控保留的细胞与基因：按 ``sample_id`` 对原始计数求和形成
pseudobulk，样本文库归一到 10000 后取 log1p，再以样本为观测单位对每个
group 做 one-vs-rest、按 group 升序两两比较（Welch t 检验、比较内 BH）。
质控后样本无细胞、某 group 有效重复少于两个或不足两个 group 都抛
:class:`cell_flow.errors.CellFlowDataError`。
"""

import hashlib
import math
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowDataError, CellFlowInputError
from .io import ExpressionMatrix, GZIP_MAGIC, gunzip_single
from .markers import GroupComparison, run_group_comparisons
from .qc import QCResult

NORMALIZATION_TARGET = 10000.0


@dataclass(frozen=True)
class ReplicateMetadata:
    """细胞 -> 样本 -> 分组 的归属关系；每个矩阵细胞恰好一行。"""

    path: str
    name: str                      # 文件名（路径末段）
    sha256: str                    # 原始字节（gzip 即压缩字节）的 SHA-256
    samples: Dict[str, str]        # cell_id -> sample_id
    groups: Dict[str, str]         # cell_id -> group
    sample_group: Dict[str, str]   # sample_id -> group（同一样本仅一个分组）
    sample_order: List[str]        # 样本按输入细胞首次出现顺序
    group_sizes: Dict[str, int]    # group -> 细胞数（按 group 升序，质控前）
    sample_sizes: Dict[str, int]   # sample_id -> 细胞数（质控前）


@dataclass(frozen=True)
class PseudobulkData:
    """质控后保留细胞与基因汇总出的样本级 pseudobulk。"""

    gene_ids: List[str]            # 保留基因（原矩阵行序）
    sample_ids: List[str]          # 有保留细胞的样本（首次出现顺序）
    sample_groups: List[str]       # 与 sample_ids 对齐的 group
    counts: List[List[int]]        # counts[gene][sample]：保留基因上的原始计数和
    sample_totals: List[int]       # 与 sample_ids 对齐的样本文库总计数
    values: List[List[float]]      # values[gene][sample] = ln(count/total*1e4 + 1)
    comparisons: List[GroupComparison]


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_replicate_metadata(
    path: str, cell_ids: List[str]
) -> ReplicateMetadata:
    """读取并校验重复元数据；``cell_ids`` 为表达矩阵的全部细胞 ID（列序）。"""
    if path is None or path == "":
        _fail("重复元数据路径为空")
    if not os.path.exists(path):
        _fail(f"重复元数据文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"重复元数据路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"重复元数据文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"重复元数据文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"重复元数据文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("重复元数据表头缺失：文件为空")
    header = lines[0].split("\t")
    if header != ["cell_id", "sample_id", "group"]:
        _fail(
            "重复元数据表头必须恰为 cell_id、sample_id、group 三列"
        )

    samples: Dict[str, str] = {}
    groups: Dict[str, str] = {}
    sample_group: Dict[str, str] = {}
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"重复元数据第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != 3:
            _fail(
                f"重复元数据第 {offset} 行列数为 {len(fields)}，"
                f"与表头三列不一致"
            )
        cell_id, sample_id, group = fields
        if cell_id == "":
            _fail(f"重复元数据第 {offset} 行细胞 ID 为空")
        if sample_id == "":
            _fail(f"重复元数据第 {offset} 行样本 ID 为空（细胞 {cell_id!r}）")
        if group == "":
            _fail(f"重复元数据第 {offset} 行分组为空（细胞 {cell_id!r}）")
        if cell_id in samples:
            _fail(f"重复元数据细胞 ID 重复：{cell_id!r}")
        if sample_id in sample_group:
            previous = sample_group[sample_id]
            if previous != group:
                _fail(
                    f"重复元数据同一样本 {sample_id!r} 归属多个分组："
                    f"{previous!r} 与 {group!r}"
                )
        else:
            sample_group[sample_id] = group
        samples[cell_id] = sample_id
        groups[cell_id] = group

    missing = [cell_id for cell_id in cell_ids if cell_id not in samples]
    if missing:
        _fail(
            f"重复元数据未覆盖矩阵全部细胞：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个"
        )
    matrix_cells = set(cell_ids)
    extra = [cell_id for cell_id in samples if cell_id not in matrix_cells]
    if extra:
        _fail(
            f"重复元数据包含矩阵之外的细胞：{extra[0]!r} 等 {len(extra)} 个"
        )

    group_sizes: Dict[str, int] = {}
    for group in groups.values():
        group_sizes[group] = group_sizes.get(group, 0) + 1
    # 样本序与样本规模均按输入细胞首次出现顺序（而非文件行序），
    # 使同归属关系的行重排不改变任何确定性输出
    sample_order: List[str] = []
    sample_sizes: Dict[str, int] = {}
    seen_samples = set()
    for cell_id in cell_ids:
        sample_id = samples[cell_id]
        sample_sizes[sample_id] = sample_sizes.get(sample_id, 0) + 1
        if sample_id not in seen_samples:
            seen_samples.add(sample_id)
            sample_order.append(sample_id)

    return ReplicateMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        samples=samples,
        groups=groups,
        sample_group=sample_group,
        sample_order=sample_order,
        group_sizes={group: group_sizes[group] for group in sorted(group_sizes)},
        sample_sizes=sample_sizes,
    )


def build_pseudobulk(
    matrix: ExpressionMatrix,
    qc: QCResult,
    replicate: ReplicateMetadata,
) -> PseudobulkData:
    """按质控后保留细胞与基因汇总 pseudobulk 并做样本级分组差异表达。

    样本按其输入细胞在矩阵列序中的首次出现排列；质控后无保留细胞的样本
    不进入汇总。原始计数只在保留基因上求和，样本文库即该样本在保留基因
    上的总计数，归一到 10000 后取 log1p。随后以样本为观测值复用分组
    比较口径（one-vs-rest + 升序两两、Welch t、比较内 BH）。
    """
    # 保留细胞 -> 样本，按保留细胞原顺序收集；样本首次出现顺序与全输入
    # 列序中的首次出现顺序一致（保留细胞是其一个子序列）
    cells_of_sample: Dict[str, List[int]] = {}
    kept_sample_order: List[str] = []
    for c in qc.kept_cells:
        sample_id = replicate.samples[matrix.cell_ids[c]]
        if sample_id not in cells_of_sample:
            cells_of_sample[sample_id] = []
            kept_sample_order.append(sample_id)
        cells_of_sample[sample_id].append(c)

    empty_samples = [
        sample_id
        for sample_id in replicate.sample_order
        if sample_id not in cells_of_sample
    ]
    if empty_samples:
        _data_fail(
            f"质控后样本 {empty_samples[0]!r} 等 {len(empty_samples)} 个"
            f"无保留细胞，无法汇总 pseudobulk"
        )

    sample_groups = [
        replicate.sample_group[sample_id] for sample_id in kept_sample_order
    ]
    group_to_samples: Dict[str, List[str]] = {}
    for sample_id, group in zip(kept_sample_order, sample_groups):
        group_to_samples.setdefault(group, []).append(sample_id)
    present_groups = sorted(group_to_samples)
    if len(present_groups) < 2:
        _data_fail(
            f"质控后非空分组仅 {len(present_groups)} 个，不足两个，"
            f"无法进行 pseudobulk 差异表达"
        )
    single_replicate_groups = [
        group
        for group in present_groups
        if len(group_to_samples[group]) < 2
    ]
    if single_replicate_groups:
        group = single_replicate_groups[0]
        _data_fail(
            f"分组 {group!r} 质控后有效重复仅 "
            f"{len(group_to_samples[group])} 个，少于两个，"
            f"无法进行 pseudobulk 差异表达"
        )

    gene_ids = [matrix.gene_ids[g] for g in qc.kept_genes]
    n_samples = len(kept_sample_order)
    counts: List[List[int]] = []
    for g in qc.kept_genes:
        row = matrix.counts[g]
        bulk_row: List[int] = []
        for sample_id in kept_sample_order:
            total = 0
            for c in cells_of_sample[sample_id]:
                total += row[c]
            bulk_row.append(total)
        counts.append(bulk_row)

    sample_totals = [0] * n_samples
    for row in counts:
        for s, value in enumerate(row):
            sample_totals[s] += value

    values: List[List[float]] = []
    for row in counts:
        normalized_row: List[float] = []
        for s, count in enumerate(row):
            total = sample_totals[s]
            scaled = count / total * NORMALIZATION_TARGET if total > 0 else 0.0
            normalized_row.append(math.log1p(scaled))
        values.append(normalized_row)

    comparisons = run_group_comparisons(gene_ids, values, sample_groups)

    return PseudobulkData(
        gene_ids=gene_ids,
        sample_ids=list(kept_sample_order),
        sample_groups=sample_groups,
        counts=counts,
        sample_totals=sample_totals,
        values=values,
        comparisons=comparisons,
    )


def _data_fail(message: str) -> None:
    raise CellFlowDataError(message)
