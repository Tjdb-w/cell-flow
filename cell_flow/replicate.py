"""按生物学重复汇总的 pseudobulk 差异表达：读取、汇总与检验。

``--replicate-metadata`` 指向 UTF-8 制表符文本，表头恰为 ``cell_id``、
``sample_id``、``group`` 三列；每个输入细胞恰好一行，``cell_id`` 唯一，
``sample_id`` 与 ``group`` 非空，同一 ``sample_id`` 只能属于一个
``group``。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名无关）。
任何不合法（表头不符、行列数不符、ID 重复或为空、覆盖不全、多出矩阵之外
的细胞、同一样本归属多个分组，或 gzip 多成员、尾随数据、截断、校验失败）
都抛 :class:`cell_flow.errors.CellFlowInputError`（退出码 2）。

汇总只使用既有 QC 保留的细胞与基因：按 ``sample_id`` 对原始计数求和，
把每个样本文库归一到 10000 后做 log1p，样本顺序按该样本细胞在输入矩阵中
首次出现的顺序。每个 group 先做 one-vs-rest（该组全部样本对比其余全部
样本），再按 group 升序两两比较；以样本为观测值做 Welch t 检验与双侧
P 值，并在每个比较内跨基因做 BH 校正。

质控后样本无细胞、某 group 有效重复少于两个，或不足两个 group，抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。
"""

import hashlib
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .errors import CellFlowDataError, CellFlowInputError
from .io import ExpressionMatrix, GZIP_MAGIC, gunzip_single
from .linalg import benjamini_hochberg, welch_ttest
from .markers import ONE_VS_REST, PAIRWISE
from .qc import QCResult

PSEUDOBULK_TARGET = 10000.0


@dataclass(frozen=True)
class ReplicateMetadata:
    """细胞 -> 样本/分组归属；每个输入细胞恰好一行。"""

    path: str
    name: str                    # 文件名（路径末段）
    sha256: str                  # 原始字节（gzip 即压缩字节）的 SHA-256
    samples: Dict[str, str]      # cell_id -> sample_id
    groups: Dict[str, str]       # cell_id -> group
    sample_group: Dict[str, str]  # sample_id -> 唯一 group
    sample_sizes: Dict[str, int]  # sample_id -> 输入细胞数（按 sample_id 升序）
    group_sizes: Dict[str, int]   # group -> 输入细胞数（按 group 升序）


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
        _fail("重复元数据表头必须恰为 cell_id、sample_id、group 三列")

    samples: Dict[str, str] = {}
    groups: Dict[str, str] = {}
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

    sample_group: Dict[str, str] = {}
    sample_sizes: Dict[str, int] = {}
    # 按输入细胞列序遍历：同一组冲突时先出现的样本归属决定报错的确定性
    for cell_id in cell_ids:
        sample_id = samples[cell_id]
        group = groups[cell_id]
        if sample_id in sample_group:
            if sample_group[sample_id] != group:
                _fail(
                    f"同一样本归属多个分组：{sample_id!r} 同时属于 "
                    f"{sample_group[sample_id]!r} 与 {group!r}"
                )
        else:
            sample_group[sample_id] = group
        sample_sizes[sample_id] = sample_sizes.get(sample_id, 0) + 1

    group_sizes: Dict[str, int] = {}
    for cell_id in cell_ids:
        group = groups[cell_id]
        group_sizes[group] = group_sizes.get(group, 0) + 1

    return ReplicateMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        samples=samples,
        groups=groups,
        sample_group=sample_group,
        sample_sizes={s: sample_sizes[s] for s in sorted(sample_sizes)},
        group_sizes={g: group_sizes[g] for g in sorted(group_sizes)},
    )


@dataclass(frozen=True)
class PseudobulkExpression:
    """质控后保留细胞/基因按样本汇总并 log 归一化的表达矩阵。

    行（基因）保持 QC 保留基因的原矩阵行序；列（样本）按该样本细胞在
    输入矩阵中首次出现的顺序。
    """

    gene_ids: List[str]
    sample_ids: List[str]
    sample_groups: List[str]           # 与 sample_ids 对齐的 group
    values: List[List[float]]          # values[gene][sample] = ln(sum/lib*1e4 + 1)
    sample_cell_counts: List[int]      # 与 sample_ids 对齐的保留细胞数
    library_sizes: List[int]           # 与 sample_ids 对齐的保留基因总计数


@dataclass(frozen=True)
class PseudobulkMarkerRecord:
    comparison_type: str
    group_a: str
    group_b: str  # one-vs-rest 时为空
    gene_id: str
    mean_in_a: float
    mean_in_b: float
    log_fc_a_vs_b: float
    t_stat: float
    p_value: float
    p_value_adj: float


# (comparison_type, group_a, group_b, records)；one-vs-rest 的 group_b 为 ""
PseudobulkComparison = Tuple[str, str, str, List[PseudobulkMarkerRecord]]


def build_pseudobulk(
    matrix: ExpressionMatrix, qc: QCResult, replicate: ReplicateMetadata
) -> PseudobulkExpression:
    """按样本汇总 QC 保留细胞的原始计数，再做文库归一与 log1p。"""
    # 样本列序：按该样本保留细胞在输入矩阵中的首次出现顺序
    sample_ids: List[str] = []
    sample_index: Dict[str, int] = {}
    for c in qc.kept_cells:
        sample_id = replicate.samples[matrix.cell_ids[c]]
        if sample_id not in sample_index:
            sample_index[sample_id] = len(sample_ids)
            sample_ids.append(sample_id)
    if not sample_ids:
        raise CellFlowDataError("质控后没有任何样本保留细胞，无法做 pseudobulk 分析")

    n_samples = len(sample_ids)
    sums: List[List[int]] = [
        [0] * n_samples for _ in qc.kept_genes
    ]
    cell_counts = [0] * n_samples
    for c in qc.kept_cells:
        s = sample_index[replicate.samples[matrix.cell_ids[c]]]
        cell_counts[s] += 1
        for gi, g in enumerate(qc.kept_genes):
            sums[gi][s] += matrix.counts[g][c]

    library_sizes = [0] * n_samples
    for row in sums:
        for s, value in enumerate(row):
            library_sizes[s] += value

    values: List[List[float]] = []
    for row in sums:
        normalized_row: List[float] = []
        for s, total in enumerate(library_sizes):
            # 与单细胞归一化同一口径：文库为零时该样本所有基因为 0.0
            scaled = row[s] / total * PSEUDOBULK_TARGET if total > 0 else 0.0
            normalized_row.append(math.log1p(scaled))
        values.append(normalized_row)

    gene_ids = [matrix.gene_ids[g] for g in qc.kept_genes]
    sample_groups = [replicate.sample_group[s] for s in sample_ids]
    return PseudobulkExpression(
        gene_ids=gene_ids,
        sample_ids=sample_ids,
        sample_groups=sample_groups,
        values=values,
        sample_cell_counts=cell_counts,
        library_sizes=library_sizes,
    )


def _comparison(
    expression: PseudobulkExpression,
    samples_a: List[int],
    samples_b: List[int],
    comparison_type: str,
    group_a: str,
    group_b: str,
) -> List[PseudobulkMarkerRecord]:
    """两组样本之间逐基因检验：均值差即 log_fc，Welch t 检验，
    并在本比较内跨基因做 BH 校正。"""
    records: List[PseudobulkMarkerRecord] = []
    pvalues: List[float] = []
    for g, gene_id in enumerate(expression.gene_ids):
        row = expression.values[g]
        values_a = [row[s] for s in samples_a]
        values_b = [row[s] for s in samples_b]
        mean_a = sum(values_a) / len(values_a)
        mean_b = sum(values_b) / len(values_b)
        log_fc = mean_a - mean_b
        t_stat, p_value = welch_ttest(values_a, values_b)
        pvalues.append(p_value)
        records.append(
            PseudobulkMarkerRecord(
                comparison_type=comparison_type,
                group_a=group_a,
                group_b=group_b,
                gene_id=gene_id,
                mean_in_a=mean_a,
                mean_in_b=mean_b,
                log_fc_a_vs_b=log_fc,
                t_stat=t_stat,
                p_value=p_value,
                p_value_adj=0.0,
            )
        )

    adjusted = benjamini_hochberg(pvalues)
    records = [
        PseudobulkMarkerRecord(
            comparison_type=r.comparison_type,
            group_a=r.group_a,
            group_b=r.group_b,
            gene_id=r.gene_id,
            mean_in_a=r.mean_in_a,
            mean_in_b=r.mean_in_b,
            log_fc_a_vs_b=r.log_fc_a_vs_b,
            t_stat=r.t_stat,
            p_value=r.p_value,
            p_value_adj=adjusted[g],
        )
        for g, r in enumerate(records)
    ]
    records.sort(key=lambda r: (r.p_value_adj, -r.log_fc_a_vs_b, r.gene_id))
    return records


def find_pseudobulk_markers(
    expression: PseudobulkExpression,
) -> List[PseudobulkComparison]:
    """按 group 的 pseudobulk 差异表达：每组 one-vs-rest，再按组升序两两比较。

    观测值为样本（pseudobulk log 归一化表达）。返回先全部 one-vs-rest
    （group 升序）、再全部 (a, b)（a < b，均按 group 升序）成对比较的
    比较列表；每个比较内记录按校正 P 值升序、log_fc_a_vs_b 降序、
    基因 ID 升序稳定排序，BH 校正以单个比较内的全部保留基因为一个家族。
    """
    group_to_samples: Dict[str, List[int]] = {}
    for s, group in enumerate(expression.sample_groups):
        group_to_samples.setdefault(group, []).append(s)

    present = sorted(group_to_samples)
    if len(present) < 2:
        raise CellFlowDataError(
            f"质控后有效分组仅 {len(present)} 个，不足两个，"
            f"无法进行 pseudobulk 差异表达"
        )
    too_few = [g for g in present if len(group_to_samples[g]) < 2]
    if too_few:
        raise CellFlowDataError(
            f"分组 {too_few[0]!r} 质控后仅有 "
            f"{len(group_to_samples[too_few[0]])} 个有效重复样本，"
            f"不足两个，无法进行 pseudobulk 差异表达"
        )

    n_samples = len(expression.sample_ids)
    comparisons: List[PseudobulkComparison] = []
    for group in present:
        in_samples = group_to_samples[group]
        out_samples = [
            s for s in range(n_samples) if expression.sample_groups[s] != group
        ]
        records = _comparison(
            expression, in_samples, out_samples, ONE_VS_REST, group, ""
        )
        comparisons.append((ONE_VS_REST, group, "", records))
    for i, group_a in enumerate(present):
        for group_b in present[i + 1:]:
            records = _comparison(
                expression,
                group_to_samples[group_a],
                group_to_samples[group_b],
                PAIRWISE,
                group_a,
                group_b,
            )
            comparisons.append((PAIRWISE, group_a, group_b, records))
    return comparisons
