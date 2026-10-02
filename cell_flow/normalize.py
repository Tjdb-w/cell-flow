"""按细胞文库大小归一化、对数变换，以及高变基因选择。

提供 ``--batch-metadata`` 时，先在质控后保留细胞上按基因做批次均值
中心化：每个值减去该细胞所属批次在该基因上的保留细胞均值，再加回
全部保留细胞的该基因总均值。校正值用于高变基因选择及后续 PCA、
聚类、差异表达与图表数据；未校正的 log 归一化值仍原样写往
``normalized_expression.tsv``。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .io import ExpressionMatrix
from .qc import QCResult

NORMALIZATION_TARGET = 10000.0
HVG_BINS = 20


@dataclass(frozen=True)
class NormalizedData:
    """质控后矩阵的 log-归一化表达。

    ``values`` 始终是未经批次校正的 log 归一化值；提供批次元数据时
    ``batch_corrected_values`` 为与其同形的按基因批次均值中心化结果，
    高变基因选择及后续分析改用校正值，``batch_labels`` 与 ``cell_ids``
    对齐记录每个保留细胞的批次。无批次时二者分别为 None。
    """

    gene_ids: List[str]          # 保留基因（原矩阵行序）
    cell_ids: List[str]          # 保留细胞（原矩阵列序）
    values: List[List[float]]    # values[gene][cell] = ln(count/total*1e4 + 1)
    cell_totals: List[int]       # 与 cell_ids 对齐的原始文库总计数
    selected_genes: List[int]    # 高变基因在 values 中的行下标（原行序）
    batch_labels: Optional[List[str]] = None
    batch_corrected_values: Optional[List[List[float]]] = None

    @property
    def analysis_values(self) -> List[List[float]]:
        """驱动高变基因之后全部分析的表达矩阵：有校正用校正值，否则原值。"""
        if self.batch_corrected_values is not None:
            return self.batch_corrected_values
        return self.values


def log_normalize(
    matrix: ExpressionMatrix,
    qc: QCResult,
) -> Tuple[List[str], List[str], List[int], List[List[float]]]:
    """对质控后保留基因 × 保留细胞做文库大小归一化与 log1p 变换。

    返回 (gene_ids, cell_ids, cell_totals, values)；values[gene][cell]。
    """
    gene_idx = qc.kept_genes
    cell_idx = qc.kept_cells
    cell_ids = [matrix.cell_ids[c] for c in cell_idx]
    gene_ids = [matrix.gene_ids[g] for g in gene_idx]

    # 文库大小按该细胞在全部输入基因上的总计数
    totals = [qc.cell_qc[c].total_counts for c in cell_idx]

    values: List[List[float]] = []
    for g in gene_idx:
        row = matrix.counts[g]
        normalized_row: List[float] = []
        for ci, c in enumerate(cell_idx):
            total = totals[ci]
            scaled = row[c] / total * NORMALIZATION_TARGET if total > 0 else 0.0
            normalized_row.append(math.log1p(scaled))
        values.append(normalized_row)

    return gene_ids, cell_ids, totals, values


def normalize_and_select_hvg(
    matrix: ExpressionMatrix,
    qc: QCResult,
    *,
    n_hvg: int,
    cell_batches: Optional[List[str]] = None,
) -> NormalizedData:
    gene_ids, cell_ids, totals, values = log_normalize(matrix, qc)

    batch_labels: Optional[List[str]] = None
    corrected: Optional[List[List[float]]] = None
    if cell_batches is not None:
        # 批次标签与保留细胞列序对齐；校正值驱动高变基因选择
        batch_labels = list(cell_batches)
        corrected = correct_batch_means(values, batch_labels)
        selected = _select_highly_variable(gene_ids, corrected, n_hvg)
    else:
        selected = _select_highly_variable(gene_ids, values, n_hvg)

    return NormalizedData(
        gene_ids=gene_ids,
        cell_ids=cell_ids,
        values=values,
        cell_totals=totals,
        selected_genes=selected,
        batch_labels=batch_labels,
        batch_corrected_values=corrected,
    )


def correct_batch_means(
    values: List[List[float]],
    cell_batches: Sequence[str],
) -> List[List[float]]:
    """按基因做批次均值中心化。

    对每个基因 g、批次 b 与该批次保留细胞集合 C_b：
    ``corrected[g][c] = values[g][c] - mean_b(values[g]) + mean_all(values[g])``，
    其中批次均值仅用该批次质控后保留细胞，总均值用全部保留细胞。
    中心化保持各基因的保留细胞总均值不变，只去除批次间的均值位移。
    """
    n_cells = len(cell_batches)
    members: Dict[str, List[int]] = {}
    for c, batch in enumerate(cell_batches):
        members.setdefault(batch, []).append(c)

    # 各批次内细胞数（与求和顺序无关，先求一次倒数避免重复除法）
    batch_inv = {batch: 1.0 / len(idx) for batch, idx in members.items()}
    inv_all = 1.0 / n_cells if n_cells > 0 else 0.0

    corrected: List[List[float]] = []
    for row in values:
        overall = sum(row) * inv_all
        batch_mean: Dict[str, float] = {}
        for batch, idx in members.items():
            batch_mean[batch] = sum(row[c] for c in idx) * batch_inv[batch]
        corrected_row = [
            row[c] - batch_mean[cell_batches[c]] + overall for c in range(n_cells)
        ]
        corrected.append(corrected_row)
    return corrected


def _gene_mean_var(row: Sequence[float]) -> Tuple[float, float]:
    n = len(row)
    mu = sum(row) / n
    if n < 2:
        return mu, 0.0
    sq = 0.0
    for x in row:
        d = x - mu
        sq += d * d
    return mu, sq / (n - 1)


def _select_highly_variable(
    gene_ids: List[str],
    values: List[List[float]],
    n_hvg: int,
) -> List[int]:
    """Seurat 风格高变基因：log 表达的均值-离散度，按均值等频分箱后
    对离散度做箱内 z 标准化，取 z 最高的 ``n_hvg`` 个。
    返回所选基因的行下标，按原行序排列。
    """
    n_genes = len(values)
    means: List[float] = [0.0] * n_genes
    dispersions: List[float] = [0.0] * n_genes
    for g, row in enumerate(values):
        mu, var = _gene_mean_var(row)
        means[g] = mu
        dispersions[g] = var / mu if mu > 0.0 else 0.0

    # 等频分箱：按 (均值, 原下标) 排序后尽量均分到各箱
    order = sorted(range(n_genes), key=lambda g: (means[g], g))
    nbins = min(HVG_BINS, n_genes)
    z_scores = [0.0] * n_genes
    base, extra = divmod(n_genes, nbins)
    start = 0
    for b in range(nbins):
        size = base + (1 if b < extra else 0)
        members = order[start:start + size]
        start += size
        if len(members) < 2:
            for g in members:
                z_scores[g] = 0.0
            continue
        bin_disp = [dispersions[g] for g in members]
        bmean = sum(bin_disp) / len(bin_disp)
        bvar = sum((d - bmean) ** 2 for d in bin_disp) / (len(bin_disp) - 1)
        bstd = math.sqrt(bvar)
        for g in members:
            z_scores[g] = (dispersions[g] - bmean) / bstd if bstd > 0.0 else 0.0

    # 按 (z 降序, 基因 ID 升序) 确定性挑选
    ranked = sorted(range(n_genes), key=lambda g: (-z_scores[g], gene_ids[g]))
    k = min(n_hvg, n_genes)
    return sorted(ranked[:k])
