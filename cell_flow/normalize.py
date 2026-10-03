"""按细胞文库大小归一化、对数变换，以及高变基因选择。"""

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from .batch import center_by_batch
from .io import ExpressionMatrix
from .qc import QCResult

NORMALIZATION_TARGET = 10000.0
HVG_BINS = 20


@dataclass(frozen=True)
class NormalizedData:
    """质控后矩阵的 log-归一化表达。"""

    gene_ids: List[str]          # 保留基因（原矩阵行序）
    cell_ids: List[str]          # 保留细胞（原矩阵列序）
    values: List[List[float]]    # values[gene][cell] = ln(count/total*1e4 + 1)
    cell_totals: List[int]       # 与 cell_ids 对齐的原始文库总计数
    selected_genes: List[int]    # 高变基因在 values 中的行下标（原行序）
    # 批次均值中心化后的表达（与 values 同形）；未做批次校正时为 None
    corrected_values: Optional[List[List[float]]] = None

    @property
    def analysis_values(self) -> List[List[float]]:
        """下游分析（HVG/PCA/聚类/差异表达/图表）使用的表达矩阵。"""
        return (
            self.corrected_values
            if self.corrected_values is not None
            else self.values
        )


def normalize_and_select_hvg(
    matrix: ExpressionMatrix,
    qc: QCResult,
    *,
    n_hvg: int,
    batch_labels: Optional[List[str]] = None,
) -> NormalizedData:
    """归一化并选择高变基因。

    ``batch_labels`` 与保留细胞（列序）对齐；提供时先对 log 归一化值做
    批次均值中心化，高变基因在校正值上选择，校正值随结果一并返回。
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

    corrected = None
    if batch_labels is not None and len(set(batch_labels)) >= 2:
        corrected = center_by_batch(values, batch_labels)
    selected = _select_highly_variable(
        gene_ids, corrected if corrected is not None else values, n_hvg
    )
    return NormalizedData(
        gene_ids=gene_ids,
        cell_ids=cell_ids,
        values=values,
        cell_totals=totals,
        selected_genes=selected,
        corrected_values=corrected,
    )


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
