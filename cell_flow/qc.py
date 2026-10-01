"""细胞与基因水平的质量控制指标与保留标记。"""

from dataclasses import dataclass
from typing import List

from .io import ExpressionMatrix


@dataclass(frozen=True)
class CellQC:
    cell_id: str
    total_counts: int
    detected_genes: int
    mitochondrial_counts: int
    mitochondrial_fraction: float
    retained: bool


@dataclass(frozen=True)
class GeneQC:
    gene_id: str
    total_counts: int
    detected_cells: int
    retained: bool


@dataclass(frozen=True)
class QCResult:
    cell_qc: List[CellQC]          # 与输入列顺序一致，含全部细胞
    gene_qc: List[GeneQC]          # 与输入行顺序一致，含全部基因
    kept_cells: List[int]          # 保留细胞在原矩阵中的列下标
    kept_genes: List[int]          # 保留基因在原矩阵中的行下标


def compute_qc(
    matrix: ExpressionMatrix,
    *,
    min_genes: int,
    max_mito_fraction: float,
    min_cells: int,
    mito_prefix: str,
) -> QCResult:
    counts = matrix.counts
    n_genes = matrix.n_genes
    n_cells = matrix.n_cells

    cell_totals = [0] * n_cells
    cell_detected = [0] * n_cells
    cell_mito = [0] * n_cells
    gene_totals = [0] * n_genes
    gene_detected = [0] * n_genes

    for g, row in enumerate(counts):
        gene_id = matrix.gene_ids[g]
        is_mito = gene_id.startswith(mito_prefix)
        gtotal = 0
        gdet = 0
        for c, value in enumerate(row):
            gtotal += value
            cell_totals[c] += value
            if value > 0:
                gdet += 1
                cell_detected[c] += 1
                if is_mito:
                    cell_mito[c] += value
        gene_totals[g] = gtotal
        gene_detected[g] = gdet

    cell_qc: List[CellQC] = []
    kept_cells: List[int] = []
    for c, cell_id in enumerate(matrix.cell_ids):
        total = cell_totals[c]
        fraction = cell_mito[c] / total if total > 0 else 0.0
        retained = (
            cell_detected[c] >= min_genes and fraction <= max_mito_fraction
        )
        cell_qc.append(
            CellQC(
                cell_id=cell_id,
                total_counts=total,
                detected_genes=cell_detected[c],
                mitochondrial_counts=cell_mito[c],
                mitochondrial_fraction=fraction,
                retained=retained,
            )
        )
        if retained:
            kept_cells.append(c)

    gene_qc: List[GeneQC] = []
    kept_genes: List[int] = []
    for g, gene_id in enumerate(matrix.gene_ids):
        retained = gene_detected[g] >= min_cells
        gene_qc.append(
            GeneQC(
                gene_id=gene_id,
                total_counts=gene_totals[g],
                detected_cells=gene_detected[g],
                retained=retained,
            )
        )
        if retained:
            kept_genes.append(g)

    return QCResult(
        cell_qc=cell_qc,
        gene_qc=gene_qc,
        kept_cells=kept_cells,
        kept_genes=kept_genes,
    )
