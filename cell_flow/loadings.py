"""PCA 高变基因载荷（可解释性）。

沿用 :func:`cell_flow.pca.run_pca` 实际使用的最终细胞、高变基因顺序与
分析表达值（批次均值中心化时取校正值，否则取 log 归一化值），把每个
主成分坐标追溯到具体高变基因。

对 PCk，基因 g 的带符号载荷定义为该基因中心化表达向量与 PCk 得分向量的
内积除以解释方差的平方根；解释方差为 0 时该 PC 全部载荷取 0。
全部 float64，纯函数、结果确定。
"""

import math
from dataclasses import dataclass
from typing import List

from .normalize import NormalizedData
from .pca import PCAResult


@dataclass(frozen=True)
class PCALoadings:
    gene_ids: List[str]            # 高变基因（沿用既有高变基因顺序）
    loadings: List[List[float]]    # loadings[gene][pc]，带符号载荷


def compute_loadings(data: NormalizedData, pca: PCAResult) -> PCALoadings:
    selected = data.selected_genes
    n_cells = len(data.cell_ids)
    n_hvg = len(selected)
    n_pcs = pca.n_pcs

    # 与 run_pca 相同的中心化表达（analysis_values：校正后取校正值）
    centered: List[List[float]] = [[0.0] * n_cells for _ in range(n_hvg)]
    for j, g in enumerate(selected):
        row = data.analysis_values[g]
        mu = sum(row) / n_cells
        for c in range(n_cells):
            centered[j][c] = row[c] - mu

    loadings: List[List[float]] = [[0.0] * n_pcs for _ in range(n_hvg)]
    for j in range(n_hvg):
        x = centered[j]
        for k in range(n_pcs):
            variance = pca.explained_variance[k]
            if variance <= 0.0:
                # 解释方差为 0（或数值上不为正）时该 PC 载荷全部取 0
                continue
            dot = 0.0
            for c in range(n_cells):
                dot += x[c] * pca.scores[c][k]
            loadings[j][k] = dot / math.sqrt(variance)

    return PCALoadings(
        gene_ids=[data.gene_ids[g] for g in selected],
        loadings=loadings,
    )
