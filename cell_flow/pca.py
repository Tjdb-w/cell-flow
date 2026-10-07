"""主成分分析。

在高变基因的细胞中心化 log 表达上，计算细胞 x 细胞协方差矩阵，
用 :mod:`cell_flow.linalg` 的 Jacobi 分解求特征向量，再投影回细胞坐标。
全部 float64、固定符号约定，结果确定。
"""

import math
from dataclasses import dataclass
from typing import List

from .linalg import jacobi_eigh
from .normalize import NormalizedData

MAX_PCS = 20


@dataclass(frozen=True)
class PCAResult:
    n_pcs: int
    scores: List[List[float]]  # scores[cell][pc]
    explained_variance: List[float]  # 长度 n_pcs，对应特征值（降序主成分）
    cell_ids: List[str]


@dataclass(frozen=True)
class PCALoadings:
    n_pcs: int
    gene_ids: List[str]          # 高变基因（与 PCA 实际使用的选择顺序一致）
    loadings: List[List[float]]  # loadings[gene][pc]，带符号载荷


def _center_selected(data: NormalizedData) -> List[List[float]]:
    """中心化的细胞 x 高变基因矩阵（批次校正时基于校正值）。"""
    selected = data.selected_genes
    n_cells = len(data.cell_ids)
    n_hvg = len(selected)
    centered: List[List[float]] = [[0.0] * n_hvg for _ in range(n_cells)]
    for j, g in enumerate(selected):
        row = data.analysis_values[g]
        mu = sum(row) / n_cells
        for c in range(n_cells):
            centered[c][j] = row[c] - mu
    return centered


def run_pca(data: NormalizedData, n_pcs: int = MAX_PCS) -> PCAResult:
    selected = data.selected_genes
    n_cells = len(data.cell_ids)
    n_hvg = len(selected)
    # 实际主成分数取请求值、细胞数减一、高变基因数三者的最小值
    n_pcs = min(n_pcs, n_cells - 1, n_hvg)
    if n_pcs < 1:
        # 交由管线统一包装为 CellFlowDataError
        raise ValueError("PCA 无法成立：主成分数为 0")

    # 中心化的细胞 x 高变基因矩阵（批次校正时基于校正值）
    centered = _center_selected(data)

    # 细胞 x 细胞协方差（1/(n-1) 缩放对特征向量无影响，仅影响特征值）
    cov = [[0.0] * n_cells for _ in range(n_cells)]
    scale = 1.0 / (n_cells - 1) if n_cells > 1 else 1.0
    for i in range(n_cells):
        xi = centered[i]
        for j in range(i, n_cells):
            xj = centered[j]
            dot = 0.0
            for k in range(n_hvg):
                dot += xi[k] * xj[k]
            value = dot * scale
            cov[i][j] = value
            cov[j][i] = value

    eigenvalues, eigenvectors = jacobi_eigh(cov)
    # jacobi_eigh 返回升序；取最大的 n_pcs 个并转为降序
    top_values = list(reversed(eigenvalues[-n_pcs:]))
    top_vectors = list(reversed(eigenvectors[-n_pcs:]))

    scores = [[0.0] * n_pcs for _ in range(n_cells)]
    explained: List[float] = []
    for k in range(n_pcs):
        eig = top_values[k]
        explained.append(eig if eig > 0.0 else 0.0)
        # Gram(XX^T) 的单位特征向量为 v，则标准 PCA 坐标 = sqrt(lambda) * v
        factor = math.sqrt(eig) if eig > 0.0 else 0.0
        vec = top_vectors[k]
        for c in range(n_cells):
            scores[c][k] = factor * vec[c]
        # 符号约定：绝对值最大的得分为正（并列取第一个）
        pivot = 0
        pivot_abs = abs(scores[0][k])
        for c in range(1, n_cells):
            value_abs = abs(scores[c][k])
            if value_abs > pivot_abs:
                pivot_abs = value_abs
                pivot = c
        if scores[pivot][k] < 0.0:
            for c in range(n_cells):
                scores[c][k] = -scores[c][k]

    return PCAResult(
        n_pcs=n_pcs,
        scores=scores,
        explained_variance=explained,
        cell_ids=list(data.cell_ids),
    )


def compute_pca_loadings(data: NormalizedData, pca: PCAResult) -> PCALoadings:
    """高变基因在各主成分上的带符号载荷。

    沿用 PCA 实际使用的最终细胞、高变基因顺序与分析表达值（批次均值
    中心化后用校正值，否则用 log 归一化值）。基因 g 在 PCk 上的载荷为
    其中心化表达向量与 PCk 细胞得分向量的内积除以解释方差的平方根；
    解释方差为 0 时该主成分载荷全部取 0。
    """
    selected = data.selected_genes
    n_cells = len(data.cell_ids)
    centered = _center_selected(data)
    gene_ids = [data.gene_ids[g] for g in selected]
    loadings: List[List[float]] = [[0.0] * pca.n_pcs for _ in selected]
    for k in range(pca.n_pcs):
        variance = pca.explained_variance[k]
        if variance <= 0.0:
            # 解释方差为 0：该 PC 载荷全部保持 0
            continue
        denom = math.sqrt(variance)
        for j in range(len(selected)):
            dot = 0.0
            for c in range(n_cells):
                dot += centered[c][j] * pca.scores[c][k]
            loadings[j][k] = dot / denom
    return PCALoadings(n_pcs=pca.n_pcs, gene_ids=gene_ids, loadings=loadings)
