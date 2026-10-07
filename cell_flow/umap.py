"""可选二维邻域嵌入（UMAP 风格的确定性实现）。

仅当用户显式启用（``--umap``）时执行，不改变默认分析路径。
以质控后保留细胞的 PCA 得分为输入，按欧氏距离为每个细胞确定
``k = min(15, n_cells - 1)`` 个近邻（距离并列按 cell_id 升序取先），
为每个细胞拟合局部尺度后将单向近邻关系对称化为模糊连接强度，
再在二维空间做确定性优化。随机初值、负样本选择与逐轮更新顺序
全部由 ``--seed`` 经 :class:`cell_flow.prng.MT19937` 派生，不依赖
字典顺序、并发或外部随机状态，同输入同参数重复运行逐字节一致。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .prng import MT19937

N_NEIGHBORS = 15
N_EPOCHS = 200
N_NEGATIVE_SAMPLES = 5
INITIAL_ALPHA = 1.0
GAMMA = 1.0
# min_dist=0.1、spread=1.0 对应的标准曲线拟合常数
A = 1.57694346040507
B = 0.89506087791087
GRAD_CLIP = 4.0
SIGMA_ITERATIONS = 64
SIGMA_TOLERANCE = 1e-5
INIT_SPAN = 20.0


@dataclass(frozen=True)
class UmapResult:
    """二维邻域嵌入结果。

    ``cell_ids``、``embedding`` 与 ``neighbors`` 均按 PCA 细胞顺序
    （即 ``pca_scatter.tsv`` 行序）对齐；``neighbors[i]`` 为源细胞 i 的
    有向邻域记录 ``(neighbor_index, weight)``，已按 weight 降序、
    neighbor_id 升序排列，weight 为对称化后的模糊连接强度。
    """

    n_neighbors: int
    cell_ids: List[str]
    embedding: List[List[float]]  # embedding[cell] = [UMAP1, UMAP2]
    neighbors: List[List[Tuple[int, float]]]
    n_epochs: int

    @property
    def n_neighbor_records(self) -> int:
        return sum(len(records) for records in self.neighbors)


def _squared_distance(a: Sequence[float], b: Sequence[float]) -> float:
    total = 0.0
    for x, y in zip(a, b):
        d = x - y
        total += d * d
    return total


def _nearest_neighbors(
    scores: List[List[float]], cell_ids: List[str], k: int
) -> List[List[Tuple[int, float]]]:
    """每个细胞的 k 个欧氏最近邻（下标, 距离）；距离并列按 cell_id 升序取先。"""
    n = len(scores)
    neighbors: List[List[Tuple[int, float]]] = []
    for i in range(n):
        dists = []
        xi = scores[i]
        for j in range(n):
            if j == i:
                continue
            d = math.sqrt(_squared_distance(xi, scores[j]))
            dists.append((d, cell_ids[j], j))
        dists.sort(key=lambda item: (item[0], item[1]))
        neighbors.append([(j, d) for d, _, j in dists[:k]])
    return neighbors


def _smooth_sigma(distances: List[float], target: float) -> Tuple[float, float]:
    """二分拟合局部尺度 sigma，使邻居权重和逼近 log2(k)；返回 (rho, sigma)。

    rho 为到最近邻的距离；k == 1（target 为 0）时权重恒为 1，sigma 取 1。
    """
    rho = min(distances)
    if target <= 0.0:
        return rho, 1.0
    lo = 0.0
    hi = float("inf")
    mid = 1.0
    for _ in range(SIGMA_ITERATIONS):
        total = 0.0
        for d in distances:
            total += math.exp(-(d - rho) / mid)
        if abs(total - target) < SIGMA_TOLERANCE:
            break
        if total > target:
            hi = mid
            mid = (lo + hi) / 2.0
        else:
            lo = mid
            mid = mid * 2.0 if hi == float("inf") else (lo + hi) / 2.0
    return rho, mid


def _symmetrize(
    directed: List[List[Tuple[int, float]]], n: int
) -> List[Dict[int, float]]:
    """把单向模糊连接强度对称化：w = a + b - a * b（缺失方向按 0 计）。"""
    out: List[Dict[int, float]] = [dict() for _ in range(n)]
    for i in range(n):
        for j, w in directed[i]:
            out[i][j] = w
    sym: List[Dict[int, float]] = [dict() for _ in range(n)]
    for i in range(n):
        partners = set(out[i])
        for j in range(n):
            if i in out[j]:
                partners.add(j)
        for j in sorted(partners):
            a = out[i].get(j, 0.0)
            b = out[j].get(i, 0.0)
            sym[i][j] = a + b - a * b
    return sym


def _clip(value: float) -> float:
    if value > GRAD_CLIP:
        return GRAD_CLIP
    if value < -GRAD_CLIP:
        return -GRAD_CLIP
    return value


def _optimize(
    edges: List[Tuple[int, int, float]], n: int, seed: int
) -> List[List[float]]:
    """确定性二维布局优化：随机初值、逐轮边序与负样本全部由 seed 派生。"""
    rng = MT19937(seed)
    embedding = [
        [rng.random() * INIT_SPAN - INIT_SPAN / 2.0,
         rng.random() * INIT_SPAN - INIT_SPAN / 2.0]
        for _ in range(n)
    ]
    n_edges = len(edges)
    for epoch in range(N_EPOCHS):
        alpha = INITIAL_ALPHA * (1.0 - epoch / N_EPOCHS)
        # Fisher-Yates 洗牌确定本轮边处理顺序
        order = list(range(n_edges))
        for idx in range(n_edges - 1, 0, -1):
            j = rng.randbelow(idx + 1)
            order[idx], order[j] = order[j], order[idx]
        for edge_index in order:
            i, j, w = edges[edge_index]
            ei = embedding[i]
            ej = embedding[j]
            dx = ei[0] - ej[0]
            dy = ei[1] - ej[1]
            dist_sq = dx * dx + dy * dy
            if dist_sq > 0.0:
                coeff = (
                    -2.0 * A * B * (dist_sq ** (B - 1.0))
                    / (1.0 + A * dist_sq ** B)
                ) * w
            else:
                coeff = 0.0
            gx = _clip(coeff * dx)
            gy = _clip(coeff * dy)
            ei[0] += gx * alpha
            ei[1] += gy * alpha
            ej[0] -= gx * alpha
            ej[1] -= gy * alpha
            # 负样本排斥：只动当前点 i
            for _ in range(N_NEGATIVE_SAMPLES):
                m = rng.randbelow(n)
                if m == i:
                    continue
                em = embedding[m]
                dx = ei[0] - em[0]
                dy = ei[1] - em[1]
                dist_sq = dx * dx + dy * dy
                coeff = w * (2.0 * GAMMA * B) / (
                    (0.001 + dist_sq) * (1.0 + A * dist_sq ** B)
                )
                ei[0] += _clip(coeff * dx) * alpha
                ei[1] += _clip(coeff * dy) * alpha
    return embedding


def compute_umap(
    scores: List[List[float]], cell_ids: List[str], seed: int
) -> UmapResult:
    """从 PCA 得分计算二维邻域嵌入与对称化邻域记录。

    邻居数取 15 与细胞数减一的较小值；无法成立时抛 ``ValueError``，
    由管线统一包装为 :class:`CellFlowDataError`。
    """
    n = len(scores)
    if n < 2:
        raise ValueError(f"仅 {n} 个细胞，无法构成邻域嵌入")
    k = min(N_NEIGHBORS, n - 1)
    knn = _nearest_neighbors(scores, cell_ids, k)
    target = math.log2(k) if k > 1 else 0.0

    directed: List[List[Tuple[int, float]]] = []
    for i in range(n):
        distances = [d for _, d in knn[i]]
        rho, sigma = _smooth_sigma(distances, target)
        directed.append(
            [(j, math.exp(-(d - rho) / sigma)) for j, d in knn[i]]
        )

    sym = _symmetrize(directed, n)

    # 无向边（i < j，按 i、j 升序）供布局优化使用
    edges: List[Tuple[int, int, float]] = []
    for i in range(n):
        for j in sorted(sym[i]):
            if j > i:
                edges.append((i, j, sym[i][j]))

    embedding = _optimize(edges, n, seed)

    neighbors: List[List[Tuple[int, float]]] = []
    for i in range(n):
        records = sorted(
            sym[i].items(), key=lambda item: (-item[1], cell_ids[item[0]])
        )
        neighbors.append(records)

    return UmapResult(
        n_neighbors=k,
        cell_ids=list(cell_ids),
        embedding=embedding,
        neighbors=neighbors,
        n_epochs=N_EPOCHS,
    )
