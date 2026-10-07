"""可选二维邻域嵌入：质控后保留细胞 PCA 得分上的确定性 UMAP 式降维。

启用 ``--umap`` 后，以本次分析实际用于聚类的 PCA 细胞得分的欧氏距离为
每个细胞确定 k 个近邻（k 取 15 与细胞数减一的较小值，距离并列按
cell_id 升序取先），把单向近邻关系按 UMAP 模糊单纯形口径对称化为
连接强度，再由 ``--seed`` 派生的确定性随机流（自带 MT19937）初始化
两维坐标并做固定轮数的吸引/排斥优化。随机初值、负样本选择与每轮
更新顺序全部来自该随机流，同输入、同参数、同种子下坐标、邻域权重与
迭代计数逐比特可复现，不依赖字典顺序、并发或外部随机状态。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .prng import MT19937

MAX_NEIGHBORS = 15
N_EPOCHS = 200
NEGATIVE_SAMPLES = 5
SIGMA_SEARCH_ITERATIONS = 64
GRAD_CLIP = 4.0
INIT_SPAN = 10.0
MIN_DIST_OFFSET = 0.001


@dataclass(frozen=True)
class UMAPNeighbor:
    """一条有向邻域记录：源细胞的某个近邻及其对称化后的连接强度。"""

    neighbor_id: str
    weight: float


@dataclass(frozen=True)
class UMAPResult:
    """二维邻域嵌入结果。"""

    cell_ids: List[str]                    # 与 PCA 细胞顺序一致
    coords: List[Tuple[float, float]]      # 与 cell_ids 对齐的二维坐标
    neighbors: List[List[UMAPNeighbor]]    # 每个源细胞的邻域（weight 降序、id 升序）
    n_neighbors: int                       # 实际邻居数 k
    n_iterations: int                      # 优化轮数


def _knn(
    scores: Sequence[Sequence[float]], cell_ids: Sequence[str], k: int
) -> List[List[Tuple[int, float]]]:
    """每个细胞的 k 个近邻（下标, 欧氏距离），距离升序、并列按 cell_id 升序。"""
    n = len(scores)
    result: List[List[Tuple[int, float]]] = []
    for i in range(n):
        xi = scores[i]
        dists: List[Tuple[float, str, int]] = []
        for j in range(n):
            if j == i:
                continue
            total = 0.0
            xj = scores[j]
            for d in range(len(xi)):
                diff = xi[d] - xj[d]
                total += diff * diff
            dists.append((math.sqrt(total), cell_ids[j], j))
        dists.sort(key=lambda item: (item[0], item[1]))
        result.append([(j, dist) for dist, _, j in dists[:k]])
    return result


def _directed_weights(distances: Sequence[float], k: int) -> List[float]:
    """单个细胞的单向模糊隶属度：exp(-(d - rho) / sigma)。

    rho 为最近邻距离；sigma 以固定轮数的二分搜索使权重和逼近 log2(k)。
    k 为 1 时唯一邻居权重取 1。
    """
    if k == 1:
        return [1.0]
    rho = distances[0]
    target = math.log2(k)
    lo = 0.0
    hi = float("inf")
    mid = 1.0
    for _ in range(SIGMA_SEARCH_ITERATIONS):
        psum = 0.0
        for d in distances:
            psum += math.exp(-(d - rho) / mid)
        if psum > target:
            # 权重和过大：sigma 偏大，收缩上界
            hi = mid
            mid = (lo + hi) / 2.0
        else:
            lo = mid
            if hi == float("inf"):
                mid = mid * 2.0
            else:
                mid = (lo + hi) / 2.0
    return [math.exp(-(d - rho) / mid) for d in distances]


def _symmetrize(
    knn: List[List[Tuple[int, float]]],
    directed: List[List[float]],
) -> Dict[Tuple[int, int], float]:
    """把单向近邻关系对称化为无向模糊连接强度：w = p + q - p * q。

    键为 (i, j) 且 i < j；只在单边出现的边保留单边权重。
    """
    edges: Dict[Tuple[int, int], float] = {}
    for i, neighbors in enumerate(knn):
        for (j, _), p in zip(neighbors, directed[i]):
            key = (i, j) if i < j else (j, i)
            existing = edges.get(key)
            if existing is None:
                edges[key] = p
            else:
                edges[key] = existing + p - existing * p
    return edges


def _optimize(
    edges: List[Tuple[int, int, float]], n: int, rng: MT19937
) -> List[Tuple[float, float]]:
    """固定轮数的吸引/排斥优化；初值、边顺序与负样本全部来自 rng。"""
    coords = [
        [rng.random() * 2.0 * INIT_SPAN - INIT_SPAN,
         rng.random() * 2.0 * INIT_SPAN - INIT_SPAN]
        for _ in range(n)
    ]
    order = list(range(len(edges)))
    for epoch in range(N_EPOCHS):
        # 学习率从 1.0 线性退火到接近 0
        lr = 1.0 - epoch / N_EPOCHS
        # 每轮以 Fisher-Yates 打乱边更新顺序（随机流只来自 rng）
        for t in range(len(order) - 1, 0, -1):
            s = rng.randbelow(t + 1)
            order[t], order[s] = order[s], order[t]
        for edge_index in order:
            a, b, w = edges[edge_index]
            pa = coords[a]
            pb = coords[b]
            dx = pa[0] - pb[0]
            dy = pa[1] - pb[1]
            d2 = dx * dx + dy * dy
            # 吸引项：沿边把两端拉近，强度随边权与距离衰减
            coeff = 2.0 * w / (1.0 + d2) * lr
            if coeff > GRAD_CLIP:
                coeff = GRAD_CLIP
            pa[0] -= coeff * dx
            pa[1] -= coeff * dy
            pb[0] += coeff * dx
            pb[1] += coeff * dy
            # 排斥项：每条边配固定数量的负样本，把源点推离
            for _ in range(NEGATIVE_SAMPLES):
                m = rng.randbelow(n)
                if m == a:
                    continue
                pm = coords[m]
                dx = pa[0] - pm[0]
                dy = pa[1] - pm[1]
                d2 = dx * dx + dy * dy
                coeff = (
                    2.0 / ((MIN_DIST_OFFSET + d2) * (1.0 + d2)) * lr
                )
                if coeff > GRAD_CLIP:
                    coeff = GRAD_CLIP
                pa[0] += coeff * dx
                pa[1] += coeff * dy
    return [(xy[0], xy[1]) for xy in coords]


def compute_umap(
    scores: Sequence[Sequence[float]],
    cell_ids: Sequence[str],
    *,
    seed: int,
) -> UMAPResult:
    """在 PCA 细胞得分上计算确定性的二维邻域嵌入。

    调用方保证细胞数 >= 2（管线在质控后已拒绝不足两个细胞的数据）。
    """
    n = len(cell_ids)
    if n < 2:
        # 交由管线统一包装为 CellFlowDataError
        raise ValueError(f"UMAP 嵌入至少需要两个细胞，得到 {n} 个")

    k = min(MAX_NEIGHBORS, n - 1)
    knn = _knn(scores, cell_ids, k)
    directed = [
        _directed_weights([dist for _, dist in neighbors], k)
        for neighbors in knn
    ]
    edges = _symmetrize(knn, directed)

    # 优化输入的边集按下标对排序，更新顺序再交给随机流打乱
    edge_list = sorted((a, b, w) for (a, b), w in edges.items())
    rng = MT19937(seed)
    coords = _optimize(edge_list, n, rng)

    # 每个源细胞的有向邻域记录：weight 为对称化后的连接强度，
    # 按 weight 降序、neighbor_id 升序排列
    neighbors: List[List[UMAPNeighbor]] = []
    for i in range(n):
        records = [
            UMAPNeighbor(
                neighbor_id=cell_ids[j],
                weight=edges[(i, j) if i < j else (j, i)],
            )
            for j, _ in knn[i]
        ]
        records.sort(key=lambda r: (-r.weight, r.neighbor_id))
        neighbors.append(records)

    return UMAPResult(
        cell_ids=list(cell_ids),
        coords=coords,
        neighbors=neighbors,
        n_neighbors=k,
        n_iterations=N_EPOCHS,
    )
