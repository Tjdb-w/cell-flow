"""随机种子固定的 k-means 聚类。

k-means++ 初始化 + Lloyd 迭代，全部使用 :class:`cell_flow.prng.MT19937`
与确定性的并列打破规则，保证同输入同种子逐位可复现。
"""

from dataclasses import dataclass
from typing import List, Sequence

from .prng import MT19937

MAX_ITER = 300


@dataclass(frozen=True)
class KMeansResult:
    labels: List[int]       # 每个点的簇编号（0..k-1，按首次出现顺序规范化）
    centers: List[List[float]]
    iterations: int


def _squared_distance(a: Sequence[float], b: Sequence[float]) -> float:
    total = 0.0
    for x, y in zip(a, b):
        d = x - y
        total += d * d
    return total


def _kmeans_plus_plus(
    points: List[List[float]], k: int, rng: MT19937
) -> List[List[float]]:
    n = len(points)
    centers: List[List[float]] = []
    first = rng.randbelow(n)
    centers.append(list(points[first]))
    chosen_indices = {first}

    nearest_sq = [_squared_distance(p, centers[0]) for p in points]
    while len(centers) < k:
        total = sum(nearest_sq)
        if total <= 0.0:
            # 所有剩余点与中心重合：确定性挑选尚未选中的最小下标
            picked = next(
                (idx for idx in range(n) if idx not in chosen_indices), None
            )
            if picked is None:  # pragma: no cover - 上游已保证 k <= n
                break
            centers.append(list(points[picked]))
            chosen_indices.add(picked)
        else:
            threshold = rng.random() * total
            cumulative = 0.0
            chosen = n - 1
            for idx in range(n):
                cumulative += nearest_sq[idx]
                if cumulative >= threshold:
                    chosen = idx
                    break
            centers.append(list(points[chosen]))
            chosen_indices.add(chosen)
        new_center = centers[-1]
        for idx in range(n):
            d_sq = _squared_distance(points[idx], new_center)
            if d_sq < nearest_sq[idx]:
                nearest_sq[idx] = d_sq
    return centers


def _assign(
    points: Sequence[Sequence[float]], centers: Sequence[Sequence[float]]
) -> List[int]:
    labels = []
    for point in points:
        best = 0
        best_d = _squared_distance(point, centers[0])
        for c in range(1, len(centers)):
            d = _squared_distance(point, centers[c])
            if d < best_d:  # 严格小于 => 等距并列取编号最小簇
                best_d = d
                best = c
        labels.append(best)
    return labels


def _reseed_empty_clusters(
    points: List[List[float]],
    labels: List[int],
    centers: List[List[float]],
) -> None:
    """对每个空簇，确定性地重播种到“距已成形中心最远”的点。"""
    k = len(centers)
    populated = {labels[i] for i in range(len(points))}
    for empty in range(k):
        if empty in populated:
            continue
        # 距任一已成形中心的最近平方距离最大者（并列取下标最小）
        best_idx = -1
        best_d = -1.0
        for idx, point in enumerate(points):
            d_sq = float("inf")
            for c in populated:
                d = _squared_distance(point, centers[c])
                if d < d_sq:
                    d_sq = d
            if d_sq > best_d:
                best_d = d_sq
                best_idx = idx
        centers[empty] = list(points[best_idx])
        labels[best_idx] = empty
        populated.add(empty)


def kmeans(
    points: List[List[float]], k: int, seed: int
) -> KMeansResult:
    n = len(points)
    if k < 1 or k > n:
        raise ValueError(f"簇数 k={k} 对 {n} 个点无法成立")

    rng = MT19937(seed)
    centers = _kmeans_plus_plus(points, k, rng)
    n_dims = len(points[0]) if points and points[0] else 0

    labels = _assign(points, centers)
    iterations = 0
    for iteration in range(1, MAX_ITER + 1):
        iterations = iteration
        # 复制旧中心用于收敛判定
        old_centers = [list(center) for center in centers]

        member_count = [0] * k
        sums = [[0.0] * n_dims for _ in range(k)]
        for idx, point in enumerate(points):
            label = labels[idx]
            member_count[label] += 1
            target = sums[label]
            for d in range(n_dims):
                target[d] += point[d]

        empty = [c for c in range(k) if member_count[c] == 0]
        if empty:
            _reseed_empty_clusters(points, labels, centers)
            # 重播种改变了归属，重新分配后进入下一轮
            labels = _assign(points, centers)
            continue

        for c in range(k):
            inv = 1.0 / member_count[c]
            centers[c] = [s * inv for s in sums[c]]

        new_labels = _assign(points, centers)
        if new_labels == labels and centers == old_centers:
            labels = new_labels
            break
        labels = new_labels

    # 将簇编号规范为按“首个成员出现位置”排序，保证标签稳定可读
    remap: dict[int, int] = {}
    for label in labels:
        if label not in remap:
            remap[label] = len(remap)
    canonical_labels = [remap[label] for label in labels]
    canonical_centers = [centers[old] for old in sorted(remap, key=lambda x: remap[x])]
    return KMeansResult(
        labels=canonical_labels,
        centers=canonical_centers,
        iterations=iterations,
    )
