"""自动簇数选择：在 PCA 坐标上扫描 k=2..min(10, 细胞数)。

每个候选 k 复用与显式模式完全相同的确定性
:func:`cell_flow.kmeans.kmeans`（同 seed、同 k-means++ 实现），
以平均轮廓系数作为选择依据，并列取较小 k。
"""

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence

from .errors import CellFlowDataError
from .kmeans import KMeansResult, kmeans

AUTO_MAX_CLUSTERS = 10


@dataclass(frozen=True)
class KCandidate:
    """单个候选 k 的评估结果；无效行两项指标为 nan。"""

    k: int
    formed_clusters: int
    valid: bool
    iterations: int
    within_cluster_sse: float
    mean_silhouette: float
    clustering: Optional[KMeansResult]


@dataclass(frozen=True)
class SelectionResult:
    candidates: List[KCandidate]
    selected: int
    clustering: KMeansResult


def _squared_distance(a: Sequence[float], b: Sequence[float]) -> float:
    total = 0.0
    for x, y in zip(a, b):
        d = x - y
        total += d * d
    return total


def _pairwise_distances(points: Sequence[Sequence[float]]) -> List[List[float]]:
    """全部点对的欧氏距离矩阵（对称、对角为 0），各候选 k 共用。"""
    n = len(points)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        row = dist[i]
        for j in range(i + 1, n):
            d = math.sqrt(_squared_distance(points[i], points[j]))
            row[j] = d
            dist[j][i] = d
    return dist


def _within_cluster_sse(
    points: Sequence[Sequence[float]],
    labels: Sequence[int],
    centers: Sequence[Sequence[float]],
) -> float:
    total = 0.0
    for idx, point in enumerate(points):
        total += _squared_distance(point, centers[labels[idx]])
    return total


def _mean_silhouette(
    labels: Sequence[int], dist: Sequence[Sequence[float]], k: int
) -> float:
    """平均轮廓系数：单细胞簇记 0；其余 (b-a)/max(a,b)，b 取最近其他簇均距。"""
    n = len(labels)
    members: List[List[int]] = [[] for _ in range(k)]
    for idx, label in enumerate(labels):
        members[label].append(idx)

    total = 0.0
    for idx in range(n):
        own = labels[idx]
        own_members = members[own]
        if len(own_members) <= 1:
            # 单细胞簇：轮廓系数记 0，仍计入分母
            continue
        a = sum(dist[idx][j] for j in own_members if j != idx) / (
            len(own_members) - 1
        )
        b = math.inf
        for other in range(k):
            if other == own:
                continue
            other_members = members[other]
            mean_d = sum(dist[idx][j] for j in other_members) / len(other_members)
            if mean_d < b:
                b = mean_d
        if a == 0.0 and b == 0.0:
            # 与本簇及最近簇完全重合时按约定记 0
            continue
        total += (b - a) / (b if b > a else a)
    return total / n


def select_clustering(
    points: List[List[float]], seed: int
) -> SelectionResult:
    """扫描候选 k 并按平均轮廓系数选簇；无有效候选时抛 CellFlowDataError。"""
    n_cells = len(points)
    upper = min(AUTO_MAX_CLUSTERS, n_cells)

    # 点对距离只依赖 PCA 坐标，一次计算供所有候选复用
    dist = _pairwise_distances(points)

    candidates: List[KCandidate] = []
    best_score = math.nan
    best_k: Optional[int] = None
    best_clustering: Optional[KMeansResult] = None

    for k in range(2, upper + 1):
        clustering = kmeans(points, k, seed)
        formed_clusters = len(set(clustering.labels))
        valid = formed_clusters == k
        if valid:
            sse = _within_cluster_sse(
                points, clustering.labels, clustering.centers
            )
            score = _mean_silhouette(clustering.labels, dist, k)
        else:
            # 最终未形成恰好 k 个非空簇：该行无效，其余 k 照常评估
            sse = math.nan
            score = math.nan
        candidates.append(
            KCandidate(
                k=k,
                formed_clusters=formed_clusters,
                valid=valid,
                iterations=clustering.iterations,
                within_cluster_sse=sse,
                mean_silhouette=score,
                clustering=clustering if valid else None,
            )
        )
        # 严格大于才替换 => 轮廓系数并列时保留较小 k（k 按升序评估）
        if valid and (best_k is None or score > best_score):
            best_score = score
            best_k = k
            best_clustering = clustering

    if best_k is None or best_clustering is None:
        raise CellFlowDataError(
            f"候选簇数范围 2..{upper} 内无法形成两个以上不同簇"
            f"（质控后 {n_cells} 个细胞的方差不足以支撑聚类）"
        )
    return SelectionResult(
        candidates=candidates,
        selected=best_k,
        clustering=best_clustering,
    )
