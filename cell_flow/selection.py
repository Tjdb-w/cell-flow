"""``--n-clusters auto`` 的簇数自动选择。

在质控后细胞的 PCA 坐标上，对候选 k=2..min(10, 细胞数) 逐个跑与显式
模式完全相同的确定性 KMeans（:func:`cell_flow.kmeans.kmeans`，同一 seed），
按平均轮廓系数选择：

- 仅当某次划分恰好形成 k 个非空簇时该行 valid=true，selected 才有资格；
- 单细胞簇中的细胞轮廓系数记 0；
- 其余细胞按到本簇与其他簇的平均欧氏距离计算 ``(b-a)/max(a,b)``；
- 并列（平均轮廓系数相同）取较小 k；
- 无效行的 iterations 与 formed_clusters 仍写实际值，两项指标写 nan。
"""

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence

from .kmeans import KMeansResult, kmeans

MAX_AUTO_K = 10


@dataclass(frozen=True)
class KSelectionRow:
    k: int
    formed_clusters: int
    valid: bool
    iterations: int
    within_cluster_sse: float   # 无效行为 nan
    mean_silhouette: float      # 无效行为 nan
    selected: bool


@dataclass(frozen=True)
class ClusterSelectionResult:
    """全部候选行（按 k 升序）、选中的 k 与对应聚类结果。"""

    rows: List[KSelectionRow]
    selected_k: int
    clustering: KMeansResult


def _squared_distance(a: Sequence[float], b: Sequence[float]) -> float:
    total = 0.0
    for x, y in zip(a, b):
        d = x - y
        total += d * d
    return total


def _within_cluster_sse(
    points: Sequence[Sequence[float]], clustering: KMeansResult
) -> float:
    centers = clustering.centers
    total = 0.0
    for idx, label in enumerate(clustering.labels):
        total += _squared_distance(points[idx], centers[label])
    return total


def _distance_table(points: Sequence[Sequence[float]]) -> List[List[float]]:
    """全部点对欧氏距离（对称、对角线为 0）；各候选 k 共享同一 PCA 坐标。"""
    n = len(points)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        point_i = points[i]
        for j in range(i + 1, n):
            d = math.sqrt(_squared_distance(point_i, points[j]))
            dist[i][j] = d
            dist[j][i] = d
    return dist


def _mean_silhouette(
    dist: Sequence[Sequence[float]], labels: Sequence[int]
) -> float:
    n = len(labels)

    members: dict[int, List[int]] = {}
    for idx, label in enumerate(labels):
        members.setdefault(label, []).append(idx)
    other_labels = sorted(members)

    total = 0.0
    for i in range(n):
        own = labels[i]
        own_members = members[own]
        if len(own_members) <= 1:
            # 单细胞簇：轮廓系数记 0
            total += 0.0
            continue
        a = sum(dist[i][j] for j in own_members if j != i) / (
            len(own_members) - 1
        )
        b = math.inf
        for other in other_labels:
            if other == own:
                continue
            other_members = members[other]
            mean_d = sum(dist[i][j] for j in other_members) / len(other_members)
            if mean_d < b:  # 严格小于 => 并列取编号最小簇，结果确定
                b = mean_d
        denom = max(a, b)
        total += (b - a) / denom if denom > 0.0 else 0.0
    return total / n


def select_cluster_count(
    points: Sequence[Sequence[float]], seed: int
) -> Optional[ClusterSelectionResult]:
    """评估全部候选 k 并返回选择结果；无有效候选时返回 ``None``。"""
    n = len(points)
    upper = min(MAX_AUTO_K, n)
    # 所有候选在同一 PCA 坐标上评估，点对距离只计算一次
    dist = _distance_table(points)

    rows: List[KSelectionRow] = []
    best_silhouette = math.inf
    best_k: Optional[int] = None
    best_clustering: Optional[KMeansResult] = None
    have_best = False

    for k in range(2, upper + 1):
        clustering = kmeans(points, k, seed)
        formed_clusters = len(set(clustering.labels))
        valid = formed_clusters == k
        if valid:
            sse = _within_cluster_sse(points, clustering)
            silhouette = _mean_silhouette(dist, clustering.labels)
        else:
            sse = float("nan")
            silhouette = float("nan")

        rows.append(
            KSelectionRow(
                k=k,
                formed_clusters=formed_clusters,
                valid=valid,
                iterations=clustering.iterations,
                within_cluster_sse=sse,
                mean_silhouette=silhouette,
                selected=False,
            )
        )

        # k 升序遍历 + 严格大于 => 平均轮廓系数并列时保留较小 k
        if valid and (not have_best or silhouette > best_silhouette):
            best_silhouette = silhouette
            best_k = k
            best_clustering = clustering
            have_best = True

    if not have_best or best_k is None or best_clustering is None:
        return None

    final_rows = [
        KSelectionRow(
            k=row.k,
            formed_clusters=row.formed_clusters,
            valid=row.valid,
            iterations=row.iterations,
            within_cluster_sse=row.within_cluster_sse,
            mean_silhouette=row.mean_silhouette,
            selected=row.valid and row.k == best_k,
        )
        for row in rows
    ]
    return ClusterSelectionResult(
        rows=final_rows,
        selected_k=best_k,
        clustering=best_clustering,
    )
