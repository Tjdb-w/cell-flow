"""可选双细胞识别：QC 候选细胞上的确定性 doublet 评分与过滤。

启用 ``--detect-doublets`` 后，以 QC 候选细胞 × 候选基因的原始计数子矩阵
计算每个候选细胞的 doublet_score：用确定性伪随机（由 ``--seed`` 派生）合成
“两个候选细胞计数叠加”的模拟双细胞，候选细胞与其最相似模拟双细胞之间的
余弦相似度均值即其评分——真实双细胞的表达谱更接近两个细胞的叠加谱，
单细胞谱则相反。同输入、同 seed、同参数下评分逐比特可复现。

目标标记数按候选细胞数与 ``--expected-doublet-rate`` 计算；分数降序、
同分按原列序依次标记；rate 为 0 不标记，且任何情况下至少保留一个细胞。
过滤后按 ``--min-cells`` 在最终细胞上重算保留基因，下游只用最终细胞与基因。
"""

import math
from dataclasses import dataclass
from typing import List

from .io import ExpressionMatrix
from .prng import MT19937
from .qc import QCResult

SCORE_BINS = 20
TOP_K_SYNTHETIC = 10
MAX_SYNTHETIC = 2000
NORMALIZATION_TARGET = 10000.0


@dataclass(frozen=True)
class DoubletResult:
    """双细胞评分与过滤结果。"""

    cell_ids: List[str]       # QC 候选细胞 ID（原列序）
    scores: List[float]       # 与 cell_ids 对齐的 doublet_score
    ranks: List[int]          # 与 cell_ids 对齐的名次（自 1 起，分数降序、同分按原列序）
    flags: List[bool]         # 与 cell_ids 对齐：是否被标记过滤
    kept_cells: List[int]     # 过滤后保留细胞在原矩阵中的列下标（原列序）
    kept_genes: List[int]     # 过滤后按 min_cells 重算的保留基因行下标（原行序）


def _log_profile(counts: List[int]) -> List[float]:
    """子矩阵内的文库大小归一化 + log1p；总计数为零时为零向量。"""
    total = sum(counts)
    if total <= 0:
        return [0.0] * len(counts)
    return [math.log1p(value / total * NORMALIZATION_TARGET) for value in counts]


def _cosine(a: List[float], b: List[float]) -> float:
    dot = 0.0
    norm_a = 0.0
    norm_b = 0.0
    for x, y in zip(a, b):
        dot += x * y
        norm_a += x * x
        norm_b += y * y
    if norm_a <= 0.0 or norm_b <= 0.0:
        return 0.0
    return dot / math.sqrt(norm_a * norm_b)


def detect_doublets(
    matrix: ExpressionMatrix,
    qc: QCResult,
    *,
    rate: float,
    min_cells: int,
    seed: int,
) -> DoubletResult:
    """对 QC 候选细胞评分并按 rate 标记双细胞。

    调用方保证候选细胞数 >= 2；候选基因即 ``qc.kept_genes``。
    """
    cell_idx = qc.kept_cells
    gene_idx = qc.kept_genes
    n = len(cell_idx)
    n_genes = len(gene_idx)

    # 候选细胞 × 候选基因的原始计数子矩阵（按列组织，列序即候选细胞原列序）
    columns = [
        [matrix.counts[g][c] for g in gene_idx] for c in cell_idx
    ]
    profiles = [_log_profile(column) for column in columns]

    # 确定性合成模拟双细胞：每次取两个不同候选细胞的原始计数叠加
    rng = MT19937(seed)
    n_synth = min(2 * n, MAX_SYNTHETIC)
    synthetics: List[List[float]] = []
    for _ in range(n_synth):
        i = rng.randbelow(n)
        j = rng.randbelow(n - 1)
        if j >= i:
            j += 1
        summed = [columns[i][g] + columns[j][g] for g in range(n_genes)]
        synthetics.append(_log_profile(summed))

    # 评分：与最相似 k 个模拟双细胞的余弦相似度均值
    k = min(TOP_K_SYNTHETIC, n_synth)
    scores: List[float] = []
    for profile in profiles:
        sims = [_cosine(profile, synth) for synth in synthetics]
        sims.sort(reverse=True)
        scores.append(sum(sims[:k]) / k)

    # 名次：分数降序，同分按原列序；rank 自 1 起
    order = sorted(range(n), key=lambda i: (-scores[i], i))
    ranks = [0] * n
    for position, i in enumerate(order, start=1):
        ranks[i] = position

    # 目标标记数按候选数与 rate 计算；rate 为 0 不标记，且至少保留一个细胞
    n_flag = int(rate * n)
    if n_flag > n - 1:
        n_flag = n - 1
    flags = [False] * n
    for i in order[:n_flag]:
        flags[i] = True

    kept_cells = [cell_idx[i] for i in range(n) if not flags[i]]

    # 过滤后在最终细胞上按 min_cells 重算保留基因（原行序）
    kept_genes: List[int] = []
    for g in range(matrix.n_genes):
        row = matrix.counts[g]
        detected = 0
        for c in kept_cells:
            if row[c] > 0:
                detected += 1
        if detected >= min_cells:
            kept_genes.append(g)

    return DoubletResult(
        cell_ids=[matrix.cell_ids[c] for c in cell_idx],
        scores=scores,
        ranks=ranks,
        flags=flags,
        kept_cells=kept_cells,
        kept_genes=kept_genes,
    )
