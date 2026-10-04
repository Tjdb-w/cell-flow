"""确定性双细胞（doublet）识别与过滤。

启用 ``--detect-doublets`` 后，以质控（QC）候选细胞的原始计数子矩阵
（全部输入基因）评分：

1. 用固定种子的 :class:`cell_flow.prng.MT19937` 抽取等量候选细胞两两配对，
   把两个细胞的计数谱相加，再按候选细胞平均文库深度做无放回降采样
   （超几何抽样，Fisher–Yates 部分洗牌实现），得到“两个细胞叠加”的合成谱；
2. 观测谱与合成谱都按各自文库大小归一到 10000 后取 ``ln(x + 1)``；
3. 每个观测细胞在全部基因的 log 归一化空间中取最近的
   ``min(15, 2*候选数 - 1)`` 个谱（欧氏距离平方，并列按原列序：
   观测细胞在前、合成谱在后），其中合成谱所占比例即 ``doublet_score``；
4. 按分数降序、同分按原列序给出 ``doublet_rank``（自 1 起）；
   目标标记数为 ``floor(候选数 * rate)``（至少保留一个细胞），前若干名
   ``doublet_flag`` 为真并被过滤。

同输入、同种子、同参数下结果确定，可跨运行逐字节复现。
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .io import ExpressionMatrix
from .prng import MT19937
from .qc import CellQC, GeneQC, QCResult

NORMALIZATION_TARGET = 10000.0
DOUBLET_NEIGHBORS = 15
SCORE_BINS = 20


@dataclass(frozen=True)
class DoubletResult:
    """QC 候选细胞上的双细胞评分、排名、标记与分数直方图。"""

    cell_ids: List[str]           # QC 候选细胞（原矩阵列序）
    scores: List[float]           # 与 cell_ids 对齐
    ranks: List[int]              # 自 1 起；分数降序、同分原列序
    flagged: List[bool]           # 与 cell_ids 对齐；是否被标记过滤
    retained: List[bool]          # 与 cell_ids 对齐；过滤后是否保留
    bins: List[Tuple[float, float, int]]  # 20 个等宽箱 (start, end, count)


def _log_normalize_row(row: List[int], total: int) -> List[float]:
    if total <= 0:
        return [0.0] * len(row)
    scale = NORMALIZATION_TARGET / total
    return [math.log1p(value * scale) for value in row]


def _downsample_counts(
    counts: List[int], draws: int, rng: MT19937
) -> List[int]:
    """从 ``counts`` 的 UMI 池中无放回抽取 ``draws`` 个（超几何）。

    Fisher–Yates 部分洗牌：位置 j 的标签默认按累积计数二分查找，被换过
    的位置从 ``swaps`` 取，复杂度 O(draws log G)。
    """
    population = sum(counts)
    if draws >= population:
        return list(counts)
    cumulative: List[int] = []
    acc = 0
    for value in counts:
        acc += value
        cumulative.append(acc)

    def gene_of(position: int) -> int:
        lo, hi = 0, len(cumulative)
        while lo < hi:
            mid = (lo + hi) // 2
            if cumulative[mid] <= position:
                lo = mid + 1
            else:
                hi = mid
        return lo

    result = [0] * len(counts)
    swaps: Dict[int, int] = {}
    for i in range(draws):
        j = i + rng.randbelow(population - i)
        label_j = swaps.pop(j, None)
        if label_j is None:
            label_j = gene_of(j)
        label_i = swaps.pop(i, None)
        if label_i is None:
            label_i = gene_of(i)
        swaps[j] = label_i
        result[label_j] += 1
    return result


def _squared_distance(a: List[float], b: List[float]) -> float:
    total = 0.0
    for x, y in zip(a, b):
        d = x - y
        total += d * d
    return total


def score_doublets(
    matrix: ExpressionMatrix,
    candidate_cells: List[int],
    expected_rate: float,
    seed: int,
) -> DoubletResult:
    """计算 QC 候选细胞的确定性双细胞评分并按目标比例标记。"""
    n = len(candidate_cells)
    if n < 2:
        raise ValueError("双细胞识别至少需要两个 QC 候选细胞")
    n_genes = matrix.n_genes
    raw = matrix.counts

    totals = [0] * n
    for ci, c in enumerate(candidate_cells):
        total = 0
        for g in range(n_genes):
            total += raw[g][c]
        totals[ci] = total

    observed_profiles: List[List[float]] = []
    for ci, c in enumerate(candidate_cells):
        observed_profiles.append(
            _log_normalize_row([raw[g][c] for g in range(n_genes)], totals[ci])
        )

    rng = MT19937(seed)
    mean_total = sum(totals) // n
    target_depth = mean_total if mean_total >= 1 else 1

    simulated_profiles: List[List[float]] = []
    for _ in range(n):
        i = rng.randbelow(n)
        j = rng.randbelow(n - 1)
        if j >= i:
            j += 1
        ci, cj = candidate_cells[i], candidate_cells[j]
        pair_row = [raw[g][ci] + raw[g][cj] for g in range(n_genes)]
        sampled = _downsample_counts(pair_row, target_depth, rng)
        simulated_profiles.append(_log_normalize_row(sampled, sum(sampled)))

    k_neighbors = min(DOUBLET_NEIGHBORS, 2 * n - 1)
    scores: List[float] = []
    for ci in range(n):
        observed = observed_profiles[ci]
        # 索引 < n：其他观测细胞（原列序）；>= n：合成谱（生成序）。
        # 并列按该索引排序即“同分按原列序、观测在前合成在后”。
        distances: List[Tuple[float, int]] = [
            (_squared_distance(observed, observed_profiles[other]), other)
            for other in range(n)
            if other != ci
        ]
        distances.extend(
            (_squared_distance(observed, simulated_profiles[s]), n + s)
            for s in range(n)
        )
        distances.sort(key=lambda item: (item[0], item[1]))
        n_simulated_neighbors = sum(
            1 for _, idx in distances[:k_neighbors] if idx >= n
        )
        scores.append(n_simulated_neighbors / k_neighbors)

    ranked = sorted(range(n), key=lambda c: (-scores[c], c))
    ranks = [0] * n
    for position, c in enumerate(ranked):
        ranks[c] = position + 1

    n_flagged = min(n - 1, math.floor(n * expected_rate))
    if n_flagged < 0:
        n_flagged = 0
    flagged_set = set(ranked[:n_flagged])
    flagged = [c in flagged_set for c in range(n)]
    retained = [not value for value in flagged]

    return DoubletResult(
        cell_ids=[matrix.cell_ids[c] for c in candidate_cells],
        scores=scores,
        ranks=ranks,
        flagged=flagged,
        retained=retained,
        bins=_score_bins(scores),
    )


def _score_bins(scores: List[float]) -> List[Tuple[float, float, int]]:
    """按分数分 20 个等宽箱；空箱保留，边界确定。"""
    lo = min(scores)
    hi = max(scores)
    width = (hi - lo) / SCORE_BINS
    counts = [0] * SCORE_BINS
    for value in scores:
        if width <= 0.0:
            idx = 0
        else:
            idx = int((value - lo) / width)
            if idx < 0:
                idx = 0
            elif idx >= SCORE_BINS:
                idx = SCORE_BINS - 1
        counts[idx] += 1
    bins: List[Tuple[float, float, int]] = []
    for b in range(SCORE_BINS):
        start = lo + b * width
        end = lo + (b + 1) * width if width > 0.0 else hi
        bins.append((start, end, counts[b]))
    return bins


def apply_doublet_filter(
    matrix: ExpressionMatrix,
    qc: QCResult,
    doublet: DoubletResult,
    min_cells: int,
) -> QCResult:
    """过滤被标记的 QC 候选细胞，并按 ``min_cells`` 在最终细胞上重算基因。

    返回新的 :class:`QCResult`：``kept_cells`` 为最终细胞（原列序），
    ``cell_qc`` 中被标记细胞的 ``retained`` 置为假；``gene_qc`` 的
    ``total_counts``/``detected_cells``/``retained`` 全部只按最终细胞统计。
    """
    candidate_set = set(qc.kept_cells)
    flagged_cell_ids = {
        cell_id
        for cell_id, is_flagged in zip(doublet.cell_ids, doublet.flagged)
        if is_flagged
    }
    final_cells = [
        c
        for c in qc.kept_cells
        if matrix.cell_ids[c] not in flagged_cell_ids
    ]

    cell_qc = [
        (
            CellQC(
                cell_id=record.cell_id,
                total_counts=record.total_counts,
                detected_genes=record.detected_genes,
                mitochondrial_counts=record.mitochondrial_counts,
                mitochondrial_fraction=record.mitochondrial_fraction,
                retained=False,
            )
            if c in candidate_set and record.cell_id in flagged_cell_ids
            else record
        )
        for c, record in enumerate(qc.cell_qc)
    ]

    gene_qc: List[GeneQC] = []
    kept_genes: List[int] = []
    for g, gene_id in enumerate(matrix.gene_ids):
        row = matrix.counts[g]
        gtotal = 0
        gdet = 0
        for c in final_cells:
            value = row[c]
            gtotal += value
            if value > 0:
                gdet += 1
        retained = gdet >= min_cells
        gene_qc.append(
            GeneQC(
                gene_id=gene_id,
                total_counts=gtotal,
                detected_cells=gdet,
                retained=retained,
            )
        )
        if retained:
            kept_genes.append(g)

    return QCResult(
        cell_qc=cell_qc,
        gene_qc=gene_qc,
        kept_cells=final_cells,
        kept_genes=kept_genes,
    )
