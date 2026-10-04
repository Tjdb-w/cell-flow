"""可选聚类稳定性分析：判断现有聚类对细胞抽样扰动是否稳定。

仅当用户显式启用（``--stability-analysis``）时执行，不改变默认分析路径。
以完整数据的一次降维聚类结果为参照，对通过质控（含可选双细胞过滤）后的
细胞做 ``n_samples`` 次不放回随机抽样，每次沿用同一预处理、降维与聚类
配置在抽样子矩阵上重新质控、归一化、选高变基因、PCA、聚类，只在抽样
细胞与参照标签的交集上计算调整兰德指数（Adjusted Rand Index）。

抽样由独立的稳定性种子派生确定性伪随机流
（:class:`cell_flow.prng.MT19937`），与主分析的 ``--seed`` 完全隔离：
改变稳定性种子只改变抽样选择与稳定性统计，绝不改变完整数据主分析的
降维坐标、聚类标签或差异表达结果。同输入、配置与种子重复执行逐次结果、
分数与汇总完全一致。
"""

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .errors import CellFlowDataError, CellFlowStabilityConfigError
from .io import ExpressionMatrix
from .kmeans import kmeans
from .normalize import normalize_and_select_hvg
from .pca import run_pca
from .prng import MT19937
from .qc import QCResult, compute_qc
from .selection import select_cluster_count

DEFAULT_STABILITY_SAMPLES = 100
DEFAULT_STABILITY_FRACTION = 0.8
DEFAULT_STABILITY_SEED = 20240617

FIELD_N_SAMPLES = "--stability-n-samples"
FIELD_SAMPLE_FRACTION = "--stability-sample-fraction"
FIELD_SEED = "--stability-seed"

STABILITY_FILE_SCORES = "stability_scores.csv"
STABILITY_FILE_SUMMARY = "summary.json"
STABILITY_FILE_DISTRIBUTION = "stability_score_distribution.json"

DISTRIBUTION_BINS = 20


def summarize_scores(scores: Sequence[float]) -> Dict[str, float]:
    """均值、中位数、最小值、最大值；偶数个时中位数取中间两者的平均。"""
    ordered = sorted(scores)
    mid = len(ordered) // 2
    median = (
        ordered[mid]
        if len(ordered) % 2 == 1
        else (ordered[mid - 1] + ordered[mid]) / 2.0
    )
    return {
        "mean": sum(scores) / len(scores),
        "median": median,
        "min": ordered[0],
        "max": ordered[-1],
    }


@dataclass(frozen=True)
class StabilityConfig:
    """稳定性分析参数。"""

    n_samples: int
    sample_fraction: float
    seed: int


@dataclass(frozen=True)
class StabilityResult:
    """逐次抽样的稳定性结果。

    ``sampled_cell_ids``、``sample_sizes``、``intersection_sizes`` 与
    ``scores`` 按抽样次序对齐；``reference_labels`` 与完整数据保留细胞
    （``analysis_qc.kept_cells`` 列序）对齐，即主分析 ``clusters.tsv`` 的行序。
    """

    config: StabilityConfig
    reference_cell_ids: List[str]
    reference_labels: List[int]
    n_available_cells: int
    scores: List[float]
    sampled_cell_ids: List[List[str]]
    sample_sizes: List[int]
    intersection_sizes: List[int]

    @property
    def n_samples(self) -> int:
        return self.config.n_samples


def validate_stability_settings(
    *,
    n_samples: Any,
    sample_fraction: Any,
    seed: Any,
) -> StabilityConfig:
    """校验稳定性分析配置；非法一律抛 :class:`ValueError` 并指出对应配置字段。

    实际抛出 :class:`cell_flow.errors.CellFlowStabilityConfigError`——它既是
    :class:`ValueError`（满足库调用方与需求的异常类型约定），又是命令行的
    配置错误（退出码 3）；消息中包含对应配置字段名。
    """
    if (
        not isinstance(n_samples, int)
        or isinstance(n_samples, bool)
        or n_samples < 2
    ):
        raise CellFlowStabilityConfigError(
            f"{FIELD_N_SAMPLES} 必须是 >= 2 的整数，得到 {n_samples!r}"
        )
    if (
        not isinstance(sample_fraction, (int, float))
        or isinstance(sample_fraction, bool)
        or not math.isfinite(float(sample_fraction))
        or not 0.0 < float(sample_fraction) < 1.0
    ):
        raise CellFlowStabilityConfigError(
            f"{FIELD_SAMPLE_FRACTION} 必须是 (0, 1) 开区间内的有限数值，"
            f"得到 {sample_fraction!r}"
        )
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise CellFlowStabilityConfigError(
            f"{FIELD_SEED} 必须是整数，得到 {seed!r}"
        )
    return StabilityConfig(
        n_samples=n_samples,
        sample_fraction=float(sample_fraction),
        seed=seed,
    )


def _sample_without_replacement(
    n_population: int, sample_size: int, rng: MT19937
) -> List[int]:
    """``[0, n_population)`` 上等概率不放回抽 ``sample_size`` 个总体下标。

    Fisher–Yates 部分洗牌，随机数全部来自确定性的 ``rng``，故同种子下
    抽样选择完全一致；返回结果按下标升序，便于与原矩阵列序对齐。
    """
    if not 0 <= sample_size <= n_population:
        raise ValueError("抽样规模超出总体范围")
    chosen = list(range(n_population))
    for i in range(sample_size):
        j = i + rng.randbelow(n_population - i)
        chosen[i], chosen[j] = chosen[j], chosen[i]
    return sorted(chosen[:sample_size])


def adjusted_rand_index(labels_a: Sequence[int], labels_b: Sequence[int]) -> float:
    """两个等长划分的调整兰德指数。

    采用列联表计数的标准定义
    ``ARI = (Index - ExpectedIndex) / (MaxIndex - ExpectedIndex)``；
    当分母为 0（两份划分各自至多一个簇，或总对数为 0）时约定为 1.0。
    """
    n = len(labels_a)
    if n != len(labels_b):
        raise ValueError("两份标签长度必须一致")
    if n < 2:
        # 交集中不足两个点时没有可比较的点对；视为完全一致
        return 1.0

    def comb2(x: int) -> int:
        return x * (x - 1) // 2

    contingency: Dict[Tuple[int, int], int] = {}
    counts_a: Dict[int, int] = {}
    counts_b: Dict[int, int] = {}
    for a, b in zip(labels_a, labels_b):
        contingency[(a, b)] = contingency.get((a, b), 0) + 1
        counts_a[a] = counts_a.get(a, 0) + 1
        counts_b[b] = counts_b.get(b, 0) + 1

    sum_comb = sum(comb2(v) for v in contingency.values())
    sum_a = sum(comb2(v) for v in counts_a.values())
    sum_b = sum(comb2(v) for v in counts_b.values())
    total = comb2(n)

    expected = sum_a * sum_b / total if total > 0 else 0.0
    maximum = 0.5 * (sum_a + sum_b)
    denominator = maximum - expected
    if denominator == 0.0:
        return 1.0
    return (sum_comb - expected) / denominator


def _relabel(values: Sequence[int]) -> List[int]:
    """把整数标签压缩成 0..k-1（按首次出现顺序）。"""
    remap: Dict[int, int] = {}
    out: List[int] = []
    for value in values:
        if value not in remap:
            remap[value] = len(remap)
        out.append(remap[value])
    return out


def _submatrix(
    matrix: ExpressionMatrix, sampled_original: Sequence[int]
) -> ExpressionMatrix:
    """构造 全部基因 × 抽中细胞列 的子矩阵（细胞按原矩阵列序升序）。

    保留全部输入基因，使子矩阵上的细胞文库总量、检出基因数与线粒体计数
    与主分析完全一致；基因是否保留则由同一 ``min_cells`` 阈值在子矩阵上
    重新判定，沿用既有质控语义。
    """
    cell_ids = [matrix.cell_ids[c] for c in sampled_original]
    counts = [[row[c] for c in sampled_original] for row in matrix.counts]
    return ExpressionMatrix(
        gene_ids=list(matrix.gene_ids),
        cell_ids=cell_ids,
        counts=counts,
        sha256="",
        path="",
        total_counts=0,
    )


def _cluster_subset(
    matrix: ExpressionMatrix,
    subset_qc: QCResult,
    *,
    n_hvg: int,
    n_pcs: int,
    n_clusters: Any,
    seed: int,
    batch_labels: Optional[List[str]],
) -> List[int]:
    """在一个抽样子集上沿用主分析配置完成归一化/HVG/PCA/聚类。

    返回与 ``subset_qc.kept_cells`` 列序对齐的簇标签。任何无法成立的
    情形（高变基因为空、PCA 为 0 维、簇数多于细胞、簇数不足）都抛
    :class:`ValueError`。
    """
    normalized = normalize_and_select_hvg(
        matrix, subset_qc, n_hvg=n_hvg, batch_labels=batch_labels
    )
    if not normalized.selected_genes:
        raise ValueError("抽样子集高变基因选择结果为空，PCA 无法成立")

    pca = run_pca(normalized, n_pcs)
    n_cells = len(normalized.cell_ids)

    if n_clusters == "auto":
        selection = select_cluster_count(pca.scores, seed)
        if selection is None:
            raise ValueError("抽样子集簇数自动选择失败（方差不足）")
        return _relabel(selection.clustering.labels)

    if n_clusters > n_cells:
        raise ValueError(
            f"抽样子集细胞数 {n_cells} 小于请求簇数 {n_clusters}，聚类无法成立"
        )
    clustering = kmeans(pca.scores, n_clusters, seed)
    if len(set(clustering.labels)) < n_clusters:
        raise ValueError("抽样子集聚类无法形成请求的簇数（方差不足）")
    return _relabel(clustering.labels)


def run_stability_analysis(
    matrix: ExpressionMatrix,
    analysis_qc: QCResult,
    *,
    config: StabilityConfig,
    min_genes: int,
    max_mito_fraction: float,
    min_cells: int,
    mito_prefix: str,
    n_hvg: int,
    n_pcs: int,
    n_clusters: Any,
    seed: int,
    batch_labels: Optional[Sequence[str]] = None,
    reference_labels: Optional[Sequence[int]] = None,
) -> StabilityResult:
    """执行逐次抽样子集重聚类并计算调整兰德指数。

    ``analysis_qc`` 是主分析实际使用的 QC 结果（启用双细胞识别时为过滤后
    结果）；其 ``kept_cells`` 即抽样总体，``reference_labels`` 必须与之
    列序对齐。``batch_labels`` 同样与该列序对齐，提供时各抽样子集沿用
    主分析的批次均值中心化语义。``seed`` 为主分析聚类种子；抽样随机流
    只来自 ``config.seed``，两者相互独立。
    """
    available = list(analysis_qc.kept_cells)
    n_available = len(available)
    if n_available < 2:
        raise CellFlowDataError(
            f"通过质量控制后仅 {n_available} 个细胞可用，不足两个，"
            f"无法进行聚类稳定性分析"
        )
    if reference_labels is None or len(reference_labels) != n_available:
        raise ValueError("参照标签缺失或与可用细胞数不一致")

    sample_size = int(math.floor(config.sample_fraction * n_available + 1e-9))
    if sample_size < 2:
        raise CellFlowDataError(
            f"{FIELD_SAMPLE_FRACTION}={config.sample_fraction!r} 在 "
            f"{n_available} 个可用细胞上每次仅抽得 {sample_size} 个细胞，"
            f"不足两个，无法重新聚类"
        )
    if sample_size >= n_available:
        raise CellFlowDataError(
            f"{FIELD_SAMPLE_FRACTION}={config.sample_fraction!r} 在 "
            f"{n_available} 个可用细胞上每次抽得 {sample_size} 个细胞，"
            f"未形成真子集，无法进行稳定性抽样"
        )
    if isinstance(n_clusters, int) and not isinstance(n_clusters, bool):
        if sample_size < n_clusters:
            raise CellFlowDataError(
                f"{FIELD_SAMPLE_FRACTION}={config.sample_fraction!r} 在 "
                f"{n_available} 个可用细胞上每次仅抽得 {sample_size} 个细胞，"
                f"少于固定簇数 {n_clusters}，抽样子集无法聚类"
            )

    reference_cell_ids = [matrix.cell_ids[c] for c in available]
    reference_by_id = dict(zip(reference_cell_ids, reference_labels))
    batch_by_position = (
        {p: batch_labels[p] for p in range(n_available)}
        if batch_labels is not None
        else None
    )

    rng = MT19937(config.seed)
    scores: List[float] = []
    sampled_cell_ids: List[List[str]] = []
    sample_sizes: List[int] = []
    intersection_sizes: List[int] = []

    for _ in range(config.n_samples):
        picked = _sample_without_replacement(n_available, sample_size, rng)
        sampled_original = [available[p] for p in picked]

        # 在 全部基因 × 抽中细胞 的子矩阵上沿用同一 QC/预处理配置：
        # 细胞级指标是细胞固有量，抽中的保留细胞必然再次通过；
        # 基因保留按同一 min_cells 在子矩阵上重算。
        submatrix = _submatrix(matrix, sampled_original)
        subset_qc = compute_qc(
            submatrix,
            min_genes=min_genes,
            max_mito_fraction=max_mito_fraction,
            min_cells=min_cells,
            mito_prefix=mito_prefix,
        )
        if len(subset_qc.kept_cells) < 2 or not subset_qc.kept_genes:
            raise ValueError(
                "抽样子集重新质控后可用细胞或基因为空，无法重新聚类"
            )
        subset_batch_labels = (
            [batch_by_position[p] for p in picked] if batch_by_position else None
        )
        subset_labels = _cluster_subset(
            submatrix,
            subset_qc,
            n_hvg=n_hvg,
            n_pcs=n_pcs,
            n_clusters=n_clusters,
            seed=seed,
            batch_labels=subset_batch_labels,
        )

        # 只在抽样细胞与参照标签都存在的交集上对齐两份划分
        kept_subset_ids = [submatrix.cell_ids[c] for c in subset_qc.kept_cells]
        ref_part: List[int] = []
        sub_part: List[int] = []
        for cell_id, label in zip(kept_subset_ids, subset_labels):
            if cell_id in reference_by_id:
                ref_part.append(reference_by_id[cell_id])
                sub_part.append(label)
        intersection = len(ref_part)
        if intersection < 2:
            raise ValueError(
                "抽样子集与参照标签的交集不足两个细胞，无法计算稳定性分数"
            )
        score = adjusted_rand_index(ref_part, _relabel(sub_part))

        sampled_cell_ids.append([reference_cell_ids[p] for p in picked])
        sample_sizes.append(sample_size)
        intersection_sizes.append(intersection)
        scores.append(score)

    return StabilityResult(
        config=config,
        reference_cell_ids=reference_cell_ids,
        reference_labels=list(reference_labels),
        n_available_cells=n_available,
        scores=scores,
        sampled_cell_ids=sampled_cell_ids,
        sample_sizes=sample_sizes,
        intersection_sizes=intersection_sizes,
    )
