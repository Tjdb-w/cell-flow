"""可选聚类稳定性分析：对质控后细胞做无放回子抽样，度量聚类可复现性。

仅当用户显式启用（``--stability``）时由管线调用，不改变默认分析路径。
完整数据的一次降维聚类结果作为参照；每次抽样沿用与主分析完全相同的
预处理（质控、log 归一化、HVG）、降维（PCA）与聚类（k-means，auto 模式下
为轮廓系数选 k）配置生成标签，只在**抽样细胞与参照标签的交集**（即本次
抽样命中的细胞）上计算调整兰德指数（Adjusted Rand Index）。

抽样为不放回抽样，由稳定性随机种子派生确定性的 :class:`MT19937` 驱动，
保证相同输入、配置与种子的逐次结果、summary 与图表数据逐字节一致；
该随机数流与主分析 ``--seed`` 相互独立，改变稳定性种子只改变抽样选择
与稳定性统计，不触及主分析的任何结果。

输出目录中已有的其他分析产物必须保留：稳定性文件以目录内暂存 +
``os.replace`` 覆盖改名方式写入既有结果目录，只覆盖同名稳定性文件。
"""

import json
import os
import tempfile
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .errors import OutputPathError
from .io import ExpressionMatrix
from .kmeans import kmeans
from .normalize import normalize_and_select_hvg
from .output import fmt_float
from .pca import run_pca
from .prng import MT19937
from .qc import QCResult
from .selection import select_cluster_count

SCORES_FILE = "stability_scores.csv"
SUMMARY_FILE = "stability_summary.json"
CHART_FILE = "stability_score_distribution.json"

HIST_BINS = 20


@dataclass(frozen=True)
class StabilityConfig:
    """聚类稳定性分析配置（仅在显式启用时存在）。"""

    n_resamples: int
    sample_fraction: float
    seed: int


@dataclass(frozen=True)
class StabilityResult:
    """逐次抽样稳定性结果与汇总；内容全部在内存中确定后再发布。"""

    n_resamples: int
    sample_fraction: float
    seed: int
    n_population: int
    sample_size: int
    scores: List[float]
    # 与 scores 对齐：每次抽样命中细胞在总体（质控保留）序列中的下标，升序
    samples: List[List[int]]

    @property
    def mean(self) -> float:
        return sum(self.scores) / len(self.scores)

    @property
    def median(self) -> float:
        return _median(self.scores)

    @property
    def minimum(self) -> float:
        return min(self.scores)

    @property
    def maximum(self) -> float:
        return max(self.scores)


def validate_stability_config(stability: StabilityConfig) -> None:
    """校验稳定性配置；任一非法一律抛 ``ValueError`` 并指出配置字段。"""
    if (
        not isinstance(stability.n_resamples, int)
        or isinstance(stability.n_resamples, bool)
        or stability.n_resamples < 2
    ):
        raise ValueError(
            "聚类稳定性配置字段 n_resamples（--stability-n-resamples）"
            "必须是 >= 2 的整数"
        )
    fraction = stability.sample_fraction
    if (
        not isinstance(fraction, (int, float))
        or isinstance(fraction, bool)
        or not 0.0 < float(fraction) < 1.0
    ):
        raise ValueError(
            "聚类稳定性配置字段 sample_fraction（--stability-sample-fraction）"
            "必须严格位于 0 与 1 之间（0 < sample_fraction < 1）"
        )
    if not isinstance(stability.seed, int) or isinstance(stability.seed, bool):
        raise ValueError(
            "聚类稳定性配置字段 seed（--stability-seed）必须是整数"
        )


def _median(values: Sequence[float]) -> float:
    """有限数值序列的中位数（偶数个取中间两个的算术平均）。"""
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def _sample_without_replacement(
    n_population: int, sample_size: int, rng: MT19937
) -> List[int]:
    """无放回均匀抽样：返回抽中下标的升序列表。

    使用部分 Fisher–Yates 洗牌（从末尾交换 ``sample_size`` 次），
    每次交换的目标位置由确定性 ``randbelow`` 给出。
    """
    indices = list(range(n_population))
    for i in range(n_population - 1, n_population - sample_size - 1, -1):
        j = rng.randbelow(i + 1)
        indices[i], indices[j] = indices[j], indices[i]
    return sorted(indices[n_population - sample_size:])


def _sub_qc(qc: QCResult, sampled: List[int]) -> QCResult:
    """以抽中的保留细胞构造子分析的 QC 视图。

    基因沿用参照分析的保留基因（``kept_genes`` 不变）；``kept_cells`` 取
    抽样命中细胞的**原矩阵列下标**，保持其在原列序中的相对顺序。
    ``cell_qc``/``gene_qc`` 仍为完整记录：归一化只经由
    ``kept_cells``/``kept_genes`` 访问数据，故与主分析语义一致。
    """
    return QCResult(
        cell_qc=qc.cell_qc,
        gene_qc=qc.gene_qc,
        kept_cells=[qc.kept_cells[s] for s in sampled],
        kept_genes=list(qc.kept_genes),
    )


def adjusted_rand_index(labels_a: Sequence[int], labels_b: Sequence[int]) -> float:
    """两个等长标签划分上的调整兰德指数（Hubert & Arabie, 1985）。

    基于列联表计数：``ARI = (Index - ExpectedIndex) / (MaxIndex -
    ExpectedIndex)``。分母为 0（任一划分只有一个簇，或少于两个对象）时
    约定返回 1.0——两划分在交集上没有可争议的结构。
    """
    n = len(labels_a)
    if n != len(labels_b):
        raise ValueError("ARI 要求两个标签序列等长")
    if n < 2:
        return 1.0

    contingency: Dict[Tuple[int, int], int] = {}
    row_totals: Dict[int, int] = {}
    col_totals: Dict[int, int] = {}
    for a, b in zip(labels_a, labels_b):
        key = (a, b)
        contingency[key] = contingency.get(key, 0) + 1
        row_totals[a] = row_totals.get(a, 0) + 1
        col_totals[b] = col_totals.get(b, 0) + 1

    def comb2(x: int) -> float:
        return x * (x - 1) / 2.0

    index = sum(comb2(v) for v in contingency.values())
    expected = (
        sum(comb2(v) for v in row_totals.values())
        * sum(comb2(v) for v in col_totals.values())
        / comb2(n)
    )
    max_index = 0.5 * (
        sum(comb2(v) for v in row_totals.values())
        + sum(comb2(v) for v in col_totals.values())
    )
    denominator = max_index - expected
    if denominator == 0.0:
        return 1.0
    return (index - expected) / denominator


def run_stability(
    matrix: ExpressionMatrix,
    qc: QCResult,
    reference_labels: Sequence[int],
    *,
    n_hvg: int,
    n_pcs: int,
    n_clusters: Any,
    seed: int,
    stability: StabilityConfig,
    batch_labels: Optional[List[str]] = None,
) -> StabilityResult:
    """执行稳定性重抽样分析。

    ``qc`` 是主分析实际使用的质控结果（启用双细胞过滤后为过滤后的细胞/
    基因集合）；``reference_labels`` 为完整数据一次降维聚类的参照标签，
    与 ``qc.kept_cells`` 同序（即主分析 ``clustering.labels``）。
    ``batch_labels`` 与 ``qc.kept_cells`` 对齐，为实际进入主分析降维聚类
    的批次标签（未实际施加校正时为 ``None``，与主分析传给归一化的取值
    一致）。``n_clusters`` 为整数或 ``"auto"``，语义与主分析完全相同；
    ``seed`` 为主分析 ``--seed``，每次聚类沿用同一配置。

    抽样随机数仅由 ``stability.seed`` 派生。
    """
    validate_stability_config(stability)

    n_population = len(qc.kept_cells)
    if matrix.n_cells == 0 or matrix.n_genes == 0:
        raise ValueError("输入矩阵为空：没有任何细胞或特征，可用细胞数为 0")
    if n_population < 2:
        raise ValueError(
            f"通过质量控制后仅有 {n_population} 个细胞，不足两个"
            f"（可用细胞数：{n_population}）"
        )
    if len(reference_labels) != n_population:
        raise ValueError("参照标签数必须与质控后细胞数一致")

    sample_size = int(stability.sample_fraction * n_population)
    if sample_size < 2:
        raise ValueError(
            f"抽样比例 {stability.sample_fraction} 在 {n_population} 个质控后"
            f"细胞上仅得到 {sample_size} 个抽样细胞，不足两个，无法聚类"
        )
    if n_clusters != "auto" and sample_size < n_clusters:
        raise ValueError(
            f"抽样比例 {stability.sample_fraction} 在 {n_population} 个质控后"
            f"细胞上仅得到 {sample_size} 个抽样细胞，少于请求簇数 "
            f"{n_clusters}，聚类无法成立"
        )

    rng = MT19937(stability.seed)
    scores: List[float] = []
    samples: List[List[int]] = []
    for _ in range(stability.n_resamples):
        sampled = _sample_without_replacement(n_population, sample_size, rng)
        sampled_labels = _cluster_subset(
            matrix,
            qc,
            sampled=sampled,
            n_hvg=n_hvg,
            n_pcs=n_pcs,
            n_clusters=n_clusters,
            seed=seed,
            batch_labels=_batch_subset(batch_labels, sampled),
        )
        # 交集即本次抽样命中的细胞（参照标签在全部质控细胞上均有定义）
        ref_subset = [reference_labels[s] for s in sampled]
        scores.append(adjusted_rand_index(ref_subset, sampled_labels))
        samples.append(sampled)

    return StabilityResult(
        n_resamples=stability.n_resamples,
        sample_fraction=stability.sample_fraction,
        seed=stability.seed,
        n_population=n_population,
        sample_size=sample_size,
        scores=scores,
        samples=samples,
    )


def _batch_subset(
    batch_labels: Optional[List[str]], sampled: List[int]
) -> Optional[List[str]]:
    if batch_labels is None:
        return None
    return [batch_labels[s] for s in sampled]


def _cluster_subset(
    matrix: ExpressionMatrix,
    qc: QCResult,
    *,
    sampled: List[int],
    n_hvg: int,
    n_pcs: int,
    n_clusters: Any,
    seed: int,
    batch_labels: Optional[List[str]],
) -> List[int]:
    """在一次抽样细胞上沿用主分析语义生成聚类标签。

    返回与 ``sampled`` 同序的标签。归一化、HVG、PCA、k-means/auto 选 k
    全部复用主分析的函数与参数；子样本上独立重算归一化与高变基因，即
    “每次抽样沿用同一预处理配置生成标签”。
    """
    sub = _sub_qc(qc, sampled)
    normalized = normalize_and_select_hvg(
        matrix, sub, n_hvg=n_hvg, batch_labels=batch_labels
    )
    if not normalized.selected_genes:
        raise ValueError("抽样细胞上的高变基因选择结果为空，聚类无法成立")
    pca = run_pca(normalized, n_pcs)

    if n_clusters == "auto":
        selection = select_cluster_count(pca.scores, seed)
        if selection is None:
            raise ValueError(
                "抽样细胞上的簇数自动选择失败：候选范围内无法形成两个以上"
                "不同簇（方差不足）"
            )
        return list(selection.clustering.labels)

    if n_clusters > len(pca.cell_ids):
        raise ValueError(
            f"簇数 {n_clusters} 大于本次抽样细胞数 {len(pca.cell_ids)}，"
            f"聚类无法成立"
        )
    clustering = kmeans(pca.scores, n_clusters, seed)
    if len(set(clustering.labels)) < n_clusters:
        raise ValueError(
            f"抽样细胞上的聚类无法成立：请求 {n_clusters} 个簇，但数据仅能"
            f"支撑 {len(set(clustering.labels))} 个不同簇（方差不足）"
        )
    return list(clustering.labels)


def _scores_csv(result: StabilityResult) -> str:
    """逐次稳定性分数 CSV：抽样序号、本次抽样细胞数、调整兰德指数。"""
    lines = ["resample,n_sampled_cells,adjusted_rand_index"]
    for i, score in enumerate(result.scores):
        lines.append(
            ",".join(
                [str(i + 1), str(result.sample_size), fmt_float(score)]
            )
        )
    return "\n".join(lines) + "\n"


def _summary_json(result: StabilityResult) -> str:
    """summary.json：抽样次数、抽样比例、随机种子及均值/中位数/最小/最大。"""
    payload = {
        "n_resamples": result.n_resamples,
        "sample_fraction": result.sample_fraction,
        "seed": result.seed,
        "mean": result.mean,
        "median": result.median,
        "min": result.minimum,
        "max": result.maximum,
    }
    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"


def _chart_json(result: StabilityResult) -> str:
    """供绘图使用的分数分布数据：逐次分数与等宽直方图计数。

    箱体边界由分数最小/最大值确定，等宽分为 20 箱（宽度为 0 时全部计入
    第 0 箱）；空箱保留。结构确定且自描述。
    """
    scores = list(result.scores)
    lo = min(scores)
    hi = max(scores)
    width = (hi - lo) / HIST_BINS
    counts = [0] * HIST_BINS
    for value in scores:
        if width <= 0.0:
            idx = 0
        else:
            idx = int((value - lo) / width)
            if idx < 0:
                idx = 0
            elif idx >= HIST_BINS:
                idx = HIST_BINS - 1
        counts[idx] += 1
    bins = []
    for b in range(HIST_BINS):
        start = lo + b * width
        end = lo + (b + 1) * width if width > 0.0 else hi
        bins.append(
            {
                "bin_index": b,
                "bin_start": start,
                "bin_end": end,
                "count": counts[b],
            }
        )
    payload = {
        "n_resamples": result.n_resamples,
        "sample_fraction": result.sample_fraction,
        "seed": result.seed,
        "n_cells_after_qc": result.n_population,
        "sample_size": result.sample_size,
        "scores": scores,
        "histogram": {
            "n_bins": HIST_BINS,
            "bin_min": lo,
            "bin_max": hi,
            "bins": bins,
        },
    }
    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"


def publish_stability(
    output_dir: str, result: StabilityResult
) -> List[str]:
    """将三个稳定性文件写入**已存在**的结果目录，保留其余产物。

    与主分析的事务性整目录发布不同：主分析发布后结果目录非空且其中
    产物必须保留，故逐文件先写入目录内隐藏临时文件再 ``os.replace``
    原子覆盖同名文件；任何失败都不删除、不改动既有产物（同名稳定性
    文件仅在替换成功的瞬间切换，再次运行即覆盖上一次的稳定性文件）。
    """
    if not os.path.isdir(output_dir) or os.path.islink(output_dir):
        raise OutputPathError(f"结果目录不存在或不是普通目录：{output_dir}")

    payload: List[Tuple[str, str]] = [
        (SCORES_FILE, _scores_csv(result)),
        (SUMMARY_FILE, _summary_json(result)),
        (CHART_FILE, _chart_json(result)),
    ]

    written: List[str] = []
    for name, content in payload:
        final_path = os.path.join(output_dir, name)
        fd, tmp_path = tempfile.mkstemp(
            prefix=f".{name}.", suffix=".tmp", dir=output_dir
        )
        try:
            try:
                with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
                    handle.write(content)
                    handle.flush()
                    os.fsync(handle.fileno())
            except OSError as exc:
                raise OutputPathError(
                    f"稳定性结果文件写出失败：{name}（{exc}）"
                ) from exc
            try:
                os.replace(tmp_path, final_path)
            except OSError as exc:
                raise OutputPathError(
                    f"稳定性结果文件发布失败：{name}（{exc}）"
                ) from exc
            tmp_path = None  # 已发布
            written.append(name)
        finally:
            if tmp_path is not None:
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass

    return written
