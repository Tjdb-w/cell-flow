"""分析管线编排：参数校验 -> 读入 -> QC -> 归一化/HVG -> PCA -> 聚类 -> 差异表达 -> 写出。"""

from dataclasses import dataclass
from typing import Any, Dict, List, Union

from . import __version__
from .errors import CellFlowConfigError, CellFlowDataError
from .io import ExpressionMatrix, read_matrix
from .kmeans import kmeans
from .markers import find_markers, find_pairwise_markers
from .normalize import normalize_and_select_hvg
from .output import Artifacts, ensure_output_dir, write_results
from .pca import MAX_PCS, run_pca
from .qc import compute_qc
from .selection import select_cluster_count

DEFAULT_SEED = 20240617
AUTO = "auto"


@dataclass(frozen=True)
class Config:
    input_path: str
    output_dir: str
    min_genes: int = 200
    max_mito_fraction: float = 0.2
    min_cells: int = 3
    mito_prefix: str = "MT-"
    n_hvg: int = 2000
    n_pcs: int = MAX_PCS
    n_clusters: Union[int, str] = 2
    seed: int = DEFAULT_SEED

    @property
    def auto_clusters(self) -> bool:
        return self.n_clusters == AUTO

    def public_parameters(self) -> Dict[str, Any]:
        return {
            "min_genes": self.min_genes,
            "max_mito_fraction": self.max_mito_fraction,
            "min_cells": self.min_cells,
            "mito_prefix": self.mito_prefix,
            "n_hvg": self.n_hvg,
            "n_pcs": self.n_pcs,
            "n_clusters": self.n_clusters,
            "seed": self.seed,
        }


def validate_config(config: Config) -> None:
    """校验参数取值；冲突或非法一律抛 CellFlowConfigError。"""
    errors: List[str] = []
    if not isinstance(config.min_genes, int) or config.min_genes < 0:
        errors.append("--min-genes 必须是非负整数")
    if not isinstance(config.min_cells, int) or config.min_cells < 1:
        errors.append("--min-cells 必须是 >= 1 的整数")
    if (
        not isinstance(config.max_mito_fraction, (int, float))
        or isinstance(config.max_mito_fraction, bool)
        or not 0.0 <= float(config.max_mito_fraction) <= 1.0
    ):
        errors.append("--max-mito-fraction 必须是 [0, 1] 内的数值")
    if not isinstance(config.mito_prefix, str) or config.mito_prefix == "":
        errors.append("--mito-prefix 必须是非空字符串")
    if not isinstance(config.n_hvg, int) or config.n_hvg < 1:
        errors.append("--n-hvg 必须是 >= 1 的整数")
    if not isinstance(config.n_pcs, int) or config.n_pcs < 1:
        errors.append("--n-pcs 必须是 >= 1 的整数")
    if config.n_clusters != AUTO:
        if (
            not isinstance(config.n_clusters, int)
            or isinstance(config.n_clusters, bool)
            or config.n_clusters < 2
        ):
            errors.append("--n-clusters 必须是 >= 2 的整数")
    if not isinstance(config.seed, int) or isinstance(config.seed, bool):
        errors.append("--seed 必须是整数")
    if errors:
        raise CellFlowConfigError("；".join(errors))


def run(config: Config) -> List[str]:
    validate_config(config)

    # 先完成读入与全部计算，最后再触碰输出目录：
    # 任何分析阶段失败都不应留下空目录或半成品结果
    matrix = read_matrix(config.input_path)

    qc = compute_qc(
        matrix,
        min_genes=config.min_genes,
        max_mito_fraction=config.max_mito_fraction,
        min_cells=config.min_cells,
        mito_prefix=config.mito_prefix,
    )

    n_kept_cells = len(qc.kept_cells)
    if n_kept_cells < 2:
        raise CellFlowDataError(
            f"质控后仅保留 {n_kept_cells} 个细胞，不足两个，无法分析"
        )
    if not qc.kept_genes:
        raise CellFlowDataError(
            f"没有基因在至少 {config.min_cells} 个细胞中检出，无可用基因"
        )

    if config.auto_clusters:
        # 自动模式：候选 k=2..min(10, 质控后细胞数)，PCA 后逐个评估
        normalized = normalize_and_select_hvg(matrix, qc, n_hvg=config.n_hvg)
        if not normalized.selected_genes:
            raise CellFlowDataError("高变基因选择结果为空，PCA 无法成立")

        try:
            pca = run_pca(normalized, config.n_pcs)
        except ValueError as exc:
            raise CellFlowDataError(f"PCA 无法成立：{exc}") from exc

        selection = select_cluster_count(pca.scores, config.seed)
        if selection is None:
            upper = min(10, n_kept_cells)
            raise CellFlowDataError(
                f"簇数自动选择失败：候选范围 k=2..{upper} 内无法形成两个以上"
                f"不同簇（方差不足）"
            )
        clustering = selection.clustering
    else:
        if config.n_clusters > n_kept_cells:
            raise CellFlowDataError(
                f"簇数 {config.n_clusters} 大于质控后细胞数 {n_kept_cells}，聚类无法成立"
            )

        normalized = normalize_and_select_hvg(matrix, qc, n_hvg=config.n_hvg)
        if not normalized.selected_genes:
            raise CellFlowDataError("高变基因选择结果为空，PCA 无法成立")

        try:
            pca = run_pca(normalized, config.n_pcs)
        except ValueError as exc:
            raise CellFlowDataError(f"PCA 无法成立：{exc}") from exc

        try:
            clustering = kmeans(pca.scores, config.n_clusters, config.seed)
        except ValueError as exc:
            raise CellFlowDataError(f"聚类无法成立：{exc}") from exc

        n_found_clusters = len(set(clustering.labels))
        if n_found_clusters < config.n_clusters:
            # 例如全部细胞表达相同：中心必然重合，无法形成 k 个不同簇
            raise CellFlowDataError(
                f"聚类无法成立：请求 {config.n_clusters} 个簇，"
                f"但数据仅能支撑 {n_found_clusters} 个不同簇（方差不足）"
            )
        selection = None

    markers = find_markers(normalized, clustering.labels)

    # 能执行到此处说明实际簇数 >= 2（上方已对不足两簇的数据错误拒绝），
    # 成对比较覆盖实际出现标签的全部 a < b 组合
    pairwise_markers = find_pairwise_markers(normalized, clustering.labels)

    cluster_sizes: Dict[int, int] = {}
    for label in clustering.labels:
        cluster_sizes[label] = cluster_sizes.get(label, 0) + 1

    run_info: Dict[str, Any] = {
        "version": __version__,
        "input": {
            "path": matrix.path,
            "sha256": matrix.sha256,
            "n_genes": matrix.n_genes,
            "n_cells": matrix.n_cells,
            "total_counts": matrix.total_counts,
        },
        "parameters": config.public_parameters(),
        "stage_counts": {
            "input_cells": matrix.n_cells,
            "input_genes": matrix.n_genes,
            "cells_after_qc": n_kept_cells,
            "genes_after_qc": len(qc.kept_genes),
            "cells_filtered_out": matrix.n_cells - n_kept_cells,
            "genes_filtered_out": matrix.n_genes - len(qc.kept_genes),
            "highly_variable_genes": len(normalized.selected_genes),
            "max_pcs_parameter": MAX_PCS,
            "pcs": pca.n_pcs,
            # auto 模式下为选中的 k；显式模式下与请求簇数一致
            "clusters": (
                selection.selected_k if selection is not None
                else len(cluster_sizes)
            ),
            "cluster_sizes": [cluster_sizes[k] for k in sorted(cluster_sizes)],
            "marker_tests": len(markers) * len(qc.kept_genes),
            "kmeans_iterations": clustering.iterations,
        },
        "random_seed": config.seed,
    }

    artifacts = Artifacts(
        qc=qc,
        data=normalized,
        pca=pca,
        clustering=clustering,
        markers=markers,
        pairwise_markers=pairwise_markers,
        run_info=run_info,
        cluster_selection=selection,
    )
    # 全部计算已完成，此处才创建/校验输出目录并原子写出
    ensure_output_dir(config.output_dir)
    return write_results(config.output_dir, artifacts)
