"""分析管线编排：参数校验 -> 读入 -> QC -> 归一化/HVG -> 批次校正 -> PCA
-> 聚类 -> 差异表达 -> 写出。"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Union

from . import __version__
from .batch import (
    MIXING_N_NEIGHBORS,
    BatchMetadata,
    BatchMixing,
    BatchSummaryRow,
    batch_mixing_score,
    read_batch_metadata,
)
from .errors import (
    CellFlowBatchError,
    CellFlowConfigError,
    CellFlowDataError,
    CellFlowInputError,
)
from .gene_sets import read_gene_sets, score_gene_sets
from .io import find_duplicate, read_matrix
from .kmeans import kmeans
from .markers import find_group_markers, find_markers, find_pairwise_markers
from .metadata import read_metadata
from .mtx import read_mtx_directory
from .normalize import normalize_and_select_hvg
from .output import Artifacts, publish_results
from .pca import MAX_PCS, run_pca
from .qc import compute_qc
from .selection import select_cluster_count

DEFAULT_SEED = 20240617
DEFAULT_BATCH_COLUMN = "batch"
DEFAULT_SAMPLE_COLUMN = "sample_id"
AUTO = "auto"
FORMAT_TSV = "tsv"
FORMAT_MTX = "mtx"
INPUT_FORMATS = (FORMAT_TSV, FORMAT_MTX)


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
    input_format: str = FORMAT_TSV
    metadata_path: Optional[str] = None
    batch_metadata_path: Optional[str] = None
    batch_column: str = DEFAULT_BATCH_COLUMN
    sample_column: str = DEFAULT_SAMPLE_COLUMN
    gene_sets_path: Optional[str] = None

    @property
    def auto_clusters(self) -> bool:
        return self.n_clusters == AUTO

    def public_parameters(self) -> Dict[str, Any]:
        parameters = {
            "min_genes": self.min_genes,
            "max_mito_fraction": self.max_mito_fraction,
            "min_cells": self.min_cells,
            "mito_prefix": self.mito_prefix,
            "n_hvg": self.n_hvg,
            "n_pcs": self.n_pcs,
            "n_clusters": self.n_clusters,
            "seed": self.seed,
        }
        return parameters


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
    if config.input_format not in INPUT_FORMATS:
        errors.append("--input-format 只能是 tsv 或 mtx")
    if config.batch_metadata_path is not None:
        if not isinstance(config.batch_column, str) or config.batch_column == "":
            errors.append("--batch-column 必须是非空字符串")
        if not isinstance(config.sample_column, str) or config.sample_column == "":
            errors.append("--sample-column 必须是非空字符串")
        if config.batch_column == "cell_id":
            errors.append("--batch-column 不能使用固定条码列名 cell_id")
        if config.sample_column == "cell_id":
            errors.append("--sample-column 不能使用固定条码列名 cell_id")
        if config.batch_column == config.sample_column:
            errors.append("--batch-column 与 --sample-column 不能同名")
    if errors:
        raise CellFlowConfigError("；".join(errors))


def _read_batch_metadata(config: Config, cell_ids: List[str]) -> BatchMetadata:
    """读取批次元数据；读入层错误（损坏 gzip、非 UTF-8 等）统一转为
    ValueError（CellFlowBatchError），与内容冲突使用同一异常类型。"""
    try:
        return read_batch_metadata(
            config.batch_metadata_path,
            cell_ids,
            sample_column=config.sample_column,
            batch_column=config.batch_column,
        )
    except CellFlowInputError as exc:
        raise CellFlowBatchError(str(exc)) from exc


def run(config: Config) -> List[str]:
    validate_config(config)

    # 先完成读入与全部计算，最后再触碰输出目录：
    # 任何分析阶段失败都不应留下空目录或半成品结果
    batch_enabled = config.batch_metadata_path is not None
    if config.input_format == FORMAT_MTX:
        matrix = read_mtx_directory(
            config.input_path, allow_duplicate_cells=batch_enabled
        )
    else:
        matrix = read_matrix(
            config.input_path, allow_duplicate_cells=batch_enabled
        )

    # 批次功能依赖条码一一关联：矩阵条码重复必须作为 ValueError 冲突
    # 在降维聚类之前失败（默认无批次路径仍由读取器直接报输入错误）
    if batch_enabled:
        duplicate_barcode = find_duplicate(matrix.cell_ids)
        if duplicate_barcode:
            raise CellFlowBatchError(
                f"表达矩阵包含重复细胞条码：{duplicate_barcode!r}，"
                f"无法与批次元数据一一关联"
            )

    # 元数据属于输入：先完成读取与校验，任何不合法都在触碰输出目录之前失败
    metadata = None
    if config.metadata_path is not None:
        metadata = read_metadata(config.metadata_path, matrix.cell_ids)
    batch_metadata = None
    if batch_enabled:
        batch_metadata = _read_batch_metadata(config, matrix.cell_ids)
    gene_sets = None
    if config.gene_sets_path is not None:
        gene_sets = read_gene_sets(config.gene_sets_path)

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

    if metadata is not None:
        kept_groups = {
            metadata.groups[matrix.cell_ids[c]] for c in qc.kept_cells
        }
        if len(kept_groups) < 2:
            raise CellFlowDataError(
                f"质控后非空分组仅 {len(kept_groups)} 个，不足两个，"
                f"无法进行分组差异表达"
            )

    # 与保留细胞（列序）对齐的批次/样本标签。单批次是合法输入：
    # 不做校正、不引入任何数值扰动，按原流程继续分析并产出完整结果。
    batch_labels: Optional[List[str]] = None
    sample_labels: Optional[List[str]] = None
    correction_labels: Optional[List[str]] = None
    n_kept_batches = 0
    n_kept_samples = 0
    if batch_metadata is not None:
        batch_labels = [
            batch_metadata.batches[matrix.cell_ids[c]] for c in qc.kept_cells
        ]
        sample_labels = [
            batch_metadata.samples[matrix.cell_ids[c]] for c in qc.kept_cells
        ]
        n_kept_batches = len(set(batch_labels))
        n_kept_samples = len(set(sample_labels))
        correction_applied = n_kept_batches >= 2
        correction_labels = batch_labels if correction_applied else None

    if not config.auto_clusters and config.n_clusters > n_kept_cells:
        raise CellFlowDataError(
            f"簇数 {config.n_clusters} 大于质控后细胞数 {n_kept_cells}，聚类无法成立"
        )

    normalized = normalize_and_select_hvg(
        matrix, qc, n_hvg=config.n_hvg, batch_labels=correction_labels
    )
    if not normalized.selected_genes:
        raise CellFlowDataError("高变基因选择结果为空，PCA 无法成立")

    # 校正前的对照低维表示：仅多批次时额外计算一次完整的未校正归一化与
    # PCA（HVG 同样从未校正值选择），用于校正前批次混合分数；
    # 单批次不校正，校正前后即同一表示，无需重复计算。
    pca_before = None
    if correction_labels is not None:
        normalized_before = normalize_and_select_hvg(
            matrix, qc, n_hvg=config.n_hvg, batch_labels=None
        )
        try:
            pca_before = run_pca(normalized_before, config.n_pcs)
        except ValueError as exc:
            raise CellFlowDataError(f"PCA 无法成立：{exc}") from exc

    try:
        pca = run_pca(normalized, config.n_pcs)
    except ValueError as exc:
        raise CellFlowDataError(f"PCA 无法成立：{exc}") from exc

    selection = None
    if config.auto_clusters:
        # 自动模式：候选 k=2..min(10, 质控后细胞数)，PCA 后逐个评估
        selection_result = select_cluster_count(pca.scores, config.seed)
        if selection_result is None:
            upper = min(10, n_kept_cells)
            raise CellFlowDataError(
                f"簇数自动选择失败：候选范围 k=2..{upper} 内无法形成两个以上"
                f"不同簇（方差不足）"
            )
        clustering = selection_result.clustering
        selection = selection_result
    else:
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

    # 批次混合分数：每个细胞 15 个最近邻中来自其他批次的占比，再对全部
    # 细胞求平均；纯距离排序、下标打破并列，不使用随机数，结果确定。
    batch_mixing: Optional[BatchMixing] = None
    if batch_metadata is not None:
        after_score = batch_mixing_score(pca.scores, batch_labels)
        if pca_before is not None:
            before_score = batch_mixing_score(pca_before.scores, batch_labels)
        else:
            # 单批次：未发生校正，校正前后分数相同
            before_score = after_score
        batch_mixing = BatchMixing(
            before_score=before_score,
            after_score=after_score,
            applied=correction_labels is not None,
            n_neighbors=min(MIXING_N_NEIGHBORS, n_kept_cells - 1),
        )

    markers = find_markers(normalized, clustering.labels)

    # 能执行到此处说明实际簇数 >= 2（上方已对不足两簇的数据错误拒绝），
    # 成对比较覆盖实际出现标签的全部 a < b 组合
    pairwise_markers = find_pairwise_markers(normalized, clustering.labels)

    # 分组差异表达：只用质控后保留细胞与基因的 log 归一化表达
    group_markers = None
    if metadata is not None:
        cell_groups = [metadata.groups[cell_id] for cell_id in normalized.cell_ids]
        group_markers = find_group_markers(normalized, cell_groups)

    # 基因集评分：用保留基因/细胞；有批次校正取中心化值，否则取 log 归一化值。
    # 任一集合与保留基因交集为空即在此处（触碰输出目录之前）报数据错误（退出码 4）
    gene_set_scores = None
    if gene_sets is not None:
        gene_set_scores = score_gene_sets(
            gene_sets,
            kept_gene_ids=normalized.gene_ids,
            cell_ids=normalized.cell_ids,
            analysis_values=normalized.analysis_values,
        )

    # 批次汇总：按 batch 升序，后三项为批次内保留细胞均值
    batch_summary: Optional[List[BatchSummaryRow]] = None
    if batch_metadata is not None:
        cells_by_batch: Dict[str, List[int]] = {}
        for c in qc.kept_cells:
            batch = batch_metadata.batches[matrix.cell_ids[c]]
            cells_by_batch.setdefault(batch, []).append(c)
        batch_summary = []
        for batch in sorted(cells_by_batch):
            members = [qc.cell_qc[c] for c in cells_by_batch[batch]]
            n = len(members)
            batch_summary.append(
                BatchSummaryRow(
                    batch_id=batch,
                    n_cells=n,
                    mean_total_counts=sum(m.total_counts for m in members) / n,
                    mean_detected_genes=sum(m.detected_genes for m in members) / n,
                    mean_mitochondrial_fraction=(
                        sum(m.mitochondrial_fraction for m in members) / n
                    ),
                )
            )

    cluster_sizes: Dict[int, int] = {}
    for label in clustering.labels:
        cluster_sizes[label] = cluster_sizes.get(label, 0) + 1

    input_info: Dict[str, Any] = {
        "path": matrix.path,
        "sha256": matrix.sha256,
        "n_genes": matrix.n_genes,
        "n_cells": matrix.n_cells,
        "total_counts": matrix.total_counts,
    }
    if matrix.input_format == FORMAT_MTX:
        # mtx 运行仅额外记录格式与三个输入文件的来源；tsv 运行字段保持不变
        input_info["input_format"] = FORMAT_MTX
        input_info["files"] = [
            {"name": f.name, "sha256": f.sha256} for f in matrix.files
        ]
    if metadata is not None:
        # 仅记录元数据来源与分组规模；既有 input 字段保持不变
        input_info["metadata"] = {
            "name": metadata.name,
            "sha256": metadata.sha256,
            "group_sizes": metadata.group_sizes,
        }
    if batch_metadata is not None:
        # 记录批次元数据来源、实际列名、质控前后各批次/样本细胞数与透传列；
        # 既有 input 字段不变
        kept_batch_sizes: Dict[str, int] = {}
        kept_sample_sizes: Dict[str, int] = {}
        for c in qc.kept_cells:
            batch = batch_metadata.batches[matrix.cell_ids[c]]
            sample = batch_metadata.samples[matrix.cell_ids[c]]
            kept_batch_sizes[batch] = kept_batch_sizes.get(batch, 0) + 1
            kept_sample_sizes[sample] = kept_sample_sizes.get(sample, 0) + 1
        sample_sizes: Dict[str, int] = {}
        for sample in batch_metadata.samples.values():
            sample_sizes[sample] = sample_sizes.get(sample, 0) + 1
        input_info["batch_metadata"] = {
            "name": batch_metadata.name,
            "sha256": batch_metadata.sha256,
            "cell_id_column": "cell_id",
            "sample_column": batch_metadata.sample_column,
            "batch_column": batch_metadata.batch_column,
            "batch_sizes": batch_metadata.batch_sizes,
            "batch_sizes_after_qc": {
                batch: kept_batch_sizes[batch]
                for batch in sorted(kept_batch_sizes)
            },
            "sample_sizes": {
                sample: sample_sizes[sample] for sample in sorted(sample_sizes)
            },
            "sample_sizes_after_qc": {
                sample: kept_sample_sizes[sample]
                for sample in sorted(kept_sample_sizes)
            },
            "extra_columns": list(batch_metadata.extra_columns),
        }
    if gene_sets is not None and gene_set_scores is not None:
        # 记录基因集来源、原始字节 SHA-256 与各集合总数/实际使用基因数；
        # 仅提供 --gene-sets 时出现，既有 input 字段不变
        input_info["gene_sets"] = {
            "name": gene_sets.name,
            "sha256": gene_sets.sha256,
            "set_sizes": {
                set_id: gene_set_scores.totals[set_id]
                for set_id in gene_set_scores.set_order
            },
            "set_sizes_used": {
                set_id: gene_set_scores.used[set_id]
                for set_id in gene_set_scores.set_order
            },
        }

    parameters = config.public_parameters()
    stage_counts: Dict[str, Any] = {
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
    }
    batch_mixing_info: Optional[Dict[str, Any]] = None
    if batch_metadata is not None:
        # 仅批次运行记录批次参数与计数；无批次运行时 parameters/stage_counts
        # 与基线完全一致
        parameters["batch_correction"] = {
            "method": "batch_mean_centering",
            "applied": correction_labels is not None,
            "n_batches": len(batch_metadata.batch_sizes),
            "n_batches_after_qc": n_kept_batches,
            "mixing_n_neighbors": MIXING_N_NEIGHBORS,
        }
        if correction_labels is not None:
            # 实际执行中心化时保留既有参数标记；单批次不校正则不出现该键
            parameters["batch_mean_centering"] = True
        stage_counts["batches"] = len(batch_metadata.batch_sizes)
        stage_counts["batches_after_qc"] = n_kept_batches
        stage_counts["samples"] = len(set(batch_metadata.samples.values()))
        stage_counts["samples_after_qc"] = n_kept_samples
        batch_mixing_info = {
            "n_neighbors": batch_mixing.n_neighbors,
            "correction_applied": batch_mixing.applied,
            "before_correction": batch_mixing.before_score,
            "after_correction": batch_mixing.after_score,
        }

    run_info: Dict[str, Any] = {
        "version": __version__,
        "input": input_info,
        "parameters": parameters,
        "stage_counts": stage_counts,
        "random_seed": config.seed,
    }
    if batch_mixing_info is not None:
        run_info["batch_mixing"] = batch_mixing_info

    artifacts = Artifacts(
        qc=qc,
        data=normalized,
        pca=pca,
        clustering=clustering,
        markers=markers,
        pairwise_markers=pairwise_markers,
        run_info=run_info,
        cluster_selection=selection,
        group_markers=group_markers,
        batch_summary=batch_summary,
        gene_set_scores=gene_set_scores,
        batch_metadata=batch_metadata,
        batch_labels=batch_labels,
        sample_labels=sample_labels,
        batch_mixing=batch_mixing,
    )
    # 全部计算已完成才触碰文件系统：预检与写出都在 publish_results 内，
    # 任一分析阶段失败时不会创建或改动目标目录；写出阶段任何文件系统故障
    # 都回滚暂存并抛 OutputPathError。
    return publish_results(config.output_dir, artifacts)
