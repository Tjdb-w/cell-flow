"""分析管线编排：参数校验 -> 读入 -> QC -> 归一化/HVG -> PCA -> 聚类 -> 差异表达 -> 写出。"""

from dataclasses import dataclass, replace
from typing import Any, Dict, List, Optional, Union

from . import __version__
from .batch import BatchSummaryRow, read_batch_metadata
from .cell_metadata import (
    DEFAULT_BATCH_COLUMN,
    DEFAULT_SAMPLE_COLUMN,
    CellBatchReport,
    batch_mixing_score,
    read_cell_metadata,
)
from .errors import CellFlowConfigError, CellFlowDataError
from .gene_sets import read_gene_sets, score_gene_sets
from .io import ExpressionMatrix, read_matrix
from .kmeans import kmeans
from .markers import find_group_markers, find_markers, find_pairwise_markers
from .metadata import read_metadata
from .mtx import read_mtx_directory
from .normalize import normalize_and_select_hvg
from .output import Artifacts, publish_results
from .pca import MAX_PCS, run_pca
from .qc import compute_qc
from .replicate import (
    PseudobulkComparison,
    PseudobulkExpression,
    build_pseudobulk,
    find_pseudobulk_markers,
    read_replicate_metadata,
)
from .selection import select_cluster_count

DEFAULT_SEED = 20240617
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
    gene_sets_path: Optional[str] = None
    cell_metadata_path: Optional[str] = None
    batch_column: str = DEFAULT_BATCH_COLUMN
    sample_column: str = DEFAULT_SAMPLE_COLUMN
    replicate_metadata_path: Optional[str] = None

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
        if self.batch_metadata_path is not None:
            # 仅批次校正运行记录该参数；无批次运行时 parameters 与基线一致
            parameters["batch_mean_centering"] = True
        if self.gene_sets_path is not None:
            # 仅基因集评分运行记录该参数；未提供时 parameters 与基线一致
            parameters["gene_set_scoring"] = True
        if self.cell_metadata_path is not None:
            # 仅跨样本批次校正运行记录列名参数；未提供时 parameters 与基线一致
            parameters["batch_column"] = self.batch_column
            parameters["sample_column"] = self.sample_column
        if self.replicate_metadata_path is not None:
            # 仅 pseudobulk 运行记录该参数；未提供时 parameters 与基线一致
            parameters["pseudobulk_de"] = True
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
    if config.cell_metadata_path is not None:
        if config.batch_metadata_path is not None:
            errors.append("--cell-metadata 与 --batch-metadata 不能同时使用")
        if not isinstance(config.batch_column, str) or config.batch_column == "":
            errors.append("--batch-column 必须是非空字符串")
        elif config.batch_column == "cell_id":
            errors.append("--batch-column 不能使用细胞条码列名 cell_id")
        if not isinstance(config.sample_column, str) or config.sample_column == "":
            errors.append("--sample-column 必须是非空字符串")
        elif config.sample_column == "cell_id":
            errors.append("--sample-column 不能使用细胞条码列名 cell_id")
        if (
            isinstance(config.batch_column, str)
            and isinstance(config.sample_column, str)
            and config.batch_column != ""
            and config.sample_column != ""
            and config.batch_column == config.sample_column
        ):
            errors.append("--batch-column 与 --sample-column 不能相同")
    if errors:
        raise CellFlowConfigError("；".join(errors))


def run(config: Config) -> List[str]:
    validate_config(config)

    # 先完成读入与全部计算，最后再触碰输出目录：
    # 任何分析阶段失败都不应留下空目录或半成品结果
    if config.input_format == FORMAT_MTX:
        matrix = read_mtx_directory(config.input_path)
    else:
        matrix = read_matrix(config.input_path)

    # 元数据属于输入：先完成读取与校验，任何不合法都在触碰输出目录之前失败
    metadata = None
    if config.metadata_path is not None:
        metadata = read_metadata(config.metadata_path, matrix.cell_ids)
    batch_metadata = None
    if config.batch_metadata_path is not None:
        batch_metadata = read_batch_metadata(
            config.batch_metadata_path, matrix.cell_ids
        )
    # 跨样本批次校正的细胞元数据同属输入：先完成读取与校验，
    # 任何不合法（统一为 ValueError）都在降维聚类之前、触碰输出目录之前失败
    cell_metadata = None
    if config.cell_metadata_path is not None:
        cell_metadata = read_cell_metadata(
            config.cell_metadata_path,
            matrix.cell_ids,
            batch_column=config.batch_column,
            sample_column=config.sample_column,
        )
    gene_sets = None
    if config.gene_sets_path is not None:
        gene_sets = read_gene_sets(config.gene_sets_path)
    replicate_metadata = None
    if config.replicate_metadata_path is not None:
        replicate_metadata = read_replicate_metadata(
            config.replicate_metadata_path, matrix.cell_ids
        )

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

    batch_labels: Optional[List[str]] = None
    if batch_metadata is not None:
        # 与保留细胞（列序）对齐的批次标签；质控后不足两个批次无法校正
        batch_labels = [
            batch_metadata.batches[matrix.cell_ids[c]] for c in qc.kept_cells
        ]
        n_kept_batches = len(set(batch_labels))
        if n_kept_batches < 2:
            raise CellFlowDataError(
                f"质控后批次仅 {n_kept_batches} 个，不足两个，"
                f"无法进行批次校正"
            )

    # 跨样本批次校正：与保留细胞对齐的批次标签。质控后仅一个批次时
    # 不施加任何扰动（校正即恒等），下游结果与无批次运行逐字节一致
    cell_batch_labels: Optional[List[str]] = None
    correction_applied = False
    if cell_metadata is not None:
        cell_batch_labels = [
            cell_metadata.batches[matrix.cell_ids[c]] for c in qc.kept_cells
        ]
        correction_applied = len(set(cell_batch_labels)) >= 2
    # 两种批次输入互斥（validate_config 已拒绝并用），这里取实际生效的标签
    effective_batch_labels = batch_labels
    if effective_batch_labels is None and correction_applied:
        effective_batch_labels = cell_batch_labels

    if config.auto_clusters:
        # 自动模式：候选 k=2..min(10, 质控后细胞数)，PCA 后逐个评估
        normalized = normalize_and_select_hvg(
            matrix, qc, n_hvg=config.n_hvg, batch_labels=effective_batch_labels
        )
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

        normalized = normalize_and_select_hvg(
            matrix, qc, n_hvg=config.n_hvg, batch_labels=effective_batch_labels
        )
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

    # 跨样本批次校正运行始终公开校正后表达：多批次为均值中心化值，
    # 单批次为与未校正值逐字节一致的拷贝（校正即恒等，不引入扰动）
    if cell_metadata is not None and normalized.corrected_values is None:
        normalized = replace(
            normalized,
            corrected_values=[list(row) for row in normalized.values],
        )

    # 批次混合分数：校正后在实际用于聚类的 PCA 坐标上计算；
    # 校正前在未校正的归一化 -> HVG -> PCA 坐标上计算（单批次时两者相同）
    cell_batch_report: Optional[CellBatchReport] = None
    if cell_metadata is not None:
        assert cell_batch_labels is not None
        if correction_applied:
            before_normalized = normalize_and_select_hvg(
                matrix, qc, n_hvg=config.n_hvg, batch_labels=None
            )
            if not before_normalized.selected_genes:
                # 未校正数据高变基因为空时退回全部保留基因，
                # 保证校正前低维表示始终可计算
                before_normalized = replace(
                    before_normalized,
                    selected_genes=list(range(len(before_normalized.gene_ids))),
                )
            try:
                before_pca = run_pca(before_normalized, config.n_pcs)
            except ValueError as exc:
                raise CellFlowDataError(
                    f"校正前 PCA 无法成立：{exc}"
                ) from exc
            before_scores = before_pca.scores
        else:
            before_scores = pca.scores
        kept_cell_ids = [matrix.cell_ids[c] for c in qc.kept_cells]
        cell_batch_report = CellBatchReport(
            columns=cell_metadata.columns,
            rows=[cell_metadata.rows[cell_id] for cell_id in kept_cell_ids],
            batches=cell_batch_labels,
            mixing_before=batch_mixing_score(before_scores, cell_batch_labels),
            mixing_after=batch_mixing_score(pca.scores, cell_batch_labels),
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

    # pseudobulk 差异表达：只用 QC 保留细胞与基因的原始计数，按 sample_id
    # 求和，样本文库归一到 10000 后 log1p，以样本为观测值做 one-vs-rest
    # 与两两 Welch t 检验。质控后无样本细胞、分组不足两个或组内有效重复
    # 不足两个都在此处（触碰输出目录之前）报数据错误（退出码 4）
    pseudobulk: Optional[PseudobulkExpression] = None
    pseudobulk_comparisons: Optional[List[PseudobulkComparison]] = None
    if replicate_metadata is not None:
        pseudobulk = build_pseudobulk(matrix, qc, replicate_metadata)
        pseudobulk_comparisons = find_pseudobulk_markers(pseudobulk)

    # 批次汇总：按 batch 升序，后三项为批次内保留细胞均值
    batch_summary: Optional[List[BatchSummaryRow]] = None
    if batch_metadata is not None or cell_metadata is not None:
        batch_of = (
            batch_metadata.batches
            if batch_metadata is not None
            else cell_metadata.batches
        )
        cells_by_batch: Dict[str, List[int]] = {}
        for c in qc.kept_cells:
            batch = batch_of[matrix.cell_ids[c]]
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
        # 记录批次元数据来源与质控前后各批次细胞数；既有 input 字段不变
        kept_batch_sizes: Dict[str, int] = {}
        for c in qc.kept_cells:
            batch = batch_metadata.batches[matrix.cell_ids[c]]
            kept_batch_sizes[batch] = kept_batch_sizes.get(batch, 0) + 1
        input_info["batch_metadata"] = {
            "name": batch_metadata.name,
            "sha256": batch_metadata.sha256,
            "batch_sizes": batch_metadata.batch_sizes,
            "batch_sizes_after_qc": {
                batch: kept_batch_sizes[batch]
                for batch in sorted(kept_batch_sizes)
            },
        }
    if cell_metadata is not None:
        # 记录细胞元数据来源、列名与质控前后各批次细胞数；既有 input 字段不变
        kept_cell_batch_sizes: Dict[str, int] = {}
        for c in qc.kept_cells:
            batch = cell_metadata.batches[matrix.cell_ids[c]]
            kept_cell_batch_sizes[batch] = kept_cell_batch_sizes.get(batch, 0) + 1
        input_info["cell_metadata"] = {
            "name": cell_metadata.name,
            "sha256": cell_metadata.sha256,
            "batch_column": config.batch_column,
            "sample_column": config.sample_column,
            "batch_sizes": cell_metadata.batch_sizes,
            "batch_sizes_after_qc": {
                batch: kept_cell_batch_sizes[batch]
                for batch in sorted(kept_cell_batch_sizes)
            },
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
    if replicate_metadata is not None:
        # 记录重复元数据来源、原始字节 SHA-256 与质控前后样本/分组规模；
        # 仅提供 --replicate-metadata 时出现，既有 input 字段不变
        assert pseudobulk is not None
        kept_sample_sizes: Dict[str, int] = {}
        for s, sample_id in enumerate(pseudobulk.sample_ids):
            kept_sample_sizes[sample_id] = pseudobulk.sample_cell_counts[s]
        kept_group_sizes: Dict[str, int] = {}
        for group, count in zip(
            pseudobulk.sample_groups, pseudobulk.sample_cell_counts
        ):
            kept_group_sizes[group] = kept_group_sizes.get(group, 0) + count
        input_info["replicate_metadata"] = {
            "name": replicate_metadata.name,
            "sha256": replicate_metadata.sha256,
            "samples_before_qc": len(replicate_metadata.sample_group),
            "samples_after_qc": len(pseudobulk.sample_ids),
            "groups_before_qc": len(replicate_metadata.group_sizes),
            "groups_after_qc": len(kept_group_sizes),
            "sample_sizes": replicate_metadata.sample_sizes,
            "sample_sizes_after_qc": {
                sample: kept_sample_sizes[sample]
                for sample in sorted(kept_sample_sizes)
            },
            "group_sizes": replicate_metadata.group_sizes,
            "group_sizes_after_qc": {
                group: kept_group_sizes[group]
                for group in sorted(kept_group_sizes)
            },
        }

    parameters = config.public_parameters()
    if cell_metadata is not None and correction_applied:
        # 实际施加批次均值中心化时记录；单批次未施加扰动则不记录
        parameters["batch_mean_centering"] = True

    run_info: Dict[str, Any] = {
        "version": __version__,
        "input": input_info,
        "parameters": parameters,
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
    if cell_batch_report is not None:
        # 跨样本批次校正汇总：方法、批次数与校正前后混合分数；既有字段不变
        assert cell_batch_labels is not None
        run_info["batch_correction"] = {
            "method": (
                "batch_mean_centering"
                if correction_applied
                else "none_single_batch"
            ),
            "n_batches": len(set(cell_batch_labels)),
            "mixing_score_before": cell_batch_report.mixing_before,
            "mixing_score_after": cell_batch_report.mixing_after,
        }
    if pseudobulk is not None and pseudobulk_comparisons is not None:
        # pseudobulk 差异表达汇总：归一目标、质控后样本/分组规模、
        # 比较数与逐基因检验数；仅提供 --replicate-metadata 时出现
        pb_groups = sorted(set(pseudobulk.sample_groups))
        group_replicate_sizes = {
            group: pseudobulk.sample_groups.count(group) for group in pb_groups
        }
        n_pb_groups = len(pb_groups)
        n_comparisons = len(pseudobulk_comparisons)
        run_info["pseudobulk_de"] = {
            "library_size_target": 10000,
            "n_samples": len(pseudobulk.sample_ids),
            "n_groups": n_pb_groups,
            "group_replicate_sizes": group_replicate_sizes,
            "n_comparisons": n_comparisons,
            "n_one_vs_rest": n_pb_groups,
            "n_pairwise": n_comparisons - n_pb_groups,
            "marker_tests": n_comparisons * len(pseudobulk.gene_ids),
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
        group_markers=group_markers,
        batch_summary=batch_summary,
        gene_set_scores=gene_set_scores,
        cell_batch_report=cell_batch_report,
        pseudobulk=pseudobulk,
        pseudobulk_comparisons=pseudobulk_comparisons,
    )
    # 全部计算已完成才触碰文件系统：预检与写出都在 publish_results 内，
    # 任一分析阶段失败时不会创建或改动目标目录；写出阶段任何文件系统故障
    # 都回滚暂存并抛 OutputPathError。
    return publish_results(config.output_dir, artifacts)
