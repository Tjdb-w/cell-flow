"""分析管线编排：参数校验 -> 读入 -> QC -> 归一化/HVG -> PCA -> 聚类 -> 差异表达 -> 写出。"""

import math
from dataclasses import dataclass, replace
from typing import Any, Dict, List, Optional, Union

from . import __version__
from .abundance import DifferentialAbundance, compute_differential_abundance
from .adjusted_cluster_pseudobulk import (
    AdjustedClusterPseudobulkDE,
    compute_adjusted_cluster_pseudobulk_de,
)
from .batch import BatchSummaryRow, read_batch_metadata
from .cell_metadata import (
    DEFAULT_BATCH_COLUMN,
    DEFAULT_SAMPLE_COLUMN,
    CellBatchReport,
    batch_mixing_score,
    read_cell_metadata,
)
from .cell_types import (
    CellTypeAnnotations,
    annotate_clusters,
    read_cell_type_reference,
)
from .cluster_pseudobulk import (
    ClusterPseudobulkDE,
    compute_cluster_pseudobulk_de,
)
from .doublets import DoubletResult, detect_doublets
from .enrichment import MarkerEnrichment, enrich_marker_gene_sets
from .errors import CellFlowConfigError, CellFlowDataError, CellFlowError
from .gene_sets import read_gene_sets, score_gene_sets
from .io import ExpressionMatrix, read_matrix
from .kmeans import kmeans
from .markers import find_group_markers, find_markers, find_pairwise_markers
from .metadata import read_metadata
from .mtx import read_mtx_directory
from .normalize import normalize_and_select_hvg
from .numeric_cluster_pseudobulk import (
    NumericClusterPseudobulkDE,
    compute_numeric_cluster_pseudobulk_de,
)
from .numeric_pseudobulk_covariates import (
    PseudobulkNumericCovariates,
    read_pseudobulk_numeric_covariates,
)
from .output import Artifacts, publish_results
from .paired_cluster_pseudobulk import (
    PairedClusterPseudobulkDE,
    compute_paired_cluster_pseudobulk_de,
)
from .paired_replicate import (
    PairedReplicateMetadata,
    read_paired_replicate_metadata,
)
from .pca import MAX_PCS, PCALoadings, compute_pca_loadings, run_pca
from .pseudobulk_covariates import PseudobulkCovariates, read_pseudobulk_covariates
from .pseudobulk_gene_sets import (
    PseudobulkGeneSetScores,
    score_pseudobulk_gene_sets,
)
from .qc import compute_qc
from .replicate import PseudobulkData, build_pseudobulk, read_replicate_metadata
from .selection import select_cluster_count
from .stability import (
    DEFAULT_STABILITY_FRACTION,
    DEFAULT_STABILITY_SAMPLES,
    DEFAULT_STABILITY_SEED,
    StabilityConfig,
    StabilityResult,
    run_stability_analysis,
    summarize_scores,
    validate_stability_settings,
)

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
    cell_type_reference_path: Optional[str] = None
    detect_doublets: bool = False
    expected_doublet_rate: float = 0.08
    stability_analysis: bool = False
    stability_n_samples: Any = DEFAULT_STABILITY_SAMPLES
    stability_sample_fraction: Any = DEFAULT_STABILITY_FRACTION
    stability_seed: Any = DEFAULT_STABILITY_SEED
    enrich_markers: bool = False
    enrichment_alpha: float = 0.05
    enrichment_min_log_fc: float = 0.0
    pseudobulk_gene_set_de: bool = False
    differential_abundance: bool = False
    cluster_pseudobulk_de: bool = False
    paired_replicate_metadata_path: Optional[str] = None
    paired_cluster_pseudobulk_de: bool = False
    pseudobulk_covariates_path: Optional[str] = None
    pseudobulk_numeric_covariates_path: Optional[str] = None
    cluster_pseudobulk_adjusted_de: bool = False
    pca_loadings: bool = False

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
        if self.pseudobulk_gene_set_de:
            # 仅显式启用 pseudobulk 基因集分组差异时记录；未启用时与基线一致
            parameters["pseudobulk_gene_set_de"] = True
        if self.detect_doublets:
            # 仅双细胞识别运行记录这两个参数；未启用时 parameters 与基线一致
            parameters["detect_doublets"] = True
            parameters["expected-doublet-rate"] = self.expected_doublet_rate
        if self.stability_analysis:
            # 仅显式启用稳定性分析时记录这三个参数；未启用时与基线逐字节一致
            parameters["stability_analysis"] = True
            parameters["stability_n_samples"] = self.stability_n_samples
            parameters["stability_sample_fraction"] = self.stability_sample_fraction
            parameters["stability_seed"] = self.stability_seed
        if self.cell_type_reference_path is not None:
            # 仅自动细胞类型注释运行记录该参数；未提供时 parameters 与基线一致
            parameters["cell_type_annotation"] = True
        if self.enrich_markers:
            # 仅显式启用 marker 基因集富集时记录这三个参数；未启用时与基线一致
            parameters["enrich_markers"] = True
            parameters["enrichment_alpha"] = self.enrichment_alpha
            parameters["enrichment_min_log_fc"] = self.enrichment_min_log_fc
        if self.differential_abundance:
            # 仅显式启用簇级样本差异丰度时记录；未启用时与基线一致
            parameters["differential_abundance"] = True
        if self.cluster_pseudobulk_de:
            # 仅显式启用簇内 pseudobulk 差异表达时记录；未启用时与基线一致
            parameters["cluster_pseudobulk_de"] = True
        if self.paired_cluster_pseudobulk_de:
            # 仅显式启用簇内配对 pseudobulk 差异表达时记录；未启用时与基线一致
            parameters["paired_cluster_pseudobulk_de"] = True
        if self.cluster_pseudobulk_adjusted_de:
            # 仅显式启用协变量校正的簇内 pseudobulk 差异表达时记录；
            # 未启用时与基线一致
            parameters["cluster_pseudobulk_adjusted_de"] = True
        if self.pca_loadings:
            # 仅显式启用 PCA 载荷输出时记录；未启用时 parameters 与基线一致
            parameters["pca_loadings"] = True
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
    if not isinstance(config.detect_doublets, bool):
        errors.append("--detect-doublets 必须是无值开关参数")
    if (
        not isinstance(config.expected_doublet_rate, (int, float))
        or isinstance(config.expected_doublet_rate, bool)
        or not 0.0 <= float(config.expected_doublet_rate) < 1.0
    ):
        errors.append("--expected-doublet-rate 必须是 [0, 1) 内的有限数值")
    if not isinstance(config.enrich_markers, bool):
        errors.append("--enrich-markers 必须是无值开关参数")
    if config.enrich_markers and config.gene_sets_path is None:
        # 富集建立在 --gene-sets 评分基线上
        errors.append("--enrich-markers 需与 --gene-sets 同时使用")
    if (
        not isinstance(config.enrichment_alpha, (int, float))
        or isinstance(config.enrichment_alpha, bool)
        or not 0.0 < float(config.enrichment_alpha) < 1.0
    ):
        errors.append("--enrichment-alpha 必须是 (0, 1) 开区间内的有限数值")
    if (
        not isinstance(config.enrichment_min_log_fc, (int, float))
        or isinstance(config.enrichment_min_log_fc, bool)
        or not math.isfinite(float(config.enrichment_min_log_fc))
    ):
        errors.append("--enrichment-min-log-fc 必须是有限数值")
    if not isinstance(config.pseudobulk_gene_set_de, bool):
        errors.append("--pseudobulk-gene-set-de 必须是无值开关参数")
    if config.pseudobulk_gene_set_de and (
        config.gene_sets_path is None or config.replicate_metadata_path is None
    ):
        # pseudobulk 基因集分组差异建立在 --gene-sets 评分与
        # --replicate-metadata pseudobulk 两条基线上
        errors.append(
            "--pseudobulk-gene-set-de 需与 --gene-sets、"
            "--replicate-metadata 同时使用"
        )
    if not isinstance(config.differential_abundance, bool):
        errors.append("--differential-abundance 必须是无值开关参数")
    if config.differential_abundance and config.replicate_metadata_path is None:
        # 簇级样本差异丰度建立在 --replicate-metadata 样本分组基线上
        errors.append("--differential-abundance 需与 --replicate-metadata 同时使用")
    if not isinstance(config.cluster_pseudobulk_de, bool):
        errors.append("--cluster-pseudobulk-de 必须是无值开关参数")
    if not isinstance(config.paired_cluster_pseudobulk_de, bool):
        errors.append("--paired-cluster-pseudobulk-de 必须是无值开关参数")
    if not isinstance(config.pca_loadings, bool):
        errors.append("--pca-loadings 必须是无值开关参数")
    if config.cluster_pseudobulk_de and config.replicate_metadata_path is None:
        # 簇内 pseudobulk 差异表达建立在 --replicate-metadata 样本分组基线上
        errors.append(
            "--cluster-pseudobulk-de 需与 --replicate-metadata 同时使用"
        )
    if config.paired_cluster_pseudobulk_de and (
        config.replicate_metadata_path is None
        or not config.cluster_pseudobulk_de
        or config.paired_replicate_metadata_path is None
    ):
        # 簇内配对 pseudobulk 差异表达建立在 --replicate-metadata 样本基线与
        # --cluster-pseudobulk-de 簇内 pseudobulk 基线之上，并须提供配对元数据
        errors.append(
            "--paired-cluster-pseudobulk-de 需与 --replicate-metadata、"
            "--cluster-pseudobulk-de、--paired-replicate-metadata 同时使用"
        )
    if not isinstance(config.cluster_pseudobulk_adjusted_de, bool):
        errors.append("--cluster-pseudobulk-adjusted-de 必须是无值开关参数")
    if config.pseudobulk_covariates_path is not None and (
        config.replicate_metadata_path is None
        or not config.cluster_pseudobulk_de
    ):
        # 样本协变量表建立在 --replicate-metadata 样本基线与
        # --cluster-pseudobulk-de 簇内 pseudobulk 基线之上
        errors.append(
            "--pseudobulk-covariates 需与 --replicate-metadata、"
            "--cluster-pseudobulk-de 同时使用"
        )
    if config.pseudobulk_numeric_covariates_path is not None and (
        config.replicate_metadata_path is None
        or not config.cluster_pseudobulk_de
        or not config.cluster_pseudobulk_adjusted_de
    ):
        # 数值型样本协变量表建立在 --replicate-metadata 样本基线与
        # --cluster-pseudobulk-de 簇内 pseudobulk 基线之上，并须启用
        # --cluster-pseudobulk-adjusted-de（可再并用分类协变量表）
        errors.append(
            "--pseudobulk-numeric-covariates 需与 --replicate-metadata、"
            "--cluster-pseudobulk-de、--cluster-pseudobulk-adjusted-de "
            "同时使用"
        )
    if config.cluster_pseudobulk_adjusted_de and (
        config.replicate_metadata_path is None
        or not config.cluster_pseudobulk_de
        or (
            config.pseudobulk_covariates_path is None
            and config.pseudobulk_numeric_covariates_path is None
        )
    ):
        # 协变量校正的簇内 pseudobulk 差异表达建立在
        # --replicate-metadata、--cluster-pseudobulk-de 基线之上，
        # 并须提供分类或数值样本协变量表（至少其一）
        errors.append(
            "--cluster-pseudobulk-adjusted-de 需与 --replicate-metadata、"
            "--cluster-pseudobulk-de 同时使用，并提供 "
            "--pseudobulk-covariates 或 --pseudobulk-numeric-covariates"
        )
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

    # 稳定性分析参数非法时统一以 ValueError 失败并指出对应配置字段；
    # 在读入矩阵与触碰输出目录之前终止，不生成任何稳定性文件
    stability_config: Optional[StabilityConfig] = None
    if config.stability_analysis:
        stability_config = validate_stability_settings(
            n_samples=config.stability_n_samples,
            sample_fraction=config.stability_sample_fraction,
            seed=config.stability_seed,
        )

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
    paired_replicate_metadata = None
    if config.paired_replicate_metadata_path is not None:
        # 配对元数据同属输入：在触碰输出目录之前完成读取与校验；配对表必须
        # 与 --replicate-metadata 的全部样本一一对应，格式/覆盖错误报输入错误
        assert replicate_metadata is not None
        paired_replicate_metadata = read_paired_replicate_metadata(
            config.paired_replicate_metadata_path,
            replicate_metadata.sample_order,
        )
    pseudobulk_covariates = None
    if config.pseudobulk_covariates_path is not None:
        # 样本协变量同属输入：在触碰输出目录之前完成读取与校验；协变量表必须
        # 与 --replicate-metadata 的全部样本一一对应，格式/覆盖错误报输入错误
        assert replicate_metadata is not None
        pseudobulk_covariates = read_pseudobulk_covariates(
            config.pseudobulk_covariates_path,
            replicate_metadata.sample_order,
        )
    pseudobulk_numeric_covariates: Optional[PseudobulkNumericCovariates] = None
    if config.pseudobulk_numeric_covariates_path is not None:
        # 数值型样本协变量同属输入：在触碰输出目录之前完成读取与校验；
        # 表必须与 --replicate-metadata 的全部样本一一对应，空值或非法数值
        # 报输入错误
        assert replicate_metadata is not None
        pseudobulk_numeric_covariates = read_pseudobulk_numeric_covariates(
            config.pseudobulk_numeric_covariates_path,
            replicate_metadata.sample_order,
        )
    # 细胞类型标记参考同属输入：先完成读取与校验，任何不合法都在触碰输出
    # 目录之前失败；标记与保留基因无交集在质控后按数据错误（退出码 4）处理
    cell_type_reference = None
    if config.cell_type_reference_path is not None:
        cell_type_reference = read_cell_type_reference(
            config.cell_type_reference_path
        )

    qc = compute_qc(
        matrix,
        min_genes=config.min_genes,
        max_mito_fraction=config.max_mito_fraction,
        min_cells=config.min_cells,
        mito_prefix=config.mito_prefix,
    )

    n_qc_cells = len(qc.kept_cells)
    if n_qc_cells < 2:
        raise CellFlowDataError(
            f"质控后仅保留 {n_qc_cells} 个细胞，不足两个，无法分析"
        )
    if not qc.kept_genes:
        raise CellFlowDataError(
            f"没有基因在至少 {config.min_cells} 个细胞中检出，无可用基因"
        )

    # 可选双细胞识别：以 QC 候选细胞 × 候选基因的原始计数子矩阵评分，
    # 按 rate 标记后在最终细胞上按 min_cells 重算基因；下游只用最终细胞与基因。
    # 未启用时 analysis_qc 即原 QC 结果，下游行为与基线逐字节一致
    doublets: Optional[DoubletResult] = None
    analysis_qc = qc
    if config.detect_doublets:
        doublets = detect_doublets(
            matrix,
            qc,
            rate=config.expected_doublet_rate,
            min_cells=config.min_cells,
            seed=config.seed,
        )
        n_after_doublet = len(doublets.kept_cells)
        if n_after_doublet < 2:
            raise CellFlowDataError(
                f"双细胞过滤后仅保留 {n_after_doublet} 个细胞，"
                f"不足两个，无法分析"
            )
        if not doublets.kept_genes:
            raise CellFlowDataError(
                f"双细胞过滤后没有基因在至少 {config.min_cells} 个细胞中检出，"
                f"无可用基因"
            )
        analysis_qc = replace(
            qc, kept_cells=doublets.kept_cells, kept_genes=doublets.kept_genes
        )

    n_kept_cells = len(analysis_qc.kept_cells)

    if metadata is not None:
        kept_groups = {
            metadata.groups[matrix.cell_ids[c]] for c in analysis_qc.kept_cells
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
            batch_metadata.batches[matrix.cell_ids[c]] for c in analysis_qc.kept_cells
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
            cell_metadata.batches[matrix.cell_ids[c]] for c in analysis_qc.kept_cells
        ]
        correction_applied = len(set(cell_batch_labels)) >= 2
    # 两种批次输入互斥（validate_config 已拒绝并用），这里取实际生效的标签
    effective_batch_labels = batch_labels
    if effective_batch_labels is None and correction_applied:
        effective_batch_labels = cell_batch_labels

    if config.auto_clusters:
        # 自动模式：候选 k=2..min(10, 质控后细胞数)，PCA 后逐个评估
        normalized = normalize_and_select_hvg(
            matrix, analysis_qc, n_hvg=config.n_hvg, batch_labels=effective_batch_labels
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
            matrix, analysis_qc, n_hvg=config.n_hvg, batch_labels=effective_batch_labels
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

    # 可选 PCA 载荷：沿用本次 PCA 的最终细胞、高变基因顺序与分析表达值，
    # 按最终（符号约定后）细胞得分与解释方差计算带符号载荷；纯增量计算，
    # 不回写任何既有中间结果，未启用时 pca_loadings 为 None
    pca_loadings: Optional[PCALoadings] = None
    if config.pca_loadings:
        pca_loadings = compute_pca_loadings(normalized, pca)

    # 可选聚类稳定性分析：以完整数据的本次聚类标签为参照，对最终保留细胞
    # （analysis_qc.kept_cells，与 clustering.labels 列序一致）不放回抽样，
    # 每次沿用同一 QC/归一化/HVG/PCA/聚类配置重算标签并在交集上求 ARI。
    # 抽样随机流只来自独立的 stability_seed，不触碰主分析的任何中间结果；
    # 未启用时 stability_result 为 None，主分析路径与结果与基线逐字节一致。
    stability_result: Optional[StabilityResult] = None
    if stability_config is not None:
        try:
            stability_result = run_stability_analysis(
                matrix,
                analysis_qc,
                config=stability_config,
                min_genes=config.min_genes,
                max_mito_fraction=config.max_mito_fraction,
                min_cells=config.min_cells,
                mito_prefix=config.mito_prefix,
                n_hvg=config.n_hvg,
                n_pcs=config.n_pcs,
                n_clusters=config.n_clusters,
                seed=config.seed,
                batch_labels=effective_batch_labels,
                reference_labels=clustering.labels,
            )
        except CellFlowError:
            # 配置/数据错误已带确定退出码，直接向上传播
            raise
        except ValueError as exc:
            # 抽样子集 PCA/聚类无法成立等意外情形，统一按数据错误（退出码 4）报告
            raise CellFlowDataError(f"聚类稳定性分析无法成立：{exc}") from exc

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
                matrix, analysis_qc, n_hvg=config.n_hvg, batch_labels=None
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
        kept_cell_ids = [matrix.cell_ids[c] for c in analysis_qc.kept_cells]
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

    # Marker 基因集富集：沿用最终细胞、保留基因、最终簇与 markers.tsv 的
    # one-versus-rest Welch t 检验结果；命中为 p_value_adj<=alpha 且
    # log_fc>=阈值 的基因，背景为保留基因，对每个 set_id 求超几何单侧 P
    # 并在每簇内做 BH 校正。仅在 --enrich-markers（且必带 --gene-sets）
    # 时执行；基因集与保留基因交集为空已在上方评分阶段报数据错误（4）。
    marker_enrichment: Optional[MarkerEnrichment] = None
    if config.enrich_markers:
        marker_enrichment = enrich_marker_gene_sets(
            gene_sets,
            markers,
            kept_gene_ids=normalized.gene_ids,
            alpha=config.enrichment_alpha,
            min_log_fc=config.enrichment_min_log_fc,
        )

    # 自动细胞类型注释：只用保留基因/细胞；有批次校正取实际聚类所用中心化值，
    # 否则取 log 归一化值。全部标记与保留基因无交集即在此处（触碰输出目录
    # 之前）报数据错误（退出码 4）。未提供参考时 cell_type_annotations 为 None
    cell_type_annotations: Optional[CellTypeAnnotations] = None
    if cell_type_reference is not None:
        cell_type_annotations = annotate_clusters(
            cell_type_reference,
            kept_gene_ids=normalized.gene_ids,
            labels=clustering.labels,
            analysis_values=normalized.analysis_values,
        )

    # 生物学重复 pseudobulk：只用质控后保留细胞与基因，按样本求和原始计数、
    # 样本文库归一到 10000 后 log1p，以样本为观测做分组差异表达。
    # 质控后样本无细胞、group 有效重复不足两个或不足两个 group 在此报数据错误
    pseudobulk: Optional[PseudobulkData] = None
    if replicate_metadata is not None:
        pseudobulk = build_pseudobulk(matrix, analysis_qc, replicate_metadata)

    # 簇级样本差异丰度：只用质控后保留细胞、最终簇标签与样本分组，
    # 以样本为观测对簇比例做分组 Welch t/BH。仅在 --differential-abundance
    # （必带 --replicate-metadata）时执行；质控后样本空、group 不足两个或
    # 任一 group 重复不足两个已由上方 build_pseudobulk 报数据错误（退出码 4）
    differential_abundance: Optional[DifferentialAbundance] = None
    if config.differential_abundance:
        assert replicate_metadata is not None and pseudobulk is not None
        cell_samples = [
            replicate_metadata.samples[matrix.cell_ids[c]]
            for c in analysis_qc.kept_cells
        ]
        differential_abundance = compute_differential_abundance(
            clustering.labels, cell_samples, pseudobulk
        )

    # 最终簇内 pseudobulk 差异表达：对每个最终簇与样本只汇总该簇内保留细胞
    # 在保留基因上的原始计数，簇内样本文库归一到 10000 后 log1p，以样本为
    # 观测做簇内分组 Welch t/比较内 BH。仅在 --cluster-pseudobulk-de
    # （必带 --replicate-metadata）时执行。某簇须有 >=2 个 group 且每组
    # >=2 个含保留细胞的有效样本才进入分析；无任何可检验簇报数据错误（4），
    # 不触碰输出目录。无细胞样本不纳入该簇也不补零。
    cluster_pseudobulk: Optional[ClusterPseudobulkDE] = None
    if config.cluster_pseudobulk_de:
        assert replicate_metadata is not None and pseudobulk is not None
        cluster_pseudobulk = compute_cluster_pseudobulk_de(
            matrix, analysis_qc, clustering.labels, replicate_metadata
        )
        if not cluster_pseudobulk.tested_clusters:
            raise CellFlowDataError(
                "没有任何最终簇具备至少两个 group 且每个 group 至少两个"
                "含保留细胞的有效样本，无法进行簇内 pseudobulk 差异表达"
            )

    # 最终簇内配对 pseudobulk 差异表达：沿用簇内 pseudobulk 的保留范围、
    # 计数汇总与文库归一化，以 a 减 b 的 pair 差值为观测做双侧配对 t 检验、
    # 簇内跨基因 BH。仅在 --paired-cluster-pseudobulk-de（必带
    # --replicate-metadata、--cluster-pseudobulk-de 与
    # --paired-replicate-metadata）时执行。质控后 group 不为两个或任一 pair
    # 不同时含两个 group 各一个样本报输入错误（退出码 2）；无任何簇具备至少
    # 两个完整 pair 报数据错误（退出码 4）；均在触碰输出目录之前失败。
    paired_cluster_pseudobulk: Optional[PairedClusterPseudobulkDE] = None
    if config.paired_cluster_pseudobulk_de:
        assert (
            replicate_metadata is not None
            and paired_replicate_metadata is not None
            and pseudobulk is not None
        )
        paired_cluster_pseudobulk = compute_paired_cluster_pseudobulk_de(
            matrix,
            analysis_qc,
            clustering.labels,
            replicate_metadata,
            paired_replicate_metadata,
            pseudobulk.sample_ids,
        )

    # 最终簇内按样本协变量校正的 pseudobulk 差异表达：沿用簇内 pseudobulk 的
    # 保留范围、计数汇总、文库归一化与簇/比较枚举次序，以样本为观测拟合
    # 截距 + 组别项 + 协变量哑变量（字典序最小水平为参照）的线性模型，
    # 组别系数即 log_fc_adjusted，簇比较内跨保留基因做 BH。该分类校正只在
    # 提供 --pseudobulk-covariates（并启用 --cluster-pseudobulk-adjusted-de，
    # 必带 --replicate-metadata、--cluster-pseudobulk-de）时执行；只给数值
    # 协变量表时分类表为 None，跳过此步而由下方数值检验独立产出。是对
    # cluster_pseudobulk_de.tsv 的补充，不改写其 Welch 结果。设计矩阵秩不足
    # 或残差自由度不大于零在此报数据错误（退出码 4），不触碰输出目录。
    adjusted_cluster_pseudobulk: Optional[AdjustedClusterPseudobulkDE] = None
    if config.cluster_pseudobulk_adjusted_de and pseudobulk_covariates is not None:
        assert (
            replicate_metadata is not None
            and cluster_pseudobulk is not None
        )
        adjusted_cluster_pseudobulk = compute_adjusted_cluster_pseudobulk_de(
            matrix,
            analysis_qc,
            clustering.labels,
            replicate_metadata,
            pseudobulk_covariates,
        )

    # 最终簇内含数值型样本协变量的 pseudobulk 差异检验：沿用簇内 pseudobulk
    # 的保留范围、计数汇总、文库归一化与簇/比较枚举次序，以样本为观测拟合
    # 截距 + 组别项 +（可并给的）分类协变量哑变量 + 数值协变量连续列的线性
    # 模型，对组别项与每个连续列给 effect、标准误、自由度、t、双侧 P，并在
    # 每个簇比较内把组别与全部连续列跨保留基因作为一个 BH 家族。仅在提供
    # --pseudobulk-numeric-covariates（必带 --replicate-metadata、
    # --cluster-pseudobulk-de 与 --cluster-pseudobulk-adjusted-de）时执行；
    # 独立于既有 Welch 与分类校正结果，不改写它们。设计矩阵秩不足或残差自由
    # 度不大于零在此报数据错误（退出码 4），不触碰输出目录。
    numeric_cluster_pseudobulk: Optional[NumericClusterPseudobulkDE] = None
    if pseudobulk_numeric_covariates is not None:
        assert (
            replicate_metadata is not None
            and cluster_pseudobulk is not None
        )
        numeric_cluster_pseudobulk = compute_numeric_cluster_pseudobulk_de(
            matrix,
            analysis_qc,
            clustering.labels,
            replicate_metadata,
            pseudobulk_numeric_covariates,
            pseudobulk_covariates,
        )

    # pseudobulk 基因集分组差异：在 pseudobulk log1p 值上按集合（与保留基因
    # 交集）给样本评分，再以样本为观测做分组 Welch t/BH。仅在
    # --pseudobulk-gene-set-de（必带 --gene-sets 与 --replicate-metadata）
    # 时执行；任一集合交集为空在此报数据错误（退出码 4），不触碰输出目录。
    pseudobulk_gene_set_scores: Optional[PseudobulkGeneSetScores] = None
    if config.pseudobulk_gene_set_de:
        assert gene_sets is not None and pseudobulk is not None
        pseudobulk_gene_set_scores = score_pseudobulk_gene_sets(gene_sets, pseudobulk)

    # 批次汇总：按 batch 升序，后三项为批次内保留细胞均值
    batch_summary: Optional[List[BatchSummaryRow]] = None
    if batch_metadata is not None or cell_metadata is not None:
        batch_of = (
            batch_metadata.batches
            if batch_metadata is not None
            else cell_metadata.batches
        )
        cells_by_batch: Dict[str, List[int]] = {}
        for c in analysis_qc.kept_cells:
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
        for c in analysis_qc.kept_cells:
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
        for c in analysis_qc.kept_cells:
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
    if replicate_metadata is not None and pseudobulk is not None:
        # 记录重复元数据来源、样本与分组规模及质控前后样本/重复数；
        # 仅提供 --replicate-metadata 时出现，既有 input 字段不变
        kept_sample_sizes: Dict[str, int] = {}
        for c in analysis_qc.kept_cells:
            sample_id = replicate_metadata.samples[matrix.cell_ids[c]]
            kept_sample_sizes[sample_id] = kept_sample_sizes.get(sample_id, 0) + 1
        kept_groups = sorted(
            {replicate_metadata.sample_group[s] for s in pseudobulk.sample_ids}
        )
        kept_group_sizes: Dict[str, int] = {}
        for sample_id in pseudobulk.sample_ids:
            group = replicate_metadata.sample_group[sample_id]
            kept_group_sizes[group] = kept_group_sizes.get(group, 0) + 1
        input_info["replicate_metadata"] = {
            "name": replicate_metadata.name,
            "sha256": replicate_metadata.sha256,
            "sample_sizes": replicate_metadata.sample_sizes,
            "sample_sizes_after_qc": {
                sample_id: kept_sample_sizes[sample_id]
                for sample_id in sorted(kept_sample_sizes)
            },
            "group_sizes": replicate_metadata.group_sizes,
            "replicates_before_qc": len(replicate_metadata.sample_order),
            "replicates_after_qc": len(pseudobulk.sample_ids),
            "group_replicates_after_qc": {
                group: kept_group_sizes[group] for group in kept_groups
            },
        }
    if (
        config.paired_cluster_pseudobulk_de
        and paired_replicate_metadata is not None
        and pseudobulk is not None
    ):
        # 配对元数据来源、原始字节 SHA-256 与配对规模；仅显式启用配对差异
        # 表达时出现，未启用时既有 input 字段逐字节不变
        present_pairs = {
            paired_replicate_metadata.sample_pair[sample_id]
            for sample_id in pseudobulk.sample_ids
        }
        input_info["paired_replicate_metadata"] = {
            "name": paired_replicate_metadata.name,
            "sha256": paired_replicate_metadata.sha256,
            "pair_count": len(paired_replicate_metadata.pair_order),
            "pairs_after_qc": len(present_pairs),
        }
    if (
        config.cluster_pseudobulk_adjusted_de
        and pseudobulk_covariates is not None
    ):
        # 样本协变量来源、原始字节 SHA-256 与协变量列；仅显式启用校正差异
        # 表达时出现，未启用时既有 input 字段逐字节不变
        input_info["pseudobulk_covariates"] = {
            "name": pseudobulk_covariates.name,
            "sha256": pseudobulk_covariates.sha256,
            "covariates": list(pseudobulk_covariates.covariate_names),
            "n_samples": len(pseudobulk_covariates.sample_levels),
        }
    if (
        pseudobulk_numeric_covariates is not None
        and numeric_cluster_pseudobulk is not None
    ):
        # 数值型样本协变量来源、原始字节 SHA-256 与协变量列；仅提供
        # --pseudobulk-numeric-covariates 时出现，未提供时既有 input 字段不变
        input_info["pseudobulk_numeric_covariates"] = {
            "name": pseudobulk_numeric_covariates.name,
            "sha256": pseudobulk_numeric_covariates.sha256,
            "covariates": list(
                pseudobulk_numeric_covariates.covariate_names
            ),
            "n_samples": len(pseudobulk_numeric_covariates.sample_values),
        }

    if cell_type_reference is not None and cell_type_annotations is not None:
        # 记录标记参考来源、原始字节 SHA-256 与各类型总标记数；
        # 仅提供 --cell-type-reference 时出现，既有 input 字段不变
        input_info["cell_type_reference"] = {
            "name": cell_type_reference.name,
            "sha256": cell_type_reference.sha256,
            "marker_counts": {
                cell_type: cell_type_annotations.totals[cell_type]
                for cell_type in cell_type_annotations.type_order
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
            "cells_after_qc": len(qc.kept_cells),
            "genes_after_qc": len(qc.kept_genes),
            "cells_filtered_out": matrix.n_cells - len(qc.kept_cells),
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
            "marker_tests": len(markers) * len(analysis_qc.kept_genes),
            "kmeans_iterations": clustering.iterations,
        },
        "random_seed": config.seed,
    }
    if doublets is not None:
        # 双细胞识别运行记录候选/标记/过滤后计数与最终细胞、基因数；
        # 未启用时 stage_counts 与基线一致
        stage_counts = run_info["stage_counts"]
        stage_counts["doublet_candidates"] = len(doublets.cell_ids)
        stage_counts["doublets_flagged"] = sum(
            1 for flag in doublets.flags if flag
        )
        stage_counts["cells_after_doublet_filter"] = len(doublets.kept_cells)
        stage_counts["genes_after_doublet_filter"] = len(doublets.kept_genes)
        stage_counts["final_cells"] = len(doublets.kept_cells)
        stage_counts["final_genes"] = len(doublets.kept_genes)
    if pseudobulk is not None:
        # pseudobulk 运行记录样本/分组数与检验计数；未提供重复元数据时字段不变
        n_pseudobulk_genes = len(pseudobulk.gene_ids)
        n_pseudobulk_comparisons = len(pseudobulk.comparisons)
        run_info["stage_counts"]["pseudobulk_samples"] = len(pseudobulk.sample_ids)
        run_info["stage_counts"]["pseudobulk_groups"] = len(
            set(pseudobulk.sample_groups)
        )
        run_info["stage_counts"]["pseudobulk_comparisons"] = n_pseudobulk_comparisons
        run_info["stage_counts"]["pseudobulk_marker_tests"] = (
            n_pseudobulk_comparisons * n_pseudobulk_genes
        )
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
    if stability_result is not None:
        # 聚类稳定性汇总：仅显式启用时出现；既有字段不变
        stats = summarize_scores(stability_result.scores)
        run_info["stability_analysis"] = {
            "n_samples": stability_result.config.n_samples,
            "sample_fraction": stability_result.config.sample_fraction,
            "seed": stability_result.config.seed,
            "n_available_cells": stability_result.n_available_cells,
            "sample_size": stability_result.sample_sizes[0],
            "mean_adjusted_rand_index": stats["mean"],
            "median_adjusted_rand_index": stats["median"],
            "min_adjusted_rand_index": stats["min"],
            "max_adjusted_rand_index": stats["max"],
        }
    if cell_type_annotations is not None:
        # 自动细胞类型注释汇总：各类型已注释簇数与最高分不大于 0 而未注释
        # （获胜候选仍为该类型）的簇数；未提供参考时顶层字段不出现
        annotation_summary: Dict[str, Any] = {
            cell_type: {"assigned_clusters": 0, "unassigned_clusters": 0}
            for cell_type in cell_type_annotations.type_order
        }
        for annotation in cell_type_annotations.clusters:
            bucket = annotation_summary[annotation.cell_type]
            if annotation.annotation_status == "assigned":
                bucket["assigned_clusters"] += 1
            else:
                bucket["unassigned_clusters"] += 1
        run_info["annotation"] = annotation_summary
    if marker_enrichment is not None:
        # marker 基因集富集汇总：逐簇命中基因数之和（命中数），以及
        # （簇, 集合）组合中 p_value_adj <= alpha 的显著集合数
        run_info["enrichment"] = {
            "marker_hits": sum(marker_enrichment.cluster_marker_counts.values()),
            "significant_sets": sum(
                1 for r in marker_enrichment.rows if r.p_value_adj <= config.enrichment_alpha
            ),
        }
    if pseudobulk_gene_set_scores is not None:
        # pseudobulk 基因集分组差异汇总：集合数、样本数、比较数、检验数
        # （比较数 × 集合数）与全部检验中的最小校正 P 值；仅显式启用时出现
        n_pbgs_comparisons = len(pseudobulk_gene_set_scores.comparisons)
        n_pbgs_sets = len(pseudobulk_gene_set_scores.set_order)
        adjusted_values = [
            record.p_value_adj
            for _, _, _, records in pseudobulk_gene_set_scores.comparisons
            for record in records
        ]
        run_info["gene_set_de"] = {
            "set_count": n_pbgs_sets,
            "sample_count": len(pseudobulk.sample_ids),
            "comparison_count": n_pbgs_comparisons,
            "test_count": n_pbgs_comparisons * n_pbgs_sets,
            "min_p_value_adj": min(adjusted_values),
        }
    if differential_abundance is not None:
        # 簇级样本差异丰度汇总：样本数、簇数、比较数与全部检验中的最小
        # 校正 P 值；仅显式启用时出现，既有字段不变
        run_info["differential_abundance"] = {
            "sample_count": len(differential_abundance.sample_ids),
            "cluster_count": len(differential_abundance.clusters),
            "comparison_count": len(differential_abundance.comparisons),
            "min_p_value_adj": min(
                record.p_value_adj
                for _, _, _, records in differential_abundance.comparisons
                for record in records
            ),
        }
    if cluster_pseudobulk is not None:
        # 簇内 pseudobulk 差异表达汇总：已检验/跳过的簇、比较数、检验数与
        # 全部检验中的最小校正 P 值；仅显式启用时出现，既有字段不变
        run_info["cluster_pseudobulk_de"] = {
            "tested_clusters": list(cluster_pseudobulk.tested_clusters),
            "skipped_clusters": list(cluster_pseudobulk.skipped_clusters),
            "comparison_count": cluster_pseudobulk.comparison_count,
            "test_count": cluster_pseudobulk.test_count,
            "min_p_value_adj": cluster_pseudobulk.min_p_value_adj,
        }
    if paired_cluster_pseudobulk is not None:
        # 簇内配对 pseudobulk 差异表达汇总：a/b 分组、已检验/跳过的簇、
        # 各簇 pair 数、检验数与全部检验中的最小校正 P 值；仅显式启用时出现
        run_info["paired_cluster_pseudobulk_de"] = {
            "group_a": paired_cluster_pseudobulk.group_a,
            "group_b": paired_cluster_pseudobulk.group_b,
            "tested_clusters": list(paired_cluster_pseudobulk.tested_clusters),
            "skipped_clusters": list(paired_cluster_pseudobulk.skipped_clusters),
            "pair_counts": [
                pair_count for _, pair_count, _ in paired_cluster_pseudobulk.results
            ],
            "test_count": paired_cluster_pseudobulk.test_count,
            "min_p_value_adj": paired_cluster_pseudobulk.min_p_value_adj,
        }
    if adjusted_cluster_pseudobulk is not None:
        # 协变量校正的簇内 pseudobulk 差异表达汇总：参与校正的协变量、
        # 已检验/跳过的簇、比较数、检验数与全部检验中的最小校正 P 值；
        # 仅显式启用时出现，既有字段不变
        run_info["cluster_pseudobulk_adjusted_de"] = {
            "covariates": list(adjusted_cluster_pseudobulk.covariate_names),
            "tested_clusters": list(adjusted_cluster_pseudobulk.tested_clusters),
            "skipped_clusters": list(adjusted_cluster_pseudobulk.skipped_clusters),
            "comparison_count": adjusted_cluster_pseudobulk.comparison_count,
            "test_count": adjusted_cluster_pseudobulk.test_count,
            "min_p_value_adj": adjusted_cluster_pseudobulk.min_p_value_adj,
        }
    if numeric_cluster_pseudobulk is not None:
        # 含数值型协变量的簇内 pseudobulk 检验汇总：数值协变量名、
        # 已检验/跳过的簇、比较数、检验数（组别 + 各连续列 × 保留基因）与
        # 全部检验中的最小校正 P 值；仅显式提供数值协变量表时出现
        run_info["cluster_pseudobulk_numeric_covariate_de"] = {
            "covariates": list(numeric_cluster_pseudobulk.covariate_names),
            "tested_clusters": list(
                numeric_cluster_pseudobulk.tested_clusters
            ),
            "skipped_clusters": list(
                numeric_cluster_pseudobulk.skipped_clusters
            ),
            "comparison_count": numeric_cluster_pseudobulk.comparison_count,
            "test_count": numeric_cluster_pseudobulk.test_count,
            "min_p_value_adj": numeric_cluster_pseudobulk.min_p_value_adj,
        }

    artifacts = Artifacts(
        qc=qc,
        data=normalized,
        pca=pca,
        clustering=clustering,
        markers=markers,
        pairwise_markers=pairwise_markers,
        run_info=run_info,
        pca_loadings=pca_loadings,
        cluster_selection=selection,
        group_markers=group_markers,
        batch_summary=batch_summary,
        gene_set_scores=gene_set_scores,
        cell_batch_report=cell_batch_report,
        pseudobulk=pseudobulk,
        doublets=doublets,
        stability=stability_result,
        cell_type_annotations=cell_type_annotations,
        marker_enrichment=marker_enrichment,
        pseudobulk_gene_set_scores=pseudobulk_gene_set_scores,
        differential_abundance=differential_abundance,
        cluster_pseudobulk=cluster_pseudobulk,
        paired_cluster_pseudobulk=paired_cluster_pseudobulk,
        adjusted_cluster_pseudobulk=adjusted_cluster_pseudobulk,
        numeric_cluster_pseudobulk=numeric_cluster_pseudobulk,
    )
    # 全部计算已完成才触碰文件系统：预检与写出都在 publish_results 内，
    # 任一分析阶段失败时不会创建或改动目标目录；写出阶段任何文件系统故障
    # 都回滚暂存并抛 OutputPathError。
    return publish_results(config.output_dir, artifacts)
