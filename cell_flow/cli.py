"""``cell-flow`` 命令行入口。

用法::

    cell-flow analyze --input <表达矩阵.tsv|mtx目录> --output-dir <结果目录> \\
        [--input-format tsv|mtx] [--metadata <细胞分组.tsv>] \\
        [--batch-metadata <细胞批次.tsv>] [--gene-sets <基因集.tsv>] \\
        [--cell-metadata <细胞元数据.tsv>] [--batch-column batch] \\
        [--sample-column sample_id] [--replicate-metadata <重复元数据.tsv>] \\
        [--min-genes 200] [--max-mito-fraction 0.2] [--min-cells 3] \\
        [--mito-prefix MT-] [--n-hvg 2000] [--n-clusters 2|auto] [--seed 20240617] \\
        [--detect-doublets] [--expected-doublet-rate 0.08] \\
        [--cell-type-reference <标记参考.tsv>] \\
        [--enrich-markers] [--enrichment-alpha F] [--enrichment-min-log-fc F] \\
        [--pseudobulk-gene-set-de] [--differential-abundance] \\
        [--cluster-pseudobulk-de] [--pca-loadings] \\
        [--stability-analysis] [--stability-n-samples 100] \\
        [--stability-sample-fraction 0.8] [--stability-seed 20240617]
"""

import sys
from typing import Optional, Sequence

from . import __version__
from .cell_metadata import DEFAULT_BATCH_COLUMN, DEFAULT_SAMPLE_COLUMN
from .errors import CellFlowError
from .pca import MAX_PCS
from .pipeline import DEFAULT_SEED, Config, INPUT_FORMATS, run
from .stability import (
    DEFAULT_STABILITY_FRACTION,
    DEFAULT_STABILITY_SAMPLES,
    DEFAULT_STABILITY_SEED,
)

USAGE = (
    "用法：cell-flow analyze --input <表达矩阵.tsv|mtx目录> --output-dir <结果目录>\n"
    "                [--input-format tsv|mtx] [--metadata <细胞分组.tsv>]\n"
    "                [--batch-metadata <细胞批次.tsv>] [--gene-sets <基因集.tsv>]\n"
    "                [--cell-metadata <细胞元数据.tsv>] [--batch-column 列名]\n"
    "                [--sample-column 列名] [--replicate-metadata <重复元数据.tsv>]\n"
    "                [--min-genes N] [--max-mito-fraction F] [--min-cells N]\n"
    "                [--mito-prefix PREFIX] [--n-hvg N] [--n-pcs N]\n"
    "                [--n-clusters N|auto] [--seed N]\n"
    "                [--detect-doublets] [--expected-doublet-rate F]\n"
    "                [--cell-type-reference <标记参考.tsv>]\n"
    "                [--enrich-markers] [--enrichment-alpha F]\n"
    "                [--enrichment-min-log-fc F]\n"
    "                [--pseudobulk-gene-set-de] [--differential-abundance]\n"
    "                [--cluster-pseudobulk-de] [--pca-loadings]\n"
    "                [--stability-analysis] [--stability-n-samples N]\n"
    "                [--stability-sample-fraction F] [--stability-seed N]\n"
    "      cell-flow --version"
)


def _fail_usage(message: str) -> "CellFlowError":
    # 延迟导入以避免循环引用
    from .errors import CellFlowConfigError

    return CellFlowConfigError(f"{message}\n{USAGE}")


def _parse_int(name: str, raw: str) -> int:
    try:
        value = int(raw, 10)
    except (TypeError, ValueError):
        raise _fail_usage(f"参数 {name} 需要整数，得到 {raw!r}")
    return value


def _parse_float(name: str, raw: str) -> float:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        raise _fail_usage(f"参数 {name} 需要数值，得到 {raw!r}")
    if value != value or value in (float("inf"), float("-inf")):
        raise _fail_usage(f"参数 {name} 必须是有限数值，得到 {raw!r}")
    return value


def _parse_n_clusters(raw: str):
    # auto（逐字小写）启用自动选择；其余取值仍按整数解析，<2 由 validate_config 拒绝
    if raw == "auto":
        return "auto"
    return _parse_int("--n-clusters", raw)


def _build_analyze_config(options: dict) -> Config:
    if options.get("input") is None:
        raise _fail_usage("缺少必需参数 --input")
    if options.get("output_dir") is None:
        raise _fail_usage("缺少必需参数 --output-dir")
    gene_sets_path = options.get("gene_sets")
    if gene_sets_path == "":
        # --gene-sets= 属调用不合法（配置错误，退出码 3）；
        # 内容问题（文件不存在、格式非法）在管线读入阶段按输入错误处理
        raise _fail_usage("--gene-sets 路径为空")
    cell_metadata_path = options.get("cell_metadata")
    if cell_metadata_path == "":
        # --cell-metadata= 同上：空路径属配置错误
        raise _fail_usage("--cell-metadata 路径为空")
    replicate_metadata_path = options.get("replicate_metadata")
    if replicate_metadata_path == "":
        # --replicate-metadata= 同上：空路径属配置错误
        raise _fail_usage("--replicate-metadata 路径为空")
    cell_type_reference_path = options.get("cell_type_reference")
    if cell_type_reference_path == "":
        # --cell-type-reference= 属调用不合法（配置错误，退出码 3）；
        # 内容问题（文件不存在、格式非法）在管线读入阶段按输入错误处理
        raise _fail_usage("--cell-type-reference 路径为空")
    return Config(
        input_path=options["input"],
        output_dir=options["output_dir"],
        min_genes=_parse_int("--min-genes", options["min_genes"]),
        max_mito_fraction=_parse_float(
            "--max-mito-fraction", options["max_mito_fraction"]
        ),
        min_cells=_parse_int("--min-cells", options["min_cells"]),
        mito_prefix=options["mito_prefix"],
        n_hvg=_parse_int("--n-hvg", options["n_hvg"]),
        n_pcs=_parse_int("--n-pcs", options["n_pcs"]),
        n_clusters=_parse_n_clusters(options["n_clusters"]),
        seed=_parse_int("--seed", options["seed"]),
        input_format=_parse_input_format(options["input_format"]),
        metadata_path=options.get("metadata"),
        batch_metadata_path=options.get("batch_metadata"),
        gene_sets_path=gene_sets_path,
        cell_metadata_path=cell_metadata_path,
        batch_column=options["batch_column"],
        sample_column=options["sample_column"],
        replicate_metadata_path=replicate_metadata_path,
        cell_type_reference_path=cell_type_reference_path,
        detect_doublets=bool(options.get("detect_doublets")),
        expected_doublet_rate=_parse_float(
            "--expected-doublet-rate", options["expected_doublet_rate"]
        ),
        stability_analysis=bool(options.get("stability_analysis")),
        stability_n_samples=_parse_int(
            "--stability-n-samples", options["stability_n_samples"]
        ),
        stability_sample_fraction=_parse_float(
            "--stability-sample-fraction", options["stability_sample_fraction"]
        ),
        stability_seed=_parse_int(
            "--stability-seed", options["stability_seed"]
        ),
        enrich_markers=bool(options.get("enrich_markers")),
        enrichment_alpha=_parse_float(
            "--enrichment-alpha", options["enrichment_alpha"]
        ),
        enrichment_min_log_fc=_parse_float(
            "--enrichment-min-log-fc", options["enrichment_min_log_fc"]
        ),
        pseudobulk_gene_set_de=bool(options.get("pseudobulk_gene_set_de")),
        differential_abundance=bool(options.get("differential_abundance")),
        cluster_pseudobulk_de=bool(options.get("cluster_pseudobulk_de")),
        pca_loadings=bool(options.get("pca_loadings")),
    )


def _parse_input_format(raw: str) -> str:
    if raw not in INPUT_FORMATS:
        raise _fail_usage(
            f"--input-format 只能是 tsv 或 mtx，得到 {raw!r}"
        )
    return raw


def _parse_analyze(argv: Sequence[str]) -> Config:
    # argparse 的退出码是 2，与本工具“配置错误退出 3”冲突，故自行做确定解析
    value_options = {
        "--input": "input",
        "--output-dir": "output_dir",
        "--input-format": "input_format",
        "--metadata": "metadata",
        "--batch-metadata": "batch_metadata",
        "--gene-sets": "gene_sets",
        "--cell-metadata": "cell_metadata",
        "--batch-column": "batch_column",
        "--sample-column": "sample_column",
        "--replicate-metadata": "replicate_metadata",
        "--cell-type-reference": "cell_type_reference",
        "--enrichment-alpha": "enrichment_alpha",
        "--enrichment-min-log-fc": "enrichment_min_log_fc",
        "--min-genes": "min_genes",
        "--max-mito-fraction": "max_mito_fraction",
        "--min-cells": "min_cells",
        "--mito-prefix": "mito_prefix",
        "--n-hvg": "n_hvg",
        "--n-pcs": "n_pcs",
        "--n-clusters": "n_clusters",
        "--seed": "seed",
        "--expected-doublet-rate": "expected_doublet_rate",
        "--stability-n-samples": "stability_n_samples",
        "--stability-sample-fraction": "stability_sample_fraction",
        "--stability-seed": "stability_seed",
    }
    # 无值开关参数；--detect-doublets=x 形式按未知参数拒绝
    flag_options = {
        "--detect-doublets": "detect_doublets",
        "--stability-analysis": "stability_analysis",
        "--enrich-markers": "enrich_markers",
        "--pseudobulk-gene-set-de": "pseudobulk_gene_set_de",
        "--differential-abundance": "differential_abundance",
        "--cluster-pseudobulk-de": "cluster_pseudobulk_de",
        "--pca-loadings": "pca_loadings",
    }
    options = {
        "input_format": "tsv",
        "min_genes": "200",
        "max_mito_fraction": "0.2",
        "min_cells": "3",
        "mito_prefix": "MT-",
        "n_hvg": "2000",
        "n_pcs": str(MAX_PCS),
        "n_clusters": "2",
        "seed": str(DEFAULT_SEED),
        "batch_column": DEFAULT_BATCH_COLUMN,
        "sample_column": DEFAULT_SAMPLE_COLUMN,
        "expected_doublet_rate": "0.08",
        "stability_n_samples": str(DEFAULT_STABILITY_SAMPLES),
        "stability_sample_fraction": str(DEFAULT_STABILITY_FRACTION),
        "stability_seed": str(DEFAULT_STABILITY_SEED),
        "enrichment_alpha": "0.05",
        "enrichment_min_log_fc": "0",
    }

    explicit = set()
    seen_flags = set()
    index = 0
    while index < len(argv):
        token = argv[index]
        if token.startswith("--") and "=" in token:
            key, raw_value = token.split("=", 1)
            if key not in value_options:
                raise _fail_usage(f"未知参数 {key}")
            options[value_options[key]] = raw_value
            explicit.add(value_options[key])
            index += 1
            continue
        if token in value_options:
            if index + 1 >= len(argv):
                raise _fail_usage(f"参数 {token} 缺少取值")
            options[value_options[token]] = argv[index + 1]
            explicit.add(value_options[token])
            index += 2
            continue
        if token in flag_options:
            flag_key = flag_options[token]
            if flag_key == "pca_loadings" and flag_key in seen_flags:
                # --pca-loadings 只接受无值形式且不可重复；重复出现属配置错误
                raise _fail_usage(f"参数 {token} 重复出现")
            seen_flags.add(flag_key)
            options[flag_key] = True
            index += 1
            continue
        raise _fail_usage(f"无法识别的参数或多余位置参数：{token}")

    if "expected_doublet_rate" in explicit and not options.get("detect_doublets"):
        # 只给 --expected-doublet-rate 而不启用双细胞识别属配置错误
        raise _fail_usage("--expected-doublet-rate 需与 --detect-doublets 同时使用")

    stability_param_keys = (
        "stability_n_samples",
        "stability_sample_fraction",
        "stability_seed",
    )
    if any(key in explicit for key in stability_param_keys) and not options.get(
        "stability_analysis"
    ):
        # 只给稳定性参数而不启用 --stability-analysis 属配置错误
        raise _fail_usage(
            "--stability-n-samples/--stability-sample-fraction/--stability-seed "
            "需与 --stability-analysis 同时使用"
        )

    if options.get("enrich_markers") and options.get("gene_sets") is None:
        # marker 基因集富集只在 --gene-sets 评分基线上启用
        raise _fail_usage("--enrich-markers 需与 --gene-sets 同时使用")

    if options.get("pseudobulk_gene_set_de") and (
        options.get("gene_sets") is None
        or options.get("replicate_metadata") is None
    ):
        # pseudobulk 基因集分组差异只在 --gene-sets 与 --replicate-metadata
        # 同时提供时可用；无值开关带值（--flag=x）已在上方按未知参数拒绝
        raise _fail_usage(
            "--pseudobulk-gene-set-de 需与 --gene-sets、"
            "--replicate-metadata 同时使用"
        )

    if options.get("differential_abundance") and options.get(
        "replicate_metadata"
    ) is None:
        # 簇级样本差异丰度只在 --replicate-metadata 样本分组基线上启用；
        # 无值开关带值（--flag=x）已在上方按未知参数拒绝
        raise _fail_usage(
            "--differential-abundance 需与 --replicate-metadata 同时使用"
        )

    if options.get("cluster_pseudobulk_de") and options.get(
        "replicate_metadata"
    ) is None:
        # 簇内 pseudobulk 差异表达只在 --replicate-metadata 样本分组基线上启用；
        # 无值开关带值（--flag=x）已在上方按未知参数拒绝
        raise _fail_usage(
            "--cluster-pseudobulk-de 需与 --replicate-metadata 同时使用"
        )

    enrichment_param_keys = ("enrichment_alpha", "enrichment_min_log_fc")
    if any(key in explicit for key in enrichment_param_keys) and not options.get(
        "enrich_markers"
    ):
        # 只给富集参数而不启用 --enrich-markers 属配置错误
        raise _fail_usage(
            "--enrichment-alpha/--enrichment-min-log-fc 需与 --enrich-markers "
            "同时使用"
        )

    return _build_analyze_config(options)


def main(argv: Optional[Sequence[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        message = _fail_usage("缺少子命令")
        print(f"cell-flow: 错误: {message}", file=sys.stderr)
        return message.exit_code
    if argv[0] in ("-h", "--help"):
        print(USAGE)
        return 0
    if argv[0] in ("-V", "--version"):
        print(f"cell-flow {__version__}")
        return 0
    if argv[0] != "analyze":
        message = _fail_usage(f"未知子命令 {argv[0]!r}")
        print(f"cell-flow: 错误: {message}", file=sys.stderr)
        return message.exit_code

    if any(token in ("-h", "--help") for token in argv[1:]):
        print(USAGE)
        return 0

    try:
        config = _parse_analyze(argv[1:])
        written = run(config)
    except CellFlowError as exc:
        print(f"cell-flow: 错误: {exc}", file=sys.stderr)
        return exc.exit_code

    print(
        f"cell-flow: 分析完成，写出 {len(written)} 个文件至 {config.output_dir}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
