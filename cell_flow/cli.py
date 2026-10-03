"""``cell-flow`` 命令行入口。

用法::

    cell-flow analyze --input <表达矩阵.tsv|mtx目录> --output-dir <结果目录> \\
        [--input-format tsv|mtx] [--metadata <细胞分组.tsv>] \\
        [--batch-metadata <细胞批次.tsv>] [--gene-sets <基因集.tsv>] \\
        [--min-genes 200] [--max-mito-fraction 0.2] [--min-cells 3] \\
        [--mito-prefix MT-] [--n-hvg 2000] [--n-clusters 2|auto] [--seed 20240617]
"""

import sys
from typing import Optional, Sequence

from . import __version__
from .errors import CellFlowError
from .pca import MAX_PCS
from .pipeline import DEFAULT_SEED, Config, INPUT_FORMATS, run

USAGE = (
    "用法：cell-flow analyze --input <表达矩阵.tsv|mtx目录> --output-dir <结果目录>\n"
    "                [--input-format tsv|mtx] [--metadata <细胞分组.tsv>]\n"
    "                [--batch-metadata <细胞批次.tsv>] [--gene-sets <基因集.tsv>]\n"
    "                [--min-genes N] [--max-mito-fraction F] [--min-cells N]\n"
    "                [--mito-prefix PREFIX] [--n-hvg N] [--n-pcs N]\n"
    "                [--n-clusters N|auto] [--seed N]\n"
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
        "--min-genes": "min_genes",
        "--max-mito-fraction": "max_mito_fraction",
        "--min-cells": "min_cells",
        "--mito-prefix": "mito_prefix",
        "--n-hvg": "n_hvg",
        "--n-pcs": "n_pcs",
        "--n-clusters": "n_clusters",
        "--seed": "seed",
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
    }

    index = 0
    while index < len(argv):
        token = argv[index]
        if token.startswith("--") and "=" in token:
            key, raw_value = token.split("=", 1)
            if key not in value_options:
                raise _fail_usage(f"未知参数 {key}")
            options[value_options[key]] = raw_value
            index += 1
            continue
        if token in value_options:
            if index + 1 >= len(argv):
                raise _fail_usage(f"参数 {token} 缺少取值")
            options[value_options[token]] = argv[index + 1]
            index += 2
            continue
        raise _fail_usage(f"无法识别的参数或多余位置参数：{token}")

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
