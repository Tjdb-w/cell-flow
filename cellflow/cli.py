"""Command line interface: cell-flow analyze."""

import argparse
import math
import sys

from .errors import CellFlowConfigError, CellFlowError
from .io import check_output_dir, read_matrix, write_outputs
from .pipeline import run_pipeline


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        raise CellFlowConfigError(f"invalid arguments: {message}")


def _build_parser():
    parser = _Parser(
        prog="cell-flow",
        description="Deterministic single-cell analysis pipeline.",
    )
    sub = parser.add_subparsers(dest="command")
    analyze = sub.add_parser(
        "analyze",
        help="run the full analysis on a TSV expression matrix",
    )
    analyze.add_argument("--input", required=True,
                         help="path to the TSV expression matrix")
    analyze.add_argument("--output-dir", required=True,
                         help="directory for result files")
    analyze.add_argument("--min-genes", default="200",
                         help="minimum detected genes per cell (default 200)")
    analyze.add_argument("--max-mito", default="0.2",
                         help="maximum mitochondrial fraction per cell "
                              "(default 0.2)")
    analyze.add_argument("--min-cells", default="3",
                         help="minimum cells detecting a gene (default 3)")
    analyze.add_argument("--mito-prefix", default="MT-",
                         help="gene ID prefix marking mitochondrial genes "
                              "(default MT-)")
    analyze.add_argument("--n-hvg", default="2000",
                         help="number of highly variable genes "
                              "(default 2000)")
    analyze.add_argument("--n-clusters", default="2",
                         help="number of k-means clusters (default 2)")
    analyze.add_argument("--seed", default="0",
                         help="random seed for clustering (default 0)")
    return parser


def _parse_int(name, raw, minimum):
    try:
        value = int(raw)
    except ValueError:
        raise CellFlowConfigError(
            f"{name} must be an integer, got {raw!r}"
        ) from None
    if value < minimum:
        raise CellFlowConfigError(
            f"{name} must be at least {minimum}, got {value}"
        )
    return value


def _validate_config(args):
    config = {}
    config["input"] = args.input
    config["output_dir"] = args.output_dir
    config["min_genes"] = _parse_int("min-genes", args.min_genes, 0)
    config["min_cells"] = _parse_int("min-cells", args.min_cells, 1)
    config["n_hvg"] = _parse_int("n-hvg", args.n_hvg, 1)
    config["n_clusters"] = _parse_int("n-clusters", args.n_clusters, 2)
    config["seed"] = _parse_int("seed", args.seed, 0)
    try:
        max_mito = float(args.max_mito)
    except ValueError:
        raise CellFlowConfigError(
            f"max-mito must be a number, got {args.max_mito!r}"
        ) from None
    if not math.isfinite(max_mito) or max_mito < 0.0 or max_mito > 1.0:
        raise CellFlowConfigError(
            f"max-mito must be between 0 and 1, got {args.max_mito!r}"
        )
    config["max_mito"] = max_mito
    config["mito_prefix"] = args.mito_prefix
    return config


def _analyze(args):
    config = _validate_config(args)
    gene_ids, cell_ids, counts, digest = read_matrix(config["input"])
    check_output_dir(config["output_dir"])
    res = run_pipeline(config, gene_ids, cell_ids, counts)
    input_info = {"path": config["input"], "sha256": digest}
    write_outputs(config["output_dir"], config, input_info, res)
    print(
        f"cell-flow: analysis complete: {len(res.kept_cell_ids)} cells, "
        f"{len(res.kept_gene_ids)} genes, {config['n_clusters']} clusters "
        f"-> {config['output_dir']}"
    )


def main(argv=None):
    parser = _build_parser()
    try:
        args = parser.parse_args(argv)
        if args.command is None:
            parser.print_help()
            return 0
        if args.command == "analyze":
            _analyze(args)
            return 0
        raise CellFlowConfigError(f"unknown command: {args.command}")
    except CellFlowError as exc:
        print(f"cell-flow: error: {exc}", file=sys.stderr)
        return exc.exit_code


if __name__ == "__main__":
    sys.exit(main())
