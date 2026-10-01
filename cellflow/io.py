"""Input parsing and deterministic output rendering."""

import hashlib
import json
import os
import re

from . import VERSION
from .errors import CellFlowInputError, OutputPathError

_INT_RE = re.compile(r"[0-9]+")
_NEG_RE = re.compile(r"-[0-9]+")


def fmt6(value):
    """Fixed 6-decimal formatting without negative zero."""
    if value == 0.0:
        value = 0.0
    text = f"{value:.6f}"
    return "0.000000" if text == "-0.000000" else text


def fmt_p(value):
    """Scientific formatting for p-values."""
    if value == 0.0:
        value = 0.0
    return f"{value:.6e}"


def read_matrix(path):
    """Read a TSV expression matrix.

    Returns (gene_ids, cell_ids, counts, sha256_hex) where counts[g][c]
    is the UMI count of gene g in cell c. Raises CellFlowInputError.
    """
    if path == "":
        raise CellFlowInputError("input path is empty")
    if not os.path.exists(path):
        raise CellFlowInputError(f"input file not found: {path}")
    if not os.path.isfile(path):
        raise CellFlowInputError(f"input path is not a file: {path}")
    try:
        with open(path, "rb") as handle:
            raw = handle.read()
    except OSError as exc:
        raise CellFlowInputError(f"cannot read input file: {path}") from exc
    digest = hashlib.sha256(raw).hexdigest()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CellFlowInputError(
            f"input file is not valid UTF-8: {path}"
        ) from exc
    if text.startswith("\ufeff"):
        text = text[1:]

    lines = [line.rstrip("\r") for line in text.split("\n")]
    lines = [line for line in lines if line != ""]
    if not lines:
        raise CellFlowInputError("input file has no header row")

    header = lines[0].split("\t")
    if len(header) < 2:
        raise CellFlowInputError("header row has no cell columns")
    cell_ids = header[1:]
    seen = set()
    for index, cell_id in enumerate(cell_ids):
        if cell_id == "":
            raise CellFlowInputError(
                f"empty cell ID at header column {index + 2}"
            )
        if cell_id in seen:
            raise CellFlowInputError(f"duplicate cell ID: {cell_id}")
        seen.add(cell_id)

    gene_ids = []
    counts = []
    seen_genes = set()
    width = len(header)
    for row_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != width:
            raise CellFlowInputError(
                f"row {row_number} has {len(fields)} fields, "
                f"expected {width}"
            )
        gene_id = fields[0]
        if gene_id == "":
            raise CellFlowInputError(f"empty gene ID at row {row_number}")
        if gene_id in seen_genes:
            raise CellFlowInputError(f"duplicate gene ID: {gene_id}")
        seen_genes.add(gene_id)
        row = []
        for col, field in enumerate(fields[1:]):
            if _INT_RE.fullmatch(field):
                row.append(int(field))
            elif _NEG_RE.fullmatch(field):
                raise CellFlowInputError(
                    f"negative count at row {row_number}, "
                    f"column {col + 2}: {field}"
                )
            else:
                raise CellFlowInputError(
                    f"non-integer count at row {row_number}, "
                    f"column {col + 2}: {field}"
                )
        gene_ids.append(gene_id)
        counts.append(row)

    if not gene_ids:
        raise CellFlowInputError("input matrix has no gene rows")

    return gene_ids, cell_ids, counts, digest


def check_output_dir(path):
    """Validate the output directory before any computation."""
    if path == "":
        raise OutputPathError("output path is empty")
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise OutputPathError(
                f"output path exists and is not a directory: {path}"
            )
        if os.listdir(path):
            raise OutputPathError(
                f"output directory exists and is not empty: {path}"
            )


def _write(path, text):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _table(header, rows):
    out = ["\t".join(header)]
    for row in rows:
        out.append("\t".join(row))
    return "\n".join(out) + "\n"


def write_outputs(outdir, config, input_info, res):
    """Write every result file. All content is computed before this call."""
    os.makedirs(outdir, exist_ok=True)

    run = {
        "version": VERSION,
        "seed": config["seed"],
        "input": {
            "path": input_info["path"],
            "sha256": input_info["sha256"],
            "n_genes": len(res.gene_ids),
            "n_cells": len(res.cell_ids),
        },
        "params": {
            "min_genes": config["min_genes"],
            "max_mito": config["max_mito"],
            "min_cells": config["min_cells"],
            "mito_prefix": config["mito_prefix"],
            "n_hvg": config["n_hvg"],
            "n_clusters": config["n_clusters"],
            "seed": config["seed"],
        },
        "stages": {
            "qc": {
                "cells_before": len(res.cell_ids),
                "cells_after": len(res.kept_cell_ids),
                "genes_before": len(res.gene_ids),
                "genes_after": len(res.kept_gene_ids),
            },
            "highly_variable_genes": {"selected": res.n_hvg},
            "pca": {"n_components": res.n_components},
            "clustering": {
                "n_clusters": config["n_clusters"],
                "cluster_sizes": res.cluster_sizes,
            },
            "markers": {"n_tests": res.n_tests},
        },
    }
    _write(os.path.join(outdir, "run.json"),
           json.dumps(run, indent=2) + "\n")

    _write(os.path.join(outdir, "cells.tsv"), _table(
        ["cell_id", "total_counts", "n_genes", "mito_counts",
         "mito_fraction", "kept"],
        ([cell_id, str(res.cell_total_counts[i]), str(res.cell_n_genes[i]),
          str(res.cell_mito_counts[i]), fmt6(res.cell_mito_fraction[i]),
          "1" if res.cell_kept[i] else "0"]
         for i, cell_id in enumerate(res.cell_ids)),
    ))

    _write(os.path.join(outdir, "genes.tsv"), _table(
        ["gene_id", "total_counts", "n_cells", "kept"],
        ([gene_id, str(res.gene_total_counts[i]), str(res.gene_n_cells[i]),
          "1" if res.gene_kept[i] else "0"]
         for i, gene_id in enumerate(res.gene_ids)),
    ))

    pc_header = ["cell_id"] + [f"PC{k + 1}" for k in range(res.n_components)]
    _write(os.path.join(outdir, "pca.tsv"), _table(
        pc_header,
        ([cell_id] + [fmt6(v) for v in res.pca_scores[i]]
         for i, cell_id in enumerate(res.kept_cell_ids)),
    ))

    _write(os.path.join(outdir, "clusters.tsv"), _table(
        ["cell_id", "cluster"],
        ([cell_id, str(res.cluster_labels[i])]
         for i, cell_id in enumerate(res.kept_cell_ids)),
    ))

    marker_rows = []
    for rows in res.markers:
        for r in rows:
            marker_rows.append([
                str(r[0]), r[1], fmt6(r[2]), fmt6(r[3]), fmt6(r[4]),
                fmt_p(r[5]), fmt_p(r[6]),
            ])
    _write(os.path.join(outdir, "markers.tsv"), _table(
        ["cluster", "gene_id", "mean_cluster", "mean_rest", "log_fc",
         "p_value", "p_adj"],
        marker_rows,
    ))

    _write(os.path.join(outdir, "qc.tsv"), _table(
        ["cell_id", "total_counts", "n_genes", "mito_fraction", "kept"],
        ([cell_id, str(res.cell_total_counts[i]), str(res.cell_n_genes[i]),
          fmt6(res.cell_mito_fraction[i]),
          "1" if res.cell_kept[i] else "0"]
         for i, cell_id in enumerate(res.cell_ids)),
    ))

    _write(os.path.join(outdir, "pca_scatter.tsv"), _table(
        ["cell_id", "PC1", "PC2", "cluster"],
        ([cell_id,
          fmt6(res.pca_scores[i][0]),
          fmt6(res.pca_scores[i][1] if res.n_components > 1 else 0.0),
          str(res.cluster_labels[i])]
         for i, cell_id in enumerate(res.kept_cell_ids)),
    ))

    _write(os.path.join(outdir, "top_markers.tsv"), _table(
        ["cluster", "rank", "gene_id", "log_fc", "p_adj"],
        ([str(cluster), str(rank), gene_id, fmt6(log_fc), fmt_p(p_adj)]
         for cluster, rank, gene_id, log_fc, p_adj in res.top_markers),
    ))

    _write(os.path.join(outdir, "top_marker_expression.tsv"), _table(
        ["gene_id", "cell_id", "cluster", "expression"],
        ([gene_id, cell_id, str(cluster), fmt6(value)]
         for gene_id, cell_id, cluster, value in res.top_expression),
    ))
