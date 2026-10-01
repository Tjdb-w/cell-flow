"""Analysis pipeline: QC, normalization, HVG, PCA, clustering, markers."""

import math

from .errors import CellFlowDataError
from .stats import bh_adjust, kmeans, pca_scores, welch_t_pvalue

SCALE_FACTOR = 10000.0
MAX_PCS = 20
TOP_MARKERS_PER_CLUSTER = 20


class Results:
    """Everything needed to render the output files."""

    def __init__(self):
        self.gene_ids = []
        self.cell_ids = []
        self.cell_total_counts = []
        self.cell_n_genes = []
        self.cell_mito_counts = []
        self.cell_mito_fraction = []
        self.cell_kept = []
        self.gene_total_counts = []
        self.gene_n_cells = []
        self.gene_kept = []
        self.kept_cell_ids = []
        self.kept_gene_ids = []
        self.n_hvg = 0
        self.n_components = 0
        self.pca_scores = []
        self.cluster_labels = []
        self.cluster_sizes = []
        self.markers = []          # list of per-cluster sorted row lists
        self.top_markers = []      # list of (cluster, rank, gene_id, log_fc, p_adj)
        self.top_expression = []   # list of (gene_id, cell_id, cluster, value)
        self.n_tests = 0


def compute_qc(gene_ids, cell_ids, counts, mito_prefix, min_genes,
               max_mito, min_cells, res):
    """QC metrics from raw counts; fill keep flags on ``res``."""
    n_genes = len(gene_ids)
    n_cells = len(cell_ids)
    res.cell_total_counts = [0] * n_cells
    res.cell_n_genes = [0] * n_cells
    res.cell_mito_counts = [0] * n_cells
    res.gene_total_counts = [0] * n_genes
    res.gene_n_cells = [0] * n_genes
    for g in range(n_genes):
        is_mito = bool(mito_prefix) and gene_ids[g].startswith(mito_prefix)
        row = counts[g]
        for c in range(n_cells):
            value = row[c]
            res.cell_total_counts[c] += value
            res.gene_total_counts[g] += value
            if value > 0:
                res.cell_n_genes[c] += 1
                res.gene_n_cells[g] += 1
                if is_mito:
                    res.cell_mito_counts[c] += value
    res.cell_mito_fraction = [
        (res.cell_mito_counts[c] / res.cell_total_counts[c]
         if res.cell_total_counts[c] > 0 else 0.0)
        for c in range(n_cells)
    ]
    res.cell_kept = [
        res.cell_n_genes[c] >= min_genes
        and res.cell_mito_fraction[c] <= max_mito
        for c in range(n_cells)
    ]
    res.gene_kept = [res.gene_n_cells[g] >= min_cells for g in range(n_genes)]


def normalize(counts, kept_genes, kept_cells, totals):
    """Log-normalized matrix (kept genes x kept cells).

    Each cell is scaled to SCALE_FACTOR total counts, then ln(1 + x).
    """
    norm = []
    for g in kept_genes:
        row = counts[g]
        out = []
        for ci, c in enumerate(kept_cells):
            total = totals[c]
            if total > 0:
                out.append(math.log1p(row[c] / total * SCALE_FACTOR))
            else:
                out.append(0.0)
        norm.append(out)
    return norm


def select_hvgs(norm, n_hvg):
    """Indices (into kept genes) of the most variable genes.

    Ranked by population variance of log-normalized values, descending;
    ties broken by kept-gene order.
    """
    n_cells = len(norm[0])
    variances = []
    for row in norm:
        mean = math.fsum(row) / n_cells
        var = math.fsum((x - mean) ** 2 for x in row) / n_cells
        variances.append(var)
    order = sorted(range(len(norm)), key=lambda g: (-variances[g], g))
    return sorted(order[:n_hvg])


def differential_expression(norm, kept_gene_ids, labels, n_clusters):
    """Per-cluster vs rest markers on log-normalized values.

    Returns a list (per cluster) of rows
    [cluster, gene_id, mean_cluster, mean_rest, log_fc, p_value, p_adj]
    sorted by adjusted p-value, descending log fold change, gene ID.
    """
    n_cells = len(labels)
    per_cluster = []
    for cluster in range(n_clusters):
        in_idx = [i for i in range(n_cells) if labels[i] == cluster]
        out_idx = [i for i in range(n_cells) if labels[i] != cluster]
        rows = []
        for g, gene_id in enumerate(kept_gene_ids):
            row = norm[g]
            xs = [row[i] for i in in_idx]
            ys = [row[i] for i in out_idx]
            mean_in = math.fsum(xs) / len(xs)
            mean_out = math.fsum(ys) / len(ys)
            p_value = welch_t_pvalue(xs, ys)
            rows.append([cluster, gene_id, mean_in, mean_out,
                         mean_in - mean_out, p_value])
        adjusted = bh_adjust([r[5] for r in rows])
        for r, q in zip(rows, adjusted):
            r.append(q)
        rows.sort(key=lambda r: (r[6], -r[4], r[1]))
        per_cluster.append(rows)
    return per_cluster


def run_pipeline(config, gene_ids, cell_ids, counts):
    """Execute the full analysis. Raises CellFlowDataError on failure."""
    res = Results()
    res.gene_ids = gene_ids
    res.cell_ids = cell_ids

    compute_qc(gene_ids, cell_ids, counts, config["mito_prefix"],
               config["min_genes"], config["max_mito"],
               config["min_cells"], res)

    kept_cells = [c for c in range(len(cell_ids)) if res.cell_kept[c]]
    kept_genes = [g for g in range(len(gene_ids)) if res.gene_kept[g]]
    if len(kept_cells) < 2:
        raise CellFlowDataError(
            "fewer than two cells remain after quality control "
            f"({len(kept_cells)} kept)"
        )
    if not kept_genes:
        raise CellFlowDataError(
            "no genes remain after quality control"
        )
    res.kept_cell_ids = [cell_ids[c] for c in kept_cells]
    res.kept_gene_ids = [gene_ids[g] for g in kept_genes]

    norm = normalize(counts, kept_genes, kept_cells, res.cell_total_counts)

    n_hvg = min(config["n_hvg"], len(kept_genes))
    hvg_idx = select_hvgs(norm, n_hvg)
    res.n_hvg = len(hvg_idx)

    n_components = min(MAX_PCS, len(kept_cells) - 1, res.n_hvg)
    if n_components < 1:
        raise CellFlowDataError(
            "PCA cannot be computed: no components available"
        )
    res.n_components = n_components
    hvg_rows = [[norm[g][c] for g in hvg_idx] for c in range(len(kept_cells))]
    res.pca_scores = pca_scores(hvg_rows, n_components)

    n_clusters = config["n_clusters"]
    if n_clusters > len(kept_cells):
        raise CellFlowDataError(
            f"cannot form {n_clusters} clusters from "
            f"{len(kept_cells)} cells"
        )
    res.cluster_labels = kmeans(res.pca_scores, n_clusters, config["seed"])
    res.cluster_sizes = [
        sum(1 for label in res.cluster_labels if label == c)
        for c in range(n_clusters)
    ]

    res.markers = differential_expression(
        norm, res.kept_gene_ids, res.cluster_labels, n_clusters
    )
    res.n_tests = n_clusters * len(res.kept_gene_ids)

    # Top markers per cluster and their expression values.
    seen = set()
    top_gene_ids = []
    for cluster, rows in enumerate(res.markers):
        for rank, r in enumerate(rows[:TOP_MARKERS_PER_CLUSTER], start=1):
            res.top_markers.append((cluster, rank, r[1], r[4], r[6]))
            if r[1] not in seen:
                seen.add(r[1])
                top_gene_ids.append(r[1])
    gene_row = {gene_id: i for i, gene_id in enumerate(res.kept_gene_ids)}
    for gene_id in top_gene_ids:
        g = gene_row[gene_id]
        for ci, cell_id in enumerate(res.kept_cell_ids):
            res.top_expression.append(
                (gene_id, cell_id, res.cluster_labels[ci], norm[g][ci])
            )

    return res
