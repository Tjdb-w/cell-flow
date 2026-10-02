"""逐簇差异表达：每个簇对比其余所有细胞，以及任意两个已得到簇之间的成对比较。

对每个（簇, 基因）输出：
- 两组在 log-归一化表达上的平均表达量；
- 对数倍数变化（簇内均值减其余均值）；
- Welch t 检验统计量与双侧 P 值；
- 每个簇内跨基因的 Benjamini-Hochberg 校正 P 值。
排序：校正 P 值升序、对数倍数变化降序、基因 ID 升序。

成对比较对实际出现的每一对簇 (a, b)（a < b，case 为 a、control 为 b）
执行同一口径的逐基因检验，并在每个比较内部跨基因做 BH 校正。
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple

from .linalg import benjamini_hochberg, welch_ttest
from .normalize import NormalizedData


@dataclass(frozen=True)
class MarkerRecord:
    cluster: int
    gene_id: str
    mean_in_cluster: float
    mean_out_cluster: float
    log_fc: float
    t_stat: float
    p_value: float
    p_value_adj: float


def find_markers(
    data: NormalizedData, labels: List[int]
) -> Dict[int, List[MarkerRecord]]:
    n_cells = len(data.cell_ids)
    cluster_to_cells: Dict[int, List[int]] = {}
    for c, label in enumerate(labels):
        cluster_to_cells.setdefault(label, []).append(c)

    results: Dict[int, List[MarkerRecord]] = {}
    for cluster in sorted(cluster_to_cells):
        in_cells = cluster_to_cells[cluster]
        out_cells = [c for c in range(n_cells) if labels[c] != cluster]
        if not in_cells or not out_cells:
            # 理论上不会发生（簇均非空且簇数 >= 2）
            continue

        records: List[MarkerRecord] = []
        pvalues: List[float] = []
        for g, gene_id in enumerate(data.gene_ids):
            row = data.analysis_values[g]
            group_in = [row[c] for c in in_cells]
            group_out = [row[c] for c in out_cells]
            mean_in = sum(group_in) / len(group_in)
            mean_out = sum(group_out) / len(group_out)
            log_fc = mean_in - mean_out
            t_stat, p_value = welch_ttest(group_in, group_out)
            pvalues.append(p_value)
            records.append(
                MarkerRecord(
                    cluster=cluster,
                    gene_id=gene_id,
                    mean_in_cluster=mean_in,
                    mean_out_cluster=mean_out,
                    log_fc=log_fc,
                    t_stat=t_stat,
                    p_value=p_value,
                    p_value_adj=0.0,
                )
            )

        adjusted = benjamini_hochberg(pvalues)
        records = [
            MarkerRecord(
                cluster=r.cluster,
                gene_id=r.gene_id,
                mean_in_cluster=r.mean_in_cluster,
                mean_out_cluster=r.mean_out_cluster,
                log_fc=r.log_fc,
                t_stat=r.t_stat,
                p_value=r.p_value,
                p_value_adj=adjusted[i],
            )
            for i, r in enumerate(records)
        ]
        records.sort(key=lambda r: (r.p_value_adj, -r.log_fc, r.gene_id))
        results[cluster] = records
    return results


@dataclass(frozen=True)
class PairwiseMarkerRecord:
    cluster_a: int
    cluster_b: int
    gene_id: str
    mean_in_a: float
    mean_in_b: float
    log_fc_a_vs_b: float
    t_stat: float
    p_value: float
    p_value_adj: float


def find_pairwise_markers(
    data: NormalizedData, labels: List[int]
) -> List[Tuple[int, int, List[PairwiseMarkerRecord]]]:
    """每一对实际出现的簇 (a, b)（a < b）之间的逐基因差异表达。

    返回按 (cluster_a, cluster_b) 升序的比较列表；每个比较内的记录按
    校正 P 值升序、log_fc_a_vs_b 降序、基因 ID 升序稳定排序。
    BH 校正以单个比较内的全部保留基因为一个校正家族。
    """
    cluster_to_cells: Dict[int, List[int]] = {}
    for c, label in enumerate(labels):
        cluster_to_cells.setdefault(label, []).append(c)

    present = sorted(cluster_to_cells)
    comparisons: List[Tuple[int, int, List[PairwiseMarkerRecord]]] = []
    for i, cluster_a in enumerate(present):
        for cluster_b in present[i + 1:]:
            cells_a = cluster_to_cells[cluster_a]
            cells_b = cluster_to_cells[cluster_b]

            records: List[PairwiseMarkerRecord] = []
            pvalues: List[float] = []
            for g, gene_id in enumerate(data.gene_ids):
                row = data.analysis_values[g]
                group_a = [row[c] for c in cells_a]
                group_b = [row[c] for c in cells_b]
                mean_a = sum(group_a) / len(group_a)
                mean_b = sum(group_b) / len(group_b)
                log_fc = mean_a - mean_b
                t_stat, p_value = welch_ttest(group_a, group_b)
                pvalues.append(p_value)
                records.append(
                    PairwiseMarkerRecord(
                        cluster_a=cluster_a,
                        cluster_b=cluster_b,
                        gene_id=gene_id,
                        mean_in_a=mean_a,
                        mean_in_b=mean_b,
                        log_fc_a_vs_b=log_fc,
                        t_stat=t_stat,
                        p_value=p_value,
                        p_value_adj=0.0,
                    )
                )

            adjusted = benjamini_hochberg(pvalues)
            records = [
                PairwiseMarkerRecord(
                    cluster_a=r.cluster_a,
                    cluster_b=r.cluster_b,
                    gene_id=r.gene_id,
                    mean_in_a=r.mean_in_a,
                    mean_in_b=r.mean_in_b,
                    log_fc_a_vs_b=r.log_fc_a_vs_b,
                    t_stat=r.t_stat,
                    p_value=r.p_value,
                    p_value_adj=adjusted[g],
                )
                for g, r in enumerate(records)
            ]
            records.sort(
                key=lambda r: (r.p_value_adj, -r.log_fc_a_vs_b, r.gene_id)
            )
            comparisons.append((cluster_a, cluster_b, records))
    return comparisons


ONE_VS_REST = "one_vs_rest"
PAIRWISE = "pairwise"


@dataclass(frozen=True)
class GroupMarkerRecord:
    comparison_type: str
    group_a: str
    group_b: str  # one-vs-rest 时为空
    gene_id: str
    mean_in_a: float
    mean_in_b: float
    log_fc_a_vs_b: float
    t_stat: float
    p_value: float
    p_value_adj: float


# (comparison_type, group_a, group_b, records)；one-vs-rest 的 group_b 为 ""
GroupComparison = Tuple[str, str, str, List[GroupMarkerRecord]]


def _group_comparison(
    data: NormalizedData,
    cells_a: List[int],
    cells_b: List[int],
    comparison_type: str,
    group_a: str,
    group_b: str,
) -> List[GroupMarkerRecord]:
    """两组保留细胞之间逐基因检验：均值差即 log_fc，Welch t 检验，
    并在本比较内跨基因做 BH 校正。"""
    records: List[GroupMarkerRecord] = []
    pvalues: List[float] = []
    for g, gene_id in enumerate(data.gene_ids):
        row = data.analysis_values[g]
        values_a = [row[c] for c in cells_a]
        values_b = [row[c] for c in cells_b]
        mean_a = sum(values_a) / len(values_a)
        mean_b = sum(values_b) / len(values_b)
        log_fc = mean_a - mean_b
        t_stat, p_value = welch_ttest(values_a, values_b)
        pvalues.append(p_value)
        records.append(
            GroupMarkerRecord(
                comparison_type=comparison_type,
                group_a=group_a,
                group_b=group_b,
                gene_id=gene_id,
                mean_in_a=mean_a,
                mean_in_b=mean_b,
                log_fc_a_vs_b=log_fc,
                t_stat=t_stat,
                p_value=p_value,
                p_value_adj=0.0,
            )
        )

    adjusted = benjamini_hochberg(pvalues)
    records = [
        GroupMarkerRecord(
            comparison_type=r.comparison_type,
            group_a=r.group_a,
            group_b=r.group_b,
            gene_id=r.gene_id,
            mean_in_a=r.mean_in_a,
            mean_in_b=r.mean_in_b,
            log_fc_a_vs_b=r.log_fc_a_vs_b,
            t_stat=r.t_stat,
            p_value=r.p_value,
            p_value_adj=adjusted[g],
        )
        for g, r in enumerate(records)
    ]
    records.sort(key=lambda r: (r.p_value_adj, -r.log_fc_a_vs_b, r.gene_id))
    return records


def find_group_markers(
    data: NormalizedData, cell_groups: List[str]
) -> List[GroupComparison]:
    """按分组的差异表达：每个分组 one-vs-rest，再按分组升序两两比较。

    ``cell_groups`` 与 ``data.cell_ids`` 对齐（仅含质控后保留细胞）。
    检验只用保留细胞的 log-归一化表达；返回先全部 one-vs-rest（分组升序）、
    再全部 (a, b)（a < b，均按分组升序）成对比较的比较列表；每个比较内
    记录按校正 P 值升序、log_fc_a_vs_b 降序、基因 ID 升序稳定排序，
    BH 校正以单个比较内的全部保留基因为一个校正家族。
    """
    group_to_cells: Dict[str, List[int]] = {}
    for c, group in enumerate(cell_groups):
        group_to_cells.setdefault(group, []).append(c)

    present = sorted(group_to_cells)
    comparisons: List[GroupComparison] = []
    for group in present:
        in_cells = group_to_cells[group]
        out_cells = [c for c in range(len(data.cell_ids)) if cell_groups[c] != group]
        records = _group_comparison(
            data, in_cells, out_cells, ONE_VS_REST, group, ""
        )
        comparisons.append((ONE_VS_REST, group, "", records))
    for i, group_a in enumerate(present):
        for group_b in present[i + 1:]:
            records = _group_comparison(
                data,
                group_to_cells[group_a],
                group_to_cells[group_b],
                PAIRWISE,
                group_a,
                group_b,
            )
            comparisons.append((PAIRWISE, group_a, group_b, records))
    return comparisons
