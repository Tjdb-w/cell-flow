"""成对簇间差异表达：任意两个实际出现的簇两两比较。

对规范化标签从 0 开始且实际出现的每一对簇 a、b（a < b），case 为 a、
control 为 b，对每个保留基因输出：
- 两组在 log-归一化表达上的平均表达量；
- log_fc_a_vs_b（a 组均值减 b 组均值）；
- Welch t 检验统计量与双侧 P 值；
- 每个比较内跨全部保留基因的 Benjamini-Hochberg 校正 P 值。
记录先按 cluster_a、cluster_b 升序，比较内再按校正 P 值升序、
log_fc_a_vs_b 降序、基因 ID 升序稳定排序。
"""

from dataclasses import dataclass
from typing import Dict, List

from .linalg import benjamini_hochberg, welch_ttest
from .normalize import NormalizedData


@dataclass(frozen=True)
class PairwiseMarkerRecord:
    cluster_a: int
    cluster_b: int
    gene_id: str
    mean_in_a: float
    mean_in_b: float
    log_fc: float
    t_stat: float
    p_value: float
    p_value_adj: float


def find_pairwise_markers(
    data: NormalizedData, labels: List[int]
) -> List[PairwiseMarkerRecord]:
    cluster_to_cells: Dict[int, List[int]] = {}
    for c, label in enumerate(labels):
        cluster_to_cells.setdefault(label, []).append(c)

    present = sorted(cluster_to_cells)
    results: List[PairwiseMarkerRecord] = []
    for i, cluster_a in enumerate(present):
        for cluster_b in present[i + 1:]:
            cells_a = cluster_to_cells[cluster_a]
            cells_b = cluster_to_cells[cluster_b]

            records: List[PairwiseMarkerRecord] = []
            pvalues: List[float] = []
            for g, gene_id in enumerate(data.gene_ids):
                row = data.values[g]
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
                        log_fc=log_fc,
                        t_stat=t_stat,
                        p_value=p_value,
                        p_value_adj=0.0,
                    )
                )

            adjusted = benjamini_hochberg(pvalues)
            for g, r in enumerate(records):
                results.append(
                    PairwiseMarkerRecord(
                        cluster_a=r.cluster_a,
                        cluster_b=r.cluster_b,
                        gene_id=r.gene_id,
                        mean_in_a=r.mean_in_a,
                        mean_in_b=r.mean_in_b,
                        log_fc=r.log_fc,
                        t_stat=r.t_stat,
                        p_value=r.p_value,
                        p_value_adj=adjusted[g],
                    )
                )
    # 先按 cluster_a、cluster_b 升序，再按 p_value_adj 升序、
    # log_fc_a_vs_b 降序、gene_id 升序稳定排序
    results.sort(
        key=lambda r: (
            r.cluster_a,
            r.cluster_b,
            r.p_value_adj,
            -r.log_fc,
            r.gene_id,
        )
    )
    return results
