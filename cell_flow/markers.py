"""逐簇差异表达：每个簇对比其余所有细胞。

对每个（簇, 基因）输出：
- 两组在 log-归一化表达上的平均表达量；
- 对数倍数变化（簇内均值减其余均值）；
- Welch t 检验统计量与双侧 P 值；
- 每个簇内跨基因的 Benjamini-Hochberg 校正 P 值。
排序：校正 P 值升序、对数倍数变化降序、基因 ID 升序。
"""

from dataclasses import dataclass
from typing import Dict, List

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
            row = data.values[g]
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
