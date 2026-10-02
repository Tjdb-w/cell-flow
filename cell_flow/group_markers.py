"""按元数据分组的差异表达：每组 one-vs-rest，再两两成对比较。

仅使用质控后保留细胞与保留基因的 log-归一化表达。对每个分组、每个基因：
- 两组在 log-归一化表达上的平均表达量（one-vs-rest 时第二组为其余全部细胞）；
- 对数倍数变化（均值差，case 均值减 control 均值）；
- Welch t 检验统计量与双侧 P 值；
- 单个比较内跨全部保留基因的 Benjamini-Hochberg 校正 P 值。

比较顺序：先按 group 升序逐组做 one-vs-rest，再按 group 升序对每个
``a < b`` 组合做成对检验。每个比较内记录按校正 P 值升序、log_fc 降序、
基因 ID 升序排序。
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple

from .linalg import benjamini_hochberg, welch_ttest
from .normalize import NormalizedData

ONE_VS_REST = "one_vs_rest"
PAIRWISE = "pairwise"


@dataclass(frozen=True)
class GroupMarkerRecord:
    comparison_type: str
    group_a: str
    group_b: str  # one-vs-rest 时为空字符串
    gene_id: str
    mean_in_a: float
    mean_in_b: float
    log_fc: float
    t_stat: float
    p_value: float
    p_value_adj: float


# （比较类型, group_a, group_b, 该比较的逐基因记录）
GroupMarkerComparison = Tuple[str, str, str, List[GroupMarkerRecord]]


def _compare(
    comparison_type: str,
    group_a: str,
    group_b: str,
    cells_a: List[int],
    cells_b: List[int],
    data: NormalizedData,
) -> GroupMarkerComparison:
    records: List[GroupMarkerRecord] = []
    pvalues: List[float] = []
    for g, gene_id in enumerate(data.gene_ids):
        row = data.values[g]
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
                log_fc=log_fc,
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
            log_fc=r.log_fc,
            t_stat=r.t_stat,
            p_value=r.p_value,
            p_value_adj=adjusted[g],
        )
        for g, r in enumerate(records)
    ]
    records.sort(key=lambda r: (r.p_value_adj, -r.log_fc, r.gene_id))
    return (comparison_type, group_a, group_b, records)


def find_group_markers(
    data: NormalizedData, cell_groups: Dict[str, str]
) -> List[GroupMarkerComparison]:
    """按 ``cell_id -> group`` 对质控后细胞做分组差异表达。

    返回的比较列表顺序为：各组（group 升序）one-vs-rest 在前，两两成对
    比较（group_a < group_b，双关键字升序）在后。``cell_groups`` 中未覆盖
    质控后细胞的情形由输入校验阶段排除，此处不再发生。
    """
    group_to_cells: Dict[str, List[int]] = {}
    cell_group: Dict[int, str] = {}
    for c, cell_id in enumerate(data.cell_ids):
        group = cell_groups[cell_id]
        group_to_cells.setdefault(group, []).append(c)
        cell_group[c] = group

    present = sorted(group_to_cells)

    comparisons: List[GroupMarkerComparison] = []

    # 1) 每组 one-vs-rest（group 升序）
    for group in present:
        cells_in = group_to_cells[group]
        cells_out = [c for c, g in cell_group.items() if g != group]
        comparisons.append(
            _compare(
                ONE_VS_REST,
                group,
                "",
                cells_in,
                cells_out,
                data,
            )
        )

    # 2) 两两成对（group 升序，a < b）
    for i, group_a in enumerate(present):
        for group_b in present[i + 1:]:
            comparisons.append(
                _compare(
                    PAIRWISE,
                    group_a,
                    group_b,
                    group_to_cells[group_a],
                    group_to_cells[group_b],
                    data,
                )
            )

    return comparisons
