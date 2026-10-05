"""pseudobulk 基因集评分与按样本（生物学重复）分组的基因集差异。

仅在同时提供 ``--gene-sets`` 与 ``--replicate-metadata`` 并启用
``--pseudobulk-gene-set-de`` 时执行。评分只使用质控后保留的基因与细胞：
每个样本取其保留细胞在保留基因上汇总出的 pseudobulk log1p 值
（见 :mod:`cell_flow.replicate`），``set_id`` 只取与保留基因 ID 的交集，
score 为交集基因在该样本 pseudobulk 值上的算术平均；矩阵外或被 QC 剔除的
成员计入 ``n_genes_total`` 但不计分。任一集合交集为空都抛
:class:`cell_flow.errors.CellFlowDataError`（退出码 4）。

差异表达以样本为观测单位，完全沿用 pseudobulk 基因差异的比较口径：
每个 group 先 one-vs-rest，再按 group 升序两两比较，执行 Welch t 双侧检验、
均值差（score_difference_a_vs_b），并在每个比较内跨集合做 BH 校正。
"""

from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowDataError
from .gene_sets import GeneSets
from .markers import GroupComparison, run_group_comparisons
from .replicate import PseudobulkData


@dataclass(frozen=True)
class PseudobulkGeneSetScore:
    """一个基因集对一个有效样本的 pseudobulk 评分结果。"""

    set_id: str
    n_total: int
    n_used: int
    sample_id: str
    score: float


@dataclass(frozen=True)
class PseudobulkGeneSetScores:
    """全部基因集在全部有效样本上的评分与逐比较差异结果。

    ``rows`` 按 set_id 升序、集合内按样本顺序（pseudobulk 样本次序，
    即输入细胞在矩阵列序中的首次出现顺序）排列；``set_order`` 为升序
    set_id 列表，供图表列序使用。``comparisons`` 先全部 one-vs-rest
    （group 升序）再全部 pairwise（group 升序）。
    """

    rows: List[PseudobulkGeneSetScore]
    set_order: List[str]
    totals: Dict[str, int]    # set_id -> n_genes_total
    used: Dict[str, int]      # set_id -> n_genes_used
    comparisons: List[GroupComparison]


def score_pseudobulk_gene_sets(
    gene_sets: GeneSets,
    bulk: PseudobulkData,
) -> PseudobulkGeneSetScores:
    """计算各基因集对每个 pseudobulk 样本的算术平均评分并做分组差异。

    评分矩阵的行是与保留基因有非空交集的集合、列是 ``bulk.sample_ids``
    对齐的有效样本，取值为该样本 pseudobulk log1p 值在交集基因上的平均。
    任一集合与保留基因交集为空即无法评分，抛数据错误（退出码 4），
    不产出任何结果。
    """
    gene_position: Dict[str, int] = {
        gene_id: g for g, gene_id in enumerate(bulk.gene_ids)
    }

    totals: Dict[str, int] = {}
    used: Dict[str, int] = {}
    set_positions: Dict[str, List[int]] = {}
    for set_id in sorted(gene_sets.members):
        entries = gene_sets.members[set_id]
        # 集合内成员在文件中已保证 (set_id, gene_id) 唯一，直接建下标映射
        positions = [
            gene_position[gene_id]
            for gene_id in entries
            if gene_id in gene_position
        ]
        totals[set_id] = len(entries)
        used[set_id] = len(positions)
        if not positions:
            raise CellFlowDataError(
                f"基因集 {set_id!r} 与质控后保留基因的交集为空，"
                f"无法进行 pseudobulk 基因集评分"
            )
        set_positions[set_id] = positions

    set_order = sorted(gene_sets.members)
    rows: List[PseudobulkGeneSetScore] = []
    # 差异检验的观测矩阵：values[s][sample]，行序与 set_order 对齐
    values: List[List[float]] = []
    for set_id in set_order:
        positions = set_positions[set_id]
        n_used = len(positions)
        n_total = totals[set_id]
        score_row: List[float] = []
        for s, sample_id in enumerate(bulk.sample_ids):
            # 沿用既有基因集评分口径（内置 sum）累加，保证浮点口径一致
            total = sum(bulk.values[g][s] for g in positions)
            score = total / n_used
            score_row.append(score)
            rows.append(
                PseudobulkGeneSetScore(
                    set_id=set_id,
                    n_total=n_total,
                    n_used=n_used,
                    sample_id=sample_id,
                    score=score,
                )
            )
        values.append(score_row)

    # 以样本为观测复用分组比较口径（one-vs-rest + 升序两两、Welch t、
    # 比较内 BH）；记录中的 gene_id 即 set_id、log_fc_a_vs_b 即评分均值差
    comparisons = run_group_comparisons(set_order, values, bulk.sample_groups)

    return PseudobulkGeneSetScores(
        rows=rows,
        set_order=set_order,
        totals=totals,
        used=used,
        comparisons=comparisons,
    )
