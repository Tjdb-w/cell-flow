"""簇内 pseudobulk 的配对生物学重复差异表达。

仅在同时提供 ``--replicate-metadata``、启用 ``--cluster-pseudobulk-de``
基线、并启用 ``--paired-cluster-pseudobulk-de`` 且提供 ``--pair-metadata``
时执行；未启用时全部行为与既有基线逐字节一致。

``--pair-metadata`` 指向 UTF-8 制表符文本或单成员 gzip（按 gzip 魔数识别，
与文件名无关），表头恰为 ``sample_id``、``pair_id`` 两列；每个样本恰好
一行，``sample_id`` 唯一且与重复元数据的样本集合一一对应（不多不少），
``pair_id`` 非空。任何格式不合法（表头不符、列数不符、空标识、样本重复、
覆盖不全、多出重复元数据之外的样本，或 gzip 多成员、尾随数据、截断、
校验失败）都抛 :class:`cell_flow.errors.CellFlowInputError`（退出码 2），
且不改动任何已有结果。

分析沿用簇内 pseudobulk 的保留范围、计数汇总与文库归一化（log1p）口径。
质控后恰有两个 group 时按字典序定 a、b；每个 pair 必须各含一个 a、b
样本，否则属 pair 归属错误，同样抛
:class:`cell_flow.errors.CellFlowInputError`。质控后 group 数不为两个时
无法定义 a、b 比较，抛 :class:`cell_flow.errors.CellFlowDataError`
（退出码 4）。

每个最终簇只用两个成员在该簇都有保留细胞的完整 pair；完整 pair 少于两个
的簇记为跳过，只要至少一个簇可检验即正常完成，没有任何簇可检验时由调用方
抛 :class:`cell_flow.errors.CellFlowDataError`（退出码 4）。每个基因对
a 减 b 的 pair 差值做双侧配对 t 检验（t 为差值均值除以差值样本标准差再乘
pair 数平方根，自由度为 pair 数减一），并在每个簇内跨保留基因做
Benjamini-Hochberg 校正。零方差且均值为零时 t_stat 为 0、p_value 为 1；
零方差且均值非零时 t_stat 为 inf 或 -inf、p_value 为 0。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .cluster_pseudobulk import _cluster_pseudobulk_values
from .errors import CellFlowDataError, CellFlowInputError
from .io import ExpressionMatrix, GZIP_MAGIC, gunzip_single
from .linalg import benjamini_hochberg, paired_ttest
from .qc import QCResult
from .replicate import ReplicateMetadata


@dataclass(frozen=True)
class PairMetadata:
    """样本 -> pair 的归属关系；每个重复元数据样本恰好一行。"""

    path: str
    name: str                             # 文件名（路径末段）
    sha256: str                           # 原始字节（gzip 即压缩字节）的 SHA-256
    pair_of_sample: Dict[str, str]        # sample_id -> pair_id
    samples_by_pair: Dict[str, List[str]]  # pair_id -> 该 pair 的样本（文件行序）

    @property
    def pair_count(self) -> int:
        return len(self.samples_by_pair)


@dataclass(frozen=True)
class PairedClusterPseudobulkRecord:
    """一个（簇, 基因）的配对簇内 pseudobulk 差异检验记录。"""

    cluster: int
    group_a: str
    group_b: str
    pair_count: int
    gene_id: str
    mean_difference: float   # a 减 b 的 pair 差值均值
    t_stat: float
    p_value: float
    p_value_adj: float


# (cluster, pair_count, records)
PairedClusterPseudobulkCluster = Tuple[int, int, List[PairedClusterPseudobulkRecord]]


@dataclass(frozen=True)
class PairedClusterPseudobulkDE:
    """全部最终簇的配对簇内 pseudobulk 差异结果。

    ``results`` 按 cluster 升序，每簇记录该簇完整 pair 数与逐基因记录；
    簇内记录按校正 P 值升序、mean_difference 降序、gene_id 升序排列。
    ``tested_clusters`` 与 ``skipped_clusters`` 互补且各自升序，合起来即
    全部最终簇。
    """

    group_a: str
    group_b: str
    gene_ids: List[str]                   # 质控保留基因（保留行序）
    tested_clusters: List[int]
    skipped_clusters: List[int]
    results: List[PairedClusterPseudobulkCluster]

    @property
    def test_count(self) -> int:
        """全部（簇, 基因）检验数。"""
        return sum(len(records) for _, _, records in self.results)

    @property
    def min_p_value_adj(self) -> float:
        return min(
            record.p_value_adj
            for _, _, records in self.results
            for record in records
        )


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_pair_metadata(path: str, sample_ids: List[str]) -> PairMetadata:
    """读取并校验配对元数据；``sample_ids`` 为重复元数据的全部样本 ID。

    样本集合必须与 ``sample_ids`` 一一对应（不多不少）。
    """
    if path is None or path == "":
        _fail("配对元数据路径为空")
    if not os.path.exists(path):
        _fail(f"配对元数据文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"配对元数据路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"配对元数据文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"配对元数据文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"配对元数据文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("配对元数据表头缺失：文件为空")
    header = lines[0].split("\t")
    if header != ["sample_id", "pair_id"]:
        _fail("配对元数据表头必须恰为 sample_id、pair_id 两列")

    pair_of_sample: Dict[str, str] = {}
    samples_by_pair: Dict[str, List[str]] = {}
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"配对元数据第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != 2:
            _fail(
                f"配对元数据第 {offset} 行列数为 {len(fields)}，"
                f"与表头两列不一致"
            )
        sample_id, pair_id = fields
        if sample_id == "":
            _fail(f"配对元数据第 {offset} 行样本 ID 为空")
        if pair_id == "":
            _fail(f"配对元数据第 {offset} 行 pair_id 为空（样本 {sample_id!r}）")
        if sample_id in pair_of_sample:
            _fail(f"配对元数据样本 ID 重复：{sample_id!r}")
        pair_of_sample[sample_id] = pair_id
        samples_by_pair.setdefault(pair_id, []).append(sample_id)

    missing = [sample_id for sample_id in sample_ids if sample_id not in pair_of_sample]
    if missing:
        _fail(
            f"配对元数据未覆盖重复元数据全部样本：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个"
        )
    known = set(sample_ids)
    extra = [sample_id for sample_id in pair_of_sample if sample_id not in known]
    if extra:
        _fail(
            f"配对元数据包含重复元数据之外的样本：{extra[0]!r} "
            f"等 {len(extra)} 个"
        )

    return PairMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        pair_of_sample=pair_of_sample,
        samples_by_pair=samples_by_pair,
    )


def _resolve_pairs(
    matrix: ExpressionMatrix,
    qc: QCResult,
    replicate: ReplicateMetadata,
    pair_metadata: PairMetadata,
) -> Tuple[str, str, List[Tuple[str, str, str]]]:
    """确定 a、b 两个 group 并校验每个 pair 的归属。

    返回（group_a, group_b, [（pair_id, a 样本, b 样本），…]），pair 按
    pair_id 字典序排列，保证差值求和顺序确定。质控后 group 数不为两个属
    数据错误；pair 归属不合法属输入错误。
    """
    kept_samples = {
        replicate.samples[matrix.cell_ids[c]] for c in qc.kept_cells
    }
    groups = sorted(
        {replicate.sample_group[sample_id] for sample_id in kept_samples}
    )
    if len(groups) != 2:
        raise CellFlowDataError(
            f"配对簇内 pseudobulk 差异表达要求质控后恰有两个 group，"
            f"实际为 {len(groups)} 个，无法确定 a、b 比较"
        )
    group_a, group_b = groups

    pairs: List[Tuple[str, str, str]] = []
    for pair_id in sorted(pair_metadata.samples_by_pair):
        members = pair_metadata.samples_by_pair[pair_id]
        member_groups = [replicate.sample_group[sample_id] for sample_id in members]
        if (
            len(members) != 2
            or member_groups.count(group_a) != 1
            or member_groups.count(group_b) != 1
        ):
            _fail(
                f"pair {pair_id!r} 归属错误：每个 pair 须各含一个 "
                f"group {group_a!r} 与 group {group_b!r} 的样本，"
                f"实际样本为 {members!r}"
            )
        sample_a = members[member_groups.index(group_a)]
        sample_b = members[member_groups.index(group_b)]
        pairs.append((pair_id, sample_a, sample_b))
    return group_a, group_b, pairs


def compute_paired_cluster_pseudobulk_de(
    matrix: ExpressionMatrix,
    qc: QCResult,
    labels: Sequence[int],
    replicate: ReplicateMetadata,
    pair_metadata: PairMetadata,
) -> PairedClusterPseudobulkDE:
    """按最终簇逐簇汇总样本 pseudobulk 并做簇内配对差异表达。

    ``labels`` 与 ``qc.kept_cells`` 对齐（每个保留细胞一个最终簇标签）。
    样本汇总、文库归一化与样本全局次序完全沿用簇内 pseudobulk 基线；
    每簇只用两个成员在该簇都有保留细胞的完整 pair，差值按 pair_id
    字典序收集以保证浮点求和顺序确定。
    """
    group_a, group_b, pairs = _resolve_pairs(
        matrix, qc, replicate, pair_metadata
    )

    # 全样本（质控后有保留细胞）首次出现顺序；与簇内 pseudobulk 基线一致
    sample_order: List[str] = []
    seen = set()
    for kept_index in qc.kept_cells:
        sample_id = replicate.samples[matrix.cell_ids[kept_index]]
        if sample_id not in seen:
            seen.add(sample_id)
            sample_order.append(sample_id)

    gene_ids = [matrix.gene_ids[g] for g in qc.kept_genes]

    tested_clusters: List[int] = []
    skipped_clusters: List[int] = []
    results: List[PairedClusterPseudobulkCluster] = []

    for cluster in sorted(set(labels)):
        active_samples, _, values = _cluster_pseudobulk_values(
            matrix, qc, labels, replicate, cluster, sample_order
        )
        position = {sample_id: i for i, sample_id in enumerate(active_samples)}
        # 完整 pair：两个成员在该簇都有保留细胞
        complete_pairs = [
            (sample_a, sample_b)
            for _, sample_a, sample_b in pairs
            if sample_a in position and sample_b in position
        ]
        pair_count = len(complete_pairs)
        if pair_count < 2:
            skipped_clusters.append(cluster)
            continue
        tested_clusters.append(cluster)

        records: List[PairedClusterPseudobulkRecord] = []
        pvalues: List[float] = []
        for g, gene_id in enumerate(gene_ids):
            row = values[g]
            differences = [
                row[position[sample_a]] - row[position[sample_b]]
                for sample_a, sample_b in complete_pairs
            ]
            mean_difference = sum(differences) / pair_count
            t_stat, p_value = paired_ttest(differences)
            pvalues.append(p_value)
            records.append(
                PairedClusterPseudobulkRecord(
                    cluster=cluster,
                    group_a=group_a,
                    group_b=group_b,
                    pair_count=pair_count,
                    gene_id=gene_id,
                    mean_difference=mean_difference,
                    t_stat=t_stat,
                    p_value=p_value,
                    p_value_adj=0.0,
                )
            )

        adjusted = benjamini_hochberg(pvalues)
        records = [
            PairedClusterPseudobulkRecord(
                cluster=r.cluster,
                group_a=r.group_a,
                group_b=r.group_b,
                pair_count=r.pair_count,
                gene_id=r.gene_id,
                mean_difference=r.mean_difference,
                t_stat=r.t_stat,
                p_value=r.p_value,
                p_value_adj=adjusted[g],
            )
            for g, r in enumerate(records)
        ]
        records.sort(
            key=lambda r: (r.p_value_adj, -r.mean_difference, r.gene_id)
        )
        results.append((cluster, pair_count, records))

    return PairedClusterPseudobulkDE(
        group_a=group_a,
        group_b=group_b,
        gene_ids=gene_ids,
        tested_clusters=tested_clusters,
        skipped_clusters=skipped_clusters,
        results=results,
    )
