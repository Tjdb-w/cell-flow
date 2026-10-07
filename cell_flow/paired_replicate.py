"""配对生物学重复元数据读取。

``--paired-replicate-metadata`` 指向 UTF-8 制表符文本，表头恰为
``sample_id``、``pair_id`` 两列；每个 ``--replicate-metadata`` 样本恰好一行，
``sample_id`` 唯一、``pair_id`` 非空，不多不少。文件可为纯文本或单成员
gzip（按 gzip 魔数识别，与文件名无关）。任何不合法（路径不可读、表头不符、
列数不符、空样本或空配对标识、``sample_id`` 重复、未覆盖全部重复样本、
出现重复元数据之外的样本，或 gzip 多成员、尾随数据、截断、校验失败）都抛
:class:`cell_flow.errors.CellFlowInputError`（退出码 2）。

pair 与质控后两个 group 的归属约束（每个 pair 恰含按字典序确定的 a、b
各一个样本）在实际分析阶段校验，同样报输入错误。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single


@dataclass(frozen=True)
class PairedReplicateMetadata:
    """样本 -> 配对标识 的归属关系；每个重复样本恰好一行。"""

    path: str
    name: str                      # 文件名（路径末段）
    sha256: str                    # 原始字节（gzip 即压缩字节）的 SHA-256
    sample_pair: Dict[str, str]    # sample_id -> pair_id
    pair_order: List[str]          # pair_id 按文件中首次出现顺序


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_paired_replicate_metadata(
    path: str, sample_ids: List[str]
) -> PairedReplicateMetadata:
    """读取并校验配对元数据。

    ``sample_ids`` 为 ``--replicate-metadata`` 的全部样本（其首次出现顺序）；
    配对表必须与之一一对应：每个样本恰好一行，且不得出现额外样本。
    """
    if path is None or path == "":
        _fail("配对重复元数据路径为空")
    if not os.path.exists(path):
        _fail(f"配对重复元数据文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"配对重复元数据路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"配对重复元数据文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"配对重复元数据文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"配对重复元数据文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("配对重复元数据表头缺失：文件为空")
    header = lines[0].split("\t")
    if header != ["sample_id", "pair_id"]:
        _fail("配对重复元数据表头必须恰为 sample_id、pair_id 两列")

    sample_pair: Dict[str, str] = {}
    pair_order: List[str] = []
    seen_pairs = set()
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"配对重复元数据第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != 2:
            _fail(
                f"配对重复元数据第 {offset} 行列数为 {len(fields)}，"
                f"与表头两列不一致"
            )
        sample_id, pair_id = fields
        if sample_id == "":
            _fail(f"配对重复元数据第 {offset} 行样本 ID 为空")
        if pair_id == "":
            _fail(f"配对重复元数据第 {offset} 行配对标识为空（样本 {sample_id!r}）")
        if sample_id in sample_pair:
            _fail(f"配对重复元数据样本 ID 重复：{sample_id!r}")
        sample_pair[sample_id] = pair_id
        if pair_id not in seen_pairs:
            seen_pairs.add(pair_id)
            pair_order.append(pair_id)

    missing = [sample_id for sample_id in sample_ids if sample_id not in sample_pair]
    if missing:
        _fail(
            f"配对重复元数据未覆盖全部重复样本：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个"
        )
    known = set(sample_ids)
    extra = [sample_id for sample_id in sample_pair if sample_id not in known]
    if extra:
        _fail(
            f"配对重复元数据包含重复元数据之外的样本：{extra[0]!r} "
            f"等 {len(extra)} 个"
        )

    return PairedReplicateMetadata(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        sample_pair=sample_pair,
        pair_order=pair_order,
    )
