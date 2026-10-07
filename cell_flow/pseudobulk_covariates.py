"""pseudobulk 样本协变量元数据读取。

``--pseudobulk-covariates`` 指向 UTF-8 制表符文本，首列恰为 ``sample_id``，
另有至少一个唯一命名的协变量列；每个 ``--replicate-metadata`` 样本恰好
一行，``sample_id`` 唯一、所有字段非空，样本集合与重复元数据完全一致
（不多不少）。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名
无关）。协变量一律按离散标签处理，取值不做数值解释。任何不合法（路径
不可读、表头不符、协变量列缺失或命名重复、列数不符、空字段、样本重复、
未覆盖全部重复样本、出现重复元数据之外的样本，或 gzip 多成员、尾随
数据、截断、校验失败）都抛 :class:`cell_flow.errors.CellFlowInputError`
（退出码 2），且不改动任何已有结果。
"""

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single


@dataclass(frozen=True)
class PseudobulkCovariates:
    """样本 -> 协变量标签 的归属关系；每个重复样本恰好一行。"""

    path: str
    name: str                          # 文件名（路径末段）
    sha256: str                        # 原始字节（gzip 即压缩字节）的 SHA-256
    covariate_columns: List[str]       # 协变量列名（文件列序，至少一个）
    values: Dict[str, Dict[str, str]]  # sample_id -> {协变量列名: 标签}


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def read_pseudobulk_covariates(
    path: str, sample_ids: List[str]
) -> PseudobulkCovariates:
    """读取并校验样本协变量表。

    ``sample_ids`` 为 ``--replicate-metadata`` 的全部样本（其首次出现顺序）；
    协变量表必须与之一一对应：每个样本恰好一行，且不得出现额外样本。
    """
    if path is None or path == "":
        _fail("样本协变量路径为空")
    if not os.path.exists(path):
        _fail(f"样本协变量文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"样本协变量路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"样本协变量文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"样本协变量文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"样本协变量文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("样本协变量表头缺失：文件为空")
    header = lines[0].split("\t")
    if len(header) < 2 or header[0] != "sample_id":
        _fail(
            "样本协变量表头首列必须为 sample_id，且至少另有一个协变量列"
        )
    covariate_columns = header[1:]
    seen_names = set()
    for name in covariate_columns:
        if name == "":
            _fail("样本协变量表头存在空的协变量列名")
        if name in seen_names or name == "sample_id":
            _fail(f"样本协变量表头列名重复：{name!r}")
        seen_names.add(name)

    n_columns = len(header)
    values: Dict[str, Dict[str, str]] = {}
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"样本协变量第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != n_columns:
            _fail(
                f"样本协变量第 {offset} 行列数为 {len(fields)}，"
                f"与表头 {n_columns} 列不一致"
            )
        sample_id = fields[0]
        if sample_id == "":
            _fail(f"样本协变量第 {offset} 行样本 ID 为空")
        for name, field in zip(covariate_columns, fields[1:]):
            if field == "":
                _fail(
                    f"样本协变量第 {offset} 行协变量 {name!r} 为空"
                    f"（样本 {sample_id!r}）"
                )
        if sample_id in values:
            _fail(f"样本协变量样本 ID 重复：{sample_id!r}")
        values[sample_id] = dict(zip(covariate_columns, fields[1:]))

    missing = [sample_id for sample_id in sample_ids if sample_id not in values]
    if missing:
        _fail(
            f"样本协变量未覆盖全部重复样本：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个"
        )
    known = set(sample_ids)
    extra = [sample_id for sample_id in values if sample_id not in known]
    if extra:
        _fail(
            f"样本协变量包含重复元数据之外的样本：{extra[0]!r} "
            f"等 {len(extra)} 个"
        )

    return PseudobulkCovariates(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        covariate_columns=covariate_columns,
        values=values,
    )
