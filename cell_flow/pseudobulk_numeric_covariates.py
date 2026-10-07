"""数值型样本协变量元数据读取。

``--pseudobulk-numeric-covariates`` 指向 UTF-8 制表符文本，首列恰为
``sample_id``，另有至少一个唯一命名的协变量列；每个
``--replicate-metadata`` 样本恰好一行，``sample_id`` 唯一，全部字段非空，
样本集合与重复元数据完全一致（不多不少）。与
``--pseudobulk-covariates`` 的离散标签不同，本表每个协变量取值都必须
是可解析的**有限**浮点数（``nan``、``inf``、``-inf`` 与不可解析文本一律
拒绝），在回归中作为连续列使用。文件可为纯文本或单成员 gzip（按 gzip
魔数识别，与文件名无关）。任何不合法（路径不可读、表头不符、协变量列
缺失或命名重复、列数不符、空字段、非法或非有限数值、``sample_id``
重复、未覆盖全部重复样本、出现重复元数据之外的样本，或 gzip 多成员、
尾随数据、截断、校验失败）都抛
:class:`cell_flow.errors.CellFlowInputError`（退出码 2）。
"""

import hashlib
import math
import os
from dataclasses import dataclass
from typing import Dict, List

from .errors import CellFlowInputError
from .io import GZIP_MAGIC, gunzip_single


@dataclass(frozen=True)
class PseudobulkNumericCovariates:
    """样本 -> 数值协变量取值；每个重复样本恰好一行，取值均为有限浮点。"""

    path: str
    name: str                                # 文件名（路径末段）
    sha256: str                              # 原始字节（gzip 即压缩字节）的 SHA-256
    covariate_names: List[str]               # 数值协变量列名（表头顺序）
    sample_values: Dict[str, Dict[str, float]]  # sample_id -> {协变量: 有限浮点取值}


def _fail(message: str) -> None:
    raise CellFlowInputError(message)


def _parse_finite_float(value: str) -> float:
    """解析有限浮点取值；不可解析或非有限（nan/inf）返回 None。"""
    try:
        number = float(value)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def read_pseudobulk_numeric_covariates(
    path: str, sample_ids: List[str]
) -> PseudobulkNumericCovariates:
    """读取并校验数值型样本协变量表。

    ``sample_ids`` 为 ``--replicate-metadata`` 的全部样本（其首次出现顺序）；
    协变量表必须与之一一对应：每个样本恰好一行，且不得出现额外样本。
    """
    if path is None or path == "":
        _fail("数值协变量路径为空")
    if not os.path.exists(path):
        _fail(f"数值协变量文件不存在：{path}")
    if not os.path.isfile(path):
        _fail(f"数值协变量路径不是普通文件：{path}")
    if not os.access(path, os.R_OK):
        _fail(f"数值协变量文件不可读：{path}")

    # SHA-256 始终按磁盘上的实际原始字节计算（gzip 即压缩字节）
    try:
        with open(path, "rb") as handle:
            raw_bytes = handle.read()
    except OSError as exc:
        _fail(f"数值协变量文件不可读：{path}（{exc}）")

    digest = hashlib.sha256(raw_bytes).hexdigest()
    if raw_bytes[:2] == GZIP_MAGIC:
        text_bytes = gunzip_single(raw_bytes, path)
    else:
        text_bytes = raw_bytes
    try:
        text = text_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(f"数值协变量文件不是合法的 UTF-8 文本：{path}")

    # splitlines 同时兼容 \n 与 \r\n，且不会因末尾换行产生空行
    lines = text.splitlines()
    if not lines:
        _fail("数值协变量表头缺失：文件为空")
    header = lines[0].split("\t")
    if header[0] != "sample_id":
        _fail("数值协变量表头首列必须为 sample_id")
    covariate_names = header[1:]
    if not covariate_names:
        _fail("数值协变量表头缺少协变量列：至少需要一列协变量")
    if any(name == "" for name in covariate_names):
        _fail("数值协变量表头存在空协变量列名")
    if len(set(covariate_names)) != len(covariate_names):
        seen = set()
        duplicate = next(
            name for name in covariate_names if name in seen or seen.add(name)
        )
        _fail(f"数值协变量列名重复：{duplicate!r}")

    n_columns = len(header)
    sample_values: Dict[str, Dict[str, float]] = {}
    for offset, line in enumerate(lines[1:], start=2):
        if line == "":
            _fail(f"数值协变量第 {offset} 行为空，无法解析")
        fields = line.split("\t")
        if len(fields) != n_columns:
            _fail(
                f"数值协变量第 {offset} 行列数为 {len(fields)}，"
                f"与表头 {n_columns} 列不一致"
            )
        sample_id = fields[0]
        if sample_id == "":
            _fail(f"数值协变量第 {offset} 行样本 ID 为空")
        parsed: Dict[str, float] = {}
        for name, value in zip(covariate_names, fields[1:]):
            if value == "":
                _fail(
                    f"数值协变量第 {offset} 行协变量 {name!r} 为空"
                    f"（样本 {sample_id!r}）"
                )
            number = _parse_finite_float(value)
            if number is None:
                _fail(
                    f"数值协变量第 {offset} 行协变量 {name!r} 不是有限数值："
                    f"{value!r}（样本 {sample_id!r}）"
                )
            parsed[name] = number
        if sample_id in sample_values:
            _fail(f"数值协变量样本 ID 重复：{sample_id!r}")
        sample_values[sample_id] = parsed

    missing = [s for s in sample_ids if s not in sample_values]
    if missing:
        _fail(
            f"数值协变量未覆盖全部重复样本：缺少 {missing[0]!r} "
            f"等 {len(missing)} 个"
        )
    known = set(sample_ids)
    extra = [s for s in sample_values if s not in known]
    if extra:
        _fail(
            f"数值协变量包含重复元数据之外的样本：{extra[0]!r} "
            f"等 {len(extra)} 个"
        )

    return PseudobulkNumericCovariates(
        path=path,
        name=os.path.basename(path),
        sha256=digest,
        covariate_names=covariate_names,
        sample_values=sample_values,
    )
