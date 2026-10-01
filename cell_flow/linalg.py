"""纯 Python 统计与线性代数基础，全部使用 float64，结果确定。

不依赖 numpy：算法简单、规模可控（细胞 x 高变基因矩阵），且能保证
跨环境、跨版本逐位一致。
"""

import math
from typing import List, Sequence


# --------------------------------------------------------------------------- #
# 基本统计
# --------------------------------------------------------------------------- #


def mean(values: Sequence[float]) -> float:
    return sum(values) / len(values)


def unbiased_variance(values: Sequence[float], mu: float) -> float:
    """样本（无偏）方差；单点或无点定义为 0。"""
    n = len(values)
    if n < 2:
        return 0.0
    total = 0.0
    for x in values:
        d = x - mu
        total += d * d
    return total / (n - 1)


def _log_beta(a: float, b: float) -> float:
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def _betacf(a: float, b: float, x: float) -> float:
    """连分式展开（Numerical Recipes betacf），用于正则化不完全 Beta 函数。"""
    max_iter = 200
    eps = 3.0e-14
    fpmin = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _betainc_regularized(x: float, a: float, b: float) -> float:
    """正则化不完全 Beta 函数 I_x(a, b)，x in [0,1]。"""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(-_log_beta(a, b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_sf_two_sided(t: float, df: float) -> float:
    """Student-t 分布的双侧生存函数 P(|T| >= |t|)。

    利用 I_{df/(df+t^2)}(df/2, 1/2) = P(|T| <= |t|) 的关系。
    df 非有限时退化为标准正态近似。
    """
    t = abs(t)
    if not math.isfinite(df) or df <= 0.0:
        # 标准正态双侧
        return 2.0 * (1.0 - 0.5 * (1.0 + math.erf(t / math.sqrt(2.0))))
    if t == 0.0:
        return 1.0
    x = df / (df + t * t)
    # 双侧生存函数 = I_{df/(df+t^2)}(df/2, 1/2)
    p = _betainc_regularized(x, df / 2.0, 0.5)
    if p < 0.0:
        p = 0.0
    if p > 1.0:
        p = 1.0
    return p


def welch_ttest(a: Sequence[float], b: Sequence[float]) -> tuple[float, float]:
    """Welch 不等方差 t 检验，返回 (t 统计量, 双侧 P 值)。

    退化情形（两组方差均为 0）：均值相同 => t=0,p=1；均值不同 => t=inf,p=0。
    """
    n1, n2 = len(a), len(b)
    m1, m2 = mean(a), mean(b)
    v1 = unbiased_variance(a, m1)
    v2 = unbiased_variance(b, m2)
    diff = m1 - m2
    se2 = v1 / n1 + v2 / n2
    if se2 <= 0.0:
        if diff == 0.0:
            return 0.0, 1.0
        return math.inf if diff > 0 else -math.inf, 0.0
    t = diff / math.sqrt(se2)
    # 单点组的无偏方差为 0，对自由度的贡献项按 0 处理（避免 0/0）
    term1 = (v1 / n1) ** 2 / (n1 - 1) if n1 > 1 else 0.0
    term2 = (v2 / n2) ** 2 / (n2 - 1) if n2 > 1 else 0.0
    den = term1 + term2
    df = (se2 * se2) / den if den > 0.0 else math.inf
    return t, t_sf_two_sided(t, df)


def benjamini_hochberg(pvalues: Sequence[float]) -> List[float]:
    """BH-FDR 校正，返回与输入等长的 q 值；并列同值、单调化处理确定。"""
    n = len(pvalues)
    if n == 0:
        return []
    # (p 值, 原下标) 排序；p 相同按原下标，保证完全确定
    order = sorted(range(n), key=lambda i: (pvalues[i], i))
    adjusted = [0.0] * n
    running = 1.0
    for rank_from_end in range(n, 0, -1):
        idx = order[rank_from_end - 1]
        value = pvalues[idx] * n / rank_from_end
        if value < running:
            running = value
        adjusted[idx] = running
    return [min(1.0, max(0.0, q)) for q in adjusted]


# --------------------------------------------------------------------------- #
# 对称矩阵特征分解（经典 Jacobi 旋转）
# --------------------------------------------------------------------------- #


def jacobi_eigh(
    matrix: Sequence[Sequence[float]],
) -> tuple[List[float], List[List[float]]]:
    """实对称矩阵的特征分解。

    返回 ``(特征值升序, 特征向量列)``：``eigenvectors[k]`` 是对应
    ``eigenvalues[k]`` 的单位特征向量（按行返回，索引与特征值对齐）。
    纯 float64 迭代，结果确定。
    """
    n = len(matrix)
    # 复制到可变二维数组
    a = [list(row) for row in matrix]
    # v 为累积旋转（特征向量按列存放：v[i][k]）
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    max_sweeps = 100
    for _ in range(max_sweeps):
        # 收敛判定：所有非对角元足够小
        off = 0.0
        for p in range(n):
            for q in range(p + 1, n):
                off += a[p][q] * a[p][q]
        if off < 1.0e-30:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if apq == 0.0:
                    continue
                app = a[p][p]
                aqq = a[q][q]
                tau = (aqq - app) / (2.0 * apq)
                if tau >= 0.0:
                    t = 1.0 / (tau + math.sqrt(1.0 + tau * tau))
                else:
                    t = -1.0 / (-tau + math.sqrt(1.0 + tau * tau))
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c

                # 旋转对称矩阵（仅更新需要的部分，利用对称性）
                for k in range(n):
                    if k == p or k == q:
                        continue
                    akp = a[k][p]
                    akq = a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[p][k] = a[k][p]
                    a[k][q] = s * akp + c * akq
                    a[q][k] = a[k][q]
                a[p][p] = app - t * apq
                a[q][q] = aqq + t * apq
                a[p][q] = 0.0
                a[q][p] = 0.0
                for k in range(n):
                    vkp = v[k][p]
                    vkq = v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq

    eigenvalues = [a[i][i] for i in range(n)]
    eigenvectors = [[v[i][k] for i in range(n)] for k in range(n)]
    # 特征值升序排序，并列按原索引，保证结果确定
    order = sorted(range(n), key=lambda k: (eigenvalues[k], k))
    sorted_values = [eigenvalues[k] for k in order]
    sorted_vectors = [eigenvectors[k] for k in order]
    return sorted_values, sorted_vectors
