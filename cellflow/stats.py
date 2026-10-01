"""Deterministic numerical routines (pure standard library).

Everything here is free of wall-clock, platform-random and hash-order
dependencies so that repeated runs produce byte-identical output.
"""

import math
import random


# ---------------------------------------------------------------------------
# Special functions
# ---------------------------------------------------------------------------

def _betacf(a, b, x):
    """Continued fraction for the incomplete beta function."""
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


def betai(a, b, x):
    """Regularized incomplete beta function I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


# ---------------------------------------------------------------------------
# Welch's t-test and multiple-testing correction
# ---------------------------------------------------------------------------

def _mean_var(xs):
    n = len(xs)
    mean = math.fsum(xs) / n
    if n < 2:
        return mean, 0.0
    var = math.fsum((x - mean) ** 2 for x in xs) / (n - 1)
    return mean, var


def welch_t_pvalue(xs, ys):
    """Two-sided Welch t-test p-value for two samples."""
    n1, n2 = len(xs), len(ys)
    m1, v1 = _mean_var(xs)
    m2, v2 = _mean_var(ys)
    a1 = v1 / n1
    a2 = v2 / n2
    s = a1 + a2
    if s <= 0.0:
        return 1.0 if m1 == m2 else 0.0
    t = (m1 - m2) / math.sqrt(s)
    d1 = a1 * a1 / (n1 - 1) if n1 > 1 else 0.0
    d2 = a2 * a2 / (n2 - 1) if n2 > 1 else 0.0
    df_den = d1 + d2
    if df_den <= 0.0:
        return 1.0 if m1 == m2 else 0.0
    df = (s * s) / df_den
    # Two-sided p = I_{df/(df+t^2)}(df/2, 1/2)
    return betai(df / 2.0, 0.5, df / (df + t * t))


def bh_adjust(pvalues):
    """Benjamini-Hochberg adjusted p-values, input order preserved."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: (pvalues[i], i))
    adjusted = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        q = pvalues[i] * m / rank
        if q > 1.0:
            q = 1.0
        if q > prev:
            q = prev
        prev = q
        adjusted[i] = q
    return adjusted


# ---------------------------------------------------------------------------
# Eigendecomposition (cyclic Jacobi, symmetric matrices)
# ---------------------------------------------------------------------------

def jacobi_eigen(matrix):
    """Eigenpairs of a real symmetric matrix.

    Returns (eigenvalues, eigenvectors) where eigenvectors[i][k] is the
    i-th component of eigenvector k (i.e. eigenvectors are columns).
    """
    n = len(matrix)
    a = [row[:] for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    scale = math.fsum(abs(a[i][i]) for i in range(n)) + 1.0
    for _ in range(100):
        off = 0.0
        for p in range(n - 1):
            row = a[p]
            for q in range(p + 1, n):
                off += row[q] * row[q]
        if off <= 1e-24 * scale * scale:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if apq == 0.0:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * apq)
                sign = 1.0 if theta >= 0.0 else -1.0
                t = sign / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp = a[k][p]
                    akq = a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk = a[p][k]
                    aqk = a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp = v[k][p]
                    vkq = v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    return [a[i][i] for i in range(n)], v


def pca_scores(rows, n_components):
    """PCA scores for a cells x genes matrix of centered-able values.

    ``rows`` is a list of cells, each a list of gene values. The matrix is
    column-centered; scores come from the eigendecomposition of the
    cells x cells Gram matrix. Each component's sign is fixed so its
    largest-magnitude loading across cells is positive.
    """
    n_cells = len(rows)
    n_genes = len(rows[0])
    means = [math.fsum(rows[c][g] for c in range(n_cells)) / n_cells
             for g in range(n_genes)]
    centered = [[rows[c][g] - means[g] for g in range(n_genes)]
                for c in range(n_cells)]
    gram = [[0.0] * n_cells for _ in range(n_cells)]
    for i in range(n_cells):
        xi = centered[i]
        gi = gram[i]
        for j in range(i, n_cells):
            xj = centered[j]
            s = 0.0
            for k in range(n_genes):
                s += xi[k] * xj[k]
            gi[j] = s
            gram[j][i] = s
    eigenvalues, eigenvectors = jacobi_eigen(gram)
    order = sorted(range(n_cells), key=lambda k: (-eigenvalues[k], k))
    scores = [[0.0] * n_components for _ in range(n_cells)]
    for comp in range(n_components):
        k = order[comp]
        lam = eigenvalues[k]
        scale = math.sqrt(lam) if lam > 0.0 else 0.0
        col = [eigenvectors[c][k] * scale for c in range(n_cells)]
        pivot = max(range(n_cells), key=lambda c: (abs(col[c]), -c))
        if col[pivot] < 0.0:
            col = [-x for x in col]
        for c in range(n_cells):
            scores[c][comp] = col[c]
    return scores


# ---------------------------------------------------------------------------
# k-means
# ---------------------------------------------------------------------------

def _sqdist(a, b):
    s = 0.0
    for x, y in zip(a, b):
        d = x - y
        s += d * d
    return s


def kmeans(points, k, seed, max_iter=100):
    """Lloyd's k-means with seeded k-means++ initialisation.

    Fully deterministic for a given (points, k, seed). Every returned
    cluster is non-empty. Requires k <= len(points).
    """
    n = len(points)
    rng = random.Random(seed)

    used = set()
    first = rng.randrange(n)
    used.add(first)
    centers = [list(points[first])]
    while len(centers) < k:
        d2 = [min(_sqdist(p, c) for c in centers) for p in points]
        total = math.fsum(d2)
        if total <= 0.0:
            idx = next(i for i in range(n) if i not in used)
        else:
            r = rng.random() * total
            acc = 0.0
            idx = n - 1
            for i in range(n):
                acc += d2[i]
                if r < acc:
                    idx = i
                    break
        used.add(idx)
        centers.append(list(points[idx]))

    assign = [-1] * n
    for _ in range(max_iter):
        changed = False
        for i in range(n):
            best = 0
            best_d = _sqdist(points[i], centers[0])
            for c in range(1, k):
                d = _sqdist(points[i], centers[c])
                if d < best_d:
                    best_d = d
                    best = c
            if assign[i] != best:
                assign[i] = best
                changed = True
        if not changed:
            break
        dims = len(points[0])
        sums = [[0.0] * dims for _ in range(k)]
        counts = [0] * k
        for i in range(n):
            c = assign[i]
            counts[c] += 1
            pi = points[i]
            sc = sums[c]
            for d in range(dims):
                sc[d] += pi[d]
        for c in range(k):
            if counts[c] == 0:
                # Reseed the empty cluster at the point farthest from its
                # current center (lowest index wins ties).
                best_i = 0
                best_d = -1.0
                for i in range(n):
                    d = _sqdist(points[i], centers[assign[i]])
                    if d > best_d:
                        best_d = d
                        best_i = i
                centers[c] = list(points[best_i])
            else:
                centers[c] = [s / counts[c] for s in sums[c]]

    # Guarantee non-empty clusters with a deterministic fix-up pass.
    while True:
        counts = [0] * k
        for a in assign:
            counts[a] += 1
        empty = next((c for c in range(k) if counts[c] == 0), None)
        if empty is None:
            break
        best_i = -1
        best_d = -1.0
        for i in range(n):
            if counts[assign[i]] > 1:
                d = _sqdist(points[i], centers[assign[i]])
                if d > best_d:
                    best_d = d
                    best_i = i
        assign[best_i] = empty
        centers[empty] = list(points[best_i])

    return assign
