"""Moderated t (limma-style empirical Bayes: Smyth fitFDist) and Benjamini-Hochberg."""
import math

import numpy as np
from scipy import special


def trigamma_inv(x):
    if x > 1e7:
        return 1 / math.sqrt(x)
    if x < 1e-6:
        return 1 / x
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y += dif
        if -dif / y < 1e-8:
            break
    return y


def fit_fdist(s2, df):
    """Prior degrees of freedom d0 and scale s0 of the variances s2 with residual degrees of freedom df."""
    ok = np.isfinite(s2) & (s2 > 1e-12) & (df > 0)
    x, d = s2[ok], df[ok]
    if len(x) < 3:
        return np.inf, float(np.nanmedian(s2))
    z = np.log(x)
    e = z - special.digamma(d / 2) + np.log(d / 2)
    em = e.mean()
    n = len(e)
    ev = ((e - em) ** 2).sum() / (n - 1) - np.mean(special.polygamma(1, d / 2))
    if ev > 0:
        d0 = 2 * trigamma_inv(ev)
        s0 = math.exp(em + special.digamma(d0 / 2) - math.log(d0 / 2))
    else:
        d0, s0 = np.inf, math.exp(em)
    return d0, s0


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    r = np.empty(n)
    r[o] = np.arange(1, n + 1)
    q = p * n / r
    qs = np.minimum.accumulate(q[o][::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(qs, 1)
    return out
