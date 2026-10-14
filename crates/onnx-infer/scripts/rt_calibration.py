"""Per-run retention-time calibration: predicted index -> observed minutes.

Isotonic regression (pool-adjacent-violators) fitted on confident target PSMs,
then applied unchanged to every row (targets and decoys). Predictions between
fitted points are interpolated linearly and clamped at the ends. With fewer
than ``min_points`` anchors a least-squares line is used instead.

Plain numpy; no scikit-learn dependency.
"""
from __future__ import annotations

import numpy as np


def pav(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted non-decreasing least-squares fit of ``y`` (already ordered by x)."""
    n = len(y)
    val = np.empty(n)
    wt = np.empty(n)
    cnt = np.empty(n, dtype=np.int64)
    k = -1
    for i in range(n):
        k += 1
        val[k], wt[k], cnt[k] = y[i], w[i], 1
        while k > 0 and val[k - 1] > val[k]:
            tw = wt[k - 1] + wt[k]
            val[k - 1] = (val[k - 1] * wt[k - 1] + val[k] * wt[k]) / tw
            wt[k - 1] = tw
            cnt[k - 1] += cnt[k]
            k -= 1
    return np.repeat(val[: k + 1], cnt[: k + 1])


def fit_isotonic(pred: np.ndarray, obs: np.ndarray, min_points: int = 20):
    """Return ``f(index) -> minutes`` fitted on anchor pairs (pred, obs)."""
    pred = np.asarray(pred, float)
    obs = np.asarray(obs, float)
    ok = np.isfinite(pred) & np.isfinite(obs)
    pred, obs = pred[ok], obs[ok]
    if len(pred) < min_points:
        if len(pred) < 2:
            raise ValueError("RT calibration needs at least two anchors")
        a, b = np.polyfit(pred, obs, 1)
        return lambda x: a * np.asarray(x, float) + b
    o = np.argsort(pred, kind="mergesort")
    xs, ys = pred[o], obs[o]
    # collapse tied x to their mean so interpolation knots are unique
    ux, inv, cnt = np.unique(xs, return_inverse=True, return_counts=True)
    uy = np.bincount(inv, weights=ys) / cnt
    fy = pav(uy, cnt.astype(float))
    return lambda x: np.interp(np.asarray(x, float), ux, fy)


if __name__ == "__main__":
    # self-test: monotone recovery and robustness to one outlier block
    rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0, 1, 2000))
    y = 10 + 100 * x + rng.normal(0, 2, x.size)
    y[600:640] = 90.0  # a cluster eluting off-trend
    f = fit_isotonic(x, y)
    assert np.all(np.diff(f(np.linspace(0, 1, 200))) >= -1e-9)
    err = np.abs(f(x) - (10 + 100 * x))
    assert np.median(err) < 2.0, np.median(err)
    assert fit_isotonic(x[:5], y[:5])(0.5) is not None
    print("ok: median |err|", round(float(np.median(err)), 3))
