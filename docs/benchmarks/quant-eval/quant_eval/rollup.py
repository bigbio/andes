"""Shared normalisation and protein roll-up, identical for every arm.

normalise  per-run median of log2 quantity over the precursors the arm quantifies in ALL runs (>= 90 % of runs if fewer than
           200 such precursors); no truth, all species.
rollup     Tukey median polish per protein over precursors x runs; only proteins with >= 2 precursors are tested later.
Precision-weight view: transfers carry a learned weight (learn_weights), the polish uses weighted row / column medians."""
import numpy as np
import pandas as pd


def normalise(d, runs):
    """d: rows with quantity. Adds lq (log2) and n (run-centred log2). Returns (table, fallback flag, number of precursors used)."""
    d = d.assign(lq=np.log2(d.quantity.values))
    flag = None
    nr = d.groupby("key").run.nunique()
    full = nr.index[nr == len(runs)]
    if len(full) < 200:
        full = nr.index[nr >= 0.9 * len(runs)]
        flag = "fb90"
    med = d[d.key.isin(full)].groupby("run").lq.median().reindex(runs)
    med = med.fillna(med.median())
    d["n"] = d.lq.values - d.run.map(med).values + float(med.median())
    return d, flag, int(len(full))


def median_polish(M, iters=10):
    """M: precursors x runs (NaN allowed). Returns the protein value per run, NaN where a run has no value."""
    X = M.copy()
    ov = 0.0
    re_ = np.zeros(X.shape[0])
    ce = np.zeros(X.shape[1])
    for _ in range(iters):
        r = np.nan_to_num(np.nanmedian(X, axis=1))
        X = X - r[:, None]
        re_ += r
        mr = np.median(ce)
        ce -= mr
        ov += mr
        c = np.nan_to_num(np.nanmedian(X, axis=0))
        X = X - c[None, :]
        ce += c
        mc = np.median(re_)
        re_ -= mc
        ov += mc
    out = ov + ce
    out[np.all(~np.isfinite(M), axis=0)] = np.nan
    return out


def _wmedian(v, w):
    """Weighted median; equal weights give numpy's median (mean of the two middle values at an even count)."""
    o = np.argsort(v, kind="stable")
    v = v[o]
    w = w[o]
    c = np.cumsum(w)
    h = 0.5 * c[-1]
    i = int(np.searchsorted(c, h, side="left"))
    if i + 1 < len(v) and abs(c[i] - h) <= 1e-12 * c[-1]:
        return 0.5 * (v[i] + v[i + 1])
    return v[min(i, len(v) - 1)]


def _wnanmedian(X, Wt, axis):
    out = np.zeros(X.shape[0] if axis == 1 else X.shape[1])
    for k in range(len(out)):
        x = X[k] if axis == 1 else X[:, k]
        w = Wt[k] if axis == 1 else Wt[:, k]
        ok = np.isfinite(x) & (w > 0)
        out[k] = _wmedian(x[ok], w[ok]) if ok.any() else 0.0
    return out


def weighted_polish(M, Wt, iters=10):
    """Median polish with cell weights (row / column medians weighted; the median of the effects stays unweighted)."""
    X = M.copy()
    ov = 0.0
    re_ = np.zeros(X.shape[0])
    ce = np.zeros(X.shape[1])
    for _ in range(iters):
        r = _wnanmedian(X, Wt, 1)
        X = X - r[:, None]
        re_ += r
        mr = np.median(ce)
        ce -= mr
        ov += mr
        c = _wnanmedian(X, Wt, 0)
        X = X - c[None, :]
        ce += c
        mc = np.median(re_)
        re_ -= mc
        ov += mc
    out = ov + ce
    out[np.all(~np.isfinite(M), axis=0)] = np.nan
    return out


def rollup(d, prot_of, runs):
    """d: rows with key, run, n (and optionally w). Returns a protein x run DataFrame (weighted polish when any w != 1)."""
    d = d.assign(p=d.key.map(prot_of)).dropna(subset=["p"])
    weighted = "w" in d.columns and (d.w.values != 1.0).any()
    out = {}
    ri = {r: i for i, r in enumerate(runs)}
    for p, g in d.groupby("p"):
        keys = sorted(g.key.unique())
        ki = {k: i for i, k in enumerate(keys)}
        M = np.full((len(keys), len(runs)), np.nan)
        rr = [ki[k] for k in g.key.values]
        cc = [ri[r] for r in g.run.values]
        M[rr, cc] = g.n.values
        if len(keys) < 2:
            out[p] = M[0]
        elif weighted:
            Wt = np.zeros(M.shape)
            Wt[rr, cc] = g.w.values
            out[p] = weighted_polish(M, Wt)
        else:
            out[p] = median_polish(M)
    return pd.DataFrame.from_dict(out, orient="index", columns=runs)


def _polish_residual(M, iters=10):
    X = M.copy()
    for _ in range(iters):
        X = X - np.nan_to_num(np.nanmedian(X, axis=1))[:, None]
        X = X - np.nan_to_num(np.nanmedian(X, axis=0))[None, :]
    return X


def learn_weights(d, prot_of, runs):
    """Truth-free learned precision of the transfers. The variance of the protein median-polish residual (proteins with >= 3
    precursors) is estimated by decile of the transfer error probability p (identified cells form their own class); the weight
    is identified variance / class variance, clipped to [0.05, 1]. No fold-change or species information is used.
    d: normalised rows with key, run, n, is_transfer, p. Returns (d with column w, variance per class)."""
    d = d.copy()
    d["pr"] = d.key.map(prot_of)
    ri = {r: i for i, r in enumerate(runs)}
    res = pd.Series(np.nan, index=d.index)
    for _, g in d.dropna(subset=["pr"]).groupby("pr"):
        keys = sorted(g.key.unique())
        if len(keys) < 3:
            continue
        ki = {k: i for i, k in enumerate(keys)}
        M = np.full((len(keys), len(runs)), np.nan)
        rr = [ki[k] for k in g.key.values]
        cc = [ri[r] for r in g.run.values]
        M[rr, cc] = g.n.values
        res.loc[g.index] = _polish_residual(M)[rr, cc]
    tr = d.is_transfer.values & d.p.notna().values
    cls = np.zeros(len(d), int)
    if tr.any():
        cls[tr] = 1 + np.minimum((pd.Series(d.p.values[tr]).rank(pct=True).values * 10).astype(int), 9)
    r = res.values
    ok = np.isfinite(r)
    var = {c: (1.4826 * np.median(np.abs(r[ok & (cls == c)] - np.median(r[ok & (cls == c)])))) ** 2
           for c in range(11) if (ok & (cls == c)).sum() >= 50}
    v0 = var.get(0, min(var.values()) if var else 1.0)
    wc = {c: float(np.clip(v0 / max(v, 1e-6), 0.05, 1.0)) for c, v in var.items()}
    d["w"] = [wc.get(c, 1.0) for c in cls]
    return d, {int(c): round(float(v), 4) for c, v in var.items()}
