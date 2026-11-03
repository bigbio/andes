"""Shared normalisation and protein roll-up, identical for every arm.

normalise  per-run median of log2 quantity over the precursors the arm quantifies in ALL runs (>= 90 % of runs if fewer than
           200 such precursors); no truth, all species.
rollup     Tukey median polish per protein over precursors x runs; only proteins with >= 2 precursors are tested later.
Precision-weight view: transfers carry a learned weight (learn_weights), the polish uses weighted row / column medians.
Null-species view: normalise(..., ref_keys=keys of the species whose amount is constant) takes the run factors from that species
only (an evaluation view that uses the design; applied to every arm the same way)."""
import numpy as np
import pandas as pd


def normalise(d, runs, ref_keys=None):
    """d: rows with quantity. Adds lq (log2) and n (run-centred log2). Returns (table, fallback flag, number of precursors used).
    ref_keys: optional set of keys the run medians are taken from (the null-species view passes the keys of the species whose
    amount is the same in every run); default = every key."""
    d = d.assign(lq=np.log2(d.quantity.values))
    flag = None if ref_keys is None else "null"
    src = d if ref_keys is None else d[d.key.isin(ref_keys)]
    nr = src.groupby("key").run.nunique()
    full = nr.index[nr == len(runs)]
    if len(full) < 200:
        full = nr.index[nr >= 0.9 * len(runs)]
        flag = "fb90" if ref_keys is None else "null-fb90"
    med = d[d.key.isin(full)].groupby("run").lq.median().reindex(runs)
    med = med.fillna(med.median())
    d["n"] = d.lq.values - d.run.map(med).values + float(med.median())
    return d, flag, int(len(full))


def normalise_rt_local(d, runs, rt_of, win=401):
    """Optional RT-local run normalisation, the same function for every arm (replaces `normalise`; design-free, no species).
    Reference precursors = those quantified in >= 90 % of runs; their residual to their own median over runs is smoothed per run
    along the precursor retention time (`rt_of`: key -> RT, one proxy for all arms) with a running median over `win` cells, and
    subtracted from every row of that run by interpolation at the row's RT. Rows whose key has no RT get the run's median shift."""
    d = d.assign(lq=np.log2(d.quantity.values))
    nr = d.groupby("key").run.nunique()
    ref = nr.index[nr >= 0.9 * len(runs)]
    mk = d[d.key.isin(ref)].groupby("key").lq.median()
    d["rt"] = d.key.map(rt_of)
    d["res"] = d.lq - d.key.map(mk)
    sm = pd.Series(0.0, index=d.index)
    for _, g in d.groupby("run"):
        x = g[g.key.isin(ref) & g.rt.notna() & g.res.notna()].sort_values("rt")
        if len(x) < 50:
            sm[g.index] = float(g.res.median()) if g.res.notna().any() else 0.0
            continue
        sx = x.res.rolling(win, center=True, min_periods=win // 4).median().bfill().ffill()
        sm[g.index] = np.where(g.rt.notna(), np.interp(g.rt.fillna(x.rt.median()).values, x.rt.values, sx.values),
                               float(x.res.median()))
    d["n"] = d.lq.values - sm.values
    return d.drop(columns=["rt", "res"]), "rtlocal", int(len(ref))


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


def _polish(M, Wt, weighted):
    return weighted_polish(M, Wt) if weighted else median_polish(M)


def censored_polish(M, Wt, T, sigma, weighted, iters=3):
    """Missing-not-at-random roll-up: transferred cells (mask T) are treated as right-censored at their observed value (the true
    value of a cell that was not identified is taken to be at most what was read). EM: polish, fitted value mu = row effect +
    protein value, then each transferred cell is replaced by E[y | y <= x] = mu - sigma phi(z) / Phi(z), z = (x - mu) / sigma;
    sigma = the arm's learned residual SD of identified cells."""
    from scipy.stats import norm
    X = M.copy()
    for _ in range(iters):
        out = _polish(X, Wt, weighted)
        re_ = np.nanmedian(X - out[None, :], axis=1)
        mu = re_[:, None] + out[None, :]
        z = (M - mu) / sigma
        lam = np.exp(norm.logpdf(z) - norm.logcdf(z))
        X = np.where(T & np.isfinite(M), mu - sigma * lam, M)
    return _polish(X, Wt, weighted)


def identified_sigma(d, prot_of, runs):
    """Robust SD of the protein median-polish residual of identified cells (proteins with >= 3 precursors), scaled by the
    two-way residual degrees of freedom of each protein's precursor x run table."""
    d = d.assign(pr=d.key.map(prot_of))
    ri = {r: i for i, r in enumerate(runs)}
    res = []
    for _, g in d.dropna(subset=["pr"]).groupby("pr"):
        keys = sorted(g.key.unique())
        if len(keys) < 3:
            continue
        ki = {k: i for i, k in enumerate(keys)}
        M = np.full((len(keys), len(runs)), np.nan)
        rr = np.array([ki[k] for k in g.key.values])
        cc = np.array([ri[r] for r in g.run.values])
        M[rr, cc] = g.n.values
        nk, nc = M.shape[0], int(np.isfinite(M).any(0).sum())
        if nc < 2:
            continue
        R = _polish_residual(M)[rr, cc] * np.sqrt(nk * nc / ((nk - 1.0) * (nc - 1.0)))   # residual degrees of freedom
        res.append(R[~g.is_transfer.values])
    return float(np.sqrt(_robust_var(np.concatenate(res)))) if res else 1.0


def rollup(d, prot_of, runs, censor_sigma=None):
    """d: rows with key, run, n (and optionally w). Returns a protein x run DataFrame (weighted polish when any w != 1).
    censor_sigma: when given, transferred cells are treated as right-censored (censored_polish) with this SD."""
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
            continue
        Wt = None
        if weighted:
            Wt = np.zeros(M.shape)
            Wt[rr, cc] = g.w.values
        if censor_sigma is not None and g.is_transfer.any():
            T = np.zeros(M.shape, bool)
            T[rr, cc] = g.is_transfer.values
            out[p] = censored_polish(M, Wt, T, censor_sigma, weighted)
        else:
            out[p] = _polish(M, Wt, weighted)
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
    if "se" in d.columns:
        # per-transfer uncertainty: the class variance scaled by the row's se**2 relative to the class mean of se**2
        s2 = d.se.values.astype(float) ** 2
        ok_se = tr & np.isfinite(s2) & (s2 > 0)
        if ok_se.any():
            w = d.w.values.copy()
            for c, v in var.items():
                k = ok_se & (cls == c)
                if c == 0 or not k.any():
                    continue
                rel = s2[k] / s2[k].mean()
                w[k] = np.clip(v0 / np.maximum(v * rel, 1e-6), 0.05, 1.0)
            d["w"] = w
    return d, {int(c): round(float(v), 4) for c, v in var.items()}


ROLLUPS = ("median-polish", "celltype", "balance", "celltype+balance", "censored")


def _robust_var(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float((1.4826 * np.median(np.abs(x - np.median(x)))) ** 2) if len(x) else np.nan


def celltype_offsets(d, cond_of):
    """Cell-type effect (selection of transfers): per precursor, the offset of its transferred cells from its identified cells,
    measured inside the conditions where the precursor has both kinds of cell, shrunk to a learned global offset.

    For precursor k and condition c with both kinds: diff = mean(n | transfer) - mean(n | identified), weight h = 1 / (1/n_T + 1/n_I).
    delta_k = sum(h diff) / H_k, sampling variance v_k = s2 / H_k, s2 = robust within-(precursor, condition, cell type) variance.
    Global offset delta0 = median of delta_k; between-precursor variance tau2 = max(0, robust var(delta_k) - median(v_k)).
    Offset_k = delta0 + tau2 / (tau2 + v_k) * (delta_k - delta0); precursors without a condition holding both kinds get delta0.
    Every scale is learned from the arm's own rows; no truth. d: normalised rows with key, run, n, is_transfer.
    Returns (offset per key of every transferred precursor, summary dict)."""
    x = d[["key", "run", "n", "is_transfer"]].copy()
    x["c"] = x.run.map(cond_of)
    x = x[x.c.notna() & np.isfinite(x.n.values)]
    tkeys = pd.Index(x.key[x.is_transfer].unique())
    info = {"transfer_precursors": int(len(tkeys))}
    if len(tkeys) == 0:
        return pd.Series(dtype=float), dict(info, delta0=0.0, tau2=0.0, s2=None, paired_precursors=0)
    grp = x.groupby(["key", "c", "is_transfer"]).n
    mu = grp.transform("mean")
    sz = grp.transform("size")
    ok = sz.values >= 2
    res = (x.n.values - mu.values)[ok] * np.sqrt(sz.values[ok] / (sz.values[ok] - 1.0))
    s2 = _robust_var(res)
    g = grp.agg(["mean", "size"]).unstack("is_transfer")
    if (True not in g["mean"].columns) or (False not in g["mean"].columns):
        return pd.Series(0.0, index=tkeys), dict(info, delta0=0.0, tau2=0.0, s2=s2, paired_precursors=0)
    mt, mi = g["mean"][True], g["mean"][False]
    nt, ni = g["size"][True], g["size"][False]
    both = mt.notna() & mi.notna()
    diff = (mt - mi)[both]
    h = (1.0 / (1.0 / nt[both] + 1.0 / ni[both]))
    kk = diff.index.get_level_values(0)
    H = h.groupby(kk).sum()
    dk = (h * diff).groupby(kk).sum() / H
    vk = s2 / H
    delta0 = float(np.median(dk.values)) if len(dk) else 0.0
    tau2 = max(0.0, _robust_var(dk.values) - float(np.median(vk.values))) if len(dk) > 1 else 0.0
    shrink = tau2 / (tau2 + vk) if tau2 > 0 else vk * 0.0
    off = pd.Series(delta0, index=tkeys)
    est = delta0 + shrink * (dk - delta0)
    off.loc[est.index.intersection(tkeys)] = est.reindex(est.index.intersection(tkeys)).values
    info.update(delta0=round(delta0, 4), tau2=round(float(tau2), 5), s2=round(float(s2), 5), paired_precursors=int(len(dk)),
                median_shrink=round(float(np.median(shrink.values)), 4) if len(dk) else None,
                offset_q10_q50_q90=[round(float(v), 4) for v in np.quantile(off.values, [0.1, 0.5, 0.9])])
    return off, info


def apply_celltype(d, off):
    """Remove the per-precursor cell-type offset from the transferred cells (column n; quantity when n is absent)."""
    d = d.copy()
    t = d.is_transfer.values
    o = d.key.map(off).fillna(0.0).values
    if "n" in d.columns:
        d.loc[t, "n"] = d.n.values[t] - o[t]
    if "lq" in d.columns:
        d.loc[t, "lq"] = d.lq.values[t] - o[t]
    return d


def adjust_raw(raw, off):
    """The same offset on a raw table (quantity / 2**offset for transferred cells), for row-level metrics of the arm."""
    raw = raw.copy()
    t = raw.is_transfer.values
    o = raw.key.map(off).fillna(0.0).values
    raw.loc[t, "quantity"] = raw.quantity.values[t] * np.exp2(-o[t])
    return raw


def balance_weights(d, prot_of, cond_of):
    """Within-condition balance: per (protein, condition) the transferred cells are re-weighted so that their share of the
    cell weight equals the protein's share pooled over all conditions (existing weights w, e.g. the learned precision of V,
    are the starting point). A condition that holds only one kind of cell is left as it is. Returns d with column w."""
    d = d.copy()
    w = d.w.values.astype(float) if "w" in d.columns else np.ones(len(d))
    p = d.key.map(prot_of)
    c = d.run.map(cond_of)
    t = d.is_transfer.values
    x = pd.DataFrame({"p": p.values, "c": c.values, "w": w, "wt": np.where(t, w, 0.0)})
    sp = x.groupby("p").wt.sum() / x.groupby("p").w.sum()
    spc = x.groupby(["p", "c"]).wt.sum() / x.groupby(["p", "c"]).w.sum()
    s_p = p.map(sp).values.astype(float)
    s_pc = spc.reindex(pd.MultiIndex.from_arrays([p.values, c.values])).values.astype(float)
    ok = t & np.isfinite(s_p) & np.isfinite(s_pc) & (s_pc > 0) & (s_pc < 1) & (s_p > 0) & (s_p < 1)
    a = np.ones(len(d))
    a[ok] = s_p[ok] * (1 - s_pc[ok]) / (s_pc[ok] * (1 - s_p[ok]))
    d["w"] = w * a
    return d
