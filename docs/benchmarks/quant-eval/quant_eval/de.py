"""Differential expression on the protein table: moderated t, BH 5 % per condition pair, sign-aware true positives.

TP = significant protein of a species with a changing amount; TPs = TP with the sign of the observed log2FC equal to the expected
sign; FP = significant protein of the null species. Counts are given without a fold-change cut (suffix '') and with
|log2FC| >= 0.5 (suffix '5')."""
import numpy as np
import pandas as pd
from scipy import stats

from .stats import bh, fit_fdist

FC_CUT = 0.5


def de_table(P, r1, r2, prot_sp, expected_of):
    """P: protein x run log2. Returns per protein: species, observed log2FC (r2 - r1), BH q, expected log2FC."""
    A = P[r1].values
    B = P[r2].values
    n1 = np.isfinite(A).sum(1)
    n2 = np.isfinite(B).sum(1)
    ok = (n1 >= 2) & (n2 >= 2)
    A, B, n1, n2 = A[ok], B[ok], n1[ok], n2[ok]
    prots = P.index[ok]
    m1 = np.nanmean(A, 1)
    m2 = np.nanmean(B, 1)
    v = np.nansum((A - m1[:, None]) ** 2, 1) + np.nansum((B - m2[:, None]) ** 2, 1)
    df = (n1 + n2 - 2).astype(float)
    s2 = v / np.maximum(df, 1)
    d0, s0 = fit_fdist(s2, df)
    if np.isinf(d0):
        s2t = np.full_like(s2, s0)
        dft = np.full_like(df, 1e6)
    else:
        s2t = (d0 * s0 + df * s2) / (d0 + df)
        dft = df + d0
    fc = m2 - m1
    t = fc / np.sqrt(s2t * (1 / n1 + 1 / n2))
    p = np.where(np.isfinite(t), 2 * stats.t.sf(np.abs(t), dft), 1.0)
    spv = np.array([prot_sp.get(x) for x in prots], dtype=object)
    ex = np.array([expected_of(s) if s is not None else np.nan for s in spv], float)
    return pd.DataFrame({"sp": spv, "fc": fc, "q": bh(p), "ex": ex})


def qstats(x):
    """Error summary of a vector: n, median |x|, share |x| > 1, median signed."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return {"n": int(len(x)), "med_abs": float(np.median(np.abs(x))) if len(x) else None,
            "share_gt1": float((np.abs(x) > 1).mean()) if len(x) else None,
            "med_signed": float(np.median(x)) if len(x) else None}


def de_all(P, design, prot_sp, classes):
    """Every condition pair (lo < hi) with >= 2 replicates on each side.
    Returns (totals, per-pair rows, log2FC error summaries by group|species)."""
    runs = list(P.columns)
    cond = {r: design.cond[r] for r in runs if r in design.cond}
    levels = design.levels
    null = design.null
    tot = {k: 0 for k in ("TP", "TPs", "FP", "TP5", "TPs5", "FP5", "tested")}
    pairs, coll = [], {}
    for i, lo in enumerate(levels):
        for hi in levels[i + 1:]:
            r1 = [r for r in runs if cond.get(r) == lo]
            r2 = [r for r in runs if cond.get(r) == hi]
            if len(r1) < 2 or len(r2) < 2:
                continue
            T = de_table(P, r1, r2, prot_sp, lambda s, lo=lo, hi=hi: design.expected_fc(s, lo, hi))
            chg = T.sp.notna() & (T.sp != null)
            sgn = np.where(np.isfinite(T.ex), np.sign(T.ex), 1.0)
            tot["tested"] += int(len(T))
            row = {"lo": str(lo), "hi": str(hi), "tested": int(len(T))}
            for tag, cut in (("", 0.0), ("5", FC_CUT)):
                sig = (T.q <= 0.05) & (T.fc.abs() >= cut)
                tp = int((sig & chg).sum())
                tps = int((sig & chg & (np.sign(T.fc) == sgn)).sum())
                fp = int((sig & (T.sp == null)).sum())
                tot["TP" + tag] += tp
                tot["TPs" + tag] += tps
                tot["FP" + tag] += fp
                row["TP" + tag] = tp
                row["FP" + tag] = fp
            row["n_chg_tested"] = int(chg.sum())
            for s in classes:
                m = (T.sp == s).values & np.isfinite(T.ex.values)
                if m.sum():
                    row["err_" + s] = qstats(T.fc.values[m] - T.ex.values[m])
                    g = "low" if design.is_low(lo) else "high"
                    if s == "UPS" or len(levels) == 2:
                        coll.setdefault((g, s), []).append((T.ex.values[m], T.fc.values[m]))
                    coll.setdefault(("all", s), []).append((T.ex.values[m], T.fc.values[m]))
            pairs.append(row)
    comp = {}
    for (g, s), L in coll.items():
        ex = np.concatenate([a for a, _ in L])
        fc = np.concatenate([b for _, b in L])
        o = qstats(fc - ex)
        o["n_prot_pairs"] = int(len(ex))
        if np.nanstd(ex) > 1e-9 and len(ex) > 10:
            o["slope_fc_on_ex"] = float(np.polyfit(ex, fc, 1)[0])
        nz = np.abs(ex) > 1e-9
        if nz.any():
            o["median_obs_over_exp"] = float(np.median(fc[nz] / ex[nz]))
        comp[f"{g}|{s}"] = o
    return tot, pairs, comp
