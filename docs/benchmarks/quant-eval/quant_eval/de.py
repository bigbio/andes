"""Differential expression on the protein table: moderated t, BH 5 % per condition pair, sign-aware true positives.

TP = significant protein of a species with a changing amount; TPs = TP with the sign of the observed log2FC equal to the expected
sign; FP = significant protein of the null species. Counts are given without a fold-change cut (suffix '') and with
|log2FC| >= 0.5 (suffix '5').

Protein families. The BH family of a pair is the set of proteins tested in it. Three families are scored:
  own    the proteins the arm itself tests (what each tool reports; family sizes differ between tools)
  inter  only proteins tested (>= 2 precursors) in every shared arm; the whole test (variance prior and BH) is run on that set,
         identical for all arms
  union  proteins tested in any shared arm; the prior is fitted on the arm's own proteins, an arm that does not test a protein
         contributes it with p = 1 to BH
The intersection family is the headline. Each family also gets the true positives at matched pooled FDP (sweep of the ranked p
values over all pairs, null-species hits counted as false) and the per-pair AUC of changed-vs-null species by p value."""
import numpy as np
import pandas as pd
from scipy import stats

from .stats import bh, fit_fdist

FC_CUT = 0.5


def de_table(P, r1, r2, prot_sp, expected_of):
    """P: protein x run log2. Returns per protein: species, observed log2FC (r2 - r1), p, BH q over the proteins of P, expected log2FC."""
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
    return pd.DataFrame({"sp": spv, "fc": fc, "p": p, "q": bh(p), "ex": ex}, index=prots)


def qstats(x):
    """Error summary of a vector: n, median |x|, share |x| > 1, median signed."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return {"n": int(len(x)), "med_abs": float(np.median(np.abs(x))) if len(x) else None,
            "share_gt1": float((np.abs(x) > 1).mean()) if len(x) else None,
            "med_signed": float(np.median(x)) if len(x) else None}


FDP_LEVELS = (0.01, 0.05, 0.10)


def _count(T, null, fam_q=None):
    """TP / TPs / FP with and without the fold-change cut from the q column of T."""
    chg = (T.sp.notna() & (T.sp != null)).values
    nul = (T.sp == null).values
    sgn = np.where(np.isfinite(T.ex.values), np.sign(T.ex.values), 1.0)
    fc = T.fc.values
    o = {}
    for tag, cut in (("", 0.0), ("5", FC_CUT)):
        sig = ((T.q.values <= 0.05) & (np.abs(np.nan_to_num(fc)) >= cut))
        o["TP" + tag] = int((sig & chg).sum())
        o["TPs" + tag] = int((sig & chg & (np.sign(fc) == sgn)).sum())
        o["FP" + tag] = int((sig & nul).sum())
    return o


def auc_changed_vs_null(T, null):
    """AUC of the p-value ranking (small p first) for separating proteins of changing species from the null species."""
    T = T[T.sp.notna()]
    lab = (T.sp != null).values
    n1, n0 = int(lab.sum()), int((~lab).sum())
    if n1 == 0 or n0 == 0:
        return None
    r = stats.rankdata(-T.p.values)
    return float((r[lab].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def tp_at_fdp(Ts, null, levels=FDP_LEVELS):
    """True positives at matched pooled FDP: rows of all pairs (species-assigned proteins) sorted by p; the largest prefix whose
    null-species share is <= the level. Returns {level: TP} and {level: FP in that prefix}."""
    T = pd.concat(Ts) if len(Ts) else pd.DataFrame({"sp": [], "p": []})
    T = T[T.sp.notna()].sort_values("p", kind="stable")
    nul = (T.sp == null).values
    cfp, ctp = np.cumsum(nul), np.cumsum(~nul)
    fdp = cfp / np.maximum(cfp + ctp, 1)
    tp, fp = {}, {}
    for lv in levels:
        k = np.where(fdp <= lv)[0]
        key = f"{int(round(lv * 100))}"
        tp[key] = int(ctp[k[-1]]) if len(k) else 0
        fp[key] = int(cfp[k[-1]]) if len(k) else 0
    return tp, fp


def union_view(T, fam, prot_sp):
    """T of one pair padded to the union family fam (missing proteins enter with p = 1); BH is recomputed on the resulting set."""
    miss = [p for p in fam if p not in T.index]
    pad = pd.DataFrame({"sp": [prot_sp.get(p) for p in miss], "fc": np.nan, "p": 1.0, "q": 1.0, "ex": np.nan}, index=miss)
    T = pd.concat([T, pad])
    T["q"] = bh(T.p.values)
    return T


def de_all(P, design, prot_sp, classes, fam_inter=None, fam_union=None):
    """Every condition pair (lo < hi) with >= 2 replicates on each side.
    Returns (totals of the own family, per-pair rows, log2FC error summaries by group|species, family summaries).
    fam_inter / fam_union: protein sets of the shared arms (see module doc); without them only the own family is scored.
    The per-pair rows and the totals are the own family; `de_family` holds all families."""
    runs = list(P.columns)
    cond = {r: design.cond[r] for r in runs if r in design.cond}
    levels = design.levels
    null = design.null
    tot = {k: 0 for k in ("TP", "TPs", "FP", "TP5", "TPs5", "FP5", "tested")}
    pairs, coll = [], {}
    P_inter = P[P.index.isin(fam_inter)] if fam_inter is not None else None
    fams = {"own": []}
    if fam_inter is not None:
        fams["inter"] = []
    if fam_union is not None:
        fams["union"] = []
    ftot = {f: {k: 0 for k in ("TP", "TPs", "FP", "TP5", "TPs5", "FP5", "tested", "tested_null")} for f in fams}
    aucs = {f: [] for f in fams}
    for i, lo in enumerate(levels):
        for hi in levels[i + 1:]:
            r1 = [r for r in runs if cond.get(r) == lo]
            r2 = [r for r in runs if cond.get(r) == hi]
            if len(r1) < 2 or len(r2) < 2:
                continue
            T = de_table(P, r1, r2, prot_sp, lambda s, lo=lo, hi=hi: design.expected_fc(s, lo, hi))
            chg = T.sp.notna() & (T.sp != null)
            tot["tested"] += int(len(T))
            row = {"lo": str(lo), "hi": str(hi), "tested": int(len(T))}
            c = _count(T, null)
            for k, v in c.items():
                tot[k] += v
                row[k] = v
            row["n_chg_tested"] = int(chg.sum())
            views = {"own": T}
            if fam_inter is not None:
                views["inter"] = de_table(P_inter, r1, r2, prot_sp, lambda s, lo=lo, hi=hi: design.expected_fc(s, lo, hi))
            if fam_union is not None:
                views["union"] = union_view(T, fam_union, prot_sp)
            for f, V in views.items():
                fams[f].append(V)
                cf = _count(V, null)
                for k, v in cf.items():
                    ftot[f][k] += v
                ftot[f]["tested"] += int(len(V))
                ftot[f]["tested_null"] += int((V.sp == null).sum())
                a = auc_changed_vs_null(V, null)
                aucs[f].append(a)
                row[f"auc_{f}"] = a
                if f != "own":
                    for k in ("TP", "FP", "TP5", "FP5"):
                        row[f"{k}_{f}"] = cf[k]
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
    fam_out = {}
    for f, Ts in fams.items():
        tp, fp = tp_at_fdp(Ts, null)
        av = [a for a in aucs[f] if a is not None]
        fam_out[f] = dict(ftot[f], tp_at_fdp=tp, fp_at_fdp=fp, auc_mean=float(np.mean(av)) if av else None,
                          auc_min=float(np.min(av)) if av else None, n_pairs_auc=len(av))
    return tot, pairs, comp, fam_out
