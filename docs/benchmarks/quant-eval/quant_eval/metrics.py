"""Precursor accuracy, counts, absence, per-level UPS1 and replicate-SD metrics."""
import numpy as np
import pandas as pd

from .de import qstats
from .keys import bare


def errors(d, ref_ident_only=True):
    """Leave-one-condition-out precursor error. d: rows with key, cond, c (run-centred log2), E (expected log2 amount),
    is_transfer. For a row of condition c the prediction is the mean over the other conditions c' of the mean of (c - E) over the
    precursor's reference rows in c' (identified rows when ref_ident_only). err = (c - E) - prediction; NaN when not evaluable."""
    adj = d.c.values - d.E.values
    base = d[~d.is_transfer.values] if ref_ident_only else d
    t = pd.DataFrame({"key": base.key.values, "cond": base.cond.values, "adj": base.c.values - base.E.values})
    t = t[np.isfinite(t.adj.values)]
    g = t.groupby(["key", "cond"]).adj.mean()
    S = g.groupby(level=0).sum()
    N = g.groupby(level=0).size()
    own = g.reindex(pd.MultiIndex.from_arrays([d.key.values, d.cond.values])).values
    s = S.reindex(d.key.values).values
    n = N.reindex(d.key.values).values.astype(float)
    have = np.isfinite(own)
    s2 = np.where(have, s - own, s)
    n2 = np.where(have, n - 1, n)
    pred = np.where(n2 > 0, s2 / np.maximum(n2, 1), np.nan)
    return adj - pred


def accuracy(d, design, sp, shared, own_rows):
    """Precursor error by species and by UPS level. d: normalised rows (column n). own_rows: reference = the arm's identified rows
    on its own rows; otherwise only the cells in `shared` (cells every shared arm quantifies) are scored and all rows are reference."""
    x = d.assign(cond=d.run.map(design.cond), sp=d.key.map(sp), c=d.n)
    x["E"] = [(lambda v: np.nan if v is None else v)(design.expected(s, c)) for s, c in zip(x.sp.values, x.cond.values)]
    if not own_rows:
        x = x[(x.key + "|" + x.run).isin(shared)]
    x = x[np.isfinite(x.E.values)].copy()
    if len(x) == 0:
        return {}
    x["e"] = errors(x, ref_ident_only=own_rows)
    o = {s: qstats(g.e.values) for s, g in x.groupby("sp")}
    o["_all"] = qstats(x.e.values)
    if "UPS" in design.classes:
        for c, g in x[x.sp == "UPS"].groupby("cond"):
            o[f"UPS@{c}"] = qstats(g.e.values)
    return o


def basic_counts(raw, prot, sp, prot_sp, classes):
    """Rows, precursors, peptides and proteins (identified and transfer-only), overall and by species."""
    o = {"rows": int(len(raw)), "ident_rows": int((~raw.is_transfer).sum()), "transfer_rows": int(raw.is_transfer.sum())}
    ik = set(raw.key[~raw.is_transfer])
    ak = set(raw.key)
    o["precursors"] = len(ak)
    o["precursors_ident"] = len(ik)
    o["precursors_transfer_only"] = len(ak - ik)
    o["peptides"] = int(pd.Series(list(ak)).map(bare).nunique())
    o["peptides_ident"] = int(pd.Series(list(ik)).map(bare).nunique())
    pr = raw.assign(p=raw.key.map(prot)).dropna(subset=["p"])
    cnt = pr.groupby("p").key.nunique()
    idp = set(pr.p[~pr.is_transfer])
    o["proteins_ge1"] = int(len(cnt))
    o["proteins_ge2"] = int((cnt >= 2).sum())
    o["proteins_ge1_ident"] = int(len(idp))
    o["proteins_ge1_transfer_only"] = int(len(cnt) - len(idp))
    o["by_species"] = {}
    for s in classes:
        kk = {k for k in ak if sp.get(k) == s}
        ps = [p for p in cnt.index if prot_sp.get(p) == s]
        o["by_species"][s] = {"precursors": len(kk), "precursors_ident": len(kk & ik), "precursors_transfer_only": len(kk - ik),
                              "proteins_ge1": len(ps), "proteins_ge2": int(sum(cnt[p] >= 2 for p in ps)),
                              "proteins_ge1_ident": len([p for p in ps if p in idp])}
    return o


def ups_by_level(raw, P, design, sp, prot, prot_sp):
    """Per UPS1 level: precursors / proteins detected (identified, transfer-only), proteins quantified in >= 2 replicates, transfers."""
    out = {}
    pr = raw.assign(p=raw.key.map(prot), cd=raw.run.map(design.cond))
    ups_p = {p for p, s in prot_sp.items() if s == "UPS"}
    Pu = P[P.index.isin(ups_p)]
    for L in design.levels:
        rl = [r for r in design.runs if design.cond.get(r) == L]
        u = pr[(pr.cd == L) & (pr.key.map(sp) == "UPS")]
        ik = set(u.key[~u.is_transfer])
        tk = set(u.key[u.is_transfer])
        gp = u.dropna(subset=["p"])
        ip = set(gp.p[~gp.is_transfer])
        tp = set(gp.p[gp.is_transfer]) - ip
        q2 = int((np.isfinite(Pu[rl].values).sum(1) >= 2).sum()) if len(Pu) else 0
        out[str(L)] = {"n_runs": len(rl), "ups_precursors_ident": len(ik), "ups_precursors_transfer_only": len(tk - ik),
                       "ups_proteins_ident": len(ip), "ups_proteins_transfer_only": len(tp), "ups_proteins_quant_ge2reps": q2,
                       "ups_rows_ident": int((~u.is_transfer).sum()), "ups_rows_transferred": int(u.is_transfer.sum()),
                       "ups_transfer_rows_without_ident_at_level": int(((u.is_transfer) & (~u.key.isin(ik))).sum())}
    return out


def design_detection(raw, design, sp, thresholds=(1.0, 1.58), ref_skip=3):
    """Design-consistent UPS1 detection per level, next to raw detection (any row of the precursor at the level).
    Per arm: run offsets = per-run median of log2 over the null-species precursors present in >= 90 % of runs; each UPS precursor
    gets its own reference from the levels above the `ref_skip` lowest ones (median over those rows of log2 - log2 level); a row is
    design-consistent when |log2 - (reference + log2 level)| <= t. Precursors without a reference row are not scored.
    thresholds: log2 ratio limits (1.0 = within 2x, 1.58 = within 3x). Returns {level: {"any", "t<thr>", "ident", "cells", ...}}."""
    if "UPS" not in design.classes or design.null is None or len(design.levels) <= ref_skip + 1:
        return {}
    T = raw[["key", "run", "quantity", "is_transfer"]].copy()
    T = T[T.quantity > 0]
    T["lq"] = np.log2(T.quantity)
    T["sp"] = T.key.map(sp)
    ny = T[T.sp == design.null].groupby("key").run.nunique()
    full = ny.index[ny >= 0.9 * len(design.runs)]
    med = T[T.key.isin(full)].groupby("run").lq.median()
    if len(med) == 0:
        return {}
    T["n"] = T.lq - T.run.map(med) + med.median()
    U = T[T.sp == "UPS"].copy()
    U["lev"] = U.run.map(lambda r: float(design.cond[r]) if r in design.cond else np.nan)
    U = U[np.isfinite(U.lev)]
    ref_levels = {float(x) for x in design.levels[ref_skip:]}
    hi = U[U.lev.isin(ref_levels)]
    a = (hi.n - np.log2(hi.lev)).groupby(hi.key).median()
    U = U[U.key.isin(a.index)].copy()
    U["err"] = U.n - (U.key.map(a) + np.log2(U.lev))
    out = {}
    for L in design.levels:
        u = U[U.lev == float(L)]
        o = {"any": int(u.key.nunique()), "ident": int(u[~u.is_transfer].key.nunique()), "cells": int(len(u)),
             "reference_levels_from": str(design.levels[ref_skip])}
        for t in thresholds:
            o[f"t{t}"] = int(u[u.err.abs() <= t].key.nunique())
        out[str(L)] = o
    return out


def absence(raw, design, sp):
    """Rows that should not exist: every row in blank runs, and rows of a species in the runs where it is absent."""
    o = {}
    if design.blank:
        b = raw[raw.run.isin(design.blank)]
        nb = raw[~raw.run.isin(design.blank)]
        o["blank"] = {"blank_runs": len(design.blank), "blank_rows": int(len(b)), "blank_transferred": int(b.is_transfer.sum()),
                      "blank_precursors": int(b.key.nunique()), "other_rows": int(len(nb)), "other_precursors": int(nb.key.nunique())}
    for s in sorted(set(design.absent.values())):
        runs_s = {r for r, v in design.absent.items() if v == s}
        h = raw[raw.run.isin(runs_s)]
        spk = h.key.map(sp)
        o.setdefault("absent_species", {})[s] = {"absent_runs": len(runs_s), "rows": int((spk == s).sum()),
                                                  "transferred": int(((spk == s) & h.is_transfer).sum()),
                                                  "other_species_rows": int((spk != s).sum())}
    return o


def pooled_sd(M, cond_of):
    """Pooled within-condition SD (ddof 1) per row of M (rows x runs, log2), from conditions with >= 2 finite replicates."""
    ss = np.zeros(len(M))
    df = np.zeros(len(M))
    for c in sorted(set(cond_of.values()), key=str):
        cols = [r for r in M.columns if cond_of.get(r) == c]
        if len(cols) < 2:
            continue
        X = M[cols].values
        n = np.isfinite(X).sum(1)
        m = np.nanmean(X, 1)
        v = np.nansum((X - m[:, None]) ** 2, 1)
        ok = n >= 2
        ss += np.where(ok, v, 0)
        df += np.where(ok, n - 1, 0)
    return pd.Series(np.where(df > 0, np.sqrt(ss / np.maximum(df, 1)), np.nan), index=M.index)


def cv_pct(sd):
    """CV in percent of a log2-scale SD."""
    return 100 * np.sqrt(np.expm1((np.log(2) * sd) ** 2))


def matched_sd(PM, ref_names, design, sp, prot_sp, classes, runs):
    """Replicate SD at matched features: proteins (polished) and precursors shared by the reference arms, in rank-matched
    intensity quartiles (mean rank across the reference arms). PM: arm -> (protein table P, normalised rows d)."""
    res = {"ref_arms": ref_names}
    for lvl in ("protein", "precursor"):
        mats = {}
        for n, (P, d) in PM.items():
            mats[n] = P if lvl == "protein" else d.pivot_table(index="key", columns="run", values="n", aggfunc="mean").reindex(columns=runs)
        common = sorted(set.intersection(*[set(mats[n].index) for n in ref_names]))
        res[lvl] = {"n_common": len(common)}
        if len(common) < 40:
            continue
        rk = pd.DataFrame({n: mats[n].loc[common].mean(axis=1).rank(pct=True) for n in ref_names}).mean(axis=1)
        qb = pd.qcut(rk, 4, labels=False, duplicates="drop")
        spx = pd.Series([(sp.get(k) if lvl == "precursor" else prot_sp.get(k)) for k in common], index=common)
        for n, M in mats.items():
            sd = pooled_sd(M.reindex(common), design.cond)
            o = {"all": {"n": int(sd.notna().sum()), "med_sd": float(sd.median()), "med_cv_pct": float(cv_pct(sd).median())}}
            for qq in sorted(qb.dropna().unique()):
                s_ = sd[(qb == qq).values]
                o[f"Q{int(qq) + 1}"] = {"n": int(s_.notna().sum()), "med_sd": float(s_.median()), "med_cv_pct": float(cv_pct(s_).median())}
            if "UPS" in classes or len(classes) >= 3:
                for s_n in (classes if len(classes) >= 3 else ("UPS", "YEAST")):
                    ss_ = sd[(spx == s_n).values]
                    o[f"sp_{s_n}"] = {"n": int(ss_.notna().sum()), "med_sd": float(ss_.median()) if ss_.notna().any() else None}
            res[lvl][n] = o
    return res


def detection_by_quartile(ROWS, PM, refs, sp, classes, all_rows_refs=()):
    """Share of the reference arm's precursors found in each arm, by species and reference-intensity quartile
    (reference intensity = per-key median log2 of the reference arm's identified rows, or of all its rows for the arms in
    all_rows_refs; Q1 = lowest)."""
    out = {}
    for refn in refs:
        d = PM[refn][1]
        di = d if refn in all_rows_refs else d[~d.is_transfer]
        med = di.groupby("key").n.median()
        spk = med.index.map(sp)
        tab = {}
        for s in classes:
            m = med[(spk == s)]
            if len(m) < 8:
                continue
            qb = pd.qcut(m.rank(method="first"), 4, labels=False)
            tab[s] = {}
            for n, raw in ROWS.items():
                det = set(raw.key)
                idk = set(raw.key[~raw.is_transfer])
                tab[s][n] = [{"q": int(q) + 1, "n_ref": int((qb == q).sum()), "detected": int(m.index[qb == q].isin(det).sum()),
                              "detected_ident": int(m.index[qb == q].isin(idk).sum())} for q in range(4)]
        out[refn] = tab
    return out
