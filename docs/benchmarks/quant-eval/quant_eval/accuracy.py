"""Identified-cell quantity accuracy of every arm on the same cells.

Cells = (precursor, run) identified (not transferred) by every arm. Per arm: run offsets from the shared null-species cells (the
species whose amount is constant), the arm's scale centred on the shared-cell median, error = (x - E) minus the mean of the
precursor's other runs (leave-run-out, at least 4 shared cells). The bright decile is the top decile of the design-predicted
intensity (leave-run-out, mean over the arms), so every arm is judged on the same cells.
Reported per arm: robust SD of the error, share |err| > 1 log2 overall and in the bright decile (gross errors), replicate SD (median
SD of the (precursor, condition) groups that have every replicate) overall and in the top predicted decile, UPS1 slope per precursor
(all levels, lower half), HYE gross ratio share, and the arm's own-row null-species replicate SD in its own top intensity decile.
E = design log2 amount of the precursor's species at the run's condition (`--accuracy-const SPECIES`: E = 0 in every run for that
species, one condition)."""
import math

import numpy as np
import pandas as pd


def _rsd(x):
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def cell_accuracy(tabs, kinds, design, sp, const, primary, min_cells=4):
    names = [n for n, k in kinds.items() if k["kind"] not in ("sage", "flashlfq")]   # context arms do not define the shared cells
    cols = {}
    for n in names:
        t = tabs[n]
        t = t[~t.is_transfer & (t.quantity > 0)]
        cols[n] = np.log2(t.set_index(["key", "run"]).quantity).rename(n)
    L = pd.concat(cols.values(), axis=1).reset_index()
    if const:
        null = const
        L["sp"] = L.key.map(sp)
        L = L[L.sp == const].copy()
        L["cond"] = "all"
        L["E"] = 0.0
    else:
        if not design.has_de:
            return {}
        null = design.null
        L["sp"] = L.key.map(sp)
        L["cond"] = L.run.map(design.cond)
        L["E"] = [(lambda v: np.nan if v is None else v)(design.expected(s, c)) for s, c in zip(L.sp, L.cond)]
        L = L[L.E.notna()].copy()
    T = list(names)
    sh = L[T].notna().all(axis=1)
    if not sh.any():
        return {}
    S = L[sh & (L.sp == null)]
    for t in T:
        pm = S.groupby("key")[t].transform("mean")
        o = (S[t] - pm).groupby(S.run).median()
        L[t] = L[t] - L.run.map(o).fillna(0.0)
    SH = L[sh].copy().reset_index(drop=True)
    gm = {t: float(SH[t].median()) for t in T}
    for t in T:
        SH[t + "_c"] = SH[t] - gm[t]
        d_ = SH[t + "_c"] - SH.E
        g = d_.groupby(SH.key)
        s_, n_ = g.transform("sum"), g.transform("count")
        SH[t + "_err"] = np.where(n_ >= min_cells, d_ - (s_ - d_) / (n_ - 1).clip(lower=1), np.nan)
    ref = [n for n in T if n == primary or kinds[n]["kind"] != "andes"]
    d_all = np.mean([SH[t + "_c"] - SH.E for t in ref], axis=0)
    g = pd.Series(d_all).groupby(SH.key)
    s_, n_ = g.transform("sum"), g.transform("count")
    SH["pred"] = np.where(n_ >= min_cells, (s_ - d_all) / (n_ - 1).clip(lower=1) + SH.E, np.nan)
    ok = np.isfinite(SH[[t + "_err" for t in T]]).all(axis=1) & np.isfinite(SH.pred)
    E_ = SH[ok]
    if len(E_) < 100:
        return {"cells": int(len(SH)), "err_cells": int(len(E_))}
    top = E_[E_.pred >= np.nanquantile(E_.pred, 0.9)]
    res = {"cells": int(len(SH)), "err_cells": int(len(E_)), "bright_cells": int(len(top)), "arms": {}}
    SH["grp"] = SH.key + "|" + SH.cond.astype(str)
    nrep = SH.groupby("grp").run.transform("nunique")
    maxrep = SH.groupby("cond").run.transform("nunique")
    G = SH[(nrep == maxrep) & (maxrep >= 3)]
    gb = G.groupby("grp")
    SH["gpred"] = SH.groupby("grp").pred.transform("mean")
    Gt = G.assign(gpred=SH.loc[G.index, "gpred"])
    if len(Gt) and Gt.gpred.notna().any():
        cut = np.nanquantile(Gt.groupby("grp").gpred.first(), 0.9)
        gbt = Gt[Gt.gpred >= cut].groupby("grp")
    else:
        gbt = None
    slopes = {}
    levels = []
    try:
        levels = sorted({float(c) for c in SH.cond.unique()})
    except ValueError:
        pass
    if "UPS" in design.classes and levels:
        U = SH[SH.sp == "UPS"].copy()
        U["lv"] = U.cond.astype(float)
        half = levels[: len(levels) // 2 + 1]
        for t in T:
            P = U.groupby(["key", "lv"])[t + "_c"].mean().unstack()
            for rn, ls in (("all", levels), ("low", half)):
                sl = []
                for k in P.index:
                    v = P.loc[k, [l for l in ls if l in P.columns]].dropna()
                    if len(v) >= max(3, len(ls) // 2 - 1):
                        sl.append(np.polyfit(np.log2(v.index.values.astype(float)), v.values, 1)[0])
                slopes[(t, rn)] = float(np.median(sl)) if sl else None
    ratio_bad = {}
    if not const and len(design.levels) >= 2 and len(design.classes) >= 3:
        for t in T:
            g_ = SH.groupby(["key", "cond"])[t].mean().unstack()
            c_ = SH.groupby(["key", "cond"]).size().unstack()
            spk = SH.groupby("key").sp.first()
            bad = []
            conds = list(g_.columns)
            for i, c1 in enumerate(conds):
                for c2 in conds[i + 1:]:
                    okk = (c_[c1] >= 2) & (c_[c2] >= 2)
                    r_ = (g_[c2] - g_[c1])[okk]
                    for s_n in design.classes:
                        e_ = design.expected_fc(s_n, c1, c2)
                        if e_ == e_:
                            bad += list(np.abs(r_[spk.reindex(r_.index) == s_n] - e_) > 1)
            ratio_bad[t] = float(100 * np.mean(bad)) if bad else None
    for t in T:
        e = E_[t + "_err"].values
        eb = top[t + "_err"].values
        r = {"rSD": _rsd(e), "bad": float(100 * np.mean(np.abs(e) > 1)), "bad_bright": float(100 * np.mean(np.abs(eb) > 1)),
             "within2x": float(100 * np.mean(np.abs(e) <= 1)),
             "rep_sd": float(gb[t].std().median()) if gb.ngroups else None, "groups": int(gb.ngroups),
             "rep_sd_top": float(gbt[t].std().median()) if gbt is not None and gbt.ngroups else None,
             "slope_all": slopes.get((t, "all")), "slope_low": slopes.get((t, "low")), "ratio_bad": ratio_bad.get(t)}
        x = L[[t, "key", "run", "cond", "sp"]].dropna(subset=[t])
        x = x[x.sp == null]
        x["grp"] = x.key + "|" + x.cond.astype(str)
        nr = x.groupby("grp").run.transform("nunique")
        mr = x.groupby("cond").run.transform("nunique")
        x = x[(nr == mr) & (mr >= 3)]
        if len(x):
            gg = x.groupby("grp")[t]
            m, sd = gg.mean(), gg.std()
            r["own_null_sd_top"] = float(sd[m >= m.quantile(0.9)].median())
            r["own_null_groups"] = int(len(m))
        res["arms"][t] = r
    return res
