"""Absence classes and blank-run split of one arm table.

Absent species (PXD014415: yeast in the 20 human-only runs). A row of the species in a run that does not contain it is false, but about
two thirds of such rows sit on precursors whose intensity is the same in the runs that do contain the species: real signal of a human
(or shared) peptide that the species map calls yeast, which is label noise and not a search or transfer error. The class of a precursor
comes from its own cross-run profile:

  d = median over the runs that contain the species of the run-normalised log2 intensity  minus  the same median over the absent runs

after a per-run normalisation by the median of the control species (the species present in every run; its precursors have d ~ 0).
Classes of the rows of the absent species in absent runs:
  yeast_enriched    d above the boundary: intensity is higher in the runs that contain the species (a genuine species peptide, absent here)
  not_specific      d at or below the boundary: the precursor looks like a control precursor (same intensity with and without the species)
  only_absent       the precursor has no row in a run that contains the species (no reference)
  not_in_reference  (reference table only) the precursor has no row in the reference table
The profile and the boundary can be derived once from a reference table (`reference_profiles`, the andes every-row table) and every arm classified
against it, so that the class of a precursor does not depend on the rows an arm reports. The enriched count is a lower bound of the errors that are
specific to the species: the boundary is wide for a precursor seen in few absent runs, and the rows of the only-absent class cannot be judged; the
board prints the enriched rate and the enriched + only-absent rate.
The boundary is learned from the control species: the 1 - alpha quantile of d over its precursors, within bins of the number of
rows the precursor has in the absent runs (d of a precursor seen in one absent run is noisier than of one seen in twenty). The fixed
boundary of 1 log2 (2x) is reported next to it.

Blank runs. Rows of a blank run split into identified and transferred. A row is carry-over-like when its raw intensity is below a
fraction of the median intensity of the same precursor in the other (non-blank) runs; the fraction is learned from the ordinary runs: the
`q` quantile of the ratio of a non-blank row to the median of the same precursor's rows in the other non-blank runs (a blank row below
the level that only q of ordinary rows fall under). Both the learning and the classification need the precursor in at least two other runs; a
blank row of a precursor in fewer runs is counted as `no_reference`. The fixed fraction 0.5 is reported next to it."""
import warnings

import numpy as np
import pandas as pd

ALPHA = 0.01
FIXED_LOG2 = 1.0
EDGES = (1, 2, 5, 10)
MIN_CONTROL = 100
BLANK_Q = 0.05
BLANK_FIXED = 0.5
MIN_REF_RUNS = 2


def profile(raw, design, sp, species, control=None):
    """Per precursor with rows in the arm: species, median normalised log2 in the runs that contain `species` (med_pres, n_pres) and in its absent runs (med_abs, n_abs), d = med_pres - med_abs.
    The control species defaults to the most frequent species other than `species`."""
    absent = {r for r, v in design.absent.items() if v == species}
    present = set(design.runs) - absent - set(design.blank)
    d = raw[raw.run.isin(absent | present) & (raw.quantity > 0)][["key", "run", "quantity"]].copy()
    d["sp"] = d.key.map(sp)
    d["lq"] = np.log2(d.quantity)
    if control is None:
        vc = d[(d.sp != species) & d.sp.notna() & (d.sp != "NONE")].sp.value_counts()
        control = vc.index[0] if len(vc) else None
    cm = d[d.sp == control].groupby("run").lq.median()
    d["lqn"] = d.lq - d.run.map(cm)
    d["ab"] = d.run.isin(absent)
    g = d.groupby(["key", "sp", "ab"]).lqn.agg(["median", "size"]).unstack("ab")
    # unstack gives columns (stat, ab): ensure both ab values exist
    for stat in ("median", "size"):
        for ab in (False, True):
            if (stat, ab) not in g.columns:
                g[(stat, ab)] = np.nan
    g = pd.DataFrame({"med_pres": g[("median", False)], "n_pres": g[("size", False)].fillna(0), "med_abs": g[("median", True)],
                      "n_abs": g[("size", True)].fillna(0)}).reset_index()
    g["d"] = g.med_pres - g.med_abs
    return g, control


def learn_cutoffs(d_control, n_abs_control, alpha=ALPHA, edges=EDGES, min_n=MIN_CONTROL):
    """Boundary of d per bin of n_abs: the 1 - alpha quantile of the control's d. Returns a list of (lower edge of the bin, cutoff, n control precursors,
    pooled_fallback). A bin with fewer than `min_n` control precursors takes the cutoff of the pooled control (pooled_fallback True; the pooled
    quantile is not conditional on the bin and is only as wide as the whole control's spread). None when the control is empty."""
    d = np.asarray(d_control, float)
    n = np.asarray(n_abs_control, float)
    ok = np.isfinite(d)
    d, n = d[ok], n[ok]
    if len(d) == 0:
        return None
    pooled = float(np.quantile(d, 1 - alpha))
    out = []
    lo = list(edges)
    for i, a in enumerate(lo):
        b = lo[i + 1] if i + 1 < len(lo) else np.inf
        m = (n >= a) & (n < b)
        enough = int(m.sum()) >= min_n
        out.append((a, float(np.quantile(d[m], 1 - alpha)) if enough else pooled, int(m.sum()), not enough))
    return out


def cutoff_for(cutoffs, n_abs):
    c = cutoffs[0][1]
    for a, v, *_ in cutoffs:
        if n_abs >= a:
            c = v
    return c


def classify(prof_s, cutoffs=None, fixed=None):
    """Class of every precursor row of the species profile: 'only_absent', 'yeast_enriched', 'not_specific' (cutoffs: learned bins; fixed: a single log2 boundary)."""
    cls = pd.Series("not_specific", index=prof_s.index, dtype=object)
    cls[prof_s.n_pres == 0] = "only_absent"
    if fixed is not None:
        enr = prof_s.d > fixed
    else:
        enr = pd.Series([bool(np.isfinite(dd)) and dd > cutoff_for(cutoffs, na) for dd, na in zip(prof_s.d, prof_s.n_abs)], index=prof_s.index)
    cls[enr & (prof_s.n_pres > 0)] = "yeast_enriched"
    return cls


CLASSES = ("yeast_enriched", "not_specific", "only_absent", "not_in_reference")


def _class_counts(rows, cls_of):
    """Rows per class; a key the class map does not hold (a precursor missing from the reference table) is `not_in_reference`."""
    c = rows.key.map(cls_of).fillna("not_in_reference")
    return {k: int((c == k).sum()) for k in CLASSES}


def _cutoff_rows(cuts):
    return [{"n_abs_from": a, "cutoff": v, "n_control": n, "pooled_fallback": f} for a, v, n, f in (cuts or [])]


def reference_profiles(ref_raw, design, sp, alpha=ALPHA, fixed=FIXED_LOG2):
    """Class maps of every absent species from ONE reference table (the andes every-row table: every row of the finalize output): d and the learned
    cutoffs are derived once, so the class of a precursor does not depend on which rows an arm reports. {species: {control, cutoffs, learned: {key: class}, fixed: {key: class}, ...}}."""
    out = {}
    for s in sorted(set(design.absent.values())):
        prof, control = profile(ref_raw, design, sp, s)
        if control is None:
            continue
        ctrl = prof[(prof.sp == control) & (prof.n_pres > 0) & (prof.n_abs > 0)]
        cuts = learn_cutoffs(ctrl.d.values, ctrl.n_abs.values, alpha)
        if cuts is None:
            continue
        ps = prof[prof.sp == s].reset_index(drop=True)
        out[s] = {"control": control, "alpha": alpha, "fixed_log2": fixed, "control_precursors": int(len(ctrl)),
                  "control_d_quantiles": {str(q): float(np.quantile(ctrl.d, q)) for q in (0.05, 0.5, 0.95, 0.99)},
                  "cutoffs": cuts, "reference_rows": int(len(ref_raw)),
                  "learned": dict(zip(ps.key, classify(ps, cuts))), "fixed": dict(zip(ps.key, classify(ps, None, fixed)))}
    return out


def absence_classes(raw, design, sp, alpha=ALPHA, fixed=FIXED_LOG2, refs=None):
    """Per absent species: identified and transferred rows of the species in the absent runs by class, with the learned and the fixed boundary.
    With `refs` (reference_profiles) the classes come from the reference table for every arm (`learned`, `fixed`) and the classes from the arm's
    own rows are kept as `own_learned` / `own_fixed`; without it `learned` / `fixed` are the arm's own (profile and cutoffs from its rows).
    The yeast-enriched count is a lower bound of the species-specific errors: a precursor seen in few absent runs needs a large d to be called."""
    out = {}
    for s in sorted(set(design.absent.values())):
        absent = {r for r, v in design.absent.items() if v == s}
        prof, control = profile(raw, design, sp, s)
        if control is None:
            continue
        ctrl = prof[(prof.sp == control) & (prof.n_pres > 0) & (prof.n_abs > 0)]
        cuts = learn_cutoffs(ctrl.d.values, ctrl.n_abs.values, alpha)
        ps = prof[prof.sp == s].reset_index(drop=True)
        rows = raw[raw.run.isin(absent) & (raw.key.map(sp) == s) & ~raw.run.isin(design.blank)]
        tr = rows.is_transfer.astype(bool)
        o = {"control": control, "alpha": alpha, "fixed_log2": fixed, "control_precursors": int(len(ctrl)),
             "control_d_quantiles": {str(q): float(np.quantile(ctrl.d, q)) for q in (0.05, 0.5, 0.95, 0.99)} if len(ctrl) else {},
             "cutoffs": _cutoff_rows(cuts)}
        if cuts is None:
            out[s] = o
            continue
        for name, cl in (("learned", classify(ps, cuts)), ("fixed", classify(ps, None, fixed))):
            cls_of = dict(zip(ps.key, cl))
            o[("own_" + name) if refs else name] = {"identified": _class_counts(rows[~tr], cls_of), "transferred": _class_counts(rows[tr], cls_of)}
        if refs and s in refs:
            r = refs[s]
            o["reference"] = {"control": r["control"], "alpha": r["alpha"], "fixed_log2": r["fixed_log2"], "control_precursors": r["control_precursors"],
                              "control_d_quantiles": r["control_d_quantiles"], "cutoffs": _cutoff_rows(r["cutoffs"]), "reference_rows": r["reference_rows"]}
            for name in ("learned", "fixed"):
                o[name] = {"identified": _class_counts(rows[~tr], r[name]), "transferred": _class_counts(rows[tr], r[name])}
        # denominators: all rows (any species) of the arm in the absent runs
        allabs = raw[raw.run.isin(absent) & ~raw.run.isin(design.blank)]
        o["denominators"] = {"identified_rows_absent_runs": int((~allabs.is_transfer.astype(bool)).sum()),
                             "transferred_rows_absent_runs": int(allabs.is_transfer.astype(bool).sum())}
        out[s] = o
    return out


def blank_split(raw, design, q=BLANK_Q, fixed=BLANK_FIXED):
    """Blank rows split into identified and transferred, each into carry-over-like / comparable / no reference (module docstring).
    {} without blank runs. The learned fraction and the reference counts are returned with it."""
    if not design.blank:
        return {}
    d = raw[raw.quantity > 0][["key", "run", "quantity", "is_transfer"]]
    nb = d[~d.run.isin(design.blank)]
    M = nb.pivot(index="key", columns="run", values="quantity")
    cols = list(M.columns)
    arr = M.values
    ratios = []
    for j in range(len(cols)):
        others = np.delete(arr, j, axis=1)
        cnt = np.isfinite(others).sum(1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(np.where(np.isfinite(others), others, np.nan), axis=1)
        ok = np.isfinite(arr[:, j]) & (cnt >= 2)
        ratios.append(arr[ok, j] / med[ok])
    ratios = np.concatenate(ratios) if ratios else np.array([])
    frac = float(np.quantile(ratios, q)) if len(ratios) else None
    med_other = nb.groupby("key").quantity.median()
    n_other = nb.groupby("key").run.nunique()
    b = d[d.run.isin(design.blank)].copy()
    b["n_ref"] = b.key.map(n_other).fillna(0)
    b["ref"] = b.key.map(med_other).where(b.n_ref >= MIN_REF_RUNS)       # a reference needs the precursor in at least MIN_REF_RUNS other runs (as the learning of the fraction)
    b["ratio"] = b.quantity / b.ref
    out = {"blank_runs": len(design.blank), "rows": int(len(b)), "learned_fraction": frac, "fixed_fraction": fixed, "quantile": q,
           "ordinary_rows": int(len(ratios))}
    for name, sel in (("identified", ~b.is_transfer.astype(bool)), ("transferred", b.is_transfer.astype(bool))):
        x = b[sel]
        has = x.ref.notna()
        o = {"rows": int(len(x)), "no_reference": int((~has).sum()),
             "median_ratio": float(x.ratio.median()) if has.any() else None}
        for tag, f in (("learned", frac), ("fixed", fixed)):
            o[f"carry_over_{tag}"] = int((has & (x.ratio < f)).sum()) if f is not None else None
        out[name] = o
    out["carry_over_learned"] = sum(out[n]["carry_over_learned"] or 0 for n in ("identified", "transferred"))
    out["carry_over_fixed"] = sum(out[n]["carry_over_fixed"] for n in ("identified", "transferred"))
    return out
