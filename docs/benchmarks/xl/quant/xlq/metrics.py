"""Scores of a cross-link quantification table against a design with a known answer.

Every arm goes through the same steps, so that tools are compared on their quantities and not
on their downstream defaults:

1. Roll-up of precursors into residue pairs: Tukey median polish of log2 quantities over the
   precursors x runs of a residue pair (the protein roll-up of docs/benchmarks/quant-eval), or the
   sum of the precursor quantities of a run.
2. No normalisation unless the dataset asks for it (datasets.json `normalise`): in a dilution
   series of cross-linked protein the change is the signal, and a global normalisation would
   remove it. When it is on, `run_shifts` centres the runs by quant-eval's rule.
3. The statistics of quant-eval: moderated t (empirical Bayes variance prior) and
   Benjamini-Hochberg.
"""
import functools
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "quant-eval"))
from quant_eval.de import de_table  # noqa: E402
from quant_eval.rollup import median_polish  # noqa: E402

DE_Q = 0.05
WITHIN = 1.0
MIN_REPLICATES = 2
QUANT_Q_LEVELS = (0.01, 0.05)
NORMALISE_MIN_COMPLETE = 200
NORMALISE_FALLBACK_SHARE = 0.9


def _quiet(fn):
    """Run `fn` without the NaN RuntimeWarnings of all-missing rows and columns, which are expected."""
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return fn(*args, **kwargs)
    return wrapper


@_quiet
def log_matrix(df, runs, rollup="median-polish"):
    """Residue pair x run log2 matrix (NaN where a run has no value)."""
    out = {}
    for rp, sub in df.groupby("residue_pair"):
        P = sub.pivot(index="precursor", columns="run", values="quantity").reindex(columns=runs)
        L = np.log2(P.to_numpy(dtype=float))
        if rollup == "sum":
            s = np.nansum(P.to_numpy(dtype=float), axis=0)
            row = np.where(np.isfinite(L).any(axis=0), np.log2(np.where(s > 0, s, np.nan)), np.nan)
        elif rollup == "median-polish":
            row = median_polish(L)
        else:
            raise ValueError(f"unknown roll-up {rollup!r}")
        out[rp] = row
    return pd.DataFrame.from_dict(out, orient="index", columns=runs).astype(float)


@_quiet
def precursor_matrix(df, runs):
    """Precursor x run log2 matrix."""
    P = df.pivot(index="precursor", columns="run", values="quantity").reindex(columns=runs)
    return np.log2(P.astype(float))


@_quiet
def run_shifts(M, runs):
    """Per-run log2 shift that centres the runs of a residue pair x run matrix.

    The rule of quant-eval's `normalise`, on a matrix: the run medians are taken over the rows with
    a value in every run (in at least 90 % of the runs when fewer than 200 rows have all), and the
    shift of a run is its median minus the median of the run medians. Subtract the result from the
    residue-pair and the precursor matrices alike.
    """
    if M.empty:
        return pd.Series(0.0, index=runs)
    n = np.isfinite(M[runs].to_numpy(dtype=float)).sum(axis=1)
    rows = n == len(runs)
    if rows.sum() < NORMALISE_MIN_COMPLETE:
        rows = n >= NORMALISE_FALLBACK_SHARE * len(runs)
    if not rows.any():
        raise ValueError("normalise: no residue pair has a value in at least 90 % of the runs")
    med = M.loc[rows, runs].median(axis=0)
    med = med.fillna(med.median())
    return med - float(med.median())


def depth(df, M):
    """Counts of precursors, residue pairs and rows, identified and transferred."""
    ident = df[df["evidence"] == "identified"]
    return {
        "precursors": int(df["precursor"].nunique()),
        "residue_pairs": int(np.isfinite(M.to_numpy(dtype=float)).any(axis=1).sum()),
        "residue_pairs_identified": int(ident["residue_pair"].nunique()),
        "rows": int(len(df)),
        "rows_identified": int(len(ident)),
        "rows_transferred": int((df["evidence"] == "transferred").sum()),
    }


def _row_mean(X):
    """Mean of the finite values of each row; NaN for a row without any."""
    fin = np.isfinite(X)
    n = fin.sum(axis=1)
    s = np.where(fin, X, 0.0).sum(axis=1)
    return np.where(n > 0, s / np.maximum(n, 1), np.nan)


def _row_cv(V):
    """Sample CV (ddof 1) of the finite values of each row; NaN for a row with fewer than two."""
    fin = np.isfinite(V)
    n = fin.sum(axis=1)
    m = _row_mean(V)
    ss = np.where(fin, (V - m[:, None]) ** 2, 0.0).sum(axis=1)
    sd = np.sqrt(ss / np.maximum(n - 1, 1))
    return np.where(n >= MIN_REPLICATES, sd / m, np.nan)


def _summary(err):
    err = np.asarray(err, float)
    if len(err) == 0:
        return {"n": 0, "median_error": None, "median_abs_error": None, "within_2x": None}
    return {
        "n": int(len(err)),
        "median_error": round(float(np.median(err)), 4),
        "median_abs_error": round(float(np.median(np.abs(err))), 4),
        "within_2x": round(float(np.mean(np.abs(err) <= WITHIN)), 4),
    }


@_quiet
def dilution(M, level_of, reference_level):
    """Accuracy of the log2 ratios of each dilution level against the reference level.

    For a residue pair with at least MIN_REPLICATES values at a level and at the reference, the
    observed log2 ratio is the difference of the replicate means and the expected one is
    log2(level / reference). Also: the least-squares slope of observed on expected (1 means no
    ratio compression; a shift common to every level, such as a reference level that reads low,
    goes into the intercept) and, per level, the share of residue pairs called changed at BH 5 %
    with the expected sign. `run_offsets` shows whether an injection differs from its replicates.
    """
    runs = list(M.columns)
    per_level_runs = {}
    for r in runs:
        per_level_runs.setdefault(level_of[r], []).append(r)
    if reference_level not in per_level_runs:
        raise ValueError(f"dilution: no run has the reference level {reference_level:g} "
                         f"(levels: {', '.join(f'{v:g}' for v in sorted(per_level_runs))})")
    thin = [f"{lvl:g}" for lvl, rs in sorted(per_level_runs.items()) if len(rs) < MIN_REPLICATES]
    if thin:
        raise ValueError(f"dilution: fewer than {MIN_REPLICATES} runs at level(s) {', '.join(thin)}")
    ref = per_level_runs[reference_level]
    levels = sorted(set(per_level_runs) - {reference_level})
    X = M.to_numpy(dtype=float)
    ref_cols = [runs.index(r) for r in ref]
    ok_ref = np.isfinite(X[:, ref_cols]).sum(axis=1) >= MIN_REPLICATES
    mean_ref = _row_mean(X[:, ref_cols])
    per_level, xs, ys, at_all = {}, [], [], ok_ref.copy()
    for lvl in levels:
        lvl_runs = per_level_runs[lvl]
        cols = [runs.index(r) for r in lvl_runs]
        ok = ok_ref & (np.isfinite(X[:, cols]).sum(axis=1) >= MIN_REPLICATES)
        at_all &= ok
        expected = math.log2(lvl / reference_level)
        obs = _row_mean(X[ok][:, cols]) - mean_ref[ok]
        s = _summary(obs - expected)
        s["expected_log2"] = round(expected, 4)
        T = de_table(M, ref, lvl_runs, {}, lambda _: np.nan)
        hit = (T["q"] <= DE_Q) & (np.sign(T["fc"]) == np.sign(expected))
        s["de_tested"] = int(len(T))
        s["de_called_expected_sign"] = int(hit.sum())
        s["de_share"] = round(float(hit.mean()), 4) if len(T) else None
        per_level[f"{lvl:g}"] = s
        xs.extend([expected] * int(ok.sum()))
        ys.extend(obs.tolist())
    slope = float(np.polyfit(xs, ys, 1)[0]) if len(set(xs)) > 1 else None
    return {
        "reference_level": reference_level,
        "levels": per_level,
        "residue_pairs_at_every_level": int(at_all.sum()),
        "slope": round(slope, 4) if slope is not None else None,
    }


@_quiet
def precision(M, groups):
    """Replicate CV (linear scale) and completeness per group of replicate runs."""
    out = {}
    for name, runs in groups.items():
        V = 2 ** M[runs].to_numpy(dtype=float)
        n = np.isfinite(V).sum(axis=1)
        present = n > 0
        complete = n == len(runs)
        cv = _row_cv(V)
        with_cv = complete & (n >= MIN_REPLICATES)
        out[name] = {
            "runs": len(runs),
            "quantified": int(present.sum()),
            "complete": int(complete.sum()),
            "completeness": round(float(complete.sum() / present.sum()), 4) if present.any() else None,
            "median_cv_complete": round(float(np.median(cv[with_cv])), 4) if with_cv.any() else None,
            "median_cv_2plus": round(float(np.median(cv[n >= MIN_REPLICATES])), 4) if (n >= MIN_REPLICATES).any() else None,
        }
    return out


@_quiet
def run_offsets(P, groups):
    """Injection-to-injection offsets inside each group of replicate runs.

    P is the precursor x run log2 matrix. The offset of a run is the median, over the precursors
    with at least MIN_REPLICATES values in its group, of the run's value minus that precursor's
    mean over the group. `spread` is the largest minus the smallest offset of the group: a large
    spread names an injection whose amount or signal differed from its replicates.
    """
    out = {}
    for name, runs in groups.items():
        X = P[runs].to_numpy(dtype=float)
        X = X[np.isfinite(X).sum(axis=1) >= MIN_REPLICATES]
        if not len(X):
            out[name] = {"offsets": {r: None for r in runs}, "spread": None}
            continue
        D = X - _row_mean(X)[:, None]
        offsets = {}
        for i, r in enumerate(runs):
            col = D[np.isfinite(D[:, i]), i]
            offsets[r] = round(float(np.median(col)), 4) if len(col) else None
        seen = [v for v in offsets.values() if v is not None]
        out[name] = {"offsets": offsets, "spread": round(max(seen) - min(seen), 4) if seen else None}
    return out


def split_halves(runs):
    """Two halves of a list of runs in injection order: 1st, 4th, 5th, 8th, ... against 2nd, 3rd,
    6th, 7th, ... (ABBA). A drift that is linear in the injection order cancels between them when
    the number of runs is a multiple of four, and is a fifth of one injection step for ten runs;
    an odd/even split would leave a whole step."""
    first = [r for i, r in enumerate(runs) if i % 4 in (0, 3)]
    second = [r for i, r in enumerate(runs) if i % 4 in (1, 2)]
    return first, second


@_quiet
def null_calls(M, groups):
    """Differential calls between two halves of each replicate group: all are false.

    The halves are those of `split_halves`, so that drift over the injection order is not a
    difference between them; the runs of a group must be listed in injection order.
    """
    out = {}
    for name, runs in groups.items():
        if len(runs) < 2 * MIN_REPLICATES:
            continue
        first, second = split_halves(runs)
        T = de_table(M, first, second, {}, lambda _: np.nan)
        out[name] = {"tested": int(len(T)), "false_calls": int((T["q"] <= DE_Q).sum())}
    return out


def transfers(df, absent):
    """False-transfer rate where the truth says a residue pair cannot be in a run.

    `absent(row)` returns True (the residue pair cannot be in that run's sample), False (it can) or
    None (no truth for that run). Rate = transferred rows of absent residue pairs / transferred rows
    in the runs with a truth. Identified rows of absent residue pairs are identification errors and
    are counted apart. With a `quant_q` column the rate is also given at fixed quant_q limits.
    """
    flags = df.apply(absent, axis=1)
    has = flags.notna()
    d = df[has].assign(absent=flags[has].astype(bool))
    tr = d[d["evidence"] == "transferred"]
    ident = d[d["evidence"] == "identified"]
    out = {
        "runs_with_truth": int(d["run"].nunique()),
        "transferred": int(len(tr)),
        "false_transfers": int(tr["absent"].sum()),
        "false_transfer_rate": round(float(tr["absent"].mean()), 4) if len(tr) else None,
        "identified": int(len(ident)),
        "identified_absent": int(ident["absent"].sum()),
    }
    if "quant_q" in d.columns and len(tr):
        q = pd.to_numeric(tr["quant_q"], errors="coerce")
        out["at_quant_q"] = {
            f"{lim:g}": {"transferred": int((q <= lim).sum()),
                         "false_transfer_rate": round(float(tr.loc[q <= lim, "absent"].mean()), 4) if (q <= lim).any() else None}
            for lim in QUANT_Q_LEVELS
        }
    return out
