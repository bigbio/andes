"""Truth-measured false-transfer rates from a species that is absent from some runs (PXD014415: yeast in the human-only runs).

A transferred row of the absent species in a run that does not contain it is a false transfer whose truth does not depend on any
model of the tool. Rates, per absent species s, with R_abs the runs without s and R_pres the other (non-blank) runs:

  rate            = (transferred rows of s in R_abs) / (all transferred rows in R_abs)       the false-transfer rate observed on foreign species
  yeast_ratio     = (transferred rows of s in R_abs) / (transferred rows of s in R_pres)     IonQuant-style: foreign hits against genuine hits of s
  rate_scaled     = rate / share_s                                                            extension: the invisible same-species errors are
                    assumed to occur in proportion to the species' share of the identified precursors (share_s); not a measured rate
  identified_in_abs = identified (MS/MS) rows of s in R_abs, reported next to the rates (search errors, not transfers)

Blank-run carry-over is reported by `carry_over` and is not part of the absence score: a transfer decoy method cannot see it.

`depth_curve` sweeps a quality column of the transfers (transfer_q or quant_q) and returns the false-transfer rate and the number of
transferred rows kept at each threshold; `matched_depth` reads the depth of a curve at a given rate (the most transferred rows among the thresholds not above that rate).
"""
import numpy as np

GRID = [0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1, 0.2, 0.5, 1.0]


def _runs(design, s):
    absent = {r for r, v in design.absent.items() if v == s}
    pres = set(design.runs) - absent - set(design.blank)
    return absent, pres


def false_transfer_rates(raw, design, sp):
    """Per absent species: counts and rates of one arm table (columns key, run, is_transfer)."""
    out = {}
    spk = raw.key.map(sp)
    for s in sorted(set(design.absent.values())):
        absent, pres = _runs(design, s)
        in_abs = raw.run.isin(absent)
        in_pres = raw.run.isin(pres)
        blank = raw.run.isin(design.blank)
        tr = raw.is_transfer.astype(bool) & ~blank
        n_tr_abs = int((in_abs & tr).sum())
        n_false = int((in_abs & tr & (spk == s)).sum())
        n_s_pres = int((in_pres & tr & (spk == s)).sum())
        ident = raw[~raw.is_transfer.astype(bool) & ~blank]
        share = float((ident.drop_duplicates("key").key.map(sp) == s).mean()) if len(ident) else float("nan")
        o = {"absent_runs": len(absent), "present_runs": len(pres), "transfers_in_absent_runs": n_tr_abs, "false_transfers": n_false,
             "genuine_transfers_of_species": n_s_pres,
             "rate": n_false / n_tr_abs if n_tr_abs else float("nan"),
             "yeast_ratio": n_false / n_s_pres if n_s_pres else float("nan"),
             "species_share_identified": share,
             "rate_scaled": (n_false / n_tr_abs / share) if (n_tr_abs and share > 0) else float("nan"),
             "identified_in_absent_runs": int((in_abs & ~raw.is_transfer.astype(bool) & (spk == s)).sum()),
             "transfers_all_runs": int(tr.sum())}
        out[s] = o
    return out


def carry_over(raw, design, sp):
    """Rows in blank runs (carry-over), reported separately from the absence score."""
    if not design.blank:
        return None
    b = raw[raw.run.isin(design.blank)]
    return {"blank_runs": len(design.blank), "rows": int(len(b)), "transferred": int(b.is_transfer.sum()), "precursors": int(b.key.nunique())}


def depth_curve(raw, design, sp, col, grid=GRID):
    """Sweep `col` (tq = transfer_q or thr = quant_q) over the transfers of one arm table that carries every transfer (rule `all`/`weights`).
    Rows: threshold, transferred rows kept (all non-blank runs), per absent species the false transfers, transfers in absent runs, rate."""
    spk = raw.key.map(sp)
    tr = raw.is_transfer.astype(bool) & ~raw.run.isin(design.blank)
    q = raw[col].where(np.isfinite(raw[col]), np.inf)
    out = []
    for t in grid:
        kept = tr & (q <= t)
        row = {"threshold": t, "transfers": int(kept.sum())}
        for s in sorted(set(design.absent.values())):
            absent, _ = _runs(design, s)
            ia = kept & raw.run.isin(absent)
            n = int(ia.sum())
            f = int((ia & (spk == s)).sum())
            row[s] = {"transfers_in_absent_runs": n, "false": f, "rate": f / n if n else float("nan")}
        out.append(row)
    return out


def matched_depth(curve, species, rate, min_n=1000):
    """Most transferred rows of the curve among the thresholds whose false-transfer rate is <= rate (a step on the grid, no interpolation);
    on equal rows the smallest threshold. Thresholds with fewer than `min_n` transferred rows in the absent runs are skipped (their rate is
    noise). Returns (threshold, transfers) or None when no such threshold has a rate <= `rate`."""
    best = None
    for p in curve:
        r = p[species]["rate"]
        if p[species]["transfers_in_absent_runs"] >= min_n and r == r and r <= rate and (best is None or p["transfers"] > best[1]):
            best = (p["threshold"], p["transfers"])
    return best
