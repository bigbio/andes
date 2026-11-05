"""Matched true error: the quant_q cut at which an arm's measured error equals a reference's.

Two truth-measured errors (neither uses a model of the tool):
  wrong-peak rate   masked.py: share of re-found masked transfers whose apex is off the identified apex (pooled, per set)
  absence rate      rates.py: share of transferred rows in runs without a species that belong to that species (PXD014415: yeast in the 20 human-only runs)
For an arm with a quant_q column, the matched cut is the loosest quant_q limit whose error is not above the reference's error at the reference's
own operating point (quant_q <= 0.01 of the reference build) plus the reference's own seed-to-seed spread of that error (`tol`: max - min of the
reference build's rate over the quant-model seeds; a candidate whose error is above the reference's by less than the reference differs from itself
does not tighten the cut). With both errors available the cut is the smaller of the two. Depth, accuracy and DE
are then read at that cut (the `M` arm of the suite: identified rows + transfers with quant_q <= cut). A tool without a quant_q column is a
reference only (its rate is a target: `absence_cut_for_tool`)."""
from . import masked as MK
from . import rates as RT
from . import readers as R
from .fasta import map_precursors

NOMINAL = 0.01
MIN_ABSENT_TRANSFERS = 1000


def absence_curve(finalize_path, design, fasta, classes, species="YEAST"):
    """Absence curve of one andes finalize table: rate of `species` among transfers in the runs without it, per quant_q limit."""
    raw = R.read_andes(finalize_path, design.runs, "weights")
    sp, _, _ = map_precursors(sorted(set(raw.key)), fasta, classes)
    return RT.depth_curve(raw, design, sp, "thr")


def rate_on_curve(curve, species, tau=NOMINAL):
    p = min(curve, key=lambda x: abs(x["threshold"] - tau))
    return p[species]["rate"]


def absence_cut(curve, species, ref_rate):
    """Loosest quant_q limit whose absence rate is <= ref_rate (thresholds with too few absent-run transfers are skipped)."""
    best = None
    for p in curve:
        r = p[species]["rate"]
        if p[species]["transfers_in_absent_runs"] >= MIN_ABSENT_TRANSFERS and r == r and r <= ref_rate and (best is None or p["threshold"] > best):
            best = p["threshold"]
    return best


def wrong_cut(curve, ref_curve, tol=0.0):
    ref_rate = MK.rate_at(ref_curve, NOMINAL) + tol
    return MK.tau_for_rate(curve, ref_rate)


def matched_cut(build, ref_build, wrong_curves, abs_curves, species="YEAST", tol=None):
    """-> {"wrong": tau or None, "absence": tau or None, "tau": min of the available ones or None}. The reference build gets the nominal cut.
    tol = {"wrong": x, "absence": y}: the reference build's seed-to-seed spread of the wrong-peak rate and of the absence rate (0 when absent)."""
    tol = tol or {}
    if build == ref_build:
        return {"wrong": NOMINAL, "absence": NOMINAL, "tau": NOMINAL}
    o = {"wrong": None, "absence": None}
    if wrong_curves and build in wrong_curves and ref_build in wrong_curves:
        o["wrong"] = wrong_cut(wrong_curves[build], wrong_curves[ref_build], tol.get("wrong", 0.0))
    if abs_curves and build in abs_curves and ref_build in abs_curves:
        o["absence"] = absence_cut(abs_curves[build], species, rate_on_curve(abs_curves[ref_build], species) + tol.get("absence", 0.0))
    have = [v for v in o.values() if v is not None]
    o["tau"] = min(have) if have else None
    return o


def absence_cut_for_tool(curve, species, tool_rate):
    """Depth reading at a tool's false-transfer rate: the cut (and transferred rows) at which an andes arm's absence rate equals a tool's."""
    return RT.matched_depth(curve, species, tool_rate)
