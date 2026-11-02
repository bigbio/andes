#!/usr/bin/env python3
"""Scorecard tables (markdown) from the result files of downstream.py.

  python -m quant_eval.scorecard --results DS1.json DS2.json ... > scorecard.md

One section per dataset: andes (the first andes arm) against MaxQuant and FragPipe, with the margin (andes minus tool; absolute
and % of the tool) and a verdict per tool.

Verdict rule for counts and errors: `equal` when |margin| <= max(3 % of the tool value, 2 units for integer counts; the
tolerance is 10 % for FP and absence counts); otherwise
`better` / `worse` by the metric direction (lower is better for FP, FDP, error and SD). FDP rows use an absolute band of 2
points. The DE summary table uses the TP / FDP rule: better = TP >= 0.98 x tool TP and FDP lower by > 2 points, or
TP > 1.03 x tool TP and FDP not higher by > 2 points; worse = the mirror; else equal.

DE rows use the intersection BH family (proteins tested by every shared arm; variance prior and BH on that set) as the headline; the family
table gives own / intersection / union. True positives at matched pooled FDP (1 / 5 / 10 %) and per-pair AUC come from the ranked
p-value list. Absence (species absent from runs), the false-transfer rates and blank-run carry-over have their own tables; they are not part of the DE or detection rows."""
import argparse
import json

import numpy as np


def g(d, *ks):
    for k in ks:
        if d is None:
            return None
        d = d.get(k) if isinstance(d, dict) else None
    return d


def fmt(x, nd=0):
    if x is None:
        return "n/a"
    if isinstance(x, float) and nd == 0 and abs(x) < 10:
        nd = 2
    return f"{x:,.{nd}f}"


def verdict(a, t, hi=True, tol=0.03, floor=0.0, abs_tol=None):
    if a is None or t is None:
        return "n/a"
    band = abs_tol if abs_tol is not None else max(abs(t) * tol, floor)
    if abs(a - t) <= band:
        return "equal"
    return "better" if ((a - t) > 0) == hi else "worse"


def margin(a, t, nd=0, pct=True):
    if a is None or t is None:
        return "n/a"
    d = a - t
    s = f"{d:+,.{nd if nd else (2 if abs(d) < 10 else 0)}f}"
    if pct and t not in (0, None) and abs(t) > 1e-12:
        s += f" ({100 * d / abs(t):+.1f}%)"
    return s


def de_verdict(tp, fp, tp0, fp0):
    f, f0 = 100 * fp / max(tp + fp, 1), 100 * fp0 / max(tp0 + fp0, 1)
    if (tp >= 0.98 * tp0 and f < f0 - 2) or (tp > 1.03 * tp0 and f <= f0 + 2):
        return "better"
    if (tp0 >= 0.98 * tp and f0 < f - 2) or (tp0 > 1.03 * tp and f0 <= f + 2):
        return "worse"
    return "equal"


class Section:
    def __init__(self, res):
        self.r = res
        self.A = res["arms"]
        roles = res["arm_roles"]
        self.C = res["primary"]
        self.MQ = next((n for n, k in roles.items() if k["kind"] == "maxquant"), None)
        self.FP = next((n for n, k in roles.items() if k["kind"] == "fragpipe"), None)
        self.views = [n for n, k in roles.items() if k["kind"] == "andes" and n != self.C]
        self.rows = []
        self.fam = "inter" if "inter" in (g(self.A.get(self.C), "de_family") or {}) else "own"

    def dv(self, arm, key):
        """DE count of an arm in the family used for the headline rows."""
        if arm is None or arm not in self.A:
            return None
        return g(self.A[arm], "de", key) if self.fam == "own" else g(self.A[arm], "de_family", self.fam, key)

    def val(self, arm, *ks):
        return g(self.A.get(arm), *ks) if arm else None

    def row(self, name, a, m, p, hi=True, tol=0.03, nd=0, abs_tol=None, pct=True):
        floor = 2.0 if (nd == 0 and abs_tol is None) else 0.0
        va, vp = verdict(a, m, hi, tol, floor, abs_tol), verdict(a, p, hi, tol, floor, abs_tol)
        self.rows.append(f"| {name} | {fmt(a, nd)} | {fmt(m, nd)} | {fmt(p, nd)} | {margin(a, m, nd, pct)} | {margin(a, p, nd, pct)} | MQ {va} / FP {vp} |")

    def r3(self, name, ks, **kw):
        self.row(name, self.val(self.C, *ks), self.val(self.MQ, *ks), self.val(self.FP, *ks), **kw)

    def r3pct(self, name, ks, **kw):
        v = [self.val(x, *ks) for x in (self.C, self.MQ, self.FP)]
        self.row(name, *[100 * x if x is not None else None for x in v], **kw)

    def head(self, text):
        self.rows.append(f"| **{text}** | | | | | | |")

    def build(self):
        A, C, r = self.A, self.C, self.r
        classes = r.get("classes", [])
        self.head("1 Identification")
        for lab, k in (("precursors (all rows)", "precursors"), ("precursors identified (MS/MS)", "precursors_ident"),
                       ("precursors transfer-only", "precursors_transfer_only"), ("rows identified", "ident_rows"),
                       ("rows transferred", "transfer_rows"), ("peptides", "peptides"), ("proteins (>= 1 precursor)", "proteins_ge1"),
                       ("proteins identified (>= 1 MS/MS precursor)", "proteins_ge1_ident"), ("proteins >= 2 precursors", "proteins_ge2")):
            self.r3(lab, (k,))
        for s in classes:
            self.r3(f"{s}: precursors", ("by_species", s, "precursors"))
            self.r3(f"{s}: proteins >= 2 precursors", ("by_species", s, "proteins_ge2"))
        if "ups_strat" in A[C]:
            self.head("2 UPS1 by level (detected = any replicate; quantified = protein with >= 2 precursors in >= 2 replicates of the level)")
            for L in r["levels"]:
                nrep = A[C]["ups_strat"][L]["n_runs"]
                tot = lambda arm, a_, b_: (None if self.val(arm, "ups_strat", L) is None else self.val(arm, "ups_strat", L, a_) + self.val(arm, "ups_strat", L, b_))
                self.row(f"L={L}: UPS proteins detected (ident + transferred)", *[tot(x, "ups_proteins_ident", "ups_proteins_transfer_only") for x in (C, self.MQ, self.FP)])
                self.r3(f"L={L}: UPS proteins quantified in >= 2 of {nrep} reps", ("ups_strat", L, "ups_proteins_quant_ge2reps"))
                self.row(f"L={L}: UPS precursors detected", *[tot(x, "ups_precursors_ident", "ups_precursors_transfer_only") for x in (C, self.MQ, self.FP)])
                self.r3(f"L={L}: UPS precursors identified", ("ups_strat", L, "ups_precursors_ident"))
        if g(A[C], "ups_design"):
            ud = A[C]["ups_design"]
            ref0 = next(iter(ud.values()))["reference_levels_from"]
            self.head(f"2b UPS1 precursor detection, design-consistent (own reference = median over levels >= {ref0}; run offsets from the null-species precursors; "
                      f"precursors without a reference row not scored; within 2x = |log2 error| <= 1, within 3x <= 1.58)")
            for L in r["levels"]:
                for lab, key in (("detected, any row", "any"), ("design-consistent, within 2x", "t1.0"), ("design-consistent, within 3x", "t1.58")):
                    self.row(f"L={L}: UPS precursors {lab}", *[g(A.get(x), "ups_design", L, key) for x in (C, self.MQ, self.FP)])
        for refn, tab in (r.get("hye_det") or {}).items():
            self.head(f"2 detection by species and intensity quartile (share of reference-arm precursors found; reference = {refn}; Q1 = lowest)")
            for s, arms in tab.items():
                for qi in range(4):
                    sh = lambda arm: None if arm not in arms else 100.0 * arms[arm][qi]["detected"] / max(arms[arm][qi]["n_ref"], 1)
                    self.row(f"{s} Q{qi + 1} (n={arms[C][qi]['n_ref']}) share detected %", sh(C), sh(self.MQ), sh(self.FP), nd=1, pct=False, abs_tol=2.0)
        if "acc_own_sp" in A[C]:
            self.head("3 Accuracy (precursor error = leave-one-condition-out log2 ratio error; protein error = log2FC - expected)")
            for s in [k for k in A[C]["acc_own_sp"] if not k.startswith("_") and not k.startswith("UPS@")]:
                self.r3pct(f"precursor own rows, {s}: share abs(err) > 1 %", ("acc_own_sp", s, "share_gt1"), hi=False, nd=2)
                self.r3(f"precursor own rows, {s}: median abs(err)", ("acc_own_sp", s, "med_abs"), hi=False, nd=3)
            self.r3pct("precursor shared cells, all species: share abs(err) > 1 %", ("acc_shared_sp", "_all", "share_gt1"), hi=False, nd=2)
            if "UPS" in classes:
                for L in r["levels"]:
                    self.r3pct(f"UPS precursors L={L}: share abs(err) > 1 % (own rows)", ("acc_own_sp", f"UPS@{L}", "share_gt1"), hi=False, nd=1)
            for key in sorted(A[C].get("prot_err", {})):
                gname, s = key.split("|")
                if gname == "all" or (gname in ("low", "high") and s == "UPS"):
                    lab = {"all": "all pairs", "low": "pairs with lower level <= low-max", "high": "pairs with lower level > low-max"}[gname]
                    self.r3(f"protein log2FC error, {s}, {lab}: median abs(err)", ("prot_err", key, "med_abs"), hi=False, nd=3)
                    self.r3pct(f"protein log2FC error, {s}, {lab}: share abs(err) > 1 %", ("prot_err", key, "share_gt1"), hi=False, nd=2)
                    if self.val(C, "prot_err", key, "slope_fc_on_ex") is not None:
                        self.r3(f"compression slope (observed on expected log2FC), {s}, {lab}", ("prot_err", key, "slope_fc_on_ex"), nd=3)
                    if len(classes) >= 3 and s != "HUMAN":
                        self.r3(f"observed / expected log2FC (median), {s}", ("prot_err", key, "median_obs_over_exp"), nd=3)
            cvd = r.get("cv", {})
            for lvl in ("protein", "precursor"):
                x = cvd.get(lvl, {})
                if x.get("n_common", 0) < 40 or C not in x:
                    continue
                self.head(f"matched-intensity SD, {lvl} level: {x['n_common']:,} shared by {', '.join(cvd['ref_arms'])}; quartiles from mean rank")
                for q in ("all", "Q1", "Q2", "Q3", "Q4"):
                    self.row(f"replicate SD (log2, median) {lvl}, {q}", g(x, C, q, "med_sd"), g(x, self.MQ, q, "med_sd"), g(x, self.FP, q, "med_sd"), hi=False, nd=3)
                self.row(f"replicate CV % (median) {lvl}, all", g(x, C, "all", "med_cv_pct"), g(x, self.MQ, "all", "med_cv_pct"), g(x, self.FP, "all", "med_cv_pct"), hi=False, nd=1)
        if "de" in A[C]:
            fam = self.fam
            famtxt = {"inter": f"intersection family, {r['families']['intersection']:,} proteins tested by every shared arm; variance prior and BH on that set",
                      "own": "own family, each arm tests its own proteins"}[fam]
            self.head(f"4 DE (BH 5 %, per level pair, summed; {famtxt}; FDP = FP / (TP + FP) %)")
            for lab, arm in [("andes", C)] + [(v, v) for v in self.views]:
                for cut, tp, fp in (("no cut", "TP", "FP"), ("FC >= 0.5", "TP5", "FP5")):
                    self.row(f"DE TP, {cut}, {lab}", self.dv(arm, tp), self.dv(self.MQ, tp), self.dv(self.FP, tp))
                    self.row(f"DE FP, {cut}, {lab}", self.dv(arm, fp), self.dv(self.MQ, fp), self.dv(self.FP, fp), hi=False, tol=0.10)

                    def fdp(t):
                        x, y = self.dv(t, tp), self.dv(t, fp)
                        return None if x is None else 100 * y / max(x + y, 1)
                    self.row(f"DE FDP %, {cut}, {lab}", fdp(arm), fdp(self.MQ), fdp(self.FP), hi=False, nd=1, abs_tol=2.0, pct=False)
            if "UPS" in classes:
                lowl = set(r["levels"][:3])
                sfx = "" if fam == "own" else "_" + fam
                self.head(f"UPS pairs (FC >= 0.5, {fam} family; lower level among the 3 lowest levels vs the rest; summed over pairs)")
                for lab, arm in [("andes", C)] + [(v, v) for v in self.views]:
                    for grp, test in (("lowest-3-level pairs", lambda p: p["lo"] in lowl), ("other pairs", lambda p: p["lo"] not in lowl)):
                        s_ = lambda t, k: None if (t is None or t not in A) else sum(p[k + sfx] for p in A[t]["pairs"] if test(p))
                        self.row(f"{lab}: TP in {grp}", s_(arm, "TP5"), s_(self.MQ, "TP5"), s_(self.FP, "TP5"))
                        self.row(f"{lab}: FP in {grp}", s_(arm, "FP5"), s_(self.MQ, "FP5"), s_(self.FP, "FP5"), hi=False, tol=0.10)
        return self.rows

    def absence_rows(self):
        """Absence table rows: rows that should not exist (species absent from runs, UPS transfers below detection). Blank-run carry-over is
        not part of this score; see carry_over_rows."""
        A, C, r = self.A, self.C, self.r
        self.rows = []
        for arm in [C] + self.views:
            for s, x in (g(r.get("absence", {}), arm, "absent_species") or {}).items():
                gs = lambda t, k: g(r["absence"], t, "absent_species", s, k)
                self.row(f"{s} rows in the {x['absent_runs']} runs without {s} [{arm}]", x["rows"], gs(self.MQ, "rows"), gs(self.FP, "rows"), hi=False, tol=0.10)
                self.row(f"  of which transferred [{arm}]", x["transferred"], gs(self.MQ, "transferred"), gs(self.FP, "transferred"), hi=False, tol=0.10)
        if "ups_strat" in A[C]:
            lv = r["levels"][:3]
            low = lambda t, key: None if (t is None or t not in A) else sum(A[t]["ups_strat"][L][key] for L in lv)
            self.row(f"UPS transferred rows at the 3 lowest levels ({', '.join(lv)})", *[low(x, "ups_rows_transferred") for x in (C, self.MQ, self.FP)], hi=False, tol=0.10)
            self.row("  of which the precursor has no identified row at that level", *[low(x, "ups_transfer_rows_without_ident_at_level") for x in (C, self.MQ, self.FP)], hi=False, tol=0.10)
        return self.rows


    def carry_over_rows(self):
        """Blank-run rows (carry-over from the previous injection): reported apart from the absence score."""
        r, C = self.r, self.C
        self.rows = []
        for arm in [C] + self.views:
            b = g(r.get("absence", {}), arm, "blank")
            if b:
                self.row(f"blank-run rows ({b['blank_runs']} runs) [{arm}]", b["blank_rows"], g(r["absence"], self.MQ, "blank", "blank_rows"), g(r["absence"], self.FP, "blank", "blank_rows"), hi=False, tol=0.10)
                self.row(f"blank-run rows transferred [{arm}]", b["blank_transferred"], g(r["absence"], self.MQ, "blank", "blank_transferred"), g(r["absence"], self.FP, "blank", "blank_transferred"), hi=False, tol=0.10)
        return self.rows


def rates_table(res):
    """Truth-measured false-transfer rates (rates.py): foreign-species transfers in the runs that lack the species, over all transfers there."""
    R = res.get("rates")
    if not R:
        return []
    out = []
    for s in sorted({s for v in R.values() for s in v}):
        out += [f"Absent species {s} (runs without {s}); rate = {s} transferred rows / all transferred rows in those runs; ratio = {s} transferred rows in those runs / "
                f"{s} transferred rows in the runs that contain it; scaled = rate / {s} share of the identified precursors (extension, not a measured rate).", "",
                f"| arm | transferred rows in runs without {s} | of which {s} | rate % | ratio % | scaled rate % | {s} share of identified precursors % | identified {s} rows in those runs | transferred rows, all runs |",
                "|---|---|---|---|---|---|---|---|---|"]
        for n, v in R.items():
            o = v.get(s)
            if not o:
                continue
            f = lambda x, nd=2: "-" if x != x else f"{100 * x:.{nd}f}"
            out.append(f"| {n} | {o['transfers_in_absent_runs']:,} | {o['false_transfers']:,} | {f(o['rate'])} | {f(o['yeast_ratio'])} | {f(o['rate_scaled'])} | "
                       f"{f(o['species_share_identified'], 1)} | {o['identified_in_absent_runs']:,} | {o['transfers_all_runs']:,} |")
        out.append("")
    return out


def depth_table(res, curve_arm=None):
    """Transferred rows of an andes arm at the false-transfer rate each tool shows (conservative step on the threshold grid)."""
    from . import rates as RT
    R, C = res.get("rates"), res.get("depth_curves")
    if not R or not C:
        return []
    curve_arm = curve_arm or next(iter(C))
    out = []
    for s in sorted({s for v in R.values() for s in v}):
        out += [f"Depth at a matched false-transfer rate ({s}): transferred rows (all non-blank runs) of `{curve_arm}` at the threshold with the most rows among those whose rate is <= the tool's rate; "
                "thresholds on the grid of rates.py (no interpolation). `-` = the arm's rate is above the tool's at every threshold of the grid.", "",
                "| reference | its rate % | its transferred rows | " + curve_arm + " by transfer_q: threshold / rows | by quant_q: threshold / rows | rows ratio (transfer_q) | rows ratio (quant_q) |",
                "|---|---|---|---|---|---|---|"]
        for n, v in R.items():
            o = v.get(s)
            if not o or n == curve_arm or o["rate"] != o["rate"]:
                continue
            cells, ratios = [], []
            for key in ("transfer_q", "quant_q"):
                m = RT.matched_depth(C[curve_arm][key], s, o["rate"])
                cells.append("-" if m is None else f"{m[0]:g} / {m[1]:,}")
                ratios.append("-" if (m is None or not o["transfers_all_runs"]) else f"{m[1] / o['transfers_all_runs']:.2f}")
            out.append(f"| {n} | {100 * o['rate']:.2f} | {o['transfers_all_runs']:,} | {cells[0]} | {cells[1]} | {ratios[0]} | {ratios[1]} |")
        out.append("")
    return out


def de_summary(res):
    """Per-dataset DE table with the TP / FDP verdict of every andes arm against each tool."""
    A = res["arms"]
    if "de" not in A.get(res["primary"], {}):
        return []
    tools = [(n, k["kind"]) for n, k in res["arm_roles"].items() if k["kind"] != "andes"]
    out = ["| arm | TP / FP (FDP %) | TP / FP FC>=0.5 | " + " | ".join(f"vs {n} (no cut / FC)" for n, _ in tools) + " |",
           "|---" * (3 + len(tools)) + "|"]
    for n, k in res["arm_roles"].items():
        de = g(A[n], "de_family", "inter") or g(A[n], "de")
        if not de:
            continue
        vs = []
        for t, _ in tools:
            b = g(A.get(t), "de_family", "inter") or g(A.get(t), "de")
            vs.append("-" if (not b or t == n) else de_verdict(de["TP"], de["FP"], b["TP"], b["FP"]) + "/" + de_verdict(de["TP5"], de["FP5"], b["TP5"], b["FP5"]))
        out.append(f"| {n} | {de['TP']} / {de['FP']} ({100 * de['FP'] / max(de['TP'] + de['FP'], 1):.1f}) | {de['TP5']} / {de['FP5']} | " + " | ".join(vs) + " |")
    return out


def family_table(res):
    """DE counts (no cut / FC >= 0.5) of every arm in the own, intersection and union BH families."""
    A = res["arms"]
    if not A.get(res["primary"], {}).get("de_family"):
        return []
    F = res["families"]
    out = [f"Family sizes (proteins tested per pair, summed): intersection {F['intersection']:,} proteins and union {F['union']:,} proteins of arms {', '.join(F['arms'])}; "
           "own = " + ", ".join(f"{n} {v:,}" for n, v in F["own"].items()) + " proteins.", "",
           "| arm | own TP / FP (FDP %) | intersection TP / FP (FDP %) | union TP / FP (FDP %) | own FC>=0.5 | intersection FC>=0.5 | union FC>=0.5 |", "|---" * 7 + "|"]
    for n in res["arm_roles"]:
        f = g(A[n], "de_family")
        if not f:
            continue
        c = lambda k, a="TP", b="FP": f"{f[k][a]:,} / {f[k][b]:,} ({100 * f[k][b] / max(f[k][a] + f[k][b], 1):.1f})"
        c5 = lambda k: f"{f[k]['TP5']:,} / {f[k]['FP5']:,}"
        out.append(f"| {n} | {c('own')} | {c('inter')} | {c('union')} | {c5('own')} | {c5('inter')} | {c5('union')} |")
    return out


def fdp_table(res):
    """True positives at matched pooled FDP and per-pair AUC (changed vs null species), intersection family (headline) and own."""
    A, C = res["arms"], res["primary"]
    if not A.get(C, {}).get("de_family"):
        return []
    roles = res["arm_roles"]
    MQ = next((n for n, k in roles.items() if k["kind"] == "maxquant"), None)
    FP = next((n for n, k in roles.items() if k["kind"] == "fragpipe"), None)
    out = ["| arm | TP at FDP 1 % (inter / own) | 5 % | 10 % | mean per-pair AUC (inter / own) | lowest pair AUC (inter) |", "|---" * 6 + "|"]
    for n in roles:
        f = g(A[n], "de_family")
        if not f:
            continue
        tpc = lambda lv: f"{f['inter']['tp_at_fdp'][lv]:,} / {f['own']['tp_at_fdp'][lv]:,}"
        au = lambda k: "n/a" if f[k]["auc_mean"] is None else f"{f[k]['auc_mean']:.3f}"
        lo = "n/a" if f["inter"]["auc_min"] is None else f"{f['inter']['auc_min']:.3f}"
        out.append(f"| {n} | {tpc('1')} | {tpc('5')} | {tpc('10')} | {au('inter')} / {au('own')} | {lo} |")
    out += ["", "TP at matched FDP, intersection family: margin and verdict of andes (first andes arm) against each tool (equal within max(3 %, 2)):", "",
            "| level | andes | MaxQuant | FragPipe | margin vs MQ | margin vs FP | verdict |", "|---|---|---|---|---|---|---|"]
    for lv in ("1", "5", "10"):
        a, m, p = (g(A.get(x), "de_family", "inter", "tp_at_fdp", lv) for x in (C, MQ, FP))
        v = lambda t: verdict(a, t, True, 0.03, 2.0)
        out.append(f"| FDP {lv} % | {fmt(a)} | {fmt(m)} | {fmt(p)} | {margin(a, m)} | {margin(a, p)} | MQ {v(m)} / FP {v(p)} |")
    return out


def auc_pairs_table(res):
    """Per-pair AUC (intersection family) for pairs of adjacent levels, and the share of pairs below 0.9, 0.95."""
    A, C = res["arms"], res["primary"]
    if not A.get(C, {}).get("pairs") or "auc_inter" not in A[C]["pairs"][0]:
        return []
    arms = [n for n in res["arm_roles"] if A[n].get("pairs")]
    lv = res.get("levels") or []
    out = ["| pair (lower / higher level) | " + " | ".join(arms) + " |", "|---" * (len(arms) + 1) + "|"]
    ix = {l: i for i, l in enumerate(lv)}
    for k, p in enumerate(A[C]["pairs"]):
        if ix.get(p["hi"], -9) - ix.get(p["lo"], 0) == 1 or len(lv) <= 3:
            vals = []
            for n in arms:
                q = A[n]["pairs"][k]
                vals.append("n/a" if q.get("auc_inter") is None else f"{q['auc_inter']:.3f}")
            out.append(f"| {p['lo']} / {p['hi']} | " + " | ".join(vals) + " |")
    out.append("| **share of pairs with AUC < 0.9 / < 0.95** | " + " | ".join(
        "{} / {}".format(*[f"{100 * np.mean([q['auc_inter'] < t for q in A[n]['pairs'] if q.get('auc_inter') is not None]):.0f} %" for t in (0.9, 0.95)])
        for n in arms) + " |")
    return out


def rt_table(res):
    """DE with and without the RT-local normalisation, intersection family of each variant."""
    A = res["arms"]
    if not any("rtlocal" in a for a in A.values()):
        return []
    F = res["rtlocal_families"]
    out = [f"RT-local family: intersection {F['intersection']:,} proteins (arms {', '.join(F['arms'])}); RT proxy for {F['rt_keys']:,} precursor keys from the andes apex RT.", "",
           "| arm | normalisation | TP / FP (FDP %) | TP / FP FC>=0.5 | TP at FDP 5 % | TP at FDP 10 % | mean per-pair AUC |", "|---" * 7 + "|"]
    for n in res["arm_roles"]:
        pair = (("global median (default)", g(A[n], "de_family", "inter")), ("RT-local", g(A[n], "rtlocal", "de_family", "inter")))
        for lab, f in pair:
            if not f:
                continue
            au = "n/a" if f["auc_mean"] is None else f"{f['auc_mean']:.3f}"
            out.append(f"| {n} | {lab} | {f['TP']:,} / {f['FP']:,} ({100 * f['FP'] / max(f['TP'] + f['FP'], 1):.1f}) | {f['TP5']:,} / {f['FP5']:,} | "
                       f"{f['tp_at_fdp']['5']:,} | {f['tp_at_fdp']['10']:,} | {au} |")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", nargs="+", required=True, help="JSON files written by quant_eval.downstream")
    a = ap.parse_args(argv)
    out = ["# Quantification scorecard", ""]
    for path in a.results:
        r = json.load(open(path))
        arms = ", ".join("%s (%s)" % (n, k["kind"]) for n, k in r["arm_roles"].items())
        out += [f"## {r['dataset']} ({r['runs']} runs)", "",
                f"Arms: {arms}. Shared cells (arms {', '.join(r['shared_arms'])}): {r['shared_cells']:,}.", ""]
        hdr = ["| metric | andes | MaxQuant | FragPipe | margin vs MQ | margin vs FP | verdict |", "|---|---|---|---|---|---|---|"]
        out += hdr + Section(r).build() + [""]
        ds = de_summary(r)
        if ds:
            out += ["### DE summary (TP / FDP verdict rule; intersection family)", ""] + ds + [""]
        for title, tab in (("DE by BH family (own / intersection / union)", family_table(r)),
                           ("True positives at matched FDP and per-pair AUC", fdp_table(r)),
                           ("Per-pair AUC (changed vs null species, intersection family), adjacent levels", auc_pairs_table(r)),
                           ("DE with and without RT-local normalisation (same step for every arm)", rt_table(r))):
            if tab:
                out += [f"### {title}", ""] + tab + [""]
        ab = Section(r)
        rows = ab.absence_rows()
        if rows:
            out += ["### Absence (rows that should not exist; lower is better)", "", *hdr, *rows, ""]
        rt = rates_table(r)
        if rt:
            out += ["### False-transfer rates (truth: foreign species in runs that lack it)", ""] + rt
        dt = depth_table(r)
        if dt:
            out += ["### Depth at a matched false-transfer rate", ""] + dt
        co = Section(r).carry_over_rows()
        if co:
            out += ["### Blank-run carry-over (reported separately; not part of the absence score)", "", *hdr, *co, ""]
    print("\n".join(out))


if __name__ == "__main__":
    main()
