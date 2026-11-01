#!/usr/bin/env python3
"""Scorecard tables (markdown) from the result files of downstream.py.

  python -m quant_eval.scorecard --results DS1.json DS2.json ... > scorecard.md

One section per dataset: andes (the first andes arm) against MaxQuant and FragPipe, with the margin (andes minus tool; absolute
and % of the tool) and a verdict per tool.

Verdict rule for counts and errors: `equal` when |margin| <= max(3 % of the tool value, 2 units for integer counts; the
tolerance is 10 % for FP and absence counts); otherwise
`better` / `worse` by the metric direction (lower is better for FP, FDP, error and SD). FDP rows use an absolute band of 2
points. The DE summary table uses the TP / FDP rule: better = TP >= 0.98 x tool TP and FDP lower by > 2 points, or
TP > 1.03 x tool TP and FDP not higher by > 2 points; worse = the mirror; else equal."""
import argparse
import json


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
            self.head("4 DE (BH 5 %, per level pair, summed; FDP = FP / (TP + FP) %)")
            for lab, arm in [("andes", C)] + [(v, v) for v in self.views]:
                for cut, tp, fp in (("no cut", "TP", "FP"), ("FC >= 0.5", "TP5", "FP5")):
                    a_tp, a_fp = g(A[arm], "de", tp), g(A[arm], "de", fp)
                    self.row(f"DE TP, {cut}, {lab}", a_tp, self.val(self.MQ, "de", tp), self.val(self.FP, "de", tp))
                    self.row(f"DE FP, {cut}, {lab}", a_fp, self.val(self.MQ, "de", fp), self.val(self.FP, "de", fp), hi=False, tol=0.10)

                    def fdp(t):
                        if t is None or t not in A:
                            return None
                        x, y = g(A[t], "de", tp), g(A[t], "de", fp)
                        return 100 * y / max(x + y, 1)
                    self.row(f"DE FDP %, {cut}, {lab}", fdp(arm), fdp(self.MQ), fdp(self.FP), hi=False, nd=1, abs_tol=2.0, pct=False)
            if "UPS" in classes:
                lowl = set(r["levels"][:3])
                self.head("UPS pairs (FC >= 0.5; lower level among the 3 lowest levels vs the rest; summed over pairs)")
                for lab, arm in [("andes", C)] + [(v, v) for v in self.views]:
                    for grp, test in (("lowest-3-level pairs", lambda p: p["lo"] in lowl), ("other pairs", lambda p: p["lo"] not in lowl)):
                        s_ = lambda t, k: None if (t is None or t not in A) else sum(p[k] for p in A[t]["pairs"] if test(p))
                        self.row(f"{lab}: TP in {grp}", s_(arm, "TP5"), s_(self.MQ, "TP5"), s_(self.FP, "TP5"))
                        self.row(f"{lab}: FP in {grp}", s_(arm, "FP5"), s_(self.MQ, "FP5"), s_(self.FP, "FP5"), hi=False, tol=0.10)
        self.head("5 Absence")
        for arm in [C] + self.views:
            b = g(r.get("absence", {}), arm, "blank")
            if b:
                self.row(f"blank-run rows ({b['blank_runs']} runs) [{arm}]", b["blank_rows"], g(r["absence"], self.MQ, "blank", "blank_rows"), g(r["absence"], self.FP, "blank", "blank_rows"), hi=False, tol=0.10)
                self.row(f"blank-run rows transferred [{arm}]", b["blank_transferred"], g(r["absence"], self.MQ, "blank", "blank_transferred"), g(r["absence"], self.FP, "blank", "blank_transferred"), hi=False, tol=0.10)
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


def de_summary(res):
    """Per-dataset DE table with the TP / FDP verdict of every andes arm against each tool."""
    A = res["arms"]
    if "de" not in A.get(res["primary"], {}):
        return []
    tools = [(n, k["kind"]) for n, k in res["arm_roles"].items() if k["kind"] != "andes"]
    out = ["| arm | TP / FP (FDP %) | TP / FP FC>=0.5 | " + " | ".join(f"vs {n} (no cut / FC)" for n, _ in tools) + " |",
           "|---" * (3 + len(tools)) + "|"]
    for n, k in res["arm_roles"].items():
        de = g(A[n], "de")
        if not de:
            continue
        vs = []
        for t, _ in tools:
            b = g(A.get(t), "de")
            vs.append("-" if (not b or t == n) else de_verdict(de["TP"], de["FP"], b["TP"], b["FP"]) + "/" + de_verdict(de["TP5"], de["FP5"], b["TP5"], b["FP5"]))
        out.append(f"| {n} | {de['TP']} / {de['FP']} ({100 * de['FP'] / max(de['TP'] + de['FP'], 1):.1f}) | {de['TP5']} / {de['FP5']} | " + " | ".join(vs) + " |")
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
            out += ["### DE summary (TP / FDP verdict rule)", ""] + ds + [""]
    print("\n".join(out))


if __name__ == "__main__":
    main()
