#!/usr/bin/env python3
"""Adoption gate: a candidate andes build against the baseline build, from the suite result JSONs that contain both builds.

  python -m quant_eval.gate --results results/ --candidate cand [--baseline main]

Criteria (pass / fail printed for each, exit status 1 when any fails):
  1 identified rows        identified rows not fewer; identified-cell error share, bright-decile gross share and top-decile replicate SD not worse
  2 wrong-peak             masked wrong-peak rate at quant_q <= 0.01 not worse in any cell (set x seed x pooled / level A / level B), tolerance --wp-tol
  3 UPS1 within 2x         not lower in any seed at any level (V rows of the UPS sets)
  4 DE (pairwise rule)       candidate vs baseline on (TP, FDP) over set x seed x {default, V, H} x {no cut, FC >= 0.5} x {all-species, null-species normalisation}, and the
                           candidate's matched-cut arm against the baseline's H arm: equal or better in >= 95 % of the cells and no worse cell on a UPS1 set
  5 matched depth          at the matched true-error cut (matched.py) the candidate keeps >= 97 % of the baseline H arm's precursors and transferred rows
  6 cost                   CPU and RSS (median over the seeds of a set) of the candidate's finalize + extract not above the baseline's + 3 %; wall is informational
                           (single timings on shared nodes vary up to 2x between identical runs: repeat a failing cost cell in a controlled run before acting on it)
Held-out datasets: results flagged held_out are refused (never part of a gate or tuning run) unless --held-out-eval is given."""
import argparse
import glob
import json
import os
import sys

import numpy as np

DE_VIEWS = (("", "de_family"), ("nullnorm", "nullnorm"))


def verd(a, t):
    """Pairwise rule. a, t = (TP, FDP %). better: TP >= +3 % with FDP not higher by > 2 points, or FDP lower by >= 2 points with TP >= -3 %; worse: the mirror."""
    tp = a[0] / max(t[0], 1) - 1
    d = a[1] - t[1]
    if (tp >= .03 and d <= 2) or (d <= -2 and tp >= -.03):
        return "better"
    if (tp <= -.03 and d >= -2) or (d >= 2 and tp <= .03):
        return "worse"
    return "equal"


def g(d, *k):
    for x in k:
        d = d.get(x) if isinstance(d, dict) else None
        if d is None:
            return None
    return d


def fdp(tp, fp):
    return 100.0 * fp / max(tp + fp, 1)


def de_cells(r, arm):
    """{(view, cut): (TP, FDP)} of one arm."""
    out = {}
    for view, key in DE_VIEWS:
        f = g(r["arms"], arm, "de_family", "inter") if view == "" else g(r["arms"], arm, "nullnorm", "de_family", "inter")
        if not f:
            continue
        out[(view or "all", "none")] = (f["TP"], fdp(f["TP"], f["FP"]))
        out[(view or "all", "FC")] = (f["TP5"], fdp(f["TP5"], f["FP5"]))
    return out


def run(R, cand, base, wp_tol=0.0002, de_share=0.95, ups_tol=0):
    res = {k: {"pass": True, "notes": []} for k in ("identified", "wrong_peak", "ups1_within2x", "de", "matched_depth", "cost")}

    def fail(k, msg):
        res[k]["pass"] = False
        res[k]["notes"].append(msg)

    de_total = de_ok = 0
    for dn, S in R.items():
        for sd, r in S.items():
            if cand not in r["builds"] or base not in r["builds"]:
                continue
            A = r["arms"]
            tag = f"{dn} s{sd}"
            # 1
            if A[cand]["ident_rows"] < A[base]["ident_rows"]:
                fail("identified", f"{tag}: identified rows {A[cand]['ident_rows']} < {A[base]['ident_rows']}")
            ac, ab = g(r, "cell_accuracy", "arms", cand), g(r, "cell_accuracy", "arms", base)
            if ac and ab:
                for k, tol in (("bad", 0.05), ("bad_bright", 0.05)):
                    if ac[k] > ab[k] + tol:
                        fail("identified", f"{tag}: {k} {ac[k]:.2f} > {ab[k]:.2f}")
                if ac.get("rep_sd_top") and ab.get("rep_sd_top") and ac["rep_sd_top"] > 1.02 * ab["rep_sd_top"]:
                    fail("identified", f"{tag}: top-decile replicate SD {ac['rep_sd_top']:.3f} > {ab['rep_sd_top']:.3f}")
            # 2
            wc, wb = g(r, "wrong_peak", cand), g(r, "wrong_peak", base)
            if wc and wb:
                pc = min(wc, key=lambda x: abs(x["tau"] - 0.01))
                pb = min(wb, key=lambda x: abs(x["tau"] - 0.01))
                for lab, a, b in (("pooled", pc, pb), ("A", pc["A"], pb["A"]), ("B", pc["B"], pb["B"])):
                    rc, rb = a["wrong"] / max(a["n"], 1), b["wrong"] / max(b["n"], 1)
                    if rc > rb + wp_tol:
                        fail("wrong_peak", f"{tag} {lab}: {100 * rc:.2f} % > {100 * rb:.2f} %")
            # 3
            uc, ub = g(A, cand + "_V", "ups_design"), g(A, base + "_V", "ups_design")
            if uc and ub:
                for L, o in ub.items():
                    if uc.get(L, {}).get("t1.0", 0) < o["t1.0"] - ups_tol:
                        fail("ups1_within2x", f"{tag} level {L}: {uc.get(L, {}).get('t1.0', 0)} < {o['t1.0']}")
            # 4
            is_ups = "UPS" in r.get("classes", [])
            pairs = [(cand, base), (cand + "_V", base + "_V"), (cand + "_H", base + "_H")]
            if cand + "_M" in A:
                pairs.append((cand + "_M", base + "_H"))
            for ca, ba in pairs:
                if ca not in A or ba not in A:
                    continue
                dc, db = de_cells(r, ca), de_cells(r, ba)
                for k in dc:
                    if k in db:
                        v = verd(dc[k], db[k])
                        de_total += 1
                        de_ok += v != "worse"
                        if v == "worse" and is_ups:
                            fail("de", f"{tag} {ca} vs {ba} {k}: worse ({dc[k][0]:,} / {dc[k][1]:.1f} vs {db[k][0]:,} / {db[k][1]:.1f})")
            # 5
            if cand + "_M" in A:
                m, h = A[cand + "_M"], A[base + "_H"]
                for k in ("precursors", "transfer_rows"):
                    if m[k] < 0.97 * h[k]:
                        fail("matched_depth", f"{tag}: {k} {m[k]:,} < 97 % of {h[k]:,}")
            elif g(r, "matched_cuts", cand, "tau") is None:
                if g(r, "wrong_peak", cand) or r.get("absence_curves"):
                    fail("matched_depth", f"{tag}: candidate error above the baseline's at every cut")
                else:
                    res["matched_depth"]["notes"].append(f"{tag}: no masked run and no absence curve, no matched cut")
    # 6: cost, median over the seeds of each set (one finalize timing varies up to 2x between identical runs on shared nodes)
    for dn, S in R.items():
        cs = [(g(r, "cost", cand), g(r, "cost", base)) for r in S.values() if cand in r["builds"] and base in r["builds"]]
        cs = [(c, b) for c, b in cs if c and b and c["files"] and b["files"]]
        if not cs:
            continue
        for k in ("cpu_s", "rss_gb"):
            mc, mb = float(np.median([c[k] for c, b in cs])), float(np.median([b[k] for c, b in cs]))
            if mc > 1.03 * mb:
                fail("cost", f"{dn}: median {k} over {len(cs)} seed(s) {mc:.1f} > 1.03 x {mb:.1f}")
        res["cost"]["notes"].append(f"{dn}: wall median {np.median([c['wall_s'] for c, b in cs]):.0f} s vs {np.median([b['wall_s'] for c, b in cs]):.0f} s (informational)")
    if de_total and de_ok / de_total < de_share:
        fail("de", f"equal-or-better in {100 * de_ok / de_total:.1f} % of {de_total} cells (< {100 * de_share:.0f} %)")
    res["de"]["notes"].append(f"{de_ok} of {de_total} cells equal or better")
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--baseline", default="main")
    ap.add_argument("--wp-tol", type=float, default=0.0002, help="wrong-peak rate tolerance (absolute), the seed noise of the r34C reproduction")
    ap.add_argument("--ups-tol", type=int, default=0, help="rows of slack in the UPS1 within-2x criterion (default 0: strict)")
    ap.add_argument("--held-out-eval", action="store_true")
    a = ap.parse_args(argv)
    R = {}
    for f in sorted(glob.glob(os.path.join(a.results, "*_s*.json"))):
        j = json.load(open(f))
        if "suite_dataset" not in j:
            continue
        if j.get("held_out") and not a.held_out_eval:
            raise SystemExit(f"{f}: held-out dataset refused in a gate (explicit --held-out-eval mode required)")
        R.setdefault(j["suite_dataset"], {})[j["seed"]] = j
    res = run(R, a.candidate, a.baseline, a.wp_tol, ups_tol=a.ups_tol)
    if a.held_out_eval:
        print("HELD-OUT EVALUATION (not a tuning or gate run)")
    bad = 0
    for k, v in res.items():
        print(f"{'PASS' if v['pass'] else 'FAIL'}  {k}")
        for n in v["notes"][:12]:
            print("      " + n)
        bad += not v["pass"]
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
