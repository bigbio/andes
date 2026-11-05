#!/usr/bin/env python3
"""ONE scoreboard table from the suite result JSONs (suite.py): per set and tool (arm).

  python -m quant_eval.board --results results/ [--seed 42] [--build main] > scoreboard.md

Columns (seed `--seed`; DE also as the range over all seeds present):
  prec / prot         precursors and proteins with >= 2 precursors; rows identified / transferred
  acc                 identified cells, every arm on the same cells: share |err| > 1 log2, within-2x share, bright-decile gross share (accuracy.py)
  prec.SD / CV        replicate SD of the top predicted-intensity decile (log2); protein replicate CV, median % (matched features)
  L1                  UPS1 precursors within 2x of their design reference, summed over the three lowest levels
  absence             PXD014415: yeast share of the transferred rows in the 20 human-only runs; single-cell set: rows in the Blank run
  wrong-peak          masked transfers off the identified apex at the arm's quant_q cut, pooled (level B) (andes hard-cut arms)
  DE                  intersection family, BH 5 % per pair: TP / FP without a fold-change cut, with |log2FC| >= 0.5, TP at pooled FDP 10 %
Arm names: `X` transfer_q <= 0.01, `X_V` all accepted transfers with learned precision weights, `X_H` hard cut quant_q <= 0.01, `X_I` identified rows only (no match between runs), `X_M` hard cut at the
matched true-error cut of the reference build (matched.py). Tools are context only."""
import argparse
import glob
import json
import os

NAME = {"MQ": "MaxQuant (MBR)", "MQ_noMBR": "MaxQuant no MBR", "FP": "FragPipe", "Sage": "Sage"}


def load(d):
    R = {}
    for f in sorted(glob.glob(os.path.join(d, "*_s*.json"))):
        j = json.load(open(f))
        if "suite_dataset" in j:
            R.setdefault(j["suite_dataset"], {})[j["seed"]] = j
    return R


def g(d, *k):
    for x in k:
        d = d.get(x) if isinstance(d, dict) else None
        if d is None:
            return None
    return d


def fmt(x, f="{:.2f}", none="-"):
    return none if x is None else f.format(x)


def de(A, n, key, seed_range=None):
    f = g(A, n, "de_family", "inter")
    if not f:
        return "-"
    return f"{f[key[0]]:,} / {f[key[1]]:,}"


def label(n, builds):
    for b in builds:
        if n == b:
            return f"andes {b}"
        for suf, nm in (("_V", "V"), ("_H", "H"), ("_I", "no MBR"), ("_M", "M (matched)")):
            if n == b + suf:
                return f"andes {b} {nm}"
    return NAME.get(n, n)


def build_rows(R, seed, builds, order):
    hdr = ("| set | arm | prec | prot | rows id / tr | err>1 % | within 2x % | bright gross % | prec.SD top | prot CV % | L1 | absence | wrong-peak % (B) | DE TP / FP | DE FC>=0.5 | TP@FDP10 | DE TP/FP range over seeds |\n"
           "|" + "---|" * 17)
    out = [hdr]
    for dn in order:
        if dn not in R:
            continue
        S = R[dn]
        sd = seed if seed in S else sorted(S)[0]
        r = S[sd]
        A = r["arms"]
        ca = r.get("cell_accuracy") or {}
        for n in A:
            a = A[n]
            plain = n in builds or r["arm_roles"][n]["kind"] != "andes"
            acc = (ca.get("arms") or {}).get(n) if plain else None
            ud = a.get("ups_design") or {}
            lv = r.get("levels") or []
            l1 = sum(g(ud, str(L), "t1.0") or 0 for L in lv[:3]) if ud else None
            ab = None
            if dn_is_p14(r):
                ab = g(r, "rates", n, "YEAST", "rate")
                ab = None if ab is None else f"{100 * ab:.2f} %"
            elif r.get("carry_over"):
                co = r["carry_over"].get(n)
                ab = None if not co else f"{co['rows']} rows"
            wp = "-"
            for b in builds:
                if n in (b + "_H", b + "_M"):
                    cu = g(r, "wrong_peak", b)
                    tau = 0.01 if n.endswith("_H") else g(r, "matched_cuts", b, "tau")
                    if cu and tau:
                        p = min(cu, key=lambda x: abs(x["tau"] - tau))
                        wp = f"{100 * p['wrong'] / max(p['n'], 1):.2f} ({100 * p['B']['wrong'] / max(p['B']['n'], 1):.2f})" + (f" @{tau:g}" if n.endswith("_M") else "")
            f = g(A, n, "de_family", "inter")
            rng = "-"
            if f:
                tps = [g(S[s], "arms", n, "de_family", "inter", "TP") for s in S]
                fps = [g(S[s], "arms", n, "de_family", "inter", "FP") for s in S]
                tps, fps = [x for x in tps if x is not None], [x for x in fps if x is not None]
                if len(tps) > 1:
                    rng = f"{min(tps):,}-{max(tps):,} / {min(fps):,}-{max(fps):,}"
            cvp = g(r, "cv", "protein", n, "all", "med_cv_pct")
            out.append("| " + " | ".join([
                r["dataset"] + (f" s{sd}" if False else ""), label(n, builds), f"{a['precursors']:,}", f"{a['proteins_ge2']:,}", f"{a['ident_rows']:,} / {a['transfer_rows']:,}",
                fmt(acc and acc["bad"]), fmt(acc and acc["within2x"]), fmt(acc and acc["bad_bright"]), fmt(acc and acc.get("rep_sd_top"), "{:.3f}"), fmt(cvp, "{:.1f}"),
                fmt(l1, "{:,}"), ab or "-", wp, de(A, n, ("TP", "FP")), de(A, n, ("TP5", "FP5")),
                (f"{f['tp_at_fdp']['10']:,}" if f else "-"), rng]) + " |")
    return "\n".join(out)


def dn_is_p14(r):
    return bool(r.get("rates"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--seed", default="42")
    ap.add_argument("--order", default="u9740,u9815,ups27,qe,astral,pxd030780,p14,sc")
    a = ap.parse_args(argv)
    R = load(a.results)
    builds = sorted({b for S in R.values() for r in S.values() for b in r["builds"]})
    order = a.order.split(",")
    print(build_rows(R, a.seed, builds, order + [d for d in sorted(R) if d not in order]))


if __name__ == "__main__":
    main()
