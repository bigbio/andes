#!/usr/bin/env python3
"""Markdown tables from the rollup_sens.py result JSONs.

  python -m quant_eval.rollup_sens_table --results res/ [--dataset u9740] > tables.md

One table per dataset: rows = roll-up x arm; the andes arms are the mean over the seeds present, with the range (max - min) in
brackets. Columns: proteins quantified (own / inside the matched-depth intersection family), protein log2-ratio error against the
design truth per species on the intersection family (mean |err|, bias), median CV %, DE true / false positives at BH 5 % (intersection
family) and true positives at pooled FDP 10 %."""
import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np


def load(d):
    R = defaultdict(dict)
    for f in sorted(glob.glob(os.path.join(d, "*_s*.json"))):
        j = json.load(open(f))
        R[j["dataset"]][j["seed"]] = j
    return R


def cell(vals, f="{:.2f}"):
    v = [x for x in vals if x is not None]
    if not v:
        return "-"
    if len(v) == 1:
        return f.format(v[0])
    return f.format(np.mean(v)) + " (" + f.format(max(v) - min(v)) + ")"


def species_of(R):
    for j in R.values():
        for ru in j["roll_ups"].values():
            if "arms" in ru:
                a = next(iter(ru["arms"].values()))
                return list(a["err_inter"])
    return []


def table(ds, runs):
    seeds = sorted(runs)
    first = runs[seeds[0]]
    sps = species_of(runs)
    lines = [f"### {first.get('label', ds)} ({first['runs']} runs; seeds {','.join(seeds)}; andes cells mean (range over seeds))", "",
             "| roll-up | arm | prot own | prot shared | " + " | ".join(f"{s} |err| / bias" for s in sps) + " | CV % | TP / FP (shared) | TP @10% FDP |",
             "|---|---|---|---|" + "---|" * len(sps) + "---|---|---|"]
    for ru in first["roll_ups"]:
        if "arms" not in first["roll_ups"][ru]:
            lines.append(f"| {ru} | - | failed: {first['roll_ups'][ru].get('error', '')[:60]} |" + " |" * (len(sps) + 4))
            continue
        for arm in first["roll_ups"][ru]["arms"]:
            per = [runs[s]["roll_ups"][ru]["arms"][arm] for s in seeds if ru in runs[s]["roll_ups"] and arm in runs[s]["roll_ups"][ru].get("arms", {})]
            if not per:
                continue
            row = [ru, arm, cell([p["proteins"] for p in per], "{:.0f}"), cell([p["proteins_inter"] for p in per], "{:.0f}")]
            for s in sps:
                e = [p["err_inter"].get(s) for p in per]
                row.append(cell([x["mean_abs"] for x in e if x], "{:.2f}") + " / " + cell([x["bias"] for x in e if x], "{:+.2f}"))
            row.append(cell([p["cv_inter"]["med_cv_pct"] for p in per], "{:.1f}"))
            row.append(cell([p["de_inter"]["TP"] for p in per], "{:.0f}") + " / " + cell([p["de_inter"]["FP"] for p in per], "{:.0f}"))
            row.append(cell([p["tp_at_fdp_inter"].get("10") for p in per], "{:.0f}"))
            lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--dataset")
    a = ap.parse_args()
    R = load(a.results)
    for ds, runs in R.items():
        if a.dataset and ds != a.dataset:
            continue
        print(table(ds, runs), "\n")


if __name__ == "__main__":
    main()
