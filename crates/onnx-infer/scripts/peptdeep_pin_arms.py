#!/usr/bin/env python3
"""Append PeptDeep features (from the `rescore_pin` example) to a PIN.

RT: predicted PeptDeep RT is calibrated to observed RT on confident target PSMs
(q <= --calib-q in a baseline Percolator run) with a binned-median curve, then
``rt_diff = |observed - calibrated|`` is computed identically for every row.

Rows without predictions keep their row: each feature is imputed with its
median over all predicted rows and ``pdeep_missing`` = 1 marks them.

Writes one PIN per arm (``ms2``, ``rt``, ``both``) plus agreement statistics.

Usage::

  python peptdeep_pin_arms.py --pin search.pin --features features.tsv \
      --baseline-psms baseline_seed42.psms --out-prefix run
"""
from __future__ import annotations

import argparse
import csv

import numpy as np
import pandas as pd

MS2 = [
    "spec_pearson_norm",
    "ionb_pearson_norm",
    "iony_pearson_norm",
    "spec_mse_norm",
    "dotprod_norm",
    "cos_norm",
    "spectral_angle",
    "pred_explained",
]
RT = ["rt_diff", "predicted_retention_time"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pin", required=True)
    p.add_argument("--features", required=True)
    p.add_argument("--baseline-psms", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--calib-q", type=float, default=0.01)
    p.add_argument("--bins", type=int, default=50)
    return p.parse_args()


def calibrate(pred, obs, bins):
    """Binned-median map predicted -> observed RT, linear between bin centres."""
    qs = np.quantile(pred, np.linspace(0, 1, bins + 1))
    xs, ys = [], []
    for lo, hi in zip(qs[:-1], qs[1:]):
        m = (pred >= lo) & (pred <= hi)
        if m.sum() >= 5:
            xs.append(np.median(pred[m]))
            ys.append(np.median(obs[m]))
    xs, ys = np.array(xs), np.maximum.accumulate(np.array(ys))
    return lambda x: np.interp(x, xs, ys)


def main():
    a = parse_args()
    with open(a.pin) as f:
        header = f.readline().rstrip("\n").split("\t")
        lines = [ln.rstrip("\n").split("\t") for ln in f]
    k = header.index("Peptide")
    ids = [x[0] for x in lines]
    labels = np.array([int(x[1]) for x in lines])
    feat = pd.read_csv(a.features, sep="\t", na_values="NA", quoting=csv.QUOTE_NONE)
    assert len(feat) == len(lines) and (feat.SpecId.values == np.array(ids)).all(), "row order"

    conf = set()
    with open(a.baseline_psms) as f:
        h = f.readline().rstrip("\n").split("\t")
        qi = h.index("q-value")
        for ln in f:
            x = ln.split("\t")
            if float(x[qi]) <= a.calib_q:
                conf.add(x[0])
    is_conf = feat.SpecId.isin(conf).values
    ok = (feat.pdeep_missing.values == 0) & np.isfinite(feat.observed_retention_time.values)
    cm = is_conf & ok & (labels == 1)
    cal = calibrate(feat.pdeep_irt.values[cm], feat.observed_retention_time.values[cm], a.bins)
    pred_rt = cal(feat.pdeep_irt.values)
    feat["predicted_retention_time"] = pred_rt
    feat["rt_diff"] = np.abs(feat.observed_retention_time.values - pred_rt)
    res = feat.rt_diff.values[cm]
    print(f"RT calibration on {cm.sum()} confident targets: median |resid| {np.median(res):.2f} min, "
          f"p90 {np.quantile(res, 0.9):.2f} min; decoys median rt_diff "
          f"{np.nanmedian(feat.rt_diff.values[(labels == -1) & ok]):.2f} min")

    for c in MS2 + RT:
        v = feat[c].values.astype(float)
        bad = ~np.isfinite(v) | (feat.pdeep_missing.values == 1)
        v[bad] = np.median(v[~bad])
        feat[c] = v
    miss = feat.pdeep_missing.values.astype(int)
    print(f"missing rows {miss.sum()} of {len(miss)}")

    # agreement stats
    dec = labels == -1
    for c in ["spectral_angle", "cos_norm", "spec_pearson_norm", "pred_explained", "rt_diff"]:
        v = feat[c].values
        print(f"{c:18s} median: confident targets {np.median(v[is_conf & ok]):.3f} | "
              f"all targets {np.median(v[(labels == 1) & ok]):.3f} | decoys {np.median(v[dec & ok]):.3f}")

    arms = {"ms2": MS2, "rt": RT, "both": MS2 + RT}
    for arm, cols in arms.items():
        cols = cols + ["pdeep_missing"]
        vals = feat[cols].values
        with open(f"{a.out_prefix}_{arm}.pin", "w") as out:
            out.write("\t".join(header[:k] + cols + header[k:]) + "\n")
            for x, row in zip(lines, vals):
                out.write("\t".join(x[:k] + [f"{v:.6g}" for v in row] + x[k:]) + "\n")
        print("wrote", f"{a.out_prefix}_{arm}.pin")


if __name__ == "__main__":
    main()
