#!/usr/bin/env python3
"""RT accuracy on a benchmark run and RT-only rescoring PINs, per predictor.

Inputs:
  --pin            andes PIN
  --rt-tsv         ``rt_pin`` example output (SpecId, observed_retention_time,
                   andes_rt, rt_missing, proforma), one row per PIN row
  --baseline-psms  Percolator PSMs of the base PIN (confident set: q <= 0.01)
  --peptdeep DIR   AlphaPeptDeep generic pretrained dir (rt.pth)
  --andes-index    optional ``frag_model_baseline`` output (rt_seed_index, rt_features)

Accuracy: confident target PSMs, peptides split in half by sequence;
isotonic calibration on one half, error on the other.

Arms (``--out-prefix``_<arm>.pin): ``rt_diff`` = |observed - calibrated| and
the calibrated RT, calibration fitted on all confident targets and applied to
every row; rows without a prediction get the median and ``rt_missing`` = 1.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from rt_calibration import fit_isotonic  # noqa: E402
from rt_eval import peptdeep_rt, ptm_class  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pin", required=True)
    p.add_argument("--rt-tsv", required=True)
    p.add_argument("--baseline-psms", required=True)
    p.add_argument("--peptdeep", required=True)
    p.add_argument("--andes-index")
    p.add_argument("--extra-pred", action="append", default=[],
                   help="NAME=TSV:COLUMN[:MISSING_COLUMN], a per-PIN-row predictor keyed by SpecId (e.g. DeepLC)")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--no-arms", action="store_true")
    a = p.parse_args()

    with open(a.pin) as f:
        header = f.readline().rstrip("\n").split("\t")
        lines = [ln.rstrip("\n").split("\t") for ln in f]
    k = header.index("Peptide")
    labels = np.array([int(x[1]) for x in lines])
    rt = pd.read_csv(a.rt_tsv, sep="\t", na_values="NA", quoting=csv.QUOTE_NONE, keep_default_na=False)
    assert len(rt) == len(lines) and (rt.SpecId.values == np.array([x[0] for x in lines])).all(), "row order"
    forms = rt.proforma[rt.rt_missing == 0].unique()
    pd_rt = dict(zip(forms, peptdeep_rt(list(forms), a.peptdeep)))
    rt["pdeep_rt"] = rt.proforma.map(pd_rt).astype(float)
    obs = rt.observed_retention_time.values
    preds = {"andes-rt": rt.andes_rt.values.astype(float), "PeptDeep": rt.pdeep_rt.values}
    for spec in a.extra_pred:
        name, rest = spec.split("=", 1)
        path, col, *mcol = rest.split(":")
        x = pd.read_csv(path, sep="\t", na_values="NA", quoting=csv.QUOTE_NONE, usecols=["SpecId", col] + mcol)
        assert len(x) == len(rt) and (x.SpecId.values == rt.SpecId.values).all(), f"{name}: row order"
        v = x[col].values.astype(float)
        if mcol:
            v[x[mcol[0]].values == 1] = np.nan
        preds[name] = v
        rt[f"extra_{name}"] = v

    conf = set()
    with open(a.baseline_psms) as f:
        h = f.readline().rstrip("\n").split("\t")
        qi = h.index("q-value")
        for ln in f:
            x = ln.split("\t")
            if float(x[qi]) <= 0.01:
                conf.add(x[0])
    cm = rt.SpecId.isin(conf).values & (labels == 1) & np.isfinite(obs)

    # ---- accuracy on confident targets, held-out half of the peptides ----
    c = rt[cm].copy()
    c["seq"] = c.proforma.fillna("").str.replace(r"\[[^\]]*\]-?|-", "", regex=True)
    if a.andes_index:
        ai = pd.read_csv(a.andes_index, sep="\t", quoting=csv.QUOTE_NONE)
        c = c.merge(ai[["SpecId", "rt_seed_index", "rt_features"]], on="SpecId", how="left")
    ok = np.isfinite(c.andes_rt.values.astype(float)) & np.isfinite(c.pdeep_rt.values)
    for name in preds:
        if name not in ("andes-rt", "PeptDeep"):
            ok &= np.isfinite(c[f"extra_{name}"].values)
    if a.andes_index:
        ok &= c.rt_seed_index.notna().values
    c = c[ok].reset_index(drop=True)
    u = c.seq.unique()
    np.random.default_rng(0).shuffle(u)
    tr = c.seq.isin(set(u[: len(u) // 2])).values
    y = c.observed_retention_time.values
    cols = {"andes-rt": c.andes_rt.values.astype(float), "PeptDeep": c.pdeep_rt.values}
    for name in preds:
        if name not in cols:
            cols[name] = c[f"extra_{name}"].values
    if a.andes_index:
        cols["andes RT index (shipped)"] = c.rt_seed_index.values
        F = np.array([np.array(s.split(","), float) for s in c.rt_features])
        X = np.c_[F, np.ones(len(F))]
        coef, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        cols["andes RT index (refit in-run)"] = X @ coef
    c["ptm"] = c.proforma.map(ptm_class)
    span = np.quantile(y, 0.99) - np.quantile(y, 0.01)
    print(f"accuracy: {len(c)} confident target PSMs, {len(u)} peptides, held-out PSMs {(~tr).sum()}, "
          f"RT span {span:.1f} min")
    for name, v in cols.items():
        f = fit_isotonic(v[tr], y[tr])
        e = np.abs(f(v[~tr]) - y[~tr])
        r2 = 1 - np.sum((f(v[~tr]) - y[~tr]) ** 2) / np.sum((y[~tr] - y[~tr].mean()) ** 2)
        by = "; ".join(f"{k} {np.mean(e[c.ptm.values[~tr] == k]):.2f} (n={np.sum(c.ptm.values[~tr] == k)})"
                       for k in sorted(set(c.ptm)))
        print(f"  {name:30s} MAE {e.mean():.2f} MedAE {np.median(e):.2f} p90 {np.quantile(e, .9):.2f} min "
              f"R2 {r2:.4f} | MAE by PTM: {by}")
    if a.no_arms:
        return

    # ---- RT-only arms ----
    for name, v in preds.items():
        good = np.isfinite(v) & np.isfinite(obs)
        f = fit_isotonic(v[cm & good], obs[cm & good])
        cal = f(np.where(np.isfinite(v), v, 0))
        diff = np.abs(obs - cal)
        miss = (~good).astype(int)
        cal[miss == 1] = np.median(cal[good])
        diff[miss == 1] = np.median(diff[good])
        dec = labels == -1
        print(f"arm {name}: calibrated on {np.sum(cm & good)} confident targets; missing {miss.sum()}; "
              f"median rt_diff confident {np.median(diff[cm & good]):.2f} decoys {np.median(diff[dec & good]):.2f} min")
        tag = name.replace(" ", "_")
        with open(f"{a.out_prefix}_{tag}.pin", "w") as out:
            out.write("\t".join(header[:k] + ["rt_diff", "predicted_retention_time", "rt_missing"] + header[k:]) + "\n")
            for x, d_, c_, m_ in zip(lines, diff, cal, miss):
                out.write("\t".join(x[:k] + [f"{d_:.6g}", f"{c_:.6g}", str(m_)] + x[k:]) + "\n")
        print("wrote", f"{a.out_prefix}_{tag}.pin")


if __name__ == "__main__":
    main()
