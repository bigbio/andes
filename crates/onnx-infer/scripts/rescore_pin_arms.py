#!/usr/bin/env python3
"""Build A/B Percolator PINs from `rescore_pin` feature TSVs of one or more models.

``--pin`` and the feature TSVs may be gzipped. Each ``--features SRC=PATH`` is a TSV written by the ``rescore_pin`` example
(one row per PIN row, same order). Its columns ``<p>_irt`` / ``<p>_missing``
carry the model's raw RT prediction and missing flag; the spectral-match
features follow (NA without an MS2 model).

Per source:

* RT: the raw prediction is calibrated to observed RT on confident target
  PSMs (q <= --calib-q in the baseline Percolator run) with a binned-median
  curve; ``<src>_rt_diff = |observed - calibrated|`` and
  ``<src>_predicted_rt`` are written for every row.
* A source is RT-calibrated only when some ``--arm`` requests its ``rt``
  block, and its MS2 block is built only when some arm requests ``ms2``.
* Rows a model cannot encode keep their row: every feature of that source is
  imputed with its median over predicted rows, and ``<src>_missing`` = 1.

``--arm NAME=SRC:KIND[,SRC:KIND...]`` (KIND = ms2 | rt) writes
``<out-prefix>_<NAME>.pin`` with those feature blocks inserted before the
Peptide column. Feature columns are prefixed with the source name.

Usage::

  python rescore_pin_arms.py --pin search.pin --baseline-psms base_seed42.psms \\
      --features ms2pip=ms2pip.tsv deeplc=deeplc.tsv --out-prefix run \\
      --arm ms2pip=ms2pip:ms2 deeplc=deeplc:rt both=ms2pip:ms2,deeplc:rt
"""
from __future__ import annotations

import argparse
import csv
import gzip

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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pin", required=True)
    p.add_argument("--features", nargs="+", required=True, help="SRC=PATH")
    p.add_argument("--baseline-psms", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--arm", nargs="+", required=True, help="NAME=SRC:KIND[,SRC:KIND...]")
    p.add_argument("--calib-q", type=float, default=0.01)
    p.add_argument("--bins", type=int, default=50)
    p.add_argument("--gzip", action="store_true", help="write <out-prefix>_<NAME>.pin.gz")
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


def col_with_suffix(df, suffix):
    cols = [c for c in df.columns if c.endswith(suffix)]
    assert len(cols) == 1, f"expected one *{suffix} column, got {cols}"
    return cols[0]


KINDS = ("ms2", "rt")


def parse_arms(specs):
    """``NAME=SRC:KIND[,SRC:KIND...]`` -> [(name, [(src, kind), ...])]."""
    arms = []
    for spec in specs:
        name, parts = spec.split("=", 1)
        pairs = []
        for part in parts.split(","):
            src, kind = part.split(":")
            if kind not in KINDS:
                raise SystemExit(f"--arm {spec}: unknown kind {kind!r} (expected one of {KINDS})")
            pairs.append((src, kind))
        arms.append((name, pairs))
    return arms


def blocks_for(src, path, ids, labels, conf, a, kinds):
    """Return the requested ``kinds`` ({"ms2", "rt"}) of imputed, prefixed
    features, plus the "missing" block."""
    feat = pd.read_csv(path, sep="\t", na_values="NA", quoting=csv.QUOTE_NONE)
    assert len(feat) == len(ids) and (feat.SpecId.values == ids).all(), f"{src}: row order"
    miss = feat[col_with_suffix(feat, "_missing")].values.astype(int)
    irt = feat[col_with_suffix(feat, "_irt")].values.astype(float)
    obs = feat.observed_retention_time.values.astype(float)
    is_conf = feat.SpecId.isin(conf).values
    out = {}
    print(f"[{src}] missing rows {miss.sum()} of {len(miss)}")

    def impute(v):
        v = v.astype(float).copy()
        bad = ~np.isfinite(v) | (miss == 1)
        v[bad] = np.median(v[~bad])
        return v

    if "rt" in kinds and np.isfinite(irt).any():
        ok = (miss == 0) & np.isfinite(irt) & np.isfinite(obs)
        cm = is_conf & ok & (labels == 1)
        cal = calibrate(irt[cm], obs[cm], a.bins)
        pred_rt = cal(np.where(np.isfinite(irt), irt, np.nanmedian(irt)))
        diff = np.abs(obs - pred_rt)
        res = diff[cm]
        print(f"[{src}] RT calibration on {cm.sum()} confident targets: median |resid| "
              f"{np.median(res):.2f} min, p90 {np.quantile(res, 0.9):.2f} min; decoys median "
              f"rt_diff {np.nanmedian(diff[(labels == -1) & ok]):.2f} min")
        out["rt"] = pd.DataFrame({f"{src}_rt_diff": impute(diff), f"{src}_predicted_rt": impute(pred_rt)})
    if "ms2" in kinds and feat[MS2].notna().any().all():
        ok = miss == 0
        dec = labels == -1
        for c in ["spectral_angle", "cos_norm", "spec_pearson_norm", "pred_explained"]:
            v = feat[c].values
            print(f"[{src}] {c:18s} median: confident targets {np.nanmedian(v[is_conf & ok]):.3f} | "
                  f"all targets {np.nanmedian(v[(labels == 1) & ok]):.3f} | decoys {np.nanmedian(v[dec & ok]):.3f}")
        out["ms2"] = pd.DataFrame({f"{src}_{c}": impute(feat[c].values) for c in MS2})
    out["missing"] = pd.DataFrame({f"{src}_missing": miss})
    return out


def main():
    a = parse_args()
    with (gzip.open(a.pin, "rt") if a.pin.endswith(".gz") else open(a.pin)) as f:
        header = f.readline().rstrip("\n").split("\t")
        lines = [ln.rstrip("\n").split("\t") for ln in f]
    k = header.index("Peptide")
    ids = np.array([x[0] for x in lines])
    labels = np.array([int(x[1]) for x in lines])

    conf = set()
    with open(a.baseline_psms) as f:
        h = f.readline().rstrip("\n").split("\t")
        qi = h.index("q-value")
        for ln in f:
            x = ln.split("\t")
            if float(x[qi]) <= a.calib_q:
                conf.add(x[0])

    arms = parse_arms(a.arm)
    requested = {}
    for _, pairs in arms:
        for src, kind in pairs:
            requested.setdefault(src, set()).add(kind)

    blocks = {}
    for spec in a.features:
        src, path = spec.split("=", 1)
        blocks[src] = blocks_for(src, path, ids, labels, conf, a, requested.get(src, set()))

    for name, pairs in arms:
        frames, srcs = [], []
        for src, kind in pairs:
            frames.append(blocks[src][kind])
            if src not in srcs:
                srcs.append(src)
        frames += [blocks[s]["missing"] for s in srcs]
        df = pd.concat(frames, axis=1)
        cols = list(df.columns)
        vals = df.values
        path = f"{a.out_prefix}_{name}.pin" + (".gz" if a.gzip else "")
        with (gzip.open(path, "wt", compresslevel=3) if a.gzip else open(path, "w")) as out:
            out.write("\t".join(header[:k] + cols + header[k:]) + "\n")
            for x, row in zip(lines, vals):
                out.write("\t".join(x[:k] + [f"{v:.6g}" for v in row] + x[k:]) + "\n")
        print("wrote", path, "with", len(cols), "columns:", " ".join(cols))


if __name__ == "__main__":
    main()
