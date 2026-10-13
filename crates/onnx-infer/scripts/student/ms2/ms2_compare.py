#!/usr/bin/env python3
"""Spectral-angle comparison of fragment-intensity predictors on confident PSMs.

Each ``NAME=PATH`` is a feature TSV written by ``examples/rescore_pin.rs`` for
the same PIN (one row per PIN row, keyed by ``SpecId``), or the output of
``model-train/examples/frag_model_baseline.rs`` for the andes GBDT. Only PSMs listed in
--psms (a TSV with ``SpecId`` and ``peptide``, e.g. Percolator q <= 0.01 targets)
and predicted by every model are compared. Prints the median spectral angle
(IQR), median Pearson and the fraction with SA > 0.7 per model, then per
precursor charge (from the ``charge`` column of --psms when present) and length.

Usage::

  python ms2_compare.py --psms conf.tsv pdeep=pdeep.tsv ams2=ams2.tsv ms2pip=ms2pip.tsv
"""
from __future__ import annotations

import argparse
import csv

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--psms", required=True)
    p.add_argument("features", nargs="+", help="NAME=rescore_pin feature TSV")
    a = p.parse_args()
    conf = pd.read_csv(a.psms, sep="\t", quoting=csv.QUOTE_NONE)
    d = conf.copy()
    names = []
    for item in a.features:
        name, path = item.split("=", 1)
        f = pd.read_csv(path, sep="\t", quoting=csv.QUOTE_NONE, na_values="NA")
        if "spectral_angle" not in f.columns and "gbdt_sa" in f.columns:
            # model-train's frag_model_baseline output (andes frag GBDT)
            f = f.rename(columns={"gbdt_sa": "spectral_angle", "gbdt_pearson": "spec_pearson_norm"})
            if "charge" not in d.columns and "charge" in f.columns:
                d = d.merge(f[["SpecId", "charge"]], on="SpecId", how="left")
        else:
            miss = [c for c in f.columns if c.endswith("_missing")][0]
            f = f[f[miss] == 0]
        f = f[["SpecId", "spectral_angle", "spec_pearson_norm"]]
        f.columns = ["SpecId", f"{name}_sa", f"{name}_pearson"]
        d = d.merge(f, on="SpecId")
        names.append(name)
    print(f"confident PSMs {len(conf)}, predicted by every model {len(d)}")
    for n in names:
        sa = d[f"{n}_sa"]
        print(f"{n:10s} SA median {sa.median():.3f} (IQR {sa.quantile(.25):.3f}-{sa.quantile(.75):.3f})  "
              f"pearson {d[f'{n}_pearson'].median():.3f}  frac SA>0.7 {np.mean(sa > 0.7):.3f}")
    if "charge" in d.columns:
        for z, g in d.groupby("charge"):
            if len(g) >= 50:
                print(f"  z={z} n={len(g):6d} " + "  ".join(f"{n} {g[f'{n}_sa'].median():.3f}" for n in names))
    if "peptide" in d.columns:
        L = d.peptide.str.replace(r"^[^.]*\.|\.[^.]*$", "", regex=True).str.replace(r"[^A-Z]", "", regex=True).str.len()
        for lo, hi in [(6, 10), (11, 15), (16, 20), (21, 60)]:
            g = d[(L >= lo) & (L <= hi)]
            if len(g):
                print(f"  len {lo}-{hi} n={len(g):6d} " + "  ".join(f"{n} {g[f'{n}_sa'].median():.3f}" for n in names))


if __name__ == "__main__":
    main()
