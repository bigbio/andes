#!/usr/bin/env python3
"""Compare OpenMS MS1LabeledWorkflow ratios (QPX feature table) with andes.

    compare_ms1labeled.py QPX_FEATURE_PARQUET LABEL_PARQUET RUN=DESIGN_LOG2 ...

One log2 heavy/light per (run, sequence, charge) on each side (the median over
its features or channel rows); medians, IQRs and the median absolute error
against the design, for all values and for the peptides both quantify.
"""
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

feat, lab, *designs = sys.argv[1:]
f = pq.read_table(feat).to_pandas()
f = f[(f.sequence != "") & ~f.is_decoy]


def ratio(ints):
    d = {x["label"]: x["intensity"] for x in ints}
    lo, hi = d.get("SILAC light", 0.0), d.get("SILAC heavy", 0.0)
    return np.log2(hi / lo) if lo > 0 and hi > 0 else np.nan


f["o"] = f.intensities.map(ratio)
f = f.dropna(subset=["o"]).groupby(["run_file_name", "sequence", "charge"]).o.median().reset_index()
a = pq.read_table(lab).to_pandas()
a = a[(a.channel == 2) & a.log2_ratio.notna()]
a = a.groupby(["filename", "sequence", "charge"]).log2_ratio.median().reset_index()
iqr = lambda x: x.quantile(0.75) - x.quantile(0.25)
for rd in designs:
    run, d = rd.split("=")
    d = float(d)
    o = f[f.run_file_name == run]
    n = a[a.filename == run]
    b = o.merge(n, left_on=["sequence", "charge"], right_on=["sequence", "charge"])
    err = lambda x: np.median(np.abs(x - d))
    print(f"{run}: MS1LabeledWorkflow n={len(o)} median {o.o.median():.3f} IQR {iqr(o.o):.3f} |error| {err(o.o):.3f} | "
          f"andes n={len(n)} median {n.log2_ratio.median():.3f} IQR {iqr(n.log2_ratio):.3f} |error| {err(n.log2_ratio):.3f}")
    print(f"{run} same {len(b)} peptides: MS1LabeledWorkflow median {b.o.median():.3f} IQR {iqr(b.o):.3f} |error| {err(b.o):.3f} | "
          f"andes median {b.log2_ratio.median():.3f} IQR {iqr(b.log2_ratio):.3f} |error| {err(b.log2_ratio):.3f}")
