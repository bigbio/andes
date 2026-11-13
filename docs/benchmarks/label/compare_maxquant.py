"""Compare andes label.parquet log-ratios with MaxQuant evidence.txt ratios.

MaxQuant: per evidence (peptide, charge, raw file) the `Ratio M/L` / `Ratio H/L`
columns (non-normalised), log2; contaminants and reverse hits left out.
andes: per (sequence, charge, run) channel log2_ratio. Both are matched on
the bare sequence and charge within a run.

usage: compare_maxquant.py LABEL_PARQUET EVIDENCE_TXT RUN DESIGN [DESIGN ...]
  DESIGN: per channel ratio (2 for M/L ... in channel order), log2 expected
  values, or 'class' to split lysine / arginine peptides (PXD003327).
"""
import sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

lab, ev, run, *design = sys.argv[1:]
t = pq.read_table(lab).to_pandas()
t = t[t.filename == run]
e = pd.read_csv(ev, sep="\t", low_memory=False)
e = e[(e["Raw file"] == run)]
for col in ("Reverse", "Potential contaminant"):
    if col in e:
        e = e[e[col] != "+"]
ratio_cols = [c for c in ("Ratio M/L", "Ratio H/L") if c in e]
iqr = lambda x: x.quantile(.75) - x.quantile(.25)
for k, rc in enumerate(ratio_cols, start=2):
    if len(ratio_cols) == 1:
        k = 2
    m = e[["Sequence", "Charge", rc]].dropna()
    m = m[m[rc] > 0].copy()
    m["mq"] = np.log2(m[rc])
    # one value per (sequence, charge): the median of its evidences
    m = m.groupby(["Sequence", "Charge"]).mq.median().reset_index()
    a = t[(t.channel == k) & t.log2_ratio.notna()][["sequence", "charge", "log2_ratio"]]
    a = a.groupby(["sequence", "charge"]).log2_ratio.median().reset_index()
    both = m.merge(a, left_on=["Sequence", "Charge"], right_on=["sequence", "charge"])
    name = rc.replace("Ratio ", "")
    if design and design[0] == "class":
        for cls, sel in (("K", lambda s: s.str.contains("K") & ~s.str.contains("R")),
                         ("R", lambda s: s.str.contains("R") & ~s.str.contains("K"))):
            mm, aa, bb = m[sel(m.Sequence)], a[sel(a.sequence)], both[sel(both.sequence)]
            print(f"{run} {name} {cls}: MaxQuant n={len(mm)} median {mm.mq.median():.2f} IQR {iqr(mm.mq):.2f} | "
                  f"andes n={len(aa)} median {aa.log2_ratio.median():.2f} IQR {iqr(aa.log2_ratio):.2f} | "
                  f"both n={len(bb)} MaxQuant {bb.mq.median():.2f}/{iqr(bb.mq):.2f} andes {bb.log2_ratio.median():.2f}/{iqr(bb.log2_ratio):.2f} "
                  f"median |diff| {np.median(np.abs(bb.mq - bb.log2_ratio)):.3f}")
    else:
        d = float(design[min(k - 2, len(design) - 1)]) if design else 0.0
        err = lambda x: np.median(np.abs(x - d))
        print(f"{run} {name}: MaxQuant n={len(m)} median {m.mq.median():.3f} IQR {iqr(m.mq):.3f} |error| {err(m.mq):.3f} | "
              f"andes n={len(a)} median {a.log2_ratio.median():.3f} IQR {iqr(a.log2_ratio):.3f} |error| {err(a.log2_ratio):.3f}")
        print(f"{run} {name} same {len(both)} peptides: MaxQuant median {both.mq.median():.3f} IQR {iqr(both.mq):.3f} |error| {err(both.mq):.3f} | "
              f"andes median {both.log2_ratio.median():.3f} IQR {iqr(both.log2_ratio):.3f} |error| {err(both.log2_ratio):.3f}")
