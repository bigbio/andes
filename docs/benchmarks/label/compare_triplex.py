"""Compare andes triplex label ratios with OpenMS FeatureFinderMultiplex.

usage: compare_triplex.py LABEL_PARQUET FFM_TSV [RUN]
RUN is the label.parquet `filename` the TSV belongs to (needed when the table
holds several runs).
"""
import argparse, sys, numpy as np, pandas as pd, pyarrow.parquet as pq
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("label_parquet")
ap.add_argument("ffm_tsv")
ap.add_argument("run", nargs="?", help="the label.parquet filename of the TSV's run")
args = ap.parse_args()
lab, ffm = args.label_parquet, args.ffm_tsv
t = pq.read_table(lab).to_pandas()
if args.run:
    t = t[t.filename == args.run]
elif t.filename.nunique() > 1:
    sys.exit(f"{lab} holds {t.filename.nunique()} runs: name the run of {ffm}")
ref = t[t.channel == 1][["peptide", "charge", "sequence", "mz", "rt_apex", "status"]]
# OpenMS multiplets: light channel (map 0) m/z, and log2 ratios of maps 1, 2 to 0.
rows = []
for _, r in pd.read_csv(ffm, sep="\t").iterrows():
    ch = {int(a): (float(b), float(c)) for a, b, c in (x.split(":") for x in r.channels.split(";"))}
    if 0 not in ch: continue
    rows.append((r.rt, r.charge, ch[0][0], *(np.log2(ch[k][1] / ch[0][1]) if k in ch and ch[0][1] > 0 and ch[k][1] > 0 else np.nan for k in (1, 2))))
f = pd.DataFrame(rows, columns=["rt", "charge", "mz", "r1", "r2"]).sort_values("mz").reset_index(drop=True)
mzs = f.mz.values
out = []
for _, g in ref.iterrows():
    lo, hi = np.searchsorted(mzs, [g.mz * (1 - 10e-6), g.mz * (1 + 10e-6)])
    c = f.iloc[lo:hi]
    c = c[(c.charge == g.charge) & ((c.rt - g.rt_apex).abs() <= 30)]
    if len(c) == 0: continue
    b = c.iloc[(c.rt - g.rt_apex).abs().argmin()]
    out.append((g.peptide, g.charge, g.sequence, b.r1, b.r2))
m = pd.DataFrame(out, columns=["peptide", "charge", "sequence", "o1", "o2"])
a = (t[t.channel > 1]
     .pivot_table(index=["peptide", "charge"], columns="channel", values="log2_ratio", dropna=False)
     .reindex(columns=[2, 3])
     .rename(columns={2: "a1", 3: "a2"})
     .reset_index())
m = m.merge(a, on=["peptide", "charge"])
nK = m.sequence.str.count("K"); nR = m.sequence.str.count("R")
m["cls"] = np.where((nK > 0) & (nR == 0), "K", np.where((nR > 0) & (nK == 0), "R", "KR"))
print(f"andes groups {len(ref)}, OpenMS multiplets {len(f)}, matched {len(m)}")
for k, (oa, aa) in {"M/L": ("o1", "a1"), "H/L": ("o2", "a2")}.items():
    for cls in ("K", "R"):
        s = m[m.cls == cls]
        both = s[[oa, aa]].dropna()
        iqr = lambda x: x.quantile(.75) - x.quantile(.25)
        print(f"{k} {cls}: OpenMS n={s[oa].notna().sum()} median {s[oa].median():.2f} IQR {iqr(s[oa]):.2f} | andes n={s[aa].notna().sum()} median {s[aa].median():.2f} IQR {iqr(s[aa]):.2f} | both n={len(both)} median |diff| {np.median(np.abs(both[oa]-both[aa])):.3f}")
