"""Fragment-intensity and RT baselines on confident PSMs: PeptDeep vs the andes
frag-intensity GBDT vs a flat prediction; PeptDeep iRT vs the andes RT index.

Inputs in the working directory, for tag T:
  T_conf_pdeep.tsv  confident target PSMs (SpecId, peptide, q) joined with the
                    per-row features written by onnx-infer's rescore_pin example
  T_gbdt.tsv        output of model-train's frag_model_baseline example

RT calibrations are fitted on half of the peptides and evaluated on the other
half. The binned-median curve is the one in peptdeep_pin_arms.py; it is not
robust to a single outlier bin (see NOTES.md).

Usage: python compare_frag_baseline.py T
"""
import sys, numpy as np, pandas as pd, csv
t = sys.argv[1]
pd_ = pd.read_csv(f"{t}_conf_pdeep.tsv", sep="\t", quoting=csv.QUOTE_NONE, na_values="NA")
gb = pd.read_csv(f"{t}_gbdt.tsv", sep="\t", quoting=csv.QUOTE_NONE)
n0 = len(pd_)
assert pd_.SpecId.is_unique and gb.SpecId.is_unique
d = pd_.merge(gb, on="SpecId")
d = d[d.pdeep_missing == 0]
print(f"== {t}: confident PSMs {n0}, with GBDT row {len(gb)}, joined & PeptDeep-predicted {len(d)}")
# sanity: the observed RT from both tools must agree (same scan)
assert np.nanmax(np.abs(d.observed_retention_time - d.obs_rt_min)) < 1e-3, "scan mismatch"
def row(name, sa, cos, pr):
    print(f"{name:10s} SA median {np.median(sa):.3f} (IQR {np.quantile(sa,.25):.3f}-{np.quantile(sa,.75):.3f})  "
          f"cos {np.median(cos):.3f}  pearson {np.median(pr):.3f}  frac SA>0.7 {np.mean(sa>0.7):.3f}")
row("PeptDeep", d.spectral_angle, d.cos_norm, d.spec_pearson_norm)
row("andesGBDT", d.gbdt_sa, d.gbdt_cos, d.gbdt_pearson)
row("flat", d.flat_sa, d.flat_cos, d.flat_pearson)
print(f"GBDT SA > PeptDeep SA on {np.mean(d.gbdt_sa > d.spectral_angle):.3f} of PSMs; corr(SA_pdeep, SA_gbdt) {np.corrcoef(d.spectral_angle, d.gbdt_sa)[0,1]:.3f}")
for z, g in d.groupby("charge"):
    if len(g) >= 100:
        print(f"  z={z} n={len(g):6d} SA PeptDeep {np.median(g.spectral_angle):.3f} GBDT {np.median(g.gbdt_sa):.3f} flat {np.median(g.flat_sa):.3f}")
d["L"] = d.peptide.str.replace(r"[^A-Z]", "", regex=True).str.len() - 2
for lo, hi in [(6, 10), (11, 15), (16, 20), (21, 60)]:
    g = d[(d.L >= lo) & (d.L <= hi)]
    print(f"  len {lo}-{hi} n={len(g):6d} SA PeptDeep {np.median(g.spectral_angle):.3f} GBDT {np.median(g.gbdt_sa):.3f}")

# ---- RT: hold out by peptide (calibrate on half the peptides, evaluate on the rest) ----
d["seq"] = d.peptide.str[2:-2]
peps = d.seq.unique(); rng = np.random.default_rng(0); rng.shuffle(peps)
tr = d.seq.isin(set(peps[: len(peps) // 2])).values; te = ~tr
obs = d.obs_rt_min.values
def binned(x, y, bins=50):
    qs = np.quantile(x, np.linspace(0, 1, bins + 1)); xs, ys = [], []
    for lo, hi in zip(qs[:-1], qs[1:]):
        m = (x >= lo) & (x <= hi)
        if m.sum() >= 5: xs.append(np.median(x[m])); ys.append(np.median(y[m]))
    xs, ys = np.array(xs), np.maximum.accumulate(np.array(ys))
    return lambda q: np.interp(q, xs, ys)
def rep(name, x):
    lin = np.polyfit(x[tr], obs[tr], 1); r_lin = np.abs(np.polyval(lin, x[te]) - obs[te])
    r_bin = np.abs(binned(x[tr], obs[tr])(x[te]) - obs[te])
    print(f"RT {name:22s} linear: MAE {r_lin.mean():.2f} median {np.median(r_lin):.2f} p90 {np.quantile(r_lin,.9):.2f} | "
          f"binned-median: MAE {r_bin.mean():.2f} median {np.median(r_bin):.2f} p90 {np.quantile(r_bin,.9):.2f} min  (r={np.corrcoef(x[te], obs[te])[0,1]:.3f})")
print(f"RT range {obs.min():.1f}-{obs.max():.1f} min, held-out PSMs {te.sum()}")
rep("PeptDeep iRT", d.pdeep_irt.values)
rep("andes RT seed index", d.rt_seed_index.values)
F = np.array([np.array(s.split(","), float) for s in d.rt_features])
X = np.c_[F, np.ones(len(F))]
coef, *_ = np.linalg.lstsq(X[tr], obs[tr], rcond=None)
rep("andes RT index refit*", X @ coef)
print("  *refit = RtIndexModel-style least squares on rt_features fitted on the calibration half (in-run, best case)")
