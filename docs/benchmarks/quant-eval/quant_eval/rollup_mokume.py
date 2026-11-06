"""Protein roll-up through mokume (mokume-py components), the same function for every arm.

The evaluation harness keeps its own median polish (rollup.py, the default of every other module). This module only adapts the
raw precursor x run tables of the arms to the peptide table that mokume's quantification methods read (columns ProteinName,
PeptideCanonical, NormIntensity, SampleID; one row per precursor and run, the precursor key is the peptide identity) and calls
mokume's own normalisation and quantification classes. Nothing is re-implemented here.

Normalisation views
  harness   the harness run normalisation (rollup.normalise: median of log2 over the precursors quantified in all runs) applied
            first, the intensities are put back on the linear scale: isolates the roll-up
  mokume    mokume's default sample normalisation `globalMedian` (GlobalMedianNormalizer with the median map the pipeline builds
            in Feature.get_median_map: per-sample median intensity / median of those medians, over all rows of the sample)
Methods (mokume.quantification): maxlfq (MaxLFQQuantification, Cox et al. Eq. 3, min_ratio_count 2), directlfq
(DirectLFQQuantification, the directlfq package, its own sample normalisation off), top3 (TopNQuantification n=3, the mean of the
three most intense precursors), sum (AllPeptidesQuantification), median (MedianQuantification).
mokume returns linear intensities; zeros (a sample with no valid MaxLFQ edge) become missing. Output: protein x run, log2."""
import numpy as np
import pandas as pd

from . import rollup as RU

METHODS = ("maxlfq", "directlfq", "top3", "sum", "median")
NORMS = ("harness", "mokume")


def versions():
    import importlib.metadata as md
    out = {}
    for p in ("mokume-py", "mokume", "directlfq", "pyopenms", "numpy", "pandas"):
        try:
            out[p] = md.version(p)
        except md.PackageNotFoundError:
            pass
    return out


def peptide_table(raw, prot, runs, norm):
    """Long table of one arm for mokume: only precursors mapped to one protein, proteins with >= 2 precursors."""
    if norm == "harness":
        d, _, _ = RU.normalise(raw, runs)
        val = np.exp2(d.n.values)
    else:
        d, val = raw, raw.quantity.values
    t = pd.DataFrame({"ProteinName": d.key.map(prot).values, "PeptideCanonical": d.key.values, "SampleID": d.run.values,
                      "NormIntensity": val})
    t = t.dropna(subset=["ProteinName"])
    t = t[np.isfinite(t.NormIntensity) & (t.NormIntensity > 0)]
    n = t.groupby("ProteinName").PeptideCanonical.nunique()
    t = t[t.ProteinName.isin(n.index[n >= 2])].copy()
    if norm == "mokume":
        from mokume.normalization.sample_normalizers import GlobalMedianNormalizer
        med = t.groupby("SampleID").NormIntensity.median()
        mm = (med / np.median(med.values)).to_dict()
        g = GlobalMedianNormalizer()
        t = pd.concat([g.normalize(x.copy(), med_map=mm, sample=s) for s, x in t.groupby("SampleID", sort=False)], ignore_index=True)
    return t


def quantifier(method, threads=4):
    from mokume import quantification as Q
    if method == "maxlfq":
        return Q.MaxLFQQuantification(min_ratio_count=2, threads=threads)
    if method == "directlfq":
        return Q.DirectLFQQuantification(min_nonan=1, num_cores=threads, deactivate_normalization=True)
    if method == "top3":
        return Q.TopNQuantification(n=3)
    if method == "sum":
        return Q.AllPeptidesQuantification()
    if method == "median":
        from mokume.quantification.median import MedianQuantification
        return MedianQuantification()
    raise SystemExit(f"method {method!r}")


def protein_table(t, method, runs, threads=4):
    """mokume protein quantification of the peptide table -> log2 protein x run (NaN where mokume reports none)."""
    r = quantifier(method, threads).quantify(t, protein_column="ProteinName", peptide_column="PeptideCanonical",
                                             intensity_column="NormIntensity", sample_column="SampleID")
    col = "Intensity" if "Intensity" in r.columns else r.select_dtypes("number").columns[-1]
    w = r.pivot_table(index="ProteinName", columns="SampleID", values=col, aggfunc="first")
    w = w.reindex(columns=runs).astype(float)
    w = w.where(w > 0)
    return np.log2(w)
