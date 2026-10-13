#!/usr/bin/env python3
"""Sampling weights and hold-out assignment for the andes-rt / andes-ms2 training corpus.

Input: the merged catalog ``datasets.parquet`` (``corpus_catalog.py merge``). Output: one row
per dataset with its stratum, split and per-PSM sampling weight, plus the stratum budget.

Plan (CORPUS.md, "Sampling plan"):
1. **Stratum** = activation x analyser x labelling x enrichment-class
   (enrichment-class: HLA-I, HLA-II, nonspecific/non-tryptic, PTM, none).
2. **Temperature sampling over strata**: stratum share p_s proportional to N_s ** alpha
   (alpha = 0.5 by default; 1 = natural mix, 0 = uniform), then multiplied by a priority boost
   (HLA and PTM strata, the user's priorities) and renormalised.
3. **Per-dataset cap inside a stratum**: no dataset contributes more than ``--max-dataset-share``
   of its stratum (a 26M-PSM dataset cannot drown the other 32).
4. **Per-precursor cap** (applied at selection time, ``ms2_select.py --max-per-precursor``):
   at most k spectra per (peptidoform, charge, activation, NCE bin); RT uses one median per
   (run, peptidoform).
5. **Hold-outs**: whole datasets (``--holdout-datasets``, default the andes-rt v1 three plus
   one per new stratum chosen by ``--holdout-per-stratum``) and whole PTM classes
   (``--holdout-ptms``: every PSM carrying one of these Unimod ids is test-only, in every
   dataset) - the zero-shot PTM evaluation. Benchmarks are dropped, never held out.

Usage::

  python corpus_sampling.py --catalog CAT/datasets.parquet --out weights.tsv \\
      --holdout-ptms UNIMOD:121 UNIMOD:1 --alpha 0.5
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_catalog import EXCLUDED, PREDICTOR_HELDOUT  # noqa: E402

BOOST = {"HLA-I": 4.0, "HLA-II": 3.0, "nonspecific": 2.0, "PTM": 2.5, "none": 1.0}
PTM_ENRICH = ("phospho", "glyco", "ubiquitin", "acetyl", "methyl", "acyl", "citrullination", "ADP-ribosyl")


def enrichment_class(row) -> str:
    e = str(row.get("enrichment") or "").lower()
    if "hla-ii" in e:
        return "HLA-II"
    if "hla" in e:
        return "HLA-I"
    if any(p.lower() in e for p in PTM_ENRICH):
        return "PTM"
    if str(row.get("enzyme_specificity") or "") == "nonspecific" or (
            row.get("tryptic_cterm_frac") is not None and not pd.isna(row.get("tryptic_cterm_frac"))
            and row.get("tryptic_cterm_frac") < 0.5):
        return "nonspecific"
    return "none"


ANALYSER_NORM = {"Orbitrap": "FTMS", "FTMS": "FTMS", "Astral": "ASTMS", "ASTMS": "ASTMS", "TOF(timsTOF)": "timsTOF",
                 "TOF": "TOF", "IonTrap": "ITMS", "ITMS": "ITMS", "FT-ICR/IT": "FTMS|ITMS?",
                 "Orbitrap/IT (tribrid)": "FTMS|ITMS?", "Orbitrap/IT (hybrid)": "FTMS|ITMS?"}


def stratum(row) -> str:
    frag = str(row.get("fragmentation") or "unknown").split("|")[0].split(" ")[0]
    an = str(row.get("analyser") or "unknown").split("|")[0]
    an = ANALYSER_NORM.get(an, an)
    lab = str(row.get("labelling") or "LFQ").split("|")[0].replace("LFQ?", "LFQ")
    return f"{frag}/{an}/{lab}/{enrichment_class(row)}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--catalog", required=True)
    p.add_argument("--runs", help="merged runs.parquet: split mixed-analyser datasets into per-analyser units")
    p.add_argument("--out", required=True)
    p.add_argument("--alpha", type=float, default=0.5)
    p.add_argument("--max-dataset-share", type=float, default=0.25)
    p.add_argument("--holdout-datasets", nargs="*", default=sorted(PREDICTOR_HELDOUT))
    p.add_argument("--holdout-per-stratum", type=int, default=1,
                   help="extra held-out datasets per stratum with >= 3 datasets (smallest-but-one by PSMs)")
    p.add_argument("--holdout-ptms", nargs="*", default=[])
    p.add_argument("--budget", type=float, default=50e6, help="training draws per epoch")
    a = p.parse_args()

    D = pd.read_parquet(a.catalog)
    D = D[~D.accession.str.split("-").str[0].isin(EXCLUDED)].copy()
    if a.runs:
        # a dataset whose runs differ in analyser (e.g. PXD009875: FTMS and ITMS runs) becomes one
        # sampling unit per analyser, so its ion-trap PSMs land in the ion-trap stratum
        R = pd.read_parquet(a.runs)
        R = R[R.analyser.notna()]
        units = []
        for acc, g in R.groupby("accession"):
            if g.analyser.nunique() < 2 or acc not in set(D.accession):
                continue
            base = D[D.accession == acc].iloc[0]
            for an, h in g.groupby("analyser"):
                u = base.copy()
                u["analyser"], u["n_psms"] = an, h.n_psms.sum()
                u["fragmentation"] = "CID?" if an == "ITMS" and not base.get("fragmentation") else base.get("fragmentation")
                u["accession"] = f"{acc}:{an}"
                units.append(u)
            D = D[D.accession != acc]
        if units:
            D = pd.concat([D, pd.DataFrame(units)], ignore_index=True)
    D["n_psms"] = D.n_psms.fillna(0).astype(float)
    D = D[D.n_psms > 0]
    D["stratum"] = D.apply(stratum, axis=1)
    D["eclass"] = D.apply(enrichment_class, axis=1)
    D["split"] = np.where(D.accession.str.split(":").str[0].isin(a.holdout_datasets), "test", "train")
    for s, g in D[D.split == "train"].groupby("stratum"):
        if len(g) >= 3 and a.holdout_per_stratum:
            pick = g.sort_values("n_psms").iloc[1:1 + a.holdout_per_stratum].index
            D.loc[pick, "split"] = "test"
    T = D[D.split == "train"]
    N = T.groupby("stratum").n_psms.sum()
    boost = pd.Series({s: BOOST[s.split("/")[-1]] for s in N.index})
    share = (N ** a.alpha) * boost
    share /= share.sum()
    budget = share * a.budget
    w = []
    for _, r in D.iterrows():
        if r.split != "train":
            w.append(0.0)
            continue
        g = T[T.stratum == r.stratum]
        capped = np.minimum(g.n_psms, a.max_dataset_share * g.n_psms.sum()) if len(g) > 1 else g.n_psms
        eff = min(r.n_psms, a.max_dataset_share * g.n_psms.sum()) if len(g) > 1 else r.n_psms
        w.append(float(budget[r.stratum] * eff / capped.sum() / r.n_psms))
    D["psm_weight"] = w
    D["expected_draws"] = D.psm_weight * D.n_psms
    cols = ["accession", "source", "stratum", "split", "n_psms", "psm_weight", "expected_draws"]
    D[cols].sort_values(["stratum", "n_psms"]).to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    summary = pd.DataFrame({"psms": N, "share_natural": N / N.sum(), "share_sampled": share})
    print(summary.sort_values("share_sampled", ascending=False).to_string(float_format=lambda x: f"{x:.3g}"))
    with open(os.path.splitext(a.out)[0] + ".holdout.json", "w") as fh:
        json.dump({"datasets": sorted(D.accession[D.split == "test"]), "ptm_classes": a.holdout_ptms,
                   "rule": "a PSM carrying any held-out PTM class is test-only in every dataset"}, fh, indent=1)


if __name__ == "__main__":
    main()
