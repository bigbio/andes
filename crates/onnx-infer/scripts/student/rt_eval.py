#!/usr/bin/env python3
"""Evaluate RT predictors on held-out peptide sequences of the qpx RT table.

Per run (up to --runs-per-dataset per dataset): calibration anchors are up to
--anchors non-test rows of that run (isotonic, ``rt_calibration``); the scored
rows are the run's "test"-split peptidoforms, whose sequences no model trained
on. Datasets named in --holdout-datasets were never trained on at all.

Predictors: each ``--model NAME=DIR`` (andes-rt checkpoints) and, with
--peptdeep DIR, the AlphaPeptDeep generic RT model on the same rows.

Writes per-row predictions (parquet) and prints MAE / median AE in minutes,
MAE as % of the run's 1-99% RT span, and R^2, by dataset group and PTM class.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from rt_calibration import fit_isotonic  # noqa: E402
from rt_model import RtCnn, parse_proforma, predict  # noqa: E402

PD_NAMES = {
    (4, "C"): "Carbamidomethyl@C", (35, "M"): "Oxidation@M", (21, "S"): "Phospho@S", (21, "T"): "Phospho@T",
    (21, "Y"): "Phospho@Y", (7, "N"): "Deamidated@N", (7, "Q"): "Deamidated@Q", (737, "K"): "TMT6plex@K",
}
PD_NTERM = {1: "Acetyl@Protein_N-term", 737: "TMT6plex@Any_N-term", 28: "Gln->pyro-Glu@Q^Any_N-term",
            27: "Glu->pyro-Glu@E^Any_N-term"}


def to_alphabase(form):
    seq, mods = parse_proforma(form)
    names, sites = [], []
    for pos, uid in mods:
        if pos == 0:
            n = PD_NTERM.get(uid)
        elif pos <= len(seq):
            n = PD_NAMES.get((uid, seq[pos - 1]))
        else:
            n = None
        if n is None:
            return None
        names.append(n)
        sites.append(str(pos))
    return seq, ";".join(names), ";".join(sites)


def peptdeep_rt(forms, pretrained_dir):
    from peptdeep.model.rt import AlphaRTModel
    dev = "gpu" if torch.cuda.is_available() else "cpu"
    m = AlphaRTModel(device=dev)
    m.load(os.path.join(pretrained_dir, "rt.pth"))
    out = np.full(len(forms), np.nan)
    conv = [to_alphabase(f) for f in forms]
    ok = [i for i, c in enumerate(conv) if c is not None and len(c[0]) <= 50]
    df = pd.DataFrame([conv[i] for i in ok], columns=["sequence", "mods", "mod_sites"])
    df["nAA"] = df.sequence.str.len()
    df["idx"] = ok
    df = m.predict(df, batch_size=4096)
    out[df["idx"].values] = df["rt_pred"].values
    return out


PTM_ORDER = [("UNIMOD:737]", "TMT"), ("UNIMOD:21]", "phospho"), ("UNIMOD:1]-", "acetyl-Nterm"),
             ("UNIMOD:28]", "pyro-Glu"), ("UNIMOD:27]", "pyro-Glu"), ("UNIMOD:7]", "deamidated"),
             ("UNIMOD:35]", "oxidation")]


def ptm_class(form):
    """First matching class in PTM_ORDER, else unmodified/CAM."""
    for tag, name in PTM_ORDER:
        if tag in form:
            return name
    return "unmodified/CAM"


def metrics(g, col):
    e = (g[col] - g.rt_min).abs()
    span = g.span.values
    ss = ((g.rt_min - g.groupby("run_key").rt_min.transform("mean")) ** 2).sum()
    return dict(n=len(g), MAE=e.mean(), MedAE=e.median(), rel_MAE_pct=100 * np.mean(e / span),
                R2=1 - ((g[col] - g.rt_min) ** 2).sum() / ss)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", required=True)
    p.add_argument("--model", action="append", default=[], help="NAME=DIR")
    p.add_argument("--peptdeep", help="pretrained generic dir with rt.pth")
    p.add_argument("--holdout-datasets", nargs="*", default=["PXD019263", "PXD005235", "PXD012255"])
    p.add_argument("--runs-per-dataset", type=int, default=12)
    p.add_argument("--anchors", type=int, default=2000)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    rng = np.random.default_rng(11)
    cats = ["dataset", "run", "peptidoform", "split"]
    t = pq.read_table(a.table, columns=cats + ["rt_min"], read_dictionary=cats).to_pandas()
    t["run_id"] = t.groupby(["dataset", "run"], observed=True, sort=True).ngroup()
    pick = []
    for ds, g in t.groupby("dataset", observed=True):
        runs = np.sort(g.run_id.unique())
        pick += list(rng.choice(runs, min(a.runs_per_dataset, len(runs)), replace=False))
    t = t[t.run_id.isin(set(pick))].copy()
    for c in cats:
        t[c] = t[c].astype(str)
    t["run_key"] = t.dataset + "/" + t.run
    parts = []
    for rk, g in t.groupby("run_key"):
        anc = g[g.split != "test"]
        anc = anc.sample(min(a.anchors, len(anc)), random_state=0)
        te = g[g.split == "test"]
        if len(te) < 10 or len(anc) < 50:
            continue
        span = np.quantile(g.rt_min, 0.99) - np.quantile(g.rt_min, 0.01)
        parts.append(pd.concat([anc.assign(role="anchor"), te.assign(role="eval")]).assign(span=span))
    d = pd.concat(parts, ignore_index=True)
    forms = d.peptidoform.unique()
    print(f"runs {d.run_key.nunique()}, anchor rows {(d.role == 'anchor').sum()}, eval rows {(d.role == 'eval').sum()}, "
          f"forms {len(forms)}", flush=True)
    preds = {}
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    for spec in a.model:
        name, mdir = spec.split("=", 1)
        m = RtCnn().to(dev)
        m.load_state_dict(torch.load(os.path.join(mdir, "model.pt"), map_location=dev))
        preds[name] = predict(m, list(forms), device=dev)
    if a.peptdeep:
        preds["PeptDeep"] = peptdeep_rt(list(forms), a.peptdeep)
    fi = {f: i for i, f in enumerate(forms)}
    ix = d.peptidoform.map(fi).values
    cols = []
    for name, v in preds.items():
        d[f"idx_{name}"] = v[ix]
        d[f"cal_{name}"] = np.nan
        cols.append(name)
    for rk, g in d.groupby("run_key"):
        an = g[g.role == "anchor"]
        for name in cols:
            ok = np.isfinite(an[f"idx_{name}"])
            if ok.sum() < 20:
                continue
            f = fit_isotonic(an[f"idx_{name}"][ok].values, an.rt_min[ok].values)
            d.loc[g.index, f"cal_{name}"] = f(g[f"idx_{name}"].values)
    ev = d[d.role == "eval"].copy()
    common = np.logical_and.reduce([np.isfinite(ev[f"cal_{n}"]) for n in cols])
    print(f"eval rows predicted by every model: {common.sum()} of {len(ev)}")
    ev = ev[common]
    ev["group"] = np.where(ev.dataset.isin(a.holdout_datasets), "held-out dataset", "trained dataset")
    ev["ptm"] = ev.peptidoform.map(ptm_class)
    rows = []
    for key, g in list(ev.groupby("group")) + list(ev.groupby("ptm")) + list(ev.groupby("dataset")) + [("ALL", ev)]:
        for name in cols:
            rows.append(dict(slice=key, model=name, **metrics(g, f"cal_{name}")))
    res = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_rows", 500)
    print(res.pivot_table(index="slice", columns="model", values=["MAE", "rel_MAE_pct", "R2"]).round(3).to_string())
    print(res.groupby("slice").n.first().to_string())
    ev.to_parquet(a.out)


if __name__ == "__main__":
    main()
