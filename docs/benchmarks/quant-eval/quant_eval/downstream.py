#!/usr/bin/env python3
"""One downstream pipeline, identical for every arm, on RAW precursor x run tables (one per arm).

  python -m quant_eval.downstream --name DS --design design.tsv [--amounts amounts.tsv] --fasta db.fasta --classes A,B --null A \
      --andes NAME:RULE:finalize.parquet ... --maxquant NAME:evidence.txt --fragpipe NAME:combined_ion.tsv --out DS.json

Arms (the only arm-specific step is the row rule of the reader, see readers.py):
  --andes NAME:RULE:PATH   RULE in default | all | hard-cut | weights; `weights` is the precision-weight view (V)
  --maxquant NAME:PATH     evidence.txt Intensity;  --fragpipe NAME:PATH   combined_ion.tsv with IonQuant normalisation off
Every arm table carries one run tag and no duplicated (key, run); no column of one arm is joined to another arm's table.

Normalisation and roll-up: rollup.py. DE: de.py (moderated t, BH 5 % per pair, sign-aware TP, with and without |log2FC| >= 0.5).
Precursor accuracy: leave-one-condition-out error, share |err| > 1, on the cells every shared arm quantifies ('shared') and on
the arm's own rows ('own'). Shared arms = the first andes arm and the tool arms."""
import argparse
import json
import sys
import time
import warnings

import numpy as np
import pandas as pd

from . import metrics as M
from . import readers as R
from . import runguard as G
from .de import de_all
from .design import Design
from .fasta import map_precursors
from .rollup import learn_weights, normalise, rollup


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help="dataset label written into the result")
    ap.add_argument("--design", required=True, help="design.tsv (see design.py)")
    ap.add_argument("--amounts", help="amounts.tsv; without it only counts and absence metrics are computed")
    ap.add_argument("--fasta", help="target database for species and protein of each precursor")
    ap.add_argument("--classes", default="", help="comma-separated species classes (HUMAN, YEAST, ECOLI, UPS)")
    ap.add_argument("--null", help="species whose amount does not change (false positives are its significant proteins)")
    ap.add_argument("--low-max", type=float, default=500.0, help="UPS pairs whose lower level is <= this count as low-level pairs")
    ap.add_argument("--andes", action="append", default=[], metavar="NAME:RULE:PATH")
    ap.add_argument("--maxquant", action="append", default=[], metavar="NAME:PATH")
    ap.add_argument("--fragpipe", action="append", default=[], metavar="NAME:PATH")
    ap.add_argument("--hard-cut", type=float, default=0.01, help="quant_q limit of the hard-cut rule")
    ap.add_argument("--out", required=True)
    return ap.parse_args(argv)


def load_arms(a, runs):
    tabs, kinds = {}, {}
    for spec in a.andes:
        name, rule, path = spec.split(":", 2)
        tabs[name] = R.read_andes(path, runs, rule, a.hard_cut)
        kinds[name] = {"kind": "andes", "rule": rule}
    for spec in a.maxquant:
        name, path = spec.split(":", 1)
        tabs[name] = R.read_maxquant(path, runs)
        kinds[name] = {"kind": "maxquant"}
    for spec in a.fragpipe:
        name, path = spec.split(":", 1)
        tabs[name] = R.read_fragpipe(path, runs)
        kinds[name] = {"kind": "fragpipe"}
    if not tabs:
        raise SystemExit("no arm given")
    for k, d in tabs.items():
        G.tag_of(d)
        assert not d.duplicated(["key", "run"]).any(), k
    return tabs, kinds


def evaluate(a):
    t0 = time.time()
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    classes = tuple(c for c in a.classes.split(",") if c)
    design = Design(a.design, a.amounts, a.null, classes, a.low_max)
    runs = design.runs
    tabs, kinds = load_arms(a, runs)
    print(a.name, "arms", {k: (len(v), G.tag_of(v)) for k, v in tabs.items()}, file=sys.stderr)

    allkeys = sorted(set().union(*[set(v.key) for v in tabs.values()]))
    if a.fasta:
        sp, prot, ent_sp = map_precursors(allkeys, a.fasta, classes)
        prot_sp = {p: ent_sp.get(p) for p in prot.dropna().unique()}
    else:
        sp = pd.Series("NONE", index=allkeys)
        prot = pd.Series(np.nan, index=allkeys, dtype=object)
        prot_sp = {}

    primary = next(n for n, k in kinds.items() if k["kind"] == "andes") if any(k["kind"] == "andes" for k in kinds.values()) else None
    tools = [n for n, k in kinds.items() if k["kind"] != "andes"]
    shared_arms = ([primary] if primary else []) + tools
    cells = lambda d: set((d.key + "|" + d.run).values)
    shared = set.intersection(*[cells(tabs[n]) for n in shared_arms])
    print("rows", {n: len(r) for n, r in tabs.items()}, "shared cells", len(shared), file=sys.stderr)

    RES = {"dataset": a.name, "runs": len(runs), "arms": {}, "arm_roles": kinds, "primary": primary, "shared_arms": shared_arms,
           "shared_cells": len(shared), "tags": {n: G.tag_of(r) for n, r in tabs.items()}, "classes": list(classes)}
    if design.has_de:
        RES["levels"] = [str(x) for x in design.levels]
    PM = {}
    for n, raw in tabs.items():
        o = M.basic_counts(raw, prot, sp, prot_sp, classes)
        if design.has_de:
            d, flag, nfull = normalise(raw, runs)
            if kinds[n].get("rule") == "weights":
                d, o["weight_variance_by_class"] = learn_weights(d, prot, runs)
            cntp = raw.groupby(raw.key.map(prot)).key.nunique()
            Pm = rollup(d, prot, runs)
            P = Pm[cntp.reindex(Pm.index).fillna(0).values >= 2]
            PM[n] = (P, d)
            o["de"], o["pairs"], o["prot_err"] = de_all(P, design, prot_sp, classes)
            o["acc_own_sp"] = M.accuracy(d, design, sp, shared, True)
            o["acc_shared_sp"] = M.accuracy(d, design, sp, shared, False)
            if "UPS" in classes:
                o["ups_strat"] = M.ups_by_level(raw, P, design, sp, prot, prot_sp)
        RES["arms"][n] = o
        print(n, {k: o[k] for k in ("rows", "precursors", "proteins_ge2")}, round(time.time() - t0), "s", file=sys.stderr)

    if design.has_de:
        ref_names = [n for n in shared_arms if n in PM]
        RES["cv"] = M.matched_sd(PM, ref_names, design, sp, prot_sp, classes, runs)
        if len(classes) >= 3:
            refs = [n for n in ([primary] + [n for n in tools if kinds[n]["kind"] == "maxquant"]) if n in PM]
            RES["hye_det"] = M.detection_by_quartile(tabs, PM, refs, sp, classes,
                                                     all_rows_refs=[n for n in refs if kinds[n]["kind"] != "andes"])
    RES["absence"] = {n: M.absence(raw, design, sp) for n, raw in tabs.items()}
    RES["elapsed_s"] = time.time() - t0
    return RES


def main(argv=None):
    a = parse_args(argv)
    RES = evaluate(a)
    with open(a.out, "w") as f:
        json.dump(RES, f, indent=1, default=float)
    print("wrote", a.out, file=sys.stderr)


if __name__ == "__main__":
    main()
