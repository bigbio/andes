#!/usr/bin/env python3
"""One downstream pipeline, identical for every arm, on RAW precursor x run tables (one per arm).

  python -m quant_eval.downstream --name DS --design design.tsv [--amounts amounts.tsv] --fasta db.fasta --classes A,B --null A \
      --andes NAME:RULE:finalize.parquet ... --maxquant NAME:evidence.txt --fragpipe NAME:combined_ion.tsv --out DS.json

Arms (the only arm-specific step is the row rule of the reader, see readers.py):
  --andes NAME:RULE:PATH   RULE in default | all | hard-cut | weights | weights-bias; `weights` is the precision-weight view (V),
                           `weights-bias` the same view with the table's per-transfer bias and its se applied (readers.py)
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
from . import rates as RT
from . import readers as R
from . import runguard as G
from .de import de_all
from .design import Design
from .fasta import map_precursors
from .rollup import (ROLLUPS, adjust_raw, apply_celltype, balance_weights, celltype_offsets, identified_sigma, learn_weights,
                     normalise, normalise_rt_local, rollup)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help="dataset label written into the result")
    ap.add_argument("--design", required=True, help="design.tsv (see design.py)")
    ap.add_argument("--amounts", help="amounts.tsv; without it only counts and absence metrics are computed")
    ap.add_argument("--fasta", help="target database for species and protein of each precursor")
    ap.add_argument("--classes", default="", help="comma-separated species classes (HUMAN, YEAST, ECOLI, UPS)")
    ap.add_argument("--null", help="species whose amount does not change (false positives are its significant proteins)")
    ap.add_argument("--low-max", type=float, default=500.0, help="UPS pairs whose lower level is <= this count as low-level pairs")
    ap.add_argument("--andes", action="append", default=[], metavar="NAME:RULE[@LIBQ]:PATH",
                    help="RULE in default | all | hard-cut | weights | weights-bias; @LIBQ sets the library_q limit of the arm (default 0.01), e.g. weights@0.05")
    ap.add_argument("--maxquant", action="append", default=[], metavar="NAME:PATH")
    ap.add_argument("--fragpipe", action="append", default=[], metavar="NAME:PATH")
    ap.add_argument("--flashlfq", action="append", default=[], metavar="NAME:PEAKS[:PIPQ]",
                    help="FlashLFQ QuantifiedPeaks.tsv (PIP-ECHO MBR FDR); optional PIPQ keeps MBR peaks with PIP q-value <= PIPQ")
    ap.add_argument("--rt-local", action="store_true",
                    help="also run DE after the RT-local run normalisation (same function for every arm; RT proxy = the first andes arm's apex RT)")
    ap.add_argument("--norm-null", action="store_true",
                    help="also run DE with run factors from the null species only (the species known constant; same step for every arm)")
    ap.add_argument("--rollup", default="median-polish", choices=ROLLUPS,
                    help="protein roll-up, the same for every arm: median-polish (default); celltype = per-precursor offset of the "
                         "transferred cells, learned inside conditions and shrunk to a global offset, removed first; balance = "
                         "transfers re-weighted per (protein, condition) to the protein's pooled transfer share; or both; censored = transferred cells treated as right-censored at their value (EM with the identified cells' learned residual SD)")
    ap.add_argument("--hard-cut", type=float, default=0.01, help="quant_q limit of the hard-cut rule")
    ap.add_argument("--out", required=True)
    return ap.parse_args(argv)


def is_weighted(kind):
    """Arms whose roll-up uses learned precision weights (the V views)."""
    return kind.get("rule") in ("weights", "weights-bias")


def prepare(d, kind, how, prot, cond_of, runs):
    """Normalised rows of one arm -> rows the roll-up takes (learned weights for V views, then the --rollup arm).
    Returns (rows, cell-type offsets or None, summary dict)."""
    info, off = {}, None
    if "celltype" in how:
        off, info["celltype"] = celltype_offsets(d, cond_of)
        d = apply_celltype(d, off)
    if is_weighted(kind):
        d, info["weight_variance_by_class"] = learn_weights(d, prot, runs)
    if "balance" in how:
        d = balance_weights(d, prot, cond_of)
    if how == "censored":
        info["censor_sigma"] = round(identified_sigma(d, prot, runs), 5)
    return d, off, info


def load_arms(a, runs):
    tabs, kinds = {}, {}
    for spec in a.andes:
        name, rule, path = spec.split(":", 2)
        rule, _, lq = rule.partition("@")
        tabs[name] = R.read_andes(path, runs, rule, a.hard_cut, float(lq) if lq else R.IDQ)
        kinds[name] = {"kind": "andes", "rule": rule, "library_q": float(lq) if lq else R.IDQ}
    for spec in a.maxquant:
        name, path = spec.split(":", 1)
        tabs[name] = R.read_maxquant(path, runs)
        kinds[name] = {"kind": "maxquant"}
    for spec in a.fragpipe:
        name, path = spec.split(":", 1)
        tabs[name] = R.read_fragpipe(path, runs)
        kinds[name] = {"kind": "fragpipe"}
    for spec in a.flashlfq:
        parts = spec.split(":")
        name, path = parts[0], parts[1]
        pq = float(parts[2]) if len(parts) > 2 else None
        tabs[name] = R.read_flashlfq(path, runs, pq)
        kinds[name] = {"kind": "flashlfq", "pip_q": pq}
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

    RES = {"dataset": a.name, "rollup": a.rollup, "runs": len(runs), "arms": {}, "arm_roles": kinds, "primary": primary, "shared_arms": shared_arms,
           "shared_cells": len(shared), "tags": {n: G.tag_of(r) for n, r in tabs.items()}, "classes": list(classes)}
    if design.has_de:
        RES["levels"] = [str(x) for x in design.levels]
    PM = {}
    for n, raw in tabs.items():
        o = M.basic_counts(raw, prot, sp, prot_sp, classes)
        if design.has_de:
            d, flag, nfull = normalise(raw, runs)
            d, off, info = prepare(d, kinds[n], a.rollup, prot, design.cond, runs)
            o.update(info)
            cntp = raw.groupby(raw.key.map(prot)).key.nunique()
            Pm = rollup(d, prot, runs, info.get("censor_sigma"))
            P = Pm[cntp.reindex(Pm.index).fillna(0).values >= 2]
            PM[n] = (P, d)
            o["acc_own_sp"] = M.accuracy(d, design, sp, shared, True)
            o["acc_shared_sp"] = M.accuracy(d, design, sp, shared, False)
            if "UPS" in classes:
                o["ups_strat"] = M.ups_by_level(raw, P, design, sp, prot, prot_sp)
                o["ups_design"] = M.design_detection(raw if off is None else adjust_raw(raw, off), design, sp)
        RES["arms"][n] = o
        print(n, {k: o[k] for k in ("rows", "precursors", "proteins_ge2")}, round(time.time() - t0), "s", file=sys.stderr)

    if design.has_de:
        # DE after every arm's protein table exists: the shared arms define the intersection / union BH families
        fam_arms = [n for n in shared_arms if n in PM]
        inter = set.intersection(*[set(PM[n][0].index) for n in fam_arms])
        union = set().union(*[set(PM[n][0].index) for n in fam_arms])
        RES["families"] = {"arms": fam_arms, "intersection": len(inter), "union": len(union),
                           "own": {n: int(len(PM[n][0])) for n in PM}}
        for n, (P, _) in PM.items():
            o = RES["arms"][n]
            o["de"], o["pairs"], o["prot_err"], o["de_family"] = de_all(P, design, prot_sp, classes, inter, union)

    def alt_pass(label, norm):
        """DE of every arm after another run normalisation `norm(raw) -> (table, flag, n_ref)`; own intersection family."""
        PA = {}
        for n, raw in tabs.items():
            d, _, nref = norm(raw)
            d, _, info = prepare(d, kinds[n], a.rollup, prot, design.cond, runs)
            cntp = raw.groupby(raw.key.map(prot)).key.nunique()
            Pm = rollup(d, prot, runs, info.get("censor_sigma"))
            PA[n] = (Pm[cntp.reindex(Pm.index).fillna(0).values >= 2], nref)
        fam_arms = [n for n in shared_arms if n in PA]
        inter = set.intersection(*[set(PA[n][0].index) for n in fam_arms])
        union = set().union(*[set(PA[n][0].index) for n in fam_arms])
        fams = {"arms": fam_arms, "intersection": len(inter), "union": len(union)}
        for n, (P, nref) in PA.items():
            de, _, _, fam = de_all(P, design, prot_sp, classes, inter, union)
            RES["arms"][n][label] = {"de": de, "de_family": fam, "proteins": int(len(P)), "reference_precursors": nref}
        print(label, "done", round(time.time() - t0), "s", file=sys.stderr)
        return fams

    if design.has_de and a.rt_local:
        first = next(s for s in a.andes).split(":", 2)[2]
        rt_of = R.read_andes_rt(first, runs)
        RES["rtlocal_families"] = alt_pass("rtlocal", lambda raw: normalise_rt_local(raw, runs, rt_of))
        RES["rtlocal_families"]["rt_keys"] = int(len(rt_of))

    if design.has_de and a.norm_null:
        null_keys = set(sp.index[sp.values == design.null]) if a.null else set()
        if not null_keys:
            raise SystemExit("--norm-null needs --null and a species map (--fasta)")
        RES["nullnorm_families"] = alt_pass("nullnorm", lambda raw: normalise(raw, runs, null_keys))
        RES["nullnorm_families"]["null"] = design.null

    if design.has_de:
        ref_names = [n for n in shared_arms if n in PM]
        RES["cv"] = M.matched_sd(PM, ref_names, design, sp, prot_sp, classes, runs)
        if len(classes) >= 3:
            refs = [n for n in ([primary] + [n for n in tools if kinds[n]["kind"] == "maxquant"]) if n in PM]
            RES["hye_det"] = M.detection_by_quartile(tabs, PM, refs, sp, classes,
                                                     all_rows_refs=[n for n in refs if kinds[n]["kind"] != "andes"])
    RES["absence"] = {n: M.absence(raw, design, sp) for n, raw in tabs.items()}
    if design.absent:
        # truth-measured false-transfer rates (rates.py); blank-run carry-over is reported apart from the absence score
        RES["rates"] = {n: RT.false_transfer_rates(raw, design, sp) for n, raw in tabs.items()}
        RES["depth_curves"] = {n: {"transfer_q": RT.depth_curve(raw, design, sp, "tq"), "quant_q": RT.depth_curve(raw, design, sp, "thr")}
                               for n, raw in tabs.items() if kinds[n].get("rule") in ("all", "weights", "weights-bias") or kinds[n]["kind"] == "flashlfq"}
    if design.blank:
        RES["carry_over"] = {n: RT.carry_over(raw, design, sp) for n, raw in tabs.items()}
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
