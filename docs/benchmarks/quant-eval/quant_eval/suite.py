#!/usr/bin/env python3
"""Scoreboard suite: one manifest -> one result JSON per (dataset, seed) -> one scoreboard (board.py).

  python -m quant_eval.suite --manifest manifest.json --out results/ [--only u9740,qe] [--seeds 42,1,7]
  python -m quant_eval.board --results results/ > scoreboard.md
  python -m quant_eval.gate  --results results/ --candidate cand [--baseline main]

The manifest (JSON; relative paths resolve against the manifest's folder) names, per dataset, the design files, the database, the tool outputs
that already exist and, per andes build and quant-model seed, the finalize parquet (and the masked-run files of the wrong-peak measurement).
The suite reads finished tables only; it launches nothing. `runguard` tags every table (a single lineage per arm).

  {"reference_build": "main", "builds": ["main", "cand"],
   "datasets": {"u9740": {"label": "PXD009740", "held_out": false, "design": "...", "amounts": "...", "fasta": "...", "classes": "UPS,YEAST", "null": "YEAST",
       "accuracy_const": null,
       "tools": {"maxquant": {"MQ": "evidence.txt"}, "fragpipe": {"FP": "combined_ion.tsv"}, "flashlfq": {"FL001": "QuantifiedPeaks.tsv:0.01"}, "sage": {"Sage": "lfq.tsv"}},
       "andes": {"main": {"42": {"finalize": "...mbr.parquet", "time_files": ["finalize.time"],
                                 "masked": {"full_dir": "...", "masked_dir": "...", "pairs": "mask_pairs.parquet", "finalize": "...mbr.parquet"}}}}}}}

Arms of a build `b`: `b` (transfer_q <= 0.01), `b_V` (precision weights), `b_H` (hard cut quant_q <= 0.01), `b_I` (identified rows only: the same run without match between runs) and, when a matched cut exists, `b_M`
(hard cut at the matched true-error cut, matched.py). Datasets flagged "held_out": true are refused unless --held-out-eval is given; a held-out
dataset never enters a gate or a tuning run (see gate.py)."""
import argparse
import glob
import json
import os
import re
import sys
import time

from . import downstream as DS
from . import masked as MK
from . import matched as MT
from .design import Design

SEEDS = ("42", "1", "7")
CUTS = (False, True)


def _p(base, x):
    return x if (x is None or os.path.isabs(x)) else os.path.join(base, x)


def parse_time(paths):
    """Sum CPU seconds, max wall seconds and max RSS (GB) over `/usr/bin/time -v` outputs."""
    cpu = wall = rss = 0.0
    n = 0
    for pth in paths:
        for f in glob.glob(pth):
            u = s = w = r = None
            for line in open(f, errors="replace"):
                if "User time (seconds)" in line:
                    u = float(line.split(":")[-1])
                elif "System time (seconds)" in line:
                    s = float(line.split(":")[-1])
                elif "Elapsed (wall clock)" in line:
                    t = line.strip().split()[-1].split(":")
                    w = sum(float(x) * 60 ** i for i, x in enumerate(reversed(t)))
                elif "Maximum resident set size" in line:
                    r = float(line.split(":")[-1]) / 1048576
            if u is not None:
                cpu += u + (s or 0)
                wall = max(wall, w or 0)
                rss = max(rss, r or 0)
                n += 1
    return {"files": n, "cpu_s": cpu, "wall_s": wall, "rss_gb": rss}


def tool_args(D, base):
    a = []
    for kind, flag in (("maxquant", "--maxquant"), ("fragpipe", "--fragpipe"), ("flashlfq", "--flashlfq"), ("sage", "--sage")):
        for name, pth in (D.get("tools") or {}).get(kind, {}).items():
            if kind == "flashlfq":
                f, _, q = pth.partition(":")
                pth = _p(base, f) + (":" + q if q else "")
            else:
                pth = _p(base, pth)
            a += [flag, f"{name}:{pth}"]
    # the match-q filtered MaxQuant arm: the manifest's `maxquant_matchq` {name: evidence.txt} (an evidence table that carries a filled Match q-value),
    # by default the evidence table of the as-reported arm `MQ`
    tools = D.get("tools") or {}
    mq_q = tools.get("maxquant_matchq")
    if mq_q is None and "MQ" in tools.get("maxquant", {}):
        mq_q = {"MQ_q01": tools["maxquant"]["MQ"]}
    for name, pth in (mq_q or {}).items():
        a += ["--maxquant-matchq", f"{name}:{_p(base, pth)}:0.01"]
    return a


def load_manifest(path):
    M = json.load(open(path))
    M["_base"] = os.path.dirname(os.path.abspath(path))
    return M


def absence_pre(M, seed, out):
    """quant_q absence curve of every build on the absence set (PXD014415) for this seed -> matched cuts of the other sets."""
    for dn, D in M["datasets"].items():
        if D.get("absent_set") and not D.get("held_out"):
            base = M["_base"]
            des = Design(_p(base, D["design"]), None, None, tuple(c for c in D["classes"].split(",") if c))
            curves = {}
            for b in M["builds"]:
                e = D["andes"].get(b, {}).get(seed)
                if e:
                    curves[b] = MT.absence_curve(_p(base, e["finalize"]), des, _p(base, D["fasta"]), des.classes)
            return curves
    return {}


_WRONG = {}
_ABS = {}


def wrong_curves_of(M, dn, seed, builds):
    """masked wrong-peak curve of each build of one dataset and seed (cached)."""
    D = M["datasets"][dn]
    base = M["_base"]
    out = {}
    for b in builds:
        mk = D["andes"].get(b, {}).get(seed, {}).get("masked")
        if mk:
            if (dn, seed, b) not in _WRONG:
                lt = MK.labelled_transfers(_p(base, mk["full_dir"]), _p(base, mk["masked_dir"]), _p(base, mk["pairs"]), _p(base, mk["finalize"]))
                _WRONG[(dn, seed, b)] = MK.wrong_curve(lt)
            out[b] = _WRONG[(dn, seed, b)]
    return out


def reference_noise(M, dn, seeds=None):
    """Seed-to-seed spread (max - min) of the reference build's masked wrong-peak rate at quant_q <= 0.01 and, for the absence set, of its absence rate: the
    noise the matched cut must not react to. Computed over every seed of the manifest's reference build; 0 with fewer than two."""
    ref = M.get("reference_build", M["builds"][0])
    D = M["datasets"][dn]
    base = M["_base"]
    tol = {"wrong": 0.0, "absence": 0.0}
    # every seed the manifest holds for the reference build, whatever --seeds selects: the noise is a property of the reference
    seeds = sorted(s for s in D.get("andes", {}).get(ref, {}) if D["andes"][ref][s])
    if len(seeds) < 2:
        return tol
    rates = [MK.rate_at(wrong_curves_of(M, dn, s, [ref])[ref], MT.NOMINAL) for s in seeds if ref in wrong_curves_of(M, dn, s, [ref])]
    rates = [r for r in rates if r == r]
    if len(rates) >= 2:
        tol["wrong"] = max(rates) - min(rates)
    if D.get("absent_set"):
        des = Design(_p(base, D["design"]), None, None, tuple(c for c in D["classes"].split(",") if c))
        ar = []
        for s in seeds:
            if (dn, s, ref) not in _ABS:
                _ABS[(dn, s, ref)] = MT.absence_curve(_p(base, D["andes"][ref][s]["finalize"]), des, _p(base, D["fasta"]), des.classes)
            ar.append(MT.rate_on_curve(_ABS[(dn, s, ref)], D.get("absent_species", "YEAST")))
        tol["absence"] = max(ar) - min(ar)
    return tol


def run_one(M, dn, seed, out, abs_curves, held_out_eval=False, tol=None):
    D = M["datasets"][dn]
    if D.get("held_out") and not held_out_eval:
        raise SystemExit(f"{dn} is held out: refused (explicit --held-out-eval is required, and it is never part of a gate or tuning run)")
    base = M["_base"]
    builds = [b for b in M["builds"] if D.get("andes", {}).get(b, {}).get(seed)]
    if not builds:
        return None
    ref = M.get("reference_build", M["builds"][0])
    wrong = wrong_curves_of(M, dn, seed, builds)
    cuts = {b: MT.matched_cut(b, ref, wrong, abs_curves, D.get("absent_species", "YEAST"), tol) for b in builds} if ref in builds else {}
    argv = ["--name", D.get("label", dn), "--design", _p(base, D["design"]), "--fasta", _p(base, D["fasta"]), "--classes", D["classes"], "--out", os.path.join(out, f"{dn}_s{seed}.json")]
    if D.get("amounts"):
        argv += ["--amounts", _p(base, D["amounts"])]
    if D.get("null"):
        argv += ["--null", D["null"], "--norm-null"]
    argv += ["--accuracy"]
    if D.get("accuracy_const"):
        argv += ["--accuracy-const", D["accuracy_const"]]
    for b in builds:
        fin = _p(base, D["andes"][b][seed]["finalize"])
        argv += ["--andes", f"{b}:default:{fin}", "--andes", f"{b}_V:weights:{fin}", "--andes", f"{b}_H:hard-cut:{fin}",
                 "--andes", f"{b}_I:ident:{fin}"]
        t = cuts.get(b, {}).get("tau")
        if t is not None and b != ref:
            argv += ["--andes", f"{b}_M:hard-cut={t}:{fin}"]
        elif b == ref:
            pass
    argv += tool_args(D, base)
    t0 = time.time()
    RES = DS.evaluate(DS.parse_args(argv))
    RES.update({"suite_dataset": dn, "seed": seed, "held_out": bool(D.get("held_out")), "builds": builds, "reference_build": ref,
                "wrong_peak": wrong, "matched_cuts": cuts, "matched_tolerance": tol, "absence_curves": abs_curves if D.get("absent_set") else None,
                "cost": {b: parse_time([_p(base, f) for f in D["andes"][b][seed].get("time_files", [])]) for b in builds},
                "suite_elapsed_s": time.time() - t0})
    with open(os.path.join(out, f"{dn}_s{seed}.json"), "w") as f:
        json.dump(RES, f, indent=1, default=float)
    return RES


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", help="comma-separated dataset keys (default: every dataset that is not held out)")
    ap.add_argument("--seeds", default=",".join(SEEDS))
    ap.add_argument("--held-out-eval", action="store_true",
                    help="explicit held-out evaluation mode: only then may held-out datasets be run; their results are marked and refused by gate.py")
    a = ap.parse_args(argv)
    M = load_manifest(a.manifest)
    os.makedirs(a.out, exist_ok=True)
    want = a.only.split(",") if a.only else [d for d, D in M["datasets"].items() if not D.get("held_out")]
    for d in want:
        if d not in M["datasets"]:
            raise SystemExit(f"unknown dataset {d}")
        if M["datasets"][d].get("held_out") and not a.held_out_eval:
            raise SystemExit(f"{d} is held out: refused without --held-out-eval")
    noise = {d: reference_noise(M, d) for d in want if M["datasets"][d].get("andes")}
    for d in want:
        noise.setdefault(d, None)
    for seed in a.seeds.split(","):
        abs_curves = absence_pre(M, seed, a.out)
        for d in want:
            D = M["datasets"][d]
            if not D.get("andes") and not D.get("tools"):
                print(f"skip {d}: no inputs in the manifest", file=sys.stderr)
                continue
            print(f"== {d} seed {seed}", file=sys.stderr)
            run_one(M, d, seed, a.out, abs_curves, a.held_out_eval, noise[d])


if __name__ == "__main__":
    main()
