#!/usr/bin/env python3
"""Protein-level sensitivity to the roll-up: every arm's RAW precursor intensities through the same set of roll-ups.

  python -m quant_eval.rollup_sens --manifest manifest.json --dataset u9740 --seed 42 --out out.json [--tool-only]

Arms (as in suite.py): the andes build `main` (transfers with transfer_q <= 0.01, i.e. with match between runs) and `main_I`
(identified rows only), MaxQuant `MQ` (match between runs) and `MQ_noMBR`, FragPipe `FP`. Roll-ups: `polish` (the harness median
polish after the harness normalisation) and the mokume methods of rollup_mokume.py under two normalisation views (harness,
mokume). Reference rows (not part of the shared BH family): the tool's own protein intensity where the run has one (MaxQuant
`LFQ intensity` of proteinGroups.txt, next to the evidence.txt of the arm).

Per roll-up x arm: proteins quantified (own, and inside the intersection family of the shared arms: matched depth), protein
log2-ratio error against the design truth per species (mean |err|, RMSE, bias, slope of observed on expected log2FC), median
CV, DE true / false positives at BH 5 % (own and intersection family) and true positives at pooled FDP 1 / 5 / 10 %."""
import argparse
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

from . import downstream as DS
from . import metrics as M
from . import readers as R
from . import rollup as RU
from . import rollup_mokume as RM
from . import suite as SU
from .de import de_all, de_table
from .design import Design
from .fasta import map_precursors


def read_mq_lfq(evidence_path, runs, ent_sp, min_peptides=2):
    """MaxQuant proteinGroups.txt (next to evidence.txt): `LFQ intensity <run>` -> log2 protein x run, indexed by the FASTA entry
    of the group (all ids of the group must map to entries of one species; reverse / contaminant groups and groups with fewer than
    `min_peptides` peptides are dropped)."""
    pg = pd.read_csv(os.path.join(os.path.dirname(evidence_path), "proteinGroups.txt"), sep="\t", low_memory=False)
    tok = {}
    for name in ent_sp:
        for t in name.split("|"):
            tok.setdefault(t, name)
        tok.setdefault(name, name)
    norm = {r.replace("-", "_").lower(): r for r in runs}
    cols = {c[len("LFQ intensity "):]: c for c in pg.columns if c.startswith("LFQ intensity ")}
    cmap = {norm[k.replace("-", "_").lower()]: c for k, c in cols.items() if k.replace("-", "_").lower() in norm}
    bad = pg.get("Reverse", pd.Series("", index=pg.index)).fillna("").eq("+") | pg.get("Potential contaminant", pd.Series("", index=pg.index)).fillna("").eq("+")
    rows, idx, drop = [], [], {"unmapped": 0, "mixed": 0, "few_peptides": 0, "dup": 0}
    for i, r in pg[~bad].iterrows():
        if float(r["Peptides"]) < min_peptides:
            drop["few_peptides"] += 1
            continue
        ents = [tok.get(x.strip()) for x in str(r["Protein IDs"]).split(";") if not x.startswith(("CON__", "REV__"))]
        if not ents or any(e is None for e in ents):
            drop["unmapped"] += 1
            continue
        if len({ent_sp[e] for e in ents}) != 1:
            drop["mixed"] += 1
            continue
        idx.append(ents[0])
        rows.append([r[cmap[x]] if x in cmap else np.nan for x in runs])
    P = pd.DataFrame(rows, index=idx, columns=runs).astype(float)
    P = P.where(P > 0)
    dup = P.index.duplicated(keep="first")
    drop["dup"] = int(dup.sum())
    print(f"mq-lfq {evidence_path}: groups {len(P)} dropped {drop}, runs mapped {len(cmap)}/{len(runs)}", file=sys.stderr)
    return np.log2(P[~dup])


def ratio_errors(P, design, prot_sp, classes):
    """Pooled over the condition pairs: per species, n, mean |err|, RMSE, bias (mean observed - expected log2FC) and the slope of
    observed on expected log2FC (species whose expected change varies)."""
    runs = list(P.columns)
    cond = {r: design.cond[r] for r in runs if r in design.cond}
    coll = {s: [] for s in classes}
    for i, lo in enumerate(design.levels):
        for hi in design.levels[i + 1:]:
            r1 = [r for r in runs if cond.get(r) == lo]
            r2 = [r for r in runs if cond.get(r) == hi]
            if len(r1) < 2 or len(r2) < 2:
                continue
            T = de_table(P, r1, r2, prot_sp, lambda s, lo=lo, hi=hi: design.expected_fc(s, lo, hi))
            for s in classes:
                m = (T.sp == s).values & np.isfinite(T.ex.values) & np.isfinite(T.fc.values)
                if m.any():
                    coll[s].append((T.ex.values[m], T.fc.values[m]))
    out = {}
    for s, L in coll.items():
        if not L:
            continue
        ex = np.concatenate([a for a, _ in L])
        fc = np.concatenate([b for _, b in L])
        e = fc - ex
        o = {"n": int(len(e)), "mean_abs": float(np.mean(np.abs(e))), "rmse": float(np.sqrt(np.mean(e ** 2))), "bias": float(np.mean(e))}
        if np.nanstd(ex) > 1e-9 and len(ex) > 10:
            o["slope"] = float(np.polyfit(ex, fc, 1)[0])
        out[s] = o
    return out


def cv_median(P, design):
    sd = M.pooled_sd(P, design.cond)
    return {"n": int(sd.notna().sum()), "med_cv_pct": float(M.cv_pct(sd).median()) if sd.notna().any() else None}


def score(PS, fam_arms, design, prot_sp, classes):
    """PS: arm -> log2 protein x run table (rows with >= 2 precursors, or the reference rows). The shared family is the
    intersection of the proteins of fam_arms; reference rows are scored on the same family."""
    inter = set.intersection(*[set(PS[n].index[PS[n].notna().any(axis=1)]) for n in fam_arms])
    union = set().union(*[set(PS[n].index) for n in fam_arms])
    out = {"family": {"arms": fam_arms, "intersection": len(inter), "union": len(union)}, "arms": {}}
    for n, P in PS.items():
        P = P[P.notna().any(axis=1)]
        Pi = P[P.index.isin(inter)]
        de, _, comp, fam = de_all(P, design, prot_sp, classes, inter, union)
        out["arms"][n] = {"proteins": int(len(P)), "proteins_inter": int(len(Pi)),
                          "de_own": {k: de[k] for k in ("TP", "FP", "TP5", "FP5", "tested")},
                          "de_inter": {k: fam["inter"][k] for k in ("TP", "FP", "TP5", "FP5", "tested")},
                          "tp_at_fdp_inter": fam["inter"]["tp_at_fdp"], "tp_at_fdp_own": fam["own"]["tp_at_fdp"],
                          "auc_inter": fam["inter"]["auc_mean"],
                          "err_inter": ratio_errors(Pi, design, prot_sp, classes), "err_own": ratio_errors(P, design, prot_sp, classes),
                          "cv_inter": cv_median(Pi, design), "cv_own": cv_median(P, design)}
    return out


def arm_specs(D, M_, seed):
    """(argv fragment, kinds) of the arms of one dataset and seed: andes main (default rule) and main_I, the tools of the manifest
    (MaxQuant, FragPipe)."""
    base = M_["_base"]
    fin = SU._p(base, D["andes"]["main"][seed]["finalize"]) if D.get("andes", {}).get("main", {}).get(seed) else None
    argv = []
    if fin:
        argv += ["--andes", f"andes:default:{fin}", "--andes", f"andes_I:ident:{fin}"]
    tools = {k: v for k, v in (D.get("tools") or {}).items() if k in ("maxquant", "fragpipe")}
    argv += SU.tool_args({"tools": tools}, base)
    return argv


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--seed", default="42")
    ap.add_argument("--out", required=True)
    ap.add_argument("--methods", default=",".join(RM.METHODS))
    ap.add_argument("--norms", default=",".join(RM.NORMS))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--held-out-eval", action="store_true")
    a = ap.parse_args(argv)
    warnings.filterwarnings("ignore")
    t0 = time.time()
    Mf = SU.load_manifest(a.manifest)
    D = Mf["datasets"][a.dataset]
    if D.get("held_out") and not a.held_out_eval:
        raise SystemExit(f"{a.dataset} is held out: refused without --held-out-eval")
    base = Mf["_base"]
    args = ["--name", D.get("label", a.dataset), "--design", SU._p(base, D["design"]), "--fasta", SU._p(base, D["fasta"]),
            "--classes", D["classes"], "--out", "/dev/null", "--amounts", SU._p(base, D["amounts"]), "--null", D["null"]]
    args += arm_specs(D, Mf, a.seed)
    A = DS.parse_args(args)
    classes = tuple(c for c in D["classes"].split(",") if c)
    design = Design(A.design, A.amounts, A.null, classes, A.low_max)
    runs = design.runs
    tabs, kinds = DS.load_arms(A, runs)
    allkeys = sorted(set().union(*[set(v.key) for v in tabs.values()]))
    sp, prot, ent_sp = map_precursors(allkeys, A.fasta, classes)
    cnt = {n: raw.assign(p=raw.key.map(prot)).dropna(subset=["p"]).groupby("p").key.nunique() for n, raw in tabs.items()}
    primary = next((n for n, k in kinds.items() if k["kind"] == "andes"), None)
    fam_arms = ([primary] if primary else []) + [n for n, k in kinds.items() if k["kind"] in ("maxquant", "fragpipe")]
    refs = {}
    for spec in A.maxquant:
        name, path = spec.split(":", 1)
        if os.path.exists(os.path.join(os.path.dirname(path), "proteinGroups.txt")):
            refs[name + "_LFQ"] = read_mq_lfq(path, runs, ent_sp)
    res = {"dataset": a.dataset, "label": D.get("label"), "seed": a.seed, "held_out": bool(D.get("held_out")), "runs": len(runs),
           "versions": RM.versions(), "arms": list(tabs), "reference_rows": list(refs), "roll_ups": {}}
    print("arms", {k: len(v) for k, v in tabs.items()}, "refs", list(refs), file=sys.stderr)

    def done(label, PS):
        res["roll_ups"][label] = score(PS, fam_arms, design, ent_sp, classes)
        print(label, round(time.time() - t0), "s", file=sys.stderr)

    PS = {}
    for n, raw in tabs.items():
        d, _, _ = RU.normalise(raw, runs)
        Pm = RU.rollup(d, prot, runs)
        PS[n] = Pm[cnt[n].reindex(Pm.index).fillna(0).values >= 2]
    PS.update(refs)
    done("harness|polish", PS)
    for norm in a.norms.split(","):
        T = {n: RM.peptide_table(raw, prot, runs, norm) for n, raw in tabs.items()}
        for m in a.methods.split(","):
            try:
                PS = {n: RM.protein_table(T[n], m, runs, a.threads) for n in tabs}
            except Exception as e:     # a method that fails is recorded, not hidden
                res["roll_ups"][f"{norm}|{m}"] = {"error": repr(e)[:400]}
                print(norm, m, "FAILED", repr(e)[:200], file=sys.stderr)
                continue
            PS.update(refs)
            done(f"{norm}|{m}", PS)
    res["elapsed_s"] = time.time() - t0
    with open(a.out, "w") as f:
        json.dump(res, f, indent=1, default=float)
    print("wrote", a.out, file=sys.stderr)


if __name__ == "__main__":
    main()
