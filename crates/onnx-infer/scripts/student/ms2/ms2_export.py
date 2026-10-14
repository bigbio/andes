#!/usr/bin/env python3
"""Write an andes-ms2 weight file and the Rust parity reference for it.

Either a trained run (--model-dir with ``model.pt`` and ``train.json``) or a
randomly initialised model of a given size (--random, used for the committed
test fixture). Every parameter of a random model, including the zero-initialised
FiLM output layer, is drawn from N(0, --random-std) so all paths are exercised.

The reference TSV has one row per (peptide, conditions): ``format`` (proforma |
andes), ``peptide``, ``charge``, ``nce`` (NA = unknown), ``activation``,
``instrument``, ``values`` (site-major, ``N_ION`` per site, comma-separated,
clipped at 0 with invalid channels 0; ``NA`` when the peptide cannot be encoded)
and ``raw`` (the unclipped network output, same layout).
Peptides come from --shards (sampled, stratified so phospho / oxidation / TMT /
acetyl / pyro-Glu forms are present) and a built-in list of edge cases, each also
written in andes PIN notation.

Usage::

  python ms2_export.py --model-dir RUN --shards SHARDS --out-weights andes_ms2.ams2 --out-reference ref.tsv
  python ms2_export.py --random --ch 8 --emb 4 --dilations 1 2 --n-heads 2 --cond-emb 2 --film-hidden 4 \\
      --out-weights tiny.ams2 --out-reference tiny_reference.tsv
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys

import numpy as np
import pyarrow.parquet as pq
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import ACTIVATIONS, parse_proforma  # noqa: E402
from ms2_model import Ms2Cnn, export_weights, predict  # noqa: E402
from mod_composition import MASSES  # noqa: E402

EDGE = [
    "PEPTIDEK", "AIK", "GK", "LIVIDLLK",
    "[UNIMOD:1]-M[UNIMOD:35]PEPS[UNIMOD:21]TK",
    "[UNIMOD:737]-K[UNIMOD:737]PEPS[UNIMOD:21]K[UNIMOD:737]",
    "[UNIMOD:28]-QPEPT[UNIMOD:21]M[UNIMOD:35]K",
    "AC[UNIMOD:4]DEFGHIKLM[UNIMOD:35]NPQRSTVWY",
    "N[UNIMOD:7]GT[UNIMOD:21]Y[UNIMOD:21]S[UNIMOD:21]R",
    "PEPTIDEK-[UNIMOD:2]",
    "ELVISLIVESK" * 4 + "R",
]
UNSUPPORTED = [("proforma", "PEPXIDEK"), ("proforma", "A" * 51), ("andes", "K.PEPXK.D"), ("andes", "K.PEP+1.23456K.D"),
               ("proforma", "K")]
INSTRUMENT_NAMES = ["", "Q Exactive HF", "Orbitrap Exploris 480", "Orbitrap Fusion Lumos", "LTQ Orbitrap Velos",
                    "Orbitrap Astral", "LTQ"]


def to_andes(form: str) -> str | None:
    p = parse_proforma(form)
    if p is None:
        return None
    seq, mods = p
    if any(pos == len(seq) + 1 for pos, _ in mods):
        return None  # andes PIN notation has no C-terminal slot
    per = {}
    for pos, uid in mods:
        per.setdefault(max(1, pos), []).append((pos != 0, MASSES[uid]))
    out = []
    for i, c in enumerate(seq, 1):
        out.append(c)
        for _, m in sorted(per.get(i, []), key=lambda t: t[0]):
            out.append(f"{m:+.5f}")
    return "-." + "".join(out) + ".-"


def andes_to_proforma(s: str) -> str | None:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from rt_model import parse_andes

    from ms2_common import to_proforma

    p = parse_andes(s)
    return to_proforma(*p) if p else None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model-dir")
    p.add_argument("--random", action="store_true")
    p.add_argument("--random-std", type=float, default=0.3)
    p.add_argument("--emb", type=int, default=32)
    p.add_argument("--ch", type=int, default=96)
    p.add_argument("--dilations", type=int, nargs="+", default=[1, 2, 4, 8, 1, 2])
    p.add_argument("--n-heads", type=int, default=4)
    p.add_argument("--cond-emb", type=int, default=8)
    p.add_argument("--film-hidden", type=int, default=64)
    p.add_argument("--shards", nargs="*", default=[])
    p.add_argument("--n", type=int, default=500, help="peptidoforms sampled from --shards")
    p.add_argument("--out-weights", required=True)
    p.add_argument("--out-reference", required=True)
    p.add_argument("--seed", type=int, default=7)
    a = p.parse_args()
    torch.manual_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    if a.model_dir:
        cfg = json.load(open(os.path.join(a.model_dir, "train.json")))["config"]
        m = Ms2Cnn(cfg["emb"], cfg["ch"], cfg["dilations"], cfg["n_heads"], cfg["cond_emb"], cfg["film_hidden"])
        m.load_state_dict(torch.load(os.path.join(a.model_dir, "model.pt"), map_location="cpu"))
    elif a.random:
        m = Ms2Cnn(a.emb, a.ch, a.dilations, a.n_heads, a.cond_emb, a.film_hidden)
        with torch.no_grad():
            for prm in m.parameters():
                prm.normal_(0, a.random_std)
    else:
        sys.exit("--model-dir or --random")
    export_weights(m, a.out_weights)

    forms = list(EDGE)
    files = []
    for s in a.shards:
        files += sorted(glob.glob(os.path.join(s, "**", "*.parquet"), recursive=True)) if os.path.isdir(s) else [s]
    if files:
        pool = sorted(set().union(*(set(pq.read_table(f, columns=["peptidoform"])["peptidoform"].to_pylist())
                                    for f in files)))
        pick = set(rng.choice(pool, min(a.n, len(pool)), replace=False).tolist())
        for tag in ["UNIMOD:21]", "UNIMOD:35]", "UNIMOD:737]", "UNIMOD:28]", "UNIMOD:1]-"]:
            sub = [f for f in pool if tag in f]
            if sub:
                pick |= set(rng.choice(sub, min(50, len(sub)), replace=False).tolist())
        forms += sorted(pick)

    rows = []
    for k, f in enumerate(forms):
        z = 1 + k % 4
        nce = float("nan") if k % 3 == 0 else 20.0 + (k % 7) * 2.5
        act = ACTIVATIONS[k % len(ACTIVATIONS)]
        instr = INSTRUMENT_NAMES[k % len(INSTRUMENT_NAMES)]
        rows.append(("proforma", f, f, z, nce, act, instr))
        an = to_andes(f)
        if an is not None:
            rows.append(("andes", an, andes_to_proforma(an), z, nce, act, instr))
    for fmt, s in UNSUPPORTED:
        rows.append((fmt, s, None if fmt == "andes" and andes_to_proforma(s) is None else
                     (andes_to_proforma(s) if fmt == "andes" else s), 2, 30.0, "HCD-FTMS", ""))
    reqs = [(pf or "", z, nce, act, instr) for _, _, pf, z, nce, act, instr in rows]
    preds = predict(m, reqs)
    raws = predict(m, reqs, raw=True)
    n_ok = 0
    with open(a.out_reference, "w") as fh:
        fh.write("format\tpeptide\tcharge\tnce\tactivation\tinstrument\tvalues\traw\n")
        for (fmt, s, _, z, nce, act, instr), y, r in zip(rows, preds, raws):
            vals = "NA" if y is None else ",".join(f"{v:.7g}" for v in y.ravel())
            rawv = "NA" if r is None else ",".join(f"{v:.7g}" for v in r.ravel())
            n_ok += y is not None
            nce_s = "NA" if not math.isfinite(nce) else f"{nce:g}"
            fh.write(f"{fmt}\t{s}\t{z}\t{nce_s}\t{act}\t{instr}\t{vals}\t{rawv}\n")
    print(f"{a.out_weights}: {sum(p.numel() for p in m.parameters())} params; "
          f"{a.out_reference}: {len(rows)} rows, {n_ok} encodable")


if __name__ == "__main__":
    main()
