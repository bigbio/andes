#!/usr/bin/env python3
"""Write the Rust parity reference for a trained andes-rt model.

Samples peptidoforms from the training table (stratified so phospho, TMT,
pyro-Glu, deamidation and acetyl forms are present), writes each in ProForma
and in andes PIN notation, adds unsupported cases (unknown residue, too long),
and records the PyTorch index (NaN when unsupported).

Usage::

  python rt_export.py --model-dir DIR --table rt_table.parquet --out rt_reference.tsv
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pyarrow.parquet as pq
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod_composition import MASSES  # noqa: E402
from rt_model import RtCnn, parse_proforma, predict  # noqa: E402


def to_andes(form: str) -> str:
    seq, mods = parse_proforma(form)
    per = {}
    for pos, uid in mods:
        per.setdefault(max(1, min(pos, len(seq))), []).append((pos != 0, MASSES[uid]))
    out = []
    for i, c in enumerate(seq, 1):
        out.append(c)
        # N-terminal deltas first, as andes prints them
        for _, m in sorted(per.get(i, []), key=lambda t: t[0]):
            out.append(f"{m:+.5f}")
    return "-." + "".join(out) + ".-"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model-dir", required=True)
    p.add_argument("--table", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--n", type=int, default=2000)
    a = p.parse_args()
    m = RtCnn()
    m.load_state_dict(torch.load(os.path.join(a.model_dir, "model.pt"), map_location="cpu"))
    forms = pq.read_table(a.table, columns=["peptidoform"])["peptidoform"].unique().to_pylist()
    rng = np.random.default_rng(7)
    pick = set(rng.choice(forms, a.n, replace=False).tolist())
    for tag in ["UNIMOD:21]", "UNIMOD:737]", "UNIMOD:28]", "UNIMOD:27]", "UNIMOD:7]", "UNIMOD:1]-"]:
        sub = [f for f in forms if tag in f]
        pick |= set(rng.choice(sub, min(100, len(sub)), replace=False).tolist()) if sub else set()
    pick = sorted(pick)
    rows = [("proforma", f) for f in pick] + [("andes", to_andes(f)) for f in pick]
    rows += [("proforma", "PEPXIDEK"), ("proforma", "A" * 51), ("andes", "K.PEPXK.D"), ("andes", "K.PEP+1.23456K.D")]
    vals = np.full(len(rows), np.nan)
    for fmt in ["proforma", "andes"]:
        idx = [i for i, r in enumerate(rows) if r[0] == fmt]
        vals[idx] = predict(m, [rows[i][1] for i in idx], fmt=fmt)
    pf = vals[: len(pick)]
    an = vals[len(pick): 2 * len(pick)]
    both = np.isfinite(pf) & np.isfinite(an)
    print(f"{len(rows)} rows; proforma vs andes-notation max diff {np.max(np.abs(pf[both] - an[both])):.2e} "
          f"on {both.sum()} forms; andes unparsed {np.sum(np.isfinite(pf) & ~np.isfinite(an))}")
    with open(a.out, "w") as f:
        f.write("format\tpeptide\tindex\n")
        for (fmt, s), v in zip(rows, vals):
            f.write(f"{fmt}\t{s}\t{'NA' if not np.isfinite(v) else repr(float(v))}\n")


if __name__ == "__main__":
    main()
