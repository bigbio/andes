#!/usr/bin/env python3
"""PeptDeep reference data for the `rescore_pin` example.

Writes two files next to the exported ONNX models (or into --out-dir):

* ``peptdeep_mods.json`` -- the alphabase ``MOD_TO_FEATURE`` vectors for the
  modifications named on the command line, plus the MS2 fragment-type column
  order and the charge / NCE input scaling that peptdeep applies before the
  network (the exported graph does not include them).
* ``peptdeep_reference.tsv`` -- RT and MS2 predictions from peptdeep's own
  ``predict`` path for a handful of peptides, with the andes PIN spelling of
  each, so a Rust consumer can check its peptide parsing and tensor packing
  against the reference implementation.

It also runs the exported ONNX graphs under onnxruntime with the same packing
the Rust example uses and reports the largest deviation from peptdeep.

Usage::

  python peptdeep_reference.py --pretrained-dir .../generic --onnx-dir .../onnx

``--mods`` defaults to every modification andes encodes for AlphaPeptDeep
(``ANDES_MODS``); a modification missing from ``peptdeep_mods.json`` makes every
peptide that carries it a missing row in ``andes --predict ms2``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import onnxruntime as ort
import pandas as pd

from peptdeep.model.featurize import MOD_TO_FEATURE, get_batch_aa_indices
from peptdeep.model.ms2 import pDeepModel
from peptdeep.model.rt import AlphaRTModel

PEPTIDES = [
    # (sequence, mods, mod_sites, charge, andes PIN peptide). Sites follow
    # alphabase: 0 = N-terminus (including X^Any_N-term mods such as pyro-Glu),
    # 1..n = residue. The Rust parity check parses the andes string and requires
    # the same names and sites.
    ("PEPTIDEK", "", "", 2, "K.PEPTIDEK.A"),
    ("LGEYGFQNALIVR", "", "", 2, "K.LGEYGFQNALIVR.A"),
    ("AMCDEFGHIKLR", "Oxidation@M;Carbamidomethyl@C", "2;3", 3, "K.AM+15.99491C+57.02146DEFGHIKLR.A"),
    ("SAEDEMKLR", "Acetyl@Protein_N-term;Oxidation@M", "0;6", 2, "K.S+42.01056AEDEM+15.99491KLR.A"),
    ("VVVGAGGVGK", "", "", 1, "K.VVVGAGGVGK.A"),
    ("YCCQWNKR", "Carbamidomethyl@C;Carbamidomethyl@C", "2;3", 2, "K.YC+57.02146C+57.02146QWNKR.A"),
    ("QPEPTIDEK", "Gln->pyro-Glu@Q^Any_N-term", "0", 2, "K.Q-17.02655PEPTIDEK.A"),
    ("EAGHYLLR", "Glu->pyro-Glu@E^Any_N-term", "0", 2, "K.E-18.01056AGHYLLR.A"),
    ("KPEPTIDER", "TMT6plex@Any_N-term", "0", 2, "R.K+229.16293PEPTIDER.A"),
]


# The alphabase names of the modifications andes parses (onnx_infer::predictors::form::Mod).
ANDES_MODS = [
    "Carbamidomethyl@C", "Oxidation@M", "Acetyl@Protein_N-term", "Phospho@S", "Phospho@T", "Phospho@Y",
    "TMT6plex@K", "TMT6plex@Any_N-term", "Deamidated@N", "Deamidated@Q",
    "Gln->pyro-Glu@Q^Any_N-term", "Glu->pyro-Glu@E^Any_N-term",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pretrained-dir", type=Path, required=True)
    p.add_argument("--onnx-dir", type=Path, required=True)
    p.add_argument("--mods", nargs="+", default=ANDES_MODS)
    p.add_argument("--out-dir", type=Path, default=None, help="default: --onnx-dir")
    p.add_argument("--nce", type=float, default=30.0)
    p.add_argument("--instrument", default="Lumos")
    return p.parse_args()


def pack(seq, mods, sites, mod_dim):
    """Same packing as the Rust example: ascii-1 indices, zero pad at both termini."""
    aa = get_batch_aa_indices([seq]).astype(np.int64)
    mod_x = np.zeros((1, len(seq) + 2, mod_dim), dtype=np.float32)
    for m, s in zip(filter(None, mods.split(";")), filter(None, sites.split(";"))):
        s = int(s)
        mod_x[0, s if s >= 0 else len(seq) + 1] += MOD_TO_FEATURE[m]
    return aa, mod_x


def main():
    a = parse_args()
    ms2 = pDeepModel()
    ms2.load(str(a.pretrained_dir / "ms2.pth"))
    rt = AlphaRTModel()
    rt.load(str(a.pretrained_dir / "rt.pth"))
    frag_types = list(ms2.charged_frag_types)
    mod_dim = len(next(iter(MOD_TO_FEATURE.values())))

    from peptdeep.model.featurize import instrument_dict, unknown_inst_index

    meta = {
        "mod_dim": mod_dim,
        "mods": {m: [float(v) for v in MOD_TO_FEATURE[m]] for m in a.mods},
        "ms2_frag_types": frag_types,
        "charge_factor": float(ms2.charge_factor),
        "nce_factor": float(ms2.NCE_factor),
        "instrument_index": dict(instrument_dict),
        "unknown_instrument_index": int(unknown_inst_index),
        "aa_index": "ascii(residue) - ascii('A') + 1, zero pad at N- and C-terminus (len+2)",
        "ms2_rows": "nAA-1 cleavage sites; row i = b(i+1) and y(nAA-i-1)",
    }
    out_dir = a.out_dir or a.onnx_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    missing = {m for _, mods, _, _, _ in PEPTIDES for m in filter(None, mods.split(";"))} - set(a.mods)
    if missing:
        raise SystemExit(f"--mods must include the panel's modifications; missing {sorted(missing)}")
    out_json = out_dir / "peptdeep_mods.json"
    out_json.write_text(json.dumps(meta, indent=1) + "\n")
    print("wrote", out_json, "frag types:", frag_types)

    df = pd.DataFrame(PEPTIDES, columns=["sequence", "mods", "mod_sites", "charge", "andes"])
    df["nAA"] = df.sequence.str.len()
    df["nce"] = a.nce
    df["instrument"] = a.instrument
    df = df.sort_values("nAA").reset_index(drop=True)
    rt_df = rt.predict(df.copy())

    s_rt = ort.InferenceSession(str(a.onnx_dir / "peptdeep_rt_dynamic.onnx"))
    s_ms2 = ort.InferenceSession(str(a.onnx_dir / "peptdeep_ms2_dynamic.onnx"))
    inst = instrument_dict.get(a.instrument.upper(), unknown_inst_index)
    rows, max_rt, max_ms2 = [], 0.0, 0.0
    for i, r in df.iterrows():
        aa, mod_x = pack(r.sequence, r.mods, r.mod_sites, mod_dim)
        o_rt = s_rt.run(None, {"input_sequences": aa, "mod_x": mod_x})[0].reshape(-1)[0]
        o_ms2 = s_ms2.run(
            None,
            {
                "aa_indices": aa,
                "mod_x": mod_x,
                "charges": np.array([[r.charge * ms2.charge_factor]], dtype=np.float32),
                "nce": np.array([[a.nce * ms2.NCE_factor]], dtype=np.float32),
                "instrument_indices": np.array([inst], dtype=np.int64),
            },
        )[0][0]
        o_ms2 = np.clip(o_ms2, 0, None)
        # peptdeep zeroes fragment charges above the precursor charge
        for j, ft in enumerate(frag_types):
            if int(ft.rsplit("z", 1)[1]) > r.charge:
                o_ms2[:, j] = 0
        ref_rt = rt_df.loc[i, "rt_pred"]
        max_rt = max(max_rt, abs(float(ref_rt) - float(o_rt)))
        fi = ms2.predict(df.iloc[[i]].copy())[frag_types]
        flat = ",".join(f"{v:.7g}" for v in fi.values.reshape(-1))
        rows.append((r.sequence, r.mods, r.mod_sites, r.charge, r.andes, float(ref_rt), float(o_rt), flat))
        max_ms2 = max(max_ms2, float(np.abs(fi.values - o_ms2).max()))
    ref_tsv = out_dir / "peptdeep_reference.tsv"
    pd.DataFrame(
        rows,
        columns=["sequence", "mods", "mod_sites", "charge", "andes", "rt_peptdeep", "rt_onnxruntime", "ms2_peptdeep"],
    ).to_csv(ref_tsv, sep="\t", index=False)
    print(f"max |rt peptdeep - onnxruntime| = {max_rt:.3g}")
    print(f"max |ms2 peptdeep - onnxruntime| = {max_ms2:.3g}")
    for r in rows:
        print(r[:7])


if __name__ == "__main__":
    main()
