#!/usr/bin/env python3
"""Prosit reference predictions from Koina for the onnx-infer parity check.

Koina (https://koina.wilhelmlab.org) serves the same SavedModels the exporter
converts, with its own pre- and post-processing. This script asks Koina for
``Prosit_2020_intensity_HCD`` and ``Prosit_2019_irt`` predictions on a panel
(CAM-C / oxidised M, charges 1-4, NCE 25-35), runs the exported ONNX graphs
under onnxruntime with Koina's encoding, and writes ``prosit_reference.tsv``
(Koina intensities: 174 values, max-normalised, -1 = impossible ion; Koina
iRT = raw * 43.39373 + 56.35363441).

Usage::

  python prosit_reference.py --onnx-dir models/prosit --out prosit_reference.tsv
"""
from __future__ import annotations

import argparse
import json
import re
import urllib.request
from pathlib import Path

import numpy as np

KOINA = "https://koina.wilhelmlab.org/v2/models"
IRT_SCALE, IRT_SHIFT = 43.39373, 56.35363441
SEQ_LEN = 30
# Koina Prosit_Preprocess_peptide/1/sequence_conversion.py
ALPHABET = {
    "A": 1, "C": 2, "D": 3, "E": 4, "F": 5, "G": 6, "H": 7, "I": 8, "K": 9, "L": 10,
    "M": 11, "N": 12, "P": 13, "Q": 14, "R": 15, "S": 16, "T": 17, "V": 18, "W": 19,
    "Y": 20, "M[UNIMOD:35]": 21, "C[UNIMOD:4]": 2,
}
# (Koina ProForma, andes PIN string, charge, NCE)
PANEL = [
    ("PEPTIDEK", "K.PEPTIDEK.A", 2, 30.0),
    ("LGEYGFQNALIVR", "K.LGEYGFQNALIVR.Y", 2, 28.0),
    ("AM[UNIMOD:35]C[UNIMOD:4]DEFGHIKLR", "K.AM+15.99491C+57.02146DEFGHIKLR.A", 3, 30.0),
    ("SAEDEM[UNIMOD:35]KLR", "K.SAEDEM+15.99491KLR.A", 2, 25.0),
    ("VVVGAGGVGK", "K.VVVGAGGVGK.S", 1, 35.0),
    ("YC[UNIMOD:4]C[UNIMOD:4]QWNKR", "K.YC+57.02146C+57.02146QWNKR.A", 4, 30.0),
    ("AC[UNIMOD:4]DEFGHIKLMNPQRSTVWY", "K.AC+57.02146DEFGHIKLMNPQRSTVWY.A", 3, 27.0),
    ("LLLLLLLLLLLLLLLLLLLLLLLLLLLLLK", "K.LLLLLLLLLLLLLLLLLLLLLLLLLLLLLK.A", 3, 30.0),
]


def encode(seq: str) -> np.ndarray:
    pat = re.compile("|".join(re.escape(k) for k in sorted(ALPHABET, key=len, reverse=True)))
    toks = pat.findall(seq)
    assert "".join(toks) == seq and len(toks) <= SEQ_LEN, seq
    out = np.zeros(SEQ_LEN, dtype=np.int64)
    out[: len(toks)] = [ALPHABET[t] for t in toks]
    return out


def koina(model: str, inputs: list[dict]) -> dict:
    body = json.dumps({"id": "0", "inputs": inputs}).encode()
    req = urllib.request.Request(f"{KOINA}/{model}/infer", body, {"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        out = json.load(r)
    return {o["name"]: np.asarray(o["data"]).reshape(o["shape"]) for o in out["outputs"]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--onnx-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    n = len(PANEL)
    seqs = [x[0] for x in PANEL]
    k_int = koina("Prosit_2020_intensity_HCD", [
        {"name": "peptide_sequences", "shape": [n, 1], "datatype": "BYTES", "data": seqs},
        {"name": "precursor_charges", "shape": [n, 1], "datatype": "INT32", "data": [x[2] for x in PANEL]},
        {"name": "collision_energies", "shape": [n, 1], "datatype": "FP32", "data": [x[3] for x in PANEL]},
    ])["intensities"].astype(np.float64)
    k_irt = koina("Prosit_2019_irt", [
        {"name": "peptide_sequences", "shape": [n, 1], "datatype": "BYTES", "data": seqs},
    ])["irt"].reshape(-1).astype(np.float64)

    import onnxruntime as ort

    enc = np.stack([encode(s) for s in seqs])
    charge = np.zeros((n, 6), dtype=np.float32)
    for i, x in enumerate(PANEL):
        charge[i, x[2] - 1] = 1
    ce = np.array([[x[3] / 100] for x in PANEL], dtype=np.float32)

    def run(path):
        s = ort.InferenceSession(str(path))
        feeds = {}
        for i in s.get_inputs():
            arr = {"peptides_in": enc, "precursor_charge_in": charge, "collision_energy_in": ce}[i.name]
            feeds[i.name] = arr.astype({"tensor(float)": np.float32, "tensor(int64)": np.int64,
                                        "tensor(int32)": np.int32}[i.type])
        return s.run(None, feeds)[0]

    raw = run(a.onnx_dir / "prosit_intensity_hcd.onnx").reshape(n, -1)
    irt = run(a.onnx_dir / "prosit_irt.onnx").reshape(-1) * IRT_SCALE + IRT_SHIFT
    valid = k_int >= 0
    norm = np.where(raw < 0, np.finfo(np.float32).eps, raw)
    norm = np.where(valid, norm, np.nan)
    norm = norm / np.nanmax(norm, axis=1, keepdims=True)
    d_int = np.nanmax(np.abs(norm - k_int)[valid])
    print(f"onnxruntime vs Koina: intensity max |d| {d_int:.3g} ({valid.sum()} values), "
          f"iRT max |d| {np.abs(irt - k_irt).max():.3g}")
    fmt = lambda v: ",".join(f"{x:.9g}" for x in np.asarray(v).reshape(-1))
    with open(a.out, "w") as f:
        f.write("proforma\tandes\tcharge\tnce\tintensities_koina\tirt_koina\n")
        for i, (pf, andes, z, nce) in enumerate(PANEL):
            f.write(f"{pf}\t{andes}\t{z}\t{nce}\t{fmt(k_int[i])}\t{k_irt[i]:.9g}\n")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
