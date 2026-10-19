#!/usr/bin/env python3
"""DeepLC reference predictions for the onnx-infer parity check.

Runs ``deeplc.predict`` (default model, uncalibrated) on a peptide panel with
the modifications andes searches and writes ``deeplc_reference.tsv``: per
peptide the andes-style sequence with mass deltas, DeepLC's encoded inputs
(``x_atom`` / ``x_global`` / ``x_one_hot``, flattened) and its prediction. The
Rust side re-encodes the andes string and runs the exported graph with tract
(``model_parity --model deeplc``). The exported ONNX graph is also run here
under onnxruntime on DeepLC's own encoding.

Usage::

  python deeplc_reference.py --onnx-dir models/deeplc --out deeplc_reference.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import onnxruntime as ort
from psm_utils import PSM, Peptidoform, PSMList

import deeplc
from deeplc import _model_ops
from deeplc.core import DEFAULT_MODEL
from deeplc.data import DeepLCDataset

# (ProForma for DeepLC, andes PIN peptide string)
PANEL = [
    ("PEPTIDEK", "K.PEPTIDEK.A"),
    ("LGEYGFQNALIVR", "K.LGEYGFQNALIVR.Y"),
    ("AM[Oxidation]C[Carbamidomethyl]DEFGHIKLR", "K.AM+15.99491C+57.02146DEFGHIKLR.A"),
    ("[Acetyl]-SAEDEM[Oxidation]KLR", "-.S+42.01056AEDEM+15.99491KLR.A"),
    ("[Acetyl]-M[Oxidation]LRK", "-.M+42.01056+15.99491LRK.A"),
    ("VVVGAGGVGK", "K.VVVGAGGVGK.S"),
    ("YC[Carbamidomethyl]C[Carbamidomethyl]QWNKR", "K.YC+57.02146C+57.02146QWNKR.A"),
    ("[TMT6plex]-PEPTK[TMT6plex]IDEK[TMT6plex]", "K.P+229.16293EPTK+229.16293IDEK+229.16293.A"),
    ("[TMT6plex]-K[TMT6plex]PEPTIDER", "R.K+229.16293+229.16293PEPTIDER.A"),
    ("Q[Gln->pyro-Glu]PEPS[Phospho]TIDEN[Deamidated]K", "K.Q-17.02655PEPS+79.96633TIDEN+0.98402K.A"),
    ("E[Glu->pyro-Glu]AGHY[Phospho]LLR", "K.E-18.01056AGHY+79.96633LLR.A"),
    ("ACDEFGHIKLMNPQRSTVWY", "K.ACDEFGHIKLMNPQRSTVWY.A"),
    ("LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLK", "K.LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLK.A"),
    ("ASDK", "K.ASDK.A"),
]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--onnx-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    model = _model_ops.load_model(DEFAULT_MODEL, device="cpu")
    spec = model.feature_spec
    psms = PSMList(psm_list=[PSM(peptidoform=Peptidoform(x), spectrum_id=i) for i, (x, _) in enumerate(PANEL)])
    pred = deeplc.predict(psms)
    ds = DeepLCDataset.from_psm_list(
        psms,
        add_terminal_composition=spec["add_terminal_composition"],
        padding_length=spec["padding_length"],
        legacy_positional_deltas=spec["legacy_positional_deltas"],
        include_rolling_sum=False,
    )
    x_atom, _, x_global, x_hc = ds.encode_batch(range(len(PANEL)))
    sess = ort.InferenceSession(str(a.onnx_dir / "deeplc_rt.onnx"))
    o = sess.run(None, {"x_atom": x_atom.numpy(), "x_global": x_global.numpy(), "x_one_hot": x_hc.numpy()})[0]
    d = np.abs(o - pred)
    print(f"onnxruntime vs deeplc.predict: max |d| {d.max():.3g} min, max rel {(d / np.abs(pred)).max():.3g}")
    fmt = lambda v: ",".join(f"{x:.9g}" for x in np.asarray(v).reshape(-1))
    with open(a.out, "w") as f:
        f.write("proforma\tandes\tx_atom\tx_global\tx_one_hot\trt_deeplc\n")
        for i, (pf, andes) in enumerate(PANEL):
            f.write("\t".join([pf, andes, fmt(x_atom[i]), fmt(x_global[i]), fmt(x_hc[i]), f"{pred[i]:.9g}"]) + "\n")
    print(f"wrote {len(PANEL)} rows to {a.out}")
    for (pf, _), r, g in zip(PANEL, pred, o):
        print(f"{pf:60s} deeplc {r:9.4f} onnxruntime {g:9.4f}")


if __name__ == "__main__":
    main()
