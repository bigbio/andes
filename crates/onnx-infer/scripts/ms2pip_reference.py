#!/usr/bin/env python3
"""MS²PIP reference predictions for the onnx-infer parity check.

Runs ms2pip's own ``predict_batch`` (ms2rescore-rs features + XGBoost) on a
peptide panel with modifications and charges 1-4 and writes
``ms2pip_reference.tsv``: per peptide and model the plain sequence, charge,
the 139-value feature rows and the b / y log2 predictions in cleavage-site
order (row i = b(i+1), y(n-1-i)). The Rust side recomputes the features and
runs the exported ONNX graphs with tract (``model_parity --model ms2pip``).

Usage::

  python ms2pip_reference.py --model-dir ~/.ms2pip --out ms2pip_reference.tsv \\
      --models HCD2021 CID TMT
"""
from __future__ import annotations

import argparse

import numpy as np
from ms2rescore_rs import ms2pip_compute_features
from psm_utils import PSM, Peptidoform, PSMList

from ms2pip.core import predict_batch

PANEL = [
    "PEPTIDEK/2",
    "LGEYGFQNALIVR/2",
    "AM[Oxidation]C[Carbamidomethyl]DEFGHIKLR/3",
    "[Acetyl]-SAEDEM[Oxidation]KLR/2",
    "VVVGAGGVGK/1",
    "YC[Carbamidomethyl]C[Carbamidomethyl]QWNKR/4",
    "[TMT6plex]-PEPTK[TMT6plex]IDEK/2",
    "Q[Gln->pyro-Glu]PEPS[Phospho]TIDEN[Deamidated]K/3",
    "ACDEFGHIKLMNPQRSTVWY/3",
    "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLK/4",
]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--models", nargs="+", default=["HCD2021", "CID", "TMT"])
    a = p.parse_args()
    psms = PSMList(psm_list=[PSM(peptidoform=Peptidoform(x), spectrum_id=i) for i, x in enumerate(PANEL)])
    rows = []
    for model in a.models:
        res = predict_batch(psms, model=model, model_dir=a.model_dir, processes=1)
        for r in res:
            pf = r.psm.peptidoform
            seq, z = pf.sequence, pf.precursor_charge
            feats = ms2pip_compute_features([f"{seq}/{z}"])[0]
            b = np.asarray(r.predicted_intensity["b"])
            y = np.asarray(r.predicted_intensity["y"])[::-1]  # back to cleavage-site order
            fmt = lambda v: ",".join(f"{x:.9g}" for x in v)
            rows.append("\t".join([model, str(pf), seq, str(z), fmt(feats), fmt(b), fmt(y)]))
    with open(a.out, "w") as f:
        f.write("model\tproforma\tsequence\tcharge\tfeatures\tb_log2\ty_log2\n")
        f.write("\n".join(rows) + "\n")
    print(f"wrote {len(rows)} rows to {a.out}")


if __name__ == "__main__":
    main()
