#!/usr/bin/env python3
"""Export Prosit 2020 intensity (HCD) and Prosit 2019 iRT to ONNX for `onnx-infer`.

The weights are Koina's TensorFlow SavedModels (Apache-2.0):

  intensity: https://zenodo.org/record/7681452/files/model.savedmodel.zip
  iRT:       https://zenodo.org/record/7671314/files/model.savedmodel.zip

Each SavedModel's ``serving_default`` signature is converted with tf2onnx
(opset 17) and its inputs / outputs renamed to the `onnx_infer::prosit`
contract:

  intensity: peptides_in [B,30], precursor_charge_in [B,6], collision_energy_in [B,1]
             -> intensities_raw [B,174]  (29 positions x y+1..3, b+1..3; before masking)
  iRT:       peptides_in [B,30] -> irt [B]  (Prosit's normalised iRT, before Koina's
             rescaling, which is affine and therefore irrelevant after RT calibration)

It then writes ``prosit_reference.tsv`` with TensorFlow's own outputs on a
peptide panel (CAM-C / oxidised M, charges 1-4, NCE 25-35) for the Rust parity
check (``model_parity --model prosit``), and reports onnxruntime vs TensorFlow.
Needs tensorflow-cpu and tf2onnx (only for the export).

Usage::

  python export_prosit_to_onnx.py --intensity-zip int.zip --irt-zip irt.zip --out-dir models/prosit
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import onnx

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
    ("ACDEFGHIKLMNPQRSTVWY", "K.AC+57.02146DEFGHIKLMNPQRSTVWY.A", 3, 27.0),
    ("LLLLLLLLLLLLLLLLLLLLLLLLLLLLLK", "K.LLLLLLLLLLLLLLLLLLLLLLLLLLLLLK.A", 3, 30.0),
]


def encode(seq: str) -> np.ndarray:
    import re

    pat = re.compile("|".join(re.escape(k) for k in sorted(ALPHABET, key=len, reverse=True)))
    toks = pat.findall(seq)
    assert "".join(toks) == seq and len(toks) <= SEQ_LEN, seq
    out = np.zeros(SEQ_LEN, dtype=np.int64)
    out[: len(toks)] = [ALPHABET[t] for t in toks]
    return out


def unzip(z: Path, dst: Path) -> Path:
    with zipfile.ZipFile(z) as f:
        f.extractall(dst)
    return next(p.parent for p in dst.rglob("saved_model.pb"))


def convert(saved: Path, out: Path, rename: dict, out_name: str, squeeze: bool):
    subprocess.run(
        [sys.executable, "-m", "tf2onnx.convert", "--saved-model", str(saved), "--opset", "17",
         "--output", str(out)],
        check=True,
    )
    m = onnx.load(str(out))
    g = m.graph
    mapping = {}
    for vi in g.input:
        new = next(v for k, v in rename.items() if k in vi.name.lower())
        mapping[vi.name] = new
        vi.name = new
    assert len(g.output) == 1, [o.name for o in g.output]
    mapping[g.output[0].name] = out_name if not squeeze else out_name + "_2d"
    g.output[0].name = mapping[g.output[0].name]
    for node in g.node:
        node.input[:] = [mapping.get(x, x) for x in node.input]
        node.output[:] = [mapping.get(x, x) for x in node.output]
    if squeeze:
        axes = onnx.helper.make_tensor("squeeze_axes", onnx.TensorProto.INT64, [1], [1])
        g.initializer.append(axes)
        g.node.append(onnx.helper.make_node("Squeeze", [out_name + "_2d", "squeeze_axes"], [out_name]))
        old = g.output[0]
        g.output.remove(old)
        g.output.append(onnx.helper.make_tensor_value_info(out_name, onnx.TensorProto.FLOAT, ["B"]))
    onnx.checker.check_model(m)
    onnx.save(m, str(out))


def tf_run(saved: Path, feeds: dict):
    import tensorflow as tf

    fn = tf.saved_model.load(str(saved)).signatures["serving_default"]
    spec = fn.structured_input_signature[1]
    args = {}
    for name, s in spec.items():
        key = next(k for k in feeds if k in name.lower())
        args[name] = tf.constant(feeds[key].astype(s.dtype.as_numpy_dtype))
    out = fn(**args)
    assert len(out) == 1
    return next(iter(out.values())).numpy()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--intensity-zip", type=Path, required=True)
    p.add_argument("--irt-zip", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    import onnxruntime as ort
    import tensorflow as tf

    seqs = np.stack([encode(x[0]) for x in PANEL])
    charge = np.zeros((len(PANEL), 6), dtype=np.float32)
    for i, x in enumerate(PANEL):
        charge[i, x[2] - 1] = 1
    ce = np.array([[x[3] / 100] for x in PANEL], dtype=np.float32)
    feeds = {"peptide": seqs, "charge": charge, "collision": ce}
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sm_int = unzip(a.intensity_zip, td / "int")
        sm_irt = unzip(a.irt_zip, td / "irt")
        ref_int = tf_run(sm_int, feeds).reshape(len(PANEL), -1)
        ref_irt = tf_run(sm_irt, {"peptide": seqs, "sequence": seqs}).reshape(-1)
        rename = {"peptide": "peptides_in", "sequence": "peptides_in", "charge": "precursor_charge_in",
                  "collision": "collision_energy_in"}
        convert(sm_int, a.out_dir / "prosit_intensity_hcd.onnx", rename, "intensities_raw", False)
        convert(sm_irt, a.out_dir / "prosit_irt.onnx", rename, "irt", True)

    def ort_run(path, names):
        s = ort.InferenceSession(str(path))
        f = {}
        for i in s.get_inputs():
            arr = {"peptides_in": seqs, "precursor_charge_in": charge, "collision_energy_in": ce}[i.name]
            f[i.name] = arr.astype({"tensor(float)": np.float32, "tensor(int64)": np.int64,
                                    "tensor(int32)": np.int32}[i.type])
        return s.run(names, f)[0]

    o_int = ort_run(a.out_dir / "prosit_intensity_hcd.onnx", ["intensities_raw"]).reshape(len(PANEL), -1)
    o_irt = ort_run(a.out_dir / "prosit_irt.onnx", ["irt"]).reshape(-1)
    e_int, e_irt = float(np.abs(o_int - ref_int).max()), float(np.abs(o_irt - ref_irt).max())
    print(f"onnxruntime vs tensorflow: intensity max |d| {e_int:.3g}, iRT max |d| {e_irt:.3g}")
    for stem, out, err, src in [
        ("prosit_intensity_hcd", "intensities_raw", e_int, "Prosit_2020_intensity_HCD (Zenodo 7681452)"),
        ("prosit_irt", "irt", e_irt, "Prosit_2019_irt (Zenodo 7671314)"),
    ]:
        manifest = {
            "id": stem,
            "description": f"{src} via Koina SavedModel + tf2onnx",
            "onnx": f"{stem}.onnx",
            "task": "ms2_intensity" if "intensity" in stem else "retention_time",
            "outputs": [{"name": out}],
            "provenance": {"origin": f"Koina {src}", "license": "Apache-2.0",
                           "notice": "https://github.com/wilhelm-lab/koina",
                           "notes": "trained on ProteomeTools synthetic peptides"},
            "extras": {"max_abs_err_vs_tensorflow": err, "tensorflow": tf.__version__},
        }
        (a.out_dir / f"{stem}.json").write_text(json.dumps(manifest, indent=1) + "\n")
    fmt = lambda v: ",".join(f"{x:.9g}" for x in np.asarray(v).reshape(-1))
    with open(a.out_dir / "prosit_reference.tsv", "w") as f:
        f.write("proforma\tandes\tcharge\tnce\tintensities_raw\tirt\n")
        for i, (pf, andes, z, nce) in enumerate(PANEL):
            f.write(f"{pf}\t{andes}\t{z}\t{nce}\t{fmt(ref_int[i])}\t{ref_irt[i]:.9g}\n")
    print("wrote", a.out_dir / "prosit_reference.tsv")


if __name__ == "__main__":
    main()
