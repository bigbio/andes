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

The SavedModels use the GPU-only ``CudnnRNN`` GRU kernel, so TensorFlow cannot
execute them on a CPU; tf2onnx still converts them (to ONNX ``GRU``). The
reference predictions therefore come from Koina's public server, which serves
these exact SavedModels: see ``prosit_reference.py``. Needs tensorflow and
tf2onnx (only for the export).

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

import onnx

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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--intensity-zip", type=Path, required=True)
    p.add_argument("--irt-zip", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    import tensorflow as tf

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sm_int = unzip(a.intensity_zip, td / "int")
        sm_irt = unzip(a.irt_zip, td / "irt")
        rename = {"peptide": "peptides_in", "sequence": "peptides_in", "charge": "precursor_charge_in",
                  "collision": "collision_energy_in"}
        convert(sm_int, a.out_dir / "prosit_intensity_hcd.onnx", rename, "intensities_raw", False)
        convert(sm_irt, a.out_dir / "prosit_irt.onnx", rename, "irt", True)
    for stem, out, src in [
        ("prosit_intensity_hcd", "intensities_raw", "Prosit_2020_intensity_HCD (Zenodo 7681452)"),
        ("prosit_irt", "irt", "Prosit_2019_irt (Zenodo 7671314)"),
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
            "extras": {"tensorflow": tf.__version__},
        }
        (a.out_dir / f"{stem}.json").write_text(json.dumps(manifest, indent=1) + "\n")
        print("wrote", a.out_dir / f"{stem}.onnx")


if __name__ == "__main__":
    main()
