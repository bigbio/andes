#!/usr/bin/env python3
"""Export PeptDeep RT / MS2 / CCS pretrained PyTorch models to ONNX.

Adapted from OpenMS `tools/scripts/export_peptdeep_models_to_onnx.py` so
quantms-rescoring and andes `onnx-infer` share one I/O contract:

  RT:  input_sequences, mod_x → rt_pred
  MS2: aa_indices, mod_x, charges, nce, instrument_indices → ms2_intensities
  CCS: aa_indices, mod_x, charges → ccs_pred

Usage (after `pip install peptdeep torch onnx onnxruntime`):

  python export_peptdeep_to_onnx.py \\
    --pretrained-dir ~/peptdeep/pretrained_models/pretrained_models_v3/generic \\
    --out-dir ./models/peptdeep

Does not run during CI; weights are large and licensed separately.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

import numpy as np
import onnx
import onnxruntime as ort
import torch

from peptdeep.model.rt import AlphaRTModel
from peptdeep.settings import model_const

DEFAULT_PRETRAINED_DIR = Path(
    "~/peptdeep/pretrained_models/pretrained_models_v3/generic"
).expanduser()
MOD_HIDDEN = len(model_const["mod_elements"])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pretrained-dir", type=Path, default=DEFAULT_PRETRAINED_DIR)
    p.add_argument("--out-dir", type=Path, default=Path("models/peptdeep"))
    p.add_argument(
        "--tasks",
        nargs="+",
        default=["rt", "ms2", "ccs"],
        choices=["rt", "ms2", "ccs"],
    )
    return p.parse_args()


def get_model_class(module_name, candidate_names):
    module = importlib.import_module(module_name)
    for name in candidate_names:
        if hasattr(module, name):
            return getattr(module, name)
    available = [n for n in dir(module) if "Model" in n or "model" in n]
    raise ImportError(
        f"Could not find any of {candidate_names} in {module_name}. "
        f"Available: {available}"
    )


def make_dummy_inputs(task, batch_size=1, seq_len=16):
    aa = torch.zeros((batch_size, seq_len), dtype=torch.long)
    mod_x = torch.zeros((batch_size, seq_len, MOD_HIDDEN), dtype=torch.float32)
    charges = torch.full((batch_size, 1), 2.0, dtype=torch.float32)
    nce = torch.full((batch_size, 1), 30.0, dtype=torch.float32)
    instrument_indices = torch.zeros((batch_size,), dtype=torch.long)
    if task == "rt":
        return aa, mod_x
    if task == "ccs":
        return aa, mod_x, charges
    if task == "ms2":
        return aa, mod_x, charges, nce, instrument_indices
    raise ValueError(task)


def warmup_and_load(model, pth_path, task):
    model.set_device("cpu")
    model.model.to(model.device)
    model.model.eval()
    with torch.no_grad():
        model.model(*make_dummy_inputs(task, seq_len=6))
    model.load(str(pth_path))
    model.model.eval()
    return model


def export_and_check(model, args, out_path, input_names, output_names, dynamic_axes):
    out_path = Path(out_path)
    torch.onnx.export(
        model.model,
        args,
        str(out_path),
        input_names=input_names,
        output_names=output_names,
        dynamic_axes=dynamic_axes,
        opset_version=17,
        do_constant_folding=True,
        dynamo=False,
        external_data=False,
    )
    onnx_model = onnx.load(str(out_path))
    onnx.checker.check_model(onnx_model)
    print(f"Wrote {out_path}")
    sess = ort.InferenceSession(str(out_path), providers=["CPUExecutionProvider"])
    for inp in sess.get_inputs():
        print(" ", inp.name, inp.type, inp.shape)
    return sess


def graph_io(onnx_path: Path):
    """Inputs and outputs of the exported graph as manifest tensor specs.

    Each shape is the one the graph declares: a fixed axis as its size, a
    dynamic (symbolic or unknown) axis as None. A tensor without a declared
    shape gets no "shape" entry. andes checks these against the loaded graph.
    """
    graph = onnx.load(str(onnx_path)).graph
    initializers = {t.name for t in graph.initializer}

    def spec(value):
        out = {"name": value.name}
        tensor = value.type.tensor_type
        if tensor.HasField("shape"):
            out["shape"] = [
                d.dim_value if d.HasField("dim_value") else None
                for d in tensor.shape.dim
            ]
        return out

    inputs = [spec(v) for v in graph.input if v.name not in initializers]
    return inputs, [spec(v) for v in graph.output]


def write_manifest(out_dir: Path, task: str, onnx_name: str):
    task_map = {"rt": "retention_time", "ms2": "ms2_intensity", "ccs": "ccs"}
    inputs, outputs = graph_io(out_dir / onnx_name)
    man = {
        "id": f"peptdeep-{task}-onnx",
        "description": f"PeptDeep {task} ONNX (quantms-rescoring contract)",
        "onnx": onnx_name,
        "inputs": inputs,
        "outputs": outputs,
        "task": task_map[task],
        "provenance": {
            "origin": "AlphaPeptDeep / peptdeep pretrained",
            "license": "see AlphaPeptDeep NOTICE",
            "notes": "Exported via onnx-infer/scripts/export_peptdeep_to_onnx.py",
        },
        "extras": {
            "mod_hidden": MOD_HIDDEN,
            "nce_default": 30.0,
            "consumer": "quantms-rescoring",
        },
    }
    path = out_dir / f"peptdeep_{task}_dynamic.json"
    path.write_text(json.dumps(man, indent=2) + "\n")
    print(f"Wrote {path}")


def export_rt(pretrained_dir, out_dir):
    model = warmup_and_load(AlphaRTModel(), pretrained_dir / "rt.pth", "rt")
    name = "peptdeep_rt_dynamic.onnx"
    export_and_check(
        model=model,
        args=make_dummy_inputs("rt"),
        out_path=out_dir / name,
        input_names=["input_sequences", "mod_x"],
        output_names=["rt_pred"],
        dynamic_axes={
            "input_sequences": {0: "batch_size", 1: "seq_length"},
            "mod_x": {0: "batch_size", 1: "seq_length"},
            "rt_pred": {0: "batch_size"},
        },
    )
    write_manifest(out_dir, "rt", name)


def export_ms2(pretrained_dir, out_dir):
    cls = get_model_class(
        "peptdeep.model.ms2", ["AlphaMS2Model", "pDeepModel", "MS2Model"]
    )
    model = warmup_and_load(cls(), pretrained_dir / "ms2.pth", "ms2")
    name = "peptdeep_ms2_dynamic.onnx"
    export_and_check(
        model=model,
        args=make_dummy_inputs("ms2"),
        out_path=out_dir / name,
        input_names=["aa_indices", "mod_x", "charges", "nce", "instrument_indices"],
        output_names=["ms2_intensities"],
        dynamic_axes={
            "aa_indices": {0: "batch_size", 1: "seq_length"},
            "mod_x": {0: "batch_size", 1: "seq_length"},
            "charges": {0: "batch_size"},
            "nce": {0: "batch_size"},
            "instrument_indices": {0: "batch_size"},
            "ms2_intensities": {0: "batch_size", 1: "frag_position"},
        },
    )
    write_manifest(out_dir, "ms2", name)


def export_ccs(pretrained_dir, out_dir):
    cls = get_model_class("peptdeep.model.ccs", ["AlphaCCSModel", "CCSModel"])
    model = warmup_and_load(cls(), pretrained_dir / "ccs.pth", "ccs")
    name = "peptdeep_ccs_dynamic.onnx"
    export_and_check(
        model=model,
        args=make_dummy_inputs("ccs"),
        out_path=out_dir / name,
        input_names=["aa_indices", "mod_x", "charges"],
        output_names=["ccs_pred"],
        dynamic_axes={
            "aa_indices": {0: "batch_size", 1: "seq_length"},
            "mod_x": {0: "batch_size", 1: "seq_length"},
            "charges": {0: "batch_size"},
            "ccs_pred": {0: "batch_size"},
        },
    )
    write_manifest(out_dir, "ccs", name)


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    exporters = {"rt": export_rt, "ms2": export_ms2, "ccs": export_ccs}
    for task in args.tasks:
        exporters[task](args.pretrained_dir, args.out_dir)
    print("mod_hidden =", MOD_HIDDEN)
    # Silence unused numpy import when ort path unused in some envs
    _ = np


if __name__ == "__main__":
    main()
