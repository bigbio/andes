#!/usr/bin/env python3
"""Export the DeepLC default RT model (PyTorch, deeplc >= 4.1) to ONNX for `onnx-infer`.

DeepLC 4.x ships ``multitask_flexcnn_model.pt``: a convolutional trunk over
per-residue atom composition and residue identity, shared by a low-rank head
with one output per LC setup (6,543 setups). ``deeplc.predict`` without
calibration reports the ``PXD005573_mcp`` setup (the DeepLC 1.x-3.x training
data); that one head is exported:

  inputs  ``x_atom``    [B, 60, 6]  f32  per-position atom counts (C H N O S P)
          ``x_global``  [B, 67]     f32  global vector incl. terminal composition
          ``x_one_hot`` [B, 60, 20] f32  one-hot residues (DeepLC order), padding rows 0
  output  ``rt``        [B]         f32  minutes on the PXD005573_mcp gradient

The graph is the model's own forward with two ops rewritten into equivalents
that export to plain ONNX and run under tract: the residue embedding becomes
a one-hot matmul (the padding row of the embedding is applied where a
position is empty), the ``scatter_add`` residue counts become a sum of the
one-hot matrix over positions, and the masked max pool subtracts 1e30 at
padding positions instead of filling -inf (tract cannot type that ``Where``). ``scripts/deeplc_reference.py`` checks the result against
``deeplc.predict``.

Usage::

  python export_deeplc_to_onnx.py --out-dir models/deeplc
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from torch import nn

import deeplc
from deeplc import _model_ops
from deeplc.core import DEFAULT_MODEL, DEFAULT_TASK_NAME


class SingleHead(nn.Module):
    def __init__(self, model, task: int):
        super().__init__()
        self.enc = model.encoder
        head = model.head
        self.proj = head.proj
        self.register_buffer("emb", head.embedding[task].detach().clone())
        self.register_buffer("scale", head.scale[task].detach().clone())
        self.register_buffer("shift", head.shift[task].detach().clone())
        pad = self.enc.PAD_INDEX
        self.register_buffer("w_res", self.enc.embed.weight[:pad].detach().clone())
        self.register_buffer("w_pad", self.enc.embed.weight[pad].detach().clone())

    def forward(self, x_atom, x_global, x_one_hot):
        enc = self.enc
        occupied = x_one_hot.sum(dim=2, keepdim=True)  # [B, L, 1], 1 = residue, 0 = padding
        emb = x_one_hot @ self.w_res + (1.0 - occupied) * self.w_pad  # [B, L, E]
        atom = x_atom.transpose(1, 2)
        if enc.stem is not None:
            atom = enc.stem(atom)
        hidden = torch.cat([atom, emb.transpose(1, 2)], dim=1)
        for block in enc.blocks:
            hidden = block(hidden)
        valid = occupied.transpose(1, 2)  # [B, 1, L], 1.0 / 0.0
        hidden = hidden * valid
        summed = hidden.sum(dim=2)
        # Masked max over residues: padding positions are pushed far below any
        # activation. Every encoded peptide has residues, so the all-padding
        # case the model maps to 0 never arises.
        maxed = (hidden - (1.0 - valid) * 1e30).amax(dim=2)
        pooled = enc.pool_norm(torch.cat([summed, maxed], dim=1))
        dense = torch.cat([x_global, x_one_hot.sum(dim=1)], dim=1)
        trunk = enc.net(torch.cat([enc.norm(dense), pooled], dim=1))
        return (self.proj(trunk) @ self.emb) * self.scale + self.shift


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--task", default=DEFAULT_TASK_NAME)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)

    model = _model_ops.load_model(DEFAULT_MODEL, device="cpu").eval()
    spec = model.feature_spec
    assert spec["global_dim"] == 67 and spec["padding_length"] == 60, spec
    assert not spec["legacy_positional_deltas"] and spec["add_terminal_composition"], spec
    task = list(model.task_names).index(a.task)
    wrapped = SingleHead(model, task).eval()

    L, B = spec["padding_length"], 3
    rng = np.random.default_rng(0)
    x_atom = torch.tensor(rng.integers(0, 6, (B, L, 6)), dtype=torch.float32)
    x_global = torch.tensor(rng.integers(0, 30, (B, 67)), dtype=torch.float32)
    idx = rng.integers(0, 20, (B, L))
    x_hc = torch.zeros(B, L, 20)
    for b, n in enumerate([7, 20, 60]):
        x_hc[b, torch.arange(n), torch.tensor(idx[b, :n])] = 1.0
        x_atom[b, n:] = 0
    with torch.no_grad():
        ref = model(x_atom, torch.zeros(B, 0, 6), x_global, x_hc, torch.tensor([task]))[:, 0]
        got = wrapped(x_atom, x_global, x_hc)
    err = float((ref - got).abs().max())
    print(f"wrapper vs model forward on random inputs: max |d| = {err:.3g} min")

    onnx_path = a.out_dir / "deeplc_rt.onnx"
    torch.onnx.export(
        wrapped,
        (x_atom, x_global, x_hc),
        str(onnx_path),
        input_names=["x_atom", "x_global", "x_one_hot"],
        output_names=["rt"],
        dynamic_axes={"x_atom": {0: "B"}, "x_global": {0: "B"}, "x_one_hot": {0: "B"}, "rt": {0: "B"}},
        opset_version=17,
        dynamo=False,
    )
    import onnxruntime as ort

    sess = ort.InferenceSession(str(onnx_path))
    o = sess.run(None, {"x_atom": x_atom.numpy(), "x_global": x_global.numpy(), "x_one_hot": x_hc.numpy()})[0]
    err_ort = float(np.abs(o - ref.numpy()).max())
    print(f"onnxruntime vs model forward: max |d| = {err_ort:.3g} min")
    manifest = {
        "id": "deeplc_rt",
        "description": f"DeepLC {deeplc.__version__} multitask_flexcnn_model, head {a.task} (minutes)",
        "onnx": onnx_path.name,
        "task": "retention_time",
        "inputs": [
            {"name": "x_atom", "shape": [None, L, 6]},
            {"name": "x_global", "shape": [None, 67]},
            {"name": "x_one_hot", "shape": [None, L, 20]},
        ],
        "outputs": [{"name": "rt", "shape": [None]}],
        "provenance": {
            "origin": f"deeplc {deeplc.__version__} package_data/models/{Path(DEFAULT_MODEL).name}",
            "license": "Apache-2.0",
            "notice": "https://github.com/compomics/DeepLC",
            "notes": f"{model.n_tasks} LC setups; exported head {a.task} (index {task}); "
            f"feature_spec {json.dumps(spec)}",
        },
        "extras": {"task": a.task, "task_index": task, "padding_length": L, "max_abs_err_vs_torch": err_ort},
    }
    (a.out_dir / "deeplc_rt.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print("wrote", onnx_path, f"{onnx_path.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
