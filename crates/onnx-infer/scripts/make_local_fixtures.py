#!/usr/bin/env python3
"""Build local ONNX/AGBD fixtures that do not need PeptDeep or torch.

Produces under fixtures/:
  - tree_leq.agbd / tree_leq.onnx     (1-feature split: x0 <= 0.5 -> 1.0 else 3.0)
  - peak_gbdt.agbd / peak_gbdt.onnx   (3-feature multi-tree, andes peak-like)
  - intensity_mlp.onnx                (tiny MatMul+ReLU; tract-runnable)
  - identity.onnx                     (no nodes; input x is also the output)

Usage:
  python scripts/make_local_fixtures.py
"""
from __future__ import annotations

import json
import struct
from pathlib import Path

import numpy as np
from onnx import TensorProto, checker, helper, numpy_helper, save

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "fixtures"


def write_agbd(
    path: Path,
    nfeat: int,
    flags: int,
    trees: list,
    iso_x: list | None = None,
    iso_y: list | None = None,
) -> None:
    """Write the layout `GbdtPeakModel::to_bytes` writes, calibration footer included."""
    iso_x = list(iso_x or [])
    iso_y = list(iso_y or [])
    if len(iso_x) != len(iso_y):
        raise ValueError("iso_x and iso_y must have the same length")
    buf = bytearray()
    buf += b"AGBD"
    buf += struct.pack("<IIII", 1, nfeat, flags, len(trees))
    for t in trees:
        n = len(t["feature"])
        buf += struct.pack("<I", n)
        buf += struct.pack(f"<{n}i", *t["feature"])
        buf += struct.pack(f"<{n}f", *t["threshold"])
        buf += struct.pack(f"<{n}i", *t["left"])
        buf += struct.pack(f"<{n}i", *t["right"])
        buf += struct.pack(f"<{n}f", *t["value"])
        buf += bytes(t["default_left"])
    buf += struct.pack("<I", len(iso_x))
    buf += struct.pack(f"<{len(iso_x)}f", *iso_x)
    buf += struct.pack(f"<{len(iso_y)}f", *iso_y)
    path.write_bytes(buf)


def _export_onnx(nfeat: int, trees: list, onnx_path: Path) -> None:
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    from export_agbd_to_onnx import trees_to_onnx as conv

    save(conv(nfeat, trees), onnx_path)


def make_tree_leq():
    """One tree, one split; no calibration points. `tree_leq.json` is hand-written."""
    trees = [
        dict(
            feature=[0, -1, -1],
            threshold=[0.5, 0.0, 0.0],
            left=[1, 0, 0],
            right=[2, 0, 0],
            value=[0.0, 1.0, 3.0],
            default_left=[1, 0, 0],
        )
    ]
    agbd = FIX / "tree_leq.agbd"
    onnx_path = FIX / "tree_leq.onnx"
    write_agbd(agbd, 1, 0, trees)
    _export_onnx(1, trees, onnx_path)
    print("wrote", agbd, onnx_path)


def make_peak_gbdt():
    """3-feature ensemble mirroring scoring::gbdt_eval packed_evaluation shape."""
    trees = [
        dict(
            feature=[0, 1, -1, -1, -1],
            threshold=[0.5, -2.0, 0.0, 0.0, 0.0],
            left=[1, 3, -1, -1, -1],
            right=[2, 4, -1, -1, -1],
            value=[0.0, 0.0, 0.37, -1.25, 2.5e-3],
            default_left=[1, 0, 0, 0, 0],
        ),
        dict(
            feature=[2, 0, -1, -1, -1],
            threshold=[1e-3, 0.5, 0.0, 0.0, 0.0],
            left=[1, 3, -1, -1, -1],
            right=[2, 4, -1, -1, -1],
            value=[0.0, 0.0, 0.11, 0.42, -0.08],
            default_left=[0, 1, 0, 0, 0],
        ),
        dict(
            feature=[-1],
            threshold=[0.0],
            left=[-1],
            right=[-1],
            value=[0.125],
            default_left=[0],
        ),
    ]
    nfeat = 3
    agbd = FIX / "peak_gbdt.agbd"
    onnx_path = FIX / "peak_gbdt.onnx"
    # Sigmoid flag plus a three-point isotonic map, so the footer is exercised by
    # the native decoder test; the ONNX graph is the raw tree sum either way.
    write_agbd(
        agbd, nfeat, 1, trees, iso_x=[0.0, 0.5, 1.0], iso_y=[0.05, 0.4, 0.95]
    )
    _export_onnx(nfeat, trees, onnx_path)
    man = {
        "id": "peak-gbdt-fixture",
        "description": "3-feature multi-tree AGBD→ONNX (andes peak-GBDT shape)",
        "onnx": "peak_gbdt.onnx",
        "inputs": [{"name": "features", "shape": [None, 3]}],
        "outputs": [{"name": "score", "shape": [None, 1]}],
        "task": "tree_ensemble",
        "provenance": {
            "origin": "onnx-infer fixture (andes AGBD layout)",
            "license": "Apache-2.0",
        },
        "extras": {
            "n_features": 3,
            "n_trees": 3,
            "agbd": "peak_gbdt.agbd",
            "output": "raw_tree_sum",
            "apply_sigmoid": True,
            "iso_x": [0.0, 0.5, 1.0],
            "iso_y": [0.05, 0.4, 0.95],
            "consumer": "quantms-rescoring",
        },
    }
    (FIX / "peak_gbdt.json").write_text(json.dumps(man, indent=2) + "\n")
    print("wrote", agbd, onnx_path)


def make_intensity_mlp():
    """Tiny [B,4] → ReLU → [B,1] intensity score; tract-friendly opset."""
    # y = ReLU(x @ W + b) with fixed weights
    W = np.array([[0.5], [-0.25], [0.75], [0.1]], dtype=np.float32)  # 4x1
    b = np.array([0.05], dtype=np.float32)
    w_init = numpy_helper.from_array(W, name="W")
    b_init = numpy_helper.from_array(b, name="B")
    nodes = [
        helper.make_node("MatMul", ["features", "W"], ["xw"]),
        helper.make_node("Add", ["xw", "B"], ["pre"]),
        helper.make_node("Relu", ["pre"], ["score"]),
    ]
    graph = helper.make_graph(
        nodes,
        "intensity_mlp",
        # Symbolic batch dim (name "batch") — tract rejects bare None/0 dims.
        [helper.make_tensor_value_info("features", TensorProto.FLOAT, ["batch", 4])],
        [helper.make_tensor_value_info("score", TensorProto.FLOAT, ["batch", 1])],
        initializer=[w_init, b_init],
    )
    model = helper.make_model(
        graph,
        producer_name="onnx-infer/make_local_fixtures",
        opset_imports=[helper.make_opsetid("", 17)],
        ir_version=8,
    )
    checker.check_model(model)
    out = FIX / "intensity_mlp.onnx"
    save(model, out)
    man = {
        "id": "intensity-mlp-fixture",
        "description": "Toy intensity MLP (MatMul+Add+ReLU); tract CI + rescoring batch smoke",
        "onnx": "intensity_mlp.onnx",
        "inputs": [{"name": "features", "shape": [None, 4]}],
        "outputs": [{"name": "score", "shape": [None, 1]}],
        "task": "rescoring",
        "provenance": {"origin": "onnx-infer fixture", "license": "Apache-2.0"},
        "extras": {"n_features": 4, "consumer": "quantms-rescoring"},
    }
    (FIX / "intensity_mlp.json").write_text(json.dumps(man, indent=2) + "\n")
    # parity with numpy
    x = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [-1, -1, -1, -1]], dtype=np.float32)
    y = np.maximum(0, x @ W + b)
    print("wrote", out, "smoke", y.reshape(-1))


def make_identity():
    """A graph without nodes whose input `x` is also its output."""
    x = helper.make_tensor_value_info("x", TensorProto.FLOAT, ["n"])
    graph = helper.make_graph([], "identity", [x], [x])
    model = helper.make_model(
        graph,
        producer_name="onnx-infer/make_local_fixtures",
        opset_imports=[helper.make_opsetid("", 17)],
        ir_version=8,
    )
    checker.check_model(model)
    out = FIX / "identity.onnx"
    save(model, out)
    print("wrote", out)


def main():
    FIX.mkdir(exist_ok=True)
    make_tree_leq()
    make_peak_gbdt()
    make_intensity_mlp()
    make_identity()


if __name__ == "__main__":
    main()
