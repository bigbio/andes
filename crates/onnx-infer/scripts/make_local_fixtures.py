#!/usr/bin/env python3
"""Build local ONNX/AGBD fixtures that do not need PeptDeep or torch.

Produces under fixtures/:
  - tree_leq.agbd                     (1-feature split; tree_leq.onnx is committed)
  - peak_gbdt.agbd / peak_gbdt.onnx   (3-feature multi-tree, andes peak-like)
  - peak_gbdt_calibrated.agbd         (peak_gbdt trees with a 3-point
                                       calibration footer; parser test input)
  - parity/agbd.json                  (rows + expected raw scores for both AGBD
                                       fixtures; checked against the native Rust
                                       decoder and against ONNX Runtime)
  - intensity_mlp.onnx                (tiny MatMul+ReLU; tract-runnable)
  - identity.onnx                     (no nodes; input x is also the output)
  - identity_named.onnx               (Identity node `copy_node`: x -> y; output-name check)

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


def write_agbd(path: Path, nfeat: int, flags: int, trees: list, iso_x=(), iso_y=()) -> None:
    """Write AGBD v1 exactly as `GbdtPeakModel::to_bytes` does.

    The blob ends with the calibration footer: `u32 n_iso`, then `n_iso` f32
    `iso_x` values and `n_iso` f32 `iso_y` values. `n_iso = 0` is an
    uncalibrated model.
    """
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


def walk_trees(trees: list, x: list) -> float:
    """Raw sum of leaf values; NaN or a missing feature follows `default_left`."""
    total = np.float32(0.0)
    for t in trees:
        node = 0
        while t["feature"][node] >= 0:
            feat = t["feature"][node]
            v = x[feat] if feat < len(x) else float("nan")
            if np.isnan(v):
                go_left = t["default_left"][node] == 1
            else:
                go_left = np.float32(v) <= np.float32(t["threshold"][node])
            node = t["left"][node] if go_left else t["right"][node]
        total = np.float32(total + np.float32(t["value"][node]))
    return float(total)


TREE_LEQ_TREES = [
    dict(
        feature=[0, -1, -1],
        threshold=[0.5, 0.0, 0.0],
        left=[1, 0, 0],
        right=[2, 0, 0],
        value=[0.0, 1.0, 3.0],
        default_left=[1, 0, 0],
    )
]
TREE_LEQ_ROWS = [[0.0], [0.5], [0.6], [float("nan")]]

PEAK_GBDT_ROWS = [
    [0.1, -3.0, 0.0],
    [0.9, 4.0, 2.0],
    [float("nan"), float("nan"), float("nan")],
    [0.5, -2.0, 1e-3],
]

PEAK_GBDT_TREES = [
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


def make_tree_leq():
    """1-feature tree: feature0 <= 0.5 -> 1.0 else 3.0 (NaN goes left)."""
    agbd = FIX / "tree_leq.agbd"
    write_agbd(agbd, 1, 0, TREE_LEQ_TREES)
    print("wrote", agbd)


def make_peak_gbdt():
    """3-feature ensemble mirroring scoring::gbdt_eval packed_evaluation shape."""
    trees = PEAK_GBDT_TREES
    nfeat = 3
    agbd = FIX / "peak_gbdt.agbd"
    onnx_path = FIX / "peak_gbdt.onnx"
    write_agbd(agbd, nfeat, 0, trees)
    # import sibling script
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    from export_agbd_to_onnx import trees_to_onnx as conv

    save(conv(nfeat, trees), onnx_path)
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
            "consumer": "quantms-rescoring",
        },
    }
    (FIX / "peak_gbdt.json").write_text(json.dumps(man, indent=2) + "\n")
    print("wrote", agbd, onnx_path)


PEAK_GBDT_ISO_X = [0.0, 0.5, 1.0]
PEAK_GBDT_ISO_Y = [0.1, 0.4, 0.9]


def make_peak_gbdt_calibrated():
    """peak_gbdt trees with a non-empty calibration footer (n_iso = 3)."""
    agbd = FIX / "peak_gbdt_calibrated.agbd"
    write_agbd(agbd, 3, 0, PEAK_GBDT_TREES, PEAK_GBDT_ISO_X, PEAK_GBDT_ISO_Y)
    print("wrote", agbd)


def make_agbd_parity():
    """Rows and expected raw scores shared by the Rust and Python parity tests.

    NaN feature values are written as JSON `null`.
    """

    def entry(agbd: str, onnx_name: str, trees: list, rows: list) -> dict:
        return {
            "agbd": agbd,
            "onnx": onnx_name,
            "rows": [[None if np.isnan(v) else v for v in r] for r in rows],
            "expected": [walk_trees(trees, r) for r in rows],
        }

    parity = [
        entry("tree_leq.agbd", "tree_leq.onnx", TREE_LEQ_TREES, TREE_LEQ_ROWS),
        entry("peak_gbdt.agbd", "peak_gbdt.onnx", PEAK_GBDT_TREES, PEAK_GBDT_ROWS),
    ]
    # Not in fixtures/ itself: every JSON there is a model manifest.
    out = FIX / "parity" / "agbd.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(parity, indent=2) + "\n")
    print("wrote", out)


def make_intensity_mlp():
    """Tiny [B,4] → ReLU → [B,1] intensity score; tract-friendly opset."""
    # y = ReLU(x @ W + b) with fixed weights
    W = np.array([[0.5], [-0.25], [0.75], [0.1]], dtype=np.float32)  # 4x1
    b = np.array([0.05], dtype=np.float32)
    w_init = numpy_helper.from_array(W, name="W")
    b_init = numpy_helper.from_array(b, name="B")
    # Node names differ from their output tensor names (as in torch exports, where
    # the output `rt_pred` is produced by a node `/Squeeze`), so backends are
    # tested to report tensor names.
    nodes = [
        helper.make_node("MatMul", ["features", "W"], ["xw"], name="/MatMul"),
        helper.make_node("Add", ["xw", "B"], ["pre"], name="/Add"),
        helper.make_node("Relu", ["pre"], ["score"], name="/Relu"),
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


def make_identity_named():
    """`Identity(x) -> y` in a node named `copy_node`.

    The node name differs from the graph output name, so a backend that reports
    node names instead of output tensor names fails the tract test on this file.
    """
    graph = helper.make_graph(
        [helper.make_node("Identity", ["x"], ["y"], name="copy_node")],
        "identity_named",
        [helper.make_tensor_value_info("x", TensorProto.FLOAT, ["batch", 1])],
        [helper.make_tensor_value_info("y", TensorProto.FLOAT, ["batch", 1])],
    )
    model = helper.make_model(
        graph,
        producer_name="onnx-infer/make_local_fixtures",
        opset_imports=[helper.make_opsetid("", 17)],
        ir_version=8,
    )
    checker.check_model(model)
    out = FIX / "identity_named.onnx"
    save(model, out)
    print("wrote", out)


def main():
    FIX.mkdir(exist_ok=True)
    make_tree_leq()
    make_peak_gbdt()
    make_peak_gbdt_calibrated()
    make_agbd_parity()
    make_intensity_mlp()
    make_identity()
    make_identity_named()


if __name__ == "__main__":
    main()
