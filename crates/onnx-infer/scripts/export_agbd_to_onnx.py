#!/usr/bin/env python3
"""Export an andes AGBD v1 GBDT blob to ONNX TreeEnsembleRegressor.

The native AGBD evaluator in `scoring::gbdt_eval` remains the in-search hot path.
This export is for interoperability with quantms-rescoring / OpenMS / ort, and
for offline parity checks against the Rust evaluator.

The ONNX graph computes the raw ensemble output (the sum of the tree leaf values),
which is `GbdtPeakModel::predict_value`. It does not apply the sigmoid flag or the
isotonic calibration footer that `predict_proba` / `predict_logit` apply on top of
that sum; both are written to the manifest `extras` (`apply_sigmoid`, `iso_x`,
`iso_y`) so a consumer that needs the calibrated probability can apply them:
  p = sigmoid(raw) if apply_sigmoid else raw
  s = piecewise-linear interpolation of (iso_x, iso_y) at p, clamped to the end
      points (identity when iso_x is empty)

Usage:
  python export_agbd_to_onnx.py --agbd path/to/model.agbd --out model.onnx

Requires: numpy, onnx
"""
from __future__ import annotations

import argparse
import io
import json
import struct
from pathlib import Path

import numpy as np
from onnx import TensorProto, checker, helper, save


FLAG_SIGMOID = 1


def _read_exact(c: io.BytesIO, n: int, what: str) -> bytes:
    data = c.read(n)
    if len(data) != n:
        raise ValueError(f"truncated AGBD blob: {what} needs {n} bytes, got {len(data)}")
    return data


def parse_agbd(blob: bytes):
    """Decode an AGBD v1 blob as written by `GbdtPeakModel::to_bytes`.

    Returns `(n_features, flags, trees, calibration)` where `calibration` is
    `{"iso_x": [...], "iso_y": [...]}` from the footer that follows the trees
    (`n_iso` u32, then `iso_x` f32 x n_iso, then `iso_y` f32 x n_iso).
    """
    c = io.BytesIO(blob)
    magic = _read_exact(c, 4, "magic")
    if magic != b"AGBD":
        raise ValueError(f"bad magic {magic!r}")
    ver, nfeat, flags, ntree = struct.unpack("<IIII", _read_exact(c, 16, "header"))
    if ver != 1:
        raise ValueError(f"unsupported AGBD version {ver}")
    trees = []
    for t in range(ntree):
        n = struct.unpack("<I", _read_exact(c, 4, f"tree {t} node count"))[0]
        feat = list(struct.unpack(f"<{n}i", _read_exact(c, 4 * n, f"tree {t} features")))
        thr = list(struct.unpack(f"<{n}f", _read_exact(c, 4 * n, f"tree {t} thresholds")))
        left = list(struct.unpack(f"<{n}i", _read_exact(c, 4 * n, f"tree {t} left children")))
        right = list(struct.unpack(f"<{n}i", _read_exact(c, 4 * n, f"tree {t} right children")))
        val = list(struct.unpack(f"<{n}f", _read_exact(c, 4 * n, f"tree {t} values")))
        dl = list(_read_exact(c, n, f"tree {t} default-left flags"))
        trees.append(
            dict(
                feature=feat,
                threshold=thr,
                left=left,
                right=right,
                value=val,
                default_left=dl,
            )
        )
    n_iso = struct.unpack("<I", _read_exact(c, 4, "calibration point count"))[0]
    iso_x = list(struct.unpack(f"<{n_iso}f", _read_exact(c, 4 * n_iso, "calibration iso_x")))
    iso_y = list(struct.unpack(f"<{n_iso}f", _read_exact(c, 4 * n_iso, "calibration iso_y")))
    trailing = len(blob) - c.tell()
    if trailing:
        raise ValueError(
            f"AGBD blob has {trailing} trailing bytes after {ntree} trees and "
            f"{n_iso} calibration points"
        )
    return nfeat, flags, trees, dict(iso_x=iso_x, iso_y=iso_y)


def trees_to_onnx(nfeat: int, trees: list, input_name: str = "features", output_name: str = "score"):
    nodes_treeids = []
    nodes_nodeids = []
    nodes_featureids = []
    nodes_modes = []
    nodes_values = []
    nodes_truenodeids = []
    nodes_falsenodeids = []
    nodes_missing = []
    target_treeids = []
    target_nodeids = []
    target_ids = []
    target_weights = []

    for ti, t in enumerate(trees):
        n = len(t["feature"])
        for ni in range(n):
            nodes_treeids.append(ti)
            nodes_nodeids.append(ni)
            feat = t["feature"][ni]
            if feat < 0:
                nodes_featureids.append(0)
                nodes_modes.append("LEAF")
                nodes_values.append(0.0)
                nodes_truenodeids.append(0)
                nodes_falsenodeids.append(0)
                nodes_missing.append(0)
                target_treeids.append(ti)
                target_nodeids.append(ni)
                target_ids.append(0)
                target_weights.append(float(t["value"][ni]))
            else:
                nodes_featureids.append(int(feat))
                # AGBD: go left when v <= threshold (or NaN and default_left)
                nodes_modes.append("BRANCH_LEQ")
                nodes_values.append(float(t["threshold"][ni]))
                nodes_truenodeids.append(int(t["left"][ni]))
                nodes_falsenodeids.append(int(t["right"][ni]))
                nodes_missing.append(1 if t["default_left"][ni] else 0)

    node = helper.make_node(
        "TreeEnsembleRegressor",
        inputs=[input_name],
        outputs=[output_name],
        domain="ai.onnx.ml",
        n_targets=1,
        aggregate_function="SUM",
        post_transform="NONE",
        nodes_treeids=nodes_treeids,
        nodes_nodeids=nodes_nodeids,
        nodes_featureids=nodes_featureids,
        nodes_modes=nodes_modes,
        nodes_values=nodes_values,
        nodes_truenodeids=nodes_truenodeids,
        nodes_falsenodeids=nodes_falsenodeids,
        nodes_missing_value_tracks_true=nodes_missing,
        target_treeids=target_treeids,
        target_nodeids=target_nodeids,
        target_ids=target_ids,
        target_weights=target_weights,
    )
    graph = helper.make_graph(
        [node],
        "agbd_tree_ensemble",
        [helper.make_tensor_value_info(input_name, TensorProto.FLOAT, ["batch", nfeat])],
        [helper.make_tensor_value_info(output_name, TensorProto.FLOAT, ["batch", 1])],
    )
    model = helper.make_model(
        graph,
        producer_name="onnx-infer/export_agbd_to_onnx",
        opset_imports=[
            helper.make_opsetid("", 17),
            helper.make_opsetid("ai.onnx.ml", 3),
        ],
        ir_version=8,
    )
    checker.check_model(model)
    return model


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--agbd", type=Path, required=True, help="Input AGBD blob")
    ap.add_argument("--out", type=Path, required=True, help="Output .onnx path")
    ap.add_argument("--manifest", type=Path, default=None, help="Optional sidecar JSON")
    ap.add_argument("--id", default="agbd-tree", help="Manifest id")
    args = ap.parse_args()

    blob = args.agbd.read_bytes()
    nfeat, flags, trees, calib = parse_agbd(blob)
    apply_sigmoid = bool(flags & FLAG_SIGMOID)
    model = trees_to_onnx(nfeat, trees)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    save(model, args.out)
    print(f"Wrote {args.out} ({len(trees)} trees, {nfeat} features, flags={flags})")
    if apply_sigmoid or calib["iso_x"]:
        print(
            f"note: the graph outputs the raw tree sum; sigmoid={apply_sigmoid} and "
            f"{len(calib['iso_x'])} isotonic points are recorded in the manifest extras, "
            "not applied in the graph"
        )

    man_path = args.manifest or args.out.with_suffix(".json")
    man = {
        "id": args.id,
        "description": f"AGBD→ONNX TreeEnsemble from {args.agbd.name}",
        "onnx": args.out.name,
        "inputs": [{"name": "features", "shape": [None, nfeat]}],
        "outputs": [{"name": "score", "shape": [None, 1]}],
        "task": "tree_ensemble",
        "provenance": {
            "origin": "andes AGBD v1",
            "license": "Apache-2.0",
            "notes": f"flags={flags}; output is the raw tree sum; prefer native AGBD in-search",
        },
        "extras": {
            "n_features": nfeat,
            "n_trees": len(trees),
            "agbd_flags": flags,
            "output": "raw_tree_sum",
            "apply_sigmoid": apply_sigmoid,
            "iso_x": calib["iso_x"],
            "iso_y": calib["iso_y"],
        },
    }
    man_path.write_text(json.dumps(man, indent=2) + "\n")
    print(f"Wrote {man_path}")

    # Optional parity smoke with onnxruntime
    try:
        import onnxruntime as ort

        sess = ort.InferenceSession(str(args.out), providers=["CPUExecutionProvider"])
        x = np.zeros((2, nfeat), dtype=np.float32)
        y = sess.run(None, {"features": x})[0]
        print("ort smoke:", y.reshape(-1)[:2])
    except Exception as e:
        print("ort smoke skipped:", e)


if __name__ == "__main__":
    main()
