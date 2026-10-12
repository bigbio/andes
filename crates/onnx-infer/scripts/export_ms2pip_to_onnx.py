#!/usr/bin/env python3
"""Export MS²PIP XGBoost models to ONNX for `onnx-infer`.

MS²PIP (ms2pip >= 4.0) predicts each ion series with one XGBoost regressor
over 139 per-cleavage-site features (``ms2rescore_rs.ms2pip_compute_features``).
Each ``.xgboost`` file becomes one ONNX graph:

  input  ``features``       [N, 139] f32  (one row per cleavage site)
  output ``log2_intensity`` [N]      f32  = log2(intensity + 0.001), unclipped

The trees are written as plain tensor ops, not ``ai.onnx.ml`` tree ensembles:
tract has no ``TreeEnsembleRegressor``, and its ``TreeEnsembleClassifier``
looks up a leaf's weight by walking the leaf table from its start, which costs
O(#leaves) per tree visit (hours for the 9 M-leaf HCD2021 y model). Instead
all trees share four flat node tables (feature, threshold, left, right) plus
a leaf-value table, indexed by a global node id. Starting from every tree's
root, one unrolled step per tree level does

  node = Where(x[feature[node]] < threshold[node], left[node], right[node])

and leaves point to themselves, so after ``max_depth`` steps every tree sits
on its leaf; the output is the sum of leaf values plus the XGBoost base score.
The features never contain NaN, so XGBoost's missing-value direction is not
encoded.

The rows of the y model are in cleavage-site order, like the features:
row ``i`` is y_(n-1-i). ms2pip itself reverses the y predictions into
y1..y(n-1) order after prediction.

Usage::

  python export_ms2pip_to_onnx.py --model-dir ~/.ms2pip --out-dir models/ms2pip \\
      --models HCD2021 CID TMT

Writes ``ms2pip_<MODEL>_<ion>.onnx`` + ``.json`` manifest per ion type, then
checks every graph against ``xgboost.Booster.predict`` on the features of a
peptide panel (onnxruntime).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import onnx
import xgboost as xgb
from onnx import TensorProto, helper

from ms2pip.constants import MODELS

N_FEATURES = 139
PANEL = [
    "PEPTIDEK/2",
    "LGEYGFQNALIVR/2",
    "AMCDEFGHIKLR/3",
    "SAEDEMKLR/2",
    "VVVGAGGVGK/1",
    "YCCQWNKR/4",
    "ACDEFGHIKLMNPQRSTVWY/3",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model-dir", type=Path, required=True, help="directory with .xgboost files")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--models", nargs="+", default=["HCD2021", "CID", "TMT"])
    return p.parse_args()


def booster_base_score(b: xgb.Booster) -> float:
    cfg = json.loads(b.save_config())
    return float(cfg["learner"]["learner_model_param"]["base_score"])


def tree_tables(b: xgb.Booster):
    """Global node tables from the per-tree JSON dump, one tree at a time (the
    HCD2021 y model has ~19 M nodes; whole-model Python objects do not fit).

    Returns (roots, feature, threshold, left, right, leaf_value, max_depth)
    where node ids are global: tree ``t``'s nodes occupy one contiguous block,
    renumbered densely (pruned XGBoost trees leave gaps in node ids)."""
    roots, feat, thr, left, right, val = [], [], [], [], [], []
    offset, max_depth = 0, 0
    for dump in b.get_dump(dump_format="json"):
        recs = []
        stack = [(json.loads(dump), 0)]
        while stack:
            n, d = stack.pop()
            recs.append((n, d))
            if "leaf" not in n:
                stack.extend((c, d + 1) for c in n["children"])
        ids = np.sort(np.fromiter((n["nodeid"] for n, _ in recs), dtype=np.int64))
        k = len(ids)
        g = lambda i: offset + int(np.searchsorted(ids, i))
        f = np.zeros(k, np.int32)
        t = np.zeros(k, np.float32)
        lf = np.zeros(k, np.int32)
        rt = np.zeros(k, np.int32)
        v = np.zeros(k, np.float32)
        for n, d in recs:
            j = g(n["nodeid"]) - offset
            max_depth = max(max_depth, d)
            if "leaf" in n:
                lf[j] = rt[j] = offset + j
                v[j] = n["leaf"]
            else:
                f[j] = int(n["split"].lstrip("f"))
                t[j] = n["split_condition"]
                lf[j] = g(n["yes"])
                rt[j] = g(n["no"])
        roots.append(g(0))
        for lst, arr in zip((feat, thr, left, right, val), (f, t, lf, rt, v)):
            lst.append(arr)
        offset += k
    cat = np.concatenate
    return np.asarray(roots, np.int32), cat(feat), cat(thr), cat(left), cat(right), cat(val), max_depth


def to_onnx(b: xgb.Booster, name: str) -> onnx.ModelProto:
    from onnx import numpy_helper

    roots, feat, thr, left, right, val, depth = tree_tables(b)
    n_trees = len(roots)
    inits = [
        numpy_helper.from_array(feat, "feature"),
        numpy_helper.from_array(thr, "threshold"),
        numpy_helper.from_array(left, "left"),
        numpy_helper.from_array(right, "right"),
        numpy_helper.from_array(val, "leaf_value"),
        numpy_helper.from_array(roots, "roots"),
        numpy_helper.from_array(np.asarray([0], np.int64), "zero1"),
        numpy_helper.from_array(np.asarray([1], np.int64), "one1"),
        numpy_helper.from_array(np.asarray([n_trees], np.int64), "n_trees"),
        numpy_helper.from_array(np.asarray([1], np.int64), "axis1"),
        numpy_helper.from_array(np.asarray(booster_base_score(b), np.float32), "base_score"),
    ]
    nodes = [
        helper.make_node("Shape", ["features"], ["shape"]),
        helper.make_node("Slice", ["shape", "zero1", "one1"], ["n_rows"]),
        helper.make_node("Concat", ["n_rows", "n_trees"], ["grid"], axis=0),
        helper.make_node("Expand", ["roots", "grid"], ["node_0"]),
    ]
    for d in range(depth):
        cur, nxt = f"node_{d}", f"node_{d + 1}"
        nodes += [
            helper.make_node("Gather", ["feature", cur], [f"f_{d}"]),
            helper.make_node("GatherElements", ["features", f"f_{d}"], [f"x_{d}"], axis=1),
            helper.make_node("Gather", ["threshold", cur], [f"t_{d}"]),
            helper.make_node("Less", [f"x_{d}", f"t_{d}"], [f"go_left_{d}"]),
            helper.make_node("Gather", ["left", cur], [f"l_{d}"]),
            helper.make_node("Gather", ["right", cur], [f"r_{d}"]),
            helper.make_node("Where", [f"go_left_{d}", f"l_{d}", f"r_{d}"], [nxt]),
        ]
    nodes += [
        helper.make_node("Gather", ["leaf_value", f"node_{depth}"], ["leaves"]),
        helper.make_node("ReduceSum", ["leaves", "axis1"], ["margin"], keepdims=0),
        helper.make_node("Add", ["margin", "base_score"], ["log2_intensity"]),
    ]
    graph = helper.make_graph(
        nodes,
        name,
        [helper.make_tensor_value_info("features", TensorProto.FLOAT, ["N", N_FEATURES])],
        [helper.make_tensor_value_info("log2_intensity", TensorProto.FLOAT, ["N"])],
        initializer=inits,
    )
    model = helper.make_model(
        graph,
        opset_imports=[helper.make_opsetid("", 17)],
        producer_name="onnx-infer export_ms2pip_to_onnx.py",
    )
    model.ir_version = 8
    print(f"{name}: {n_trees} trees, {len(feat)} nodes, max depth {depth}")
    return model


def panel_features():
    from ms2rescore_rs import ms2pip_compute_features

    feats = ms2pip_compute_features(PANEL)
    return np.concatenate([f.reshape(-1, N_FEATURES) for f in feats]).astype(np.float32)


def main():
    a = parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    import onnxruntime as ort
    import ms2pip

    x = panel_features()
    for model in a.models:
        files = MODELS[model]["xgboost_model_files"]
        for ion, fname in files.items():
            b = xgb.Booster()
            b.load_model(str(a.model_dir / fname))
            stem = f"ms2pip_{model}_{ion}"
            onnx_path = a.out_dir / f"{stem}.onnx"
            m = to_onnx(b, stem)
            onnx.save(m, str(onnx_path))
            ref = b.predict(xgb.DMatrix(x))
            sess = ort.InferenceSession(str(onnx_path))
            got = sess.run(["log2_intensity"], {"features": x})[0]
            err = float(np.abs(ref - got).max())
            manifest = {
                "id": stem,
                "description": f"MS²PIP {model} {ion}-ion XGBoost model as an ONNX tree-traversal graph "
                "(log2(intensity + 0.001), one row per cleavage site)",
                "onnx": onnx_path.name,
                "task": "ms2_intensity",
                "inputs": [{"name": "features", "shape": [None, N_FEATURES]}],
                "outputs": [{"name": "log2_intensity", "shape": [None]}],
                "provenance": {
                    "origin": f"ms2pip {ms2pip.__version__} {fname}",
                    "license": "Apache-2.0",
                    "notice": "https://github.com/compomics/ms2pip",
                    "notes": f"sha1 {MODELS[model]['model_hash'][fname]}",
                },
                "extras": {"ion": ion, "model": model, "max_abs_err_vs_xgboost": err},
            }
            (a.out_dir / f"{stem}.json").write_text(json.dumps(manifest, indent=1) + "\n")
            print(f"{stem}: {b.num_boosted_rounds()} trees, {onnx_path.stat().st_size / 1e6:.1f} MB, "
                  f"max |onnxruntime - xgboost| = {err:.3g}")


if __name__ == "__main__":
    main()
