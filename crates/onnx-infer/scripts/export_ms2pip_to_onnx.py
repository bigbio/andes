#!/usr/bin/env python3
"""Export MS²PIP XGBoost models to ONNX for `onnx-infer`.

MS²PIP (ms2pip >= 4.0) predicts each ion series with one XGBoost regressor
over 139 per-cleavage-site features (``ms2rescore_rs.ms2pip_compute_features``).
Each ``.xgboost`` file becomes one ONNX graph:

  input  ``features`` [N, 139] f32  (one row per cleavage site)
  output ``label``    [N]      i64  (ignored)
  output ``scores``   [N, 3]   f32  column 0 = log2(intensity + 0.001)

tract (the pure-Rust backend) implements ``ai.onnx.ml.TreeEnsembleClassifier``
but not ``TreeEnsembleRegressor``, so the regressor is written as a classifier
with three classes whose leaves all add to class 0, ``post_transform = NONE``
and ``base_values[0]`` = the XGBoost base score. Column 0 of ``scores`` is then
exactly the regressor margin; columns 1 and 2 stay at zero. (Two classes would
trigger tract's binary-classifier output layout, three do not.) Node ids are
renumbered densely per tree: pruned XGBoost trees leave gaps, and tract
addresses nodes by id.

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


def tree_arrays(b: xgb.Booster):
    """Node arrays from the per-tree JSON dump, one tree at a time (the
    HCD2021 y model has ~19 M nodes; whole-model Python objects do not fit)."""
    keys = ["tree", "node", "feat", "thr", "yes", "no", "miss_left", "leaf_val", "is_leaf"]
    dtypes = [np.int32, np.int32, np.int16, np.float32, np.int32, np.int32, bool, np.float32, bool]
    parts = {k: [] for k in keys}
    for t, dump in enumerate(b.get_dump(dump_format="json")):
        cols = {k: [] for k in keys}
        stack = [json.loads(dump)]
        while stack:
            n = stack.pop()
            leaf = "leaf" in n
            cols["tree"].append(t)
            cols["node"].append(n["nodeid"])
            cols["is_leaf"].append(leaf)
            cols["leaf_val"].append(n["leaf"] if leaf else 0.0)
            cols["feat"].append(0 if leaf else int(n["split"].lstrip("f")))
            cols["thr"].append(0.0 if leaf else n["split_condition"])
            cols["yes"].append(0 if leaf else n["yes"])
            cols["no"].append(0 if leaf else n["no"])
            cols["miss_left"].append(False if leaf else n["missing"] == n["yes"])
            if not leaf:
                stack.extend(n["children"])
        # Pruned XGBoost trees leave gaps in node ids; tract indexes nodes by
        # id within a tree, so renumber them densely in id order.
        ids = np.asarray(cols["node"])
        dense = np.sort(ids)
        for k in ("node", "yes", "no"):
            v = np.asarray(cols[k])
            cols[k] = np.where(np.asarray(cols["is_leaf"]) & (k != "node"), 0, np.searchsorted(dense, v))
        for k, dt in zip(keys, dtypes):
            parts[k].append(np.asarray(cols[k], dtype=dt))
    return {k: np.concatenate(v) for k, v in parts.items()}


def to_onnx(b: xgb.Booster, name: str) -> onnx.ModelProto:
    a = tree_arrays(b)
    order = np.lexsort((a["node"], a["tree"]))
    a = {k: v[order] for k, v in a.items()}
    leaf = a["is_leaf"]
    modes = np.where(leaf, "LEAF", "BRANCH_LT")
    node = helper.make_node(
        "TreeEnsembleClassifier",
        ["features"],
        ["label", "scores"],
        domain="ai.onnx.ml",
        classlabels_int64s=[0, 1, 2],
        base_values=[booster_base_score(b), 0.0, 0.0],
        post_transform="NONE",
        nodes_treeids=a["tree"].tolist(),
        nodes_nodeids=a["node"].tolist(),
        nodes_featureids=a["feat"].tolist(),
        nodes_values=a["thr"].tolist(),
        nodes_modes=modes.tolist(),
        nodes_truenodeids=a["yes"].tolist(),
        nodes_falsenodeids=a["no"].tolist(),
        nodes_missing_value_tracks_true=a["miss_left"].astype(np.int64).tolist(),
        class_treeids=a["tree"][leaf].tolist(),
        class_nodeids=a["node"][leaf].tolist(),
        class_ids=[0] * int(leaf.sum()),
        class_weights=a["leaf_val"][leaf].tolist(),
    )
    graph = helper.make_graph(
        [node],
        name,
        [helper.make_tensor_value_info("features", TensorProto.FLOAT, ["N", N_FEATURES])],
        [
            helper.make_tensor_value_info("label", TensorProto.INT64, ["N"]),
            helper.make_tensor_value_info("scores", TensorProto.FLOAT, ["N", 3]),
        ],
    )
    model = helper.make_model(
        graph,
        opset_imports=[helper.make_opsetid("", 17), helper.make_opsetid("ai.onnx.ml", 3)],
        producer_name="onnx-infer export_ms2pip_to_onnx.py",
    )
    model.ir_version = 8
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
            got = sess.run(["scores"], {"features": x})[0][:, 0]
            err = float(np.abs(ref - got).max())
            manifest = {
                "id": stem,
                "description": f"MS²PIP {model} {ion}-ion XGBoost model as ONNX TreeEnsembleClassifier "
                "(scores[:, 0] = log2(intensity + 0.001), one row per cleavage site)",
                "onnx": onnx_path.name,
                "task": "ms2_intensity",
                "inputs": [{"name": "features", "shape": [None, N_FEATURES]}],
                "outputs": [{"name": "label"}, {"name": "scores", "shape": [None, 3]}],
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
