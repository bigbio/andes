#!/usr/bin/env python3
"""Parity tests for committed onnx-infer fixtures (no torch required).

Checks:
  1. AGBD decode → tree walk matches ONNX Runtime TreeEnsemble, and both match
     the expected scores in fixtures/parity/agbd.json, which the native Rust
     decoder is checked against (tests/agbd_native.rs)
  2. intensity_mlp MatMul+ReLU matches numpy reference
  3. export_agbd_to_onnx round-trip still loads

Exit 0 on success. Requires: numpy, onnx, onnxruntime
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import onnxruntime as ort
from onnx import save

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "fixtures"
sys.path.insert(0, str(ROOT / "scripts"))
from export_agbd_to_onnx import parse_agbd, parse_agbd_full, trees_to_onnx  # noqa: E402


def eval_agbd_trees(trees, x: np.ndarray) -> float:
    """Mirror scoring::gbdt_eval::Tree::eval (NaN → default_left)."""
    total = 0.0
    for t in trees:
        node = 0
        while True:
            feat = t["feature"][node]
            if feat < 0:
                total += t["value"][node]
                break
            v = x[feat] if feat < len(x) else np.nan
            if np.isnan(v):
                go_left = t["default_left"][node] == 1
            else:
                go_left = v <= t["threshold"][node]
            node = t["left"][node] if go_left else t["right"][node]
    return total


def check_tree_fixture(agbd_name: str, onnx_name: str, rows: np.ndarray, expected: np.ndarray):
    nfeat, flags, trees = parse_agbd((FIX / agbd_name).read_bytes())
    assert rows.shape[1] == nfeat, (rows.shape, nfeat)
    sess = ort.InferenceSession(str(FIX / onnx_name), providers=["CPUExecutionProvider"])
    y_ort = sess.run(None, {"features": rows.astype(np.float32)})[0].reshape(-1)
    y_nat = np.array([eval_agbd_trees(trees, r) for r in rows], dtype=np.float32)
    if not np.allclose(y_ort, y_nat, atol=1e-5):
        raise AssertionError(f"{onnx_name}: ort={y_ort} native={y_nat} flags={flags}")
    if not np.allclose(y_ort, expected, atol=1e-5):
        raise AssertionError(f"{onnx_name}: ort={y_ort} expected (native Rust)={expected}")
    # round-trip export
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "rt.onnx"
        save(trees_to_onnx(nfeat, trees), out)
        y2 = ort.InferenceSession(str(out), providers=["CPUExecutionProvider"]).run(
            None, {"features": rows.astype(np.float32)}
        )[0].reshape(-1)
        if not np.allclose(y_ort, y2, atol=1e-5):
            raise AssertionError(f"{onnx_name} round-trip drift: {y_ort} vs {y2}")
    print(f"OK {onnx_name}: {len(rows)} rows, max|Δ|={np.max(np.abs(y_ort - y_nat)):.2e}")


def test_intensity_mlp():
    W = np.array([[0.5], [-0.25], [0.75], [0.1]], dtype=np.float32)
    b = np.array([0.05], dtype=np.float32)
    x = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [-1, -1, -1, -1]], dtype=np.float32)
    expect = np.maximum(0, x @ W + b).reshape(-1)
    sess = ort.InferenceSession(
        str(FIX / "intensity_mlp.onnx"), providers=["CPUExecutionProvider"]
    )
    got = sess.run(None, {"features": x})[0].reshape(-1)
    if not np.allclose(got, expect, atol=1e-5):
        raise AssertionError(f"mlp: got={got} expect={expect}")
    print(f"OK intensity_mlp: {got}")


def load_parity():
    entries = json.loads((FIX / "parity" / "agbd.json").read_text())
    out = []
    for e in entries:
        rows = np.array(
            [[np.nan if v is None else v for v in r] for r in e["rows"]], dtype=np.float32
        )
        out.append((e["agbd"], e["onnx"], rows, np.array(e["expected"], dtype=np.float32)))
    return out


def test_tree_fixtures():
    for agbd_name, onnx_name, rows, expected in load_parity():
        check_tree_fixture(agbd_name, onnx_name, rows, expected)


def test_parser_accepts_native_blobs():
    """Every AGBD fixture is byte-identical to `GbdtPeakModel::to_bytes`
    (tests/agbd_native.rs), so these are natively serialized models."""
    for name in ("tree_leq.agbd", "peak_gbdt.agbd"):
        _, _, trees, iso_x, iso_y = parse_agbd_full((FIX / name).read_bytes())
        assert trees and iso_x == [] and iso_y == [], name
    nfeat, _, trees, iso_x, iso_y = parse_agbd_full((FIX / "peak_gbdt_calibrated.agbd").read_bytes())
    assert (nfeat, len(trees)) == (3, 3)
    assert np.allclose(iso_x, [0.0, 0.5, 1.0]) and np.allclose(iso_y, [0.1, 0.4, 0.9])
    assert parse_agbd((FIX / "peak_gbdt_calibrated.agbd").read_bytes())[2] == trees
    print("OK parser accepts native uncalibrated and calibrated AGBD blobs")


def _expect_error(blob: bytes, needle: str):
    try:
        parse_agbd(blob)
    except ValueError as e:
        assert needle in str(e), (needle, str(e))
        return
    raise AssertionError(f"parse_agbd accepted a malformed blob (expected {needle!r})")


def test_parser_rejects_malformed_footer():
    blob = (FIX / "peak_gbdt_calibrated.agbd").read_bytes()
    _expect_error(blob + b"\0", "1 trailing bytes")
    _expect_error(blob[:-4], "truncated AGBD blob")
    _expect_error((FIX / "tree_leq.agbd").read_bytes()[:-4], "calibration count")
    print("OK parser rejects trailing and truncated footers")


def main():
    test_parser_accepts_native_blobs()
    test_parser_rejects_malformed_footer()
    test_tree_fixtures()
    test_intensity_mlp()
    print("all fixture tests passed")


if __name__ == "__main__":
    main()
