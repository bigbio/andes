"""Tests for the per-dataset training cap, the training cell and the cross-run hard rule (pytest)."""
from __future__ import annotations

import os
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import row_cell  # noqa: E402
from ms2_mix import plan  # noqa: E402
from ms2_select import explained_cuts, is_hard  # noqa: E402


def shard(path, ds, run, n, activation="HCD-FTMS", instrument="Q Exactive HF", label="none"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pq.write_table(pa.table({
        "dataset": [ds] * n, "run": [run] * n, "label": [label] * n, "activation": [activation] * n,
        "instrument": [instrument] * n, "n_matched": [6] * n, "explained": np.linspace(0, 1, n),
        "peptidoform": [f"PEPTIDE{i % 7}K" for i in range(n)]}), path)


def test_row_cell():
    assert row_cell("none", "HCD-FTMS", "Q Exactive HF") == "QE-Exploris"
    assert row_cell("none", "HCD-FTMS", "Orbitrap Fusion Lumos") == "Tribrid"
    assert row_cell("none", "CID-ITMS", "Orbitrap Fusion Lumos") == "CID-ITMS"
    assert row_cell("TMT", "HCD-FTMS", "Q Exactive") == "TMT"
    assert row_cell("TMT", "CID-ITMS", "Orbitrap Fusion Lumos") == "TMT-CID-ITMS"
    assert row_cell("iTRAQ", "CID-ITMS", "Orbitrap Fusion Lumos") == "iTRAQ"
    assert row_cell("none", "HCD-TOF", "timsTOF Pro") == "timsTOF"
    assert row_cell("none", "HCD-FTMS", "Orbitrap Fusion Lumos", hla=True) == "HLA"


def test_cap_per_dataset_and_cell(tmp_path):
    big = str(tmp_path / "raw" / "PXDA" / "part-00001.parquet")
    hla = str(tmp_path / "raw" / "PXDA" / "part-00002.parquet")
    small = str(tmp_path / "raw" / "PXDB" / "part-00001.parquet")
    uvpd = str(tmp_path / "raw" / "PXDB" / "part-00002.parquet")
    shard(big, "PXDA", "r1", 500)
    shard(hla, "PXDA", "r2.raw", 300)
    shard(small, "PXDB", "s1", 50)
    shard(uvpd, "PXDB", "s2", 40, activation="UVPD")
    files = [big, hla, small, uvpd]
    keep, mix = plan(files, {("PXDA", "r2")}, 4, 0.0, 100, 200, seed=1)
    got = {(d, c): (n, e, k) for d, c, n, e, k, _, _ in mix}
    distinct = {(d, c): (de, dk) for d, c, _, _, _, de, dk in mix}
    assert distinct[("PXDA", "QE-Exploris")] == (7, 7)
    assert got[("PXDA", "QE-Exploris")] == (500, 500, 100)
    assert got[("PXDA", "HLA")] == (300, 300, 200)
    assert got[("PXDB", "QE-Exploris")] == (90, 50, 50)
    assert keep[big].sum() == 100 and keep[hla].sum() == 200 and keep[small].all() and not keep[uvpd].any()
    again, _ = plan(files, {("PXDA", "r2")}, 4, 0.0, 100, 200, seed=1)
    assert all((again[f] == keep[f]).all() for f in files)
    nocap, _ = plan(files, set(), 4, 0.0, 0, 0, seed=1)
    assert nocap[big].all() and nocap[hla].all()


def test_cross_run_cut_differs_from_per_run():
    acts = ["HCD-FTMS"] * 8
    low_run, high_run = [0.1, 0.2, 0.3, 0.4], [0.6, 0.7, 0.8, 0.9]
    per_run = explained_cuts(acts[:4], low_run)["HCD-FTMS"]
    pooled = explained_cuts([("PXDA", "QE-Exploris", a) for a in acts], low_run + high_run)
    cut = pooled[("PXDA", "QE-Exploris", "HCD-FTMS")]
    assert per_run < cut
    assert sum(is_hard("PEPTIDEKAAR", e, [], cut) for e in low_run + high_run) == 2
    assert sum(is_hard("PEPTIDEKAAR", e, [], per_run) for e in low_run) == 1
    assert is_hard("PEPTIDEK", 0.9, [], cut)
    assert is_hard("PEPTIDEKAAR", 0.9, [(3, 213)], cut)
