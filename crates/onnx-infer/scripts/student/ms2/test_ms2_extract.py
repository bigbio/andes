"""Regression tests for the andes-ms2 extraction scripts (pytest; numpy and pyarrow only)."""
from __future__ import annotations

import glob
import json
import os
import subprocess
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ms2_common import fragment_mz, precursor_mz, valid_mask  # noqa: E402
from mod_composition import MASSES  # noqa: E402

AMIDATED = 2


def run_script(name, *args):
    r = subprocess.run([sys.executable, os.path.join(HERE, name), *map(str, args)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return r.stdout


def perfect_spectrum(seq, mods, z):
    theo = fragment_mz(seq, mods)[valid_mask(seq, mods, z)]
    mz = np.sort(theo[np.isfinite(theo)])
    return mz.tolist(), np.full(len(mz), 1000.0).tolist()


def test_flats_keep_c_terminal_amidation(tmp_path):
    seq, z = "PEPTIDEK", 2
    rows = {k: [] for k in ("seq", "charge", "prec_mz", "res_mod_pos", "res_mod_delta",
                            "nterm_delta", "cterm_delta", "mz", "intensity")}
    for mods, cterm in (([], 0.0), ([(len(seq) + 1, AMIDATED)], MASSES[AMIDATED])):
        mz, inten = perfect_spectrum(seq, mods, z)
        for k, v in (("seq", seq), ("charge", z), ("prec_mz", precursor_mz(seq, mods, z)), ("res_mod_pos", []),
                     ("res_mod_delta", []), ("nterm_delta", 0.0), ("cterm_delta", cterm), ("mz", mz),
                     ("intensity", inten)):
            rows[k].append(v)
    schema = pa.schema([("seq", pa.string()), ("charge", pa.int32()), ("prec_mz", pa.float64()),
                        ("res_mod_pos", pa.list_(pa.int32())), ("res_mod_delta", pa.list_(pa.float64())),
                        ("nterm_delta", pa.float64()), ("cterm_delta", pa.float64()),
                        ("mz", pa.list_(pa.float64())), ("intensity", pa.list_(pa.float64()))])
    flat = tmp_path / "flat.parquet"
    pq.write_table(pa.table(rows, schema=schema), flat)
    out = tmp_path / "out"
    run_script("ms2_extract_flats.py", "--dataset", "PXD000000", "--out", out, flat)
    st = json.load(open(out / "flats.json"))
    assert st["unmapped"] == 0 and st["kept"] == 2, st
    kept = pq.read_table(glob.glob(str(out / "part-flat-*.parquet"))[0]).column("peptidoform").to_pylist()
    assert kept == ["PEPTIDEK", "PEPTIDEK-[UNIMOD:2]"], kept

