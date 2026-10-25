"""Tests for ms2_extract_msnet.py (pytest; numpy and pyarrow only)."""
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
from ms2_extract_msnet import conditions, parse_named  # noqa: E402


def test_named_modifications_parse():
    assert parse_named("PEPM[Oxidation]K") == ("PEPMK", [(4, 35)])
    assert parse_named("[Acetyl]-PEPC[Carbamidomethyl]K") == ("PEPCK", [(0, 1), (4, 4)])
    assert parse_named("PEPK-[Amidated]") == ("PEPK", [(4 + 1, 2)])
    assert parse_named("NN[Deamidated]GAGYK") == ("NNGAGYK", [(2, 7)])
    assert parse_named("PEP[Mystery]K") is None
    assert parse_named("PEP[+15.9949]K") is None


def test_conditions_for_both_cv_layouts():
    a = [{"cv_name": "instrument", "cv_value": "Q Exactive"}, {"cv_name": "dissociation method", "cv_value": "HCD"},
         {"cv_name": "normalized collision energy", "cv_value": "27.00"}, {"cv_name": "mass analyzer", "cv_value": "FTMS"}]
    assert conditions(a) == ("HCD-FTMS", "FTMS", 27.0, "Q Exactive")
    b = {"Instrument": "LTQ Orbitrap Velos", "Fragmentation Method": "CID", "Collision Energy": "35.00",
         "Mass Analyzer": "ITMS"}
    assert conditions(b) == ("CID-ITMS", "ITMS", 35.0, "LTQ Orbitrap Velos")
    act, an, nce, inst = conditions({"Instrument": "timsTOF Pro", "Fragmentation Method": "CID", "Collision Energy": "30",
                                     "Mass Analyzer": "TOF"})
    assert act == "unknown" and np.isnan(nce)


def spectrum(seq, mods, z):
    theo = fragment_mz(seq, mods)[valid_mask(seq, mods, z)]
    mz = np.sort(theo)
    return mz.tolist(), np.full(len(mz), 1000.0).tolist()


def test_both_layouts_extract(tmp_path):
    seq, z = "PEPTIDEKR", 2
    mz, inten = spectrum(seq, [], z)
    pmz = precursor_mz(seq, [], z)
    cv_a = [{"cv_name": "instrument", "cv_value": "Q Exactive HF"}, {"cv_name": "dissociation method", "cv_value": "HCD"},
            {"cv_name": "normalized collision energy", "cv_value": "28.00"}, {"cv_name": "mass analyzer", "cv_value": "FTMS"}]
    cv_b = {"Instrument": "LTQ Orbitrap Velos", "Fragmentation Method": "CID", "Collision Energy": "35.00",
            "Mass Analyzer": "ITMS"}
    base = {"peptidoform": [seq, seq, seq], "scan": [[1], [2], [3]], "mz_array": [mz] * 3, "intensity_array": [inten] * 3,
            "posterior_error_probability": [1e-5, 1e-5, 0.5], "is_decoy": [False, False, False]}
    a = pa.table({**base, "charge": [z] * 3, "observed_mz": [pmz] * 3, "run_file_name": ["r"] * 3, "cv_params": [cv_a] * 3})
    b = pa.table({**base, "precursor_charge": [z] * 3, "exp_mass_to_charge": [pmz] * 3, "reference_file_name": ["r"] * 3,
                  "cv_params": [cv_b] * 3})
    for name, t in (("a", a), ("b", b)):
        f = tmp_path / f"{name}.parquet"
        pq.write_table(t, f)
        out = tmp_path / f"out_{name}"
        r = subprocess.run([sys.executable, os.path.join(HERE, "ms2_extract_msnet.py"), "--dataset", "PXD000000",
                            "--out", str(out), str(f)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        st = json.load(open(out / "msnet.json"))
        # the two PEP <= 0.01 rows share one key: one kept, one key_full is only hit beyond --per-key
        assert st["kept"] == 2 and st["decoy_or_pep"] == 1, st
        rows = pq.read_table(glob.glob(str(out / "part-*.parquet"))[0]).to_pylist()
        assert rows[0]["activation"] == ("HCD-FTMS" if name == "a" else "CID-ITMS")
        assert rows[0]["nce"] == (28.0 if name == "a" else 35.0)
