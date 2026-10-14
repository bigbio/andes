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
from ms2_selftest import write_mzml  # noqa: E402
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


def test_flats_decode_summed_nterm_plus_residue_delta():
    import re
    from pathlib import Path
    from ms2_extract_flats import parse_flat
    fixture = Path(HERE).parents[2] / "fixtures" / "andes_writer_peptides.tsv"
    andes = fixture.read_text().splitlines()[1].split("\t")[0]  # _.C+99.03203PEPTIDEK.R
    delta = float(re.search(r"C\+([0-9.]+)", andes).group(1))
    assert parse_flat("CPEPTIDEK", [1], [delta], 0.0, 0.0) == ("CPEPTIDEK", [(0, 1), (1, 4)])
    # C-terminal amidation still decodes alongside
    assert parse_flat("CPEPTIDEK", [1], [delta], 0.0, MASSES[AMIDATED]) == \
        ("CPEPTIDEK", [(0, 1), (1, 4), (10, AMIDATED)])
    assert parse_flat("CPEPTIDEK", [1], [1.23456], 0.0, 0.0) is None


def test_keep_mzml_preserves_the_converted_file(tmp_path):
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    (raw_dir / "run_a.raw").write_bytes(b"")
    source = tmp_path / "source.mzML"
    write_mzml(source, [])
    # Stand-in converter: copies the prepared mzML to <outdir>/<stem>.mzML.
    trfp = tmp_path / "trfp.sh"
    trfp.write_text('#!/bin/sh\n'
                    'while [ $# -gt 0 ]; do case "$1" in -i) I="$2"; shift 2;; -o) O="$2"; shift 2;; *) shift;; esac; done\n'
                    f'cp "{source}" "$O/$(basename "$I" .raw).mzML"\n')
    trfp.chmod(0o755)
    selected = tmp_path / "selected.parquet"
    pq.write_table(pa.table({"run": ["run_a"], "scan": [1], "peptidoform": ["PEPTIDEK"], "charge": [2],
                             "pep": [0.001]}), selected)
    tasks = tmp_path / "tasks.tsv"
    tasks.write_text(f"task\tdataset\truns\tselected\traw_dir\n1\tPXD000000\trun_a\t{selected}\t{raw_dir}\n")
    work = tmp_path / "work"
    common = ["--tasks", tasks, "--task-index", 1, "--out-dir", tmp_path / "shards", "--work-dir", work,
              "--trfp", trfp, "--overwrite"]

    run_script("ms2_extract_raw.py", *common, "--keep-mzml")
    kept = glob.glob(str(work / "**" / "run_a.mzML"), recursive=True)
    assert len(kept) == 1, kept
    assert open(kept[0]).read() == open(source).read()

    # Without the flag the scratch directory is removed.
    for d in glob.glob(str(work / "*")):
        subprocess.run(["rm", "-rf", d], check=True)
    run_script("ms2_extract_raw.py", *common)
    assert glob.glob(str(work / "**" / "*.mzML"), recursive=True) == []


def test_failed_run_does_not_complete_the_task(tmp_path):
    seq, z = "PEPTIDEK", 2
    mz, inten = perfect_spectrum(seq, [], z)
    pmz = precursor_mz(seq, [], z)
    mzml_dir = tmp_path / "mzml"
    mzml_dir.mkdir()
    write_mzml(mzml_dir / "run_a.mzML",
               [(1, f"FTMS + p NSI d Full ms2 {pmz:.2f}@hcd27.00 [100.00-2000.00]", pmz, z, mz, inten)])
    selected = tmp_path / "selected.parquet"
    pq.write_table(pa.table({"run": ["run_a", "run_b"], "scan": [1, 1], "peptidoform": [seq, seq],
                             "charge": [z, z], "pep": [0.001, 0.001]}), selected)
    tasks = tmp_path / "tasks.tsv"
    tasks.write_text(f"task\tdataset\truns\tselected\traw_dir\n1\tPXD000000\trun_a|run_b\t{selected}\t{tmp_path}\n")
    shards = tmp_path / "shards"
    cmd = [sys.executable, os.path.join(HERE, "ms2_extract_raw.py"), "--tasks", str(tasks), "--task-index", "1",
           "--out-dir", str(shards), "--work-dir", str(tmp_path / "work"), "--mzml-dir", str(mzml_dir)]
    part = shards / "PXD000000" / "part-00001.parquet"
    summary = shards / "PXD000000" / "part-00001.json"

    # run_b has no mzML: the task fails, writes no part and is not complete.
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 1, r.stdout + r.stderr
    assert not part.exists()
    st = json.load(open(summary))
    assert st["complete"] is False and st["total"]["runs_failed"] == 1, st

    # The resubmission is not skipped; once run_b is readable the task completes.
    write_mzml(mzml_dir / "run_b.mzML",
               [(1, f"FTMS + p NSI d Full ms2 {pmz:.2f}@hcd27.00 [100.00-2000.00]", pmz, z, mz, inten)])
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert json.load(open(summary))["complete"] is True
    assert pq.read_table(part).column("run").to_pylist() == ["run_a", "run_b"]

    # A complete task is skipped.
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0 and "skipping" in r.stdout, r.stdout + r.stderr
