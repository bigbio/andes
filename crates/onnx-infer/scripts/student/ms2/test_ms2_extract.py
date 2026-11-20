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
    seq, z = "PEPTIDEKAR", 2
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
    # both rows are annotated; amidation is outside the common mods, so that row is a hard example and,
    # with no submitter identification in a flat, it is dropped
    assert st["unmapped"] == 0 and st["kept"] == 2 and st["hard_dropped"] == 1, st
    kept = pq.read_table(glob.glob(str(out / "part-flat-*.parquet"))[0]).column("peptidoform").to_pylist()
    assert kept == ["PEPTIDEKAR"], kept


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
    converter = tmp_path / "converter.sh"
    converter.write_text('#!/bin/sh\n'
                    'while [ $# -gt 0 ]; do case "$1" in -i) I="$2"; shift 2;; -o) O="$2"; shift 2;; *) shift;; esac; done\n'
                    f'cp "{source}" "$O/$(basename "$I" .raw).mzML"\n')
    converter.chmod(0o755)
    selected = tmp_path / "selected.parquet"
    pq.write_table(pa.table({"run": ["run_a"], "scan": [1], "peptidoform": ["PEPTIDEK"], "charge": [2],
                             "pep": [0.001]}), selected)
    tasks = tmp_path / "tasks.tsv"
    tasks.write_text(f"task\tdataset\truns\tselected\traw_dir\n1\tPXD000000\trun_a\t{selected}\t{raw_dir}\n")
    work = tmp_path / "work"
    common = ["--tasks", tasks, "--task-index", 1, "--out-dir", tmp_path / "shards", "--work-dir", work,
              "--converter", converter, "--overwrite"]

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
    seq, z = "PEPTIDEKAR", 2
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


def test_v2_masses_equal_model_mass_rs():
    """Every fragment constant equals the f64 that ``crates/model/src/mass.rs`` computes, bit for bit."""
    import re
    from pathlib import Path
    import ms2_common as mc
    src = (Path(HERE).parents[4] / "crates" / "model" / "src" / "mass.rs").read_text()
    env = {}
    for name, expr in re.findall(r"pub const (\w+): f64 = ([^;]+);", src):
        env[name] = eval(expr.replace("_f64", ""), {}, dict(env))
    for py, rs in (("PROTON", "PROTON"), ("H2O", "H2O"), ("NH3", "NH3"), ("CO", "CO"), ("H3PO4", "H3PO4"),
                   ("CH4SO", "CH4SO"), ("Z_DOT_OFFSET", "Z_DOT_OFFSET")):
        assert getattr(mc, py) == env[rs], (py, getattr(mc, py), env[rs])
    g = 2.0 * env["C"] + 3.0 * env["H"] + 1.0 * env["N"] + 1.0 * env["O"] + 0.0 * env["S"]
    assert mc.RESIDUE_MASS["G"] == g
    m = 5.0 * env["C"] + 9.0 * env["H"] + 1.0 * env["N"] + 1.0 * env["O"] + 1.0 * env["S"]
    assert mc.RESIDUE_MASS["M"] == m


def ethcd_spectrum(seq, mods, z):
    """Every valid v2 channel of an EThcD scan at intensity 1000, c1 of the first site at 4000."""
    from ms2_common import ION_TYPES_V2, activation_index
    act = activation_index("EThcD")
    theo = fragment_mz(seq, mods, ION_TYPES_V2)
    valid = valid_mask(seq, mods, z, ION_TYPES_V2, act)
    mz = np.sort(theo[valid])
    inten = np.full(len(mz), 1000.0)
    inten[np.argmin(np.abs(mz - theo[0, 15]))] = 4000.0
    return mz.tolist(), inten.tolist()


def test_v2_channels_annotated(tmp_path):
    """An EThcD scan read from mzML annotates c and z-dot channels and writes a 2.0 shard."""
    from ms2_common import N_ION_V2
    seq, z = "PEPTIDEKAR", 3
    mz, inten = ethcd_spectrum(seq, [], z)
    pmz = precursor_mz(seq, [], z)
    mzml_dir = tmp_path / "mzml"
    mzml_dir.mkdir()
    write_mzml(mzml_dir / "run_a.mzML",
               [(1, f"FTMS + p NSI d Full ms2 {pmz:.2f}@etd50.00@hcd20.00 [100.00-2000.00]", pmz, z, mz, inten)])
    selected = tmp_path / "selected.parquet"
    pq.write_table(pa.table({"run": ["run_a"], "scan": [1], "peptidoform": [seq], "charge": [z], "pep": [0.001]}),
                   selected)
    tasks = tmp_path / "tasks.tsv"
    tasks.write_text(f"task\tdataset\truns\tselected\traw_dir\n1\tPXD000000\trun_a\t{selected}\t{tmp_path}\n")
    shards = tmp_path / "shards"
    run_script("ms2_extract_raw.py", "--tasks", tasks, "--task-index", 1, "--out-dir", shards,
               "--work-dir", tmp_path / "work", "--mzml-dir", mzml_dir)
    t = pq.read_table(shards / "PXD000000" / "part-00001.parquet")
    row = t.to_pylist()[0]
    assert row["activation"] == "EThcD" and row["nce"] == 20.0
    v = np.asarray(row["intensities"]).reshape(len(seq) - 1, N_ION_V2)
    assert v[0, 15] == 1.0
    assert (v[:, 15:19][v[:, 15:19] >= 0] > 0).all() and (v[:, 15:19] > 0).sum() > 10
    assert (v[1, 15:17] == -1).all(), "no c ion N-terminal to proline"
    assert (v[:, 12:14] > 0).any(), "3+ fragments of a 3+ precursor"


def test_shard_schema_version(tmp_path):
    from ms2_common import SHARD_COLUMNS, shard_table
    row = {"dataset": "D", "run": "r", "scan": 1, "peptidoform": "PEPTIDEK", "sequence": "PEPTIDEK", "charge": 2,
           "nce": 27.0, "activation": "HCD-FTMS", "instrument": "Q Exactive", "pep": 0.001, "n_matched": 6,
           "explained": 0.5, "intensities": np.zeros((7, 19), np.float32), "label": "TMT", "hard": True,
           "agree": None}
    t = shard_table([row])
    assert t.schema.metadata[b"andes.ms2_shard.schema"] == b"2.0"
    assert t.column_names == SHARD_COLUMNS
    assert {"label", "hard", "agree"} <= set(t.column_names)
    pq.write_table(t, tmp_path / "s.parquet")
    back = pq.read_table(tmp_path / "s.parquet")
    assert back.schema.metadata[b"andes.ms2_shard.schema"] == b"2.0"
    assert back.column("agree").to_pylist() == [None] and back.column("hard").to_pylist() == [True]
    assert back.column("label").to_pylist() == ["TMT"]


def _row(seq, explained, pep_form=None, run="r", scan=1, act="HCD-FTMS"):
    from ms2_common import parse_proforma
    form = pep_form or seq
    s, mods = parse_proforma(form)
    return {"dataset": "D", "run": run, "scan": scan, "peptidoform": form, "sequence": s, "mods": mods,
            "charge": 2, "nce": 27.0, "activation": act, "instrument": "Q Exactive", "pep": 0.001, "n_matched": 6,
            "explained": explained, "intensities": np.zeros((len(s) - 1, 19), np.float32)}


def test_hard_requires_agree():
    """Hard rows (length <= 9, bottom-quartile matched share, uncommon mod) need submitter agreement."""
    from ms2_select import finalize_rows
    rows = [
        _row("SHORTPEK", 0.6, scan=1),                       # hard (length 8), no submitter ID
        _row("SHORTPEK", 0.6, scan=2),                       # hard, submitter agrees
        _row("SHORTPEK", 0.6, scan=3),                       # hard, submitter disagrees
        _row("LONGPEPTIDEK", 0.6, scan=4),                   # easy, no submitter ID: kept, agree null
        _row("LONGPEPTIDEK", 0.6, "LONGPEPTIDEK-[UNIMOD:2]", scan=5),  # hard (amidation), no ID
        _row("LONGPEPTIDEK", 0.6, "LONGPEPTIDEM[UNIMOD:35]K", scan=6),  # easy: oxidation
        _row("LONGPEPTIDEK", 0.6, "[UNIMOD:737]-LONGPEPTIDEK[UNIMOD:737]", scan=7),  # easy: TMT
    ] + [_row("ANOTHERLONGPEPK", 0.7 + 0.01 * k, scan=10 + k) for k in range(5)] + [
        _row("ANOTHERLONGPEPK", 0.01, scan=20),              # hard: bottom quartile of explained
    ]
    sub = {("r", 2): "SHORTPEK", ("r", 3): "SHORTPEKK", ("r", 4): "LONGPEPTIDEK"}
    kept, st = finalize_rows(rows, sub)
    by_scan = {r["scan"]: r for r in kept}
    assert set(by_scan) == {2, 4, 6, 7, 10, 11, 12, 13, 14}, sorted(by_scan)
    assert by_scan[2]["hard"] and by_scan[2]["agree"] is True
    assert not by_scan[4]["hard"] and by_scan[4]["agree"] is True
    assert by_scan[6]["agree"] is None and not by_scan[6]["hard"] and not by_scan[7]["hard"]
    assert st["hard"] == 5 and st["hard_dropped"] == 4 and st["agree_true"] == 2 and st["agree_false"] == 1
    # without any submitter table every hard row is dropped and agree stays null
    kept, st = finalize_rows([dict(r) for r in rows], None)
    assert all(r["agree"] is None and not r["hard"] for r in kept)


def test_submitter_agreement_ignores_i_l_and_carbamidomethyl():
    from ms2_select import same_identification
    assert same_identification("PEPIC[UNIMOD:4]K", "PEPLCK")
    assert not same_identification("PEPICK", "PEPICR")
    assert not same_identification("PEPM[UNIMOD:35]K", "PEPMK")


def test_run_label_from_tags():
    from ms2_select import finalize_rows
    tmt = "[UNIMOD:737]-LONGPEPTIDEK[UNIMOD:737]"
    rows = [_row("LONGPEPTIDEK", 0.6, tmt, run="t", scan=1), _row("LONGPEPTIDEK", 0.6, tmt, run="t", scan=2),
            _row("LONGPEPTIDEK", 0.6, run="t", scan=3),
            _row("LONGPEPTIDEK", 0.6, "[UNIMOD:214]-LONGPEPTIDEK[UNIMOD:214]", run="i", scan=4),
            _row("LONGPEPTIDEK", 0.6, run="n", scan=5)]
    kept, _ = finalize_rows(rows, None)
    assert {r["scan"]: r["label"] for r in kept} == {1: "TMT", 2: "TMT", 3: "TMT", 4: "iTRAQ", 5: "none"}


def test_filter_v2_electron_activations():
    from ms2_common import parse_filter_v2
    assert parse_filter_v2("FTMS + p NSI d Full ms2 445.12@etd50.00@hcd20.00 [110.00-905.00]") == ("EThcD", 20.0)
    act, nce = parse_filter_v2("ITMS + c NSI d Full ms2 445.12@etd70.00 [110.00-905.00]")
    assert act == "ETD" and np.isnan(nce)
    assert parse_filter_v2("FTMS + p NSI d Full ms2 500.00@hcd28.00 [100.00-1000.00]") == ("HCD-FTMS", 28.0)
