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
    # amidation is outside the common mods: a hard example, kept in the candidate shard
    assert st["unmapped"] == 0 and st["kept"] == 2 and st["hard"] == 1, st
    kept = pq.read_table(glob.glob(str(out / "part-flat-*.parquet"))[0]).column("peptidoform").to_pylist()
    assert kept == ["PEPTIDEKAR", "PEPTIDEKAR-[UNIMOD:2]"], kept


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
           "source_engine": "comet", "agree_andes": None, "agree_comet": True, "agree_msfragger": False,
           "agree_sage": None, "agree_maxquant": None, "agree_submitter": None, "msnet_agree": None,
           "n_engines": 1}
    t = shard_table([row])
    assert t.schema.metadata[b"andes.ms2_shard.schema"] == b"2.0"
    assert t.schema.metadata[b"andes.ms2_shard.finalized"] == b"false"
    assert shard_table([row], finalized=True).schema.metadata[b"andes.ms2_shard.finalized"] == b"true"
    assert t.column_names == SHARD_COLUMNS
    assert {"label", "hard", "agree_andes", "agree_comet", "agree_msfragger", "agree_sage", "agree_maxquant",
            "agree_submitter", "n_engines"} <= set(t.column_names)
    pq.write_table(t, tmp_path / "s.parquet")
    back = pq.read_table(tmp_path / "s.parquet")
    assert back.schema.metadata[b"andes.ms2_shard.schema"] == b"2.0"
    assert back.column("agree_andes").to_pylist() == [None] and back.column("agree_comet").to_pylist() == [True]
    assert back.column("hard").to_pylist() == [True] and back.column("label").to_pylist() == ["TMT"]


def _row(seq, explained, pep_form=None, run="r", scan=1, act="HCD-FTMS"):
    from ms2_common import parse_proforma
    form = pep_form or seq
    s, mods = parse_proforma(form)
    return {"dataset": "D", "run": run, "scan": scan, "peptidoform": form, "sequence": s, "mods": mods,
            "charge": 2, "nce": 27.0, "activation": act, "instrument": "Q Exactive", "pep": 0.001, "n_matched": 6,
            "explained": explained, "intensities": np.zeros((len(s) - 1, 19), np.float32)}


def test_hard_flags_at_the_candidate_stage():
    """Hard = length <= 9, bottom-quartile matched share, or an uncommon mod; candidates keep every row."""
    from ms2_select import annotate_rows
    rows = [
        _row("SHORTPEK", 0.6, scan=1),
        _row("LONGPEPTIDEK", 0.6, "LONGPEPTIDEK-[UNIMOD:2]", scan=2),
        _row("LONGPEPTIDEK", 0.6, "LONGPEPTIDEM[UNIMOD:35]K", scan=3),
        _row("LONGPEPTIDEK", 0.6, "[UNIMOD:737]-LONGPEPTIDEK[UNIMOD:737]", scan=4),
    ] + [_row("ANOTHERLONGPEPK", 0.7 + 0.01 * k, scan=10 + k) for k in range(5)] + [
        _row("ANOTHERLONGPEPK", 0.01, scan=20)]
    out = annotate_rows(rows, "comet")
    assert len(out) == len(rows)
    assert {r["scan"] for r in out if r["hard"]} == {1, 2, 20}
    assert all(r["source_engine"] == "comet" and r["agree_comet"] is None and r["n_engines"] is None for r in out)


def test_hard_requires_agree():
    """A label needs >= 2 agreeing non-andes votes, a hard row >= 3; andes never counts, MS-GF+ sources neither."""
    from ms2_select import annotate_rows, finalize_rows
    rows = annotate_rows([
        _row("SHORTPEK", 0.6, scan=1),            # hard; Comet, MSFragger, Sage -> kept (3)
        _row("SHORTPEK", 0.6, scan=2),            # hard; Comet, Sage (+ andes) -> dropped (2 < 3)
        _row("LONGPEPTIDEK", 0.6, scan=3),        # easy; MSFragger + MaxQuant -> kept (2)
        _row("LONGPEPTIDEK", 0.6, scan=4),        # easy; andes + Comet -> dropped (1)
        _row("LONGPEPTIDEK", 0.6, scan=5),        # easy; Comet + submitter -> kept (2)
    ] + [_row("ANOTHERLONGPEPK", 0.7, scan=10 + k) for k in range(3)], "msgf+")
    votes = {"andes": {("r", 2): "SHORTPEK", ("r", 4): "LONGPEPTIDEK"},
             "comet": {("r", 1): "SHORTPEK", ("r", 2): "SHORTPEK", ("r", 4): "LONGPEPTIDEK", ("r", 5): "LONGPEPTIDEK"},
             "msfragger": {("r", 1): "SHORTPEK", ("r", 2): "SHORTPEKK", ("r", 3): "LONGPEPTIDEK"},
             "sage": {("r", 1): "SHORTPEK", ("r", 2): "SHORTPEK"},
             "maxquant": {("r", 3): "LONGPEPTIDEK"},
             "submitter": {("r", 5): "LONGPEPTIDEK"}}
    kept, st = finalize_rows([dict(r) for r in rows], votes)
    by = {r["scan"]: r for r in kept}
    assert set(by) == {1, 3, 5}, sorted(by)
    assert by[1]["n_engines"] == 3 and by[3]["n_engines"] == 2 and by[5]["n_engines"] == 2
    assert by[3]["agree_andes"] is False and by[1]["agree_maxquant"] is False
    assert st["dropped_hard"] == 1 and st["dropped_votes"] == 4
    # the source never votes: MsNet's own identification is only the informational msnet_agree column
    one = [_row("LONGPEPTIDEK", 0.7, scan=4), _row("ANOTHERLONGPEPK", 0.7, scan=9)]
    rows_m = annotate_rows([dict(r) for r in one], "msnet-consensus")
    kept, _ = finalize_rows([dict(r) for r in rows_m], {"comet": {("r", 4): "LONGPEPTIDEK"}})
    assert kept == []
    kept, _ = finalize_rows([dict(r) for r in rows_m], {"comet": {("r", 4): "LONGPEPTIDEK"},
                                                       "sage": {("r", 4): "LONGPEPTIDEK"}})
    assert [r["scan"] for r in kept] == [4] and kept[0]["n_engines"] == 2 and kept[0]["msnet_agree"] is True
    rows_c = annotate_rows([dict(r) for r in one], "comet")
    kept, _ = finalize_rows(rows_c, {"comet": {("r", 4): "LONGPEPTIDEK"}, "sage": {("r", 4): "LONGPEPTIDEK"},
                                     "msnet": {("r", 4): "LONGPEPTIDEK"}})
    assert kept[0]["msnet_agree"] is True and kept[0]["n_engines"] == 2
    assert all(r["msnet_agree"] is None for r in finalize_rows([dict(r) for r in rows], votes)[0])
    # engines that did not run leave their column null
    kept, _ = finalize_rows([dict(r) for r in rows], {"comet": {}, "msfragger": {}})
    assert kept == []


def test_finalize_script(tmp_path):
    from ms2_common import shard_table
    from ms2_select import annotate_rows
    rows = annotate_rows([_row("LONGPEPTIDEK", 0.6, scan=k) for k in (1, 2)], "comet")
    cand = tmp_path / "cand" / "cellA" / "D"
    cand.mkdir(parents=True)
    pq.write_table(shard_table(rows), cand / "part-msnet-0000.parquet")
    for e, scans in (("sage", [1, 2]), ("msfragger", [1])):
        pq.write_table(pa.table({"dataset": ["D"] * len(scans), "run": ["r"] * len(scans), "scan": scans,
                                 "peptidoform": ["LONGPEPTIDEK"] * len(scans)}), tmp_path / f"{e}.parquet")
    out = tmp_path / "final"
    run_script("ms2_finalize.py", "--candidates", tmp_path / "cand", "--out", out,
               "--votes", f"sage={tmp_path / 'sage.parquet'}", "--votes", f"msfragger={tmp_path / 'msfragger.parquet'}")
    t = pq.read_table(out / "cellA" / "D" / "part-msnet-0000.parquet")
    assert t.schema.metadata[b"andes.ms2_shard.finalized"] == b"true"
    assert t.column("scan").to_pylist() == [1] and t.column("agree_comet").to_pylist() == [None]
    assert t.column("n_engines").to_pylist() == [2]


def test_submitter_agreement_ignores_i_l_and_carbamidomethyl():
    from ms2_select import same_identification
    assert same_identification("PEPIC[UNIMOD:4]K", "PEPLCK")
    assert not same_identification("PEPICK", "PEPICR")
    assert not same_identification("PEPM[UNIMOD:35]K", "PEPMK")


def test_run_label_from_tags():
    from ms2_select import annotate_rows
    tmt = "[UNIMOD:737]-LONGPEPTIDEK[UNIMOD:737]"
    rows = [_row("LONGPEPTIDEK", 0.6, tmt, run="t", scan=1), _row("LONGPEPTIDEK", 0.6, tmt, run="t", scan=2),
            _row("LONGPEPTIDEK", 0.6, run="t", scan=3),
            _row("LONGPEPTIDEK", 0.6, "[UNIMOD:214]-LONGPEPTIDEK[UNIMOD:214]", run="i", scan=4),
            _row("LONGPEPTIDEK", 0.6, run="n", scan=5)]
    out = annotate_rows(rows, "comet")
    assert {r["scan"]: r["label"] for r in out} == {1: "TMT", 2: "TMT", 3: "TMT", 4: "iTRAQ", 5: "none"}


def test_filter_v2_electron_activations():
    from ms2_common import parse_filter_v2
    assert parse_filter_v2("FTMS + p NSI d Full ms2 445.12@etd50.00@hcd20.00 [110.00-905.00]") == ("EThcD", 20.0)
    act, nce = parse_filter_v2("ITMS + c NSI d Full ms2 445.12@etd70.00 [110.00-905.00]")
    assert act == "ETD" and np.isnan(nce)
    assert parse_filter_v2("FTMS + p NSI d Full ms2 500.00@hcd28.00 [100.00-1000.00]") == ("HCD-FTMS", 28.0)


def test_tof_mzml_without_filter_is_hcd_tof(tmp_path):
    """A timsconvert-style mzML (no filter string, timsTOF header) is HCD-TOF with NCE nan; a Q Exactive
    scan without a filter string is still skipped (Thermo behaviour unchanged)."""
    seq, z = "PEPTIDEKAR", 2
    mz, inten = perfect_spectrum(seq, [], z)
    pmz = precursor_mz(seq, [], z)
    mzml_dir = tmp_path / "mzml"
    mzml_dir.mkdir()
    write_mzml(mzml_dir / "run_t.mzML", [(7, "", pmz, z, mz, inten)])
    p = mzml_dir / "run_t.mzML"
    p.write_text(p.read_text().replace('accession="MS:1001911" name="Q Exactive"', 'accession="MS:1003005" name="timsTOF Pro"'))
    write_mzml(mzml_dir / "run_q.mzML", [(7, "", pmz, z, mz, inten)])
    selected = tmp_path / "selected.parquet"
    pq.write_table(pa.table({"run": ["run_t", "run_q"], "scan": [7, 7], "peptidoform": [seq, seq], "charge": [z, z],
                             "pep": [0.001, 0.001]}), selected)
    tasks = tmp_path / "tasks.tsv"
    tasks.write_text(f"task\tdataset\truns\tselected\traw_dir\n1\tPXD000000\trun_t|run_q\t{selected}\t{tmp_path}\n")
    shards = tmp_path / "shards"
    run_script("ms2_extract_raw.py", "--tasks", tasks, "--task-index", 1, "--out-dir", shards,
               "--work-dir", tmp_path / "work", "--mzml-dir", mzml_dir)
    rows = pq.read_table(shards / "PXD000000" / "part-00001.parquet").to_pylist()
    assert [r["run"] for r in rows] == ["run_t"]
    assert rows[0]["activation"] == "HCD-TOF" and np.isnan(rows[0]["nce"]) and rows[0]["instrument"] == "timsTOF Pro"
    runs = {r["run"]: r for r in json.load(open(shards / "PXD000000" / "part-00001.json"))["runs"]}
    assert runs["run_q"]["activation_skipped"] == 1 and runs["run_t"]["no_filter"] == 1


def _extract(tmp_path, spectra_by_run, seq, z, extra=()):
    """Write one mzML per run (``spectra_by_run``: run -> list of write_mzml tuples), select scan 1 of every run
    for ``seq``, run ms2_extract_raw.py and return (rows by run, per-run json stats)."""
    mzml_dir = tmp_path / "mzml"
    mzml_dir.mkdir()
    for run, spectra in spectra_by_run.items():
        write_mzml(mzml_dir / f"{run}.mzML", spectra)
    runs = sorted(spectra_by_run)
    selected = tmp_path / "selected.parquet"
    pq.write_table(pa.table({"run": runs, "scan": [1] * len(runs), "peptidoform": [seq] * len(runs),
                             "charge": [z] * len(runs), "pep": [0.001] * len(runs)}), selected)
    tasks = tmp_path / "tasks.tsv"
    tasks.write_text(f"task\tdataset\truns\tselected\traw_dir\n1\tPXD000000\t{'|'.join(runs)}\t{selected}\t{tmp_path}\n")
    shards = tmp_path / "shards"
    run_script("ms2_extract_raw.py", "--tasks", tasks, "--task-index", 1, "--out-dir", shards,
               "--work-dir", tmp_path / "work", "--mzml-dir", mzml_dir, *extra)
    rows = {r["run"]: r for r in pq.read_table(shards / "PXD000000" / "part-00001.parquet").to_pylist()}
    stats = {r["run"]: r for r in json.load(open(shards / "PXD000000" / "part-00001.json"))["runs"]}
    return rows, stats


def test_uvpd_declared_by_sdrf(tmp_path):
    """PXD003109: a UVPD laser in a Q Exactive HCD cell records '@hcd1.00' with the beam-type CID term. The SDRF's
    dissociation method makes that run UVPD (NCE = the filter energy); an HCD run of the same SDRF stays HCD-FTMS."""
    seq, z = "PEPTIDEKAR", 2
    mz, inten = perfect_spectrum(seq, [], z)
    pmz = precursor_mz(seq, [], z)
    filt = f"FTMS + p NSI d Full ms2 {pmz:.2f}@hcd1.00 [100.00-2000.00]"
    sdrf = tmp_path / "x.sdrf.tsv"
    sdrf.write_text("source name\tcomment[data file]\tcomment[dissociation method]\n"
                    "s1\trun_u.raw\tNT=ultraviolet photodissociation;AC=MS:1003246\n"
                    "s2\trun_h.raw\tNT=beam-type collision-induced dissociation;AC=MS:1000422\n")
    rows, stats = _extract(tmp_path, {"run_u": [(1, filt, pmz, z, mz, inten)],
                                      "run_h": [(1, filt.replace("@hcd1.00", "@hcd25.00"), pmz, z, mz, inten)]},
                           seq, z, ["--sdrf", sdrf])
    assert rows["run_u"]["activation"] == "UVPD" and rows["run_u"]["nce"] == 1.0
    assert rows["run_h"]["activation"] == "HCD-FTMS" and rows["run_h"]["nce"] == 25.0
    assert stats["run_u"]["uvpd"] == 1 and stats["run_h"]["uvpd"] == 0
    # Without the SDRF the same scan is HCD: the mzML alone cannot show this UVPD.
    (tmp_path / "nosdrf").mkdir()
    rows2, _ = _extract(tmp_path / "nosdrf", {"run_u": [(1, filt, pmz, z, mz, inten)]}, seq, z)
    assert rows2["run_u"]["activation"] == "HCD-FTMS"


def test_uvpd_from_scan_metadata(tmp_path):
    """A '@uvpd' filter step or the MS:1003246 activation term is UVPD on its own."""
    from ms2_common import ACTIVATIONS, SHARD_ONLY_ACTIVATIONS, activation_index, scan_activation, v2_activation
    from mzml_reader import read_scans
    act, nce = scan_activation("FTMS + p NSI d Full ms2 500.00@uvpd25.00 [100.00-1000.00]")
    assert act == "UVPD" and np.isnan(nce)
    assert scan_activation("FTMS + p NSI d Full ms2 500.00@hcd28.00 [100.00-1000.00]", ["UVPD"]) == ("UVPD", 28.0)
    assert scan_activation("FTMS + p NSI d Full ms2 500.00@hcd28.00 [100.00-1000.00]", ["HCD"]) == ("HCD-FTMS", 28.0)
    assert v2_activation("ultraviolet photodissociation", "FTMS", "Tribrid") == "UVPD"
    assert "UVPD" in SHARD_ONLY_ACTIVATIONS and "UVPD" not in ACTIVATIONS and activation_index("UVPD") == 0
    seq, z = "PEPTIDEKAR", 2
    mz, inten = perfect_spectrum(seq, [], z)
    p = tmp_path / "u.mzML"
    write_mzml(p, [(1, "FTMS + p NSI d Full ms2 500.00@hcd1.00 [100.00-1000.00]", 500.0, z, mz, inten)])
    assert next(read_scans(str(p)))["dissociation"] == ["HCD"]
    p.write_text(p.read_text().replace('accession="MS:1000422" name="beam-type collision-induced dissociation"',
                                       'accession="MS:1003246" name="ultraviolet photodissociation"'))
    assert next(read_scans(str(p)))["dissociation"] == ["UVPD"]


def test_adp_ribosyl_supported_and_annotated(tmp_path):
    """PXD023835: ADP-Ribosyl (UNIMOD:213, +541.06111 on S/R) parses, is supported (ms2_select / the corpus
    extract.py keep it) and an EThcD scan annotates c / z-dot ions carrying the intact ADP-ribose, as the
    engines' label peptidoforms do (no ADP-ribose loss channel)."""
    from ms2_common import MOD_MASS, N_ION_V2, parse_proforma, supported
    from ms2_select import peptidoform_key
    assert abs(MOD_MASS[213] - 541.06111) < 1e-5 and abs(MASSES[213] - 541.06111) < 1e-5
    pf = "SGAALS[UNIMOD:213]KR"
    seq, mods = parse_proforma(pf)
    assert (seq, mods) == ("SGAALSKR", [(6, 213)]) and supported(seq, mods) and peptidoform_key(pf) is not None
    z = 3
    mz, inten = ethcd_spectrum(seq, mods, z)
    pmz = precursor_mz(seq, mods, z)
    assert abs(pmz - precursor_mz(seq, [], z) - 541.06111 / z) < 1e-5
    rows, stats = _extract(tmp_path, {"run_a": [(1, f"FTMS + p NSI d Full ms2 {pmz:.2f}@etd50.00@hcd20.00 "
                                                    "[100.00-2000.00]", pmz, z, mz, inten)]}, pf, z)
    row = rows["run_a"]
    assert row["activation"] == "EThcD" and row["peptidoform"] == pf
    v = np.asarray(row["intensities"]).reshape(len(seq) - 1, N_ION_V2)
    # c ions after the modified S (site index >= 5) and z-dot ions holding it (site index <= 4) carry +541.06.
    assert (v[5:, 15] > 0).all() and (v[:5, 17] > 0).all()
    assert (v[:, 4:8] == -1).all(), "ADP-ribosyl S is not a phospho site: no H3PO4 loss channel"
    assert stats["run_a"]["precursor_mismatch"] == 0
