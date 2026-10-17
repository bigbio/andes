#!/usr/bin/env python3
"""Annotate the observed b/y intensities of selected PSMs, one array task at a time.

For each run of the task (``ms2_select.py`` ``tasks.tsv``):
1. find ``<run>.raw`` under the dataset's raw directory (or ``<run>.mzML[.gz]``
   under --mzml-dir, which skips conversion);
2. convert it with a raw-to-mzML converter (``--converter``: the command prefix,
   e.g. a container exec; ``--converter-args``: the argument template) into
   --work-dir;
3. stream the mzML, keeping only the selected scans; read NCE and activation
   from the scan's filter string, the instrument model from the file header;
4. check the precursor m/z against the peptidoform (isotope offsets -1..3 within
   --prec-tol-ppm), annotate the ``ms2_common.ION_TYPES`` channels at the
   analyser's tolerance and normalise;
5. delete the mzML. With --keep-mzml the task's scratch directory under
   --work-dir, with the converted mzML files, is kept and its path printed.

Spectra that are not HCD/CID, fail the precursor check or match fewer than
--min-matched plain b/y ions are dropped and counted. The task writes
``<out-dir>/<dataset>/part-<task>.json`` with the counts and, when every run of
the task succeeded, ``part-<task>.parquet`` (schema ``ms2_common.SHARD_COLUMNS``;
not written when no spectrum was kept). The JSON's ``complete`` field marks the
task done: a complete task is skipped unless --overwrite, so a re-submitted
array resumes. When any run fails (missing file, conversion error, unreadable
mzML), no parquet is written, a stale one is removed, ``complete`` is false and
the script exits with status 1; the next submission redoes the whole task.

Usage::

  python ms2_extract_raw.py --tasks SEL/tasks.tsv --task-index N --out-dir SHARDS \\
      --work-dir SCRATCH --converter "CONTAINER_EXEC CONVERTER"
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import annotate, parse_filter, parse_proforma, precursor_mz, shard_table  # noqa: E402
from mzml_reader import instrument_name, read_scans  # noqa: E402

DEFAULT_CONVERTER_ARGS = "-i {raw} -o {outdir} -f 1"


def tolerance(activation: str, a) -> tuple[float, bool]:
    """(tolerance, is_ppm) for the MS2 analyser of an activation name."""
    if activation.endswith("-ITMS"):
        return a.itms_tol_da, False
    return a.ftms_tol_ppm, True


_INDEX: dict = {}


def find_file(root: str, stem: str, exts) -> str | None:
    """``<root>/**/<stem><ext>`` for the first matching extension (case-insensitive)."""
    key = (root, tuple(exts))
    if key not in _INDEX:
        idx = {}
        for dirpath, _, files in os.walk(root):
            for f in files:
                low = f.lower()
                for e in exts:
                    if low.endswith(e.lower()):
                        idx.setdefault((f[: -len(e)].lower(), e.lower()), os.path.join(dirpath, f))
        _INDEX[key] = idx
    for e in exts:
        hit = _INDEX[key].get((stem.lower(), e.lower()))
        if hit:
            return hit
    return None


def convert(raw: str, work: str, a) -> str:
    outdir = os.path.join(work, "mzml")
    os.makedirs(outdir, exist_ok=True)
    cmd = shlex.split(a.converter) + shlex.split(a.converter_args.format(raw=raw, outdir=outdir))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, timeout=a.converter_timeout)
    stem = os.path.splitext(os.path.basename(raw))[0]
    for ext in (".mzML", ".mzML.gz"):
        p = os.path.join(outdir, stem + ext)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"the converter wrote no mzML for {raw} in {outdir}")


def prec_ok(seq, mods, z, mz_obs, tol_ppm) -> bool:
    if not np.isfinite(mz_obs):
        return True
    theo = precursor_mz(seq, mods, z)
    return any(abs(mz_obs - (theo + k * 1.0033548378 / z)) <= tol_ppm * 1e-6 * theo for k in range(-1, 4))


def process_run(mzml: str, psms, dataset: str, run: str, a, stats) -> list[dict]:
    want = {int(s): r for s, r in zip(psms.scan, psms.itertuples(index=False))}
    instr = instrument_name(mzml) or ""
    out = []
    for spec in read_scans(mzml, set(want)):
        r = want[spec["scan"]]
        stats["found"] += 1
        act, nce = parse_filter(spec["filter"])
        if act == "unknown":
            stats["activation_skipped"] += 1
            continue
        if not np.isfinite(nce) and np.isfinite(spec["collision_energy"]):
            nce = spec["collision_energy"]
        seq, mods = parse_proforma(r.peptidoform)
        if not prec_ok(seq, mods, int(r.charge), spec["precursor_mz"], a.prec_tol_ppm):
            stats["precursor_mismatch"] += 1
            continue
        tol, ppm = tolerance(act, a)
        vals, n_base, explained = annotate(seq, mods, int(r.charge), spec["mz"], spec["intensity"], tol, ppm)
        if vals is None or n_base < a.min_matched:
            stats["too_few_ions"] += 1
            continue
        stats["kept"] += 1
        out.append({
            "dataset": dataset, "run": run, "scan": spec["scan"], "peptidoform": r.peptidoform,
            "sequence": seq, "charge": int(r.charge), "nce": nce, "activation": act,
            "instrument": instr, "pep": float(r.pep), "n_matched": n_base, "explained": explained,
            "intensities": vals,
        })
    stats["not_in_file"] = len(want) - stats["found"]
    return out


def task_complete(summary: str) -> bool:
    """Whether ``part-<task>.json`` records a task whose runs all succeeded."""
    try:
        with open(summary) as fh:
            return json.load(fh).get("complete") is True
    except (OSError, ValueError):
        return False


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tasks", required=True)
    p.add_argument("--task-index", type=int, required=True, help="1-based row of tasks.tsv (SLURM_ARRAY_TASK_ID)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--work-dir", required=True, help="scratch for mzML files")
    p.add_argument("--converter", help="command prefix that runs the raw-to-mzML converter (required unless --mzml-dir)")
    p.add_argument("--converter-args", default=DEFAULT_CONVERTER_ARGS, help="argument template ({raw}, {outdir})")
    p.add_argument("--converter-timeout", type=int, default=3600)
    p.add_argument("--mzml-dir", help="use <run>.mzML[.gz] from here instead of converting raw files")
    p.add_argument("--keep-mzml", action="store_true")
    p.add_argument("--ftms-tol-ppm", type=float, default=20.0)
    p.add_argument("--itms-tol-da", type=float, default=0.5)
    p.add_argument("--prec-tol-ppm", type=float, default=20.0)
    p.add_argument("--min-matched", type=int, default=4)
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args()
    if not a.converter and not a.mzml_dir:
        p.error("--converter is required unless --mzml-dir is given")

    task = next((r for r in csv.DictReader(open(a.tasks), delimiter="\t") if int(r["task"]) == a.task_index), None)
    if task is None:
        sys.exit(f"task {a.task_index} not in {a.tasks}")
    ds = task["dataset"]
    out_dir = os.path.join(a.out_dir, ds)
    os.makedirs(out_dir, exist_ok=True)
    part = os.path.join(out_dir, f"part-{a.task_index:05d}.parquet")
    summary = part[:-8] + ".json"
    if task_complete(summary) and not a.overwrite:
        print(f"{summary} marks the task complete; skipping")
        return
    runs = task["runs"].split("|")
    sel = pq.read_table(task["selected"], filters=[("run", "in", runs)]).to_pandas()
    work = os.path.join(a.work_dir, f"ms2-{ds}-{a.task_index}-{os.getpid()}")
    os.makedirs(work, exist_ok=True)
    t0 = time.time()
    rows, log = [], []
    try:
        for run, psms in sel.groupby("run", sort=True):
            st = {"run": run, "selected": len(psms), "found": 0, "activation_skipped": 0,
                  "precursor_mismatch": 0, "too_few_ions": 0, "kept": 0, "not_in_file": 0}
            t1 = time.time()
            try:
                if a.mzml_dir:
                    mzml = find_file(a.mzml_dir, run, [".mzML.gz", ".mzML"])
                    if mzml is None:
                        raise FileNotFoundError(f"no mzML for {run} under {a.mzml_dir}")
                else:
                    raw = find_file(task["raw_dir"], run, [".raw"])
                    if raw is None:
                        raise FileNotFoundError(f"no raw file for {run} under {task['raw_dir']}")
                    mzml = convert(raw, work, a)
                rows += process_run(mzml, psms, ds, run, a, st)
                if not a.mzml_dir and not a.keep_mzml:
                    os.remove(mzml)
            except (FileNotFoundError, subprocess.SubprocessError, ValueError, OSError) as e:
                st["error"] = f"{type(e).__name__}: {e}"
            st["seconds"] = round(time.time() - t1, 1)
            log.append(st)
            print(json.dumps(st), flush=True)
    finally:
        if a.keep_mzml and not a.mzml_dir:
            print(f"kept the converted mzML files in {work}", flush=True)
        else:
            shutil.rmtree(work, ignore_errors=True)
    tot = {k: sum(s.get(k, 0) for s in log) for k in
           ("selected", "found", "activation_skipped", "precursor_mismatch", "too_few_ions", "kept")}
    tot["runs"], tot["runs_failed"] = len(log), sum("error" in s for s in log)
    complete = tot["runs_failed"] == 0
    if complete and rows:
        pq.write_table(shard_table(rows), part + ".tmp", compression="zstd")
        os.replace(part + ".tmp", part)
    elif os.path.exists(part):
        # Either a failed task (its part would be missing the failed runs) or a
        # complete task that kept nothing: no part may outlive the summary.
        os.remove(part)
    tot["seconds"] = round(time.time() - t0, 1)
    with open(summary + ".tmp", "w") as fh:
        json.dump({"task": a.task_index, "dataset": ds, "complete": complete, "total": tot, "runs": log},
                  fh, indent=1)
    os.replace(summary + ".tmp", summary)
    print(f"task {a.task_index} {ds}: {json.dumps(tot)}", flush=True)
    if not complete:
        sys.exit(f"task {a.task_index} {ds}: {tot['runs_failed']} of {tot['runs']} runs failed; "
                 "no part written, resubmit to retry")


if __name__ == "__main__":
    main()
