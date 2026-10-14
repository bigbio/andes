#!/usr/bin/env python3
"""Choose the PSMs whose spectra andes-ms2 extracts, and cut them into array tasks.

For every dataset of ``ms2_manifest.py``'s TSV: keep target PSMs with
PEP <= --max-pep whose peptidoform parses and carries only modifications of the
built-in Unimod table (``mod_composition``); keep at most --max-per-precursor
PSMs per (peptidoform, charge) with the lowest PEP (ties broken at random, seeded);
write ``selected/<dataset>.parquet`` (dataset, run, scan, peptidoform, charge, pep).

Runs are then grouped into tasks of about --psms-per-task selected PSMs (whole
runs only), written to ``tasks.tsv`` (task, dataset, raw_dir, selected, runs
separated by "|"); one task is one SLURM array element of ``ms2_extract_raw.py``.

Usage::

  python ms2_select.py --manifest datasets.tsv --out-dir SEL [--datasets PXD...]
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import parse_proforma, supported  # noqa: E402

CHARGE_COLS = ["precursor_charge", "charge"]
SCAN_COLS = ["scan", "scan_number"]
PEP_COLS = ["posterior_error_probability", "pep", "PEP"]
_SCAN = re.compile(r"(\d+)\s*$")
RAW_EXT = re.compile(r"\.(raw|mzml|mzml\.gz|mgf|d)$", re.I)


def run_stem(name: str) -> str:
    return RAW_EXT.sub("", os.path.basename(str(name)))


def pick(names, cands, what):
    for c in cands:
        if c in names:
            return c
    raise SystemExit(f"no {what} column among {cands} (have {names})")


def scan_number(v) -> int:
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (list, tuple)):
        v = v[0] if v else ""
    m = _SCAN.search(str(v)) or re.search(r"scan=(\d+)", str(v))
    return int(m.group(1)) if m else -1


def load(paths, max_pep):
    parts = []
    for path in paths:
        names = pq.ParquetFile(path).schema_arrow.names
        cc, sc, pc_ = pick(names, CHARGE_COLS, "charge"), pick(names, SCAN_COLS, "scan"), pick(names, PEP_COLS, "PEP")
        t = pq.read_table(path, columns=["peptidoform", "is_decoy", "run_file_name", cc, sc, pc_])
        keep = pc.and_(pc.invert(pc.fill_null(t["is_decoy"], False)), pc.less_equal(t[pc_], max_pep))
        t = t.filter(keep)
        scan = t[sc]
        if not pa.types.is_integer(scan.type):
            scan = pa.array([scan_number(v) for v in scan.to_pylist()], pa.int64())
        parts.append(pa.table({
            "run": t["run_file_name"].cast(pa.string()),
            "scan": pc.cast(scan, pa.int64()),
            "peptidoform": t["peptidoform"].cast(pa.string()),
            "charge": pc.cast(t[cc], pa.int8()),
            "pep": pc.cast(t[pc_], pa.float32()),
        }))
    return pa.concat_tables(parts).to_pandas()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--datasets", nargs="*", help="subset of the manifest (default: all)")
    p.add_argument("--max-pep", type=float, default=0.01)
    p.add_argument("--max-per-precursor", type=int, default=4)
    p.add_argument("--psms-per-task", type=int, default=150000)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    os.makedirs(os.path.join(a.out_dir, "selected"), exist_ok=True)
    rows = list(csv.DictReader(open(a.manifest), delimiter="\t"))
    if a.datasets:
        rows = [r for r in rows if r["dataset"] in set(a.datasets)]
    tasks = []
    for r in rows:
        ds = r["dataset"]
        df = load(r["psm_tables"].split(","), a.max_pep)
        n0 = len(df)
        forms = df.peptidoform.unique()
        ok = {f for f in forms if (q := parse_proforma(f)) is not None and supported(*q)}
        df = df[df.peptidoform.isin(ok) & (df.charge > 0) & (df.scan > 0)]
        n1 = len(df)
        rng = np.random.default_rng(a.seed)
        df = df.assign(_r=rng.random(len(df))).sort_values(["pep", "_r"])
        df = df.groupby(["peptidoform", "charge"], sort=False).head(a.max_per_precursor)
        df = df.drop(columns="_r").drop_duplicates(["run", "scan"])
        df["run"] = df.run.map(run_stem)
        df.insert(0, "dataset", ds)
        df = df.sort_values(["run", "scan"]).reset_index(drop=True)
        sel = os.path.join(a.out_dir, "selected", f"{ds}.parquet")
        pq.write_table(pa.Table.from_pandas(df, preserve_index=False), sel)
        per_run = df.groupby("run", sort=True).size()
        chunk, size = [], 0
        for run, n in per_run.items():
            chunk.append(run)
            size += n
            if size >= a.psms_per_task:
                tasks.append((ds, r["raw_dir"], sel, chunk))
                chunk, size = [], 0
        if chunk:
            tasks.append((ds, r["raw_dir"], sel, chunk))
        print(f"{ds}: {n0} confident target PSMs, {n1} with supported mods, {len(df)} selected "
              f"({df.peptidoform.nunique()} peptidoforms, {len(per_run)} runs)", flush=True)
    with open(os.path.join(a.out_dir, "tasks.tsv"), "w") as fh:
        fh.write("task\tdataset\traw_dir\tselected\truns\n")
        for i, (ds, raw, sel, runs) in enumerate(tasks, 1):
            fh.write(f"{i}\t{ds}\t{raw}\t{sel}\t{'|'.join(runs)}\n")
    print(f"{len(tasks)} tasks -> {os.path.join(a.out_dir, 'tasks.tsv')} (array 1-{len(tasks)})")


if __name__ == "__main__":
    main()
