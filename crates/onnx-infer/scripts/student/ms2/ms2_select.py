#!/usr/bin/env python3
"""Choose the PSMs whose spectra andes-ms2 extracts, and cut them into array tasks.

For every dataset of ``ms2_manifest.py``'s TSV: keep target PSMs with
PEP <= --max-pep whose peptidoform parses and carries only modifications of the
built-in Unimod table (``mod_composition``); keep at most --max-per-precursor
PSMs per (peptidoform, charge) with the lowest PEP (ties broken at random, seeded);
write ``selected/<dataset>.parquet`` (dataset, run, scan, peptidoform, charge, pep).

--runs (TSV with ``dataset`` and ``run`` columns, e.g. the corpus ``manifest_v2.tsv``)
keeps only the listed runs; --max-psms-per-dataset keeps whole runs, in a seeded random
order, until the dataset reaches that many selected PSMs (fewer raw conversions).

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
from ms2_common import CARBAMIDOMETHYL, ITRAQ_MODS, TMT_MODS, label_of_mods, parse_proforma, supported  # noqa: E402

# Training-label rules (plan Global Constraints): a hard example is a peptide of at most
# HARD_MAX_LENGTH residues, a matched-intensity share in the bottom quartile of its extraction unit
# and activation, or a modification outside COMMON_MODS. A hard example is kept only when the
# submitter's identification of the same spectrum agrees.
HARD_MAX_LENGTH = 9
HARD_EXPLAINED_QUANTILE = 0.25
COMMON_MODS = frozenset({CARBAMIDOMETHYL, 35, 21}) | TMT_MODS | ITRAQ_MODS
# A run is labelled TMT / iTRAQ when at least this share of its PSMs carries the tag.
RUN_LABEL_SHARE = 0.5

CHARGE_COLS = ["precursor_charge", "charge"]
SCAN_COLS = ["scan", "scan_number"]
PEP_COLS = ["posterior_error_probability", "pep", "PEP"]
_SCAN = re.compile(r"(\d+)\s*$")
RAW_EXT = re.compile(r"\.(raw|mzml|mzml\.gz|mgf|d)$", re.I)


def run_stem(name: str) -> str:
    return RAW_EXT.sub("", os.path.basename(str(name)))


def _ident_key(peptidoform: str):
    p = parse_proforma(peptidoform)
    if p is None:
        return None
    seq, mods = p
    return seq.replace("I", "L"), tuple(sorted(m for m in mods if m[1] != CARBAMIDOMETHYL))


def same_identification(a: str, b: str) -> bool:
    """Whether two ProForma peptidoforms name the same identification: I and L are equal, and
    carbamidomethylation (a fixed mod some submitters leave out) is ignored."""
    ka = _ident_key(a)
    return ka is not None and ka == _ident_key(b)


def load_submitter_ids(path: str | None, dataset: str | None = None) -> dict | None:
    """``{(run, scan): peptidoform}`` from a parquet with ``run``, ``scan``, ``peptidoform`` (ProForma
    with Unimod accessions; optional ``dataset`` column filtered to ``dataset``); None without a path."""
    if not path:
        return None
    t = pq.read_table(path)
    if dataset is not None and "dataset" in t.column_names:
        t = t.filter(pc.equal(t["dataset"], dataset))
    runs = [run_stem(r) for r in t["run"].to_pylist()]
    return {(r, int(sc)): pf for r, sc, pf in zip(runs, t["scan"].to_pylist(), t["peptidoform"].to_pylist())}


def finalize_rows(rows: list[dict], submitter: dict | None):
    """Set ``label``, ``hard`` and ``agree`` on annotated shard rows and drop hard rows without agreement.

    ``label`` is the run's isobaric label (``RUN_LABEL_SHARE`` of its PSMs carry the tag). ``agree`` is
    True / False when ``submitter`` (``{(run, scan): peptidoform}``) identifies the spectrum, None otherwise.
    The explained-share quartile is taken per activation over ``rows``. Returns (kept rows, counts)."""
    st = {"rows": len(rows), "hard": 0, "hard_dropped": 0, "agree_true": 0, "agree_false": 0, "agree_null": 0}
    if not rows:
        return [], st
    by_run: dict = {}
    for r in rows:
        if "mods" not in r:
            r["mods"] = parse_proforma(r["peptidoform"])[1]
        by_run.setdefault(r["run"], []).append(label_of_mods(r["mods"]))
    run_label = {}
    for run, labs in by_run.items():
        top = max(("TMT", "iTRAQ"), key=labs.count)
        run_label[run] = top if labs.count(top) >= RUN_LABEL_SHARE * len(labs) else "none"
    cut = {}
    for act in {r["activation"] for r in rows}:
        ex = np.array([r["explained"] for r in rows if r["activation"] == act], np.float64)
        cut[act] = float(np.quantile(ex, HARD_EXPLAINED_QUANTILE))
    kept = []
    for r in rows:
        r["label"] = run_label[r["run"]]
        r["hard"] = (len(r["sequence"]) <= HARD_MAX_LENGTH or r["explained"] < cut[r["activation"]]
                     or any(uid not in COMMON_MODS for _, uid in r["mods"]))
        theirs = submitter.get((r["run"], int(r["scan"]))) if submitter else None
        r["agree"] = None if theirs is None else same_identification(r["peptidoform"], theirs)
        st["agree_null" if r["agree"] is None else "agree_true" if r["agree"] else "agree_false"] += 1
        st["hard"] += r["hard"]
        if r["hard"] and r["agree"] is not True:
            st["hard_dropped"] += 1
            continue
        del r["mods"]
        kept.append(r)
    return kept, st


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
    p.add_argument("--runs", help="TSV with dataset and run columns: keep only these runs")
    p.add_argument("--max-psms-per-dataset", type=int, default=0, help="0 = no cap; else whole runs at random")
    a = p.parse_args()
    keep_runs = None
    if a.runs:
        keep_runs = {}
        for r in csv.DictReader(open(a.runs), delimiter="\t"):
            keep_runs.setdefault(r["dataset"], set()).add(run_stem(r["run"]))
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
        if keep_runs is not None:
            df = df[df.run.isin(keep_runs.get(ds, set()))]
        if a.max_psms_per_dataset and len(df) > a.max_psms_per_dataset:
            sizes = df.groupby("run").size()
            order = sizes.index.to_numpy()[np.random.default_rng(a.seed).permutation(len(sizes))]
            total = np.cumsum(sizes.loc[order].to_numpy())
            n_keep = int(np.searchsorted(total, a.max_psms_per_dataset) + 1)
            df = df[df.run.isin(set(order[:n_keep]))]
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
