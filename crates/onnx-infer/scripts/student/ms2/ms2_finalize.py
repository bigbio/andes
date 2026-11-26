#!/usr/bin/env python3
"""Turn candidate andes-ms2 v2 shards into training shards with the multi-engine label rule.

Candidate shards (``andes.ms2_shard.finalized`` = ``false``) come from the extractors with the agreement
columns null. Each ``--votes ENGINE=PATH`` gives one source's confident identifications (ENGINE is andes,
comet, msfragger, sage, maxquant or submitter, or msnet for the informational ``msnet_agree`` column, which
never counts; PATH a parquet with ``run``, ``scan``, ``peptidoform`` as ProForma with Unimod
accessions, and optionally ``dataset``; for an engine, its own Percolator at q <= 0.01 on the same FASTA,
modifications and enzyme). A source without ``--votes`` did not run: its column stays null and gives no vote.
``ms2_select.finalize_rows`` sets the columns and keeps a row when at least two independent non-andes
identifications agree (three for a hard row).

``hard`` is taken from the candidate part (``--hard-scope part``: the extractor's quartile over the runs of its task),
or set again before the label rule over every candidate part of the dataset (``--hard-scope unit``: the
matched-intensity quartile per extraction unit = (dataset, training cell ``ms2_common.row_cell``) and activation,
the cross-run rule of the multi-run extraction tasks; needed when each part holds one run). ``--hla-runs`` (TSV
``dataset<TAB>run``) names the HLA runs for the unit cell.

Every ``part-*.parquet`` under --candidates is written to the same relative path under --out
(``andes.ms2_shard.finalized`` = ``true``; a part that keeps no row is not written), with the counts in
``<out>/finalize.json``.

Usage::

  python ms2_finalize.py --candidates CAND --out SHARDS --votes comet=comet.parquet --votes msfragger=msf.parquet
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import parse_proforma, row_cell, shard_table  # noqa: E402
from ms2_mix import load_hla_runs  # noqa: E402
from ms2_select import INFO_SOURCES, VOTE_SOURCES, explained_cuts, finalize_rows, is_hard, load_votes, run_stem  # noqa: E402

SOURCES = VOTE_SOURCES + INFO_SOURCES
UNIT_COLS = ["dataset", "run", "label", "activation", "instrument", "explained"]


def unit_key(ds, run, label, activation, instrument, hla_runs):
    return ds, row_cell(label, activation, instrument, (ds, run_stem(run)) in hla_runs), activation


def unit_cuts(files, hla_runs) -> dict:
    """``{(dataset, cell, activation): explained cut}`` over the candidate rows of ``files``."""
    keys, ex = [], []
    for f in files:
        d = pq.read_table(f, columns=UNIT_COLS).to_pydict()
        memo: dict = {}
        for ds, run, lab, act, ins, e in zip(*(d[c] for c in UNIT_COLS)):
            k = (ds, run, lab, act, ins)
            if k not in memo:
                memo[k] = unit_key(ds, run, lab, act, ins, hla_runs)
            keys.append(memo[k])
            ex.append(e)
    return explained_cuts(keys, ex)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--candidates", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--votes", action="append", default=[], help="ENGINE=PATH (repeatable)")
    p.add_argument("--hard-scope", choices=("part", "unit"), default="part")
    p.add_argument("--hla-runs", help="TSV dataset, run of the HLA runs (cell of --hard-scope unit)")
    a = p.parse_args()
    paths = {}
    for v in a.votes:
        e, _, path = v.partition("=")
        if e not in SOURCES or not path:
            p.error(f"--votes {v!r}: expected one of {SOURCES}=PATH")
        paths[e] = path
    total: dict = {}
    cache: dict = {}
    files = sorted(glob.glob(os.path.join(a.candidates, "**", "part-*.parquet"), recursive=True))
    hla_runs = load_hla_runs(a.hla_runs)
    cuts = unit_cuts(files, hla_runs) if a.hard_scope == "unit" else None
    for f in files:
        t = pq.read_table(f)
        rows = t.to_pylist()
        if not rows:
            continue
        ds = rows[0]["dataset"]
        if ds not in cache:
            cache = {ds: {e: load_votes(paths.get(e), ds) for e in SOURCES}}
        for r in rows:
            r["intensities"] = np.asarray(r["intensities"], np.float32)
            if cuts is not None:
                k = unit_key(r["dataset"], r["run"], r["label"], r["activation"], r["instrument"], hla_runs)
                r["hard"] = is_hard(r["sequence"], r["explained"], parse_proforma(r["peptidoform"])[1], cuts[k])
        kept, st = finalize_rows(rows, cache[ds])
        for k, v in st.items():
            total[k] = total.get(k, 0) + v
        if kept:
            out = os.path.join(a.out, os.path.relpath(f, a.candidates))
            os.makedirs(os.path.dirname(out), exist_ok=True)
            pq.write_table(shard_table(kept, finalized=True), out + ".tmp", compression="zstd")
            os.replace(out + ".tmp", out)
        print(f"{f}: {json.dumps(st)}", flush=True)
    total["votes"] = paths
    total["hard_scope"] = a.hard_scope
    if cuts is not None:
        total["hard_cuts"] = {"|".join(map(str, k)): round(v, 6) for k, v in sorted(cuts.items())}
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "finalize.json"), "w") as fh:
        json.dump(total, fh, indent=1)
    print(json.dumps(total))


if __name__ == "__main__":
    main()
