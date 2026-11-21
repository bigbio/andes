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
from ms2_common import shard_table  # noqa: E402
from ms2_select import INFO_SOURCES, VOTE_SOURCES, finalize_rows, load_votes  # noqa: E402

SOURCES = VOTE_SOURCES + INFO_SOURCES


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--candidates", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--votes", action="append", default=[], help="ENGINE=PATH (repeatable)")
    a = p.parse_args()
    paths = {}
    for v in a.votes:
        e, _, path = v.partition("=")
        if e not in SOURCES or not path:
            p.error(f"--votes {v!r}: expected one of {SOURCES}=PATH")
        paths[e] = path
    total: dict = {}
    cache: dict = {}
    for f in sorted(glob.glob(os.path.join(a.candidates, "**", "part-*.parquet"), recursive=True)):
        t = pq.read_table(f)
        rows = t.to_pylist()
        if not rows:
            continue
        ds = rows[0]["dataset"]
        if ds not in cache:
            cache = {ds: {e: load_votes(paths.get(e), ds) for e in SOURCES}}
        for r in rows:
            r["intensities"] = np.asarray(r["intensities"], np.float32)
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
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "finalize.json"), "w") as fh:
        json.dump(total, fh, indent=1)
    print(json.dumps(total))


if __name__ == "__main__":
    main()
