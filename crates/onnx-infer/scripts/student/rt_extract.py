#!/usr/bin/env python3
"""Build the andes-rt training table from quantms qpx PSM tables.

For each input ``psm.parquet`` keep target PSMs with PEP <= --max-pep, drop
peptidoforms carrying a modification outside ``mod_composition.UNIMOD``, and
reduce to one row per (dataset, run, peptidoform) with the median RT in minutes
and the PSM count.

Output parquet columns: dataset, run, peptidoform, sequence, rt_min, n_psm,
split (train/val/test by a hash of the stripped sequence, so every modified
form of a sequence lands in the same split), has_phospho, has_tmt.

Usage::

  python rt_extract.py --out rt_table.parquet DATASET=path/psm.parquet [...]
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
import os

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod_composition import UNIMOD  # noqa: E402

TOKEN = re.compile(r"\[UNIMOD:(\d+)\]")


def split_of(seq: str) -> str:
    h = int(hashlib.md5(seq.encode()).hexdigest()[:8], 16) % 100
    return "test" if h < 5 else "val" if h < 10 else "train"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True)
    p.add_argument("--max-pep", type=float, default=0.01)
    p.add_argument("inputs", nargs="+", help="DATASET=path")
    a = p.parse_args()
    parts = []
    for item in a.inputs:
        ds, path = item.split("=", 1)
        t = pq.read_table(path, columns=["peptidoform", "sequence", "posterior_error_probability", "is_decoy", "rt", "run_file_name"])
        n0 = t.num_rows
        keep = pc.and_(pc.invert(t["is_decoy"]), pc.less_equal(t["posterior_error_probability"], a.max_pep))
        keep = pc.and_(keep, pc.is_valid(t["rt"]))
        t = t.filter(keep)
        g = t.group_by(["run_file_name", "peptidoform", "sequence"]).aggregate([("rt", "approximate_median"), ("rt", "count")])
        forms = g["peptidoform"].to_pylist()
        ok_form = {}
        for f in set(forms):
            ok_form[f] = all(int(m) in UNIMOD for m in TOKEN.findall(f))
        mask = np.array([ok_form[f] for f in forms], dtype=bool)
        dropped = len(forms) - int(mask.sum())
        g = g.filter(pa.array(mask))
        seqs = g["sequence"].to_pylist()
        fs = g["peptidoform"].to_pylist()
        sp = {s: split_of(s) for s in set(seqs)}
        out = pa.table(
            {
                "dataset": pa.array([ds] * g.num_rows),
                "run": g["run_file_name"],
                "peptidoform": g["peptidoform"],
                "sequence": g["sequence"],
                "rt_min": pc.divide(pc.cast(g["rt_approximate_median"], pa.float64()), 60.0),
                "n_psm": g["rt_count"],
                "split": pa.array([sp[s] for s in seqs]),
                "has_phospho": pa.array(["UNIMOD:21]" in f for f in fs]),
                "has_tmt": pa.array(["UNIMOD:737]" in f or "UNIMOD:2016]" in f for f in fs]),
            }
        )
        print(f"{ds}: rows {n0} kept PSMs {t.num_rows} -> {out.num_rows} (run, peptidoform); "
              f"dropped {dropped} with unmapped mods; runs {len(set(out['run'].to_pylist()))}", flush=True)
        parts.append(out)
    allt = pa.concat_tables(parts)
    pq.write_table(allt, a.out)
    print("total", allt.num_rows, "rows;", len(set(allt["peptidoform"].to_pylist())), "peptidoforms")


if __name__ == "__main__":
    main()
