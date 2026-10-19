#!/usr/bin/env python3
"""Build the andes-rt training table from quantms qpx PSM tables.

For each input ``psm.parquet`` keep target PSMs with PEP <= --max-pep, drop
peptidoforms carrying a modification outside ``mod_composition.UNIMOD``, and
reduce to one row per (dataset, run, peptidoform) with the exact median RT in minutes
and the PSM count.

Output: a directory of parquet parts, columns: dataset, run, peptidoform, sequence, rt_min, n_psm,
split (train/val/test by a hash of the stripped sequence, so every modified
form of a sequence lands in the same split), has_phospho, has_tmt.

Usage::

  python rt_extract.py --out rt_table/ DATASET=path/psm.parquet [...]
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import os

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod_composition import UNIMOD  # noqa: E402

def split_of(seq: str) -> str:
    """train / val / test by an md5 of the stripped sequence, I as L (same rule as ms2_common)."""
    h = int(hashlib.md5(seq.replace("I", "L").encode()).hexdigest()[:8], 16) % 100
    return "test" if h < 5 else "val" if h < 10 else "train"


def group_median(t: pa.Table) -> pa.Table:
    """Exact median RT and PSM count per (run, peptidoform), via dictionary codes and one sort."""
    run = pc.dictionary_encode(t["run_file_name"]).combine_chunks()
    form = pc.dictionary_encode(t["peptidoform"]).combine_chunks()
    seq = pc.dictionary_encode(t["sequence"]).combine_chunks()
    key = run.indices.to_numpy().astype(np.int64) * len(form.dictionary) + form.indices.to_numpy()
    rt = t["rt"].combine_chunks().to_numpy(zero_copy_only=False).astype(np.float64)
    order = np.lexsort((rt, key))
    key, rt = key[order], rt[order]
    seq_idx = seq.indices.to_numpy()[order]
    start = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    count = np.diff(np.r_[start, len(key)])
    lo = start + (count - 1) // 2
    hi = start + count // 2
    med = (rt[lo] + rt[hi]) / 2
    k = key[start]
    return pa.table({
        "run_file_name": pc.take(run.dictionary, pa.array(k // len(form.dictionary))),
        "peptidoform": pc.take(form.dictionary, pa.array(k % len(form.dictionary))),
        "sequence": pc.take(seq.dictionary, pa.array(seq_idx[start])),
        "rt_median": pa.array(med),
        "rt_count": pa.array(count),
    })


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, help="output directory (one parquet part per input)")
    p.add_argument("--max-pep", type=float, default=0.01)
    p.add_argument("inputs", nargs="+", help="DATASET=path")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    known = "|".join(str(k) for k in sorted(UNIMOD, reverse=True))
    total = 0
    for item in a.inputs:
        ds, path = item.split("=", 1)
        t = pq.read_table(path, columns=["peptidoform", "sequence", "posterior_error_probability", "is_decoy", "rt", "run_file_name"])
        n0 = t.num_rows
        keep = pc.and_(pc.invert(t["is_decoy"]), pc.less_equal(t["posterior_error_probability"], a.max_pep))
        keep = pc.and_(keep, pc.is_valid(t["rt"]))
        t = t.filter(keep).select(["run_file_name", "peptidoform", "sequence", "rt"])
        kept = t.num_rows
        g = group_median(t)
        del t
        stripped = pc.replace_substring_regex(g["peptidoform"], rf"\[UNIMOD:({known})\]", "")
        ok = pc.invert(pc.match_substring(stripped, "UNIMOD"))
        dropped = g.num_rows - pc.sum(ok).as_py()
        g = g.filter(ok)
        enc = pc.dictionary_encode(g["sequence"]).combine_chunks()
        splits = pa.array([split_of(s) for s in enc.dictionary.to_pylist()])
        out = pa.table(
            {
                "dataset": pa.array([ds] * g.num_rows),
                "run": g["run_file_name"],
                "peptidoform": g["peptidoform"],
                "sequence": g["sequence"],
                "rt_min": pc.divide(g["rt_median"], 60.0),
                "n_psm": g["rt_count"],
                "split": pc.take(splits, enc.indices),
                "has_phospho": pc.match_substring(g["peptidoform"], "UNIMOD:21]"),
                "has_tmt": pc.or_(pc.match_substring(g["peptidoform"], "UNIMOD:737]"),
                                  pc.match_substring(g["peptidoform"], "UNIMOD:2016]")),
            }
        )
        tag = hashlib.md5(path.encode()).hexdigest()[:8]
        pq.write_table(out, os.path.join(a.out, f"part-{ds}-{tag}.parquet"))
        total += out.num_rows
        print(f"{ds}: rows {n0} kept PSMs {kept} -> {out.num_rows} (run, peptidoform); "
              f"dropped {dropped} with unmapped mods; runs {pc.count_distinct(out['run']).as_py()}", flush=True)
    print("total", total, "rows")


if __name__ == "__main__":
    main()
