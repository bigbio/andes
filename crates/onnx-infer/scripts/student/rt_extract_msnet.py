#!/usr/bin/env python3
"""Build andes-rt training rows from MsNet PSM tables (the layout of ``rt_extract.py`` output).

Per input file: targets with PEP <= --max-pep, one row per (run, peptidoform) with the exact median RT in
minutes and the PSM count. Peptidoforms are written with named modifications in MsNet; they are mapped to
the built-in Unimod table and re-written as ProForma with ``[UNIMOD:id]`` tags (an unknown name drops the
row). Layouts: ``rt`` / ``run_file_name`` or ``retention_time`` / ``reference_file_name`` (seconds). A file
without a retention-time column is skipped and reported. The extra column ``instrument`` holds the
instrument model of the file's first PSM (from ``cv_params``), ``analyser`` its mass analyser.

Usage::

  python rt_extract_msnet.py --out rt_table_msnet/ DATASET=path/to/file.parquet [...]
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "ms2"))
from ms2_common import split_of, to_proforma  # noqa: E402
from ms2_extract_msnet import conditions, parse_named  # noqa: E402
from rt_extract import group_median  # noqa: E402


def first_present(names, have):
    for n in names:
        if n in have:
            return n
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True)
    p.add_argument("--max-pep", type=float, default=0.01)
    p.add_argument("inputs", nargs="+", help="DATASET=path")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for item in a.inputs:
        ds, path = item.split("=", 1)
        pf = pq.ParquetFile(path)
        names = set(pf.schema_arrow.names)
        c_rt = first_present(["rt", "retention_time"], names)
        c_run = first_present(["run_file_name", "reference_file_name"], names)
        if c_rt is None or c_run is None or not {"peptidoform", "sequence", "posterior_error_probability", "is_decoy"} <= names:
            print(f"{ds}: skipped, no usable retention time ({sorted(names)})", flush=True)
            continue
        inst = analyser = ""
        if "cv_params" in names:
            first = next(pf.iter_batches(batch_size=50, columns=["cv_params"])).to_pydict()["cv_params"]
            for cv in first:
                _, analyser, _, inst = conditions(cv)
                if inst:
                    break
        parts = []
        n0 = 0
        for batch in pf.iter_batches(batch_size=1 << 20, columns=["peptidoform", "sequence", "posterior_error_probability",
                                                                 "is_decoy", c_rt, c_run]):
            t = pa.Table.from_batches([batch])
            n0 += t.num_rows
            keep = pc.and_(pc.invert(pc.fill_null(t["is_decoy"], False)),
                           pc.fill_null(pc.less_equal(t["posterior_error_probability"], a.max_pep), False))
            keep = pc.and_(keep, pc.is_valid(t[c_rt]))
            t = t.filter(keep).select([c_run, "peptidoform", "sequence", c_rt]).rename_columns(
                ["run_file_name", "peptidoform", "sequence", "rt"])
            parts.append(t.cast(pa.schema([("run_file_name", pa.string()), ("peptidoform", pa.string()),
                                           ("sequence", pa.string()), ("rt", pa.float64())])))
        t = pa.concat_tables(parts)
        kept = t.num_rows
        g = group_median(t)
        del t, parts
        forms = pc.dictionary_encode(g["peptidoform"]).combine_chunks()
        mapped = []
        for name in forms.dictionary.to_pylist():
            r = parse_named(name)
            mapped.append(to_proforma(*r) if r else None)
        mapped_arr = pa.array(mapped, pa.string())
        pro = pc.take(mapped_arr, forms.indices)
        ok = pc.is_valid(pro)
        dropped = g.num_rows - pc.sum(ok).as_py()
        g = g.filter(ok)
        pro = pro.filter(ok) if hasattr(pro, "filter") else pc.filter(pro, ok)
        enc = pc.dictionary_encode(g["sequence"]).combine_chunks()
        splits = pa.array([split_of(s) for s in enc.dictionary.to_pylist()])
        out = pa.table({
            "dataset": pa.array([ds] * g.num_rows),
            "run": g["run_file_name"],
            "peptidoform": pro,
            "sequence": g["sequence"],
            "rt_min": pc.divide(g["rt_median"], 60.0),
            "n_psm": g["rt_count"],
            "split": pc.take(splits, enc.indices),
            "has_phospho": pc.match_substring(pro, "UNIMOD:21]"),
            "has_tmt": pc.or_(pc.match_substring(pro, "UNIMOD:737]"), pc.match_substring(pro, "UNIMOD:2016]")),
            "instrument": pa.array([inst] * g.num_rows),
            "analyser": pa.array([analyser] * g.num_rows),
        })
        tag = hashlib.md5(path.encode()).hexdigest()[:8]
        pq.write_table(out, os.path.join(a.out, f"part-{ds}-{tag}.parquet"))
        print(f"{ds}: rows {n0} kept PSMs {kept} -> {out.num_rows} (run, peptidoform); dropped {dropped} with unmapped mods; "
              f"runs {pc.count_distinct(out['run']).as_py()}; instrument {inst!r} {analyser!r}", flush=True)


if __name__ == "__main__":
    main()
