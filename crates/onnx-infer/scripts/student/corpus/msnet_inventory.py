#!/usr/bin/env python3
"""Inventory an MsNet collection (or any directory tree of PSM files) for the corpus catalog.

Walks ROOT and writes one row per file: relative path, size, format, dataset accession
(PXD / MSV / IPX / other id found in the path), and for parquet files the schema
fingerprint, row count, whether peak lists are attached (non-empty ``mz``/``intensity``
lists in a sample), which identification columns exist (score / q-value / PEP / decoy /
rt / scan / instrument / activation / NCE) and the MsNet-flat compatibility verdict
(``read_msnet_parquet`` needs: seq, charge, prec_mz, res_mod_pos, res_mod_delta,
nterm_delta, cterm_delta, mz, intensity).

A per-dataset roll-up (files, rows, GB, flat-compatible files, peak lists yes/no, the
union of informative columns) is written next to it, plus ``qc_commands.sh`` with one
``corpus_qc.py flat`` command per flat-compatible dataset, so the QC step is a plain
array job over that file.

Usage::

  python msnet_inventory.py --root MSNET_ROOT --out-dir INV [--sample-rows 200] [--max-files 0]
"""
from __future__ import annotations

import argparse
import collections
import csv
import os
import re
import sys

_ACC = re.compile(r"(PXD\d{6}|MSV\d{9}|IPX\d{10}|IPX\d+|JPST\d{6}|PDC\d{6}|PMID\d+)")
FLAT_COLS = {"seq", "charge", "prec_mz", "res_mod_pos", "res_mod_delta", "nterm_delta", "cterm_delta", "mz", "intensity"}
INFO = {
    "score": ("score", "hyperscore", "expect", "andromeda", "xcorr", "evalue", "e_value"),
    "qvalue": ("q_value", "qvalue", "q-value", "global_qvalue"),
    "pep": ("pep", "posterior_error_probability"),
    "decoy": ("is_decoy", "decoy", "label"),
    "rt": ("rt", "retention_time", "rt_sec", "rt_min"),
    "scan": ("scan", "scan_number", "spectrum_ref", "native_id", "title"),
    "run": ("run", "run_file_name", "raw_file", "file"),
    "instrument": ("instrument", "instrument_model"),
    "activation": ("activation", "fragmentation", "dissociation", "frag_method"),
    "nce": ("nce", "collision_energy"),
    "peptidoform": ("peptidoform", "modified_sequence", "proforma"),
}
FORMATS = [("parquet", r"\.parquet$"), ("mzML", r"\.mzml(\.gz)?$"), ("mgf", r"\.mgf(\.gz)?$"), ("raw", r"\.raw$"),
           ("mzIdentML", r"\.mzid(\.gz)?$"), ("mzTab", r"\.mztab(\.gz)?$"), ("pepXML", r"\.pep\.xml(\.gz)?$"),
           ("tsv", r"\.(tsv|txt|csv)(\.gz)?$"), ("msp", r"\.(msp|sptxt|blib|dlib|speclib)$"),
           ("archive", r"\.(zip|tar|tgz|gz|bz2|xz)$"), ("json", r"\.json$"), ("hdf5", r"\.(h5|hdf5)$")]


def fmt_of(name):
    n = name.lower()
    for k, rx in FORMATS:
        if re.search(rx, n):
            return k
    return "other"


def parquet_info(path, sample_rows):
    import pyarrow.parquet as pq

    pf = pq.ParquetFile(path)
    names = pf.schema_arrow.names
    lower = {c.lower(): c for c in names}
    info = {"rows": pf.metadata.num_rows, "columns": ",".join(names)}
    info["flat_compatible"] = int(FLAT_COLS <= set(names))
    for k, cands in INFO.items():
        info[f"has_{k}"] = next((lower[c] for c in cands if c in lower), "")
    peaks = ""
    mzc = next((c for c in names if c.lower() in ("mz", "mz_array", "mzs")), None)
    if mzc and pf.metadata.num_rows:
        try:
            b = next(pf.iter_batches(batch_size=sample_rows, columns=[mzc]))
            lens = [len(x) if x is not None else 0 for x in b.column(0).to_pylist()]
            peaks = f"{sum(1 for x in lens if x > 0)}/{len(lens)} non-empty, median {sorted(lens)[len(lens) // 2]}"
        except Exception as e:  # noqa: BLE001
            peaks = f"unreadable: {e}"
    info["peaks_sample"] = peaks
    # instrument / activation values if the columns exist (first batch only)
    for k in ("instrument", "activation", "nce"):
        c = info[f"has_{k}"]
        if c:
            try:
                b = next(pf.iter_batches(batch_size=sample_rows, columns=[c]))
                vals = collections.Counter(map(str, b.column(0).to_pylist()))
                info[f"{k}_values"] = ";".join(f"{v}:{n}" for v, n in vals.most_common(4))
            except Exception:  # noqa: BLE001
                pass
    return info


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample-rows", type=int, default=200)
    p.add_argument("--max-files", type=int, default=0)
    p.add_argument("--qc-script", default="corpus_qc.py", help="path written into qc_commands.sh")
    p.add_argument("--qc-out", default="QC", help="QC output root written into qc_commands.sh")
    a = p.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    cols = ["path", "dataset", "format", "bytes", "rows", "flat_compatible", "peaks_sample", "columns"] + \
        [f"has_{k}" for k in INFO] + ["instrument_values", "activation_values", "nce_values", "error"]
    per_ds = collections.defaultdict(lambda: collections.Counter())
    ds_cols = collections.defaultdict(set)
    ds_flats = collections.defaultdict(list)
    n = 0
    with open(os.path.join(a.out_dir, "files.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for dirpath, _, files in os.walk(a.root):
            for f in sorted(files):
                path = os.path.join(dirpath, f)
                rel = os.path.relpath(path, a.root)
                hits = _ACC.findall(rel)
                ds = hits[0] if hits else rel.split(os.sep)[0]
                row = {"path": rel, "dataset": ds, "format": fmt_of(f), "bytes": os.path.getsize(path)}
                if row["format"] == "parquet":
                    try:
                        row.update(parquet_info(path, a.sample_rows))
                    except Exception as e:  # noqa: BLE001
                        row["error"] = str(e)[:200]
                w.writerow(row)
                c = per_ds[ds]
                c["files"] += 1
                c[f"fmt_{row['format']}"] += 1
                c["bytes"] += row["bytes"]
                c["rows"] += int(row.get("rows") or 0)
                if row.get("flat_compatible"):
                    c["flat_files"] += 1
                    ds_flats[ds].append(path)
                if str(row.get("peaks_sample", "")).split("/")[0] not in ("", "0"):
                    c["files_with_peaks"] += 1
                for k in INFO:
                    if row.get(f"has_{k}"):
                        ds_cols[ds].add(k)
                n += 1
                if a.max_files and n >= a.max_files:
                    break
            if a.max_files and n >= a.max_files:
                break
    with open(os.path.join(a.out_dir, "datasets.tsv"), "w") as fh:
        fh.write("dataset\tfiles\tGB\trows\tflat_files\tfiles_with_peaks\tformats\tinformative_columns\n")
        for ds, c in sorted(per_ds.items()):
            fm = ",".join(f"{k[4:]}:{v}" for k, v in sorted(c.items()) if k.startswith("fmt_"))
            fh.write(f"{ds}\t{c['files']}\t{c['bytes'] / 1e9:.2f}\t{c['rows']}\t{c['flat_files']}\t"
                     f"{c['files_with_peaks']}\t{fm}\t{','.join(sorted(ds_cols[ds]))}\n")
    with open(os.path.join(a.out_dir, "qc_commands.sh"), "w") as fh:
        for ds, paths in sorted(ds_flats.items()):
            fh.write(f"python {a.qc_script} flat --dataset {ds} --out-dir {a.qc_out}/{ds} --workers 8 "
                     + " ".join(paths) + "\n")
    print(f"{n} files, {len(per_ds)} datasets, {sum(1 for d in ds_flats)} with flat-compatible parquet -> {a.out_dir}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
