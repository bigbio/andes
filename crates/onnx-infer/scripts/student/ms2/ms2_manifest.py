#!/usr/bin/env python3
"""List the PSM tables available for andes-ms2 training, grouped by fragmentation.

Inputs are qpx ``psm.parquet`` tables given as ``DATASET=path`` or as bare paths
(the accession is then taken from the path). When a dataset has several tables,
the most recently modified one is kept unless ``--all-tables``.

Each dataset gets a group:
- ``HCD-LFQ``, ``HCD-TMT``, ``CID-TMT``, ``CID-ITMS-LFQ`` from ``--metadata``
  (TSV with ``dataset`` and ``group`` columns) or the built-in table below;
- otherwise ``HCD-TMT`` / ``unknown`` from a sample of its peptidoforms (TMT or not).
The group only drives reporting and the training mix: the activation of every
spectrum is read from its own scan header during extraction.

Benchmark accessions are never listed; ``--exclude`` adds more. The raw-file
directory of a dataset is ``--raw-template`` with ``{dataset}`` substituted.

Output TSV columns: dataset, group, phospho, psm_tables (comma-separated),
raw_dir, n_rows.

Usage::

  python ms2_manifest.py --out datasets.tsv --raw-template RAWROOT/{dataset} PSM_TABLES...
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import pyarrow.parquet as pq

# Never trained on: the andes benchmark datasets.
BENCHMARK = {"PXD070049", "PXD001819", "PXD007683", "PXD007653", "PXD005553", "PXD002047"}

KNOWN_GROUPS = {
    **{d: "HCD-LFQ" for d in (
        "PXD042233", "PXD005354", "PXD000612", "PXD002395", "PXD030406", "PXD009875",
        "PXD000396", "PXD036643", "PXD004940", "PXD014058", "PXD018450", "PMID23933261-FTMS",
        "PXD000269", "PXD013658", "PXD001305", "PXD005481", "PXD010150", "PXD019263",
        "PXD012255", "PXD003209")},
    **{d: "HCD-TMT" for d in ("PXD023662", "PXD011896", "PXD005235")},
    **{d: "CID-TMT" for d in ("MSV000085836", "PXD026581")},
    "PMID23933261-ITMS": "CID-ITMS-LFQ",
}
PHOSPHO = {"PXD000612", "PXD018450", "PXD012255"}

_ACC = re.compile(r"(PXD\d{6}|MSV\d{9}|PMID\d+(?:-[A-Z]+)?)")


def accession(path: str) -> str | None:
    hits = _ACC.findall(path)
    return hits[-1] if hits else None


def base_accession(ds: str) -> str:
    return ds.split("-", 1)[0]


def sample_tmt(path: str, n: int = 20000) -> float:
    pf = pq.ParquetFile(path)
    rows = []
    for b in pf.iter_batches(batch_size=n, columns=["peptidoform"]):
        rows = b.column(0).to_pylist()
        break
    if not rows:
        return 0.0
    return sum(("UNIMOD:737]" in r or "UNIMOD:2016]" in r) for r in rows) / len(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--raw-template", required=True, help="raw directory per dataset, e.g. ROOT/{dataset}/raw")
    p.add_argument("--metadata", help="TSV with dataset and group columns (overrides the built-in table)")
    p.add_argument("--exclude", nargs="*", default=[])
    p.add_argument("--all-tables", action="store_true")
    p.add_argument("tables", nargs="+", help="DATASET=psm.parquet or psm.parquet")
    a = p.parse_args()

    groups = dict(KNOWN_GROUPS)
    if a.metadata:
        import csv

        with open(a.metadata) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                groups[r["dataset"]] = r["group"]
    excluded = BENCHMARK | set(a.exclude)

    by_ds: dict[str, list[str]] = {}
    for item in a.tables:
        ds, path = item.split("=", 1) if "=" in item and not os.path.exists(item) else (accession(item), item)
        if ds is None:
            sys.exit(f"cannot infer the dataset accession of {item}; pass DATASET=path")
        if ds in excluded or base_accession(ds) in excluded:
            print(f"excluded: {ds} ({path})", file=sys.stderr)
            continue
        by_ds.setdefault(ds, []).append(path)

    with open(a.out, "w") as out:
        out.write("dataset\tgroup\tphospho\tpsm_tables\traw_dir\tn_rows\n")
        for ds in sorted(by_ds):
            paths = sorted(by_ds[ds], key=os.path.getmtime)
            if not a.all_tables:
                paths = paths[-1:]
            n = sum(pq.ParquetFile(x).metadata.num_rows for x in paths)
            g = groups.get(ds) or groups.get(base_accession(ds))
            if g is None:
                g = "HCD-TMT" if sample_tmt(paths[0]) > 0.5 else "unknown"
            ph = int(base_accession(ds) in PHOSPHO)
            raw = a.raw_template.format(dataset=ds, accession=base_accession(ds))
            out.write(f"{ds}\t{g}\t{ph}\t{','.join(paths)}\t{raw}\t{n}\n")
            print(f"{ds}\t{g}\tphospho={ph}\trows={n}\ttables={len(paths)}\traw={raw}"
                  f"{'' if os.path.isdir(raw) else ' (MISSING)'}")


if __name__ == "__main__":
    main()
