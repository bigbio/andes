#!/usr/bin/env python3
"""Per-dataset balance of andes-ms2 training shards.

A raw back-fill of one large dataset can hold millions of rows where an extraction unit of the older corpus kept at most
``DATASET_CAP`` (``DATASET_CAP_SCARCE`` in a scarce cell), so training caps the rows of each (dataset, training cell)
unit at load time; every row stays on disk. A unit above its cap keeps a uniform random sample (``seed``) of its
eligible rows (``n_matched`` and ``explained`` at least the training minimums, an activation the model carries). A v1
shard has no ``run`` or ``label`` column: its rows are unlabelled and HLA only by their directory.

The training cell of a row is ``ms2_common.row_cell``. HLA comes from the run: a file under ``msnet/HLA/`` (an MsNet
HLA extraction unit), or a (dataset, run) listed in ``--hla-runs`` (a TSV ``dataset<TAB>run``, run without
extension).

Usage (the mix only, no training)::

  python ms2_mix.py --shards SHARDS/ --hla-runs SHARDS/hla_runs.tsv --out mix.tsv
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import ACTIVATIONS, CELLS, row_cell  # noqa: E402
from ms2_select import run_stem  # noqa: E402

DATASET_CAP = 60000
DATASET_CAP_SCARCE = 400000
SCARCE_CELLS = frozenset({"CID-ITMS", "TMT-CID-ITMS", "timsTOF", "EThcD-ETD", "HLA", "Astral"})
MIX_COLS = ["dataset", "run", "label", "activation", "instrument", "n_matched", "explained", "peptidoform"]


def shard_files(paths):
    out = []
    for p in paths:
        out += sorted(glob.glob(os.path.join(p, "**", "*.parquet"), recursive=True)) if os.path.isdir(p) else [p]
    return out


def load_hla_runs(path: str | None) -> set:
    if not path:
        return set()
    out = set()
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2 and c[0] and c[1] and c[0] != "dataset":
            out.add((c[0], run_stem(c[1])))
    return out


def in_hla_unit(path: str) -> bool:
    return f"{os.sep}msnet{os.sep}HLA{os.sep}" in os.path.abspath(path)


def row_cells(t, path: str, hla_runs: set) -> np.ndarray:
    """``CELLS`` index of every row of a shard table with the ``MIX_COLS`` columns."""
    unit_hla = in_hla_unit(path)
    keys = [t[c].to_pylist() for c in ("dataset", "run", "label", "activation", "instrument")]
    memo: dict = {}
    out = np.empty(t.num_rows, np.int16)
    for i, (ds, run, lab, act, ins) in enumerate(zip(*keys)):
        k = (ds, run, lab, act, ins)
        c = memo.get(k)
        if c is None:
            hla = unit_hla or (ds, run_stem(run)) in hla_runs
            c = memo[k] = CELLS.index(row_cell(lab, act, ins, hla))
        out[i] = c
    return out


def cap_of(cell: str, cap: int, cap_scarce: int) -> int:
    return cap_scarce if cell in SCARCE_CELLS else cap


def plan(files, hla_runs: set, min_matched: int, min_explained: float, cap: int, cap_scarce: int, seed: int):
    """Rows to train per file and the mix. Returns ({file: bool mask over its rows}, mix rows
    ``(dataset, cell, rows, eligible, kept, distinct eligible, distinct kept)``), the distinct counts
    being peptidoforms. A cap of 0 keeps every eligible row of its units."""
    masks, codes, datasets, forms, form_ids = {}, [], {}, {}, []
    for f in files:
        names = pq.ParquetFile(f).schema_arrow.names
        t = pq.read_table(f, columns=[c for c in MIX_COLS if c in names])
        for c in MIX_COLS:
            if c not in names:
                t = t.append_column(c, pa.nulls(t.num_rows, pa.string()))
        ok = pc.and_(pc.greater_equal(t["n_matched"], min_matched), pc.greater_equal(t["explained"], min_explained))
        ok = ok.to_numpy(zero_copy_only=False) & np.isin(np.asarray(t["activation"].to_pylist(), object), ACTIVATIONS)
        cells = row_cells(t, f, hla_runs)
        ds = np.asarray([datasets.setdefault(d, len(datasets)) for d in t["dataset"].to_pylist()], np.int64)
        codes.append(ds * len(CELLS) + cells)
        enc = pc.dictionary_encode(t["peptidoform"]).combine_chunks()
        gid = np.asarray([forms.setdefault(x, len(forms)) for x in enc.dictionary.to_pylist()], np.int64)
        form_ids.append(gid[enc.indices.to_numpy(zero_copy_only=False).astype(np.int64)] if len(gid) else
                        np.full(t.num_rows, -1, np.int64))
        masks[f] = ok
    code = np.concatenate(codes) if codes else np.zeros(0, np.int64)
    elig = np.concatenate([masks[f] for f in files]) if files else np.zeros(0, bool)
    form = np.concatenate(form_ids) if form_ids else np.zeros(0, np.int64)
    keep = elig.copy()
    rng = np.random.default_rng(seed)
    names = {v: k for k, v in datasets.items()}
    mix = []
    for u in np.unique(code):
        ds, cell = names[int(u) // len(CELLS)], CELLS[int(u) % len(CELLS)]
        rows = np.flatnonzero(code == u)
        e = rows[elig[rows]]
        lim = cap_of(cell, cap, cap_scarce)
        if lim and len(e) > lim:
            keep[e] = False
            keep[rng.choice(e, lim, replace=False)] = True
        kept = rows[keep[rows]]
        mix.append((ds, cell, len(rows), len(e), len(kept), len(np.unique(form[e])), len(np.unique(form[kept]))))
    off = 0
    for f in files:
        n = len(masks[f])
        masks[f] = keep[off: off + n]
        off += n
    return masks, mix


def write_mix(mix, path: str):
    with open(path, "w") as fh:
        fh.write("dataset\tcell\trows\teligible\tkept\tdistinct_eligible\tdistinct_kept\n")
        for r in sorted(mix):
            fh.write("\t".join(map(str, r)) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--shards", nargs="+", required=True)
    p.add_argument("--hla-runs")
    p.add_argument("--min-matched", type=int, default=4)
    p.add_argument("--min-explained", type=float, default=0.0)
    p.add_argument("--dataset-cap", type=int, default=DATASET_CAP)
    p.add_argument("--dataset-cap-scarce", type=int, default=DATASET_CAP_SCARCE)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", required=True, help="mix TSV (dataset, cell, rows, eligible, kept)")
    a = p.parse_args()
    _, mix = plan(shard_files(a.shards), load_hla_runs(a.hla_runs), a.min_matched, a.min_explained, a.dataset_cap,
                  a.dataset_cap_scarce, a.seed)
    write_mix(mix, a.out)
    tot = lambda i: sum(r[i] for r in mix)
    print(f"units {len(mix)} rows {tot(2)} eligible {tot(3)} kept {tot(4)}")


if __name__ == "__main__":
    main()
