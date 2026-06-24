#!/usr/bin/env python3
"""Split a single-file model store into per-protocol Parquet partitions.

This is the *shipping* writer for `resources/models/protocol=<P>/models.parquet`.
It mirrors the logic of the Rust `split_store_by_protocol` (in
`crates/model-train/src/store/split.rs`), but writes the partitions through
pyarrow (parquet-cpp) so the high-entropy GBDT blob columns
(`rich_ion_model_bytes`, `frag_intensity_model_bytes`) compress as tightly as
the original single-file store. The Rust splitter remains the dev/parity tool;
its parquet-rs writer encodes these sparse Binary blob columns several times
larger than parquet-cpp, which is the only reason this script exists.

The split is lossless: every row that references a model (manifest, table,
source and stat rows) is routed into that model's protocol partition by
`model_id`, and the union of all partitions is row-for-row identical to the
input file. The Rust `ModelStore::open` on the resulting directory yields the
same in-memory store as opening the single file (the
`partition_parity` test gate proves this).

Usage:
    python scripts/split_store_by_protocol.py <single.parquet> <out_dir> \
        [--compression {zstd,snappy}]

Default compression is zstd (smallest); snappy matches the historical
single-file codec. Dictionary encoding is disabled for parity with the Rust
splitter (every model's GBDT blob is distinct, so a dictionary only adds an
index page).
"""

import argparse
import os
import shutil
import sys

import pyarrow.compute as pc
import pyarrow.parquet as pq


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("src", help="single-file model store (models_v1_own.parquet)")
    ap.add_argument("out_dir", help="output partitioned directory (e.g. resources/models)")
    ap.add_argument(
        "--compression",
        choices=["zstd", "snappy"],
        default="zstd",
        help="parquet codec for partition files (default: zstd)",
    )
    args = ap.parse_args()

    full = pq.ParquetFile(args.src).read()
    n = full.num_rows

    # 1. model_id -> protocol from manifest rows.
    manifest_mask = pc.equal(full["record_kind"], "manifest")
    m_mid = pc.filter(full["model_id"], manifest_mask).to_pylist()
    m_proto = pc.filter(full["protocol"], manifest_mask).to_pylist()
    model_protocol = {}
    for mid, proto in zip(m_mid, m_proto):
        model_protocol[mid] = proto if proto is not None else ""

    # 2. assign a protocol to every row by its model_id.
    all_mid = full["model_id"].to_pylist()
    missing = set()
    row_proto = []
    for mid in all_mid:
        if mid is None:
            row_proto.append(None)
            continue
        p = model_protocol.get(mid)
        if p is None:
            missing.add(mid)
            row_proto.append(None)
        else:
            row_proto.append(p)
    if missing:
        print(
            f"ERROR: model_id(s) {sorted(missing)} have no manifest row "
            "(cannot assign a protocol partition)",
            file=sys.stderr,
        )
        return 2

    import pyarrow as pa

    rp = pa.array(row_proto, type=pa.string())
    protocols = sorted({p for p in row_proto if p is not None})

    # 3. write one partition file per protocol.
    if os.path.isdir(args.out_dir):
        shutil.rmtree(args.out_dir)
    os.makedirs(args.out_dir, exist_ok=True)

    total = 0
    rows_written = 0
    print(f"Splitting {args.src} ({n} rows) -> {args.out_dir} [{args.compression}]")
    for p in protocols:
        sub = full.filter(pc.equal(rp, p))
        part_dir = os.path.join(args.out_dir, f"protocol={p}")
        os.makedirs(part_dir, exist_ok=True)
        part_path = os.path.join(part_dir, "models.parquet")
        pq.write_table(
            sub,
            part_path,
            compression=args.compression,
            use_dictionary=False,
        )
        sz = os.path.getsize(part_path)
        total += sz
        rows_written += sub.num_rows
        print(f"  protocol={p:<16} {sz:>12,} bytes  ({sub.num_rows} rows)")

    assert rows_written == n, f"row count mismatch: wrote {rows_written}, input {n}"
    print(f"TOTAL: {total / 1e6:.2f} MB ({rows_written} rows, lossless)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
