#!/usr/bin/env python3
"""Pool per-fraction glyco PIN files into one PIN for a single Percolator run.

Pooling is not an optimisation, it is a correctness requirement. One fraction of a
typical glyco run yields on the order of 0-2 decoy glycopeptides, so a per-fraction 1%
q-value is estimated from almost no data and swings between runs; differences measured
that way are noise. Search each file separately, pool here, run Percolator once.

Scan numbers collide across fractions, so each SpecId is prefixed with a per-file tag
(`f1_`, `f2_`, ...) derived from input order. Downstream evaluators recover the fraction
from that prefix.

Usage:
    pool_pins.py frac1.glyco.pin frac2.glyco.pin ... > pooled.pin

All inputs must share a header; a mismatch is an error rather than a silent
column-misalignment in the pooled file.
"""
import sys

pins = sys.argv[1:]
if not pins:
    sys.exit(__doc__)

header = None
rows = 0
for i, path in enumerate(pins, start=1):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n")
        if not hdr:
            sys.exit(f"error: {path} is empty")
        if header is None:
            header = hdr
            print(header)
        elif hdr != header:
            sys.exit(
                f"error: {path} has a different header from {pins[0]}.\n"
                "Pooling PINs with different feature columns would misalign them; "
                "re-run every fraction with the same andes version and flags."
            )
        for line in fh:
            if not line.strip():
                continue
            print(f"f{i}_{line.rstrip()}")
            rows += 1

print(f"pooled {rows} rows from {len(pins)} files", file=sys.stderr)
