#!/usr/bin/env python3
"""Extract benchmark metrics from GNU time output + MS-GF+ Percolator-pin output.

PR #23 removed mzIdentML output entirely; .pin is the only modern format. This
script counts native target / decoy rows directly from the .pin (column 2 = label,
{1, -1}). These counts are deterministic across runs (search produces the same
PSMs given the same inputs), so they form a stable correctness gate. Wall-time
and RSS come from GNU time -v.

For 1 % FDR PSM counts, run Percolator on the .pin separately — that's a
sensitivity gate, not a search-correctness gate, and Percolator's SVM has its
own stochasticity (seed 42 stabilises it). Keep the two gates separate.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_gnu_time(path: Path) -> tuple[str, str]:
    text = path.read_text(errors="replace")
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    cpu = re.search(r"Percent of CPU this job got: (\d+)", text)
    return (rss.group(1) if rss else "NA", cpu.group(1) if cpu else "NA")


def parse_pin(path: Path) -> tuple[int, int]:
    """Return (native_target_count, native_decoy_count) by counting label rows.

    A Percolator .pin is TSV with the second column being the label (1 = target,
    -1 = decoy). Header row is skipped. The file can be tens of millions of rows;
    streaming line-at-a-time keeps memory bounded.
    """
    targets = 0
    decoys = 0
    with path.open("r", encoding="utf-8", errors="replace") as f:
        next(f, None)  # header
        for line in f:
            cols = line.split("\t", 2)
            if len(cols) < 2:
                continue
            label = cols[1].strip()
            if label == "1":
                targets += 1
            elif label == "-1":
                decoys += 1
    return targets, decoys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--time", type=Path, required=True, help="GNU time -v output")
    ap.add_argument("--pin", type=Path, required=True, help="MS-GF+ Percolator .pin output")
    ap.add_argument("--wall", type=int, required=True, help="Wall-clock seconds (int)")
    ap.add_argument("--output", type=Path, required=True, help="Destination key=value file")
    args = ap.parse_args()

    rss_kb, cpu_pct = parse_gnu_time(args.time)
    targets, decoys = parse_pin(args.pin)

    lines = [
        "dataset=PXD001819",
        f"wall_time_sec={args.wall}",
        f"native_target_count={targets}",
        f"native_decoy_count={decoys}",
        f"peak_rss_kb={rss_kb}",
        f"cpu_percent={cpu_pct}",
    ]
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.output.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
