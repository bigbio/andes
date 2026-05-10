#!/usr/bin/env python3
"""Parity analysis: stratified PIN diff for Java MS-GF+ vs Rust msgf-rust.

Reads java.pin and rust.pin from --java-pin / --rust-pin, emits a markdown
report covering ΔRawScore decomposition (Track A) and stratified flip
analysis (Track B) per the design doc at
docs/parity-analysis/2026-05-09-parity-analysis-design.md.

Usage:
  python3 benchmark/parity/analyze_parity.py \\
      --java-pin benchmark/results/PXD001819-parity/java.pin \\
      --rust-pin benchmark/results/PXD001819-parity/rust.pin \\
      --output docs/parity-analysis/reports/2026-05-09-parity-report.md

  python3 benchmark/parity/analyze_parity.py --self-test
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean, median, pstdev
from typing import Callable


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--java-pin", type=Path, help="Path to Java MS-GF+ PIN file")
    parser.add_argument("--rust-pin", type=Path, help="Path to Rust msgf-rust PIN file")
    parser.add_argument("--output", type=Path, help="Output markdown report path")
    parser.add_argument("--self-test", action="store_true", help="Run embedded unit tests and exit")
    args = parser.parse_args()

    if args.self_test:
        return run_self_tests()

    if not (args.java_pin and args.rust_pin and args.output):
        parser.error("--java-pin, --rust-pin, --output are required (or use --self-test)")

    print(f"Reading {args.java_pin}...")
    print(f"Reading {args.rust_pin}...")
    print(f"Will write to {args.output}")
    return 0


def parse_pin(path: Path) -> list[dict[str, str]]:
    """Parse a PIN file into a list of dicts keyed by column name.

    PIN files are tab-separated with a header row. Blank lines are skipped.
    All values are returned as raw strings; downstream callers cast as needed.
    """
    rows: list[dict[str, str]] = []
    with path.open() as f:
        header_line = f.readline().rstrip("\n")
        header = header_line.split("\t")
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                # Tolerate short rows (some PIN writers emit fewer columns for noise rows).
                # Pad with empty strings so the dict still has all keys.
                parts = parts + [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))
    return rows


def run_self_tests() -> int:
    """Stub for now; tests get appended as components are implemented."""
    print("Self-tests: 0 ran, 0 failed")
    return 0


# ── Tests for parse_pin ─────────────────────────────────────────────────

def _test_parse_pin_basic():
    import tempfile
    pin_text = (
        "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tmass\tRawScore\tDeNovoScore\t"
        "lnSpecEValue\tlnEValue\tisotope_error\tpeplen\tPeptide\n"
        "scan=5_5_1\t-1\t5\t1014.68\t1014.68\t1014.68\t-34\t12\t-8.77\t6.78\t0\t11\tK.SLKKISVIK.D\n"
        "scan=5_5_1\t-1\t5\t1014.68\t1015.69\t1014.68\t-34\t12\t-8.77\t6.78\t1\t11\tK.KPFIKIIR.D\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".pin", delete=False) as f:
        f.write(pin_text)
        path = Path(f.name)
    rows = parse_pin(path)
    assert len(rows) == 2, f"expected 2 rows, got {len(rows)}"
    assert rows[0]["ScanNr"] == "5"
    assert rows[0]["Label"] == "-1"
    assert rows[0]["RawScore"] == "-34"
    assert rows[0]["Peptide"] == "K.SLKKISVIK.D"
    assert rows[1]["isotope_error"] == "1"

def _test_parse_pin_skips_blank():
    import tempfile
    pin_text = (
        "SpecId\tLabel\tScanNr\tRawScore\tPeptide\n"
        "\n"
        "scan=1\t1\t1\t10\tK.AAA.B\n"
        "\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".pin", delete=False) as f:
        f.write(pin_text)
        path = Path(f.name)
    rows = parse_pin(path)
    assert len(rows) == 1
    assert rows[0]["RawScore"] == "10"


def run_self_tests() -> int:
    tests = [
        ("parse_pin basic", _test_parse_pin_basic),
        ("parse_pin skips blank lines", _test_parse_pin_skips_blank),
    ]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS: {name}")
        except AssertionError as e:
            print(f"  FAIL: {name}: {e}")
            failed += 1
    print(f"Self-tests: {len(tests)} ran, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
