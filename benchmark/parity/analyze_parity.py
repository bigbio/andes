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


def run_self_tests() -> int:
    """Stub for now; tests get appended as components are implemented."""
    print("Self-tests: 0 ran, 0 failed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
