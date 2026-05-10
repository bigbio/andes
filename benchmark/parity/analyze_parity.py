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


# Pre-compiled regexes for mod counting. Tolerant of trailing decimal
# variants (Java may write +15.9949 or +15.99491).
_OXIDATION_RE = re.compile(r"M\+15\.994\d*")
_CARBAMIDOMETHYL_RE = re.compile(r"C\+57\.0214?\d*")


def peptide_features(row: dict[str, str]) -> dict[str, object]:
    """Extract diagnostic features from a single PIN row.

    Returns a dict with keys: length, charge, n_oxidation,
    n_carbamidomethyl, iso_off, last_aa, pre_aa, is_decoy, score_bucket.
    All purely from the row — no spectrum or FASTA needed.
    """
    pep = row.get("Peptide", "")
    # Peptide is "X.RESIDUES.Y" form. Strip flanking residues.
    # Handle mods with decimals (e.g., +15.99491) by finding first/last
    # period that's followed by exactly one letter (the flanking AA).
    first_dot = pep.find(".")
    last_dot = pep.rfind(".")

    if first_dot >= 0 and last_dot > first_dot:
        pre_aa = pep[first_dot - 1:first_dot] if first_dot > 0 else "?"
        body = pep[first_dot + 1:last_dot]
        # last_aa = last A-Z letter in body (mods are non-letters between residues)
        body_letters = [c for c in body if c.isalpha()]
        last_aa = body_letters[-1] if body_letters else "?"
        length = len(body_letters)
    else:
        pre_aa = "?"
        last_aa = "?"
        length = 0

    # Charge from one-hot columns (PIN convention).
    charge = 2  # default fallback
    for z, key in [(2, "charge2"), (3, "charge3"), (4, "charge4")]:
        if row.get(key, "0") == "1":
            charge = z
            break

    raw_score = int(row.get("RawScore", "0"))
    if raw_score <= -10:
        score_bucket = "very_weak"
    elif raw_score <= 0:
        score_bucket = "weak"
    elif raw_score <= 49:
        score_bucket = "medium"
    elif raw_score <= 199:
        score_bucket = "strong"
    else:
        score_bucket = "very_strong"

    return {
        "length": length,
        "charge": charge,
        "n_oxidation": len(_OXIDATION_RE.findall(pep)),
        "n_carbamidomethyl": len(_CARBAMIDOMETHYL_RE.findall(pep)),
        "iso_off": int(row.get("isotope_error", "0")),
        "last_aa": last_aa,
        "pre_aa": pre_aa,
        "is_decoy": row.get("Label", "1") == "-1",
        "score_bucket": score_bucket,
    }


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


# ── Tests for peptide_features ──────────────────────────────────────────

def _test_peptide_features_basic_tryptic():
    row = {
        "Label": "1", "isotope_error": "0", "peplen": "11",
        "RawScore": "65", "charge2": "1", "charge3": "0", "charge4": "0",
        "Peptide": "K.SLKKISVIK.D",
    }
    f = peptide_features(row)
    assert f["length"] == 9, f"expected length 9, got {f['length']}"
    assert f["charge"] == 2
    assert f["n_oxidation"] == 0
    assert f["n_carbamidomethyl"] == 0
    assert f["iso_off"] == 0
    assert f["last_aa"] == "K"
    assert f["pre_aa"] == "K"
    assert f["is_decoy"] is False
    assert f["score_bucket"] == "strong"

def _test_peptide_features_with_mods():
    row = {
        "Label": "-1", "isotope_error": "1", "peplen": "12",
        "RawScore": "-15", "charge2": "0", "charge3": "1", "charge4": "0",
        "Peptide": "R.AC+57.0214DEM+15.9949FKPSQ.G",
    }
    f = peptide_features(row)
    assert f["length"] == 10
    assert f["charge"] == 3
    assert f["n_oxidation"] == 1
    assert f["n_carbamidomethyl"] == 1
    assert f["iso_off"] == 1
    assert f["last_aa"] == "Q"
    assert f["pre_aa"] == "R"
    assert f["is_decoy"] is True
    assert f["score_bucket"] == "very_weak"

def _test_peptide_features_protein_n_term():
    row = {
        "Label": "1", "isotope_error": "0", "peplen": "8",
        "RawScore": "5", "charge2": "1", "charge3": "0", "charge4": "0",
        "Peptide": "_.MSEAQR.K",
    }
    f = peptide_features(row)
    assert f["pre_aa"] == "_", f"expected pre_aa='_' for protein N-term, got {f['pre_aa']}"
    assert f["length"] == 6
    assert f["score_bucket"] == "medium"

def _test_peptide_features_score_buckets():
    cases = [
        (-30, "very_weak"), (-10, "very_weak"), (-5, "weak"),
        (0, "weak"), (49, "medium"),
        (50, "strong"), (199, "strong"),
        (200, "very_strong"), (500, "very_strong"),
    ]
    base_row = {
        "Label": "1", "isotope_error": "0", "peplen": "8",
        "charge2": "1", "charge3": "0", "charge4": "0", "Peptide": "K.AAA.B",
    }
    for raw, expected in cases:
        row = {**base_row, "RawScore": str(raw)}
        f = peptide_features(row)
        assert f["score_bucket"] == expected, f"raw={raw}: expected {expected}, got {f['score_bucket']}"


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
        ("peptide_features basic tryptic", _test_peptide_features_basic_tryptic),
        ("peptide_features with mods", _test_peptide_features_with_mods),
        ("peptide_features protein N-term", _test_peptide_features_protein_n_term),
        ("peptide_features score buckets", _test_peptide_features_score_buckets),
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
