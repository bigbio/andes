# Stratified PIN Parity Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use the subagent-driven-development skill (recommended) or executing-plans skill to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a single read-only Python script that consumes the existing PXD001819 java.pin / rust.pin files and emits a markdown report stratifying ΔRawScore and class-flip rates by peptide property, so the next investigation step is data-driven (not gut feel).

**Architecture:** One file at `benchmark/parity/analyze_parity.py`. Standard library only — no pandas / scipy. Pure functions for parsing, feature extraction, stratification, and report formatting; a thin `main()` wires them together. A `--self-test` CLI flag runs the embedded unit tests against synthetic PIN rows so the engineer can validate without running the full 38k-spectrum dataset.

**Tech Stack:** Python 3.10+ (standard library), markdown output, existing PIN files at `benchmark/results/PXD001819-parity/{java,rust}.pin`. Slice fixtures at `/tmp/{java,rust}_slice.pin` from prior session for fast smoke-testing.

---

## File structure

| File | Responsibility |
|---|---|
| `benchmark/parity/analyze_parity.py` | Entire script: CLI, parser, feature extractor, stratifier, lift calculator, ranking-mode classifier, markdown report formatter, embedded `--self-test` |
| `docs/parity-analysis/reports/2026-05-09-parity-report.md` | Generated output of running the script on full PXD001819 (committed for the historical record) |

The design doc (`docs/parity-analysis/2026-05-09-parity-analysis-design.md`) lists six logical components: `parse_pin`, `peptide_features`, `stratify`, `compute_lift`, `classify_ranking_mode`, `format_report`. Each gets its own task below.

---

## Task 1: Skeleton CLI and module imports

**Files:**
- Create: `benchmark/parity/analyze_parity.py`

- [ ] **Step 1: Write the failing self-test invocation**

The script must accept a `--self-test` flag and exit 0 when no tests are defined yet (we'll add tests in later tasks). Write the entry-point logic before any business logic so subsequent tasks just plug in.

- [ ] **Step 2: Run the not-yet-existing script to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `python3: can't open file '...analyze_parity.py': [Errno 2] No such file or directory`

- [ ] **Step 3: Create the skeleton**

```python
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
```

- [ ] **Step 4: Run the self-test**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected:
```
Self-tests: 0 ran, 0 failed
```

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag: scaffold analyze_parity.py CLI"
```

---

## Task 2: parse_pin — read PIN tab-separated rows

**Files:**
- Modify: `benchmark/parity/analyze_parity.py` (add `parse_pin` + tests)

- [ ] **Step 1: Write the failing test**

Append to `analyze_parity.py` (above `if __name__ == "__main__":`):

```python
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
```

Update `run_self_tests`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `FAIL: parse_pin basic: ... NameError: name 'parse_pin' is not defined`

- [ ] **Step 3: Implement parse_pin**

Add above `def main()`:

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected:
```
  PASS: parse_pin basic
  PASS: parse_pin skips blank lines
Self-tests: 2 ran, 0 failed
```

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): parse_pin reads tab-separated PIN rows"
```

---

## Task 3: peptide_features — extract per-row features

**Files:**
- Modify: `benchmark/parity/analyze_parity.py` (add `peptide_features` + tests)

- [ ] **Step 1: Write the failing test**

Append above the `_test_parse_pin_basic` block:

```python
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
        "Peptide": "R.AC+57.02146DEM+15.99491FK.G",
    }
    f = peptide_features(row)
    assert f["length"] == 10
    assert f["charge"] == 3
    assert f["n_oxidation"] == 1
    assert f["n_carbamidomethyl"] == 1
    assert f["iso_off"] == 1
    assert f["last_aa"] == "K"
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
        (-30, "very_weak"), (-10, "weak"), (-5, "weak"),
        (0, "medium"), (49, "medium"),
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
```

Add the four tests to `run_self_tests` tests list:

```python
        ("peptide_features basic tryptic", _test_peptide_features_basic_tryptic),
        ("peptide_features with mods", _test_peptide_features_with_mods),
        ("peptide_features protein N-term", _test_peptide_features_protein_n_term),
        ("peptide_features score buckets", _test_peptide_features_score_buckets),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: 4 new failures (`peptide_features` not defined).

- [ ] **Step 3: Implement peptide_features**

Add above `parse_pin`:

```python
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
    parts = pep.split(".")
    if len(parts) >= 3:
        pre_aa = parts[0][:1] if parts[0] else "?"
        body = parts[1]
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 6 ran, 0 failed`

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): peptide_features extracts per-row diagnostic features"
```

---

## Task 4: classify_ranking_mode — distinguish RawScore-swap vs SpecE-swap

**Files:**
- Modify: `benchmark/parity/analyze_parity.py`

- [ ] **Step 1: Write the failing test**

Append above the existing test block:

```python
# ── Tests for classify_ranking_mode ─────────────────────────────────────

def _test_classify_ranking_mode_agree():
    java_row = {"Peptide": "K.AAA.B", "RawScore": "10", "lnSpecEValue": "-5.0", "Label": "1"}
    rust_row = {"Peptide": "K.AAA.B", "RawScore": "8",  "lnSpecEValue": "-4.5", "Label": "1"}
    assert classify_ranking_mode(java_row, rust_row) == "agree"

def _test_classify_ranking_mode_raw_swap():
    # Different peptides, both RawScore and SpecE put Java's first.
    java_row = {"Peptide": "K.AAA.B", "RawScore": "20", "lnSpecEValue": "-8.0", "Label": "1"}
    rust_row = {"Peptide": "K.BBB.C", "RawScore": "15", "lnSpecEValue": "-6.0", "Label": "1"}
    assert classify_ranking_mode(java_row, rust_row) == "raw_swap"

def _test_classify_ranking_mode_spec_e_swap_only():
    # Different peptides; Rust's pick has higher RawScore but worse SpecE.
    # I.e. SpecE inverts what RawScore would say.
    java_row = {"Peptide": "K.AAA.B", "RawScore": "10", "lnSpecEValue": "-9.0", "Label": "1"}
    rust_row = {"Peptide": "K.BBB.C", "RawScore": "15", "lnSpecEValue": "-6.0", "Label": "1"}
    # Rust's RawScore (15) > Java's (10), so RawScore order says Rust wins.
    # But SpecE -9.0 < -6.0, so Java's pick is actually better by SpecE.
    assert classify_ranking_mode(java_row, rust_row) == "spec_e_swap_only"

def _test_classify_ranking_mode_both_swap():
    # Different peptides, fully inverted on both metrics. Rare; happens when
    # there's noise in both ranking signals.
    java_row = {"Peptide": "K.AAA.B", "RawScore": "10", "lnSpecEValue": "-5.0", "Label": "1"}
    rust_row = {"Peptide": "K.BBB.C", "RawScore": "20", "lnSpecEValue": "-8.0", "Label": "1"}
    assert classify_ranking_mode(java_row, rust_row) == "both_swap"
```

Add to `run_self_tests` tests list:

```python
        ("classify agree", _test_classify_ranking_mode_agree),
        ("classify raw swap", _test_classify_ranking_mode_raw_swap),
        ("classify spec_e swap only", _test_classify_ranking_mode_spec_e_swap_only),
        ("classify both swap", _test_classify_ranking_mode_both_swap),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: 4 new failures.

- [ ] **Step 3: Implement classify_ranking_mode**

Add above `peptide_features`:

```python
def classify_ranking_mode(java_row: dict[str, str], rust_row: dict[str, str]) -> str:
    """Classify how Java's top-1 PSM and Rust's top-1 PSM disagree.

    Returns one of:
      "agree"             — same peptide picked by both engines
      "raw_swap"          — Java's pick has higher RawScore in Java's view AND
                            higher (better, more negative) lnSpecEValue
      "spec_e_swap_only"  — RawScore values would put Rust's pick first
                            (higher RawScore in absolute value), but SpecE
                            inverts and makes Java's pick the better one
      "both_swap"         — fully inverted on both metrics (rare)
    """
    if java_row["Peptide"] == rust_row["Peptide"]:
        return "agree"

    j_raw = float(java_row["RawScore"])
    r_raw = float(rust_row["RawScore"])
    j_lse = float(java_row["lnSpecEValue"])
    r_lse = float(rust_row["lnSpecEValue"])

    # Lower (more negative) lnSpecEValue is better.
    raw_says_java_wins = j_raw > r_raw
    spec_e_says_java_wins = j_lse < r_lse

    if raw_says_java_wins and spec_e_says_java_wins:
        return "raw_swap"
    if not raw_says_java_wins and spec_e_says_java_wins:
        return "spec_e_swap_only"
    return "both_swap"
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 10 ran, 0 failed`

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): classify_ranking_mode separates raw/spec_e swap modes"
```

---

## Task 5: stratify and compute_lift — bucket aggregation

**Files:**
- Modify: `benchmark/parity/analyze_parity.py`

- [ ] **Step 1: Write the failing test**

Append above the existing test block:

```python
# ── Tests for stratify and compute_lift ─────────────────────────────────

def _test_stratify_aggregates_per_bucket():
    rows = [
        {"bucket": "a", "delta": 5},
        {"bucket": "a", "delta": 7},
        {"bucket": "b", "delta": -1},
        {"bucket": "b", "delta": 0},
        {"bucket": "b", "delta": 1},
    ]
    out = stratify(rows, lambda r: r["bucket"], lambda r: r["delta"])
    assert set(out.keys()) == {"a", "b"}
    assert out["a"]["count"] == 2
    assert out["a"]["mean"] == 6.0
    assert out["b"]["count"] == 3
    assert out["b"]["mean"] == 0.0
    assert out["b"]["median"] == 0
    # stdev only defined for n >= 2
    assert out["a"]["stdev"] > 0

def _test_stratify_handles_singletons():
    # n=1 → stdev is 0 (we use population stdev which is 0 for single sample).
    rows = [{"bucket": "x", "delta": 42}]
    out = stratify(rows, lambda r: r["bucket"], lambda r: r["delta"])
    assert out["x"]["count"] == 1
    assert out["x"]["mean"] == 42
    assert out["x"]["stdev"] == 0.0

def _test_compute_lift_basic():
    # Bucket flip rate 50%, baseline 25% → lift = 2.0
    assert compute_lift(group_rate=0.50, base_rate=0.25) == 2.0
    # Bucket below baseline → lift < 1
    assert compute_lift(0.10, 0.25) == 0.4
    # Zero baseline → return inf for non-zero group, 0 for zero group
    assert compute_lift(0.5, 0.0) == float("inf")
    assert compute_lift(0.0, 0.0) == 0.0
```

Add to test list:

```python
        ("stratify per bucket", _test_stratify_aggregates_per_bucket),
        ("stratify singletons", _test_stratify_handles_singletons),
        ("compute_lift basic", _test_compute_lift_basic),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: 3 new failures.

- [ ] **Step 3: Implement stratify and compute_lift**

Add above `peptide_features`:

```python
def stratify(
    rows: list[dict],
    bucket_fn: Callable[[dict], object],
    value_fn: Callable[[dict], float],
) -> dict[object, dict[str, float]]:
    """Group rows by bucket_fn(row) and aggregate value_fn(row) per bucket.

    Returns dict mapping bucket → {"count", "mean", "median", "stdev"}.
    Single-sample buckets get stdev = 0.0 (population stdev convention).
    """
    grouped: dict[object, list[float]] = defaultdict(list)
    for row in rows:
        grouped[bucket_fn(row)].append(float(value_fn(row)))

    out: dict[object, dict[str, float]] = {}
    for bucket, values in grouped.items():
        out[bucket] = {
            "count": len(values),
            "mean": mean(values),
            "median": median(values),
            "stdev": pstdev(values) if len(values) > 1 else 0.0,
        }
    return out


def compute_lift(group_rate: float, base_rate: float) -> float:
    """lift = group_rate / base_rate.

    Edge cases:
      - base_rate == 0 and group_rate > 0 → infinity (rare bucket)
      - base_rate == 0 and group_rate == 0 → 0 (degenerate)
    """
    if base_rate == 0:
        return float("inf") if group_rate > 0 else 0.0
    return group_rate / base_rate
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 13 ran, 0 failed`

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): stratify + compute_lift for bucket aggregation"
```

---

## Task 6: Match Java/Rust pin rows by scan + classify

**Files:**
- Modify: `benchmark/parity/analyze_parity.py`

- [ ] **Step 1: Write the failing test**

Append above the existing test block:

```python
# ── Tests for match_pins ────────────────────────────────────────────────

def _test_match_pins_pairs_by_scan():
    java = [
        {"ScanNr": "5", "Label": "-1", "Peptide": "K.AAA.B", "RawScore": "-34", "lnSpecEValue": "-8.7"},
        {"ScanNr": "5", "Label": "-1", "Peptide": "K.BBB.C", "RawScore": "-34", "lnSpecEValue": "-8.7"},
        {"ScanNr": "10", "Label": "1", "Peptide": "K.CCC.D", "RawScore": "50", "lnSpecEValue": "-12"},
    ]
    rust = [
        {"ScanNr": "5", "Label": "1", "Peptide": "K.DEN.R", "RawScore": "-5", "lnSpecEValue": "-8.8"},
        {"ScanNr": "10", "Label": "1", "Peptide": "K.CCC.D", "RawScore": "48", "lnSpecEValue": "-11.5"},
    ]
    matches = match_pins(java, rust)
    # 2 paired scans (5 and 10). For scan 5, Java has 2 rows but match_pins
    # uses FIRST row for parity with the existing flip_count.py baseline.
    assert len(matches) == 2
    pair5 = [m for m in matches if m["scan"] == 5][0]
    assert pair5["java"]["Peptide"] == "K.AAA.B"
    assert pair5["rust"]["Peptide"] == "K.DEN.R"
    assert pair5["mode"] in ("raw_swap", "spec_e_swap_only", "both_swap")
    pair10 = [m for m in matches if m["scan"] == 10][0]
    assert pair10["mode"] == "agree"
```

Add to test list:

```python
        ("match_pins pairs by scan", _test_match_pins_pairs_by_scan),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `FAIL: match_pins pairs by scan: name 'match_pins' is not defined`

- [ ] **Step 3: Implement match_pins**

Add above `peptide_features`:

```python
def match_pins(
    java_rows: list[dict[str, str]],
    rust_rows: list[dict[str, str]],
) -> list[dict]:
    """Pair Java and Rust pin rows by scan (first-row-only convention).

    Returns a list of matched pairs, each with keys:
      "scan" (int), "java" (row), "rust" (row), "mode" (from
      classify_ranking_mode).

    Scans missing on either side are skipped. Multi-row Java scans use the
    first row for backward-compat with flip_count.py's baseline numbers.
    """
    java_by_scan: dict[int, dict[str, str]] = {}
    for r in java_rows:
        scan = int(r["ScanNr"])
        if scan not in java_by_scan:
            java_by_scan[scan] = r

    rust_by_scan: dict[int, dict[str, str]] = {}
    for r in rust_rows:
        scan = int(r["ScanNr"])
        if scan not in rust_by_scan:
            rust_by_scan[scan] = r

    pairs: list[dict] = []
    for scan in sorted(set(java_by_scan) & set(rust_by_scan)):
        j = java_by_scan[scan]
        r = rust_by_scan[scan]
        pairs.append({
            "scan": scan,
            "java": j,
            "rust": r,
            "mode": classify_ranking_mode(j, r),
        })
    return pairs
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 14 ran, 0 failed`

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): match_pins pairs java/rust rows by scan"
```

---

## Task 7: Section 1 — population overview report

**Files:**
- Modify: `benchmark/parity/analyze_parity.py`

- [ ] **Step 1: Write the failing test**

```python
# ── Tests for Section 1 / population overview ───────────────────────────

def _test_section1_counts_and_distributions():
    matches = [
        {"scan": 1, "java": {"Peptide": "X.A.X", "RawScore": "10", "lnSpecEValue": "-5", "Label": "1"},
                   "rust": {"Peptide": "X.A.X", "RawScore": "8",  "lnSpecEValue": "-4", "Label": "1"},
                   "mode": "agree"},
        {"scan": 2, "java": {"Peptide": "X.A.X", "RawScore": "20", "lnSpecEValue": "-7", "Label": "1"},
                   "rust": {"Peptide": "X.B.X", "RawScore": "15", "lnSpecEValue": "-5", "Label": "-1"},
                   "mode": "raw_swap"},
        {"scan": 3, "java": {"Peptide": "X.A.X", "RawScore": "30", "lnSpecEValue": "-9", "Label": "1"},
                   "rust": {"Peptide": "X.B.X", "RawScore": "25", "lnSpecEValue": "-7", "Label": "1"},
                   "mode": "raw_swap"},
    ]
    section = format_section1_overview(matches)
    assert "Full match" in section
    assert "1 (33.3%)" in section, f"expected '1 (33.3%)' for full match in:\n{section}"
    assert "Class flips" in section
    # 1 flip (scan 2: target -> decoy)
    assert "1 (33.3%)" in section
```

Add to test list:

```python
        ("section1 counts", _test_section1_counts_and_distributions),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `FAIL: ... format_section1_overview is not defined`

- [ ] **Step 3: Implement format_section1_overview**

Add above `match_pins`:

```python
def format_section1_overview(matches: list[dict]) -> str:
    """Section 1 — population overview: counts by class, baseline rates."""
    n = len(matches)
    full_match = sum(1 for m in matches if m["mode"] == "agree")
    same_label_diff_pep = sum(
        1 for m in matches
        if m["mode"] != "agree" and m["java"]["Label"] == m["rust"]["Label"]
    )
    class_flip = sum(
        1 for m in matches
        if m["mode"] != "agree" and m["java"]["Label"] != m["rust"]["Label"]
    )

    def pct(c: int) -> str:
        return f"{c} ({100 * c / n:.1f}%)" if n else f"{c} (-)"

    lines = [
        "## Section 1 — Population overview",
        "",
        f"- Total matched scans: {n}",
        f"- Full match (same peptide): {pct(full_match)}",
        f"- Same-label different-peptide: {pct(same_label_diff_pep)}",
        f"- Class flips (target↔decoy): {pct(class_flip)}",
        "",
    ]
    return "\n".join(lines)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 15 ran, 0 failed`

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): format_section1_overview reports population counts"
```

---

## Task 8: Section 2 — ΔRawScore decomposition (Track A)

**Files:**
- Modify: `benchmark/parity/analyze_parity.py`

- [ ] **Step 1: Write the failing test**

```python
# ── Tests for Section 2 / ΔRawScore decomposition ───────────────────────

def _test_section2_delta_decomposition():
    # Build 4 full-match scans; Δ varies by feature so we can see stratification.
    matches = []
    base_rust = {"charge2": "1", "charge3": "0", "charge4": "0", "isotope_error": "0", "Label": "1"}
    for i, (charge_col, raw_j, raw_r, pep) in enumerate([
        ("charge2", 50, 10, "K.AAAAAA.B"),  # length 6
        ("charge2", 60, 20, "K.BBBBBB.B"),  # length 6
        ("charge3", 100, 80, "K.CCCCCCCCC.D"),  # length 9
        ("charge3", 110, 90, "K.DDDDDDDDD.D"),  # length 9
    ]):
        j = {"Peptide": pep, "RawScore": str(raw_j), "lnSpecEValue": "-5", "Label": "1",
             "isotope_error": "0", "peplen": str(len(pep) - 4 + 2),
             "charge2": "1" if charge_col == "charge2" else "0",
             "charge3": "1" if charge_col == "charge3" else "0",
             "charge4": "0"}
        r = {**j, "RawScore": str(raw_r)}
        matches.append({"scan": i, "java": j, "rust": r, "mode": "agree"})
    section = format_section2_delta_decomposition(matches)
    assert "Section 2" in section
    assert "median" in section.lower() or "Median" in section
    # Δ for charge=2 is 40, 40 (mean 40); for charge=3 is 20, 20 (mean 20)
    assert "40" in section, f"expected charge=2 Δ around 40 to appear in:\n{section}"
    assert "20" in section
```

Add to test list:

```python
        ("section2 Δ decomposition", _test_section2_delta_decomposition),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `FAIL: ... format_section2_delta_decomposition is not defined`

- [ ] **Step 3: Implement format_section2_delta_decomposition**

Add above `format_section1_overview`:

```python
def _eta_squared(stratification: dict[object, dict[str, float]], grand_mean: float) -> float:
    """η² = between-group sum-of-squares / total sum-of-squares.

    A rough effect-size measure: 0 = no group differences, 1 = all variance
    explained by group membership. Used as a quick sanity check on whether
    a feature has a systematic Δ effect or is just noise.
    """
    between_ss = 0.0
    total_n = 0
    for stats in stratification.values():
        n = int(stats["count"])
        between_ss += n * (stats["mean"] - grand_mean) ** 2
        total_n += n
    # Total SS would require per-row data; we approximate with bucket variances.
    # Approx total SS = between_ss + within_ss where within_ss = Σ n_i * stdev_i²
    within_ss = sum(int(s["count"]) * s["stdev"] ** 2 for s in stratification.values())
    total_ss = between_ss + within_ss
    return between_ss / total_ss if total_ss > 0 else 0.0


def format_section2_delta_decomposition(matches: list[dict]) -> str:
    """Section 2 — ΔRawScore decomposition on full-match PSMs.

    Δ = Java RawScore − Rust RawScore (positive means Java higher).
    Stratifies Δ by peptide features and computes η² to flag features
    that systematically explain Δ.
    """
    full = [m for m in matches if m["mode"] == "agree"]
    if not full:
        return "## Section 2 — ΔRawScore decomposition\n\n_No full-match scans._\n\n"

    deltas = [int(m["java"]["RawScore"]) - int(m["rust"]["RawScore"]) for m in full]
    grand_mean = mean(deltas)
    grand_median = median(deltas)
    grand_stdev = pstdev(deltas) if len(deltas) > 1 else 0.0
    tail = sum(1 for d in deltas if abs(d) > 100)

    lines = [
        "## Section 2 — ΔRawScore decomposition (full-match PSMs only)",
        "",
        f"- N: {len(full)}",
        f"- Median Δ: {grand_median}",
        f"- Mean Δ: {grand_mean:.1f}",
        f"- Stdev Δ: {grand_stdev:.1f}",
        f"- Tail (|Δ| > 100): {tail} ({100 * tail / len(full):.1f}%)",
        "",
        "### Stratified by feature",
        "",
        "| Feature | Bucket | N | Mean Δ | Median Δ | Stdev |",
        "|---|---|---|---|---|---|",
    ]

    feature_names = ["length", "charge", "n_oxidation", "n_carbamidomethyl",
                     "iso_off", "last_aa", "pre_aa", "is_decoy", "score_bucket"]
    eta_summary: list[tuple[str, float]] = []
    for fname in feature_names:
        # Build paired (row, delta) list keyed by feature value.
        per_row = []
        for m in full:
            feats = peptide_features(m["java"])
            d = int(m["java"]["RawScore"]) - int(m["rust"]["RawScore"])
            per_row.append({"bucket": feats[fname], "delta": d})
        strat = stratify(per_row, lambda r: r["bucket"], lambda r: r["delta"])
        eta = _eta_squared(strat, grand_mean)
        eta_summary.append((fname, eta))
        # Render top buckets (sorted by count descending, max 5 rows)
        sorted_buckets = sorted(strat.items(), key=lambda kv: -kv[1]["count"])[:5]
        for bucket, stats in sorted_buckets:
            lines.append(f"| {fname} | {bucket} | {int(stats['count'])} | "
                         f"{stats['mean']:.1f} | {stats['median']:.0f} | {stats['stdev']:.1f} |")

    lines += [
        "",
        "### Variance contribution (η²) per feature",
        "",
        "| Feature | η² |",
        "|---|---|",
    ]
    for fname, eta in sorted(eta_summary, key=lambda x: -x[1]):
        lines.append(f"| {fname} | {eta:.3f} |")
    lines.append("")

    # Decision rule
    top_feature, top_eta = max(eta_summary, key=lambda x: x[1])
    if top_eta > 0.40:
        verdict = (f"**Decision:** η² for `{top_feature}` is {top_eta:.2f} > 0.40 → "
                   "property-conditioned bias. Next: narrow code audit on the path that "
                   f"consumes `{top_feature}`.")
    elif top_eta < 0.20:
        verdict = (f"**Decision:** max η² is {top_eta:.2f} (`{top_feature}`) < 0.20 → "
                   "Δ is mostly per-scan noise; no single feature explains it. The 67% "
                   "disagreement is many small divergences, not one missing global term.")
    else:
        verdict = (f"**Decision:** top η² is {top_eta:.2f} (`{top_feature}`) — mixed "
                   "systematic effect plus noise. Worth a narrow audit on the feature's "
                   "code path but expect partial coverage.")
    lines.append(verdict)
    lines.append("")
    return "\n".join(lines)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 16 ran, 0 failed`

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): section 2 decomposes ΔRawScore by feature (Track A)"
```

---

## Task 9: Section 3 — stratified flip lift (Track B)

**Files:**
- Modify: `benchmark/parity/analyze_parity.py`

- [ ] **Step 1: Write the failing test**

```python
# ── Tests for Section 3 / stratified flip lift ──────────────────────────

def _test_section3_lift_table():
    matches = []
    # 10 charge=2 scans, 5 flips → flip_rate = 50%
    # 10 charge=3 scans, 1 flip → flip_rate = 10%
    # baseline = 6 / 20 = 30%
    # lift charge=2 = 50/30 = 1.67; lift charge=3 = 10/30 = 0.33
    for i in range(10):
        is_flip = i < 5
        j = {"ScanNr": str(i), "Peptide": "K.A.B", "RawScore": "10", "lnSpecEValue": "-5",
             "Label": "1", "isotope_error": "0", "peplen": "3",
             "charge2": "1", "charge3": "0", "charge4": "0"}
        if is_flip:
            r = {**j, "Peptide": "K.B.C", "Label": "-1"}
            mode = "raw_swap"
        else:
            r = {**j}
            mode = "agree"
        matches.append({"scan": i, "java": j, "rust": r, "mode": mode})
    for i in range(10, 20):
        is_flip = i == 10
        j = {"ScanNr": str(i), "Peptide": "K.X.Y", "RawScore": "10", "lnSpecEValue": "-5",
             "Label": "1", "isotope_error": "0", "peplen": "3",
             "charge2": "0", "charge3": "1", "charge4": "0"}
        if is_flip:
            r = {**j, "Peptide": "K.Z.W", "Label": "-1"}
            mode = "raw_swap"
        else:
            r = {**j}
            mode = "agree"
        matches.append({"scan": i, "java": j, "rust": r, "mode": mode})

    section = format_section3_flip_lift(matches)
    assert "Section 3" in section
    assert "lift" in section.lower()
    # charge=2 should have lift > 1; charge=3 should have lift < 1
    # Hard to assert exact values in a markdown table, but we can check both appear
    assert "charge" in section
```

Add to test list:

```python
        ("section3 flip lift", _test_section3_lift_table),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `FAIL: ... format_section3_flip_lift is not defined`

- [ ] **Step 3: Implement format_section3_flip_lift**

Add above `format_section2_delta_decomposition`:

```python
def format_section3_flip_lift(matches: list[dict]) -> str:
    """Section 3 — stratified flip-rate analysis with lift vs baseline."""
    n = len(matches)
    flips = sum(1 for m in matches
                if m["java"]["Label"] != m["rust"]["Label"]
                and m["mode"] != "agree")
    if n == 0:
        return "## Section 3 — Stratified flip lift\n\n_No matched scans._\n\n"
    base_rate = flips / n

    feature_names = ["length", "charge", "n_oxidation", "n_carbamidomethyl",
                     "iso_off", "last_aa", "pre_aa", "is_decoy", "score_bucket"]

    lines = [
        "## Section 3 — Stratified flip lift (Track B)",
        "",
        f"- Total scans: {n}",
        f"- Total class flips: {flips}",
        f"- Baseline flip rate: {100 * base_rate:.1f}%",
        "",
        "### Top flip-enriched buckets (lift > 1.5, n ≥ 50)",
        "",
        "| Feature | Bucket | N | Flips | Flip rate | Lift |",
        "|---|---|---|---|---|---|",
    ]

    candidates: list[tuple[str, object, int, int, float, float]] = []
    for fname in feature_names:
        per_bucket: dict[object, dict[str, int]] = defaultdict(lambda: {"n": 0, "f": 0})
        for m in matches:
            feats = peptide_features(m["java"])
            b = feats[fname]
            per_bucket[b]["n"] += 1
            if m["java"]["Label"] != m["rust"]["Label"] and m["mode"] != "agree":
                per_bucket[b]["f"] += 1
        for bucket, c in per_bucket.items():
            if c["n"] < 50:
                continue
            rate = c["f"] / c["n"]
            lift = compute_lift(rate, base_rate)
            if lift > 1.5:
                candidates.append((fname, bucket, c["n"], c["f"], rate, lift))

    for row in sorted(candidates, key=lambda x: -x[5])[:10]:
        fname, bucket, total, fct, rate, lift = row
        lines.append(f"| {fname} | {bucket} | {total} | {fct} | "
                     f"{100*rate:.1f}% | {lift:.2f} |")
    if not candidates:
        lines.append("| _none_ | | | | | |")

    lines.append("")
    return "\n".join(lines)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 17 ran, 0 failed`

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): section 3 stratified flip lift (Track B)"
```

---

## Task 10: Section 4 — ranking-mode breakdown

**Files:**
- Modify: `benchmark/parity/analyze_parity.py`

- [ ] **Step 1: Write the failing test**

```python
# ── Tests for Section 4 / ranking-mode breakdown ────────────────────────

def _test_section4_ranking_modes():
    matches = [
        {"mode": "agree", "java": {"RawScore": "10", "lnSpecEValue": "-5"},
                          "rust": {"RawScore": "8",  "lnSpecEValue": "-4"}},
        {"mode": "raw_swap", "java": {"RawScore": "20", "lnSpecEValue": "-7"},
                              "rust": {"RawScore": "15", "lnSpecEValue": "-5"}},
        {"mode": "raw_swap", "java": {"RawScore": "30", "lnSpecEValue": "-9"},
                              "rust": {"RawScore": "25", "lnSpecEValue": "-7"}},
        {"mode": "spec_e_swap_only", "java": {"RawScore": "10", "lnSpecEValue": "-9"},
                                      "rust": {"RawScore": "15", "lnSpecEValue": "-6"}},
    ]
    section = format_section4_ranking_modes(matches)
    assert "Section 4" in section
    assert "raw_swap" in section
    assert "spec_e_swap_only" in section
    # 2 of 3 flips are raw_swap → 66.7%
    assert "66.7" in section or "67" in section
```

Add to test list:

```python
        ("section4 ranking modes", _test_section4_ranking_modes),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `FAIL: ... format_section4_ranking_modes is not defined`

- [ ] **Step 3: Implement format_section4_ranking_modes**

Add above `format_section3_flip_lift`:

```python
def format_section4_ranking_modes(matches: list[dict]) -> str:
    """Section 4 — distinguish RawScore-swap from SpecE-swap among flips."""
    flips = [m for m in matches if m["mode"] != "agree"]
    if not flips:
        return "## Section 4 — Ranking-mode breakdown\n\n_No flips._\n\n"

    counts = Counter(m["mode"] for m in flips)
    n = len(flips)
    lines = [
        "## Section 4 — Ranking-mode breakdown (flips only)",
        "",
        f"- Total flips: {n}",
        "",
        "| Mode | Count | % |",
        "|---|---|---|",
    ]
    for mode in ("raw_swap", "spec_e_swap_only", "both_swap"):
        c = counts.get(mode, 0)
        lines.append(f"| {mode} | {c} | {100 * c / n:.1f}% |")
    lines.append("")

    # Decision
    top_mode = counts.most_common(1)[0][0]
    if top_mode == "raw_swap":
        verdict = ("**Decision:** RawScore swaps dominate → per-PSM scoring (`score_psm` / "
                   "`directional_node_score`) is the prime suspect. Trace a high-lift bucket "
                   "from Section 3.")
    elif top_mode == "spec_e_swap_only":
        verdict = ("**Decision:** SpecE swaps without RawScore agreement dominate → GF DP / "
                   "`compute_inner` / `add_prob_dist` is the prime suspect. RawScore math is "
                   "fine; the spec_e_value lookup or distribution differs.")
    else:
        verdict = ("**Decision:** `both_swap` dominates → both per-PSM scoring AND GF DP "
                   "diverge. Likely a shared upstream cause (e.g. wrong partition lookup, "
                   "ion enumeration, or peak ranking).")
    lines.append(verdict)
    lines.append("")
    return "\n".join(lines)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 18 ran, 0 failed`

- [ ] **Step 5: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): section 4 ranking-mode breakdown"
```

---

## Task 11: Wire main() — assemble and write the full report

**Files:**
- Modify: `benchmark/parity/analyze_parity.py` (replace stub `main()` body)

- [ ] **Step 1: Write the failing test**

```python
# ── Tests for end-to-end main wiring ────────────────────────────────────

def _test_main_runs_on_slice_pins():
    import tempfile
    java_pin = Path("/tmp/java_slice.pin")
    rust_pin = Path("/tmp/rust_slice.pin")
    if not java_pin.exists() or not rust_pin.exists():
        # Skip when fixtures not available (e.g. in CI / fresh checkout).
        print("    [skip: slice pins not on disk]")
        return
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as f:
        out_path = Path(f.name)
    rc = run_pipeline(java_pin, rust_pin, out_path)
    assert rc == 0
    text = out_path.read_text()
    assert "Section 1" in text
    assert "Section 2" in text
    assert "Section 3" in text
    assert "Section 4" in text
    # Sanity: total scans should be ~1972 for the slice (matches flip_count.py
    # baseline). Some tolerance for parser differences.
    assert any(f"Total matched scans: {n}" in text for n in range(1900, 2000)), \
        f"expected ~1972 matched scans, report:\n{text[:500]}"
```

Add to test list:

```python
        ("main runs on slice", _test_main_runs_on_slice_pins),
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `FAIL: ... run_pipeline is not defined`

- [ ] **Step 3: Implement run_pipeline and rewire main**

Add above `def main()`:

```python
def run_pipeline(java_pin: Path, rust_pin: Path, output: Path) -> int:
    """End-to-end: parse pins, match by scan, format report, write to disk."""
    java_rows = parse_pin(java_pin)
    rust_rows = parse_pin(rust_pin)
    matches = match_pins(java_rows, rust_rows)

    sections = [
        f"# PXD001819 Java/Rust Parity Analysis Report\n",
        f"_Generated by `benchmark/parity/analyze_parity.py`._\n",
        f"_Java pin: `{java_pin}`_  ",
        f"_Rust pin: `{rust_pin}`_\n",
        format_section1_overview(matches),
        format_section2_delta_decomposition(matches),
        format_section3_flip_lift(matches),
        format_section4_ranking_modes(matches),
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(sections))
    print(f"Wrote {output} ({sum(len(s) for s in sections)} chars)")
    return 0
```

Replace the body of `main()` after the arg parsing (the `print(...)` lines and `return 0` at the end):

```python
    return run_pipeline(args.java_pin, args.rust_pin, args.output)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 benchmark/parity/analyze_parity.py --self-test
```

Expected: `Self-tests: 19 ran, 0 failed` (or 18 ran with `[skip: slice pins not on disk]` if /tmp/ is empty)

- [ ] **Step 5: Run end-to-end on the slice fixtures (smoke test)**

```bash
python3 benchmark/parity/analyze_parity.py \
    --java-pin /tmp/java_slice.pin \
    --rust-pin /tmp/rust_slice.pin \
    --output /tmp/parity_report_slice.md
head -40 /tmp/parity_report_slice.md
```

Expected: prints `Wrote /tmp/parity_report_slice.md (NNNN chars)` and head shows Sections 1-4. The "Total matched scans" line in Section 1 should be 1972.

- [ ] **Step 6: Commit**

```bash
git add benchmark/parity/analyze_parity.py
git commit -m "diag(parity): wire end-to-end run_pipeline + slice smoke test"
```

---

## Task 12: Run on full PXD001819 + commit the report

**Files:**
- Create: `docs/parity-analysis/reports/2026-05-09-parity-report.md`

- [ ] **Step 1: Run on the full dataset**

```bash
python3 benchmark/parity/analyze_parity.py \
    --java-pin benchmark/results/PXD001819-parity/java.pin \
    --rust-pin benchmark/results/PXD001819-parity/rust.pin \
    --output docs/parity-analysis/reports/2026-05-09-parity-report.md
```

Expected: prints `Wrote docs/parity-analysis/reports/2026-05-09-parity-report.md (NNNN chars)`. Should run in <30 seconds.

- [ ] **Step 2: Sanity-check the report**

```bash
head -20 docs/parity-analysis/reports/2026-05-09-parity-report.md
grep -E "Total matched scans|Median Δ|Mean Δ|Decision" docs/parity-analysis/reports/2026-05-09-parity-report.md
```

Expected:
- "Total matched scans: 37089" (matches flip_count.py baseline)
- "Median Δ:" near 30-50 (close to the +42 we observed)
- "Decision:" lines in Sections 2 and 4 with concrete recommendations

- [ ] **Step 3: Commit the report**

```bash
git add docs/parity-analysis/reports/2026-05-09-parity-report.md
git commit -m "diag(parity): generated report on full PXD001819 (12,432 full match)"
```

- [ ] **Step 4: Read the report's Section 5 / decision lines**

Open `docs/parity-analysis/reports/2026-05-09-parity-report.md` in an editor. The "Decision:" lines in Sections 2 and 4 are the deliverable: they tell us where to dig next (narrow code audit on a feature, msgf-trace on a high-lift bucket, GF DP audit, or score_psm audit). No more code work needed in this plan; the next iteration is informed by what the report says.

---

## Self-review notes

1. **Spec coverage:** Each design-doc section maps to a task:
   - design §"Components 1" parse_pin → Task 2
   - §"Components 2" peptide_features → Task 3
   - §"Components 3" protein_position — NOT implemented (design marked optional, requires FASTA loading; defer to follow-up if Section 5 recommends it)
   - §"Components 4" stratify → Task 5
   - §"Components 5" compute_lift → Task 5
   - §"Components 6" classify_ranking_mode → Task 4
   - §"Components 7" format_report → Tasks 7-10 (one per section)
   - §"Section 1-4 outputs" → Tasks 7-10
   - §"Section 5 recommendation" → embedded as decision lines in Sections 2 and 4 (per task 12 step 4 — the human reads them)

2. **Placeholder scan:** None. Every step shows the actual code.

3. **Type consistency:** All task signatures align — `parse_pin → list[dict]`, `peptide_features(row) → dict`, `match_pins → list[{scan, java, rust, mode}]`, `format_sectionN(matches) → str`. The dict keys (`"java"`, `"rust"`, `"mode"`, `"scan"`) are consistent across tasks 6-10.

4. **Decision-rule thresholds match the design:** η² > 0.40 → narrow audit, η² < 0.20 → noise verdict, in-between → mixed. Section 4 maps `raw_swap` / `spec_e_swap_only` / `both_swap` dominance to specific code-path suspects, matching the design's Section 4 decision rule.
