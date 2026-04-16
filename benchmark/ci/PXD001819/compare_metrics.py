#!/usr/bin/env python3
"""Compare key=value metrics against a baseline TSV (metric, min, max, optional)."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


_TRUTHY = {"1", "y", "yes", "true", "t"}


def load_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _is_missing(val: str) -> bool:
    return val == "" or val.upper() == "NA"


def _check_row(metrics: dict[str, str], row: dict[str, str]) -> tuple[str | None, str | None]:
    """Return (failure, warning). At most one is non-None per row."""
    name = (row.get("metric") or "").strip()
    if not name:
        return None, None
    lo_raw = (row.get("min") or "").strip()
    hi_raw = (row.get("max") or "").strip()
    if lo_raw == "" and hi_raw == "":
        return None, None
    optional = (row.get("optional") or "").strip().lower() in _TRUTHY

    val_raw = metrics.get(name, "").strip()
    if name not in metrics or _is_missing(val_raw):
        msg = f"{name!r} missing or NA"
        return (None, f"{msg} (optional; skipped)") if optional else (f"{msg} (required by baseline)", None)

    try:
        val = float(val_raw)
    except ValueError:
        return f"{name}={val_raw!r} is not numeric", None
    lo = float(lo_raw) if lo_raw else float("-inf")
    hi = float(hi_raw) if hi_raw else float("inf")
    if val < lo or val > hi:
        return f"{name}={val} outside [{lo}, {hi}]", None
    return None, None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("metrics", type=Path, help="key=value metrics file from CI")
    ap.add_argument("baseline", type=Path, help="TSV with columns metric, min, max, optional")
    args = ap.parse_args()

    metrics = load_kv(args.metrics)
    failures: list[str] = []
    warnings: list[str] = []

    with args.baseline.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            failure, warning = _check_row(metrics, row)
            if failure:
                failures.append(failure)
            if warning:
                warnings.append(warning)

    for line in warnings:
        print(f"warning: {line}", file=sys.stderr)

    if failures:
        print("Benchmark baseline comparison failed:", file=sys.stderr)
        for line in failures:
            print(f"  - {line}", file=sys.stderr)
        return 1
    print("All checked metrics within baseline ranges.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
