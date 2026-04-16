#!/usr/bin/env python3
"""Compare key=value metrics against a baseline TSV (metric, min, max)."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def load_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _check_row(metrics: dict[str, str], row: dict[str, str]) -> str | None:
    name = (row.get("metric") or "").strip()
    if not name:
        return None
    lo_raw = (row.get("min") or "").strip()
    hi_raw = (row.get("max") or "").strip()
    if lo_raw == "" and hi_raw == "":
        return None
    if name not in metrics:
        return f"missing metric {name!r} (expected for baseline row)"
    val_raw = metrics[name].strip()
    if val_raw == "" or val_raw.upper() == "NA":
        return f"{name} is empty (cannot compare to baseline)"
    try:
        val = float(val_raw)
    except ValueError:
        return f"{name}={val_raw!r} is not numeric"
    lo = float(lo_raw) if lo_raw else float("-inf")
    hi = float(hi_raw) if hi_raw else float("inf")
    if val < lo or val > hi:
        return f"{name}={val} outside [{lo}, {hi}]"
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("metrics", type=Path, help="key=value metrics file from CI")
    ap.add_argument("baseline", type=Path, help="TSV with columns metric, min, max")
    args = ap.parse_args()

    metrics = load_kv(args.metrics)
    failures: list[str] = []

    with args.baseline.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            msg = _check_row(metrics, row)
            if msg:
                failures.append(msg)

    if failures:
        print("Benchmark baseline comparison failed:", file=sys.stderr)
        for line in failures:
            print(f"  - {line}", file=sys.stderr)
        return 1
    print("All checked metrics within baseline ranges.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
