#!/usr/bin/env python3
"""Compare whether baseline top-ranked peptides are recovered in an experimental top-K pin."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


RANK_RE = re.compile(r"_(\d+)$")


def peptide_core(value: str) -> str:
    parts = value.split(".")
    if len(parts) == 3:
        return parts[1]
    return value


def pin_rank(row: dict[str, str]) -> int:
    spec_id = row.get("SpecId", "")
    match = RANK_RE.search(spec_id)
    if match:
        return int(match.group(1))
    return 1


def load_pin_groups(path: Path, max_rank: int, label_filter: str) -> dict[str, set[str]]:
    groups: dict[str, set[str]] = defaultdict(set)
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ScanNr", "Peptide", "Label"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            if pin_rank(row) > max_rank:
                continue
            if label_filter != "any" and row.get("Label") != label_filter:
                continue
            scan = row["ScanNr"]
            groups[scan].add(peptide_core(row["Peptide"]))
    return groups


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Measure whether the baseline top-ranked peptide for each spectrum "
            "is recovered in an experimental top-K .pin file."
        )
    )
    parser.add_argument("baseline_pin", type=Path)
    parser.add_argument("experimental_pin", type=Path)
    parser.add_argument("--baseline-rank", type=int, default=1, help="Maximum baseline rank to consider")
    parser.add_argument("--experimental-rank", type=int, default=10, help="Maximum experimental rank to consider")
    parser.add_argument(
        "--only-label",
        choices=["1", "-1", "any"],
        default="1",
        help="Restrict baseline/experimental rows by label; default keeps only target PSMs",
    )
    parser.add_argument(
        "--show-missing",
        type=int,
        default=10,
        help="How many missing baseline examples to print",
    )
    args = parser.parse_args()

    baseline_groups = load_pin_groups(args.baseline_pin, args.baseline_rank, args.only_label)
    experimental_groups = load_pin_groups(args.experimental_pin, args.experimental_rank, args.only_label)

    considered = 0
    recovered = 0
    missing_examples: list[tuple[str, list[str]]] = []

    for scan, baseline_peptides in sorted(baseline_groups.items(), key=lambda item: int(item[0])):
        if not baseline_peptides:
            continue
        considered += 1
        experimental_peptides = experimental_groups.get(scan, set())
        if baseline_peptides.intersection(experimental_peptides):
            recovered += 1
        elif len(missing_examples) < args.show_missing:
            missing_examples.append((scan, sorted(baseline_peptides)))

    recall = (recovered / considered * 100.0) if considered else 0.0

    print(f"baseline_pin={args.baseline_pin}")
    print(f"experimental_pin={args.experimental_pin}")
    print(f"baseline_rank<={args.baseline_rank}")
    print(f"experimental_rank<={args.experimental_rank}")
    print(f"label_filter={args.only_label}")
    print(f"spectra_considered={considered}")
    print(f"recovered={recovered}")
    print(f"missing={considered - recovered}")
    print(f"recall_pct={recall:.2f}")

    if missing_examples:
        print("")
        print("missing_examples:")
        for scan, peptides in missing_examples:
            print(f"  scan={scan} baseline_peptides={','.join(peptides)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
