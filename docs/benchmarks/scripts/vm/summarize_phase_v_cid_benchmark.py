#!/usr/bin/env python3
"""Summarize Phase V CID benchmark: MSFragger vs andes rank/strong."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_summary(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    lines = path.read_text().splitlines()
    if not lines:
        return rows
    header = lines[0].split("\t")
    for line in lines[1:]:
        if not line.strip():
            continue
        vals = line.split("\t")
        rows.append(dict(zip(header, vals)))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", type=Path, help="benchmark output directory")
    args = ap.parse_args()
    summary = args.out_dir / "summary.tsv"
    if not summary.exists():
        print(f"missing {summary}", file=sys.stderr)
        return 1

    rows = read_summary(summary)
    if not rows:
        print("empty summary")
        return 1

    print(f"# CID benchmark — {args.out_dir.name}\n")
    print("| Dataset | Engine | PSMs@1% | Peptides@1% |")
    print("|---------|--------|---------|-------------|")
    for r in sorted(rows, key=lambda x: (x.get("dataset", ""), x.get("engine", ""))):
        print(
            f"| {r.get('dataset','')} | {r.get('engine','')} | "
            f"{r.get('psms_1pct','')} | {r.get('peptides_1pct','')} |"
        )

    # Head-to-head: andes strong vs best external per dataset.
    by_ds: dict[str, list[dict[str, str]]] = {}
    for r in rows:
        by_ds.setdefault(r.get("dataset", ""), []).append(r)

    print("\n## Decision hints (andes vs MSFragger)\n")
    for ds, grp in sorted(by_ds.items()):
        def psms(engine: str) -> int:
            for r in grp:
                if r.get("engine") == engine:
                    return int(r.get("psms_1pct") or 0)
            return 0

        fragger = psms("msfragger")
        rank = psms("andes_rank")
        strong = psms("andes_strong")
        delta_strong = strong - fragger
        delta_rank = rank - fragger
        print(
            f"- **{ds}**: MSFragger {fragger:,}; "
            f"andes rank {rank:,} ({delta_rank:+,}); "
            f"andes strong {strong:,} ({delta_strong:+,})"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
