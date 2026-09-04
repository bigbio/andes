#!/usr/bin/env python3
"""Summarize Phase V gate results and emit PASS/FAIL.

Usage: summarize_phase_v_gate.py <phase-v-output-dir>

Reads summary.tsv produced by phase_v_strong_score_gate.sh.
"""
from __future__ import annotations

import sys
from pathlib import Path


def load_summary(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        vals = line.split("\t")
        rows.append(dict(zip(header, vals)))
    return rows


def f(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def main() -> int:
    out_dir = Path(sys.argv[1])
    summary_path = out_dir / "summary.tsv"
    if not summary_path.exists():
        print(f"missing {summary_path}", file=sys.stderr)
        return 2

    rows = load_summary(summary_path)
    by_ds: dict[str, dict[str, dict[str, dict[str, str]]]] = {}
    for r in rows:
        ds = r["dataset"]
        mode = r["score_mode"]
        db = r["db"]
        by_ds.setdefault(ds, {}).setdefault(mode, {})[db] = r

    print(f"Phase V gate summary: {out_dir}")
    print()
    all_pass = True
    speed_pass = True

    for ds in sorted(by_ds):
        rank = by_ds[ds].get("rank", {})
        strong = by_ds[ds].get("strong", {})
        if not rank or not strong:
            print(f"[{ds}] INCOMPLETE — missing rank or strong runs")
            all_pass = False
            continue

        rn = rank.get("normal", {})
        sn = strong.get("normal", {})
        re_ = rank.get("entrap", {})
        se = strong.get("entrap", {})

        psms_rank = f(rn.get("psms_1pct", ""))
        psms_strong = f(sn.get("psms_1pct", ""))
        fdp_rank = f(re_.get("entrap_combined_fdp", ""))
        fdp_strong = f(se.get("entrap_combined_fdp", ""))
        wall_rank = f(rn.get("wall_s", ""))
        wall_strong = f(sn.get("wall_s", ""))

        psms_ok = psms_strong >= psms_rank
        fdp_ok = True
        if re_ and se:
            fdp_ok = fdp_strong <= fdp_rank + 1e-9
        speed_ok = True
        if wall_rank > 0:
            speed_ok = wall_strong <= 1.10 * wall_rank
            speed_pass = speed_pass and speed_ok

        ds_pass = psms_ok and fdp_ok and speed_ok
        all_pass = all_pass and ds_pass

        status = "PASS" if ds_pass else "FAIL"
        print(f"[{ds}] {status}")
        print(f"  PSMs@1%:  rank={psms_rank:.0f}  strong={psms_strong:.0f}  "
              f"({'≥' if psms_ok else '<'} rank)")
        if re_ and se:
            print(f"  FDP:      rank={fdp_rank:.4f}  strong={fdp_strong:.4f}  "
                  f"({'≤' if fdp_ok else '>'} rank)")
        else:
            print("  FDP:      entrap run missing")
            all_pass = False
        if wall_rank > 0:
            pct = 100.0 * wall_strong / wall_rank
            print(f"  Wall (s): rank={wall_rank:.1f}  strong={wall_strong:.1f}  "
                  f"({pct:.1f}% of rank, gate ≤110%)")
        print()

    if "ptm" not in by_ds:
        print("WARN: PTM dataset not run — gate requires Astral + TMT + PTM for full sign-off")
        all_pass = False

    print("=" * 60)
    if all_pass and speed_pass:
        print("OVERALL: PASS — strong beats rank on all completed datasets")
        return 0
    print("OVERALL: FAIL — do NOT flip default --score; investigate deltas")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
