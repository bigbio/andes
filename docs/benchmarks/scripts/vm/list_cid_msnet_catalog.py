#!/usr/bin/env python3
"""List MSNet catalog entries useful for CID rank training / Lumos hold-out."""
import json
import sys
from pathlib import Path

CATALOG = Path(__file__).resolve().parents[3] / "internal-docs/msnet-model-catalog/msnet_model_catalog.jsonl"


def has_cid(frag: str) -> bool:
    return "CID" in [t.strip().upper() for t in (frag or "").split(";") if t.strip()]


def has_tmt(mods: str, ptm: str) -> bool:
    return "TMT" in f"{mods} {ptm}".upper()


def main() -> None:
    rows = [json.loads(line) for line in CATALOG.read_text().splitlines() if line.strip()]
    cid = [d for d in rows if has_cid(d.get("fragmentation_methods"))]
    print("# CID MSNet training candidates (true CID token)\n")
    print("| dataset | PSMs | instrument | TMT | url |")
    print("|---------|-----:|------------|-----|-----|")
    for d in sorted(cid, key=lambda x: -(x.get("psm_rows") or 0)):
        print(
            f"| {d['dataset_id']} | {d.get('psm_rows', 0):,} | "
            f"{(d.get('instruments') or '')[:28]} | "
            f"{'yes' if has_tmt(d.get('fixed_mods', ''), d.get('ptm_classes', '')) else 'no'} | "
            f"{d.get('msnet_parquet_url', '')} |"
        )
    lumos = [
        d for d in rows
        if "Lumos" in (d.get("instruments") or "") and has_tmt(d.get("fixed_mods", ""), d.get("ptm_classes", ""))
    ]
    print("\n# Lumos + TMT in MSNet (no true CID in catalog today)\n")
    if not lumos:
        print("(none — PXD007683 has no MSNet parquet; use PRIDE hold-out only)")
    else:
        for d in lumos:
            print(f"- {d['dataset_id']}: {d.get('fragmentation_methods')} {d.get('msnet_parquet_url')}")
    print("\n# PRIDE hold-out (benchmark, not MSNet train)\n")
    print("- PXD007683 — Lumos TMT CID (primary gap)")
    print("- PXD016999 — Fusion TMT CID (regression guard)")


if __name__ == "__main__":
    main()
