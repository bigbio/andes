#!/usr/bin/env python3
"""Extract benchmark metrics from GNU time output + MS-GF+ mzIdentML.

Uses xml.etree.ElementTree.iterparse to stream mzIdentML (files can be
hundreds of MB) and count SpectrumIdentificationItem elements and the
subset with PSM-level Q-value (MS:1002054) <= 0.01.
"""
from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path

PSM_QVALUE_ACCESSION = "MS:1002054"
PSM_QVALUE_THRESHOLD = 0.01

_NS_RE = re.compile(r"^\{[^}]+\}")


def _localname(tag: str) -> str:
    return _NS_RE.sub("", tag)


def parse_gnu_time(path: Path) -> tuple[str, str]:
    text = path.read_text(errors="replace")
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    cpu = re.search(r"Percent of CPU this job got: (\d+)", text)
    return (rss.group(1) if rss else "NA", cpu.group(1) if cpu else "NA")


def parse_mzid(path: Path) -> tuple[int, int]:
    """Return (sii_count, psm_1pct_fdr_count) via streaming iterparse."""
    sii_count = 0
    psm_1pct = 0

    context = ET.iterparse(str(path), events=("end",))
    for _, elem in context:
        if _localname(elem.tag) != "SpectrumIdentificationItem":
            continue
        sii_count += 1
        for child in elem:
            if _localname(child.tag) != "cvParam":
                continue
            if child.get("accession") != PSM_QVALUE_ACCESSION:
                continue
            value = child.get("value", "")
            try:
                if float(value) <= PSM_QVALUE_THRESHOLD:
                    psm_1pct += 1
            except ValueError:
                pass
            break
        elem.clear()
    return sii_count, psm_1pct


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--time", type=Path, required=True, help="GNU time -v output")
    ap.add_argument("--mzid", type=Path, required=True, help="MS-GF+ mzIdentML output")
    ap.add_argument("--wall", type=int, required=True, help="Wall-clock seconds (int)")
    ap.add_argument("--output", type=Path, required=True, help="Destination key=value file")
    args = ap.parse_args()

    rss_kb, cpu_pct = parse_gnu_time(args.time)
    sii_count, psm_1pct = parse_mzid(args.mzid)

    lines = [
        "dataset=PXD001819",
        f"wall_time_sec={args.wall}",
        f"sii_count={sii_count}",
        f"psm_1pct_fdr={psm_1pct}",
        f"peak_rss_kb={rss_kb}",
        f"cpu_percent={cpu_pct}",
    ]
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.output.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
