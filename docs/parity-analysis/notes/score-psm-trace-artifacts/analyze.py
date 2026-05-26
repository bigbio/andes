#!/usr/bin/env python3
"""One-shot I5 analysis: align Rust per-PSM JSON trace against Java per-scan
TRACE log for the 5 PXD001819 label-flip PSMs. Java trace has no peptide
attribution (it's per-(scan, nominal_mass, isPrefix, ion, theo_mz) — one
record per ion within a getNodeScore call). Rust JSON has per-PSM per-ion
records keyed by theo_mz.

For each Rust PSM ion, find Java's matching (ion_kind, theo_mz) within the
same scan with a 1e-3 Da tolerance. Tally divergences.
"""

import collections
import json
import os
import re
import struct
import sys

ART = "."
SCANS = [41522, 34685, 23272, 23082, 16629]


def normalize_rust_ion(s):
    """Rust IonType Debug -> 'b/<c>+<off:.5f>' or 'y/<c>+<off:.5f>' or 'Noise'."""
    s = s.strip()
    if "Noise" in s:
        return "Noise"
    m = re.match(r"(Prefix|Suffix)\s*\{\s*charge:\s*(\d+),\s*offset_bits:\s*(\d+)\s*\}", s)
    if m:
        kind = "b" if m.group(1) == "Prefix" else "y"
        c = int(m.group(2))
        off_bits = int(m.group(3))
        off = struct.unpack(">f", struct.pack(">I", off_bits))[0]
        return f"{kind}/{c}+{off:.5f}"
    return s


def normalize_java_ion(s):
    """Java 'b/<c>+<float>' -> 'b/<c>+<off:.5f>'."""
    m = re.match(r"([by])/(\d+)\+(-?[\d.]+)", s)
    if m:
        kind = m.group(1)
        c = int(m.group(2))
        off = float(m.group(3))
        return f"{kind}/{c}+{off:.5f}"
    return s


def load_rust(scan):
    path = f"{ART}/rust-trace-scan-{scan}.json"
    with open(path) as fh:
        data = json.load(fh)
    return data  # list of PSMs


def load_java(scan):
    """Return list of dicts per ion. Handles both .log and .log.gz."""
    import gzip
    base = f"{ART}/java-trace-scan-{scan}.log"
    if os.path.exists(base):
        fh = open(base)
    elif os.path.exists(base + ".gz"):
        fh = gzip.open(base + ".gz", "rt")
    else:
        raise FileNotFoundError(f"neither {base} nor {base}.gz")
    out = []
    with fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.startswith("TRACE"):
                continue
            fields = {}
            for part in line.split("\t")[1:]:
                if "=" in part:
                    k, v = part.split("=", 1)
                    fields[k] = v
            try:
                rec = {
                    "scan": int(fields["scan"]),
                    "nominalMass": int(fields["nominalMass"]),
                    "isPrefix": fields["isPrefix"] == "true",
                    "ion_kind": normalize_java_ion(fields["ion"]),
                    "theo_mz": float(fields["theo_mz"]),
                    "rank": int(fields["rank"]) if fields["rank"] != "-1" else None,
                    "log_prob": float(fields["log_prob"]),
                    "contribution": float(fields["contribution"]),
                }
            except (KeyError, ValueError):
                continue
            out.append(rec)
    return out


def index_java(java_ions, mz_tol=1e-3):
    """Index by (ion_kind, theo_mz_rounded). Multiple entries possible if
    Java emits the same nominal_mass repeatedly during scoring of different
    candidate peptides (values should be identical)."""
    idx = collections.defaultdict(list)
    for r in java_ions:
        key = (r["ion_kind"], round(r["theo_mz"] / mz_tol))
        idx[key].append(r)
    return idx


def compare_psm(psm, java_idx, mz_tol=1e-3):
    """Yields (ion_kind, theo_mz, rust, java_or_None, flags)."""
    rows = []
    for rust_ion in psm["ions"]:
        rkind = normalize_rust_ion(rust_ion["ion_type"])
        rkey = (rkind, round(rust_ion["theo_mz"] / mz_tol))
        candidates = java_idx.get(rkey, [])
        # Pick the first matching Java ion. (All should have the same numeric
        # values since they're per-(scan, nominal_mass, ion).)
        java_ion = candidates[0] if candidates else None
        flags = []
        if java_ion is None:
            flags.append("RUST_ONLY")
        else:
            if rust_ion.get("rank") != java_ion.get("rank"):
                flags.append("RANK_DIFF")
            if abs(rust_ion["log_prob"] - java_ion["log_prob"]) > 1e-3:
                flags.append("LOGPROB_DIFF")
            if abs(rust_ion["contribution"] - java_ion["contribution"]) > 1e-3:
                flags.append("CONTRIB_DIFF")
        rows.append((rkind, rust_ion["theo_mz"], rust_ion, java_ion, flags))
    return rows


def fmt_num(v, prec):
    return f"{v:>{8+prec}.{prec}f}" if v is not None else "-" * (8 + prec)


def main():
    summary = []
    for scan in SCANS:
        rust_psms = load_rust(scan)
        java_ions = load_java(scan)
        java_idx = index_java(java_ions)
        print(f"\n{'=' * 78}\nSCAN {scan}  |  Rust PSMs traced: {len(rust_psms)}  |  Java ions: {len(java_ions)}")
        for psm in rust_psms:
            pep = psm["peptide"]
            rscore = psm["rust_rank_score"]
            print(f"\n  PSM: peptide={pep}  charge={psm['charge']}  rust_rank_score={rscore}")
            rows = compare_psm(psm, java_idx)
            rust_total = sum(r[2]["contribution"] for r in rows)
            java_matched = sum(r[3]["contribution"] for r in rows if r[3] is not None)
            divergences = collections.Counter()
            for kind, mz, rust, java, flags in rows:
                for f in flags:
                    divergences[f] += 1
            print(f"    ions: {len(rows)} (rust-only: {divergences.get('RUST_ONLY', 0)})")
            print(f"    rust contribution sum:  {rust_total:>10.4f}")
            print(f"    java contribution sum:  {java_matched:>10.4f}  (matched ions only)")
            print(f"    delta (rust - java):    {rust_total - java_matched:>+10.4f}")
            print(f"    divergence counts: {dict(divergences)}")
            summary.append((scan, pep, rscore, len(rows), divergences))

    # Aggregate across all 5 scans / 10 PSMs
    print("\n" + "=" * 78)
    print("AGGREGATE (5 scans x ~2 PSMs each):")
    total_div = collections.Counter()
    for scan, pep, rscore, nions, divs in summary:
        total_div.update(divs)
    print(f"  Total divergences across all traced PSMs:")
    for cat, count in total_div.most_common():
        print(f"    {cat}: {count}")


if __name__ == "__main__":
    main()
