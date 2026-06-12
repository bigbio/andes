#!/usr/bin/env python3
"""
Diff per-PSM per-ion trace outputs from Rust (msgf-trace --trace-json) and
Java (instrumented java-legacy stderr). For each (scan, peptide) PSM, align
records by (ion_kind, theo_mz tolerance 1e-3 Da) and emit a side-by-side
table.

Usage:
    diff_score_psm_traces.py --rust rust-trace.json --java java-trace.log \\
        [--mz-tol 1e-3] [--scan SCAN] [--peptide PEP]

Outputs to stdout. Exit code 0 = success.

Rust JSON shape (per PSM):
    {
      "scan": int,
      "peptide": str,
      "charge": int,
      "rust_rank_score": int,
      "ions": [
        {"ion_type": str, "theo_mz": float, "rank": int|null,
         "max_rank": int, "log_prob": float, "contribution": float},
        ...
      ]
    }

Java log shape (one line per ion):
    TRACE\\tscan=<int>\\tpeptide=<str>\\tion=<str>\\ttheo_mz=<float>\\trank=<int>\\tlog_prob=<float>\\tcontribution=<float>

Java represents a missing rank as rank=-1 (Rust uses null).
"""

import argparse
import collections
import json
import re
import struct
import sys


def normalize_ion_kind(s: str) -> str:
    """Map both Rust and Java ion-type representations to a normalized key.

    Rust format: `Prefix { charge: 1, offset_bits: 0 }`
    Java format: `b/1+0.00000`
    Normalize to: `b/<charge>+<offset:.5f>` or `y/<charge>+<offset:.5f>` or `Noise`.
    """
    s = s.strip()
    if "Noise" in s:
        return "Noise"
    # Rust format
    rust_match = re.match(
        r"(Prefix|Suffix)\s*\{\s*charge:\s*(\d+),\s*offset_bits:\s*(\d+)\s*\}",
        s,
    )
    if rust_match:
        kind = "b" if rust_match.group(1) == "Prefix" else "y"
        charge = int(rust_match.group(2))
        off_bits = int(rust_match.group(3))
        off = struct.unpack(">f", struct.pack(">I", off_bits))[0]
        return f"{kind}/{charge}+{off:.5f}"
    # Java format
    java_match = re.match(r"([by])/(\d+)\+([\d.+\-eE]+)", s)
    if java_match:
        kind = java_match.group(1)
        charge = int(java_match.group(2))
        off = float(java_match.group(3))
        return f"{kind}/{charge}+{off:.5f}"
    return s


def parse_rust_json(path: str) -> dict:
    """Returns {(scan, peptide): [{ion fields}, ...]}."""
    out = {}
    with open(path) as fh:
        data = json.load(fh)
    for psm in data:
        key = (psm["scan"], psm["peptide"])
        out[key] = psm["ions"]
    return out


def parse_java_log(path: str) -> dict:
    """Returns {(scan, peptide): [{ion fields}, ...]}."""
    out = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.startswith("TRACE\t"):
                continue
            fields = {}
            for part in line.split("\t")[1:]:
                if "=" not in part:
                    continue
                k, v = part.split("=", 1)
                fields[k] = v
            try:
                scan = int(fields["scan"])
                peptide = fields["peptide"]
                raw_rank = fields.get("rank", "")
                rank = None if raw_rank in ("", "-1", "null") else int(raw_rank)
                ion = {
                    "ion_type": fields.get("ion", "?"),
                    "theo_mz": float(fields.get("theo_mz", "nan")),
                    "rank": rank,
                    "log_prob": float(fields.get("log_prob", "nan")),
                    "contribution": float(fields.get("contribution", "nan")),
                }
            except (KeyError, ValueError) as e:
                print(
                    f"WARN: skipping malformed Java TRACE line: {line[:80]}... ({e})",
                    file=sys.stderr,
                )
                continue
            out[(scan, peptide)].append(ion)
    return out


def align_and_diff(rust_ions, java_ions, mz_tol):
    """Yields (key, rust_ion_or_None, java_ion_or_None, flags) per ion."""
    java_by_key = collections.defaultdict(list)
    for ion in java_ions:
        key = (normalize_ion_kind(ion["ion_type"]), round(ion["theo_mz"] / mz_tol))
        java_by_key[key].append(ion)

    matched_java_ids = set()
    for rust_ion in rust_ions:
        rust_key = (
            normalize_ion_kind(rust_ion["ion_type"]),
            round(rust_ion["theo_mz"] / mz_tol),
        )
        candidates = java_by_key.get(rust_key, [])
        java_ion = candidates.pop(0) if candidates else None
        if java_ion is not None:
            matched_java_ids.add(id(java_ion))
        flags = []
        if java_ion is None:
            flags.append("RUST_ONLY")
        else:
            if rust_ion.get("rank") != java_ion.get("rank"):
                flags.append("RANK_DIFF")
            if abs(rust_ion["log_prob"] - java_ion["log_prob"]) > 1e-4:
                flags.append("LOGPROB_DIFF")
            if abs(rust_ion["contribution"] - java_ion["contribution"]) > 1e-4:
                flags.append("CONTRIB_DIFF")
        yield (rust_key, rust_ion, java_ion, flags)

    for ion in java_ions:
        if id(ion) in matched_java_ids:
            continue
        key = (normalize_ion_kind(ion["ion_type"]), round(ion["theo_mz"] / mz_tol))
        yield (key, None, ion, ["JAVA_ONLY"])


def format_row(key, rust_ion, java_ion, flags):
    def fmt(v, w, prec=None):
        if v is None:
            return "-" * w
        if isinstance(v, float) and prec is not None:
            return f"{v:>{w}.{prec}f}"
        return f"{str(v):>{w}}"

    theo_mz = (rust_ion or java_ion)["theo_mz"]
    return "  ".join([
        fmt(key[0], 22),
        fmt(theo_mz, 10, prec=4),
        fmt(rust_ion.get("rank") if rust_ion else None, 5),
        fmt(java_ion.get("rank") if java_ion else None, 5),
        fmt(rust_ion["log_prob"] if rust_ion else None, 9, prec=4),
        fmt(java_ion["log_prob"] if java_ion else None, 9, prec=4),
        fmt(rust_ion["contribution"] if rust_ion else None, 9, prec=4),
        fmt(java_ion["contribution"] if java_ion else None, 9, prec=4),
        ",".join(flags) if flags else "",
    ])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--rust",
        required=True,
        help="Rust trace JSON from msgf-trace --trace-json",
    )
    ap.add_argument(
        "--java",
        required=True,
        help="Java instrumented trace log (TRACE lines)",
    )
    ap.add_argument(
        "--mz-tol",
        type=float,
        default=1e-3,
        help="m/z alignment tolerance (Da, default 1e-3)",
    )
    ap.add_argument(
        "--scan",
        type=int,
        default=None,
        help="Restrict to one scan",
    )
    ap.add_argument(
        "--peptide",
        default=None,
        help="Restrict to one peptide",
    )
    args = ap.parse_args()

    rust = parse_rust_json(args.rust)
    java = parse_java_log(args.java)

    all_keys = sorted(set(rust.keys()) | set(java.keys()))
    for key in all_keys:
        scan, pep = key
        if args.scan is not None and scan != args.scan:
            continue
        if args.peptide is not None and pep != args.peptide:
            continue
        print(f"\n=== scan={scan} peptide={pep} ===")
        rust_ions = rust.get(key, [])
        java_ions = java.get(key, [])
        if not rust_ions and not java_ions:
            print("  (no data on either side)")
            continue
        print(
            "  ion_type                theo_mz     R_rk   J_rk    R_logP    J_logP    R_ctrb     J_ctrb    flags"
        )
        rust_total = 0.0
        java_total = 0.0
        category_counts = collections.Counter()
        for row in align_and_diff(rust_ions, java_ions, args.mz_tol):
            print("  " + format_row(*row))
            if row[1] is not None:
                rust_total += row[1]["contribution"]
            if row[2] is not None:
                java_total += row[2]["contribution"]
            for f in row[3]:
                category_counts[f] += 1
        print(
            f"  TOTAL contribution: rust={rust_total:.4f}  java={java_total:.4f}  "
            f"delta={rust_total - java_total:+.4f}"
        )
        if category_counts:
            print(f"  DIVERGENCES: {dict(category_counts)}")


if __name__ == "__main__":
    main()
