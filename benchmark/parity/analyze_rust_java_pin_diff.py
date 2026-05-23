#!/usr/bin/env python3
"""
Per-PSM Rust↔Java PIN diff harness.

Reads two .pin files (one from Java MS-GF+, one from msgf-rust) produced
on the same input mzML + FASTA with matching CLI flags, and emits a
report localizing where the two engines diverge:

  - Top-1-per-scan agreement / ranking-flip / Java-only / Rust-only buckets
  - Per-feature distributional diff on the agreement bucket
  - CSV of per-PSM diffs for further drill-down

Usage:
    python3 analyze_rust_java_pin_diff.py \\
        --java path/to/astral-java.pin \\
        --rust path/to/astral-rust.pin \\
        --out-dir path/to/output-dir

Python stdlib only; no pandas dependency.
"""

import argparse
import csv
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


# Feature columns to compare. Everything between RawScore (index 6) and
# matchedIonRatio (index 34) is numeric and a candidate for diff.
NUMERIC_COLS = [
    "RawScore", "DeNovoScore", "lnSpecEValue", "lnEValue", "isotope_error",
    "peplen", "dm", "absdm",
    "charge2", "charge3", "charge4",
    "enzN", "enzC", "enzInt",
    "NumMatchedMainIons", "longest_b", "longest_y", "longest_y_pct",
    "ExplainedIonCurrentRatio", "NTermIonCurrentRatio", "CTermIonCurrentRatio",
    "MS2IonCurrent", "IsolationWindowEfficiency",
    "MeanErrorTop7", "StdevErrorTop7", "MeanRelErrorTop7", "StdevRelErrorTop7",
    "lnDeltaSpecEValue", "matchedIonRatio",
]


def strip_flanking_and_mods(pep: str) -> str:
    """
    Strip Percolator flanking (`X.PEPTIDE.Y`) and mod-mass tokens like
    `+57.021` from a .pin-format peptide string. Returns the residue-only
    sequence in uppercase. Matches the Rust common::strip_flanking_and_mods
    behaviour at crates/search/tests/common/mod.rs.
    """
    if len(pep) < 5 or pep[1] != "." or pep[-2] != ".":
        return ""
    middle = pep[2:-2]
    out = []
    i = 0
    while i < len(middle):
        c = middle[i]
        if c == "+" or c == "-":
            # Consume digits and at most one decimal.
            i += 1
            while i < len(middle) and (middle[i].isdigit() or middle[i] == "."):
                i += 1
        elif c.isascii() and c.isupper():
            out.append(c)
            i += 1
        else:
            i += 1
    return "".join(out)


def parse_pin(path: Path):
    """
    Parse a .pin file. Returns:
        header (list[str])
        rows (list[dict]): each row has the column->str map plus
            `_scan`, `_peptide_residues`, `_label` precomputed.

    The Proteins column may have multiple tab-separated entries (one
    accession per matching protein); they are collected into a single
    list at row["Proteins"].
    """
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        # Locate fixed-position columns by name.
        scan_idx = header.index("ScanNr")
        label_idx = header.index("Label")
        pep_idx = header.index("Peptide")
        prot_idx = header.index("Proteins")
        # Anything from prot_idx onward is the multi-accession Proteins tail.
        rows = []
        for fields in reader:
            if len(fields) <= prot_idx:
                continue
            row = {h: fields[i] for i, h in enumerate(header) if i < prot_idx}
            row["Peptide"] = fields[pep_idx]
            row["Proteins"] = fields[prot_idx:]
            try:
                row["_scan"] = int(fields[scan_idx])
            except ValueError:
                continue
            try:
                row["_label"] = int(fields[label_idx])
            except ValueError:
                continue
            row["_peptide_residues"] = strip_flanking_and_mods(fields[pep_idx])
            rows.append(row)
    return header, rows


def best_by_lnspec(rows):
    """Pick the row with the smallest lnSpecEValue (best PSM). Returns None for empty."""
    best = None
    best_val = None
    for r in rows:
        try:
            v = float(r["lnSpecEValue"])
        except ValueError:
            continue
        if best_val is None or v < best_val:
            best_val = v
            best = r
    return best


def quantile(values, q):
    if not values:
        return float("nan")
    sv = sorted(values)
    pos = q * (len(sv) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sv) - 1)
    frac = pos - lo
    return sv[lo] * (1 - frac) + sv[hi] * frac


def feature_diff_stats(java_rows_by_psm, rust_rows_by_psm, agreement_keys):
    """
    For each (scan, peptide) key in the agreement bucket, compute the
    per-feature (rust - java) deltas. Return:
        stats: dict[col] -> dict of summary stats (n, mean, median, p5, p95,
                                                    mean_abs, mean_rel)
        per_row: list of (scan, peptide, {col: delta})
    """
    # Sentinel filter: drop rows where either side carries an out-of-range
    # placeholder (Rust uses i32::MIN ~ -2.1e9 for DeNovoScore when GF didn't
    # compute; Java filters those rows pre-emit so it never carries them).
    # Including them skews the mean Δ but tells us nothing about feature
    # formation — they're a "PSM was never enriched" signal, not a divergence.
    SENTINEL_THRESHOLD = -1e8

    per_row = []
    diffs_by_col = defaultdict(list)
    rel_diffs_by_col = defaultdict(list)
    sentinels_dropped = Counter()
    for key in agreement_keys:
        j = java_rows_by_psm[key]
        r = rust_rows_by_psm[key]
        deltas = {}
        for col in NUMERIC_COLS:
            try:
                jv = float(j[col])
                rv = float(r[col])
            except (KeyError, ValueError):
                continue
            if jv < SENTINEL_THRESHOLD or rv < SENTINEL_THRESHOLD:
                sentinels_dropped[col] += 1
                continue
            d = rv - jv
            deltas[col] = d
            diffs_by_col[col].append(d)
            if abs(jv) > 1e-12:
                rel_diffs_by_col[col].append(d / abs(jv))
        per_row.append((key[0], key[1], deltas))
    stats = {}
    for col in NUMERIC_COLS:
        vals = diffs_by_col[col]
        rvals = rel_diffs_by_col[col]
        if not vals:
            continue
        abs_vals = [abs(v) for v in vals]
        stats[col] = {
            "n": len(vals),
            "mean": statistics.mean(vals),
            "median": statistics.median(vals),
            "stdev": statistics.pstdev(vals) if len(vals) > 1 else 0.0,
            "p5": quantile(vals, 0.05),
            "p95": quantile(vals, 0.95),
            "mean_abs": statistics.mean(abs_vals),
            "mean_rel": statistics.mean(rvals) if rvals else float("nan"),
            "frac_diff_gt_1pct": (
                sum(1 for v in rvals if abs(v) > 0.01) / len(rvals)
                if rvals else float("nan")
            ),
        }
    return stats, per_row


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--java", required=True, type=Path)
    ap.add_argument("--rust", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Parsing {args.java}...", file=sys.stderr)
    _, java_rows = parse_pin(args.java)
    print(f"  Java rows: {len(java_rows)}", file=sys.stderr)

    print(f"Parsing {args.rust}...", file=sys.stderr)
    _, rust_rows = parse_pin(args.rust)
    print(f"  Rust rows: {len(rust_rows)}", file=sys.stderr)

    # ── group rows by scan ─────────────────────────────────────────────────
    java_by_scan = defaultdict(list)
    rust_by_scan = defaultdict(list)
    for r in java_rows:
        java_by_scan[r["_scan"]].append(r)
    for r in rust_rows:
        rust_by_scan[r["_scan"]].append(r)

    all_scans = sorted(set(java_by_scan) | set(rust_by_scan))

    # ── per-scan top-1 buckets ─────────────────────────────────────────────
    buckets = Counter()
    bucket_examples = defaultdict(list)
    # iter38 (diagnostic upgrade): collect ALL scans in the 3 non-converging
    # buckets along with both engines' top-1 rows. Written to non_converging.csv
    # to enable per-feature analysis of the residual PSM-divergence buckets.
    non_converging: list = []  # (scan, bucket, j_row, r_row) tuples
    for scan in all_scans:
        j = best_by_lnspec(java_by_scan.get(scan, []))
        r = best_by_lnspec(rust_by_scan.get(scan, []))
        if j is None and r is None:
            buckets["both_missing"] += 1
            continue
        if j is None:
            if r["_label"] == 1:
                buckets["rust_only_target"] += 1
                if len(bucket_examples["rust_only_target"]) < 5:
                    bucket_examples["rust_only_target"].append((scan, r["_peptide_residues"]))
            else:
                buckets["rust_only_decoy"] += 1
            continue
        if r is None:
            if j["_label"] == 1:
                buckets["java_only_target"] += 1
                if len(bucket_examples["java_only_target"]) < 5:
                    bucket_examples["java_only_target"].append((scan, j["_peptide_residues"]))
            else:
                buckets["java_only_decoy"] += 1
            continue
        # Both present.
        j_lab, r_lab = j["_label"], r["_label"]
        if j_lab == 1 and r_lab == 1:
            if j["_peptide_residues"] == r["_peptide_residues"]:
                buckets["both_target_same_peptide"] += 1
            else:
                buckets["both_target_diff_peptide"] += 1
                if len(bucket_examples["both_target_diff_peptide"]) < 5:
                    bucket_examples["both_target_diff_peptide"].append(
                        (scan, j["_peptide_residues"], r["_peptide_residues"])
                    )
                non_converging.append((scan, "both_target_diff_peptide", j, r))
        elif j_lab == 1 and r_lab == -1:
            buckets["java_target_rust_decoy"] += 1
            if len(bucket_examples["java_target_rust_decoy"]) < 5:
                bucket_examples["java_target_rust_decoy"].append(
                    (scan, j["_peptide_residues"], r["_peptide_residues"])
                )
            non_converging.append((scan, "java_target_rust_decoy", j, r))
        elif j_lab == -1 and r_lab == 1:
            buckets["rust_target_java_decoy"] += 1
            if len(bucket_examples["rust_target_java_decoy"]) < 5:
                bucket_examples["rust_target_java_decoy"].append(
                    (scan, j["_peptide_residues"], r["_peptide_residues"])
                )
            non_converging.append((scan, "rust_target_java_decoy", j, r))
        else:
            buckets["both_decoy"] += 1

    # ── per-feature diffs on (scan, peptide) full-join agreement set ───────
    java_by_psm = {}
    for r in java_rows:
        key = (r["_scan"], r["_peptide_residues"])
        # Keep best PSM per (scan, peptide) — same peptide can appear at multiple
        # ranks; we want the best score for the diff.
        prev = java_by_psm.get(key)
        if prev is None:
            java_by_psm[key] = r
        else:
            try:
                if float(r["lnSpecEValue"]) < float(prev["lnSpecEValue"]):
                    java_by_psm[key] = r
            except ValueError:
                pass

    rust_by_psm = {}
    for r in rust_rows:
        key = (r["_scan"], r["_peptide_residues"])
        prev = rust_by_psm.get(key)
        if prev is None:
            rust_by_psm[key] = r
        else:
            try:
                if float(r["lnSpecEValue"]) < float(prev["lnSpecEValue"]):
                    rust_by_psm[key] = r
            except ValueError:
                pass

    # Agreement subset: same (scan, peptide) AND same label on both sides
    # (typically Label=1; we restrict to target-target here).
    agreement_keys = []
    for key, j in java_by_psm.items():
        r = rust_by_psm.get(key)
        if r is None:
            continue
        if j["_label"] != 1 or r["_label"] != 1:
            continue
        agreement_keys.append(key)

    print(f"Agreement (scan, peptide, both target) PSMs: {len(agreement_keys)}",
          file=sys.stderr)

    stats, per_row = feature_diff_stats(java_by_psm, rust_by_psm, agreement_keys)

    # ── write CSV ──────────────────────────────────────────────────────────
    csv_path = args.out_dir / "per_psm_diff.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scan", "peptide"] + NUMERIC_COLS)
        for scan, pep, deltas in per_row:
            w.writerow([scan, pep] + [
                f"{deltas.get(c, ''):.6g}" if c in deltas else ""
                for c in NUMERIC_COLS
            ])
    print(f"Wrote {csv_path} ({len(per_row)} rows)", file=sys.stderr)

    # iter38 diagnostic upgrade: dump per-feature values for the
    # 3 non-converging buckets so future audits can characterize which
    # features are driving the divergence. Each row has both engines'
    # top-1 PSM values side-by-side for every NUMERIC_COL feature.
    nc_path = args.out_dir / "non_converging.csv"
    with open(nc_path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["scan", "bucket", "java_peptide", "rust_peptide"]
        for c in NUMERIC_COLS:
            header.append(f"java_{c}")
            header.append(f"rust_{c}")
        w.writerow(header)
        for scan, bucket, j, r in non_converging:
            row = [
                scan,
                bucket,
                j.get("_peptide_residues", ""),
                r.get("_peptide_residues", ""),
            ]
            for c in NUMERIC_COLS:
                jv = j.get(c, "")
                rv = r.get(c, "")
                row.append(jv)
                row.append(rv)
            w.writerow(row)
    print(f"Wrote {nc_path} ({len(non_converging)} non-converging PSMs)",
          file=sys.stderr)

    # ── write markdown report ──────────────────────────────────────────────
    md_path = args.out_dir / "report.md"
    with open(md_path, "w") as f:
        f.write("# Rust↔Java PIN diff report\n\n")
        f.write(f"- Java:  `{args.java}` ({len(java_rows)} rows, "
                f"{len(java_by_scan)} scans)\n")
        f.write(f"- Rust:  `{args.rust}` ({len(rust_rows)} rows, "
                f"{len(rust_by_scan)} scans)\n")
        f.write(f"- Total unique scans: {len(all_scans)}\n\n")

        f.write("## Top-1-per-scan buckets\n\n")
        total = sum(buckets.values())
        f.write("| Bucket | Count | % of total |\n")
        f.write("|---|---:|---:|\n")
        # Print buckets in a meaningful order.
        bucket_order = [
            "both_target_same_peptide",
            "both_target_diff_peptide",
            "java_target_rust_decoy",
            "rust_target_java_decoy",
            "java_only_target",
            "rust_only_target",
            "both_decoy",
            "java_only_decoy",
            "rust_only_decoy",
            "both_missing",
        ]
        for b in bucket_order:
            n = buckets.get(b, 0)
            pct = (n / total * 100.0) if total else 0.0
            f.write(f"| {b} | {n:,} | {pct:.2f}% |\n")
        f.write(f"| **total** | **{total:,}** | 100.00% |\n\n")

        # Bucket examples — useful for spot-checking
        if bucket_examples:
            f.write("### Sample disagreements\n\n")
            for b, samples in bucket_examples.items():
                f.write(f"**{b}** (first {len(samples)}):\n\n")
                for s in samples:
                    f.write(f"- `{s}`\n")
                f.write("\n")

        f.write("## Per-feature diff (agreement bucket: same scan + peptide, both target)\n\n")
        f.write(f"_{len(agreement_keys):,} PSMs in agreement bucket._\n\n")

        # Rank features by mean absolute delta. Skip those with n=0.
        ranked = sorted(stats.items(), key=lambda kv: -kv[1]["mean_abs"])
        f.write("Sorted by mean |Δ| (Rust - Java), descending:\n\n")
        f.write("| Feature | n | mean Δ | median Δ | stdev | p5 | p95 | mean \\|Δ\\| | mean rel Δ | %frac \\|relΔ\\|>1% |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for col, s in ranked:
            frac = s["frac_diff_gt_1pct"]
            frac_str = f"{frac*100:.1f}%" if not math.isnan(frac) else "—"
            f.write(
                f"| {col} | {s['n']:,} | "
                f"{s['mean']:+.4g} | {s['median']:+.4g} | {s['stdev']:.4g} | "
                f"{s['p5']:+.4g} | {s['p95']:+.4g} | "
                f"{s['mean_abs']:.4g} | "
                f"{s['mean_rel']:+.4g} | "
                f"{frac_str} |\n"
            )
        f.write("\n")

        f.write("## Notes\n\n")
        f.write("- Δ = (Rust value) - (Java value).\n")
        f.write("- `mean rel Δ` = mean of (Δ / |java|) over PSMs with |java| > 1e-12.\n")
        f.write("- `%frac |relΔ|>1%` = fraction of PSMs where the relative diff exceeds 1%.\n")
        f.write("- Agreement bucket restricts to scans + peptides present as target on BOTH sides; "
                "this strips ranking-flip + retention-only effects so the table measures FEATURE divergence specifically.\n")

    print(f"Wrote {md_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
