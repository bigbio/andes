#!/usr/bin/env python3
"""Per-scan top-1 PSM agreement between two .pin files.

Usage: diff_pin_psms.py pre.pin post.pin
Exits 0 if agreement >= 99.9% (<= 37 disagreements out of 37,112).
"""
import csv, sys

def load(p):
    with open(p) as f:
        return {r['ScanNr']: r for r in csv.DictReader(f, delimiter='\t')}

def main(pre, post):
    a, b = load(pre), load(post)
    common = set(a) & set(b)
    disagree = [s for s in common if a[s]['Peptide'] != b[s]['Peptide']]
    only_a = set(a) - set(b); only_b = set(b) - set(a)
    n_total = len(common)
    n_dis = len(disagree)
    pct_agree = 100 * (1 - n_dis / n_total) if n_total else 0
    print(f"Scans in both: {n_total}")
    print(f"Top-1 peptide disagreements: {n_dis} ({100-pct_agree:.3f}%)")
    print(f"Pre-only scans: {len(only_a)}")
    print(f"Post-only scans: {len(only_b)}")
    print(f"Agreement: {pct_agree:.3f}% (gate: >= 99.9%)")
    for s in disagree[:5]:
        print(f"  scan={s}: pre={a[s]['Peptide']!r}  post={b[s]['Peptide']!r}")
    return 0 if pct_agree >= 99.9 else 1

if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
