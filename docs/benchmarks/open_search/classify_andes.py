#!/usr/bin/env python3
"""Classify an andes open-search PSM list into substitution-containing / -free.

The paper's comparison (Yu et al., Nat. Commun. 2020, Fig. 3) splits the PSMs of
each tool into those whose true peptide carries a substitution and those that do
not. This script does the same for an andes run against the substituted
database: a PSM is *containing* when its peptide equals a span of the mutated
protein that differs from the plain protein at exactly one residue, *free* when
it equals the plain span, and *other* otherwise.

    python3 -I classify_andes.py --dir OUT --plain plain.fasta --mutated mut.fasta \
        --out report.json
"""
import argparse
import json
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze as A  # noqa: E402


def read_fasta(path):
    out = {}
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def classify(plain, mutated, peptide, accessions):
    """'containing', 'free' or None for one PSM."""
    best = None
    for acc in accessions:
        p, m = plain.get(acc), mutated.get(acc)
        if not p or not m or len(p) != len(m):
            continue
        start = m.find(peptide)
        while start >= 0:
            subs = sum(1 for i in range(len(peptide)) if p[start + i] != m[start + i])
            if subs == 0:
                return "free"
            if subs == 1:
                best = "containing"
            start = m.find(peptide, start + 1)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--plain", required=True)
    ap.add_argument("--mutated", required=True)
    ap.add_argument("--q", type=float, default=0.01)
    ap.add_argument(
        "--truth",
        default=None,
        help="truth_by_scan.json; adds the Type-1 (matches the reference truth) fraction",
    )
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    plain = read_fasta(a.plain)
    mutated = read_fasta(a.mutated)
    psms = A.read_percolator(os.path.join(a.dir, "open.psms"))
    ok = {k: v for k, v in psms.items() if v[0] <= a.q}
    truth = json.load(open(a.truth)) if a.truth else None
    truth_peps = {t["mutant"] for t in truth.values()} if truth else set()
    counts = {"containing": 0, "free": 0, "other": 0}
    type1 = {"containing": 0, "free": 0}
    type2 = {"containing": 0, "free": 0}
    for spec, (_q, peptide, prots) in ok.items():
        bare, _mods = A.parse_peptide(peptide)
        accs = [x for cell in prots for x in cell.split(";") if x]
        kind = classify(plain, mutated, bare, accs)
        counts[kind or "other"] += 1
        if truth is not None and kind in type1:
            import re
            m = re.search(r"scan=(\d+)", spec)
            t = truth.get(m.group(1)) if m else None
            # Type 1: the PSM's peptide is the reference truth's peptide.
            if t is not None and t["mutant"] == bare:
                type1[kind] += 1
            elif bare in truth_peps:
                # Type 2: a different sequence, but one the reference found
                # elsewhere.
                type2[kind] += 1
    report = {"q": a.q, "psms": len(ok), **counts}
    if truth is not None:
        report["type1"] = type1
        report["type2"] = type2
        report["type1_fraction_containing"] = (
            type1["containing"] / counts["containing"] if counts["containing"] else 0.0
        )
        report["precision_containing"] = (
            (type1["containing"] + type2["containing"]) / counts["containing"]
            if counts["containing"]
            else 0.0
        )
    with open(a.out, "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
