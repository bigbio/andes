#!/usr/bin/env python3
"""Score the simulated-substitution benchmark (Yu et al., Nat. Commun. 2020).

The reference search (plain database) gives the truth peptides; the modified
database carries random amino-acid substitutions. For every reference PSM whose
span has exactly ONE substituted residue, the open search should report the
substituted peptide with an offset equal to the substitution's mass difference
(original - substituted) and the substituted residue as the site.

    python3 -I analyze_substitution.py --dir OUT --plain plain.fasta \
        --mutated modified.fasta --out OUT/report.json
"""
import argparse
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze as A  # noqa: E402

H2O = 18.01056
CAM = 57.021464
RES_MASS = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406,
    "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
    "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841,
    "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}


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


def span_mass(seq):
    """Neutral monoisotopic mass of a peptide span with carbamidomethyl C."""
    return sum(RES_MASS[a] for a in seq) + H2O + seq.count("C") * CAM


def find_truth(plain, mutated, peptide, accessions):
    """(mutant span, delta, 0-based position) for a span with one substitution.

    `delta` is the offset the open search sees: the reference peptide's measured
    neutral mass (CalcMass, variable mods included) minus the substituted span's
    mass. Returns None when the span has zero or several substitutions, or the
    peptide does not occur in the named proteins.
    """
    for acc in accessions:
        p, m = plain.get(acc), mutated.get(acc)
        if not p or not m or len(p) != len(m):
            continue
        start = p.find(peptide)
        while start >= 0:
            sub = [i for i in range(len(peptide)) if p[start + i] != m[start + i]]
            if len(sub) == 1:
                mut = m[start:start + len(peptide)]
                return mut, sub[0]
            start = p.find(peptide, start + 1)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--plain", required=True)
    ap.add_argument("--mutated", required=True)
    ap.add_argument("--q", type=float, default=0.01)
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument(
        "--truth-json",
        default=None,
        help="also write the truth keyed by mzML scan number, for external-tool comparisons",
    )
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    d = a.dir
    plain = read_fasta(a.plain)
    mutated = read_fasta(a.mutated)
    ref_pin = A.read_pin(os.path.join(d, "reference.pin"))
    open_pin = A.read_pin(os.path.join(d, "open.pin"))
    ref_psms = A.read_percolator(os.path.join(d, "reference.psms"))
    open_psms = A.read_percolator(os.path.join(d, "open.psms"))
    ref_ok = {k: v for k, v in ref_psms.items() if v[0] <= a.q}
    open_ok = {k: v for k, v in open_psms.items() if v[0] <= a.q}

    truth = {}
    multi = 0
    variable = 0
    for spec_id, (_q, peptide, prots) in ref_ok.items():
        bare, mods = A.parse_peptide(peptide)
        # The simulation is on unmodified spectra: keep only peptides whose
        # only modification is the fixed carbamidomethyl on C.
        if any(not (bare[pos] == "C" and abs(delta - CAM) < 0.01) for pos, delta in mods):
            variable += 1
            continue
        row = A.pick(ref_pin, spec_id, peptide)
        if row is None:
            continue
        accs = [x for cell in prots for x in cell.split(";") if x]
        hit = find_truth(plain, mutated, bare, accs)
        if hit is None:
            multi += 1
            continue
        mut, pos = hit
        truth[spec_id] = {
            "mutant": mut,
            "delta": row["calc"] - span_mass(mut),
            "site": pos + 1,
        }

    recovered_pep = recovered = localized_called = localized_ok = 0
    errors = []
    rows = []
    for spec_id, t in truth.items():
        hit = open_ok.get(spec_id)
        if hit is None:
            continue
        bare, _ = A.parse_peptide(hit[1])
        if bare != t["mutant"]:
            continue
        recovered_pep += 1
        row = A.pick(open_pin, spec_id, hit[1])
        if row is None:
            continue
        delta = A.pin_delta(row)
        if abs(delta - t["delta"]) > a.tol:
            continue
        recovered += 1
        errors.append(abs(delta - t["delta"]))
        site = row["open_site"]
        if site != 0:
            localized_called += 1
            if site == t["site"]:
                localized_ok += 1
        rows.append({"spec": spec_id, "delta": delta, "truth": t["delta"],
                     "site": site, "truth_site": t["site"]})

    if a.truth_json:
        import re
        by_scan = {}
        for spec_id, t in truth.items():
            m = re.search(r"scan=(\d+)", spec_id)
            if m:
                by_scan[m.group(1)] = t
        with open(a.truth_json, "w") as fh:
            json.dump(by_scan, fh)

    report = {
        "q": a.q,
        "reference_psms": len(ref_ok),
        "truth_cases": len(truth),
        "skipped_variable_mod_reference": variable,
        "skipped_zero_or_multiple_substitutions": multi,
        "open_psms": len(open_ok),
        "recovered_peptide": recovered_pep,
        "recovered_offset": recovered,
        "recovered_offset_rate": recovered / len(truth) if truth else 0.0,
        "localized_called": localized_called,
        "localized_correct": localized_ok,
        "localized_accuracy": (localized_ok / localized_called) if localized_called else 0.0,
        "offset_abs_error_median": statistics.median(errors) if errors else 0.0,
        "rows": rows[:100],
    }
    with open(a.out, "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
