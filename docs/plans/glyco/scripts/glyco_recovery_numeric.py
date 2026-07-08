#!/usr/bin/env python3
"""Parser-bug-free @1% backbone-correct recovery.

glyco_recovery_fdr.py matches truth by re-parsing the Percolator peptide string
(residue_mass), which mis-parses ~46% of multi-mod/long peptides. This joins each
@q<=thr surviving PSM back to its PIN row and matches numerically via
CalcMass - GlycanMass - H2O (the authoritative backbone convention).

Usage: glyco_recovery_numeric.py <truth.tsv> <percolator_psms> <pin> [q=0.01] [tol=0.05]
"""
import csv, re, sys
H2O = 18.010565
truth_f, psms_f, pin_f = sys.argv[1], sys.argv[2], sys.argv[3]
Q = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
TOL = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05

truth = {int(float(r["scan"])): float(r["backbone_mass"]) + H2O
         for r in csv.DictReader(open(truth_f), delimiter="\t")}

# PIN: SpecId -> backbone_neutral (CalcMass - GlycanMass). Target rows only.
H = open(pin_f).readline().rstrip("\n").split("\t"); idx = {h: i for i, h in enumerate(H)}
def g(p, n):
    try: return float(p[idx[n]])
    except (KeyError, ValueError): return 0.0
bb_by_spec = {}
for line in open(pin_f):
    p = line.rstrip("\n").split("\t")
    if len(p) < len(H) or p[0] == "SpecId": continue
    if int(g(p, "Label")) < 0: continue  # skip decoys
    bb_by_spec[p[0]] = g(p, "CalcMass") - g(p, "GlycanMass")

# Percolator psms (tab): PSMId score q-value posterior_error_prob peptide proteinIds
scan_re = re.compile(r"scan=(\d+)")
survived_scans = set(); correct_scans = set(); total_surv = 0; matched_pin = 0
with open(psms_f) as fh:
    rd = csv.reader(fh, delimiter="\t"); header = next(rd)
    hi = {h: i for i, h in enumerate(header)}
    qcol = hi.get("q-value", 2); idcol = hi.get("PSMId", 0)
    for row in rd:
        if len(row) <= qcol: continue
        try: q = float(row[qcol])
        except ValueError: continue
        if q > Q: continue
        total_surv += 1
        spec = row[idcol]; m = scan_re.search(spec)
        if not m: continue
        scan = int(m.group(1)); survived_scans.add(scan)
        bbn = bb_by_spec.get(spec)
        if bbn is None: continue
        matched_pin += 1
        tn = truth.get(scan)
        if tn is not None and abs(bbn - tn) <= TOL:
            correct_scans.add(scan)

print(f"total target PSMs @ q<={Q}          : {total_surv}")
print(f"  joined to a PIN target row        : {matched_pin}")
print(f"truth scans with a PSM @ q<={Q}     : {len(survived_scans & set(truth))}")
print(f"  backbone-correct (NUMERIC, tol {TOL}): {len(correct_scans)}")
