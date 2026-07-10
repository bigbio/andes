#!/usr/bin/env python3
"""Parser-bug-free @1% backbone-correct recovery.

glyco_recovery_fdr.py matches truth by re-parsing the Percolator peptide string
(residue_mass), which mis-parses ~46% of multi-mod/long peptides. This joins each
@q<=thr surviving PSM back to its PIN row and matches numerically via
CalcMass - GlycanMass - H2O (the authoritative backbone convention).

Usage: glyco_recovery_numeric.py <truth.tsv> <percolator_psms> <pin> [q=0.01] [tol=0.05]
"""
import csv, math, re, sys
H2O = 18.010565

USAGE = ("usage: glyco_recovery_numeric.py <truth.tsv> <percolator_psms> <pin> "
         "[q=0.01] [tol=0.05]")
if len(sys.argv) < 4:
    sys.exit(USAGE)
truth_f, psms_f, pin_f = sys.argv[1], sys.argv[2], sys.argv[3]
try:
    Q = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
    TOL = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05
except ValueError:
    sys.exit(USAGE)

# Truth backbone_mass MUST be a RESIDUE mass (no water); we add one H2O to
# compare against the PIN neutral peptide mass (CalcMass - GlycanMass).
truth = {}
with open(truth_f, newline="") as fh:
    trd = csv.DictReader(fh, delimiter="\t")
    if trd.fieldnames is None or "scan" not in trd.fieldnames or "backbone_mass" not in trd.fieldnames:
        sys.exit(f"{truth_f}: truth TSV needs 'scan' and 'backbone_mass' columns")
    for r in trd:
        s, m = r.get("scan"), r.get("backbone_mass")
        if not s or not m:
            continue
        try:
            scan = int(float(s)); mass = float(m)
        except ValueError:
            continue
        if not math.isfinite(mass):
            continue
        truth[scan] = mass + H2O

# PIN: SpecId -> backbone_neutral (CalcMass - GlycanMass). Target rows only.
with open(pin_f) as fh:
    H = fh.readline().rstrip("\n").split("\t")
idx = {h: i for i, h in enumerate(H)}
_required = ("SpecId", "Label", "CalcMass", "GlycanMass")
_missing = [c for c in _required if c not in idx]
if _missing:
    sys.exit(f"{pin_f}: PIN missing required columns: {', '.join(_missing)}")

def gnum(p, n):
    """Required numeric PIN field; None if missing/malformed/non-finite."""
    try:
        v = float(p[idx[n]])
    except (IndexError, ValueError):
        return None
    return v if math.isfinite(v) else None

bb_by_spec = {}
with open(pin_f) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(H) or p[idx["SpecId"]] == "SpecId":
            continue
        label = gnum(p, "Label")
        if label is None or label < 0:  # skip decoys and malformed-label rows
            continue
        # CalcMass and GlycanMass are an atomic pair; require both (never 0-default).
        calc = gnum(p, "CalcMass")
        gly = gnum(p, "GlycanMass")
        if calc is None or gly is None:
            continue
        bb_by_spec[p[idx["SpecId"]]] = calc - gly

# Percolator psms (tab): PSMId score q-value posterior_error_prob peptide proteinIds
scan_re = re.compile(r"scan=(\d+)")
survived_scans = set(); correct_scans = set(); total_surv = 0; matched_pin = 0
with open(psms_f) as fh:
    rd = csv.reader(fh, delimiter="\t"); header = next(rd)
    hi = {h: i for i, h in enumerate(header)}
    if "q-value" not in hi or "PSMId" not in hi:
        sys.exit(f"{psms_f}: percolator output needs 'q-value' and 'PSMId' columns")
    qcol = hi["q-value"]; idcol = hi["PSMId"]
    for row in rd:
        if len(row) <= max(qcol, idcol): continue
        try: q = float(row[qcol])
        except ValueError: continue
        if not math.isfinite(q): continue
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
