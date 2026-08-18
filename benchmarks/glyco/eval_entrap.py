#!/usr/bin/env python3
"""Measure the entrapment false-discovery proportion of a glyco search.

`build_entrap.py` appends an unrelated proteome to the search database as TARGETS.
Those proteins are not in the sample, and prokaryotes/yeast lack the mammalian
N-glycosylation machinery, so any accepted glycopeptide on an `ENTRAP_` protein is
false by construction. Counting them gives an absolute error estimate that does not
depend on a reference identification set, and that no amount of search-space expansion
can flatter.

This is the check that must accompany `eval_yield.py`. Yield alone rises whenever the
candidate space grows, real or not: the full 4034-composition glycan list looked like a
gain of 59 compositions on yield and turned out to inflate this number 5.4x.

FDP is reported as entrap / (entrap + sample), which is the fraction of accepted IDs
that are demonstrably wrong. Compare it against the nominal q-value threshold: if FDP is
at or below nominal, the reported FDR is conservative; well above it, the FDR is
optimistic and the IDs should not be trusted.

Usage:
    eval_entrap.py <pooled.pin> <percolator.psms> [q-threshold, default 0.01]
"""
import sys
from collections import defaultdict

if len(sys.argv) < 3:
    sys.exit(__doc__)
pin, psms = sys.argv[1], sys.argv[2]
qmax = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01

# SpecId -> accessions. Proteins occupy every column after "Peptide".
prot_of = {}
header = None
with open(pin) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = f
            i_pep = header.index("Peptide")
            continue
        prot_of[f[0]] = f[i_pep + 1:]

n_entrap = n_sample = 0
per_charge = defaultdict(lambda: [0, 0])
with open(psms) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] == "PSMId" or len(f) < 3:
            continue
        try:
            if float(f[2]) > qmax:
                continue
        except ValueError:
            continue
        accs = prot_of.get(f[0])
        if accs is None:
            continue
        # An ID is entrapment only if EVERY accession is an entrapment protein. A
        # peptide shared with a real sample protein is not evidence of an error.
        if all("ENTRAP_" in a for a in accs):
            n_entrap += 1
        else:
            n_sample += 1

total = n_entrap + n_sample
fdp = (n_entrap / total) if total else 0.0
print(f"q threshold      : {qmax:.4f}")
print(f"sample glycoPSMs : {n_sample}")
print(f"entrapment hits  : {n_entrap}")
print(f"entrapment FDP   : {fdp * 100:.2f}%")
if total == 0:
    print("verdict          : NO ACCEPTED IDs — nothing to judge")
elif fdp <= qmax:
    print(f"verdict          : CONSERVATIVE (FDP <= nominal {qmax * 100:.1f}%)")
else:
    print(f"verdict          : OPTIMISTIC — reported FDR understates the real error")
