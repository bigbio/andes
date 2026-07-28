#!/usr/bin/env python3
"""Build an entrapment FASTA: mouse (target+decoy, unchanged) + yeast/E.coli entrapment
(target+decoy). Any glycopeptide identified on an ENTRAP_ protein at 1% FDR is a
false positive, giving an absolute error estimate that does not depend on the reference
identification set. Yeast/E.coli are near-orthogonal to mouse brain, and prokaryotes lack the
mammalian N-glycosylation machinery, so a glyco ID there is unambiguously wrong."""
import re
import os
import sys

# Usage: build_entrap.py <target.fasta> <entrapment.fasta> <out.fasta>
# Defaults kept only so the script still runs bare on the benchmark host.
MOUSE = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ENTRAP_TARGET", "target-decoy.fasta")
HYE = sys.argv[2] if len(sys.argv) > 2 else os.environ.get("ENTRAP_SOURCE", "entrapment.fasta")
OUT = sys.argv[3] if len(sys.argv) > 3 else os.environ.get("ENTRAP_OUT", "target-entrap.fasta")

def read_fasta(p):
    acc, seq = None, []
    for ln in open(p):
        if ln.startswith(">"):
            if acc: yield acc, "".join(seq)
            acc, seq = ln.rstrip("\n"), []
        else:
            seq.append(ln.strip())
    if acc: yield acc, "".join(seq)

n_m = n_e = 0
with open(OUT, "w") as out:
    for acc, seq in read_fasta(MOUSE):          # mouse target+decoy verbatim
        out.write(f"{acc}\n{seq}\n"); n_m += 1
    for acc, seq in read_fasta(HYE):
        if "Cont_" in acc: continue
        if not re.search(r"_(YEAST|ECOLI)\b", acc): continue
        ident = acc[1:].split()[0]
        out.write(f">ENTRAP_{ident}\n{seq}\n")
        out.write(f">DECOY_ENTRAP_{ident}\n{seq[::-1]}\n")
        n_e += 1
print(f"mouse entries (target+decoy): {n_m}")
print(f"entrapment proteins: {n_e}  (+{n_e} reversed decoys)")
