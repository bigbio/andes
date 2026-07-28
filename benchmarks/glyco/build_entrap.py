#!/usr/bin/env python3
"""Build an entrapment FASTA: mouse (target+decoy, unchanged) + yeast/E.coli entrapment
(target+decoy). Any glycopeptide identified on an ENTRAP_ protein at 1% FDR is a
false positive, giving an absolute error estimate that does not depend on the Byonic
reference. Yeast/E.coli are near-orthogonal to mouse brain, and prokaryotes lack the
mammalian N-glycosylation machinery, so a glyco ID there is unambiguously wrong."""
import re
HYE = "/srv/data/msgf-bench/astral-data/ProteoBenchFASTA_MixedSpecies_HYE.fasta"
MOUSE = "/srv/data/msgf-bench/ethcd/mouse-decoy.fasta"
OUT = "/srv/data/msgf-bench/ethcd/mouse-entrap.fasta"

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
