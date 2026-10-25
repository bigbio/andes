#!/usr/bin/env python3
"""Hex-containing O-glycan rate — the pGlyco3 Fig. 3b accuracy metric, for andes.

pGlyco3's O-glycan benchmark (IHMO HEK-293, Nat. Methods 2021, Fig. 3b) validates
O-glycopeptides with "Hex-suggested" glycan-level error rates: the IHMO inhibitor
(GalNAc-O-bn) truncates O-glycans to GalNAc (Tn) and sialyl-Tn, so any accepted
O-glycopeptide whose glycan still contains Hex (galactose) is a false positive unless
its Hex diagnostic ions (m/z 163.060 / 366.139) are present. pGlyco3 scored 2/484 =
~0.4%; MetaMorpheus 8/273 = 2.9%, Byonic 48/488 = 10%, MSFragger 64/320 = 20%.

This script computes the same denominator for andes: the fraction of accepted target
O-GPSMs (@1% q-value) whose glycan composition carries Hex > 0. It is an UPPER bound on
the error (it counts every Hex-containing ID as wrong, without the diagnostic-ion
rescue pGlyco3 applies), so andes only has to clear ~0.4% to match pGlyco3.

Usage:
    eval_hex_entrap.py <pooled.glyco.pin> <percolator.psms>
"""
import re
import sys

if len(sys.argv) < 3:
    sys.exit(__doc__)
pin, psms = sys.argv[1], sys.argv[2]

COMP = re.compile(r"HexNAc(\d+)Hex(\d+)Fuc(\d+)NeuAc(\d+)NeuGc(\d+)")


def comp_of(pep):
    m = COMP.search(pep)
    if not m:
        return None
    return tuple(int(x) for x in m.groups())  # hexnac, hex, fuc, neuac, neugc


# SpecId -> (label, glycan composition) for every PIN row.
rows = {}
hdr = None
with open(pin) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            i_pep, i_lab = hdr.index("Peptide"), hdr.index("Label")
            continue
        try:
            if int(f[i_lab]) != 1:
                continue
        except (ValueError, IndexError):
            continue
        c = comp_of(f[i_pep])
        if c:
            rows[f[0]] = c

total = hex_pos = 0
fuc_pos = neugc_pos = 0
comps = {}
QI = open(psms).readline().rstrip("\n").split("\t").index("q-value")  # Percolator column order varies (filename column)
with open(psms) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if f[0] == "PSMId" or len(f) < 3:
            continue
        try:
            if float(f[QI]) > 0.01:
                continue
        except ValueError:
            continue
        c = rows.get(f[0])
        if not c:
            continue
        total += 1
        hexnac, hex, fuc, neuac, neugc = c
        comps[c] = comps.get(c, 0) + 1
        if hex > 0:
            hex_pos += 1
        if fuc > 0:
            fuc_pos += 1
        if neugc > 0:
            neugc_pos += 1

print(f"  accepted O-GPSMs @1%      : {total}")
print(f"  Hex-containing (Hex>0)    : {hex_pos}  ({100 * hex_pos / total:.2f}%)" if total else
      f"  Hex-containing (Hex>0)    : {hex_pos}")
print(f"  Fuc-containing            : {fuc_pos}  ({100 * fuc_pos / total:.2f}%)" if total else
      f"  Fuc-containing            : {fuc_pos}")
print(f"  NeuGc-containing          : {neugc_pos}  ({100 * neugc_pos / total:.2f}%)" if total else
      f"  NeuGc-containing          : {neugc_pos}")
print(f"  pGlyco3 Hex error (Fig 3b): 2/484 = 0.41%   (MetaMorpheus 2.9%, Byonic 10%, MSFragger 20%)")
print("\n  top accepted compositions:")
for c, n in sorted(comps.items(), key=lambda kv: -kv[1])[:12]:
    print(f"    HexNAc{c[0]}Hex{c[1]}Fuc{c[2]}NeuAc{c[3]}NeuGc{c[4]}  x{n}")
