#!/usr/bin/env python3
"""ABSOLUTE yield metrics - what andes actually finds, independent of the reference.

The benchmark scores agreement with the reference identification set. That
structurally cannot reward a real
glycopeptide the reference never reported, and entrapment puts andes's true error at
~0.45%, so the uncredited IDs are overwhelmingly genuine. These are the numbers that
answer "more IDs / more glycans":
   glycoPSMs        - PSMs @1%
   glycopeptides    - unique (peptide, glycan composition)
   glycan comps     - unique compositions observed
   glycosites       - unique (protein, peptide) carrying a glycan
Usage: eval_yield.py <pooled.pin> <psms>
"""
import sys, re
pin, psms = sys.argv[1], sys.argv[2]

def bare(s):
    s = re.sub(r"\[[^\]]*\]$", "", s.strip())
    if len(s) > 4 and s[1] == "." and s[-2] == ".": s = s[2:-2]
    return re.sub(r"[^A-Za-z]", "", s).upper()

hdr = None; rows = {}
for ln in open(pin):
    f = ln.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f; iP = hdr.index("Peptide"); iL = hdr.index("Label")
        iPr = len(hdr) - 1
        continue
    try:
        if int(f[iL]) != 1: continue
        pep_raw = f[iP]
        m = re.search(r"\[([^\]]*)\]$", pep_raw.strip())
        rows[f[0]] = (bare(pep_raw), m.group(1) if m else "", f[iPr] if len(f) > iPr else "")
    except (ValueError, IndexError):
        continue

npsm = 0
gp, comps, sites = set(), set(), set()
for ln in open(psms):
    f = ln.rstrip("\n").split("\t")
    if f[0] == "PSMId" or len(f) < 3: continue
    try:
        if float(f[2]) > 0.01: continue
    except ValueError: continue
    r = rows.get(f[0])
    if not r: continue
    pep, comp, prot = r
    npsm += 1
    gp.add((pep, comp)); comps.add(comp)
    sites.add((prot.split(";")[0] if prot else "", pep))
print(f"  glycoPSMs @1%        : {npsm}")
print(f"  unique glycopeptides : {len(gp)}")
print(f"  unique glycan comps  : {len(comps)}")
print(f"  unique glycosites    : {len(sites)}")
