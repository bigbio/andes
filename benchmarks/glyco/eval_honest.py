#!/usr/bin/env python3
"""Honest evaluator for the andes glyco benchmark.

Fixes three defects in eth_bench_eval.py:
  1. DENOMINATOR: it divided by all-6-fraction truth (5088) even when only a
     subset of fractions was searched. We derive the searched fractions from the
     PIN's `f{N}_` SpecId prefixes and score against THAT subset.
  2. CORRECTNESS: it declared a hit "backbone-correct" on a 0.1 Da MASS match
     (~40 ppm) and never compared the peptide SEQUENCE, which both files carry.
     We do I/L-equivalent sequence matching and report the mass-only overcount.
  3. UNCOUNTED IDs: it silently ignored every PSM that did not land on a truth
     scan. We report them, since the reference (Byonic, no glycan FDR) is known
     to miss a large fraction of what other engines find on these files.

Usage: eval_honest.py <pooled.pin> <psms> [dpsms]
"""
import sys, os, re
from collections import defaultdict

H2O = 18.010565
pin, psms = sys.argv[1], sys.argv[2]
dpsms = sys.argv[3] if len(sys.argv) > 3 else None

def bare(seq):
    """PIN peptide 'X.SEQ.Y[glycan]' -> bare uppercase sequence, I->L."""
    s = re.sub(r"\[[^\]]*\]$", "", seq.strip())   # drop trailing glycan tag
    if len(s) > 4 and s[1] == "." and s[-2] == ".":
        s = s[2:-2]                                # drop single-char flanks
    return re.sub(r"[^A-Za-z]", "", s).upper().replace("I", "L")

# ---- truth, per fraction -------------------------------------------------
truth = {}
per_frac = defaultdict(int)
for fr in range(1, 7):
    p = f"ethcd_truth_frac{fr}.tsv"
    if not os.path.exists(p):
        continue
    for ln in open(p):
        if ln.startswith("scan"):
            continue
        c = ln.rstrip("\n").split("\t")
        truth[(fr, int(c[0]))] = (float(c[4]), int(c[2]), bare(c[1]))
        per_frac[fr] += 1

# ---- PIN ------------------------------------------------------------------
hdr = None
info = {}
searched = set()
for ln in open(pin):
    f = ln.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f
        iS, iC, iG = hdr.index("ScanNr"), hdr.index("CalcMass"), hdr.index("GlycanMass")
        iP, iL = hdr.index("Peptide"), hdr.index("Label")
        continue
    sid = f[0]
    m = re.match(r"f(\d+)_", sid)
    if not m:
        continue
    fr = int(m.group(1))
    searched.add(fr)
    try:
        info[sid] = (fr, int(float(f[iS])), float(f[iC]), float(f[iG]), bare(f[iP]), int(f[iL]))
    except (ValueError, IndexError):
        continue

denom = sum(per_frac[fr] for fr in searched)

# ---- score PSMs @1% -------------------------------------------------------
seq_ok, mass_only, wrong_on_truth, off_truth = set(), set(), set(), 0
ntgt = 0
by_charge_ok = defaultdict(set)
for ln in open(psms):
    f = ln.rstrip("\n").split("\t")
    if f[0] == "PSMId" or len(f) < 3:
        continue
    try:
        if float(f[2]) > 0.01:
            continue
    except ValueError:
        continue
    ntgt += 1
    rec = info.get(f[0])
    if rec is None:
        continue
    fr, scan, calc, gly, pep, label = rec
    k = (fr, scan)
    t = truth.get(k)
    if t is None:
        off_truth += 1
        continue
    tb, ch, tpep = t
    massv = abs(calc - gly - H2O - tb) < 0.1 or abs(calc - gly - tb) < 0.1
    if pep and pep == tpep:
        seq_ok.add(k); by_charge_ok[ch].add(k)
    elif massv:
        mass_only.add(k)
    else:
        wrong_on_truth.add(k)

ndec = 0
if dpsms and os.path.exists(dpsms):
    for ln in open(dpsms):
        f = ln.rstrip("\n").split("\t")
        if f[0] == "PSMId" or len(f) < 3:
            continue
        try:
            if float(f[2]) <= 0.01:
                ndec += 1
        except ValueError:
            pass

tt = defaultdict(int)
for (fr, s), (tb, ch, tp) in truth.items():
    if fr in searched:
        tt[ch] += 1

legacy = len(seq_ok) + len(mass_only)
print(f"  fractions searched : {sorted(searched)}   truth in those = {denom}  (all-6 = {len(truth)})")
print(f"  targets @1%        : {ntgt}   decoys @1% : {ndec}")
print(f"  SEQUENCE-correct   : {len(seq_ok)}/{denom} = {100*len(seq_ok)/max(denom,1):.1f}%")
print(f"  mass-only (WRONG peptide, right mass) : {len(mass_only)}  <- legacy metric counted these as correct")
print(f"  wrong on truth scan: {len(wrong_on_truth)}")
print(f"  OFF-truth PSMs     : {off_truth}  ({100*off_truth/max(ntgt,1):.0f}% of targets; uncredited by the reference)")
print(f"  legacy-equivalent  : {legacy}/{denom} = {100*legacy/max(denom,1):.1f}%   (legacy printed /{len(truth)} = {100*legacy/max(len(truth),1):.0f}%)")
for ch in sorted(tt):
    if tt[ch] >= 3:
        print(f"    z{ch}: {len(by_charge_ok[ch])}/{tt[ch]} = {100*len(by_charge_ok[ch])/tt[ch]:.0f}%")

# ---- WHERE DO THE UNANSWERED TRUTH SCANS DIE? -----------------------------
# For every truth scan in a searched fraction, classify:
#   A emitted-correct & passed 1%      -> already won
#   B emitted-correct but FAILED 1%    -> FDR/separation-limited
#   C emitted, but the WRONG peptide   -> selection-limited
#   D no PIN row at all                -> generation/gating-limited
pass_ids = set()
for ln in open(psms):
    f = ln.rstrip("\n").split("\t")
    if f[0] == "PSMId" or len(f) < 3:
        continue
    try:
        if float(f[2]) <= 0.01:
            pass_ids.add(f[0])
    except ValueError:
        pass
best = {}
for sid, (fr, scan, calc, gly, pep, label) in info.items():
    if label != 1:
        continue
    k = (fr, scan)
    if k in truth:
        best.setdefault(k, []).append((pep, sid in pass_ids))
A = B = C = D = 0
for k, (tb, ch, tpep) in truth.items():
    if k[0] not in searched:
        continue
    rows = best.get(k)
    if not rows:
        D += 1
    elif any(p == tpep and ok for p, ok in rows):
        A += 1
    elif any(p == tpep for p, ok in rows):
        B += 1
    else:
        C += 1
tot = A + B + C + D
print(f"\n  --- where the {denom} truth scans end up ---")
print(f"  A won (correct, @1%)          : {A:5d}  {100*A/tot:.1f}%")
print(f"  B correct but BELOW 1% FDR    : {B:5d}  {100*B/tot:.1f}%   <- FDR/separation-limited")
print(f"  C emitted the WRONG peptide   : {C:5d}  {100*C/tot:.1f}%   <- selection-limited")
print(f"  D no PIN row at all           : {D:5d}  {100*D/tot:.1f}%   <- generation/gating-limited")
