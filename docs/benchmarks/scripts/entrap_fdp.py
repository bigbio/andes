#!/usr/bin/env python3
# True entrapment-FDP for andes Astral, independent of target-decoy q-values.
# DB = real HYE + 1:1 ENT_ paired entrapment (+ XXX_ decoys, already in d.psms).
# A target PSM whose peptide maps ONLY to ENT_ proteins is a known false positive.
# True FDP = 2 * E / T  (1:1 entrapment: each ENT hit implies ~1 undetected false
# hit among real proteins). q-values come from the XXX_ decoys, as published.
import os

RES = "/srv/data/msgf-bench/astral-entrapment"
PUBLISHED = {"top1": 36782, "chim": 69968}

def analyze(tag):
    path = f"{RES}/{tag}.t.psms"
    if not os.path.exists(path):
        print(f"### {tag}: NO OUTPUT (search/percolator failed)"); return
    rows = []  # (q, is_entrapment_unique)
    with open(path) as f:
        hdr = next(f).rstrip("\n").split("\t")
        qi = hdr.index("q-value")
        pi = hdr.index("proteinIds")
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= pi:
                continue
            try:
                q = float(p[qi])
            except ValueError:
                continue
            accs = [a for a in p[pi:] if a]           # proteinIds = rest-of-line
            is_ent = bool(accs) and all(a.startswith("ENT_") for a in accs)
            rows.append((q, is_ent))
    rows.sort(key=lambda r: r[0])
    print(f"### {tag}  (published @1%% = {PUBLISHED.get(tag,'?')}) ###")
    print(f"    {'q<=':>7} {'targets':>9} {'ENT':>6} {'realIDs':>8} {'trueFDP':>9}")
    for thr in (0.005, 0.01, 0.02, 0.05):
        T = sum(1 for q, e in rows if q <= thr)
        E = sum(1 for q, e in rows if q <= thr and e)
        fdp = (2.0 * E / T) if T else 0.0
        flag = "  <-- 1% line" if abs(thr - 0.01) < 1e-9 else ""
        print(f"    {thr:>7.3f} {T:>9d} {E:>6d} {T-E:>8d} {100*fdp:>8.2f}%{flag}")
    # Score-sorted deepest-prefix-with-FDP<=1% (the W1b/W1c metric) as a cross-check.
    print()

for tag in ("top1", "chim"):
    analyze(tag)
print("VERDICT: q-values are HONEST iff trueFDP tracks q (FDP ~ 1% at the q<=0.01 line).")
