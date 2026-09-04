#!/usr/bin/env python3
"""Corrected glyco gap decomposition (v2). Rewritten after a code review found four
defects in v1; each is now handled explicitly and asserted where possible.

v1 DEFECTS FIXED HERE
  D1 DUPLICATE TRUTH ROWS. Byonic's PQMs table emits one row per (PSM, protein), so a
     peptide mapping to N proteins appeared N times. v1 counted each row as a separate
     truth PSM, inflating totals (1068 rows -> 629 real spectra) and, worst, manufacturing
     a fake "one backbone = 33.8% of all misses" headline (213 rows = 11 spectra).
     FIX: dedupe to one truth entry per (fraction, scan).
  D2 WRONG DIGESTION RULE. v1 applied trypsin's "do not cleave before proline" rule.
     andes does NOT: crates/model/src/enzyme.rs has Trypsin = {after: "KR", before: ""},
     with no proline restriction. v1 therefore marked reachable peptides unreachable.
     FIX: cleave after K/R unconditionally, matching the engine under test.
  D3 WRONG LENGTH BOUNDS. v1 used 5..50; andes defaults are min_length=6, max_length=50.
  D4 STALE LABEL. v1 printed "MSFragger" while reading Byonic-derived truth.

Truth = the depositors' own Byonic results for the pure-HCD acquisitions, IsDecoy=0,
FalseDiscoveryRate1<=0.01, carrying an 'ngly' modification. Byonic's own digestion was
"KR, C-terminal cutter, <=2 missed cleavages" (read from its SearchParameters), so the
specificity comparison against andes (KR, <=3 MC for glyco) is fair.
"""
import csv, re, sys
from collections import defaultdict

B = "/srv/data/andes-bench/plasma"
FRACS = [1, 2, 3]
MIN_LEN, MAX_LEN, MAX_MC = 6, 50, 3   # andes defaults (glyco MC default = 3)

def bare(pep):
    """Strip PIN flanking residues, bracketed glycan tags and inline mod masses.
    NOTE: a naive split(".") is wrong -- inline masses ("C+57.02146") contain dots."""
    m = re.match(r"^[A-Za-z_]\.(.*)\.[A-Za-z_](\[[^\]]*\])?$", pep, re.DOTALL)
    core = m.group(1) if m else pep
    core = re.sub(r"\[[^\]]*\]", "", core)
    core = re.sub(r"[+-]\d+\.\d+", "", core)
    return re.sub(r"[^A-Z]", "", core.upper())

assert bare("R.HNSTGC+57.02146LR.M[HexNAc5Hex6@N2]") == "HNSTGCLR"
assert bare("_.MANSTGK.A[HexNAc4Hex5@N3]") == "MANSTGK"
assert bare("PREEQYNSTYR") == "PREEQYNSTYR"

def digest(fasta_path):
    """andes-equivalent trypsin: cleave after K/R with NO proline restriction."""
    seqs, cur = [], []
    for line in open(fasta_path):
        if line.startswith(">"):
            if cur: seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    peps = set()
    for seq in seqs:
        sites = [0] + [i + 1 for i in range(len(seq) - 1) if seq[i] in "KR"] + [len(seq)]
        sites = sorted(set(sites))
        for i in range(len(sites) - 1):
            for j in range(i + 1, min(i + 2 + MAX_MC, len(sites))):
                p = seq[sites[i]:sites[j]]
                if MIN_LEN <= len(p) <= MAX_LEN:
                    peps.add(p)
    return peps

print("digesting (trypsin KR, no proline rule, 6-50 aa, <=3 MC)...", file=sys.stderr)
tryptic = digest(f"{B}/human_entrap.fasta")
print(f"  {len(tryptic)} distinct peptides", file=sys.stderr)
# D2 regression check: these were mis-bucketed by v1's proline rule
for p in ("PREEQYNSTYR", "PREEQFNSTFR"):
    assert p in tryptic, f"{p} should be reachable once the proline rule is dropped"

# --- truth, deduped to one entry per (fraction, scan) [D1] ---
truth, dup_rows, conflicts = {}, 0, 0
for i in FRACS:
    with open(f"{B}/byonic_truth_frac{i}.tsv") as fh:
        rd = csv.reader(fh, delimiter="\t"); next(rd)
        for scan, pep, charge in rd:
            key = (i, int(scan)); pb = bare(pep)
            if key in truth:
                dup_rows += 1
                if truth[key] != pb: conflicts += 1
                continue
            truth[key] = pb
print(f"truth: {len(truth)} distinct (fraction, scan) entries "
      f"[{dup_rows} duplicate rows collapsed, {conflicts} with a conflicting peptide]")

# --- andes pre-Percolator emission: PIN is one row per scan (asserted) ---
emitted = {}
for i in FRACS:
    with open(f"{B}/gap_R{i}.glyco.pin") as fh:
        rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
        si, li, pi_ = hdr.index("ScanNr"), hdr.index("Label"), hdr.index("Peptide")
        for r in rd:
            key = (i, int(r[si]))
            assert key not in emitted, f"PIN not one-row-per-scan at {key}"
            emitted[key] = (bare(r[pi_]), r[li])

# --- Percolator acceptance, majority of 5 seeds ---
accepted = defaultdict(dict)
for s in range(1, 6):
    with open(f"{B}/gap_pooled_{s}.psms") as fh:
        rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
        qi, pi_ = hdr.index("q-value"), hdr.index("peptide")
        for r in rd:
            if len(r) <= pi_: continue
            mf = re.search(r"_R(\d)$", r[0]); ms = re.search(r"scan=(\d+)", r[0])
            if not (mf and ms): continue
            try: q = float(r[qi])
            except ValueError: continue
            accepted[(int(mf.group(1)), int(ms.group(1)))][s] = (bare(r[pi_]), q)

buckets = defaultdict(int)
missed_scans_by_pep = defaultdict(int)
for key, pep_b in truth.items():
    if pep_b not in tryptic:
        b = "NOT_DIGESTIBLE"
    elif key not in emitted or emitted[key][1] != "1":
        b = "NOT_EMITTED"
    elif emitted[key][0] != pep_b:
        b = "WRONG_WINNER"
    else:
        wins = sum(1 for (p, q) in accepted.get(key, {}).values() if p == pep_b and q <= 0.01)
        b = "CONFIRMED" if wins >= 3 else "FDR_REJECTED"
    buckets[b] += 1
    if b != "CONFIRMED":
        missed_scans_by_pep[pep_b] += 1

n = len(truth)
print(f"\ndecomposition over {n} truth spectra:")
for b in ("NOT_DIGESTIBLE", "NOT_EMITTED", "WRONG_WINNER", "FDR_REJECTED", "CONFIRMED"):
    print(f"  {b:16s} {buckets[b]:5d} ({100.0*buckets[b]/n:5.1f}%)")

total_missed = sum(missed_scans_by_pep.values())
top = sorted(missed_scans_by_pep.items(), key=lambda x: -x[1])[:8]
print(f"\ntop missed backbones (of {len(missed_scans_by_pep)}), by DISTINCT SPECTRA missed:")
for p, c in top:
    print(f"  {c:4d} ({100.0*c/total_missed:4.1f}% of {total_missed} missed)  {p[:50]}")

# --- cross-check: confirmed must not exceed andes's own pooled yield ---
pooled = [sum(1 for (p, q) in accepted[k].values() if q <= 0.01 and len(accepted[k]) >= 1)
          for k in accepted]
seed1 = sum(1 for k in accepted if 1 in accepted[k] and accepted[k][1][1] <= 0.01)
print(f"\ncross-check: andes accepted {seed1} PSMs @1% (seed 1, all scans); "
      f"CONFIRMED={buckets['CONFIRMED']} must be <= that: "
      f"{'OK' if buckets['CONFIRMED'] <= seed1 else 'INCONSISTENT'}")
