#!/usr/bin/env python3
"""Phase A validation: best-PSM-per-scan collapse + peptide-level entrapment-FDP.

Given a UNIFIED refine PIN's Percolator target.psms.txt (KEEPS-BOTH mode: many
rows per scan), collapse each scan to its single highest-`score` PSM (the unified
per-scan winner across unmodified ∪ modified), then:
  - count modified winners accepted @ q<=qthr (peptide carries a +/-delta);
  - compute the TRUE peptide-level entrapment-FDP on the accepted winners
    (backbone tryptic-in-real vs only-in-ENT_, 1:1 DB → combined = 2×fraction);
  - report the double-count the collapse removes (scans whose pre-collapse rows
    contained BOTH a modified and an unmodified PSM at q<=qthr).

Usage: phaseA_collapse_fdp.py <entrap.fasta> <target.psms.txt> [qthr=0.01]
"""
import sys, re, collections

fasta, psms = sys.argv[1], sys.argv[2]
qthr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01


def read_fasta(p):
    n, s = None, []
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if n is not None:
                yield n, "".join(s)
            n, s = l[1:].split()[0], []
        else:
            s.append(l.strip())
    if n is not None:
        yield n, "".join(s)


def tryptic(seq, mc=2):
    sites = [0] + [i + 1 for i in range(len(seq) - 1) if seq[i] in "KR" and seq[i + 1] != "P"] + [len(seq)]
    sites = sorted(set(sites))
    out = set()
    for a in range(len(sites) - 1):
        for b in range(a + 1, min(a + 2 + mc, len(sites))):
            x = seq[sites[a]:sites[b]]
            if 6 <= len(x) <= 50:
                out.add(x)
    return out


real, ent = set(), set()
for n, s in read_fasta(fasta):
    (ent if n.startswith("ENT_") else real).__ior__(tryptic(s))


def backbone(pep):
    pep = re.sub(r"^[A-Z_\-]\.", "", pep)
    pep = re.sub(r"\.[A-Z_\-]$", "", pep)
    return re.sub(r"[^A-Z]", "", pep)


def scan_of(psm_id):
    # SpecId = {spec}_{scan}_{rank}_{rowidx}; use the full prefix up to the last
    # two underscore fields as the scan key (robust to spec ids containing '_').
    parts = psm_id.rsplit("_", 2)
    return parts[0] if len(parts) == 3 else psm_id


def is_modified(pep):
    return re.search(r"[+\-]\d+\.\d+", pep) is not None


hdr = open(psms).readline().rstrip("\n").split("\t")
si, qi, pepi = hdr.index("score"), hdr.index("q-value"), hdr.index("peptide")

# Gather rows per scan; track q<=thr mod/unmod presence for the double-count stat.
best = {}                       # scan -> (score, q, peptide)
had_mod_q = collections.Counter()
had_unmod_q = collections.Counter()
rows_total = 0
with open(psms) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) <= pepi:
            continue
        try:
            sc, q = float(f[si]), float(f[qi])
        except ValueError:
            continue
        rows_total += 1
        scan = scan_of(f[0])
        pep = f[pepi]
        if q <= qthr:
            if is_modified(pep):
                had_mod_q[scan] += 1
            else:
                had_unmod_q[scan] += 1
        cur = best.get(scan)
        if cur is None or sc > cur[0]:
            best[scan] = (sc, q, pep)

# Best-per-scan winners accepted at q<=thr.
acc = [(s, p) for s, (sc, q, p) in best.items() if q <= qthr]
mod_winners = [(s, p) for s, p in acc if is_modified(p)]

# Entrapment-FDP on the accepted best-per-scan winners (all, and mod-only).
def fdp(rows):
    n = real_hits = ent_hits = amb = 0
    for _s, p in rows:
        b = backbone(p)
        n += 1
        if b in real:
            real_hits += 1
        elif b in ent:
            ent_hits += 1
        else:
            amb += 1
    frac = ent_hits / n if n else 0.0
    return n, real_hits, ent_hits, amb, frac

# Double-count removed: scans that had BOTH a mod and an unmod PSM at q<=thr
# pre-collapse (the disjoint-union would count such a scan twice).
both = sum(1 for s in had_mod_q if had_unmod_q.get(s, 0) > 0)

print(f"rows={rows_total}  distinct_scans={len(best)}")
print(f"best-per-scan winners accepted @q<={qthr}: {len(acc)}  (of which MODIFIED: {len(mod_winners)})")
n, rh, eh, amb, frac = fdp(acc)
print(f"[ALL winners]  n={n} real={rh} ENT={eh} amb={amb}  ENT_frac={frac:.4f}  combined-FDP(1:1)={2*frac:.4f}")
n, rh, eh, amb, frac = fdp(mod_winners)
print(f"[MOD winners]  n={n} real={rh} ENT={eh} amb={amb}  ENT_frac={frac:.4f}  combined-FDP(1:1)={2*frac:.4f}")
print(f"double-count REMOVED by collapse (scans with both mod & unmod @q<={qthr} pre-collapse): {both}")
