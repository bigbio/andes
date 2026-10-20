#!/usr/bin/env python3
"""O-glycan overlap: pGlyco3 vs andes on the IHMO HEK293 (CWQ-O-293) benchmark.

    overlap_pglyco3_o.py <pglyco3-293.txt> <andes.psms> [q]

pGlyco3's published result for CWQ-O-293 (`pGlyco3-293.txt`, 484 GPSMs, already
target-only and TotalFDR <= 0.01) is normalised to the harness's canonical form and
compared against andes's accepted PSMs at `q` (default 1%).

The two sides name the same spectra the same way (`CWQ-O-293-N` + raw scan number), so
the comparison is done BOTH as a set overlap (which peptidoforms are shared) and as a
same-scan agreement (of the spectra both engines identified, how often they call the
same thing).

Two pGlyco3 specifics that a naive diff gets wrong:

  * glycan is `Glycan(H,N,F,A)` — a 4-vector in the order Hex, HexNAc, Fuc, NeuAc
    (no NeuGc: human). Mapped to the canonical (HexNAc, Hex, Fuc, NeuAc, NeuGc).
  * `J` in the peptide is pGlyco's placeholder for the Asn (N) of an N-X-S/T sequon,
    a fixed peptide annotation it applies in O-glycan mode too — it is NOT the
    O-glycosylation site (that is `GlySite`/`LocalizedSiteGroups`, a Ser/Thr). So the
    backbone is recovered with `J -> N`, the same convention as `make_truth.py engine_l`.
    (Verified against UniProt: every `...J...` peptide is a real `...N...` in the FASTA.)
"""
import csv, re, sys

if len(sys.argv) < 3:
    sys.exit(__doc__)
pg_path, andes_path = sys.argv[1], sys.argv[2]
qcut = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01

GLY_PG = re.compile(r"HexNAc(\d+)Hex(\d+)Fuc(\d+)NeuAc(\d+)NeuGc(\d+)")


def andes_bare(pep):
    """PIN peptide `K.SEQ[+glycan].K[@site]` -> bare uppercase backbone."""
    m = re.match(r"^[A-Za-z_]\.(.*)\.[A-Za-z_](\[[^\]]*\])?$", pep, re.DOTALL)
    core = m.group(1) if m else pep
    core = re.sub(r"\[[^\]]*\]", "", core)
    core = re.sub(r"[+-]\d+\.\d+", "", core)
    return re.sub(r"[^A-Z]", "", core.upper())


def andes_glycan(pep):
    m = GLY_PG.search(pep)
    return tuple(int(g) for g in m.groups()) if m else None


def pg_backbone(pep):
    """pGlyco3 peptide -> bare backbone; `J` (the N-X-S/T sequon Asn) becomes `N`."""
    return re.sub(r"[^A-Z]", "", pep.upper().replace("J", "N"))


def pg_glycan(vec):
    """`Glycan(H,N,F,A)` `H N F A` -> (HexNAc, Hex, Fuc, NeuAc, NeuGc)."""
    p = vec.split()
    if len(p) != 4:
        return None
    try:
        h, n, f, a = (int(float(x)) for x in p)
    except ValueError:
        return None
    return (n, h, f, a, 0)


def scan_key_pg(r):
    return (r["RawName"], int(r["Scan"]))


def scan_key_andes(specid):
    run = specid.split("_controllerType")[0]
    m = re.search(r"scan=(\d+)", specid)
    return (run, int(m.group(1))) if m else None


# ---- pGlyco3 reference ----
pg_rows = []
for r in csv.DictReader(open(pg_path), delimiter="\t"):
    gly = pg_glycan(r["Glycan(H,N,F,A)"])
    if gly is None:
        continue
    pg_rows.append((scan_key_pg(r), pg_backbone(r["Peptide"]), gly))

pg_scans = {k for k, _, _ in pg_rows}
pg_pep = {p for _, p, _ in pg_rows}                      # distinct backbone patterns
pg_form = {(p, g) for _, p, g in pg_rows}                # distinct peptidoforms

# ---- andes accepted ----
andes_scan, andes_pep, andes_form = {}, set(), set()
for r in csv.DictReader(open(andes_path), delimiter="\t"):
    try:
        q = float(r["q-value"])
    except (KeyError, ValueError):
        continue
    if q > qcut:
        continue
    k = scan_key_andes(r["PSMId"])
    if k is None:
        continue
    b = andes_bare(r["peptide"])
    g = andes_glycan(r["peptide"])
    andes_scan[k] = (b, g)
    andes_pep.add(b)
    if g is not None:
        andes_form.add((b, g))

# ---- agreement ----
# set-level: a pGlyco3 peptidoform is "shared" if andes emitted the same backbone+glycan
# anywhere (not necessarily on the same scan). Matching is exact equality now that both
# sides are canonical, so this is plain set intersection.
shared_pep = pg_pep & andes_pep
shared_form = pg_form & andes_form

# same-scan: spectra identified by BOTH engines.
both = set(pg_scans) & set(andes_scan)

# per-scan pGlyco3 backbone map (rank 1 only, file is already rank 1)
pg_by_scan = {}
for k, p, g in pg_rows:
    pg_by_scan[k] = (p, g)

same_backbone = 0
same_form = 0
for k in sorted(both):
    p, g = pg_by_scan[k]
    b, ag = andes_scan[k]
    if p == b:
        same_backbone += 1
        if g == ag:
            same_form += 1

n_pep = len(pg_pep)
n_form = len(pg_form)
print(f"reference (pGlyco3) : {pg_path}")
print(f"andes accepted       : {andes_path} @ q<={qcut}\n")
print(f"{'':24} {'pGlyco3':>10} {'andes':>10}")
print(f"{'glycoPSMs':24} {len(pg_rows):10,} {len(andes_scan):10,}")
print(f"{'glycopeptides (backbone)':24} {n_pep:10,} {len(andes_pep):10,}")
print(f"{'glycopeptidoforms':24} {n_form:10,} {len(andes_form):10,}\n")
print(f"  backbone agreement    {100.0*len(shared_pep)/max(n_pep,1):5.1f}%   "
      f"({len(shared_pep):,} of {n_pep:,} pGlyco3 backbones also found by andes)")
print(f"  peptidoform agreement {100.0*len(shared_form)/max(n_form,1):5.1f}%   "
      f"({len(shared_form):,} of {n_form:,} pGlyco3 peptidoforms also found by andes)")
print(f"  co-identified scans   {len(both):,}   (spectra identified by BOTH)")
print(f"  same-scan BACKBONE    {100.0*same_backbone/max(len(both),1):5.1f}%   "
      f"({same_backbone:,} of {len(both):,})")
print(f"  same-scan PEPTIDOFORM {100.0*same_form/max(len(both),1):5.1f}%   "
      f"({same_form:,} of {len(both):,})")
print(f"  scans only pGlyco3    {len(pg_scans - set(andes_scan)):,}")
print(f"  scans only andes      {len(set(andes_scan) - pg_scans):,}")

# write the actual peptidoform lists next to the andes results
import os
outdir = os.path.dirname(os.path.abspath(andes_path))


def fmt_gly(g):
    return f"HexNAc{g[0]}Hex{g[1]}Fuc{g[2]}NeuAc{g[3]}NeuGc{g[4]}"


def write(name, forms):
    p = os.path.join(outdir, name)
    with open(p, "w") as fh:
        fh.write("backbone\tglycan(HexNAc,Hex,Fuc,NeuAc,NeuGc)\n")
        for bb, g in sorted(forms):
            fh.write(f"{bb}\t{fmt_gly(g)}\n")
    return p


print(f"\n  written: {write('overlap_shared.tsv', pg_form & andes_form)}")
print(f"  written: {write('overlap_pglyco3_only.tsv', pg_form - andes_form)}")
print(f"  written: {write('overlap_andes_only.tsv', andes_form - pg_form)}")
