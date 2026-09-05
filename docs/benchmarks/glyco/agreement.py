#!/usr/bin/env python3
"""Engine-vs-reference agreement at PSM and peptide level.

    agreement.py <truth.tsv[.gz]> <psms> [q]

Reports four things that are routinely conflated:

  glycoPSMs            accepted spectra (one PSM per spectrum after the top-1 collapse)
  glycopeptides        DISTINCT backbones behind those PSMs
  glycopeptidoforms    DISTINCT (backbone, glycan composition) pairs
  peptide agreement    share of the reference's backbones andes also finds ANYWHERE
  same-scan BACKBONE   of spectra BOTH identified, share calling the same backbone
  same-scan PEPTIDOFORM ... and the same GLYCAN COMPOSITION too

BACKBONE AGREEMENT IS NOT AGREEMENT. A glycopeptide identification is a backbone AND a
glycan; calling `HNSTGCLR` with two different compositions is two different answers. The
backbone number alone flatters both engines, and this codebase already has evidence that
glycan composition is the inconsistent part (mass largely right, composition not). Quote
the peptidoform number as the agreement figure; the backbone number only isolates WHERE a
disagreement comes from.
"""
import csv, gzip, re, sys


def open_maybe_gz(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


GLY_RE = re.compile(r"HexNAc(\d+)Hex(\d+)Fuc(\d+)NeuAc(\d+)NeuGc(\d+)")


def norm_glycan(text):
    """Canonical (HexNAc, Hex, Fuc, NeuAc, NeuGc) tuple, or None if absent/unparseable.
    Both the reference tables and the PIN peptide string carry this same layout, so one
    parser serves both sides."""
    m = GLY_RE.search(text or "")
    return tuple(int(g) for g in m.groups()) if m else None


def bare(pep):
    m = re.match(r"^[A-Za-z_]\.(.*)\.[A-Za-z_](\[[^\]]*\])?$", pep, re.DOTALL)
    core = m.group(1) if m else pep
    core = re.sub(r"\[[^\]]*\]", "", core)
    core = re.sub(r"[+-]\d+\.\d+", "", core)
    return re.sub(r"[^A-Z]", "", core.upper())


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    truth_p, psms_p = sys.argv[1], sys.argv[2]
    qcut = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01

    with open_maybe_gz(truth_p) as fh:
        trows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    truth = {(r["run"], int(r["scan"])): bare(r["peptide"]) for r in trows}
    truth_gly = {(r["run"], int(r["scan"])): norm_glycan(r.get("glycan", "")) for r in trows}
    runs = {k[0] for k in truth}
    # Same tail-token fallback as score_vs_truth.py, including its refusal to guess: two
    # runs sharing a tail (sample_R1, control_R1) would otherwise both bind to whichever
    # dict order happened to reach `next()` first, and the agreement numbers would be
    # silently wrong rather than absent.
    tails = {r: r.rsplit("_", 1)[-1] for r in runs}
    if len(set(tails.values())) != len(tails):
        dupes = [r for r in runs if list(tails.values()).count(tails[r]) > 1]
        sys.exit(f"run names are ambiguous after trimming: {sorted(dupes)}; rename the "
                 f"truth `run` column or pool with matching tags")

    def key_of(specid):
        m = re.search(r"scan=(\d+)", specid)
        if not m:
            return None
        run = next((r for r in runs if r in specid), None)
        if run is None:
            run = next((r for r, t in tails.items()
                        if re.search(rf"(^|[^A-Za-z0-9]){re.escape(t)}([^A-Za-z0-9]|$)", specid)), None)
        return (run, int(m.group(1))) if run else None

    andes, andes_gly = {}, {}
    with open(psms_p) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        qi, pi = hdr.index("q-value"), hdr.index("peptide")
        for r in rd:
            if len(r) <= pi:
                continue
            try:
                q = float(r[qi])
            except ValueError:
                continue
            if q > qcut:
                continue
            k = key_of(r[0])
            if k:
                andes[k] = bare(r[pi])
                andes_gly[k] = norm_glycan(r[pi])

    t_pep, a_pep = set(truth.values()), set(andes.values())
    both = set(truth) & set(andes)
    same = sum(1 for k in both if truth[k] == andes[k])
    # peptidoform agreement is only meaningful where BOTH sides carry a composition
    comparable = [k for k in both if truth_gly.get(k) and andes_gly.get(k)]
    same_form = sum(1 for k in comparable
                    if truth[k] == andes[k] and truth_gly[k] == andes_gly[k])
    t_form = {(truth[k], truth_gly[k]) for k in truth if truth_gly.get(k)}
    a_form = {(andes[k], andes_gly[k]) for k in andes if andes_gly.get(k)}

    print(f"reference          : {truth_p}")
    print(f"q cutoff           : {qcut}\n")
    print(f"{'':22} {'reference':>12} {'andes':>12}")
    print(f"{'glycoPSMs':22} {len(truth):12,} {len(andes):12,}")
    print(f"{'glycopeptides':22} {len(t_pep):12,} {len(a_pep):12,}")
    print(f"{'glycopeptidoforms':22} {len(t_form):12,} {len(a_form):12,}\n")
    print(f"  peptide agreement    {100.0*len(t_pep & a_pep)/max(len(t_pep),1):5.1f}%   "
          f"({len(t_pep & a_pep):,} of {len(t_pep):,} reference backbones also found by andes)")
    print(f"  co-identified scans  {len(both):,}   (spectra identified by BOTH)")
    print(f"  same-scan BACKBONE   {100.0*same/max(len(both),1):5.1f}%   "
          f"({same:,} of {len(both):,})")
    if comparable:
        print(f"  same-scan PEPTIDOFORM{100.0*same_form/len(comparable):5.1f}%   "
              f"({same_form:,} of {len(comparable):,} with a composition on both sides)"
              f"   <-- the agreement figure")
    else:
        print("  same-scan PEPTIDOFORM  n/a   (reference carries no glycan composition)")
    print(f"  scans only reference {len(set(truth)-set(andes)):,}")
    print(f"  scans only andes     {len(set(andes)-set(truth)):,}")


if __name__ == "__main__":
    main()
