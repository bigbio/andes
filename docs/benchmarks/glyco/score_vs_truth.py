#!/usr/bin/env python3
"""Score a glyco run against any reference table in `truth/`, and say WHERE it loses.

    score_vs_truth.py <truth.tsv[.gz]> <pooled.pin> <psms> [q]

Replaces the plasma-only `gap_decompose.py`, which hardcoded one dataset's paths and
fraction numbering and so could not run on anything else, or on any machine but the
benchmark host.

Every miss is attributed to exactly one stage, in this order — the order matters, because
these are different problems with different fixes:

  NOT_EMITTED       no PIN row for that spectrum at all -> generation/retrieval
  DECOY_WON         a decoy beat the true peptide       -> scoring cannot separate
  WRONG_TARGET_WON  a different target won the collapse -> selection
  FDR_REJECTED      right answer emitted, lost at q     -> separation/statistics
  CONFIRMED         right answer emitted and accepted

DECOY_WON is reported separately from NOT_EMITTED on purpose: lumping them hid the single
largest cause on plasma (15.6% of all truth) behind a label that implied a generation
problem.

Matching is on the BARE backbone, because glycan-composition disagreement is a separate,
known defect and folding it in here would conflate two independent failures.
"""
import csv, gzip, re, sys
from collections import defaultdict


def open_maybe_gz(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def bare(pep):
    """PIN peptides look like `K.SEQ[+glycan].K[tag]` with INLINE mod masses that contain
    dots, so splitting on '.' truncates them. Anchor on the flanking-residue structure."""
    m = re.match(r"^[A-Za-z_]\.(.*)\.[A-Za-z_](\[[^\]]*\])?$", pep, re.DOTALL)
    core = m.group(1) if m else pep
    core = re.sub(r"\[[^\]]*\]", "", core)
    core = re.sub(r"[+-]\d+\.\d+", "", core)
    return re.sub(r"[^A-Z]", "", core.upper())


assert bare("R.HNSTGC+57.02146LR.M[HexNAc5Hex6@N2]") == "HNSTGCLR"


def load_truth(path):
    with open_maybe_gz(path) as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    return {(r["run"], int(r["scan"])): bare(r["peptide"]) for r in rows}


def run_matchers(runs):
    """Map each truth run name to the token that identifies it inside a PIN SpecId.

    The two sides are named independently: a reference table names runs after the source
    file (`plasma_sceHCD_R1`), while a pooled PIN tags SpecIds however `pool_pins.py`
    tagged them (`...scan=2302_glyco_2302_1_R1`). Requiring the full name to appear
    matched NOTHING and silently reported 100% NOT_EMITTED. So try the full name first,
    then fall back to the run's last underscore-token, and refuse to guess if that
    fallback is ambiguous."""
    full = {r: r for r in runs}
    tails = {r: r.rsplit("_", 1)[-1] for r in runs}
    if len(set(tails.values())) != len(tails):
        dupes = [r for r in runs if list(tails.values()).count(tails[r]) > 1]
        sys.exit(f"run names are ambiguous after trimming: {sorted(dupes)}; rename the "
                 f"truth `run` column or pool with matching tags")
    return full, tails


def key_of(specid, full, tails):
    """PIN SpecId carries `scan=N`; identify the run by full name, else by its tail token."""
    m = re.search(r"scan=(\d+)", specid)
    if not m:
        return None
    # Full-name match on token boundaries: a run called `R1` must not claim `R10`.
    hits = [r for r, v in full.items() if re.search(rf"(^|[^A-Za-z0-9]){re.escape(v)}([^A-Za-z0-9]|$)", specid)]
    if len(hits) > 1:
        sys.exit(f"SpecId {specid!r} matches several runs by full name: {sorted(hits)}")
    run = hits[0] if hits else None
    if run is None:
        run = next((r for r, t in tails.items() if re.search(rf"(^|[^A-Za-z0-9]){re.escape(t)}([^A-Za-z0-9]|$)", specid)), None)
    return (run, int(m.group(1))) if run else None


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    truth_p, pin_p, psms_p = sys.argv[1:4]
    qcut = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
    truth = load_truth(truth_p)
    full, tails = run_matchers({k[0] for k in truth})

    emitted = {}
    with open(pin_p) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        si, li, pi = hdr.index("SpecId"), hdr.index("Label"), hdr.index("Peptide")
        for r in rd:
            k = key_of(r[si], full, tails)
            if k:
                emitted[k] = (bare(r[pi]), r[li])

    accepted = {}
    with open(psms_p) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        qi, pi = hdr.index("q-value"), hdr.index("peptide")
        for r in rd:
            if len(r) <= pi:
                continue
            k = key_of(r[0], full, tails)
            try:
                q = float(r[qi])
            except ValueError:
                continue
            if k and q <= qcut:
                accepted[k] = bare(r[pi])

    buckets = defaultdict(int)
    for k, pep in truth.items():
        if k not in emitted:
            b = "NOT_EMITTED"
        elif emitted[k][1] != "1":
            b = "DECOY_WON"
        elif emitted[k][0] != pep:
            b = "WRONG_TARGET_WON"
        elif accepted.get(k) == pep:
            b = "CONFIRMED"
        else:
            b = "FDR_REJECTED"
        buckets[b] += 1

    n = len(truth)
    print(f"reference : {truth_p}")
    print(f"spectra   : {n}   (runs: {len(full)})")
    print(f"q cutoff  : {qcut}\n")
    for b in ("CONFIRMED", "WRONG_TARGET_WON", "DECOY_WON", "FDR_REJECTED", "NOT_EMITTED"):
        print(f"  {b:18s} {buckets[b]:6d}  ({100.0*buckets[b]/n:5.1f}%)")
    sel = buckets["WRONG_TARGET_WON"] + buckets["DECOY_WON"]
    print(f"\n  selection loss (wrong target + decoy): {sel} ({100.0*sel/n:.1f}%)")
    own = len(accepted)
    print(f"  andes accepted {own} spectra @q<={qcut}; "
          f"{own - buckets['CONFIRMED']} of those are NOT in this reference")


if __name__ == "__main__":
    main()
