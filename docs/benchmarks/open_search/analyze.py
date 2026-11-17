#!/usr/bin/env python3
"""Score the open-search benchmark produced by run_open_bench.sh.

Ground truth is the narrow (default-mods) search's oxidised-M PSMs. A truth PSM
is *recovered* when the open search reports the same scan and the same bare
peptide with an offset within 0.02 Da of +15.9949. The site is compared when
the open row localised the offset (OpenSite != 0).

The same script reports the entrapment FDP of both runs at q <= 0.01, using the
shared rule of ../entrap_common.py, so the nominal q can be checked against a
measured one.

    python3 -I analyze.py --dir results-open --fasta mouse_entrap.fasta --out report.json
"""
import argparse
import collections
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import entrap_common  # noqa: E402

RESIDUES = set("ACDEFGHIKLMNPQRSTVWYUO")
OXIDATION = 15.994915
OX_TOL_DA = 0.02
# andes' `model::mass::ISOTOPE` (C13 - C), used to turn a reported offset in the
# PIN's `isotope_error` into a neutral-mass correction.
ISOTOPE_DA = 1.00335483507


def parse_peptide(cell):
    """`pre.SEQ{+delta}.post` -> (bare sequence, [(0-based pos, delta), ...]).

    The flanking residues are single characters, so strip them positionally: a
    split on `.` would also cut a mod delta such as `+15.99491` in two.
    """
    middle = cell[2:-2] if len(cell) > 4 and cell[1] == "." and cell[-2] == "." else cell
    seq = []
    mods = []
    i = 0
    while i < len(middle):
        c = middle[i]
        if c not in RESIDUES:
            i += 1
            continue
        seq.append(c)
        i += 1
        while i < len(middle) and middle[i] in "+-0123456789.":
            j = i + 1 if middle[i] in "+-" else i
            while j < len(middle) and (middle[j].isdigit() or middle[j] == "."):
                j += 1
            if j == i + 1 and middle[i] in "+-":
                break
            mods.append((len(seq) - 1, float(middle[i:j])))
            i = j
    return "".join(seq), mods


def read_pin(path):
    """SpecId -> list of rows (the PIN holds the top-N rows per scan)."""
    rows = collections.defaultdict(list)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            rows[f[idx["SpecId"]]].append(
                {
                    "peptide": f[idx["Peptide"]],
                    "calc": float(f[idx["CalcMass"]]),
                    "exp": float(f[idx["ExpMass"]]),
                    "iso": int(f[idx["isotope_error"]]),
                    "mass": float(f[idx["DeltaMass"]]) if "DeltaMass" in idx else None,
                    "open_site": int(f[idx["OpenSite"]]) if "OpenSite" in idx else 0,
                }
            )
    return rows


def read_percolator(path):
    """PSMId -> (q-value, peptide, protein cells)."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_q = header.index("q-value")
        i_p = header.index("peptide")
        i_prot = header.index("proteinIds")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_p:
                continue
            try:
                q = float(f[i_q])
            except ValueError:
                continue
            out[f[0]] = (q, f[i_p], f[i_prot:])
    return out


def pin_delta(row):
    """Modification mass of one PIN row, preferring the reported DeltaMass."""
    if row.get("mass") is not None:
        return row["mass"]
    return row["exp"] - row["iso"] * ISOTOPE_DA - row["calc"]


def pick(pin_rows, spec_id, peptide):
    for r in pin_rows.get(spec_id, []):
        if r["peptide"] == peptide:
            return r
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--q", type=float, default=0.01)
    ap.add_argument(
        "--te",
        type=float,
        default=None,
        help="target/entrapment peptide ratio; computed from --fasta (slow) when absent",
    )
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    d = a.dir
    truth_pin = read_pin(os.path.join(d, "truth.pin"))
    open_pin = read_pin(os.path.join(d, "open.pin"))
    truth_psms = read_percolator(os.path.join(d, "truth.psms"))
    open_psms = read_percolator(os.path.join(d, "open.psms"))

    truth_ok = {k: v for k, v in truth_psms.items() if v[0] <= a.q}
    open_ok = {k: v for k, v in open_psms.items() if v[0] <= a.q}

    # Ground truth: oxidised-M PSMs of the narrow search, with the M position.
    truth_ox = {}
    for spec_id, (_q, peptide, _prots) in truth_ok.items():
        bare, mods = parse_peptide(peptide)
        sites = [pos for pos, delta in mods if abs(delta - OXIDATION) <= 0.01]
        if sites:
            truth_ox[spec_id] = (bare, sites)

    recovered = []
    for spec_id, (bare, sites) in truth_ox.items():
        hit = open_ok.get(spec_id)
        if hit is None:
            continue
        open_bare, _ = parse_peptide(hit[1])
        if open_bare != bare:
            continue
        row = pick(open_pin, spec_id, hit[1])
        if row is None:
            continue
        delta = pin_delta(row)  # modification mass (observed - peptide)
        if abs(delta - OXIDATION) > OX_TOL_DA:
            continue
        site_match = row["open_site"] != 0 and (row["open_site"] - 1) in sites
        recovered.append(
            {
                "spec": spec_id,
                "delta": delta,
                "open_site": row["open_site"],
                "truth_sites": [p + 1 for p in sites],
                "site_match": site_match,
            }
        )

    # Delta histogram of every accepted open-search PSM (target and decoy alike;
    # the decoys show the random background).
    hist = collections.Counter()
    localized = 0
    for spec_id, (q, peptide, _prots) in open_ok.items():
        row = pick(open_pin, spec_id, peptide)
        if row is None:
            continue
        delta = pin_delta(row)  # modification mass (observed - peptide)
        hist[round(delta)] += 1
        if row["open_site"] != 0:
            localized += 1

    if a.te is None:
        t, e, moved = entrap_common.te_ratio(a.fasta)
        te = t / max(e, 1)
    else:
        te, t, e, moved = a.te, None, None, None

    def fdp(path):
        n, ent, res = entrap_common.count_percolator_psms(
            path, a.q, entrap_common.SampleSequences(a.fasta)
        )
        return {"accepted": n, "entrapment_hits": ent, "rescued": res,
                "fdp": (ent / n * (1 + te)) if n else 0.0}

    # The open-search accepted set is dominated by zero-offset rows. Split the
    # entrapment count by |delta| so the offset-bearing rows get their own error
    # estimate: a single global q does not condition on delta.
    sample = entrap_common.SampleSequences(a.fasta)
    by_offset = {"zero": {"accepted": 0, "entrapment_hits": 0},
                 "offset": {"accepted": 0, "entrapment_hits": 0}}
    for spec_id, (_q, peptide, prots) in open_ok.items():
        row = pick(open_pin, spec_id, peptide)
        if row is None:
            continue
        delta = pin_delta(row)
        acc = [x for cell in prots for x in cell.split(";") if x]
        kind = entrap_common.classify(acc, peptide, sample)
        bucket = by_offset["zero" if abs(delta) <= 0.5 else "offset"]
        bucket["accepted"] += 1
        bucket["entrapment_hits"] += kind == "entrap"
    for bucket in by_offset.values():
        n = bucket["accepted"]
        bucket["fdp"] = (bucket["entrapment_hits"] / n * (1 + te)) if n else 0.0

    report = {
        "q": a.q,
        "te": {"target_peptides": t, "entrapment_peptides": e, "moved": moved, "ratio": te},
        "te_given": a.te is not None,
        "truth": {"psms": len(truth_ok), "oxidised_m": len(truth_ox),
                  **fdp(os.path.join(d, "truth.psms"))},
        "open": {"psms": len(open_ok), "localized": localized,
                 **fdp(os.path.join(d, "open.psms")), "by_offset": by_offset},
        "recovery": {
            "recovered": len(recovered),
            "of": len(truth_ox),
            "rate": len(recovered) / len(truth_ox) if truth_ox else 0.0,
            "site_called": sum(1 for r in recovered if r["open_site"] != 0),
            "site_matches": sum(1 for r in recovered if r["site_match"]),
            "delta_abs_error_median": (
                statistics.median(abs(r["delta"] - OXIDATION) for r in recovered)
                if recovered else 0.0
            ),
        },
        "delta_hist_top": hist.most_common(12),
        "recovered_rows": recovered[:50],
    }
    with open(a.out, "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps({k: report[k] for k in ("truth", "open", "recovery")}, indent=1))


if __name__ == "__main__":
    main()
