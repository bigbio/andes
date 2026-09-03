#!/usr/bin/env python3
"""Offline re-ranking gate for the glyco selector redesign.

WHY THIS EXISTS. Building a selector term and then measuring it end to end costs a
full search plus Percolator, and the measurement floor is coarse: a five-seed design
on the pooled benchmark resolves only about 117 PSMs. Most candidate terms can be
killed far earlier and far more cheaply, because the question "does this term rank the
true candidate first among the candidates the engine already generated" is answerable
from a `--debug-glyco` dump alone, with no engine changes and no FDR tool.

So this script re-ranks an EXISTING dump under alternative scoring rules and reports
truth-top-1 rate and margin for each, overall AND on the low-margin stratum (the scans
the shipped selector decides by the smallest margin, where the wrong winners live). A
rule that does not lift top-1 on the low-margin stratum does not get built.

WHAT IT IS NOT. Top-1 rate is not identifications at 1% FDR. A rule can rank truth
first more often and still lose identifications, because Percolator's reportable count
depends on target-decoy separation over the whole run, not on per-scan ordering. Use
this to REJECT terms cheaply, never to claim a yield.

INPUT
  --pin        one or more PINs written by `andes --glyco --debug-glyco` (all candidate
               rows per scan, not the top-1 collapse). Rows are keyed by (file, scan):
               a pooled multi-fraction run re-uses scan numbers across fractions, so
               scan alone would merge unrelated scans.
  --truth      TSV with columns `scan` and `peptide` (and `file` when more than one
               PIN is given) holding the REFERENCE answer per scan. Peptides are
               compared after stripping modification brackets and flanking residues,
               because the engine writes inline modifications and reference tools do
               not; comparing them raw counts every carbamidomethyl peptide as a
               disagreement, which has already produced a wrong number in this project.

               The truth MUST come from an independent engine (MSFragger-Glyco,
               Byonic, pGlyco). Feeding this script the engine's own collapsed PIN is
               circular: the shipped rule then scores 100% by construction. The script
               refuses a truth file that looks like an andes PIN and warns when the
               shipped rule is suspiciously perfect.

The `shipped` rule is the emitted ROW ORDER within a scan: the dump is already sorted
by the engine's fused selector score (which includes terms the PIN does not carry, such
as the ETD c/z weight), so re-deriving it from columns would be a worse copy. Every
other rule is defined over columns the PIN already carries, so a rule that scores well
here is implementable without guessing.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict

# Inline modifications (`C+57.02146`), flanking residues (`K.PEPTIDE.R`) and the glyco
# suffix (`[HexNAc4...@N2]`) all have to come off before two engines' peptide strings
# can be compared. Skipping this step is a measured source of a wrong answer here: an
# earlier comparison counted 81 agreements where there were 130.
_MOD = re.compile(r"[+-]\d+\.\d+")
_GLYCO = re.compile(r"\[[^\]]*\]")


def bare_peptide(p: str) -> str:
    p = _GLYCO.sub("", p)
    p = _MOD.sub("", p)
    if len(p) > 4 and p[1] == "." and p[-2] == ".":
        p = p[2:-2]
    return p.strip().upper()


def f(row: dict, col: str, default: float = 0.0) -> float:
    try:
        return float(row.get(col, default))
    except (TypeError, ValueError):
        return default


def backbone_mass(row: dict) -> float:
    """Neutral PEPTIDE mass of a candidate. `CalcMass` is the INTACT glycopeptide mass,
    which every candidate of a scan shares to within ppm, so keying a split on it puts
    all ~600 candidates in one bucket and makes any multiplicity term a constant."""
    return f(row, "CalcMass") - f(row, "GlycanMass")


def split_ids(rows: list[dict], tol_ppm: float) -> list[int]:
    """Single-linkage clustering of backbone masses within `tol_ppm`, mirroring
    `andes_glyco::glyco_psm::split_ids_by_clustering`. A rounding grid was tried first
    and put masses 5 ppm apart in different splits whenever they straddled a boundary."""
    masses = [backbone_mass(r) for r in rows]
    order = sorted(range(len(rows)), key=lambda i: masses[i])
    ids = [0] * len(rows)
    current = 0
    prev = None
    for i in order:
        m = masses[i]
        joins = prev is not None and m > 0 and abs(m - prev) <= m * tol_ppm * 1e-6
        if prev is not None and not joins:
            current += 1
        ids[i] = current
        prev = m
    return ids


# ---------------------------------------------------------------------------
# Scoring rules. Each takes the candidate rows of ONE scan (in emitted order) and
# returns a score per row. Rules see only PIN columns, so a rule that scores well here
# is implementable.
# ---------------------------------------------------------------------------

def rule_shipped(rows):
    """Emitted row order = the engine's own fused selector (rank + K*ladder + J*core_y
    + H*hyper + CZ*cz). Its two heaviest terms are functions of the backbone mass and
    composition, not of the peptide, so within a mass split they are shared constants."""
    return [float(len(rows) - i) for i in range(len(rows))]


def rule_rank_only(rows):
    return [f(r, "RankScoreFloat") for r in rows]


def rule_rawscore(rows):
    """RawScore alone: measured to rank truth at median 2 while the shipped selector
    ranks it at median 15, yet it is computed after the collapse and never consulted."""
    return [f(r, "RawScore") for r in rows]


def rule_rank_plus_raw(rows):
    return [f(r, "RankScoreFloat") + f(r, "RawScore") for r in rows]


def rule_plus_anchor(rows):
    """Adds the shipped peptide-mass-conditioned Y0/Y1 anchor (`Y0Y1Anchor`). This is
    NOT the doc's exclusive-assignment Y0/Y1/Y2 anchor, which has no PIN column yet;
    it is the closest quantity the dump carries."""
    return [
        f(r, "RankScoreFloat") + f(r, "RawScore") + 5.0 * f(r, "Y0Y1Anchor") for r in rows
    ]


def rule_plus_ytree(rows):
    """Adds the composition-specific Y-tree LLR. Requires --glyco-y-tree in the dump."""
    return [
        f(r, "RankScoreFloat") + f(r, "RawScore") + f(r, "YTreeLLR") for r in rows
    ]


def rule_plus_oxonium(rows):
    """Adds per-candidate oxonium-composition consistency (`OxoniumCompLLR`, requires
    --glyco-oxonium-llr in the dump). `SialicConsistency` is not used: it only flips
    the sign of observed intensity and carries no absence penalty."""
    return [
        f(r, "RankScoreFloat")
        + f(r, "RawScore")
        + f(r, "YTreeLLR")
        + f(r, "OxoniumCompLLR")
        for r in rows
    ]


def rule_plus_chance_masked(rows):
    """Adds the peptide-channel chance LLR (step 6; requires --glyco-chance-llr-masked
    in the dump)."""
    return [
        f(r, "RankScoreFloat")
        + f(r, "RawScore")
        + f(r, "YTreeLLR")
        + f(r, "OxoniumCompLLR")
        + f(r, "ChanceLlrMasked")
        for r in rows
    ]


def make_split_multiplicity(tol_ppm: float):
    def rule_split_multiplicity(rows):
        """The full redesign rule: score as above, then penalise each candidate by the
        log of how many DISTINCT PEPTIDES share its backbone-mass split.

        Without that penalty a split holding 40 sequon peptides draws 40 samples of the
        score distribution while a split holding 2 draws 2, so the argmax rewards how
        dense the peptide window is rather than how good the evidence is. This is the
        mechanism behind the measurement that 96.9% of decoy winners sit at a different
        backbone mass than the truth.
        """
        ids = split_ids(rows, tol_ppm)
        peptides_per_split: dict[int, set] = defaultdict(set)
        for r, sid in zip(rows, ids):
            peptides_per_split[sid].add(bare_peptide(r.get("Peptide", "")))
        out = []
        for r, sid in zip(rows, ids):
            n = len(peptides_per_split[sid])
            out.append(
                f(r, "RankScoreFloat")
                + f(r, "RawScore")
                + f(r, "YTreeLLR")
                + f(r, "OxoniumCompLLR")
                - math.log(max(n, 1))
            )
        return out

    return rule_split_multiplicity


def build_rules(tol_ppm: float):
    return {
        "shipped": rule_shipped,
        "rank_only": rule_rank_only,
        "rawscore_only": rule_rawscore,
        "rank+raw": rule_rank_plus_raw,
        "rank+raw+anchor": rule_plus_anchor,
        "rank+raw+ytree": rule_plus_ytree,
        "rank+raw+ytree+ox": rule_plus_oxonium,
        "...+chance_masked": rule_plus_chance_masked,
        "...+multiplicity": make_split_multiplicity(tol_ppm),
    }


def scan_key(row: dict, file_tag: str) -> tuple[str, str]:
    """(file, scan). `ScanNr` is 0 on MGF input (the TITLE-fallback defect), so the
    SpecId's `scan=` token is used when ScanNr carries nothing."""
    scan = str(row.get("ScanNr", "")).strip()
    if scan in ("", "0"):
        m = re.search(r"scan=(\d+)", row.get("SpecId", ""))
        scan = m.group(1) if m else row.get("SpecId", "")
    return (file_tag, scan)


def evaluate(fn, present: dict, truth: dict):
    """Per-scan truth rank and (for top-1 scans) margin over the best DIFFERENT
    peptide. Returns {scan_key: (rank, margin_or_None)}."""
    out = {}
    for key, rows in present.items():
        scores = fn(rows)
        order = sorted(range(len(rows)), key=lambda i: (-scores[i], rows[i].get("Peptide", "")))
        truth_pos = next(
            (k for k, i in enumerate(order)
             if bare_peptide(rows[i].get("Peptide", "")) == truth[key]),
            None,
        )
        if truth_pos is None:
            continue
        margin = None
        winner_pep = bare_peptide(rows[order[0]].get("Peptide", ""))
        runner = next(
            (i for i in order[1:] if bare_peptide(rows[i].get("Peptide", "")) != winner_pep),
            None,
        )
        if runner is not None:
            margin = scores[order[0]] - scores[runner]
        out[key] = (truth_pos + 1, margin)
    return out


def summarise(res: dict, keys, n_total: int):
    ranks = [res[k][0] for k in keys if k in res]
    top1 = sum(1 for r in ranks if r == 1)
    top3 = sum(1 for r in ranks if r <= 3)
    margins = sorted(res[k][1] for k in keys if k in res and res[k][0] == 1 and res[k][1] is not None)
    med_rank = sorted(ranks)[len(ranks) // 2] if ranks else 0
    med_margin = margins[len(margins) // 2] if margins else 0.0
    n = max(n_total, 1)
    return top1, 100.0 * top1 / n, 100.0 * top3 / n, med_rank, med_margin


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pin", required=True, nargs="+",
                    help="--debug-glyco PIN(s) (all candidate rows). With several, the "
                         "truth file needs a `file` column matching each PIN's basename stem.")
    ap.add_argument("--truth", required=True, help="TSV with columns scan, peptide[, file]")
    ap.add_argument("--tol-ppm", type=float, default=20.0,
                    help="backbone-mass split tolerance for the multiplicity rule")
    ap.add_argument("--min-candidates", type=int, default=2,
                    help="skip uncontested scans; they cannot discriminate between rules")
    ap.add_argument("--low-margin-frac", type=float, default=1.0 / 3.0,
                    help="fraction of scans (lowest shipped-rule margin) forming the "
                         "low-margin stratum")
    ap.add_argument("--rules", default=None, help="comma-separated subset")
    args = ap.parse_args()

    rules = build_rules(args.tol_ppm)
    multi = len(args.pin) > 1

    # Circularity guard 1: an andes PIN as truth.
    with open(args.truth, newline="") as fh:
        first = fh.readline()
    if "SpecId" in first.split("\t") or "RankScore" in first:
        print("truth file looks like an andes PIN (SpecId/RankScore header): that is the "
              "engine grading itself. Use an independent engine's identifications.",
              file=sys.stderr)
        return 2

    truth: dict[tuple[str, str], str] = {}
    with open(args.truth, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if multi and "file" not in (rd.fieldnames or []):
            print("several PINs given but the truth file has no `file` column", file=sys.stderr)
            return 2
        for row in rd:
            tag = str(row.get("file", "")).strip() if multi else ""
            truth[(tag, str(row["scan"]).strip())] = bare_peptide(row["peptide"])
    if not truth:
        print("truth file has no rows", file=sys.stderr)
        return 2

    by_scan: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for pin in args.pin:
        tag = os.path.basename(pin).split(".")[0] if multi else ""
        with open(pin, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = scan_key(row, tag)
                if key in truth:
                    by_scan[key].append(row)

    scans = {s: rs for s, rs in by_scan.items() if len(rs) >= args.min_candidates}
    if not scans:
        print("no scan in the dump matched the truth file with enough candidates", file=sys.stderr)
        return 2

    # A scan only counts toward top-1 if the true peptide is PRESENT among the
    # candidates. A rule cannot be blamed for a candidate that was never generated, and
    # conflating the two is how a selection problem gets misread as a coverage problem.
    present = {
        s: rs for s, rs in scans.items()
        if any(bare_peptide(r.get("Peptide", "")) == truth[s] for r in rs)
    }

    print(f"scans in dump matching truth : {len(scans)}")
    print(f"  ... with truth generated   : {len(present)} "
          f"({100.0 * len(present) / max(len(scans), 1):.1f}%)  <- the ceiling for every rule")
    print(f"  median candidates/scan     : "
          f"{sorted(len(r) for r in present.values())[len(present) // 2] if present else 0}")

    # Low-margin stratum: the scans the SHIPPED rule decides by the smallest margin.
    shipped = evaluate(rules["shipped"], present, truth)
    with_margin = [(k, v[1]) for k, v in shipped.items() if v[1] is not None]
    with_margin.sort(key=lambda kv: kv[1])
    n_low = max(1, int(round(len(with_margin) * args.low_margin_frac)))
    low_keys = [k for k, _ in with_margin[:n_low]]
    print(f"  low-margin stratum         : {len(low_keys)} scans (lowest "
          f"{100.0 * args.low_margin_frac:.0f}% of shipped-rule margins)")

    # Circularity guard 2: a perfect shipped rule is a sign of self-grading.
    shipped_top1 = sum(1 for v in shipped.values() if v[0] == 1)
    if present and shipped_top1 >= 0.99 * len(present):
        print("WARNING: the shipped rule ranks truth first on >=99% of scans. Either the "
              "truth is this engine's own output (circular) or the scans are not "
              "contested; nothing below can discriminate between rules.", file=sys.stderr)
    print()
    header = (f"{'rule':22s} {'top1':>6s} {'top1%':>7s} {'top3%':>7s} {'med.rank':>9s} "
              f"{'med.margin':>11s} | {'low.top1':>8s} {'low.top1%':>9s} {'low.med.rank':>12s}")
    print(header)
    print("-" * len(header))

    names = [n.strip() for n in args.rules.split(",")] if args.rules else list(rules)
    for name in names:
        fn = rules.get(name)
        if fn is None:
            print(f"unknown rule {name!r}; known: {', '.join(rules)}", file=sys.stderr)
            return 2
        res = shipped if name == "shipped" else evaluate(fn, present, truth)
        t1, p1, p3, mr, mm = summarise(res, present.keys(), len(present))
        lt1, lp1, _, lmr, _ = summarise(res, low_keys, len(low_keys))
        print(f"{name:22s} {t1:6d} {p1:6.1f}% {p3:6.1f}% {mr:9d} {mm:11.3f} | "
              f"{lt1:8d} {lp1:8.1f}% {lmr:12d}")

    print()
    print("Read this as a REJECTION filter, on the LOW-MARGIN columns first. A rule that")
    print("does not raise low-margin top-1 will not raise identifications at 1% FDR. A rule")
    print("that does still has to be measured end to end, pooled, over five seeds, with")
    print("entrapment FDP, on BOTH regimes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
