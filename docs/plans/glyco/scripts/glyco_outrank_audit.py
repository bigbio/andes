#!/usr/bin/env python3
"""Audit glyco truth scans where the true backbone is generated but out-ranked.

The input PIN should be a diagnostic all-hit file, not the normal collapsed
top-1 glyco PIN. The script compares the best candidate for each truth scan
under the andes collapse selector with the best candidate whose peptide
backbone residue mass matches the truth backbone. It is intended for the glyco
benchmark artifacts on the VM, but it only depends on TSV/PIN inputs.
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


AA_MASS = {
    "G": 57.02146,
    "A": 71.03711,
    "S": 87.03203,
    "P": 97.05276,
    "V": 99.06841,
    "T": 101.04768,
    "C": 103.00919,
    "L": 113.08406,
    "I": 113.08406,
    "N": 114.04293,
    "D": 115.02694,
    "Q": 128.05858,
    "K": 128.09496,
    "E": 129.04259,
    "M": 131.04049,
    "H": 137.05891,
    "F": 147.06841,
    "R": 156.10111,
    "Y": 163.06333,
    "W": 186.07931,
}

CAM_MASS = 57.02146
H2O_MASS = 18.010565
SCAN_RE = re.compile(r"scan=(\d+)")
BRACKET_RE = re.compile(r"\[[^\]]*\]")
PAREN_RE = re.compile(r"\([^)]*\)")


@dataclass(frozen=True)
class Truth:
    scan: int
    backbone_mass: float


@dataclass(frozen=True)
class Candidate:
    scan: int
    spec_id: str
    label: int
    peptide: str
    proteins: str
    backbone_mass: float
    rank_score: float
    rank_score_float: float
    y_ladder: float
    core_y_hits: float
    oxonium: float
    n_core_oxonium: float
    glycan_mass: float
    y0y1_anchor: float
    sialic_consistency: float
    calc_mass: float
    exp_mass: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit generated-but-outranked glyco truth scans.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--truth", required=True, type=Path, help="truth TSV with scan/backbone_mass columns")
    parser.add_argument("--pin", required=True, type=Path, help="all-hit glyco PIN/TSV candidate dump")
    parser.add_argument("--out", type=Path, help="optional per-scan audit TSV")
    parser.add_argument("--tol-da", type=float, default=0.05, help="absolute backbone mass tolerance")
    parser.add_argument(
        "--truth-mass-kind",
        choices=("residue", "neutral", "auto"),
        default="residue",
        help="truth backbone mass convention",
    )
    parser.add_argument(
        "--selector",
        choices=("yladder", "rank"),
        default="yladder",
        help="collapse selector to audit; andes default is yladder",
    )
    parser.add_argument("--include-decoys", action="store_true", help="include Label=-1 candidates")
    parser.add_argument("--examples", type=int, default=20, help="number of worst outranked examples to print")
    return parser.parse_args()


def require_column(header: list[str], names: tuple[str, ...], path: Path) -> str:
    by_lower = {h.lower(): h for h in header}
    for name in names:
        found = by_lower.get(name.lower())
        if found is not None:
            return found
    raise SystemExit(f"{path}: missing any of columns {', '.join(names)}")


def load_truth(path: Path, mass_kind: str) -> dict[int, Truth]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"{path}: empty truth file")
        scan_col = require_column(reader.fieldnames, ("scan", "ScanNr", "scan_nr", "scan_number"), path)
        mass_col = require_column(
            reader.fieldnames,
            ("backbone_mass", "backbone_residue_mass", "peptide_residue_mass", "mass"),
            path,
        )
        truth = {}
        for row in reader:
            if not row.get(scan_col) or not row.get(mass_col):
                continue
            scan = int(float(row[scan_col]))
            mass = float(row[mass_col])
            if mass_kind == "neutral":
                mass -= H2O_MASS
            truth[scan] = Truth(scan=scan, backbone_mass=mass)
    return truth


def strip_flanks(peptide: str) -> str:
    no_glycan = BRACKET_RE.sub("", peptide)
    no_paren = PAREN_RE.sub("", no_glycan)
    parts = no_paren.split(".")
    if len(parts) >= 3:
        return parts[1]
    return parts[0]


def residue_mass(peptide: str) -> float | None:
    """Return residue-sum mass, including inline numeric mods and fixed CAM."""
    core = strip_flanks(peptide)
    mass = 0.0
    i = 0
    while i < len(core):
        aa = core[i].upper()
        if aa not in AA_MASS:
            i += 1
            continue

        mass += AA_MASS[aa]
        i += 1

        explicit_mods = []
        while i < len(core) and core[i] in "+-":
            sign = 1.0 if core[i] == "+" else -1.0
            i += 1
            start = i
            while i < len(core) and (core[i].isdigit() or core[i] == "."):
                i += 1
            if start != i:
                mod = sign * float(core[start:i])
                mass += mod
                explicit_mods.append(mod)

        if aa == "C" and not any(abs(mod - CAM_MASS) <= 0.01 for mod in explicit_mods):
            mass += CAM_MASS

    return mass if mass > 0.0 else None


def parse_scan(spec_id: str) -> int | None:
    match = SCAN_RE.search(spec_id)
    return int(match.group(1)) if match else None


def to_float(row: dict[str, str], name: str, default: float = 0.0) -> float:
    value = row.get(name)
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def to_int(row: dict[str, str], name: str, default: int = 0) -> int:
    value = row.get(name)
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


def load_candidates(path: Path, truth_scans: set[int], include_decoys: bool) -> dict[int, list[Candidate]]:
    by_scan: dict[int, list[Candidate]] = defaultdict(list)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"{path}: empty PIN file")
        pep_col = require_column(reader.fieldnames, ("Peptide", "peptide"), path)
        protein_col = next((c for c in ("Proteins", "proteinIds", "proteinids") if c in reader.fieldnames), "")

        for row in reader:
            spec_id = row.get("SpecId") or row.get("PSMId") or ""
            scan = parse_scan(spec_id)
            if scan is None or scan not in truth_scans:
                continue
            label = to_int(row, "Label", 1)
            if label < 0 and not include_decoys:
                continue
            peptide = row.get(pep_col, "")
            # Backbone mass numerically as CalcMass - GlycanMass - H2O (the peptide
            # residue mass), which is the authoritative convention. Re-parsing the
            # peptide string with residue_mass() is unreliable — it mis-parses ~46%
            # of multi-mod / long peptides (e.g. a 3947 Da backbone read as 146 Da),
            # which inflates the truth_absent bucket with false absents. Fall back to
            # string parsing only when the numeric columns are missing.
            calc = to_float(row, "CalcMass")
            gly = to_float(row, "GlycanMass")
            backbone = calc - gly - H2O_MASS if calc > 0.0 else residue_mass(peptide)
            if backbone is None:
                continue
            by_scan[scan].append(
                Candidate(
                    scan=scan,
                    spec_id=spec_id,
                    label=label,
                    peptide=peptide,
                    proteins=row.get(protein_col, "") if protein_col else "",
                    backbone_mass=backbone,
                    rank_score=to_float(row, "RankScore"),
                    rank_score_float=to_float(row, "RankScoreFloat"),
                    y_ladder=to_float(row, "YLadderScore"),
                    core_y_hits=to_float(row, "CoreYHits"),
                    oxonium=to_float(row, "OxoniumScore"),
                    n_core_oxonium=to_float(row, "NCoreOxoniumIons"),
                    glycan_mass=to_float(row, "GlycanMass"),
                    y0y1_anchor=to_float(row, "Y0Y1Anchor"),
                    sialic_consistency=to_float(row, "SialicConsistency"),
                    calc_mass=to_float(row, "CalcMass"),
                    exp_mass=to_float(row, "ExpMass"),
                )
            )
    return by_scan


def selector_key(candidate: Candidate, selector: str) -> tuple[float, float, str]:
    if selector == "yladder":
        return (candidate.y_ladder, candidate.rank_score, candidate.spec_id)
    return (candidate.rank_score, candidate.y_ladder, candidate.spec_id)


def mass_matches(candidate: Candidate, truth: Truth, tol_da: float) -> bool:
    return abs(candidate.backbone_mass - truth.backbone_mass) <= tol_da


def best(candidates: list[Candidate], selector: str) -> Candidate:
    return max(candidates, key=lambda c: selector_key(c, selector))


def reason_for(winner: Candidate, true_best: Candidate, selector: str) -> str:
    rank_gap = winner.rank_score - true_best.rank_score
    y_gap = winner.y_ladder - true_best.y_ladder
    core_gap = winner.core_y_hits - true_best.core_y_hits
    anchor_gap = winner.y0y1_anchor - true_best.y0y1_anchor
    eps = 1e-9
    if selector == "yladder":
        if y_gap > eps:
            return "truth_loses_y_ladder"
        if abs(y_gap) <= eps and rank_gap > eps:
            return "yladder_tie_loses_rank"
    else:
        if rank_gap > eps:
            return "truth_loses_rank"
        if abs(rank_gap) <= eps and y_gap > eps:
            return "rank_tie_loses_y_ladder"
    if core_gap > eps:
        return "truth_loses_core_y_hits"
    if anchor_gap > eps:
        return "truth_loses_y0y1_anchor"
    return "tie_or_unmodeled_tiebreak"


def fmt(value: float) -> str:
    return f"{value:.6g}"


def row_for(
    truth: Truth,
    candidates: list[Candidate],
    selector: str,
    tol_da: float,
) -> dict[str, str]:
    truth_candidates = [c for c in candidates if mass_matches(c, truth, tol_da)]
    winner = best(candidates, selector) if candidates else None
    rank_winner = best(candidates, "rank") if candidates else None
    y_winner = best(candidates, "yladder") if candidates else None

    if not candidates:
        return {
            "scan": str(truth.scan),
            "status": "no_candidates",
            "truth_backbone_mass": fmt(truth.backbone_mass),
            "n_candidates": "0",
            "n_truth_candidates": "0",
        }
    if not truth_candidates:
        delta = winner.backbone_mass - truth.backbone_mass
        return {
            "scan": str(truth.scan),
            "status": "truth_absent",
            "truth_backbone_mass": fmt(truth.backbone_mass),
            "n_candidates": str(len(candidates)),
            "n_truth_candidates": "0",
            "winner_spec_id": winner.spec_id,
            "winner_peptide": winner.peptide,
            "winner_proteins": winner.proteins,
            "winner_backbone_mass": fmt(winner.backbone_mass),
            "winner_backbone_delta": fmt(delta),
            "winner_rank_score": fmt(winner.rank_score),
            "winner_y_ladder": fmt(winner.y_ladder),
            "winner_core_y_hits": fmt(winner.core_y_hits),
            "winner_glycan_mass": fmt(winner.glycan_mass),
            "rank_top_correct": str(False),
            "yladder_top_correct": str(False),
        }

    true_best = best(truth_candidates, selector)
    winner_is_true = mass_matches(winner, truth, tol_da)
    status = "top1_correct" if winner_is_true else "truth_outranked"
    reason = "" if winner_is_true else reason_for(winner, true_best, selector)
    return {
        "scan": str(truth.scan),
        "status": status,
        "reason": reason,
        "truth_backbone_mass": fmt(truth.backbone_mass),
        "n_candidates": str(len(candidates)),
        "n_truth_candidates": str(len(truth_candidates)),
        "winner_spec_id": winner.spec_id,
        "winner_peptide": winner.peptide,
        "winner_proteins": winner.proteins,
        "winner_backbone_mass": fmt(winner.backbone_mass),
        "winner_backbone_delta": fmt(winner.backbone_mass - truth.backbone_mass),
        "winner_rank_score": fmt(winner.rank_score),
        "winner_y_ladder": fmt(winner.y_ladder),
        "winner_core_y_hits": fmt(winner.core_y_hits),
        "winner_glycan_mass": fmt(winner.glycan_mass),
        "winner_y0y1_anchor": fmt(winner.y0y1_anchor),
        "winner_sialic_consistency": fmt(winner.sialic_consistency),
        "truth_spec_id": true_best.spec_id,
        "truth_peptide": true_best.peptide,
        "truth_proteins": true_best.proteins,
        "truth_backbone_delta": fmt(true_best.backbone_mass - truth.backbone_mass),
        "truth_rank_score": fmt(true_best.rank_score),
        "truth_y_ladder": fmt(true_best.y_ladder),
        "truth_core_y_hits": fmt(true_best.core_y_hits),
        "truth_glycan_mass": fmt(true_best.glycan_mass),
        "truth_y0y1_anchor": fmt(true_best.y0y1_anchor),
        "truth_sialic_consistency": fmt(true_best.sialic_consistency),
        "rank_gap": fmt(winner.rank_score - true_best.rank_score),
        "y_ladder_gap": fmt(winner.y_ladder - true_best.y_ladder),
        "core_y_gap": fmt(winner.core_y_hits - true_best.core_y_hits),
        "y0y1_gap": fmt(winner.y0y1_anchor - true_best.y0y1_anchor),
        "rank_top_correct": str(mass_matches(rank_winner, truth, tol_da)),
        "yladder_top_correct": str(mass_matches(y_winner, truth, tol_da)),
    }


def print_summary(rows: list[dict[str, str]], selector: str, tol_da: float) -> None:
    total = len(rows)
    counts = Counter(row.get("status", "") for row in rows)
    reasons = Counter(row.get("reason", "") for row in rows if row.get("reason"))
    rank_correct = sum(1 for row in rows if row.get("rank_top_correct") == "True")
    y_correct = sum(1 for row in rows if row.get("yladder_top_correct") == "True")

    print(f"truth scans              : {total}")
    print(f"selector                 : {selector}")
    print(f"backbone tol             : {tol_da:g} Da")
    for key in ("no_candidates", "truth_absent", "top1_correct", "truth_outranked"):
        n = counts[key]
        print(f"{key:24s}: {n:5d}  ({100 * n / total:5.1f}%)")
    print(f"rank-top correct         : {rank_correct:5d}  ({100 * rank_correct / total:5.1f}%)")
    print(f"yladder-top correct      : {y_correct:5d}  ({100 * y_correct / total:5.1f}%)")

    if reasons:
        print("outrank reasons:")
        for reason, n in reasons.most_common():
            print(f"  {reason:22s}: {n}")

    outranked = [row for row in rows if row.get("status") == "truth_outranked"]
    for field in ("rank_gap", "y_ladder_gap", "core_y_gap", "y0y1_gap"):
        vals = [float(row[field]) for row in outranked if row.get(field)]
        if vals:
            print(
                f"{field:24s}: mean={statistics.fmean(vals):.4g}  "
                f"median={statistics.median(vals):.4g}"
            )


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "scan",
        "status",
        "reason",
        "truth_backbone_mass",
        "n_candidates",
        "n_truth_candidates",
        "winner_spec_id",
        "winner_peptide",
        "winner_proteins",
        "winner_backbone_mass",
        "winner_backbone_delta",
        "winner_rank_score",
        "winner_y_ladder",
        "winner_core_y_hits",
        "winner_glycan_mass",
        "winner_y0y1_anchor",
        "winner_sialic_consistency",
        "truth_spec_id",
        "truth_peptide",
        "truth_proteins",
        "truth_backbone_delta",
        "truth_rank_score",
        "truth_y_ladder",
        "truth_core_y_hits",
        "truth_glycan_mass",
        "truth_y0y1_anchor",
        "truth_sialic_consistency",
        "rank_gap",
        "y_ladder_gap",
        "core_y_gap",
        "y0y1_gap",
        "rank_top_correct",
        "yladder_top_correct",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def print_examples(rows: list[dict[str, str]], limit: int) -> None:
    if limit <= 0:
        return
    outranked = [row for row in rows if row.get("status") == "truth_outranked"]
    outranked.sort(key=lambda row: (float(row.get("y_ladder_gap") or 0.0), float(row.get("rank_gap") or 0.0)), reverse=True)
    if not outranked:
        return
    print("worst outranked examples:")
    for row in outranked[:limit]:
        print(
            "  scan {scan}: {reason}; "
            "winner(rank={winner_rank_score}, y={winner_y_ladder}, coreY={winner_core_y_hits}, dm={winner_backbone_delta}) "
            "truth(rank={truth_rank_score}, y={truth_y_ladder}, coreY={truth_core_y_hits}, dm={truth_backbone_delta})".format(
                **row
            )
        )


def main() -> None:
    args = parse_args()
    initial_truth = load_truth(args.truth, "residue" if args.truth_mass_kind == "auto" else args.truth_mass_kind)
    candidates = load_candidates(args.pin, set(initial_truth), args.include_decoys)

    if args.truth_mass_kind == "auto":
        residue_rows = [row_for(t, candidates.get(t.scan, []), args.selector, args.tol_da) for t in initial_truth.values()]
        neutral_truth = {
            scan: Truth(scan=t.scan, backbone_mass=t.backbone_mass - H2O_MASS) for scan, t in initial_truth.items()
        }
        neutral_rows = [row_for(t, candidates.get(t.scan, []), args.selector, args.tol_da) for t in neutral_truth.values()]
        residue_present = sum(1 for row in residue_rows if row.get("status") in {"top1_correct", "truth_outranked"})
        neutral_present = sum(1 for row in neutral_rows if row.get("status") in {"top1_correct", "truth_outranked"})
        truth = initial_truth if residue_present >= neutral_present else neutral_truth
        convention = "residue" if truth is initial_truth else "neutral(-H2O)"
        print(f"auto truth mass convention: {convention}")
    else:
        truth = initial_truth

    rows = [row_for(t, candidates.get(t.scan, []), args.selector, args.tol_da) for t in truth.values()]
    rows.sort(key=lambda row: int(row["scan"]))
    print_summary(rows, args.selector, args.tol_da)
    print_examples(rows, args.examples)
    if args.out:
        write_rows(args.out, rows)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
