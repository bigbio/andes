#!/usr/bin/env python3
"""Decompose the glyco `truth_absent` bucket — the P0 half that
`glyco_outrank_audit.py` does not cover.

`glyco_outrank_audit.py` tells you a truth backbone was *absent* from the
scoreable candidate surface, but not WHY. This tool distinguishes the two causes
the handoff (SESSION-HANDOFF §5-P0/§6.1) asks for, and characterises the miss by
precursor charge and backbone mass (the miss skews large / high-charge).

Cause attribution needs TWO all-hit PINs of the SAME run set:
  --pin           baseline all-hit PIN (e.g. the honest 253/97 candidate surface)
  --pin-expanded  a max-retention / precursor-mass-on all-hit PIN (optional)

For each truth scan absent in --pin:
  * present in --pin-expanded  -> "retention_recoverable"
        (top-k / SHORTLIST_K=24 / retention dropped it; the DB *can* enumerate it)
  * absent in both             -> "never_enumerated"
        (no source proposed the backbone: no precursor-mass DB source, or truly
         unsearchable — needs leg-1 precursor-mass enumeration to even appear)

Without --pin-expanded it still reports absent-vs-present and the charge / mass /
mass-split characterisation, but cannot split the two causes.

Backbone matching is numeric on RESIDUE mass (no water), matching andes'
convention and `glyco_outrank_audit.py`. PIN peptides carry Cam-C explicitly as
`C+57.02146`; the parser adds Cam only when not already present.
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
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
    "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
CAM_MASS = 57.02146
H2O_MASS = 18.010565
SCAN_RE = re.compile(r"scan=(\d+)")
BRACKET_RE = re.compile(r"\[[^\]]*\]")
PAREN_RE = re.compile(r"\([^)]*\)")
CHARGE_COL_RE = re.compile(r"^charge(\d+)$", re.IGNORECASE)


@dataclass(frozen=True)
class Truth:
    scan: int
    backbone_mass: float


@dataclass
class Cand:
    backbone_mass: float
    glycan_mass: float
    charge: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Decompose glyco truth_absent into retention-loss vs never-enumerated.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--truth", required=True, type=Path, help="truth TSV (scan + backbone_mass)")
    p.add_argument("--pin", required=True, type=Path, help="baseline all-hit glyco PIN")
    p.add_argument("--pin-expanded", type=Path, help="max-retention / precursor-mass-on all-hit PIN")
    p.add_argument("--out", type=Path, help="optional per-absent-scan TSV")
    p.add_argument("--tol-da", type=float, default=0.05, help="absolute backbone mass tolerance")
    p.add_argument(
        "--truth-mass-kind", choices=("residue", "neutral"), default="residue",
        help="truth backbone mass convention",
    )
    p.add_argument("--include-decoys", action="store_true", help="include Label=-1 rows")
    return p.parse_args()


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
            ("backbone_mass", "backbone_residue_mass", "peptide_residue_mass", "mass"), path,
        )
        truth: dict[int, Truth] = {}
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
    return parts[1] if len(parts) >= 3 else parts[0]


def residue_mass(peptide: str) -> float | None:
    """Residue-sum mass, including inline numeric mods and fixed Cam-C."""
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
        explicit_mods: list[float] = []
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
    if not value:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def parse_charge(row: dict[str, str], charge_cols: list[tuple[str, int]]) -> int:
    for col, z in charge_cols:
        if row.get(col, "0") == "1":
            return z
    return 0


def load_candidates(
    path: Path, truth_scans: set[int], include_decoys: bool
) -> dict[int, list[Cand]]:
    by_scan: dict[int, list[Cand]] = defaultdict(list)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"{path}: empty PIN file")
        pep_col = require_column(reader.fieldnames, ("Peptide", "peptide"), path)
        charge_cols = sorted(
            ((h, int(CHARGE_COL_RE.match(h).group(1))) for h in reader.fieldnames if CHARGE_COL_RE.match(h)),
            key=lambda t: t[1],
        )
        for row in reader:
            spec_id = row.get("SpecId") or row.get("PSMId") or ""
            scan = parse_scan(spec_id)
            if scan is None or scan not in truth_scans:
                continue
            if not include_decoys and int(to_float(row, "Label", 1)) < 0:
                continue
            # Numeric backbone = CalcMass - GlycanMass - H2O (peptide residue mass),
            # the authoritative convention. The residue_mass() string parser mis-parses
            # ~46% of multi-mod/long peptides and produces false absents; use it only
            # as a fallback when the numeric columns are missing.
            gly = to_float(row, "GlycanMass")
            calc = to_float(row, "CalcMass")
            backbone = calc - gly - H2O_MASS if calc > 0.0 else residue_mass(row.get(pep_col, ""))
            if backbone is None:
                continue
            by_scan[scan].append(
                Cand(
                    backbone_mass=backbone,
                    glycan_mass=gly,
                    charge=parse_charge(row, charge_cols),
                )
            )
    return by_scan


def truth_present(cands: list[Cand], truth: Truth, tol_da: float) -> bool:
    return any(abs(c.backbone_mass - truth.backbone_mass) <= tol_da for c in cands)


def mass_bucket(mass: float) -> str:
    lo = int(mass // 250) * 250
    return f"{lo}-{lo + 250}"


def main() -> None:
    args = parse_args()
    truth = load_truth(args.truth, args.truth_mass_kind)
    base = load_candidates(args.pin, set(truth), args.include_decoys)
    expanded = (
        load_candidates(args.pin_expanded, set(truth), args.include_decoys)
        if args.pin_expanded else None
    )

    rows: list[dict[str, str]] = []
    for scan, t in sorted(truth.items()):
        base_cands = base.get(scan, [])
        present = truth_present(base_cands, t, args.tol_da)
        if present:
            continue  # covered by glyco_outrank_audit.py (top1_correct / truth_outranked)

        # Best charge/mass-split context from the baseline winner-ish surface:
        # use the strongest-glycan competitor as a stand-in for the wrong winner.
        wrong = max(base_cands, key=lambda c: c.glycan_mass, default=None)
        cause = "unknown"
        if expanded is not None:
            cause = (
                "retention_recoverable"
                if truth_present(expanded.get(scan, []), t, args.tol_da)
                else "never_enumerated"
            )
        rows.append({
            "scan": str(scan),
            "cause": cause,
            "truth_backbone_mass": f"{t.backbone_mass:.4f}",
            "mass_bucket": mass_bucket(t.backbone_mass),
            "n_base_candidates": str(len(base_cands)),
            "wrong_charge": str(wrong.charge if wrong else 0),
            "wrong_backbone_mass": f"{wrong.backbone_mass:.4f}" if wrong else "",
            "wrong_glycan_mass": f"{wrong.glycan_mass:.4f}" if wrong else "",
            "short_backbone_big_glycan": str(
                bool(wrong and wrong.backbone_mass < t.backbone_mass - args.tol_da
                     and wrong.glycan_mass > 0.0)
            ),
        })

    total_truth = len(truth)
    absent = len(rows)
    print(f"truth scans                 : {total_truth}")
    print(f"absent in baseline --pin    : {absent}  ({100 * absent / max(total_truth, 1):.1f}%)")
    if expanded is not None:
        causes = Counter(r["cause"] for r in rows)
        for key in ("retention_recoverable", "never_enumerated"):
            n = causes[key]
            print(f"  {key:22s}: {n:4d}  ({100 * n / max(absent, 1):.1f}% of absent)")
        print("  (retention_recoverable = leg-2 retention/shortlist fix; "
              "never_enumerated = needs leg-1 precursor-mass source)")
    else:
        print("  (pass --pin-expanded to split retention-loss vs never-enumerated)")

    charges = Counter(r["wrong_charge"] for r in rows)
    print("absent by wrong-winner charge:")
    for z, n in sorted(charges.items(), key=lambda t: int(t[0])):
        print(f"  z={z or '?'}: {n}")
    buckets = Counter(r["mass_bucket"] for r in rows)
    print("absent by truth backbone mass (Da):")
    for b, n in sorted(buckets.items(), key=lambda t: int(t[0].split('-')[0])):
        print(f"  {b:>10s}: {n}")
    split = sum(1 for r in rows if r["short_backbone_big_glycan"] == "True")
    print(f"short-backbone/big-glycan wrong winner: {split}/{absent}")
    masses = [float(r["truth_backbone_mass"]) for r in rows]
    if masses:
        print(f"absent truth backbone mass: median={statistics.median(masses):.1f} "
              f"mean={statistics.fmean(masses):.1f} Da")

    if args.out:
        fields = ["scan", "cause", "truth_backbone_mass", "mass_bucket", "n_base_candidates",
                  "wrong_charge", "wrong_backbone_mass", "wrong_glycan_mass",
                  "short_backbone_big_glycan"]
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
