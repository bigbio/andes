"""Readers for published Skyline results of cross-links.

Skyline before its native cross-link support (Skyline 4.x) held a cross-link as one pseudo-peptide:
peptide A, a spacer residue that carries the linker mass (a modification such as `Kxlink`), then
peptide B, with the two linked residues marked by another modification (such as `LinkSite`).
`read_document` turns such a document's peak areas into the long table of table.py; the names of
the two modifications are read from datasets.json, because each document defines its own.
`site_isomers` counts the precursors of such a document that differ only in a link site.

`read_cv_report` reads a Skyline report of per-transition CVs (`File Name`, `Peptide Modified
Sequence`, `Product Charge`, `Fragment Ion`, `Cv Area`), which is all some datasets publish. The
CV is written either as a percentage (`11.80%`) or as a fraction (`0.118`), depending on the
report.
"""
import csv
import xml.etree.ElementTree as ET

import pandas as pd

from .keys import acc_of, pair_key, side
from .table import validate

IDENTIFIED = "true"


def _split(peptide, spacer_mod, site_mod):
    """The two peptides of a pseudo-peptide, each as (sequence, link site, modifications), or None
    when it does not have exactly one spacer with one link site on each side of it.

    Link sites and modification positions are 0-based inside each peptide; the modifications are
    the variable ones other than the spacer and the link sites, as (position, name)."""
    seq = peptide.get("sequence")
    spacer, sites, other = [], [], []
    for m in peptide.findall("./variable_modifications/variable_modification"):
        name, at = m.get("modification_name", ""), int(m.get("index_aa"))
        if name == spacer_mod:
            spacer.append(at)
        elif name.startswith(site_mod):
            sites.append(at)
        else:
            other.append((at, name))
    sites.sort()
    if len(spacer) != 1 or len(sites) != 2 or not (sites[0] < spacer[0] < sites[1]):
        return None
    s = spacer[0]
    alpha = (seq[:s], sites[0], tuple(sorted(m for m in other if m[0] < s)))
    beta = (seq[s + 1:], sites[1] - s - 1, tuple(sorted((at - s - 1, name) for at, name in other if at > s)))
    return alpha, beta


def _peaks(precursor):
    """(peak element, area) of every peak of a precursor with an area above zero, and the counts
    of the peaks left out."""
    kept, no_area, zero_area = [], 0, 0
    for peak in precursor.findall("./precursor_results/precursor_peak"):
        area = peak.get("area")
        if area is None:
            no_area += 1
        elif float(area) <= 0:
            zero_area += 1
        else:
            kept.append((peak, float(area)))
    return kept, no_area, zero_area


def read_document(path, proteins, spacer_mod, site_mod):
    """Long table of a Skyline document whose cross-links are pseudo-peptides.

    Returns (table, counts). A peak counts as identified when Skyline marks it `identified="true"`
    (an MS2 identification in that run); `aligned` (an identification placed by retention-time
    alignment) and `false` count as transferred. Decoys, peptides without exactly one spacer and
    one link site on each side of it, peaks without an area and peaks with an area of zero are
    left out and counted: a zero area has no log2 value, so it is a missing value here.
    """
    root = ET.parse(path).getroot()
    rows, counts = [], {"peptides": 0, "decoy": 0, "not_crosslink": 0, "unmapped": 0, "peaks": 0,
                        "no_area": 0, "zero_area": 0}
    for pep in root.iter("peptide"):
        counts["peptides"] += 1
        if pep.get("decoy") == "true":
            counts["decoy"] += 1
            continue
        parts = _split(pep, spacer_mod, site_mod)
        if parts is None:
            counts["not_crosslink"] += 1
            continue
        (alpha, a_site, _), (beta, b_site, _) = parts
        a_side, b_side = side(proteins, alpha, a_site), side(proteins, beta, b_site)
        if not a_side or not b_side:
            counts["unmapped"] += 1
            continue
        rp = pair_key(a_side, b_side)
        for prec in pep.findall("./precursor"):
            key = f"{prec.get('modified_sequence') or pep.get('modified_sequence')}/{prec.get('charge')}"
            kept, no_area, zero_area = _peaks(prec)
            counts["peaks"] += len(kept) + no_area + zero_area
            counts["no_area"] += no_area
            counts["zero_area"] += zero_area
            for peak, area in kept:
                rows.append({
                    "run": peak.get("replicate"),
                    "precursor": key,
                    "residue_pair": rp,
                    "quantity": area,
                    "evidence": "identified" if peak.get("identified") == IDENTIFIED else "transferred",
                    "isotope_dotp": float(peak.get("isotope_dotp", "nan")),
                    "alpha_peptide": alpha, "alpha_site": a_site + 1,
                    "beta_peptide": beta, "beta_site": b_site + 1,
                })
    return validate(pd.DataFrame(rows), str(path)), counts


def site_isomers(path, spacer_mod, site_mod):
    """Count the precursors of a document that differ from another one only in a link site.

    A site-isomer set is the precursors with the same two peptides, the same modifications and the
    same charge but different link sites: they are isobaric. `overlapping` counts the precursors
    whose integration window overlaps, in at least one run, the window of another site of their
    set; MS1 cannot tell the sites apart there. Only peaks with an area above zero are read.
    """
    root = ET.parse(path).getroot()
    groups, precursors = {}, 0
    for pep in root.iter("peptide"):
        if pep.get("decoy") == "true":
            continue
        parts = _split(pep, spacer_mod, site_mod)
        if parts is None:
            continue
        first, second = sorted(parts, key=lambda x: (x[0], x[2], x[1]))
        sites = (first[1], second[1])
        for prec in pep.findall("./precursor"):
            precursors += 1
            windows = {}
            for peak, _ in _peaks(prec)[0]:
                if peak.get("start_time") and peak.get("end_time"):
                    windows[peak.get("replicate")] = (float(peak.get("start_time")), float(peak.get("end_time")))
            key = (first[0], first[2], second[0], second[2], prec.get("charge"))
            groups.setdefault(key, []).append((sites, windows))
    sets = [members for members in groups.values() if len({s for s, _ in members}) > 1]

    def overlaps(a, b):
        return any(run in b and a[run][0] < b[run][1] and b[run][0] < a[run][1] for run in a)

    overlapping = sum(
        any(other_sites != sites and overlaps(windows, other) for other_sites, other in members)
        for members in sets for sites, windows in members)
    return {"precursors": precursors, "site_isomer_sets": len(sets),
            "precursors_in_sets": sum(len(m) for m in sets), "overlapping": int(overlapping)}


def xifdr_residue_pairs(path):
    """Residue-pair keys of the target-target matches of an xiFDR PSM file (columns `isTT`,
    `Protein1`, `ProteinLinkPos1`, `Protein2`, `ProteinLinkPos2`; several proteins of one peptide
    are joined by `;` in both columns)."""
    x = pd.read_csv(path)
    tt = x[x["isTT"].astype(str).str.lower() == "true"]

    def one(accs, positions):
        return ";".join(sorted({f"{acc_of(a)}:{int(float(p))}"
                                for a, p in zip(str(accs).split(";"), str(positions).split(";"))}))

    return {pair_key(one(a, p), one(b, q))
            for a, p, b, q in zip(tt["Protein1"], tt["ProteinLinkPos1"], tt["Protein2"], tt["ProteinLinkPos2"])}


def read_cv_report(path):
    """One CV (fraction) per (precursor, fragment ion) of a Skyline CV report."""
    out = {}
    with open(path, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            text = r["Cv Area"].strip()
            if not text:
                continue
            cv = float(text[:-1]) / 100 if text.endswith("%") else float(text)
            out[(f"{r['Peptide Modified Sequence']}/{r['Product Charge']}", r["Fragment Ion"])] = cv
    return pd.DataFrame([{"precursor": k[0], "ion": k[1], "cv": v} for k, v in out.items()])
