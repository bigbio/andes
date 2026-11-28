"""Readers for published Skyline results of cross-links.

Skyline before its native cross-link support (Skyline 4.x) held a cross-link as one pseudo-peptide:
peptide A, a spacer residue that carries the linker mass (a modification such as `Kxlink`), then
peptide B, with the two linked residues marked by another modification (such as `LinkSite`).
`read_document` turns such a document's peak areas into the long table of table.py; the names of
the two modifications are read from datasets.json, because each document defines its own.

`read_cv_report` reads a Skyline report of per-transition CVs (`File Name`, `Peptide Modified
Sequence`, `Product Charge`, `Fragment Ion`, `Cv Area`), which is all some datasets publish. The
CV is written either as a percentage (`11.80%`) or as a fraction (`0.118`), depending on the
report.
"""
import csv
import xml.etree.ElementTree as ET

import pandas as pd

from .keys import pair_key, side
from .table import validate

IDENTIFIED = "true"


def _sites(peptide, spacer_mod, site_mod):
    spacer, sites = [], []
    for m in peptide.findall("./variable_modifications/variable_modification"):
        name, at = m.get("modification_name", ""), int(m.get("index_aa"))
        if name == spacer_mod:
            spacer.append(at)
        elif name.startswith(site_mod):
            sites.append(at)
    return spacer, sorted(sites)


def read_document(path, proteins, spacer_mod, site_mod):
    """Long table of a Skyline document whose cross-links are pseudo-peptides.

    Returns (table, counts). A peak counts as identified when Skyline marks it `identified="true"`
    (an MS2 identification in that run); `aligned` (an identification placed by retention-time
    alignment) and `false` count as transferred. Decoys, peptides without exactly one spacer and
    one link site on each side of it, and peaks without an area are left out and counted.
    """
    root = ET.parse(path).getroot()
    rows, counts = [], {"peptides": 0, "decoy": 0, "not_crosslink": 0, "unmapped": 0, "peaks": 0, "no_area": 0}
    for pep in root.iter("peptide"):
        counts["peptides"] += 1
        if pep.get("decoy") == "true":
            counts["decoy"] += 1
            continue
        seq = pep.get("sequence")
        spacer, sites = _sites(pep, spacer_mod, site_mod)
        if len(spacer) != 1 or len(sites) != 2 or not (sites[0] < spacer[0] < sites[1]):
            counts["not_crosslink"] += 1
            continue
        s = spacer[0]
        alpha, beta = seq[:s], seq[s + 1:]
        a_site, b_site = sites[0], sites[1] - s - 1
        a_side, b_side = side(proteins, alpha, a_site), side(proteins, beta, b_site)
        if not a_side or not b_side:
            counts["unmapped"] += 1
            continue
        rp = pair_key(a_side, b_side)
        for prec in pep.findall("./precursor"):
            key = f"{prec.get('modified_sequence') or pep.get('modified_sequence')}/{prec.get('charge')}"
            for peak in prec.findall("./precursor_results/precursor_peak"):
                counts["peaks"] += 1
                area = peak.get("area")
                if area is None:
                    counts["no_area"] += 1
                    continue
                rows.append({
                    "run": peak.get("replicate"),
                    "precursor": key,
                    "residue_pair": rp,
                    "quantity": float(area),
                    "evidence": "identified" if peak.get("identified") == IDENTIFIED else "transferred",
                    "isotope_dotp": float(peak.get("isotope_dotp", "nan")),
                    "alpha_peptide": alpha, "alpha_site": a_site + 1,
                    "beta_peptide": beta, "beta_site": b_site + 1,
                })
    return validate(pd.DataFrame(rows), str(path)), counts


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
