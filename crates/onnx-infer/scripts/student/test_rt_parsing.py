#!/usr/bin/env python3
"""Python andes-peptide parsing matches the Rust parser on andes writer output.

``fixtures/andes_writer_peptides.tsv`` holds strings produced by andes's
peptide writer (checked by ``tests/andes_writer_peptides.rs``) and their
expected ProForma. Requires: numpy, torch. Run with pytest or as a script.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod_composition import match_mass  # noqa: E402
from rt_model import parse_andes, parse_proforma  # noqa: E402

FIXTURE = Path(__file__).resolve().parents[2] / "fixtures" / "andes_writer_peptides.tsv"


def writer_cases():
    lines = FIXTURE.read_text().splitlines()[1:]
    return [tuple(line.split("\t")) for line in lines]


def test_tmtpro_is_not_itraq8plex():
    assert match_mass(304.20715, "K", False) == 2016
    assert match_mass(304.20715, None, True) == 2016
    assert match_mass(304.20536, "K", False) == 730
    assert parse_andes("K.P+304.20715EPTIDEK+304.20715.R") == ("PEPTIDEK", [(0, 2016), (8, 2016)])


def test_summed_terminal_and_residue_deltas_from_the_andes_writer():
    cases = writer_cases()
    assert cases
    for andes, proforma in cases:
        seq, mods = parse_proforma(proforma)
        assert parse_andes(andes) == (seq, sorted(mods)), andes


def test_split_deltas_still_parse():
    assert parse_andes("K.M+42.01056+15.99491PEPS+79.96633K.D") == ("MPEPSK", [(0, 1), (1, 35), (5, 21)])
    assert parse_andes("K.PEP+1.23456K.D") is None


if __name__ == "__main__":
    test_tmtpro_is_not_itraq8plex()
    test_summed_terminal_and_residue_deltas_from_the_andes_writer()
    test_split_deltas_still_parse()
    print("rt parsing tests passed")
