"""Tests for the peptidoform identity used by the multi-engine vote check (pytest)."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_select import finalize_rows, peptidoform_key, same_identification  # noqa: E402

SAME = [
    ("[UNIMOD:737]-SHAQK[UNIMOD:737]NENAR", "S[UNIMOD:737]HAQK[UNIMOD:737]NENAR"),
    ("[UNIMOD:737]-K[UNIMOD:737]ILQ", "K[UNIMOD:737][UNIMOD:737]ILQ"),
    ("[UNIMOD:1]-AHEH", "A[UNIMOD:1]HEH"),
    ("M[UNIMOD:35][UNIMOD:737]PDR", "[UNIMOD:737]-M[UNIMOD:35]PDR"),
    ("[UNIMOD:1][UNIMOD:737]-PEPK", "P[UNIMOD:737][UNIMOD:1]EPK"),
]

DIFFERENT = [
    ("[UNIMOD:737]-SHAQK[UNIMOD:737]NENAR", "[UNIMOD:737]-SHAQKNENAR"),
    ("[UNIMOD:737]-SHAQK[UNIMOD:737]NENAR", "[UNIMOD:737]-SHAQK[UNIMOD:737]NENAK"),
    ("[UNIMOD:737]-K[UNIMOD:737]ILQ", "K[UNIMOD:737]ILQ"),
    ("[UNIMOD:1]-AHEH", "AHEH"),
    ("[UNIMOD:1]-AHEH", "[UNIMOD:737]-AHEH"),
    ("PEM[UNIMOD:35]K", "PEMK[UNIMOD:35]"),
    ("PEM[UNIMOD:35]K", "M[UNIMOD:35]PEK"),
    ("PEPS[UNIMOD:21]K", "PEPSK-[UNIMOD:2]"),
    ("[UNIMOD:1]-APEPK", "AP[UNIMOD:1]EPK"),
]


def test_nterm_spellings_are_one_identification():
    for a, b in SAME:
        assert peptidoform_key(a) == peptidoform_key(b), (a, b)
        assert same_identification(a, b) and same_identification(b, a), (a, b)


def test_different_peptidoforms_stay_different():
    for a, b in DIFFERENT:
        assert peptidoform_key(a) != peptidoform_key(b), (a, b)
        assert not same_identification(a, b), (a, b)


def test_carbamidomethyl_is_ignored_on_c_only():
    assert same_identification("PEPC[UNIMOD:4]TIDEK", "PEPCTIDEK"), "fixed CAM left out by a submitter"
    assert not same_identification("PEPTIDEK[UNIMOD:4]", "PEPTIDEK"), "CAM on K is an alkylation side product"
    assert not same_identification("[UNIMOD:4]-PEPTIDEK", "PEPTIDEK"), "N-terminal CAM is a side product"
    assert same_identification("C[UNIMOD:4]PEPTIDEK", "CPEPTIDEK"), "CAM on a first-residue C"


def test_key_keeps_other_positions():
    assert peptidoform_key("[UNIMOD:737]-PEPS[UNIMOD:21]K[UNIMOD:737]") == (
        "PEPSK", ((0, 737), (4, 21), (5, 737)))
    assert peptidoform_key("[UNIMOD:1]-") is None


def test_finalize_counts_votes_written_on_first_residue():
    row = {"run": "r1", "scan": 7, "peptidoform": "[UNIMOD:737]-SHAQK[UNIMOD:737]NENAR", "hard": False,
           "source_engine": "msnet", "intensities": np.zeros(1, np.float32)}
    vote = {("r1", 7): "S[UNIMOD:737]HAQK[UNIMOD:737]NENAR"}
    kept, st = finalize_rows([row], {"comet": vote, "sage": dict(vote), "msfragger": {}})
    assert st["kept"] == 1 and kept[0]["n_engines"] == 2
    assert kept[0]["agree_comet"] and kept[0]["agree_sage"] and kept[0]["agree_msfragger"] is False


def test_finalize_joins_candidate_runs_with_extension():
    row = {"run": "20160804_QE1_proteome_02.mzML", "scan": 11, "peptidoform": "PEPTIDEK", "hard": False,
           "source_engine": "msnet", "intensities": np.zeros(1, np.float32)}
    vote = {("20160804_QE1_proteome_02", 11): "PEPTIDEK"}
    kept, st = finalize_rows([row], {"comet": vote, "sage": dict(vote)})
    assert st["kept"] == 1 and kept[0]["n_engines"] == 2
    assert kept[0]["run"] == "20160804_QE1_proteome_02.mzML"
