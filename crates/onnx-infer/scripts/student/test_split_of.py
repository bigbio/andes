#!/usr/bin/env python3
"""The RT and MS2 student corpora assign train / val / test with the same rule.

Both hash the stripped sequence with I folded to L, so an I/L-isobaric pair never
straddles splits and a peptide held out of one model is held out of the other.
Requires: numpy, pyarrow. Run with pytest or as a script.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "ms2"))
from ms2_common import split_of as ms2_split_of  # noqa: E402
from rt_extract import split_of as rt_split_of  # noqa: E402

SEQS = [
    "PEPTIDEK",
    "PEPTLDEK",
    "IIIIK",
    "LLLLK",
    "AGIVLSQK",
    "SAMPLERIK",
    "ELVISLIVESK",
] + ["".join("ACDEFGHIKLMNPQRSTVWY"[(i * 7 + j * 3) % 20] for j in range(12)) for i in range(500)]


def test_rt_and_ms2_splits_agree():
    for seq in SEQS:
        assert rt_split_of(seq) == ms2_split_of(seq), seq


def test_split_folds_isoleucine_to_leucine():
    for seq in SEQS:
        assert rt_split_of(seq) == rt_split_of(seq.replace("I", "L")), seq


if __name__ == "__main__":
    test_rt_and_ms2_splits_agree()
    test_split_folds_isoleucine_to_leucine()
    print("split tests passed")
