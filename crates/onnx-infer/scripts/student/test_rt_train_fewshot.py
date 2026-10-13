#!/usr/bin/env python3
"""rt_train.select_rows: few-shot PTM training data is bounded by --fewshot-n.

Requires: numpy, pandas, pyarrow, torch. Run with pytest or as a script.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt_train import select_rows  # noqa: E402


def table():
    """One run; 5 TMT and 5 phospho peptidoforms plus 20 unmodified, all in train."""
    forms = ([f"[UNIMOD:737]-PEPTIDE{i}K[UNIMOD:737]" for i in range(5)]
             + [f"PEPS[UNIMOD:21]IDE{i}K" for i in range(5)]
             + [f"PEPTIDE{i}R" for i in range(20)])
    rows = []
    for rep in range(2):
        for f in forms:
            rows.append(dict(dataset="PXD1", run=f"r{rep}", peptidoform=f, split="train",
                             rt_min=10.0, has_tmt="UNIMOD:737" in f, has_phospho="UNIMOD:21" in f))
    for f in forms[:3] + forms[5:7] + forms[10:12]:
        rows.append(dict(dataset="PXD1", run="r0", peptidoform=f + "V", split="val", rt_min=10.0,
                         has_tmt="UNIMOD:737" in f, has_phospho="UNIMOD:21" in f))
    t = pd.DataFrame(rows)
    for c in ["dataset", "run", "peptidoform", "split"]:
        t[c] = t[c].astype("category")
    t["run_id"] = t.groupby(["dataset", "run"], observed=True, sort=True).ngroup()
    t["holdout"] = False
    return t


def args(ptm):
    return argparse.Namespace(exclude_phospho=False, exclude_tmt=False, fewshot_ptm=ptm, fewshot_n=2)


def check(ptm, col):
    train, val, few = select_rows(table(), args(ptm), np.random.default_rng(0))
    assert few[few[col]].peptidoform.nunique() == 2
    assert not train[col].any(), f"train holds {train[train[col]].peptidoform.nunique()} {ptm} forms"
    assert not val[col].any()
    other = "has_phospho" if col == "has_tmt" else "has_tmt"
    assert train[other].any(), "the other PTM stays in train"


def test_fewshot_tmt_keeps_unselected_tmt_forms_out_of_train():
    check("tmt", "has_tmt")


def test_fewshot_phospho_keeps_unselected_phospho_forms_out_of_train():
    check("phospho", "has_phospho")


def test_without_fewshot_nothing_is_excluded():
    train, _, few = select_rows(table(), args(None), np.random.default_rng(0))
    assert few is None and train.has_tmt.any() and train.has_phospho.any()


if __name__ == "__main__":
    test_fewshot_tmt_keeps_unselected_tmt_forms_out_of_train()
    test_fewshot_phospho_keeps_unselected_phospho_forms_out_of_train()
    test_without_fewshot_nothing_is_excluded()
    print("rt_train few-shot selection tests passed")
