#!/usr/bin/env python3
"""Smoke test on synthetic data (a few seconds):
  1. the lineage guard refuses a cross-run join and accepts a same-run one;
  2. the quant_q -> p recovery inverts the running mean;
  3. the whole downstream (andes finalize parquet + MaxQuant evidence, design and amounts files) finds the changing species
     (sign-aware TP > 0), no false positive in the null species, and a small precursor error.

    python3 tests/test_smoke.py
"""
import contextlib
import io
import json
import os
import random
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from quant_eval import downstream, readers, scorecard  # noqa: E402
from quant_eval import runguard as G  # noqa: E402

RUNS = ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"]
AA = "ADEFGHKLMNPQRSTVWY"   # no C: MaxQuant keys add the fixed carbamidomethyl


def make_dataset(tmp):
    rnd = random.Random(1)
    rng = np.random.default_rng(1)
    fasta, peptides = [], []
    for sp, n in (("HUMAN", 150), ("YEAST", 20)):   # human majority keeps the all-species run median on the null species
        for i in range(n):
            peps = ["".join(rnd.choice(AA) for _ in range(9)) + "K" for _ in range(3)]
            fasta.append(f">sp|P{sp[0]}{i}|P{sp[0]}{i}_{sp} protein\n{''.join(peps)}\n")
            peptides += [(sp, p, 6 + rng.normal(0, 2)) for p in peps]
    with open(f"{tmp}/db.fasta", "w") as f:
        f.write("".join(fasta))
    with open(f"{tmp}/design.tsv", "w") as f:
        f.write("run\tcondition\n" + "".join(f"{r}\t{r[0]}\n" for r in RUNS))
    with open(f"{tmp}/amounts.tsv", "w") as f:
        f.write("species\tcondition\tamount\nHUMAN\t*\t1\nYEAST\tA\t1\nYEAST\tB\t2\n")
    rows, ev = [], []
    for sp, pep, base in peptides:
        for r in RUNS:
            lq = base + (1.0 if (sp == "YEAST" and r[0] == "B") else 0.0) + rng.normal(0, 0.2)
            rows.append((pep + "/2", r, "identified_ms2", 2.0 ** lq, 0.0, 0.0, 0.0))
            ev.append((pep, r, "MULTI-MSMS", 2.0 ** lq, ""))
    pd.DataFrame(rows, columns=["precursor_key", "run", "evidence_type", "quantity", "quant_q", "transfer_q", "library_q"]).to_parquet(f"{tmp}/final.parquet")
    with open(f"{tmp}/evidence.txt", "w") as f:
        f.write("Modified sequence\tCharge\tRaw file\tType\tIntensity\tReverse\tPotential contaminant\n")
        for pep, r, t, q, _ in ev:
            f.write(f"_{pep}_\t2\t{r}\t{t}\t{q}\t\t\n")
    return tmp


class Smoke(unittest.TestCase):
    def test_guard(self):
        a = G.tagged(pd.DataFrame({"k": [1, 2], "x": [1.0, 2.0]}), "A:1")
        b = G.tagged(pd.DataFrame({"k": [1, 2], "y": [3.0, 4.0]}), "B:1")
        with self.assertRaises(AssertionError):
            G.join(a, b, "k")
        self.assertEqual(G.join(a, G.tagged(b, "A:1"), "k").shape[0], 2)
        with self.assertRaises(AssertionError):
            G.assert_cols_agree([1.0, 2.0], [1.0, 3.0], tol=1e-6, what="demo")

    def test_recover_p(self):
        p = np.array([0.1, 0.3, 0.5])
        tr = pd.DataFrame({"run": "r", "quant_q": np.cumsum(p) / np.arange(1, 4)})
        self.assertTrue(np.allclose(readers.recover_p(tr).values, p))
        # rows in any order: p must stay attached to its own row
        p2 = np.array([0.05, 0.2, 0.9])
        tr2 = pd.concat([tr.assign(p=p), pd.DataFrame({"run": "r2", "quant_q": np.cumsum(p2) / np.arange(1, 4), "p": p2})])
        tr2 = tr2.sample(frac=1, random_state=3).reset_index(drop=True)
        self.assertTrue(np.allclose(readers.recover_p(tr2).values, tr2.p.values))

    def test_downstream(self):
        with tempfile.TemporaryDirectory() as tmp:
            make_dataset(tmp)
            out = f"{tmp}/res.json"
            downstream.main(["--name", "synthetic", "--design", f"{tmp}/design.tsv", "--amounts", f"{tmp}/amounts.tsv",
                             "--fasta", f"{tmp}/db.fasta", "--classes", "HUMAN,YEAST", "--null", "HUMAN",
                             "--andes", f"andes:default:{tmp}/final.parquet", "--andes", f"andes_V:weights:{tmp}/final.parquet",
                             "--maxquant", f"MQ:{tmp}/evidence.txt", "--out", out])
            with open(out) as f:
                res = json.load(f)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                scorecard.main(["--results", out])
        self.assertIn("DE TP, no cut, andes", buf.getvalue())
        for arm in ("andes", "andes_V", "MQ"):
            de = res["arms"][arm]["de"]
            self.assertGreaterEqual(de["TPs"], 15, arm)      # 20 changing proteins, sign-correct
            self.assertEqual(de["TPs"], de["TP"], arm)
            self.assertLessEqual(de["FP"], 3, arm)           # null species: at most chance-level calls
            self.assertEqual(de["FP5"], 0, arm)
            acc = res["arms"][arm]["acc_own_sp"]["_all"]
            self.assertLess(acc["share_gt1"], 0.01, arm)
            self.assertEqual(res["arms"][arm]["proteins_ge2"], 170, arm)
        self.assertEqual(res["shared_cells"], 170 * 3 * len(RUNS))


if __name__ == "__main__":
    unittest.main()
