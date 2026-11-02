#!/usr/bin/env python3
"""Smoke test on synthetic data (a few seconds):
  1. the lineage guard refuses a cross-run join and accepts a same-run one;
  2. the quant_q -> p recovery inverts the running mean;
  3. scorecard fixes: same-family BH, TP at matched FDP, per-pair AUC, design-consistent detection, absence table, RT-local
     normalisation (unit tests on small synthetic inputs);
  4. the whole downstream (andes finalize parquet + MaxQuant evidence, design and amounts files) finds the changing species
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
from quant_eval import de as DE  # noqa: E402
from quant_eval import downstream, metrics, readers, rollup, scorecard  # noqa: E402
from quant_eval.design import Design  # noqa: E402
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


def ups_design(tmp):
    """UPS1-like design: 6 levels x 2 replicates, amounts = level for UPS, constant for YEAST."""
    levels = [1, 2, 4, 8, 16, 32]
    runs = [f"L{lv}_{i}" for lv in levels for i in (1, 2)]
    with open(f"{tmp}/d.tsv", "w") as f:
        f.write("run\tcondition\tblank\tabsent_species\n" + "".join(f"{r}\t{r[1:].split('_')[0]}\t\t\n" for r in runs) + "BLANK\t\ttrue\t\n")
    with open(f"{tmp}/a.tsv", "w") as f:
        f.write("species\tcondition\tamount\nUPS\t*\t=condition\nYEAST\t*\t1\n")
    return Design(f"{tmp}/d.tsv", f"{tmp}/a.tsv", "YEAST", ("UPS", "YEAST"), low_max=8.0), runs


class Fixes(unittest.TestCase):
    def test_union_view(self):
        # padding adds p = 1 tests, keeps the species and never lowers a q value of the original proteins
        T = pd.DataFrame({"sp": ["UPS"] * 4 + ["YEAST"] * 4, "fc": 1.0, "p": [1e-6, 1e-4, 1e-3, 0.5, 0.9, 0.04, 0.6, 0.8],
                          "ex": [1.0] * 4 + [0.0] * 4}, index=[f"P{i}" for i in range(8)])
        T["q"] = DE.bh(T.p.values)
        uni = DE.union_view(T, {f"P{i}" for i in range(8)} | {"X1", "X2"}, {"X1": "UPS", "X2": "YEAST"})
        self.assertEqual(len(uni), 10)
        self.assertEqual(float(uni.loc["X1", "p"]), 1.0)
        self.assertEqual(uni.loc["X1", "sp"], "UPS")
        self.assertTrue((uni.q.loc[T.index].values >= T.q.values - 1e-12).all())

    def test_matched_fdp_and_auc(self):
        sp = ["UPS"] * 10 + ["YEAST"] * 10
        p = np.r_[np.linspace(1e-6, 1e-3, 9), 0.2, 0.01, np.linspace(0.3, 0.9, 9)]   # one YEAST hit among the top 11
        T = pd.DataFrame({"sp": sp, "p": p, "fc": 1.0, "q": 1.0, "ex": 1.0})
        tp, fp = DE.tp_at_fdp([T], "YEAST", (0.01, 0.10, 0.5))
        self.assertEqual((tp["1"], fp["1"]), (9, 0))                   # ranked: 9 UPS, the YEAST hit, then the 10th UPS
        self.assertEqual((tp["10"], fp["10"]), (10, 1))                # FDP 1/11 <= 10 % keeps the hit and the next UPS
        self.assertGreaterEqual(tp["50"], tp["10"])
        a = DE.auc_changed_vs_null(T, "YEAST")
        self.assertTrue(0.9 < a < 1.0)
        T2 = T.assign(p=np.r_[np.linspace(1e-6, 1e-3, 10), np.linspace(0.3, 0.9, 10)])
        self.assertEqual(DE.auc_changed_vs_null(T2, "YEAST"), 1.0)

    def test_design_detection_counts_consistent_precursors_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            design, runs = ups_design(tmp)
            rng = np.random.default_rng(5)
            rows, sp = [], {}
            for k in range(300):                                           # yeast background fixes the run offsets
                sp[f"Y{k}"] = "YEAST"
                for r in runs:
                    rows.append((f"Y{k}", r, 2.0 ** (20 + rng.normal(0, 0.05)), False))
            for k in range(20):                                            # 20 consistent UPS precursors (quantity follows the level)
                sp[f"U{k}"] = "UPS"
                for r, lv in zip(runs, [float(r[1:].split("_")[0]) for r in runs]):
                    rows.append((f"U{k}", r, 2.0 ** (10 + np.log2(lv) + rng.normal(0, 0.1)), lv < 8 and k % 2 == 0))
            for k in range(20, 30):                                        # 10 flat precursors: seen at every level, not design-consistent
                sp[f"U{k}"] = "UPS"
                for r in runs:
                    rows.append((f"U{k}", r, 2.0 ** (14 + rng.normal(0, 0.1)), False))
            raw = pd.DataFrame(rows, columns=["key", "run", "quantity", "is_transfer"])
            out = metrics.design_detection(raw, design, pd.Series(sp), ref_skip=3)
        lowest = out["1"]
        self.assertEqual(lowest["any"], 30)                                # raw detection counts every UPS precursor with a reference
        self.assertEqual(lowest["t1.0"], 20)                               # only the design-consistent ones remain
        self.assertGreaterEqual(out["32"]["t1.0"], 20)                     # at the top level the flat ones may also fall within 2x
        self.assertLess(out["1"]["t1.0"], out["1"]["any"])

    def test_rt_local_normalisation_removes_an_rt_trend(self):
        rng = np.random.default_rng(2)
        runs = ["A_1", "A_2", "B_1", "B_2"]
        keys = [f"K{i}/2" for i in range(1500)]
        rt = pd.Series(rng.uniform(0, 100, len(keys)), index=keys)
        base = pd.Series(rng.normal(20, 1, len(keys)), index=keys)
        rows = []
        for r in runs:
            trend = 0.8 * np.sin(rt.values / 15) if r == "B_2" else 0.0 * rt.values     # one run with an RT-dependent response
            lq = base.values + trend + rng.normal(0, 0.05, len(keys))
            rows += [(k, r, 2.0 ** v, False) for k, v in zip(keys, lq)]
        raw = pd.DataFrame(rows, columns=["key", "run", "quantity", "is_transfer"])
        g, _, _ = rollup.normalise(raw, runs)
        loc, flag, nref = rollup.normalise_rt_local(raw, runs, rt)
        self.assertEqual(flag, "rtlocal")
        resid = lambda d: (d.n - d.groupby("key").n.transform("median"))[d.run == "B_2"].abs().median()
        self.assertGreater(resid(g), 0.3)
        self.assertLess(resid(loc), 0.1)
        # keys without an RT proxy still get a (global) shift and no NaN
        part = rollup.normalise_rt_local(raw, runs, rt.iloc[:500])[0]
        self.assertTrue(np.isfinite(part.n).all())

    def test_absence_table_is_separate(self):
        with tempfile.TemporaryDirectory() as tmp:
            design, runs = ups_design(tmp)
            raw = pd.DataFrame({"key": ["Y1", "Y2"], "run": ["BLANK", "L1_1"], "quantity": [1.0, 2.0], "is_transfer": [True, False]})
            o = metrics.absence(raw, design, pd.Series({"Y1": "YEAST", "Y2": "YEAST"}))
        self.assertEqual(o["blank"]["blank_rows"], 1)
        self.assertEqual(o["blank"]["blank_transferred"], 1)
        res = {"arms": {"andes": {"absence_dummy": 0}}, "arm_roles": {"andes": {"kind": "andes"}}, "primary": "andes",
               "absence": {"andes": o}, "classes": []}
        sec = scorecard.Section(res)
        rows = sec.absence_rows()
        self.assertFalse(any("blank-run" in r for r in rows))                              # carry-over is not part of the absence score
        self.assertTrue(any("blank-run rows" in r for r in scorecard.Section(res).carry_over_rows()))
        self.assertFalse(any("blank-run" in r for r in scorecard.Section(res).build()))   # nor mixed into the main table


def p14_design(tmp):
    """Six runs: H1-H3 contain no YEAST (absent_species), M1-M2 contain it, BLANK is a blank run."""
    with open(f"{tmp}/d.tsv", "w") as f:
        f.write("run\tabsent_species\tblank\nH1\tYEAST\t\nH2\tYEAST\t\nH3\tYEAST\t\nM1\t\t\nM2\t\t\nBLANK\t\ttrue\n")
    return Design(f"{tmp}/d.tsv", None, None, ("HUMAN", "YEAST"))


class Rates(unittest.TestCase):
    SP = pd.Series({"h1": "HUMAN", "h2": "HUMAN", "h3": "HUMAN", "y1": "YEAST", "y2": "YEAST"})

    def table(self):
        rows = [  # key, run, is_transfer, tq, thr
            ("h1", "H1", True, 0.001, 0.002), ("h2", "H1", True, 0.004, 0.02), ("h3", "H2", True, 0.002, 0.001), ("y1", "H2", True, 0.008, 0.03),
            ("h1", "H3", False, 0.0, 0.0), ("y2", "H3", False, 0.0, 0.0),                               # identified rows: not transfers
            ("y1", "M1", True, 0.003, 0.004), ("y2", "M2", True, 0.002, 0.001), ("h1", "M1", False, 0.0, 0.0),
            ("y2", "BLANK", True, 0.001, 0.001), ("h2", "BLANK", False, 0.0, 0.0)]
        d = pd.DataFrame(rows, columns=["key", "run", "is_transfer", "tq", "thr"])
        d["quantity"] = 1.0
        return d

    def test_false_transfer_rates(self):
        from quant_eval import rates
        with tempfile.TemporaryDirectory() as tmp:
            design = p14_design(tmp)
        o = rates.false_transfer_rates(self.table(), design, self.SP)["YEAST"]
        self.assertEqual((o["transfers_in_absent_runs"], o["false_transfers"]), (4, 1))
        self.assertAlmostEqual(o["rate"], 0.25)
        self.assertEqual(o["genuine_transfers_of_species"], 2)
        self.assertAlmostEqual(o["yeast_ratio"], 0.5)
        self.assertEqual(o["identified_in_absent_runs"], 1)                      # y2 identified in H3: reported, not a transfer
        self.assertAlmostEqual(o["species_share_identified"], 0.5)             # identified keys outside the blank run: h1, y2
        self.assertEqual(o["transfers_all_runs"], 6)                              # the blank-run transfer is carry-over, not counted
        self.assertEqual(rates.carry_over(self.table(), design, self.SP), {"blank_runs": 1, "rows": 2, "transferred": 1, "precursors": 2})

    def test_depth_curve_and_matched_depth(self):
        from quant_eval import rates
        with tempfile.TemporaryDirectory() as tmp:
            design = p14_design(tmp)
        d = self.table()
        cur = rates.depth_curve(d, design, self.SP, "tq", grid=[0.001, 0.003, 0.01])
        self.assertEqual([c["transfers"] for c in cur], [1, 4, 6])             # blank-run transfer never counts
        self.assertEqual([c["YEAST"]["false"] for c in cur], [0, 0, 1])
        self.assertAlmostEqual(cur[2]["YEAST"]["rate"], 0.25)
        self.assertEqual(rates.matched_depth(cur, "YEAST", 0.10), (0.003, 4))   # loosest threshold with rate <= 10 %
        self.assertEqual(rates.matched_depth(cur, "YEAST", 0.25), (0.01, 6))
        self.assertIsNone(rates.matched_depth(rates.depth_curve(d, design, self.SP, "tq", grid=[0.0001]), "YEAST", 0.5))   # nothing kept: no rate, no match

    def test_scorecard_tables(self):
        from quant_eval import rates
        with tempfile.TemporaryDirectory() as tmp:
            design = p14_design(tmp)
        d = self.table()
        res = {"rates": {"V": rates.false_transfer_rates(d, design, self.SP), "FP": rates.false_transfer_rates(d.iloc[:7], design, self.SP)},
               "depth_curves": {"V": {k: rates.depth_curve(d, design, self.SP, c) for k, c in (("transfer_q", "tq"), ("quant_q", "thr"))}}}
        rt = "\n".join(scorecard.rates_table(res))
        self.assertIn("| V | 4 | 1 | 25.00 |", rt)
        dt = "\n".join(scorecard.depth_table(res))
        self.assertIn("| FP |", dt)
        self.assertNotIn("| V |", dt)


class FlashReader(unittest.TestCase):
    def test_read_flashlfq(self):
        hdr = ["File Name", "Full Sequence", "Precursor Charge", "Peak intensity", "Peak Detection Type", "PIP Q-Value", "Decoy Peptide", "Random RT"]
        rows = [["r_1", "PEPTIDEK", "2", "100", "MSMS", "", "False", "False"], ["r_1", "AAAAK", "2", "50", "MBR", "0.004", "False", "False"],
                ["r_2", "AAAAK", "2", "60", "MBR", "0.03", "False", "False"], ["r_2", "AAAAK", "3", "40", "MBR", "0.001", "True", "False"],
                ["r_2", "CCCCK", "2", "30", "MBR", "0.002", "False", "True"], ["r_2", "PEPTIDEK", "2", "0", "MSMS", "", "False", "False"]]
        with tempfile.TemporaryDirectory() as tmp:
            p = f"{tmp}/QuantifiedPeaks.tsv"
            with open(p, "w") as f:
                f.write("\t".join(hdr) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
            d = readers.read_flashlfq(p, ["r_1", "r_2"])
            d1 = readers.read_flashlfq(p, ["r_1", "r_2"], pip_q=0.01)
        self.assertEqual(sorted(zip(d.key, d.run, d.is_transfer)), [("AAAAK/2", "r_1", True), ("AAAAK/2", "r_2", True), ("PEPTIDEK/2", "r_1", False)])
        self.assertEqual(sorted(d.loc[d.is_transfer, "tq"]), [0.004, 0.03])                      # PIP q-value kept for transfers; decoy / random-RT / zero rows dropped
        self.assertEqual(sorted(d1.key + d1.run), ["AAAAK/2r_1", "PEPTIDEK/2r_1"])               # pip_q = 0.01 drops the 0.03 peak
        self.assertTrue((d.loc[~d.is_transfer, "tq"] == 0).all())


class LibraryGate(unittest.TestCase):
    def test_library_q_limit_per_arm(self):
        d = pd.DataFrame({"precursor_key": ["AAAAK/2", "CCCCK/2"], "run": ["r_1", "r_1"], "evidence_type": ["identified_ms2", "identified_ms2"],
                          "quantity": [1.0, 2.0], "quant_q": [0.0, 0.0], "transfer_q": [0.0, 0.0], "library_q": [0.004, 0.03]})
        with tempfile.TemporaryDirectory() as tmp:
            p = f"{tmp}/f.parquet"
            d.to_parquet(p)
            a = readers.read_andes(p, ["r_1"], "default")
            b = readers.read_andes(p, ["r_1"], "default", library_q=0.05)
            args = downstream.parse_args(["--name", "x", "--design", "d", "--andes", f"wide:default@0.05:{p}", "--out", "o"])
            tabs, kinds = downstream.load_arms(args, ["r_1"])
        self.assertEqual(list(a.key), ["AAAAK/2"])
        self.assertEqual(sorted(b.key), ["AAAAK/2", "CCCCK/2"])
        self.assertEqual(len(tabs["wide"]), 2)
        self.assertEqual(kinds["wide"], {"kind": "andes", "rule": "default", "library_q": 0.05})


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
        for title in ("DE by BH family", "True positives at matched FDP", "intersection family"):
            self.assertIn(title, buf.getvalue())
        self.assertEqual(res["families"]["intersection"], 170)
        for arm in ("andes", "MQ"):
            fam = res["arms"][arm]["de_family"]
            self.assertEqual(fam["inter"]["tested"], fam["own"]["tested"])          # identical protein sets here
            self.assertGreater(fam["inter"]["auc_mean"], 0.9)
            self.assertLessEqual(fam["inter"]["tp_at_fdp"]["1"], fam["inter"]["tp_at_fdp"]["10"])
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
