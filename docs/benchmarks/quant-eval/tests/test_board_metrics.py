#!/usr/bin/env python3
"""Tests of the board metrics (lowlevel.py, absence.py, the board columns) on small synthetic inputs (seconds):

    python3 tests/test_board_metrics.py
"""
import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from quant_eval import absence as AB  # noqa: E402
from quant_eval import board, lowlevel as LL, metrics  # noqa: E402
from quant_eval.design import Design  # noqa: E402

LEVELS = [1, 2, 4, 8, 16, 32]
RUNS = [f"L{lv}_{i}" for lv in LEVELS for i in (1, 2)]


def ups_design(tmp):
    with open(f"{tmp}/d.tsv", "w") as f:
        f.write("run\tcondition\n" + "".join(f"{r}\t{r[1:].split('_')[0]}\n" for r in RUNS))
    with open(f"{tmp}/a.tsv", "w") as f:
        f.write("species\tcondition\tamount\nUPS\t*\t=condition\nYEAST\t*\t1\n")
    return Design(f"{tmp}/d.tsv", f"{tmp}/a.tsv", "YEAST", ("UPS", "YEAST"), low_max=8.0)


def ups_raw(seed=5, n_good=20, n_flat=10):
    """Yeast background; n_good consistent UPS precursors (identified), n_flat precursors flat at 2^14 (transferred, quant_q rising with their index)."""
    rng = np.random.default_rng(seed)
    rows, sp = [], {}
    for k in range(300):
        sp[f"Y{k}"] = "YEAST"
        for r in RUNS:
            rows.append((f"Y{k}", r, 2.0 ** (20 + rng.normal(0, 0.05)), False, 0.0))
    for k in range(n_good):
        sp[f"U{k}"] = "UPS"
        for r in RUNS:
            rows.append((f"U{k}", r, 2.0 ** (10 + np.log2(float(r[1:].split("_")[0])) + rng.normal(0, 0.1)), False, 0.0))
    for k in range(n_good, n_good + n_flat):
        sp[f"U{k}"] = "UPS"
        for r in RUNS:
            rows.append((f"U{k}", r, 2.0 ** (14 + rng.normal(0, 0.1)), True, 0.001 * (k + 1)))
    return pd.DataFrame(rows, columns=["key", "run", "quantity", "is_transfer", "thr"]), pd.Series(sp)


class L1Detail(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.design = ups_design(self.tmp.name)
        self.raw, self.sp = ups_raw()

    def tearDown(self):
        self.tmp.cleanup()

    def test_levels_roles_and_agreement_with_design_detection(self):
        d = LL.l1_detail(self.raw, self.design, self.sp)
        self.assertEqual(d["l1_levels"], ["1", "2", "4"])
        self.assertEqual(d["ref_levels"], ["8", "16", "32"])
        self.assertEqual(d["self_reference"], [])
        self.assertEqual(d["per_level"]["1"]["role"], "L1")
        self.assertEqual(d["per_level"]["16"]["role"], "ref")
        old = metrics.design_detection(self.raw, self.design, self.sp)
        for L in d["levels"]:
            self.assertEqual(d["per_level"][L]["keys_within2x"], old[L]["t1.0"])     # the L1 count of the board is unchanged
            self.assertEqual(d["per_level"][L]["rows"], old[L]["cells"])

    def test_self_reference_is_flagged_when_levels_overlap(self):
        d = LL.l1_detail(self.raw, self.design, self.sp, ref_skip=2, n_l1=3)
        self.assertEqual(d["self_reference"], ["4"])

    def test_precision_and_error(self):
        d = LL.l1_detail(self.raw, self.design, self.sp)["per_level"]["1"]
        self.assertEqual(d["rows"], 60)                                  # 30 precursors x 2 runs
        self.assertEqual(d["keys_within2x"], 20)
        self.assertEqual(d["rows_within2x"] if "rows_within2x" in d else d["within2x"], 40)
        self.assertAlmostEqual(d["share_within2x"], 40 / 60)
        self.assertGreater(d["median_abs_err"], 0)
        self.assertEqual(d["transfer"]["rows"], 20)                      # the flat precursors
        self.assertEqual(d["transfer"]["within2x"], 0)                   # flat precursors are 2^14 / 2^10 above the design at level 1
        self.assertGreater(d["transfer"]["median_err"], 3)               # over-quantified (positive)
        self.assertEqual(d["identified"]["share_within2x"], 1.0)

    def test_order_and_matched_counts(self):
        self.assertTrue(LL.has_quality_order(self.raw))
        self.assertFalse(LL.has_quality_order(self.raw.assign(thr=0.0)))
        d = LL.l1_detail(self.raw, self.design, self.sp)
        pl = d["per_level"]["1"]
        # quality order: the 40 identified rows first (all within 2x), then the transfers (none within 2x)
        self.assertEqual(pl["row_ranks"], list(range(1, 41)))
        c = LL.counts_at(pl, 10)
        self.assertEqual((c["rows"], c["rows_within2x"]), (10, 10))
        self.assertLessEqual(c["keys"], 10)
        self.assertEqual(LL.counts_at(pl, 1000)["keys"], 20)

    def test_intensity_order_for_arms_without_quality(self):
        d = LL.l1_detail(self.raw.assign(thr=0.0), self.design, self.sp)
        self.assertEqual(d["order"], "intensity")
        # brightest first: the flat precursors (2^14) precede the consistent ones (2^10 x level) at level 1
        self.assertGreater(d["per_level"]["1"]["row_ranks"][0], 20)

    def test_matched_reading_limits_arms_to_the_same_rows(self):
        a = LL.l1_detail(self.raw, self.design, self.sp)
        b = LL.l1_detail(self.raw[self.raw.key.map(self.sp) != "UPS"].pipe(lambda x: x), self.design, self.sp)   # no UPS rows: {}
        self.assertEqual(b, {})
        half = self.raw[~self.raw.key.isin([f"U{k}" for k in range(10, 20)] + [f"U{k}" for k in range(25, 30)])]
        h = LL.l1_detail(half, self.design, self.sp)
        k, out = LL.matched_reading({"full": a, "half": h}, ["full", "half"], "1")
        self.assertEqual(k, h["per_level"]["1"]["rows"])
        self.assertEqual(out["half"]["rows"], k)
        self.assertEqual(out["full"]["rows"], k)                           # the full arm is cut to the half arm's row count
        self.assertEqual(out["half"]["keys"], 10)
        self.assertLessEqual(out["full"]["keys"], 20)
        self.assertEqual(LL.matched_reading({"full": a}, ["nobody"], "1"), (0, {}))

    def test_expected_keys_random(self):
        # 4 rows, one precursor with 1 correct row: drawing 2 rows finds it with probability 2/4; drawing all rows always
        lvl = {"rows": 4, "key_correct_rows": {"1": 1}}
        self.assertAlmostEqual(LL.expected_keys_random(lvl, 2), 0.5)
        self.assertAlmostEqual(LL.expected_keys_random(lvl, 4), 1.0)
        # two correct rows of 4, k = 3: C(2,3) = 0 -> always found
        self.assertAlmostEqual(LL.expected_keys_random({"rows": 4, "key_correct_rows": {"2": 3}}, 3), 3.0)
        # k = 1 of 4 rows with c = 2: 1 - C(2,1)/C(4,1) = 0.5
        self.assertAlmostEqual(LL.expected_keys_random({"rows": 4, "key_correct_rows": {"2": 1}}, 1), 0.5)
        d = LL.l1_detail(self.raw, self.design, self.sp)["per_level"]["1"]
        self.assertEqual(sum(d["key_correct_rows"].values()), d["keys_within2x"])
        self.assertAlmostEqual(LL.expected_keys_random(d, d["rows"]), d["keys_within2x"])
        self.assertLess(LL.expected_keys_random(d, 10), d["keys_within2x"])

    def test_count_only_metric_rewards_rows(self):
        """Adding wrong rows raises L1 although precision falls: the reason for the companion columns."""
        base = LL.l1_detail(self.raw, self.design, self.sp)["per_level"]["1"]
        rng = np.random.default_rng(1)
        extra = []
        for k in range(30, 330):
            for r in RUNS[:2]:
                extra.append((f"U{k}", r, 2.0 ** rng.uniform(6, 14), True, 0.5))
        sp = pd.concat([self.sp, pd.Series({f"U{k}": "UPS" for k in range(30, 330)})])
        raw = pd.concat([self.raw, pd.DataFrame(extra, columns=self.raw.columns)])
        # give the extra precursors a reference at the high levels so that they are scored
        hi = [(f"U{k}", r, 2.0 ** (10 + np.log2(float(r[1:].split("_")[0]))), False, 0.0) for k in range(30, 330) for r in RUNS if int(r[1:].split("_")[0]) >= 8]
        raw = pd.concat([raw, pd.DataFrame(hi, columns=self.raw.columns)])
        more = LL.l1_detail(raw, self.design, sp)["per_level"]["1"]
        self.assertGreater(more["keys_within2x"], base["keys_within2x"])
        self.assertLess(more["share_within2x"], base["share_within2x"])


class Absence(unittest.TestCase):
    def make(self, tmp):
        runs_abs = [f"H{i}" for i in range(1, 7)]
        runs_pres = [f"M{i}" for i in range(1, 7)]
        with open(f"{tmp}/d.tsv", "w") as f:
            f.write("run\tabsent_species\tblank\n" + "".join(f"{r}\tYEAST\t\n" for r in runs_abs) + "".join(f"{r}\t\t\n" for r in runs_pres) + "BLANK\t\ttrue\n")
        design = Design(f"{tmp}/d.tsv", None, None, ("HUMAN", "YEAST"))
        rng = np.random.default_rng(3)
        rows, sp = [], {}
        for k in range(600):                                   # control: same intensity everywhere
            sp[f"H{k}"] = "HUMAN"
            base = 20 + rng.normal(0, 1)
            for r in runs_abs + runs_pres:
                if r in runs_abs and rng.random() < 0.6:               # a control precursor is seen in a varying number of absent runs
                    continue
                rows.append((f"H{k}", r, 2.0 ** (base + rng.normal(0, 0.3)), False))
        for k in range(40):                                    # genuine yeast: 4 log2 higher where present; two false rows in absent runs
            sp[f"Y{k}"] = "YEAST"
            base = 18 + rng.normal(0, 0.5)
            for r in runs_pres:
                rows.append((f"Y{k}", r, 2.0 ** (base + 4 + rng.normal(0, 0.3)), False))
            for r in runs_abs[:2]:
                rows.append((f"Y{k}", r, 2.0 ** (base + rng.normal(0, 0.3)), k % 2 == 0))
        for k in range(40, 70):                                # mis-labelled: same intensity with and without
            sp[f"Y{k}"] = "YEAST"
            base = 19 + rng.normal(0, 0.5)
            for r in runs_pres + runs_abs[:2]:
                rows.append((f"Y{k}", r, 2.0 ** (base + rng.normal(0, 0.3)), k % 2 == 0))
        for k in range(70, 75):                                # only in absent runs
            sp[f"Y{k}"] = "YEAST"
            rows.append((f"Y{k}", runs_abs[0], 2.0 ** 17, True))
        return design, pd.DataFrame(rows, columns=["key", "run", "quantity", "is_transfer"]), pd.Series(sp)

    def test_profile_and_classes(self):
        with tempfile.TemporaryDirectory() as tmp:
            design, raw, sp = self.make(tmp)
            prof, control = AB.profile(raw, design, sp, "YEAST")
            self.assertEqual(control, "HUMAN")
            ctrl = prof[prof.sp == "HUMAN"]
            self.assertLess(abs(ctrl.d.median()), 0.1)             # control centred at 0 after the per-run control normalisation
            y = prof[prof.sp == "YEAST"].set_index("key")
            self.assertGreater(y.loc["Y3", "d"], 3.5)
            self.assertLess(abs(y.loc["Y45", "d"]), 1.0)
            self.assertEqual(y.loc["Y72", "n_pres"], 0)
            cl = AB.classify(prof[prof.sp == "YEAST"].reset_index(drop=True), None, 1.0)
            by = dict(zip(prof[prof.sp == "YEAST"].key, cl))
            self.assertEqual(by["Y3"], "yeast_enriched")
            self.assertEqual(by["Y45"], "not_specific")
            self.assertEqual(by["Y72"], "only_absent")

    def test_learned_cutoff_follows_the_control_spread(self):
        rng = np.random.default_rng(0)
        d = rng.normal(0, 0.2, 2000)
        n = np.full(2000, 6)
        cuts = AB.learn_cutoffs(d, n, alpha=0.01)
        self.assertAlmostEqual(AB.cutoff_for(cuts, 6), float(np.quantile(d, 0.99)), places=6)
        wide = AB.learn_cutoffs(rng.normal(0, 1.0, 2000), n)
        self.assertGreater(AB.cutoff_for(wide, 6), AB.cutoff_for(cuts, 6))     # a noisier control gives a wider boundary
        self.assertIsNone(AB.learn_cutoffs([], []))
        # a sparse bin falls back to the pooled cutoff
        d2 = np.r_[rng.normal(0, 0.2, 500), rng.normal(0, 2.0, 5)]
        n2 = np.r_[np.full(500, 6), np.full(5, 1)]
        c2 = AB.learn_cutoffs(d2, n2)
        self.assertEqual(AB.cutoff_for(c2, 1), float(np.quantile(d2, 0.99)))

    def test_absence_classes_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            design, raw, sp = self.make(tmp)
            o = AB.absence_classes(raw, design, sp)["YEAST"]
        for name in ("learned", "fixed"):
            c = o[name]
            enr = c["identified"]["yeast_enriched"] + c["transferred"]["yeast_enriched"]
            ns = c["identified"]["not_specific"] + c["transferred"]["not_specific"]
            self.assertEqual(enr + ns, 70 * 2)
            self.assertGreaterEqual(enr, 40 * 2)                           # every row of the 40 genuine precursors; a mislabelled one may cross the boundary by chance
            self.assertLessEqual(enr, 40 * 2 + 4)
            self.assertEqual(c["transferred"]["only_absent"] + c["identified"]["only_absent"], 5)
        self.assertEqual(o["fixed"]["transferred"]["yeast_enriched"], 20 * 2)
        self.assertEqual(len(o["cutoffs"]), 4)
        d = o["denominators"]
        self.assertEqual(d["transferred_rows_absent_runs"] + d["identified_rows_absent_runs"], len(raw[raw.run.str.startswith("H")]))

    def test_blank_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            design, raw, sp = self.make(tmp)
            blank = pd.DataFrame({"key": ["H0", "H1", "H2", "ZZ"], "run": "BLANK", "quantity": [2.0 ** 20 * 0.001, 2.0 ** 20, 2.0 ** 20 * 0.2, 5.0],
                                  "is_transfer": [True, False, False, True]})
            # make H0..H2 base ~ 2^20: take their medians as the reference
            raw = pd.concat([raw, blank])
            o = AB.blank_split(raw, design)
        self.assertEqual(o["rows"], 4)
        self.assertEqual((o["identified"]["rows"], o["transferred"]["rows"]), (2, 2))
        self.assertEqual(o["transferred"]["no_reference"], 1)              # ZZ has no other run
        self.assertEqual(o["carry_over_fixed"], 2)                         # H0 (0.001x) and H2 (0.2x); H1 is at the level of the other runs
        self.assertLess(o["learned_fraction"], 1.0)


class BoardColumns(unittest.TestCase):
    def test_matched_l1_and_new_cells(self):
        with tempfile.TemporaryDirectory() as tmp:
            design = ups_design(tmp)
            raw, sp = ups_raw()
            d_full = LL.l1_detail(raw, design, sp)
            d_int = LL.l1_detail(raw.assign(thr=0.0), design, sp)
        r = {"arms": {"main": {"l1_detail": d_full}, "MQ": {"l1_detail": d_int}, "FP": {"l1_detail": d_int}}}
        m = board.matched_l1(r, ["main"])
        self.assertEqual(set(m), {"main", "MQ", "FP"})
        self.assertLessEqual(m["main"][0], d_full["per_level"]["1"]["keys_within2x"] + d_full["per_level"]["2"]["keys_within2x"] + d_full["per_level"]["4"]["keys_within2x"])
        cells = board.new_cells(r, "main", ["main"], m)
        self.assertEqual(cells[0], "1/2/4")
        self.assertEqual(len(cells), 8)
        self.assertEqual(cells[4], f"{m['main'][0]:,}")
        self.assertEqual(board.new_cells(r, "none", ["main"], m)[0], "-")


if __name__ == "__main__":
    unittest.main()
