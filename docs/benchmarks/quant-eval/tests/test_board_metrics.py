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
from quant_eval import board, de as DE, lowlevel as LL, metrics  # noqa: E402
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
        self.assertEqual(b["l1_pooled"]["rows"], 0)                      # no UPS rows: zero counts, as in the legacy L1
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


class DEFamilies(unittest.TestCase):
    """Sign-aware counts and the common-tested family on a small protein table: conditions A (A1..A3) and B (B1..B3)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        with open(f"{self.tmp.name}/d.tsv", "w") as f:
            f.write("run\tcondition\n" + "".join(f"{r}\t{r[0]}\n" for r in ["A1", "A2", "A3", "B1", "B2", "B3"]))
        with open(f"{self.tmp.name}/a.tsv", "w") as f:
            f.write("species\tcondition\tamount\nUPS\tA\t1\nUPS\tB\t2\nYEAST\t*\t1\n")
        self.design = Design(f"{self.tmp.name}/d.tsv", f"{self.tmp.name}/a.tsv", "YEAST", ("UPS", "YEAST"))
        rng = np.random.default_rng(0)
        runs = ["A1", "A2", "A3", "B1", "B2", "B3"]
        rows, sp = {}, {}
        for i in range(40):                 # UPS expected +1 log2: 30 with the right sign, 10 with the opposite sign (-2)
            fc = 1.0 if i < 30 else -2.0
            rows[f"U{i}"] = list(10 + rng.normal(0, 0.05, 3)) + list(10 + fc + rng.normal(0, 0.05, 3))
            sp[f"U{i}"] = "UPS"
        for i in range(60):                 # null species, one of them strongly changed
            fc = 3.0 if i == 0 else 0.0
            rows[f"Y{i}"] = list(12 + rng.normal(0, 0.05, 3)) + list(12 + fc + rng.normal(0, 0.05, 3))
            sp[f"Y{i}"] = "YEAST"
        self.P = pd.DataFrame(rows, index=runs).T
        self.sp = sp

    def tearDown(self):
        self.tmp.cleanup()

    def test_wrong_sign_calls_are_counted_apart_and_as_false(self):
        tot, _, _, fam = DE.de_all(self.P, self.design, self.sp, ("UPS", "YEAST"))
        self.assertEqual(tot["TP"], tot["TPs"] + tot["WS"])
        self.assertEqual(tot["WS"], 10)
        self.assertEqual(tot["TPs"], 30)
        self.assertGreaterEqual(tot["FP"], 1)                            # the strongly changed null protein (and perhaps a chance call)
        # at pooled FDP 50 % the unsigned reading keeps the wrong-sign rows as true positives, the signed one counts them as false
        tp_u, fp_u = DE.tp_at_fdp([DE.de_table(self.P, ["A1", "A2", "A3"], ["B1", "B2", "B3"], self.sp, lambda s: 1.0 if s == "UPS" else 0.0)], "YEAST", (0.5,))
        tp_s, fp_s = DE.tp_at_fdp([DE.de_table(self.P, ["A1", "A2", "A3"], ["B1", "B2", "B3"], self.sp, lambda s: 1.0 if s == "UPS" else 0.0)], "YEAST", (0.5,), signed=True)
        self.assertGreater(tp_u["50"], tp_s["50"])
        self.assertIn("tp_at_fdp_signed", fam["own"])

    def test_common_family_has_identical_tested_rows(self):
        P2 = self.P.copy()
        P2.loc[[f"U{i}" for i in range(10)], ["B1", "B2"]] = np.nan      # arm 2 cannot test 10 proteins (one finite value on one side)
        P2.loc[[f"U{i}" for i in range(10)], ["B3"]] = 11.0
        P3 = self.P.copy()
        com = DE.common_sets([self.P, P2, P3], self.design)
        self.assertEqual(len(com[("A", "B")]), 100 - 10)
        for P in (self.P, P2, P3):
            tot, _, _, fam = DE.de_all(P, self.design, self.sp, ("UPS", "YEAST"), fam_common=com)
            self.assertEqual(fam["common"]["tested"], 90)                # the same rows for every arm
        _, _, _, own = DE.de_all(self.P, self.design, self.sp, ("UPS", "YEAST"), fam_common=com)
        self.assertEqual(own["own"]["tested"], 100)                      # the arm's own testable rows are more
        self.assertEqual(own["common"]["tested_null"], 60)

    def test_board_columns_for_de(self):
        _, _, _, fam = DE.de_all(self.P, self.design, self.sp, ("UPS", "YEAST"), fam_inter=set(self.P.index),
                                 fam_common=DE.common_sets([self.P], self.design))
        r = {"arms": {"main": {"de_family": fam}}}
        c = board.new_cells(r, "main", ["main"], None)
        self.assertEqual(len(c), board.N_NEW)
        fp = fam["inter"]["FP"]
        self.assertEqual(c[9], f"30 / {fp + 10} (10)")                   # signed TP / FP + wrong-sign (wrong-sign)
        self.assertEqual(c[10], f"40 / {fp}")
        self.assertEqual(c[12], "100 / 100 / 100")


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
        self.assertEqual(len(cells), board.N_NEW)
        self.assertEqual(cells[5], f"{m['main'][1]:.1f}")                # headline: random rows
        self.assertEqual(cells[6], f"{m['main'][0]:,}")                  # best first
        self.assertEqual(board.new_cells(r, "none", ["main"], m)[0], "-")


class Reviewed(unittest.TestCase):
    """Fixes of the independent review of the board metrics."""

    def test_binned_cutoffs_change_the_classification(self):
        rng = np.random.default_rng(1)
        n_abs = np.r_[np.full(300, 1), np.full(300, 12)]
        d = np.r_[rng.normal(0, 1.5, 300), rng.normal(0, 0.15, 300)]     # precursors seen in one absent run are noisy, in twelve are not
        cuts = AB.learn_cutoffs(d, n_abs)
        self.assertEqual([c[2] for c in cuts], [300, 0, 0, 300])           # n control per bin is reported
        self.assertEqual([c[3] for c in cuts], [False, True, True, False])  # empty bins fall back to the pooled cutoff
        self.assertGreater(AB.cutoff_for(cuts, 1), 2.0)
        self.assertLess(AB.cutoff_for(cuts, 12), 0.6)
        ps = pd.DataFrame({"key": ["a", "b"], "d": [1.5, 1.5], "n_abs": [1, 12], "n_pres": [20, 20]})
        cl = AB.classify(ps, cuts)
        self.assertEqual(list(cl), ["not_specific", "yeast_enriched"])    # the same d is enriched with twelve absent rows, not with one
        self.assertEqual(list(AB.classify(ps, None, 1.0)), ["yeast_enriched", "yeast_enriched"])

    def test_reference_classes_are_the_same_for_every_arm(self):
        with tempfile.TemporaryDirectory() as tmp:
            design, raw, sp = Absence().make(tmp)
            refs = AB.reference_profiles(raw, design, sp)
            sub = raw[~raw.key.isin(["Y3", "Y45"])]                          # an arm that does not report two precursors
            o_full = AB.absence_classes(raw, design, sp, refs=refs)["YEAST"]
            o_sub = AB.absence_classes(sub, design, sp, refs=refs)["YEAST"]
            o_own = AB.absence_classes(sub, design, sp)["YEAST"]
            extra = pd.concat([raw, pd.DataFrame({"key": ["Q1"], "run": ["H1"], "quantity": [2.0 ** 17], "is_transfer": [True]})])
            o_x = AB.absence_classes(extra, design, pd.concat([sp, pd.Series({"Q1": "YEAST"})]), refs=refs)["YEAST"]
        self.assertIn("own_learned", o_sub)
        self.assertEqual(o_sub["reference"]["cutoffs"], o_full["reference"]["cutoffs"])      # cutoffs derived once, not per arm
        enr_full = o_full["learned"]["identified"]["yeast_enriched"] + o_full["learned"]["transferred"]["yeast_enriched"]
        enr_sub = o_sub["learned"]["identified"]["yeast_enriched"] + o_sub["learned"]["transferred"]["yeast_enriched"]
        self.assertEqual(enr_full - enr_sub, 2)                              # Y3 lost its two rows; its class did not change for the others
        self.assertEqual(o_x["learned"]["transferred"]["only_absent"], o_full["learned"]["transferred"]["only_absent"])
        self.assertEqual(o_x["learned"]["transferred"]["not_in_reference"], 1)  # a precursor absent from the reference table
        self.assertNotIn("reference", o_own)                                  # without refs the classes come from the arm's own rows
        self.assertGreater(sum(x["n_control"] for x in o_own["cutoffs"]), 0)

    def test_unscored_rows_and_survivorship(self):
        with tempfile.TemporaryDirectory() as tmp:
            design = ups_design(tmp)
            raw, sp = ups_raw()
            junk = [(f"J{k}", r, 2.0 ** 10, True, 0.5) for k in range(50) for r in RUNS[:4]]   # only at the two lowest levels: no reference row
            raw2 = pd.concat([raw, pd.DataFrame(junk, columns=raw.columns)])
            sp2 = pd.concat([sp, pd.Series({f"J{k}": "UPS" for k in range(50)})])
            a, b = LL.l1_detail(raw, design, sp), LL.l1_detail(raw2, design, sp2)
        self.assertEqual(a["l1_pooled"]["unscored_rows"], 0)
        self.assertEqual(b["l1_pooled"]["unscored_rows"], 200)
        self.assertEqual(b["l1_pooled"]["rows"], a["l1_pooled"]["rows"])      # junk leaves the scored rows: the L1 count does not see it
        self.assertEqual(b["per_level"]["1"]["keys_within2x"], a["per_level"]["1"]["keys_within2x"])
        self.assertLess(b["l1_pooled"]["within2x_incl_unscored"], a["l1_pooled"]["within2x_incl_unscored"])
        self.assertAlmostEqual(b["l1_pooled"]["within2x_incl_unscored"], b["l1_pooled"]["within2x"] / (b["l1_pooled"]["rows"] + 200))
        cells = board.new_cells({"arms": {"main": {"l1_detail": b}}}, "main", ["main"], None)
        self.assertEqual(cells[1], f"{b['l1_pooled']['rows']:,} (+200)")

    def test_empty_ups_set_keeps_the_legacy_zero_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            design = ups_design(tmp)
            raw, sp = ups_raw()
            out = metrics.design_detection(raw[raw.key.str.startswith("Y")], design, sp)
        self.assertEqual(out["1"]["t1.0"], 0)
        self.assertEqual(sum(out[L]["t1.0"] for L in design.levels[:3]), 0)

    def test_blank_reference_needs_two_other_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            design, raw, sp = Absence().make(tmp)
            one = raw[raw.key == "H0"].iloc[:1]
            only1 = pd.DataFrame({"key": ["ONE"], "run": ["M1"], "quantity": [2.0 ** 20], "is_transfer": [False]})
            blank = pd.DataFrame({"key": ["H0", "ONE"], "run": "BLANK", "quantity": [2.0 ** 10, 2.0 ** 10], "is_transfer": [False, False]})
            o = AB.blank_split(pd.concat([raw, only1, blank]), design)
        self.assertEqual(o["identified"]["no_reference"], 1)               # ONE is in one other run only
        self.assertEqual(o["identified"]["carry_over_fixed"], 1)

    def test_de_counts_without_an_expected_sign(self):
        T = pd.DataFrame({"sp": ["UPS", "UPS", "UPS", "UPS", "YEAST"], "fc": [1.0, -1.0, 1.0, 1.0, 2.0], "q": 0.001,
                          "ex": [1.0, 1.0, 0.0, np.nan, 0.0], "p": 0.001})
        c = DE._count(T, "YEAST")
        self.assertEqual((c["TP"], c["TPs"], c["WS"], c["NS"], c["FP"]), (4, 1, 1, 2, 1))
        self.assertEqual(c["TP"], c["TPs"] + c["WS"] + c["NS"])
        tp, fp = DE.tp_at_fdp([T], "YEAST", (1.0,), signed=True)
        self.assertEqual((tp["100"], fp["100"]), (3, 2))                   # NaN and zero expected sign are not wrong-sign; the wrong-sign row and the null row are false

    def test_board_absence_and_blank_cells(self):
        with tempfile.TemporaryDirectory() as tmp:
            design, raw, sp = Absence().make(tmp)
            refs = AB.reference_profiles(raw, design, sp)
            ac = AB.absence_classes(raw, design, sp, refs=refs)
            bl = AB.blank_split(pd.concat([raw, pd.DataFrame({"key": ["H0", "H1", "H2"], "run": "BLANK", "quantity": [2.0 ** 20 * 0.001] * 3, "is_transfer": [True, False, False]})]), design)
        r = {"arms": {"main": {}}, "absence_classes": {"main": ac}, "blank_split": {"main": bl}}
        c = board.new_cells(r, "main", ["main"], None)
        o = ac["YEAST"]
        den_t, den_i = o["denominators"]["transferred_rows_absent_runs"], o["denominators"]["identified_rows_absent_runs"]
        e_t = o["learned"]["transferred"]["yeast_enriched"]
        a_t = o["learned"]["transferred"]["only_absent"] + o["learned"]["transferred"]["not_in_reference"]
        self.assertEqual(c[7].split(" / ")[0], f"{100 * e_t / den_t:.3f} ({100 * (e_t + a_t) / den_t:.3f})")   # both readings
        self.assertTrue(c[7].endswith("%"))
        self.assertEqual(c[8], f"2 / 1 / {bl['carry_over_learned']}")      # two identified rows (H1, H2), one transferred (H0), all three carry-over level
        self.assertEqual(board.rate_pct(1, 0, "{:.3f}"), "-")
        self.assertEqual(board.rate_pct(0, 10, "{:.3f}"), "0.000")

    def test_legacy_board_is_byte_identical_to_the_old_board(self):
        import importlib.util
        import subprocess
        here = os.path.dirname(os.path.abspath(__file__))
        try:
            top = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=here, stderr=subprocess.DEVNULL).decode().strip()
            src = subprocess.check_output(["git", "show", "6c2b72c3:docs/benchmarks/quant-eval/quant_eval/board.py"], cwd=top, stderr=subprocess.DEVNULL).decode()
        except Exception:
            self.skipTest("the base commit is not available")
        with tempfile.TemporaryDirectory() as tmp:
            f = os.path.join(tmp, "old_board.py")
            with open(f, "w") as fh:
                fh.write(src)
            spec = importlib.util.spec_from_file_location("old_board", f)
            old = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(old)

        def arm(l1, tr=100):
            return {"precursors": 1000, "proteins_ge2": 100, "ident_rows": 5000, "transfer_rows": tr,
                    "ups_design": {L: {"t1.0": v} for L, v in zip(["1", "2", "4", "8"], l1)},
                    "de_family": {"inter": {"TP": 50, "FP": 5, "TP5": 40, "FP5": 2, "tp_at_fdp": {"10": 55}}}}
        res = {"u": {"42": {"dataset": "SET", "suite_dataset": "u", "seed": "42", "builds": ["main"], "levels": ["1", "2", "4", "8"], "arm_roles": {"main": {"kind": "andes"}, "MQ": {"kind": "maxquant"}, "FL001": {"kind": "flashlfq"}},
                            "arms": {"main": arm([3, 4, 5, 6]), "MQ": arm([0, 0, 0, 0], 0), "FL001": arm([1, 1, 1, 1])}, "cv": {}, "cell_accuracy": {"arms": {"main": {"bad": 2.5, "within2x": 97.5, "bad_bright": 3.0, "rep_sd_top": 0.1}}}}},
               "p": {"42": {"dataset": "P14", "suite_dataset": "p", "seed": "42", "builds": ["main"], "arm_roles": {"main": {"kind": "andes"}},
                            "arms": {"main": {"precursors": 10, "proteins_ge2": 5, "ident_rows": 20, "transfer_rows": 0}}, "rates": {"main": {"YEAST": {"rate": float("nan")}}}}},
               "s": {"42": {"dataset": "SC", "suite_dataset": "s", "seed": "42", "builds": ["main"], "arm_roles": {"main": {"kind": "andes"}},
                            "arms": {"main": {"precursors": 10, "proteins_ge2": 5, "ident_rows": 20, "transfer_rows": 3}}, "carry_over": {"main": {"rows": 23}}}}}
        order = ["u", "p", "s"]
        self.assertEqual(old.build_rows(res, "42", ["main"], order), board.build_rows(res, "42", ["main"], order, legacy=True))
        self.assertIn("nan %", board.build_rows(res, "42", ["main"], order, legacy=True))
        new = board.build_rows(res, "42", ["main"], order, legacy=False)
        self.assertNotIn("nan %", new)
        self.assertNotIn("within 2x %", new.split("\n")[0].replace("L1 within-2x %", ""))


if __name__ == "__main__":
    unittest.main()
