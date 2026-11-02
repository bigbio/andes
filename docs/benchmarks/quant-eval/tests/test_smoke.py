#!/usr/bin/env python3
"""Smoke test on synthetic data (a few seconds):
  1. the lineage guard refuses a cross-run join and accepts a same-run one;
  2. the quant_q -> p recovery inverts the running mean;
  3. scorecard fixes: same-family BH, TP at matched FDP, per-pair AUC, design-consistent detection, absence table, RT-local
     normalisation (unit tests on small synthetic inputs);
  4. the whole downstream (andes finalize parquet + MaxQuant evidence, design and amounts files) finds the changing species
     (sign-aware TP > 0), no false positive in the null species, and a small precursor error;
  5. false-transfer rate and its f_S scaling on a design with an absent species; blank carry-over with an `order` column (and
     the unchanged counts plus a note without it); the FlashLFQ QuantifiedPeaks reader (keys, transfers, decoys, MBR q cut,
     duplicate peaks); an absence-only downstream run with andes, MaxQuant and FlashLFQ arms through the scorecard.

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
from quant_eval import downstream, keys, metrics, readers, rollup, scorecard  # noqa: E402
from quant_eval.fasta import map_precursors  # noqa: E402
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
        self.assertTrue(any("blank-run rows" in r for r in rows))
        self.assertFalse(any("blank-run" in r for r in scorecard.Section(res).build()))   # not mixed into the main table


FLQ_HEADER = ("File Name\tBase Sequence\tFull Sequence\tProtein Group\tOrganism\tPeptide Monoisotopic Mass\tMS2 Retention Time\t"
              "Precursor Charge\tTheoretical MZ\tPeak intensity\tPeak RT Start\tPeak RT Apex\tPeak RT End\tPeak FWHM\tPeak FWHM Status\t"
              "Peak MZ\tPeak Charge\tNum Charge States Observed\tPeak Detection Type\tPIP Q-Value\tPIP PEP\tPSMs Mapped\t"
              "Base Sequences Mapped\tFull Sequences Mapped\tPeak Split Valley RT\tPeak Apex Mass Error (ppm)\tDecoy Peptide\tRandom RT")


def write_flashlfq(path, rows):
    """rows: (file, full sequence, charge, intensity, detection type, PIP q, decoy, random RT) -> QuantifiedPeaks.tsv (mzLib header)."""
    cols = FLQ_HEADER.split("\t")
    with open(path, "w") as f:
        f.write(FLQ_HEADER + "\n")
        for fn, seq, z, q, t, pq, dec, rrt in rows:
            v = dict.fromkeys(cols, "")
            v.update({"File Name": fn, "Full Sequence": seq, "Base Sequence": seq, "Precursor Charge": str(z), "Peak intensity": str(q),
                      "Peak Detection Type": t, "PIP Q-Value": pq, "Decoy Peptide": dec, "Random RT": rrt})
            f.write("\t".join(v[c] for c in cols) + "\n")


class Truth(unittest.TestCase):
    """False-transfer rate (absent species), blank carry-over (design `order`), FlashLFQ reader."""

    def test_false_transfer_rate_and_scaling(self):
        with tempfile.TemporaryDirectory() as tmp:
            with open(f"{tmp}/d.tsv", "w") as f:
                f.write("run\tabsent_species\nH1\tYEAST\nH2\tYEAST\nM1\t\nM2\t\n")
            design = Design(f"{tmp}/d.tsv")
        H = [f"H{i}" for i in range(8)]
        rows = [(k, "M1", False) for k in H + ["Y0", "Y1"]]           # 10 transfer-eligible precursors, 2 of them YEAST: f = 0.2
        rows += [(k, "H2", True) for k in H] + [("T0", "H1", True)]   # 9 HUMAN transfers in the runs without YEAST (T0: transfer only)
        rows += [("Y0", "H1", True), ("Y1", "H1", False), ("Y1", "M2", True)]   # 1 visible false transfer, 1 identified YEAST row
        raw = pd.DataFrame(rows, columns=["key", "run", "is_transfer"]).assign(quantity=1.0)
        sp = pd.Series({**{k: "HUMAN" for k in H + ["T0"]}, "Y0": "YEAST", "Y1": "YEAST"})
        x = metrics.false_transfer(raw, design, sp)["YEAST"]
        self.assertEqual((x["transferred_s_rows"], x["transferred_rows"], x["identified_s_rows"]), (1, 10, 1))
        self.assertEqual((x["eligible_precursors"], x["eligible_s_precursors"]), (10, 2))
        self.assertAlmostEqual(x["visible_rate"], 0.1)
        self.assertAlmostEqual(x["f_s"], 0.2)
        self.assertAlmostEqual(x["scaled_rate"], 0.5)
        o = metrics.absence(raw, design, sp)
        self.assertEqual(o["absent_species"]["YEAST"]["transferred"], 1)          # existing count unchanged
        self.assertEqual(o["false_transfer"]["YEAST"], x)
        none = metrics.false_transfer(raw[~raw.is_transfer], design, sp)["YEAST"]
        self.assertIsNone(none["visible_rate"])
        self.assertIsNone(none["scaled_rate"])

    def test_blank_carry_over_with_order(self):
        rows = [("a", "S1"), ("b", "S1"), ("c", "S1"), ("a", "B1"), ("b", "B1"), ("x", "B1"),
                ("d", "S2"), ("d", "B2"), ("a", "B2"), ("d", "B3"), ("y", "B3")]
        raw = pd.DataFrame(rows, columns=["key", "run"]).assign(quantity=1.0, is_transfer=lambda d: d.key.isin(["x", "y"]))
        sp = pd.Series(dict.fromkeys("abcdxy", "HUMAN"))
        with tempfile.TemporaryDirectory() as tmp:
            with open(f"{tmp}/d.tsv", "w") as f:
                f.write("run\tblank\torder\nS1\t\t1\nB1\ttrue\t2\nS2\t\t3\nB2\ttrue\t4\nB3\ttrue\t5\n")
            design = Design(f"{tmp}/d.tsv")
            with open(f"{tmp}/n.tsv", "w") as f:
                f.write("run\tblank\nS1\t\nB1\ttrue\nS2\t\nB2\ttrue\nB3\ttrue\n")
            plain = Design(f"{tmp}/n.tsv")
        self.assertEqual([design.run_before(r) for r in ("S1", "B1", "B3")], [None, "S1", "B2"])
        b = metrics.absence(raw, design, sp)["blank"]
        self.assertEqual((b["blank_rows"], b["carry_over"], b["blank_rows_excl_carry_over"]), (7, 4, 3))
        self.assertEqual(b["blank_transferred_excl_carry_over"], 2)
        self.assertEqual({k: (v["previous_run"], v["carry_over"]) for k, v in b["carry_over_by_blank"].items()},
                         {"B1": ("S1", 2), "B2": ("S2", 1), "B3": ("B2", 1)})   # B3 follows a blank: carry-over chains through B2
        p = metrics.absence(raw, plain, sp)["blank"]
        self.assertIsNone(p["carry_over"])
        self.assertIn("no 'order' column", p["carry_over_note"])
        self.assertEqual({k: p[k] for k in ("blank_rows", "blank_transferred", "blank_precursors")}, {k: b[k] for k in ("blank_rows", "blank_transferred", "blank_precursors")})
        res = {"arms": {"andes": {}}, "arm_roles": {"andes": {"kind": "andes", "rule": "default"}}, "primary": "andes", "classes": [],
               "absence": {"andes": {"blank": b}}}
        self.assertTrue(any("carry-over excluded" in r for r in scorecard.Section(res).absence_rows()))
        self.assertIn("| andes | 3 | 7 | 4 | 3 | 2 |", "\n".join(scorecard.carry_over_table(res)))
        res["absence"]["andes"]["blank"] = p
        self.assertIn("not computed", scorecard.carry_over_table(res)[0])

    def test_flashlfq_reader(self):
        acm = "[Common Biological:Acetylation on X]M[Common Variable:Oxidation on M]PEC[Common Fixed:Carbamidomethyl on C]K"
        k1 = keys.andes_key("[UNIMOD:1]M[UNIMOD:35]PEC[UNIMOD:4]K/2")
        self.assertEqual(keys.flashlfq_key(acm, 2), k1)
        self.assertEqual(keys.flashlfq_key("PEPTIDEK-[+15.9949]", 2), "PEPTIDEK[+15.99491]/2")   # C-terminal tag on the last residue
        with self.assertRaises(ValueError):
            keys.flashlfq_key("PEPS[Common Biological:Phosphorylation on S]IDEK", 2)
        rows = [("run-a", acm, 2, 1000, "MSMS", "", "False", "False"),
                ("Run_B", acm, 2, 500, "MBR", "0.001", "False", "False"),
                ("run-a", "PEPTIDEK", 2, 300, "MSMS", "", "False", "False"),
                ("run-a", "PEPTIDEK", 2, 900, "MSMS", "", "False", "False"),             # second peak of the same precursor x run
                ("Run_B", "PEPTIDEK", 2, 100, "MBR", "0.2", "False", "False"),           # above an MBR q limit of 0.05
                ("Run_B", "PEPTIDEK", 3, 100, "MBR", "0.001", "False", "True"),          # random-RT decoy peak
                ("run-a", "DECOYK", 2, 100, "MSMS", "", "True", "False"),                # decoy peptide
                ("run-a", "PEPS[Common Biological:Phosphorylation on S]IDEK", 2, 100, "MSMS", "", "False", "False"),
                ("run-a", "PEPTIDEK|PEPTLDEK", 2, 100, "MSMSAmbiguousPeakfinding", "", "False", "False"),
                ("run-a", "ZEROK", 2, 0, "MSMS", "", "False", "False"),
                ("Run_B", "ISOK", 2, 50, "MBR_IsoTrack", "", "False", "False"),          # transfer without a q value: kept
                ("other", "PEPTIDEK", 2, 100, "MSMS", "", "False", "False")]             # run not in the design
        with tempfile.TemporaryDirectory() as tmp:
            write_flashlfq(f"{tmp}/QuantifiedPeaks.tsv", rows)
            d = readers.read_flashlfq(f"{tmp}/QuantifiedPeaks.tsv", ["RUN_A", "RUN_B"])
            dq = readers.read_flashlfq(f"{tmp}/QuantifiedPeaks.tsv", ["RUN_A", "RUN_B"], mbr_q=0.05)
            with open(f"{tmp}/bad.tsv", "w") as f:
                f.write("File Name\tFull Sequence\tPrecursor Charge\tPeak intensity\n")
            with self.assertRaises(SystemExit):
                readers.read_flashlfq(f"{tmp}/bad.tsv", ["RUN_A"])
            with open(f"{tmp}/noq.tsv", "w") as f:
                f.write("File Name\tFull Sequence\tPrecursor Charge\tPeak intensity\tPeak Detection Type\nrun-a\tPEPTIDEK\t2\t5\tMSMS\n")
            self.assertEqual(len(readers.read_flashlfq(f"{tmp}/noq.tsv", ["RUN_A"])), 1)
            with self.assertRaises(SystemExit):
                readers.read_flashlfq(f"{tmp}/noq.tsv", ["RUN_A"], mbr_q=0.05)
        got = {(k, r): (q, t) for k, r, q, t in zip(d.key, d.run, d.quantity, d.is_transfer)}
        self.assertEqual(got, {(k1, "RUN_A"): (1000.0, False), (k1, "RUN_B"): (500.0, True), ("PEPTIDEK/2", "RUN_A"): (900.0, False),
                               ("PEPTIDEK/2", "RUN_B"): (100.0, True), ("ISOK/2", "RUN_B"): (50.0, True)})
        self.assertEqual(set(zip(dq.key, dq.run)), set(got) - {("PEPTIDEK/2", "RUN_B")})
        self.assertTrue(G.tag_of(d).startswith("flashlfq-quantifiedpeaks:mbr_q=None:"))
        self.assertNotEqual(G.tag_of(d), G.tag_of(dq))
        self.assertEqual(list(d.columns), ["key", "run", "quantity", "is_transfer", "library_q", "tq", "thr", "run_tag"])

    def test_downstream_absence_with_flashlfq(self):
        """Absence-only dataset: andes (default rule) + MaxQuant + FlashLFQ; YEAST absent from A runs, one blank after B_1."""
        with tempfile.TemporaryDirectory() as tmp:
            make_dataset(tmp)
            peps = sorted(set(k.split("/")[0] for k in pd.read_parquet(f"{tmp}/final.parquet").precursor_key))
            spm = map_precursors([p + "/2" for p in peps], f"{tmp}/db.fasta", ("HUMAN", "YEAST"))[0]
            Hs = [k[:-2] for k in spm.index if spm[k] == "HUMAN"]
            Ys = [k[:-2] for k in spm.index if spm[k] == "YEAST"]
            # (peptide, run, transfer, transfer_q): HUMAN 0-59 and YEAST 0-19 identified in B runs, HUMAN 0-59 in A runs;
            # transfers into A_1: YEAST 0-1 (q 0.001), YEAST 2-3 (q 0.5), HUMAN 60-67 (q 0.001, identified in B_1);
            # blank after B_1: HUMAN 0-4 (carry-over) and HUMAN 70 (not in B_1)
            ev = [(p, r, False, 0.0) for p in Hs[:60] for r in ("A_1", "A_2", "A_3", "B_1", "B_2")]
            ev += [(p, r, False, 0.0) for p in Ys[:20] for r in ("B_1", "B_2")] + [(p, "B_1", False, 0.0) for p in Hs[60:68]]
            ev += [(p, "A_1", True, 0.001) for p in Ys[:2] + Hs[60:68]] + [(p, "A_1", True, 0.5) for p in Ys[2:4]]
            ev += [(p, "BLANK", False, 0.0) for p in Hs[:5] + [Hs[70]]]
            pd.DataFrame({"precursor_key": [p + "/2" for p, *_ in ev], "run": [e[1] for e in ev],
                          "evidence_type": ["transfer_mbr" if e[2] else "identified_ms2" for e in ev], "quantity": 1e6,
                          "quant_q": [e[3] for e in ev], "transfer_q": [e[3] for e in ev], "library_q": 0.0}).to_parquet(f"{tmp}/abs.parquet")
            with open(f"{tmp}/evidence_abs.txt", "w") as f:
                f.write("Modified sequence\tCharge\tRaw file\tType\tIntensity\tReverse\tPotential contaminant\n")
                f.writelines(f"_{p}_\t2\t{r}\t{'MULTI-MATCH' if t else 'MULTI-MSMS'}\t1000000\t\t\n" for p, r, t, _ in ev)
            write_flashlfq(f"{tmp}/QuantifiedPeaks.tsv", [(r, p, 2, 1e6, "MBR" if t else "MSMS", ("0.2" if q > 0.1 else "0.001") if t else "",
                                                           "False", "False") for p, r, t, q in ev])
            with open(f"{tmp}/design_abs.tsv", "w") as f:
                f.write("run\tblank\tabsent_species\torder\nB_1\t\t\t1\nBLANK\ttrue\t\t2\nA_1\t\tYEAST\t3\nA_2\t\tYEAST\t4\nA_3\t\tYEAST\t5\nB_2\t\t\t6\n")
            out = f"{tmp}/abs.json"
            downstream.main(["--name", "absence", "--design", f"{tmp}/design_abs.tsv", "--fasta", f"{tmp}/db.fasta", "--classes", "HUMAN,YEAST",
                             "--andes", f"andes:default:{tmp}/abs.parquet", "--maxquant", f"MQ:{tmp}/evidence_abs.txt",
                             "--flashlfq", f"FLQ:{tmp}/QuantifiedPeaks.tsv", "--flashlfq-mbr-q", "0.05", "--out", out])
            with open(out) as f:
                res = json.load(f)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                scorecard.main(["--results", out])
        AB = res["absence"]
        self.assertEqual(res["arm_roles"]["FLQ"], {"kind": "flashlfq", "mbr_q": 0.05})
        f_s = 20 / 89                                       # eligible: HUMAN 0-67 and 70 (69), YEAST 0-19 (20)
        for arm, key, n_s, n in (("andes", "false_transfer", 2, 10), ("andes", "false_transfer_all", 4, 12),
                                 ("MQ", "false_transfer", 4, 12), ("FLQ", "false_transfer", 2, 10)):
            x = AB[arm][key]["YEAST"]
            self.assertEqual((x["transferred_s_rows"], x["transferred_rows"], x["identified_s_rows"]), (n_s, n, 0), (arm, key))
            self.assertAlmostEqual(x["f_s"], f_s, msg=arm)
            self.assertAlmostEqual(x["scaled_rate"], n_s / n / f_s, msg=arm)
        for arm in ("andes", "MQ", "FLQ"):
            b = AB[arm]["blank"]
            self.assertEqual((b["blank_rows"], b["carry_over"], b["blank_rows_excl_carry_over"]), (6, 5, 1), arm)
            self.assertEqual(b["carry_over_by_blank"]["BLANK"]["previous_run"], "B_1")
        md = buf.getvalue()
        for text in ("False-transfer rate per arm", "| FLQ | MBR q <= 0.05 | YEAST | 3 | 2 | 10 | 0 | 20.00 | 22.47 | 89.00 |",
                     "| andes | rule all (every finalize transfer) | YEAST | 3 | 4 | 12 |", "carry-over excluded [andes]",
                     "false-transfer rate, scaled % (visible / f_YEAST) [andes]", "| FLQ | 1 | 6 | 5 | 1 | 0 | BLANK after B_1 |"):
            self.assertIn(text, md)


class AndesColumns(unittest.TestCase):
    """The rules that read andes finalize schema 1.1 columns."""

    def test_fixed_and_no_shared_rules(self):
        from quant_eval.readers import read_andes
        with tempfile.TemporaryDirectory() as tmp:
            runs = ["A_1", "A_2"]
            rows = pd.DataFrame({
                "precursor_key": ["AAAK/2", "AAAK/2", "CCCK/2", "CCCK/2", "DDDK/2"],
                "run": ["A_1", "A_2", "A_1", "A_2", "A_2"],
                "evidence_type": ["identified", "transfer_ms1", "identified", "transfer_ms1", "transfer_ms1"],
                "quantity": [1e6, 2e6, 1e6, 3e6, 4e6],
                "quant_q": [0.0, 0.001, 0.0, 0.002, 0.003],
                "transfer_q": [np.nan, 0.001, np.nan, 0.002, 0.003],
                "library_q": 0.0,
                "shared_peak_with": [None, None, None, "AAAK/2", None],
                "quantity_fixed": [np.nan, 1.5e6, np.nan, -2e5, 3.9e6],
            })
            p = f"{tmp}/final.parquet"
            rows.to_parquet(p)
            q = lambda t: {(k, r): v for k, r, v in zip(t.key, t.run, t.quantity)}
            d = q(read_andes(p, runs, "default"))
            f = q(read_andes(p, runs, "fixed"))
            n = q(read_andes(p, runs, "no-shared"))
            assert len(d) == 5
            # Identified rows keep their quantity; a non-positive fixed value is missing.
            assert len(f) == 4 and min(f.values()) > 0
            assert sorted(f.values()) == [1e6, 1e6, 1.5e6, 3.9e6]
            # The shared transfer is left out.
            assert len(n) == 4 and 3e6 not in n.values()
            # An older table without the columns is refused with a message.
            rows.drop(columns=["shared_peak_with", "quantity_fixed"]).to_parquet(p)
            with self.assertRaises(SystemExit):
                read_andes(p, runs, "fixed")


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
