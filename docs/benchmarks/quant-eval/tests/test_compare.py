#!/usr/bin/env python3
"""Tests of the short comparison (compare.py) and of the suite's worker processes on synthetic results (seconds):

    python3 tests/test_compare.py
"""
import contextlib
import io
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from quant_eval import compare as C, suite  # noqa: E402


def result(seed, tp, fp):
    def arm(prec, tp_, fp_):
        return {"precursors": prec, "proteins_ge2": prec // 10, "de_family": {"common": {"TP": tp_, "FP": fp_}}}
    return {"suite_dataset": "d", "seed": seed, "dataset": "SET", "builds": ["main", "cand"], "levels": [],
            "arm_roles": {"main": {"kind": "andes"}, "cand": {"kind": "andes"}, "MQ": {"kind": "maxquant"}, "MQ_noMBR": {"kind": "maxquant"},
                          "MQ_q01": {"kind": "maxquant_q"}, "FP": {"kind": "fragpipe"}},
            "arms": {"main": arm(900, 1, 1), "cand": arm(1000, tp, fp), "MQ": arm(800, 100, 20), "MQ_noMBR": arm(700, 1, 1), "MQ_q01": arm(700, 1, 1),
                     "FP": arm(1000, 100, 50)},
            "cell_accuracy": {"arms": {"cand": {"bad": 2.0}, "MQ": {"bad": 1.0}, "FP": {"bad": 2.0}}}}


class Compare(unittest.TestCase):
    def test_marker_and_difference(self):
        self.assertEqual(C.marker(10, (10, 10), 8, (8, 8), True), "✅")
        self.assertEqual(C.marker(10, (10, 10), 8, (8, 8), False), "❌")
        self.assertEqual(C.marker(10, (10, 10), 10, (10, 10), True), "≈")
        self.assertEqual(C.marker(10, (9, 12), 11, (11, 11), False), "≈")      # the tool lies inside the seed range of andes
        self.assertEqual(C.marker(10, (9, 10), 11, (11, 12), False), "✅")
        self.assertEqual(C.marker(None, None, 1, None, True), "-")
        self.assertEqual(C.difference(110, 100), "+10.0 %")
        self.assertEqual(C.difference(45, 100), "-55.0 %")
        self.assertEqual(C.difference(330, 100), "3.3x")
        self.assertEqual(C.difference(3, 0), "+3")
        self.assertEqual(C.difference(None, 1), "-")

    def test_table(self):
        S = {"42": result("42", 110, 150), "1": result("1", 120, 180)}
        rows = {ln.split("|")[1].strip(): [c.strip() for c in ln.split("|")[2:-1]] for ln in C.table(S, "42").split("\n") if ln.startswith("| ")}
        # default arm: the last build; default tools: the MaxQuant and FragPipe arms with match between runs
        self.assertEqual(rows["metric"], ["andes cand", "MaxQuant (MBR)", "vs MaxQuant (MBR)", "FragPipe", "vs FragPipe"])
        self.assertEqual(rows["precursors"], ["1,000", "800", "✅ +25.0 %", "1,000", "≈ +0.0 %"])
        self.assertEqual(rows["cells off by > 2x, % (shared cells)"], ["2.00", "1.00", "❌ 2.0x", "2.00", "≈ +0.0 %"])
        self.assertEqual(rows["DE true calls, common family"], ["110 (110 to 120)", "100", "✅ +10.0 %", "100", "✅ +10.0 %"])
        self.assertEqual(rows["DE false calls, common family"], ["150 (150 to 180)", "20", "❌ 7.5x", "50", "❌ 3.0x"])
        self.assertNotIn("protein CV, %", rows)                                   # a metric no arm has is not printed
        one = [ln for ln in C.table(S, "42", build="main", tools=["FP"]).split("\n") if ln.startswith("| metric")]
        self.assertEqual([c.strip() for c in one[0].split("|")[1:-1]], ["metric", "andes main", "FragPipe", "vs FragPipe"])
        self.assertIn("no arm nope", C.table(S, "42", build="nope"))

    def test_main_reads_a_results_folder(self):
        with tempfile.TemporaryDirectory() as d:
            for s in ("42", "1"):
                json.dump(result(s, 110, 150), open(os.path.join(d, f"d_s{s}.json"), "w"))
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                C.main(["--results", d])
            self.assertIn("### SET: andes cand against MaxQuant (MBR), FragPipe", buf.getvalue())
            with self.assertRaises(SystemExit):
                C.main(["--results", os.path.join(d, "empty")])


class SuiteJobs(unittest.TestCase):
    def test_units_run_once_each_with_and_without_workers(self):
        """--jobs changes how the (dataset, seed) units are run, not which ones: the same units, one result file each."""
        M = {"reference_build": "main", "builds": ["main"],
             "datasets": {k: {"design": "x", "fasta": "x", "classes": "A,B", "tools": {"sage": {"S": "x"}}} for k in ("a", "b")}}
        seen = []
        orig = suite._unit
        try:
            for jobs in ("1", "3"):
                with tempfile.TemporaryDirectory() as d:
                    json.dump(M, open(os.path.join(d, "m.json"), "w"))
                    suite._unit = _touch
                    suite.main(["--manifest", os.path.join(d, "m.json"), "--out", os.path.join(d, "o"), "--seeds", "42,1", "--jobs", jobs])
                    seen.append(sorted(os.listdir(os.path.join(d, "o"))))
        finally:
            suite._unit = orig
        self.assertEqual(seen[0], ["a_s1", "a_s42", "b_s1", "b_s42"])
        self.assertEqual(seen[0], seen[1])


def _touch(u):
    _, d, seed, out = u[:4]
    open(os.path.join(out, f"{d}_s{seed}"), "x").close()


if __name__ == "__main__":
    unittest.main()
