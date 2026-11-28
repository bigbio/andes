#!/usr/bin/env python3
"""Unit tests of the cross-link quantification harness on small synthetic inputs (seconds):
  1. residue-pair keys: protein positions, accession reduction, canonical order;
  2. the long-table contract refuses duplicates and unknown evidence;
  3. the Skyline readers: a pseudo-peptide document and both CV report formats;
  4. the metrics: dilution accuracy and slope, replicate CV and completeness, null calls,
     false-transfer rate;
  5. arm filters and run matching;
  6. the whole evaluate command on a synthetic dilution design.

    python3 tests/test_xlq.py
"""
import json
import os
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from xlq import design, evaluate, keys, metrics, skyline, table  # noqa: E402

PROTEIN = "MAKDEFGKHIKLMNPKQR"

SKY = """<?xml version="1.0" encoding="utf-8"?>
<srm_settings format_version="4.2">
 <peptide_list>
  <peptide sequence="AKDEFGKKHIKLMNPK" modified_sequence="AK[+1]DEFGKK[+26]HIK[+1]LMNPK">
   <variable_modifications>
    <variable_modification index_aa="1" modification_name="LinkSite(K,S,T,Y,D)" />
    <variable_modification index_aa="7" modification_name="Kxlink" />
    <variable_modification index_aa="10" modification_name="LinkSite(K,S,T,Y,D)" />
   </variable_modifications>
   <precursor charge="4" modified_sequence="AK[+1.0]DEFGKK[+26.0]HIK[+1.0]LMNPK">
    <precursor_results>
     <precursor_peak replicate="r1" area="1000" identified="true" isotope_dotp="0.95" />
     <precursor_peak replicate="r2" area="500" identified="aligned" isotope_dotp="0.80" />
     <precursor_peak replicate="r3" identified="false" />
    </precursor_results>
   </precursor>
  </peptide>
  <peptide sequence="DEFGKHIK" modified_sequence="DEFGKHIK">
   <precursor charge="2" modified_sequence="DEFGKHIK">
    <precursor_results><precursor_peak replicate="r1" area="7" identified="true" /></precursor_results>
   </precursor>
  </peptide>
 </peptide_list>
</srm_settings>
"""


def long_table(values, runs, rp="P1:3 -- P1:11", precursor="x/4", evidence="identified"):
    return pd.DataFrame([{"run": r, "precursor": precursor, "residue_pair": rp, "quantity": v, "evidence": evidence}
                         for r, v in zip(runs, values) if v is not None])


class Keys(unittest.TestCase):
    def test_side_and_pair(self):
        prot = {"P1": "AAKBBKAAK", "P2": "CCKAAK"}
        self.assertEqual(keys.side(prot, "AAK", 2), "P1:3;P1:9;P2:6")
        self.assertEqual(keys.pair_key("P2:5", "P1:3"), "P1:3 -- P2:5")

    def test_normalise(self):
        self.assertEqual(keys.normalise_key("sp|Q9|B_HUMAN:5 -- sp|P1|A_HUMAN:9;sp|P1|A_HUMAN:3"),
                         "P1:3;P1:9 -- Q9:5")
        self.assertEqual(keys.accessions("P1:3;Q2:4"), {"P1", "Q2"})
        with self.assertRaises(ValueError):
            keys.normalise_key("P1:3")


class Table(unittest.TestCase):
    def test_contract(self):
        df = long_table([1.0, 2.0], ["a", "b"])
        self.assertEqual(len(table.validate(df)), 2)
        with self.assertRaises(ValueError):
            table.validate(pd.concat([df, df]))
        with self.assertRaises(ValueError):
            table.validate(df.assign(evidence="mbr"))
        self.assertEqual(len(table.validate(df.assign(quantity=[0.0, 3.0]))), 1)


class Skyline(unittest.TestCase):
    def test_document(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "doc.sky"
            path.write_text(SKY, encoding="utf-8")
            df, counts = skyline.read_document(path, {"P1": PROTEIN}, "Kxlink", "LinkSite")
        self.assertEqual(counts["not_crosslink"], 1)
        self.assertEqual(counts["no_area"], 1)
        self.assertEqual(list(df["run"]), ["r1", "r2"])
        self.assertEqual(set(df["residue_pair"]), {keys.pair_key("P1:3", "P1:11")})
        self.assertEqual(list(df["evidence"]), ["identified", "transferred"])
        row = df.iloc[0]
        self.assertEqual((row["alpha_peptide"], row["alpha_site"], row["beta_peptide"], row["beta_site"]),
                         ("AKDEFGK", 2, "HIKLMNPK", 3))

    def test_cv_report(self):
        text = ("File Name,Peptide Modified Sequence,Product Charge,Fragment Ion,Cv Area\n"
                "a.raw,PEPK,3,precursor,11.80%\nb.raw,PEPK,3,precursor,11.80%\n"
                "a.raw,PEPR,4,precursor,0.25\na.raw,PEPS,4,precursor,\n")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cv.csv"
            path.write_text(text, encoding="utf-8")
            cv = skyline.read_cv_report(path).set_index("precursor")["cv"]
        self.assertAlmostEqual(cv["PEPK/3"], 0.118)
        self.assertAlmostEqual(cv["PEPR/4"], 0.25)
        self.assertNotIn("PEPS/4", cv.index)


class Metrics(unittest.TestCase):
    runs = ["l1a", "l1b", "l1c", "h1a", "h1b", "h1c"]
    level_of = {"l1a": 0.25, "l1b": 0.25, "l1c": 0.25, "h1a": 1.0, "h1b": 1.0, "h1c": 1.0}

    def dilution_matrix(self, compression):
        rng = np.random.default_rng(7)
        rows = {}
        for i in range(40):
            base = 20 + i * 0.1
            rows[f"P1:{i} -- P1:{i + 100}"] = [base + compression * np.log2(self.level_of[r]) + rng.normal(0, 0.05)
                                               for r in self.runs]
        return pd.DataFrame.from_dict(rows, orient="index", columns=self.runs)

    def test_dilution(self):
        d = metrics.dilution(self.dilution_matrix(1.0), self.level_of, 1.0)
        lvl = d["levels"]["0.25"]
        self.assertEqual(d["residue_pairs_at_every_level"], 40)
        self.assertLess(abs(lvl["median_error"]), 0.05)
        self.assertEqual(lvl["within_2x"], 1.0)
        self.assertEqual(lvl["de_called_expected_sign"], 40)
        self.assertIsNone(d["slope"])
        half = metrics.dilution(self.dilution_matrix(0.5), self.level_of, 1.0)["levels"]["0.25"]
        self.assertAlmostEqual(half["median_error"], 1.0, delta=0.05)

    def test_slope_two_levels(self):
        level_of = dict(self.level_of, l1c=0.5)
        runs = self.runs + ["m1a"]
        level_of["m1a"] = 0.5
        M = pd.DataFrame({r: [20 + 0.5 * np.log2(level_of[r]) + 0.01 * k for k in range(10)] for r in runs},
                         index=[f"P1:{k} -- P1:{k + 50}" for k in range(10)])
        self.assertAlmostEqual(metrics.dilution(M, level_of, 1.0)["slope"], 0.5, places=3)

    def test_precision(self):
        df = pd.concat([long_table([100, 120, 80], ["a", "b", "c"], precursor="x/3"),
                        long_table([5, None, None], ["a", "b", "c"], rp="P1:1 -- P1:2", precursor="y/3")])
        M = metrics.log_matrix(df, ["a", "b", "c"])
        p = metrics.precision(M, {"g": ["a", "b", "c"]})["g"]
        self.assertEqual((p["quantified"], p["complete"], p["completeness"]), (2, 1, 0.5))
        self.assertAlmostEqual(p["median_cv_complete"], 0.2, places=6)

    def test_rollups(self):
        df = pd.concat([long_table([100, 120, 80], ["a", "b", "c"], precursor="x/3"),
                        long_table([200, 240, None], ["a", "b", "c"], precursor="x/4")])
        S = metrics.log_matrix(df, ["a", "b", "c"], rollup="sum")
        self.assertAlmostEqual(2 ** S.loc["P1:3 -- P1:11", "a"], 300)
        self.assertAlmostEqual(2 ** S.loc["P1:3 -- P1:11", "c"], 80)
        P = metrics.log_matrix(df, ["a", "b", "c"])
        self.assertTrue(np.isfinite(P.to_numpy()).all())
        self.assertGreater(P.loc["P1:3 -- P1:11", "b"], P.loc["P1:3 -- P1:11", "a"])
        with self.assertRaises(ValueError):
            metrics.log_matrix(df, ["a"], rollup="mean")

    def test_null_calls(self):
        rng = np.random.default_rng(1)
        runs = [f"r{i}" for i in range(10)]
        M = pd.DataFrame(20 + rng.normal(0, 0.3, (200, 10)), columns=runs, index=[f"k{i}" for i in range(200)])
        out = metrics.null_calls(M, {"rep": runs})["rep"]
        self.assertEqual(out["tested"], 200)
        self.assertEqual(out["false_calls"], 0)
        shifted = M.copy()
        shifted[runs[0::2]] += 2.0
        self.assertEqual(metrics.null_calls(shifted, {"rep": runs})["rep"]["false_calls"], 200)

    def test_transfers(self):
        df = pd.DataFrame({
            "run": ["lib", "lib", "lib", "lib", "all"],
            "residue_pair": ["A:1 -- A:2", "B:1 -- B:2", "B:3 -- B:4", "B:1 -- B:2", "B:1 -- B:2"],
            "evidence": ["transferred", "transferred", "transferred", "identified", "transferred"],
            "quant_q": [0.001, 0.02, 0.001, np.nan, 0.001],
        })
        absent = lambda r: None if r["run"] == "all" else r["residue_pair"].startswith("B")  # noqa: E731
        t = metrics.transfers(df, absent)
        self.assertEqual((t["transferred"], t["false_transfers"], t["identified_absent"]), (3, 2, 1))
        self.assertAlmostEqual(t["false_transfer_rate"], 0.6667, places=4)
        self.assertEqual(t["at_quant_q"]["0.01"], {"transferred": 2, "false_transfer_rate": 0.5})


class Arms(unittest.TestCase):
    def test_filters(self):
        base, f = evaluate.parse_kind("skyline?min_isotope_dotp=0.9&evidence=identified")
        self.assertEqual((base, f), ("skyline", {"min_isotope_dotp": "0.9", "evidence": "identified"}))
        df = pd.DataFrame({"evidence": ["identified", "transferred", "transferred"],
                           "quant_q": [np.nan, 0.001, 0.2]})
        self.assertEqual(len(evaluate.apply_filters(df, {"max_quant_q": "0.01"})), 2)
        with self.assertRaises(ValueError):
            evaluate.apply_filters(df, {"min_idotp": "0.9"})

    def test_match_runs(self):
        dataset = {"runs": [{"run": "HSA_BS3_1x_DDA_001"}]}
        df = pd.DataFrame({"run": ["/x/hsa-bs3-1x-dda-001.mzML", "other"]})
        kept, dropped = design.match_runs(df, dataset)
        self.assertEqual(list(kept["run"]), ["HSA_BS3_1x_DDA_001"])
        self.assertEqual(dropped, ["other"])


class EndToEnd(unittest.TestCase):
    def test_evaluate(self):
        runs = [f"d{lvl}_{k}" for lvl in ("05", "1") for k in range(3)]
        level = {r: 0.5 if r.startswith("d05") else 1.0 for r in runs}
        rows = []
        for i in range(12):
            for r in runs:
                rows.append({"run": r + ".raw", "precursor": f"p{i}/3", "residue_pair": f"sp|P1|X:{i} -- sp|P1|X:{i + 40}",
                             "quantity": 2 ** (20 + np.log2(level[r]) + 0.01 * (i % 3)), "evidence": "identified"})
        manifest = [{"name": "toy", "design": "dilution", "reference_level": 1.0, "fasta": {"file": "db/x.fasta"},
                     "runs": [{"run": r, "level": level[r]} for r in runs]}]
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "m.json").write_text(json.dumps(manifest))
            pd.DataFrame(rows).to_csv(Path(tmp, "t.tsv"), sep="\t", index=False)
            out = Path(tmp, "r.json")
            evaluate.main(["--dataset", "toy", "--data", tmp, "--manifest", str(Path(tmp, "m.json")),
                           "--arm", "a:table:t.tsv", "--arm", "b:table?evidence=identified:t.tsv", "--out", str(out)])
            r = json.loads(out.read_text())
        own = r["arms"]["a"]["own"]
        self.assertEqual(own["depth"]["residue_pairs"], 12)
        self.assertLess(abs(own["dilution"]["levels"]["0.5"]["median_error"]), 1e-9)
        self.assertEqual(r["common_residue_pairs"], 12)
        self.assertEqual(r["arms"]["a"]["runs_dropped"], [])


if __name__ == "__main__":
    warnings.simplefilter("ignore", RuntimeWarning)
    unittest.main()
