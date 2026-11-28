#!/usr/bin/env python3
"""Unit tests of the cross-link quantification harness on small synthetic inputs (seconds):
  1. residue-pair keys: protein positions, accession reduction, canonical order, and the
     matching of keys between arms for the common view (never by renaming, never by chains);
  2. the long-table contract refuses duplicates and unknown evidence;
  3. the Skyline readers: a pseudo-peptide document (zero areas counted), its site isomers, an
     xiFDR PSM file and both CV report formats;
  4. the metrics: dilution accuracy and slope, the checks of the levels, run shifts, replicate CV
     and completeness, the offsets of single injections, null calls under drift,
     false-transfer rate;
  5. arm filters and run matching, with the cases that must be errors;
  6. the evaluate and published commands on synthetic designs; result files hold no NaN; an
     arm's own scores, its false-transfer rate included, do not depend on the other arms.

    python3 tests/test_xlq.py
"""
import contextlib
import io
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
from xlq import design, evaluate, keys, metrics, published, skyline, table  # noqa: E402

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
     <precursor_peak replicate="r1" area="1000" identified="true" isotope_dotp="0.95" start_time="10.0" end_time="10.5" />
     <precursor_peak replicate="r2" area="500" identified="aligned" isotope_dotp="0.80" start_time="10.1" end_time="10.6" />
     <precursor_peak replicate="r3" identified="false" />
     <precursor_peak replicate="r4" area="0" identified="false" start_time="10.0" end_time="10.5" />
    </precursor_results>
   </precursor>
  </peptide>
  <peptide sequence="HIKLMNPKKAKDEFGK" modified_sequence="HIKLMNPK[+1]K[+26]AK[+1]DEFGK">
   <variable_modifications>
    <variable_modification index_aa="7" modification_name="LinkSite(K,S,T,Y,D)" />
    <variable_modification index_aa="8" modification_name="Kxlink" />
    <variable_modification index_aa="10" modification_name="LinkSite(K,S,T,Y,D)" />
   </variable_modifications>
   <precursor charge="4" modified_sequence="HIKLMNPK[+1.0]K[+26.0]AK[+1.0]DEFGK">
    <precursor_results>
     <precursor_peak replicate="r1" area="800" identified="false" isotope_dotp="0.90" start_time="10.2" end_time="10.7" />
     <precursor_peak replicate="r2" area="400" identified="false" isotope_dotp="0.70" start_time="30.0" end_time="30.5" />
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

    def test_comparison_ids(self):
        short, long, other = "A:10 -- X:5", "A:10;B:57 -- X:5", "B:57 -- X:5"
        ids, report = keys.comparison_ids({"andes": {short, "Q:1 -- Q:9"}, "other": {long, "Q:1 -- Q:9"}})
        self.assertEqual((ids[short], ids[long], ids["Q:1 -- Q:9"]), (short, short, "Q:1 -- Q:9"))
        self.assertEqual((report["matched_sets"], report["ambiguous_sets"]), (1, 0))
        ids, report = keys.comparison_ids({"andes": {short, other}, "other": {long}})
        self.assertEqual((ids[short], ids[long], ids[other]), (short, long, other))
        self.assertEqual((report["matched_sets"], report["ambiguous_sets"]), (0, 1))
        self.assertEqual(keys.comparison_ids({"one": {"A:1 -- A:2"}})[0], {"A:1 -- A:2": "A:1 -- A:2"})

    def test_overlap_that_only_chains_is_no_match(self):
        short, long, other = "A:10 -- X:5", "A:10;B:57 -- X:5", "B:57 -- X:5"
        ids, report = keys.comparison_ids({"one": {short}, "two": {long}, "three": {other}})
        self.assertEqual(len({ids[short], ids[long], ids[other]}), 3)
        self.assertEqual((report["matched_sets"], report["ambiguous_sets"]), (0, 1))


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
            isomers = skyline.site_isomers(path, "Kxlink", "LinkSite")
        self.assertEqual(counts["not_crosslink"], 1)
        self.assertEqual((counts["peaks"], counts["no_area"], counts["zero_area"]), (6, 1, 1))
        self.assertEqual(len(df), counts["peaks"] - counts["no_area"] - counts["zero_area"])
        first = df[df["alpha_peptide"] == "AKDEFGK"]
        self.assertEqual(list(first["run"]), ["r1", "r2"])
        self.assertEqual(set(first["residue_pair"]), {keys.pair_key("P1:3", "P1:11")})
        self.assertEqual(list(first["evidence"]), ["identified", "transferred"])
        row = first.iloc[0]
        self.assertEqual((row["alpha_peptide"], row["alpha_site"], row["beta_peptide"], row["beta_site"]),
                         ("AKDEFGK", 2, "HIKLMNPK", 3))
        self.assertEqual(isomers, {"precursors": 2, "site_isomer_sets": 1, "precursors_in_sets": 2, "overlapping": 2})

    def test_site_isomers_apart(self):
        apart = SKY.replace('start_time="10.2" end_time="10.7"', 'start_time="20.0" end_time="20.5"')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "doc.sky"
            path.write_text(apart, encoding="utf-8")
            isomers = skyline.site_isomers(path, "Kxlink", "LinkSite")
        self.assertEqual((isomers["site_isomer_sets"], isomers["overlapping"]), (1, 0))

    def test_xifdr(self):
        text = ("isTT,Protein1,ProteinLinkPos1,Protein2,ProteinLinkPos2\n"
                "true,P1,11,P1,3\nTRUE,sp|P1|X;sp|Q2|Y,3;40,P1,11\nfalse,P1,5,P1,6\n")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "psm.csv"
            path.write_text(text, encoding="utf-8")
            pairs = skyline.xifdr_residue_pairs(path)
        self.assertEqual(pairs, {"P1:11 -- P1:3", "P1:11 -- P1:3;Q2:40"})

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

    def test_run_offsets(self):
        P = self.dilution_matrix(1.0)
        P["h1b"] -= 2.7
        groups = {"low": ["l1a", "l1b", "l1c"], "high": ["h1a", "h1b", "h1c"]}
        offsets = metrics.run_offsets(P, groups)
        self.assertAlmostEqual(offsets["high"]["offsets"]["h1b"], -1.8, delta=0.1)
        self.assertAlmostEqual(offsets["high"]["spread"], 2.7, delta=0.15)
        self.assertLess(offsets["low"]["spread"], 0.2)
        empty = metrics.run_offsets(P.iloc[0:0], groups)
        self.assertEqual(empty["low"]["spread"], None)

    def test_dilution_level_checks(self):
        M = self.dilution_matrix(1.0)
        with self.assertRaisesRegex(ValueError, "reference level 2"):
            metrics.dilution(M, self.level_of, 2.0)
        with self.assertRaisesRegex(ValueError, "fewer than 2 runs"):
            metrics.dilution(M, dict(self.level_of, l1c=0.5), 1.0)

    def test_run_shifts(self):
        runs = ["a", "b", "c"]
        M = pd.DataFrame({"a": [20.0, 21.0, 22.0], "b": [21.0, 22.0, 23.0], "c": [19.0, 20.0, np.nan]},
                         index=["k1", "k2", "k3"])
        shift = metrics.run_shifts(M, runs)
        self.assertEqual([round(float(shift[r]), 6) for r in runs], [0.0, 1.0, -1.0])
        centred = M - shift
        self.assertAlmostEqual(centred.loc["k1", "b"], centred.loc["k1", "a"])
        self.assertAlmostEqual(centred.loc["k1", "c"], centred.loc["k1", "a"])
        with self.assertRaises(ValueError):
            metrics.run_shifts(pd.DataFrame({"a": [1.0], "b": [np.nan], "c": [np.nan]}), runs)

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
        solo = metrics.precision(M, {"solo": ["a"]})["solo"]
        self.assertEqual((solo["complete"], solo["median_cv_complete"], solo["median_cv_2plus"]), (2, None, None))
        json.dumps(solo, allow_nan=False)

    def test_empty_table(self):
        empty = long_table([], [])
        M = metrics.log_matrix(empty.reindex(columns=["run", "precursor", "residue_pair", "quantity", "evidence"]), ["a", "b"])
        self.assertEqual(M.shape, (0, 2))
        self.assertEqual(metrics.depth(pd.DataFrame(columns=["run", "precursor", "residue_pair", "quantity", "evidence"]), M)["residue_pairs"], 0)

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
        shifted[metrics.split_halves(runs)[0]] += 2.0
        self.assertEqual(metrics.null_calls(shifted, {"rep": runs})["rep"]["false_calls"], 200)

    def test_null_calls_under_drift(self):
        rng = np.random.default_rng(3)
        runs = [f"r{i}" for i in range(8)]
        self.assertEqual(metrics.split_halves(runs), (["r0", "r3", "r4", "r7"], ["r1", "r2", "r5", "r6"]))
        drift = 0.3 * np.arange(8)
        M = pd.DataFrame(20 + drift + rng.normal(0, 0.01, (200, 8)), columns=runs, index=[f"k{i}" for i in range(200)])
        self.assertEqual(metrics.null_calls(M, {"rep": runs})["rep"]["false_calls"], 0)

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
        with self.assertRaisesRegex(ValueError, "identifed"):
            evaluate.apply_filters(df, {"evidence": "identifed"})

    def test_match_runs(self):
        dataset = {"runs": [{"run": "HSA_BS3_1x_DDA_001"}]}
        df = pd.DataFrame({"run": ["/x/hsa-bs3-1x-dda-001.mzML", "other"]})
        kept, dropped = design.match_runs(df, dataset)
        self.assertEqual(list(kept["run"]), ["HSA_BS3_1x_DDA_001"])
        self.assertEqual(dropped, ["other"])
        windows = pd.DataFrame({"run": [r"C:\\data\\HSA_BS3_1x_DDA_001.raw"]})
        self.assertEqual(list(design.match_runs(windows, dataset)[0]["run"]), ["HSA_BS3_1x_DDA_001"])
        twice = pd.DataFrame({"run": ["HSA_BS3_1x_DDA_001.raw", "HSA-BS3-1x-DDA-001.mzML"]})
        with self.assertRaisesRegex(ValueError, "more than once"):
            design.match_runs(twice, dataset)


class ClasenTransfers(unittest.TestCase):
    """The truth of a transfer design reads the proteins an arm reported, whatever other arms are
    scored next to it. The known interaction group is {B, X}; `A:10 -- X:5` is outside it, while
    `A:10;B:57 -- X:5` and `B:57 -- X:5` are inside."""
    SHORT, LONG, OTHER = "A:10 -- X:5", "A:10;B:57 -- X:5", "B:57 -- X:5"
    MANIFEST = [{"name": "toy", "design": "absence", "fasta": {"file": "x.fasta"},
                 "runs": [{"dataset": "cla_b1", "truth": "clasen:batch1"}, {"dataset": "cla_b4", "truth": "clasen:batch4"}]}]

    def setUp(self):
        ev = design._xl_eval(".")
        truth = ev.ClasenTruth.__new__(ev.ClasenTruth)
        truth.groups = {"g1_batch1": {"B", "X"}}
        truth.batch_prots = {"B", "X"}
        ev.truth_for = lambda spec: truth
        self.loader, design._xl_eval = design._xl_eval, lambda data_dir: ev

    def tearDown(self):
        design._xl_eval = self.loader

    def rates(self, arms):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp, "m.json")
            manifest.write_text(json.dumps(self.MANIFEST))
            run = design.load("toy", manifest)["runs"][0]["run"]
            argv = ["--dataset", "toy", "--data", tmp, "--manifest", str(manifest), "--out", str(Path(tmp, "r.json"))]
            for name, key in arms.items():
                pd.DataFrame([{"run": run, "precursor": "p/4", "residue_pair": key, "quantity": 5.0,
                               "evidence": "transferred"}]).to_csv(Path(tmp, f"{name}.tsv"), sep="\t", index=False)
                argv += ["--arm", f"{name}:table:{name}.tsv"]
            with contextlib.redirect_stderr(io.StringIO()):
                evaluate.main(argv)
            r = json.loads(Path(tmp, "r.json").read_text())
        own = {name: a["own"]["transfers"]["false_transfer_rate"] for name, a in r["arms"].items()}
        common = {name: a["common"]["transfers"]["false_transfer_rate"] for name, a in r["arms"].items() if "common" in a}
        return own, common, r

    def test_rate_is_invariant_to_other_arms(self):
        self.assertEqual(self.rates({"two": self.LONG})[0], {"two": 0.0})
        own, common, r = self.rates({"one": self.SHORT, "two": self.LONG})
        self.assertEqual(own, {"one": 1.0, "two": 0.0})
        self.assertEqual(common, {"one": 1.0, "two": 0.0})
        self.assertEqual(r["common_residue_pairs"], 1)
        own, _, _ = self.rates({"one": self.SHORT, "two": self.LONG, "three": self.OTHER})
        self.assertEqual(own, {"one": 1.0, "two": 0.0, "three": 0.0})

    def test_chained_keys_are_not_one_common_residue_pair(self):
        _, _, r = self.rates({"one": self.SHORT, "two": self.LONG, "three": self.OTHER})
        self.assertEqual(r["common_residue_pairs"], 0)
        self.assertEqual((r["key_matching"]["matched_sets"], r["key_matching"]["ambiguous_sets"]), (0, 1))


def strict_json(text):
    """Parse a result file as strict JSON: NaN and Infinity are refused."""
    def refuse(token):
        raise ValueError(f"{token} in a result file")
    return json.loads(text, parse_constant=refuse)


class EndToEnd(unittest.TestCase):
    runs = [f"d{lvl}_{k}" for lvl in ("05", "1") for k in range(3)]
    level = {r: 0.5 if r.startswith("d05") else 1.0 for r in runs}

    def rows(self, offsets=None, key="sp|P1|X:{i} -- sp|P1|X:{j}", suffix=".raw"):
        offsets = offsets or {}
        return [{"run": r + suffix, "precursor": f"p{i}/3", "residue_pair": key.format(i=i, j=i + 40),
                 "quantity": 2 ** (20 + np.log2(self.level[r]) + offsets.get(r, 0.0) + 0.01 * (i % 3)),
                 "evidence": "identified"} for i in range(12) for r in self.runs]

    def evaluate(self, tmp, arms, tables, **dataset):
        manifest = [dict({"name": "toy", "design": "dilution", "reference_level": 1.0, "fasta": {"file": "db/x.fasta"},
                          "runs": [{"run": r, "level": self.level[r]} for r in self.runs]}, **dataset)]
        Path(tmp, "m.json").write_text(json.dumps(manifest))
        for name, rows in tables.items():
            pd.DataFrame(rows).to_csv(Path(tmp, name), sep="\t", index=False)
        out = Path(tmp, "r.json")
        argv = ["--dataset", "toy", "--data", tmp, "--manifest", str(Path(tmp, "m.json")), "--out", str(out)]
        for arm in arms:
            argv += ["--arm", arm]
        with contextlib.redirect_stderr(io.StringIO()) as err:
            evaluate.main(argv)
        return strict_json(out.read_text()), err.getvalue()

    def test_evaluate(self):
        with tempfile.TemporaryDirectory() as tmp:
            r, _ = self.evaluate(tmp, ["a:table:t.tsv", "b:table?evidence=identified:t.tsv"], {"t.tsv": self.rows()})
        own = r["arms"]["a"]["own"]
        self.assertEqual(own["depth"]["residue_pairs"], 12)
        self.assertLess(abs(own["dilution"]["levels"]["0.5"]["median_error"]), 1e-9)
        self.assertEqual(r["common_residue_pairs"], 12)
        self.assertEqual((r["arms"]["a"]["runs_dropped"], r["arms"]["a"]["runs_without_rows"]), ([], []))
        self.assertEqual(r["key_matching"]["matched_sets"], 0)
        self.assertLess(own["run_offsets"]["1"]["spread"], 1e-9)

    def test_normalise(self):
        offsets = {"d1_1": 1.5}
        with tempfile.TemporaryDirectory() as tmp:
            raw, _ = self.evaluate(tmp, ["a:table:t.tsv"], {"t.tsv": self.rows(offsets)})
        with tempfile.TemporaryDirectory() as tmp:
            norm, _ = self.evaluate(tmp, ["a:table:t.tsv"], {"t.tsv": self.rows(offsets)}, normalise=True)
        self.assertAlmostEqual(raw["arms"]["a"]["own"]["run_offsets"]["1"]["spread"], 1.5, places=3)
        self.assertLess(norm["arms"]["a"]["own"]["run_offsets"]["1"]["spread"], 1e-6)
        self.assertGreater(raw["arms"]["a"]["own"]["precision"]["1"]["median_cv_complete"], 0.5)
        self.assertLess(norm["arms"]["a"]["own"]["precision"]["1"]["median_cv_complete"], 1e-6)

    def test_keys_join_across_arms(self):
        other = self.rows(key="P1:{i};Q9:{i} -- P1:{j}")
        with tempfile.TemporaryDirectory() as tmp:
            r, _ = self.evaluate(tmp, ["andes:table:a.tsv", "other:table:b.tsv"], {"a.tsv": self.rows(), "b.tsv": other})
        self.assertEqual(r["key_matching"]["matched_sets"], 12)
        self.assertEqual(r["common_residue_pairs"], 12)
        self.assertEqual(r["arms"]["other"]["common"]["depth"]["residue_pairs"], 12)

    def test_own_scores_do_not_depend_on_other_arms(self):
        other = self.rows(key="P1:{i};Q9:{i} -- P1:{j}")
        with tempfile.TemporaryDirectory() as tmp:
            alone, _ = self.evaluate(tmp, ["other:table:b.tsv"], {"b.tsv": other})
        with tempfile.TemporaryDirectory() as tmp:
            both, _ = self.evaluate(tmp, ["andes:table:a.tsv", "other:table:b.tsv"], {"a.tsv": self.rows(), "b.tsv": other})
        self.assertEqual(both["arms"]["other"]["own"], alone["arms"]["other"]["own"])

    def test_arm_without_rows_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit) as stop:
                self.evaluate(tmp, ["a:table:t.tsv"], {"t.tsv": self.rows(suffix="_other.raw")})
            self.assertIn("none of its runs", str(stop.exception))
            with self.assertRaises(SystemExit) as stop:
                self.evaluate(tmp, ["a:table?evidence=transferred:t.tsv"], {"t.tsv": self.rows()})
            self.assertIn("filters keep no row", str(stop.exception))

    def test_missing_runs_are_named(self):
        rows = [r for r in self.rows() if not r["run"].startswith("d1_2")]
        with tempfile.TemporaryDirectory() as tmp:
            r, err = self.evaluate(tmp, ["a:table:t.tsv"], {"t.tsv": rows})
        self.assertEqual(r["arms"]["a"]["runs_without_rows"], ["d1_2"])
        self.assertIn("no row for 1 of the 6 runs", err)

    def test_published(self):
        cv = ("File Name,Peptide Modified Sequence,Product Charge,Fragment Ion,Cv Area\n"
              "a.raw,PEPK,3,precursor,10%\na.raw,PEPR,3,precursor,30%\n")
        psms = "isTT,Protein1,ProteinLinkPos1,Protein2,ProteinLinkPos2\ntrue,P1,3,P1,11\ntrue,P1,5,P1,9\n"
        manifest = [{"name": "toy", "design": "replicates", "fasta": {"file": "x.fasta"}, "runs": [],
                     "published": {"skyline_cv_injection": {"file": "cv.csv"},
                                   "skyline_document": {"file": "doc.sky", "spacer_mod": "Kxlink", "site_mod": "LinkSite"},
                                   "xifdr_psms": {"file": "psm.csv"}}}]
        with tempfile.TemporaryDirectory() as tmp:
            for name, text in (("m.json", json.dumps(manifest)), ("cv.csv", cv), ("psm.csv", psms),
                               ("doc.sky", SKY), ("x.fasta", f">sp|P1|X\n{PROTEIN}\n")):
                Path(tmp, name).write_text(text, encoding="utf-8")
            out = Path(tmp, "p.json")
            with contextlib.redirect_stderr(io.StringIO()):
                published.main(["--dataset", "toy", "--data", tmp, "--manifest", str(Path(tmp, "m.json")), "--out", str(out)])
            r = strict_json(out.read_text())
        self.assertEqual(r["cv_reports"]["injection"]["median_cv"], {"precursor": 0.2})
        check = r["skyline_document"]["xifdr_psms"]
        self.assertEqual((check["residue_pairs_document"], check["residue_pairs_xifdr"], check["shared"]), (2, 2, 1))
        self.assertEqual(r["skyline_document"]["site_isomers"]["site_isomer_sets"], 1)


if __name__ == "__main__":
    warnings.simplefilter("ignore", RuntimeWarning)
    unittest.main()
