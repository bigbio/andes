#!/usr/bin/env python3
"""Tests for ../proteobench/andes_to_proteobench.py (andes MBR output -> ProteoBench Custom table).

The peptidoform tests need only the standard library; the end-to-end test writes a small
mbr.parquet and quant directory and needs pyarrow (skipped without it).

    python3 -I docs/benchmarks/tests/test_proteobench_convert.py
"""
import csv
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "proteobench"))
import andes_to_proteobench as conv  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover
    pa = None


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Peptidoform(unittest.TestCase):
    def test_unimod_accessions(self):
        self.assertEqual(conv.to_proforma("PEPC[UNIMOD:4]M[UNIMOD:35]K"),
                         ("PEPCMK", "PEPC[Carbamidomethyl]M[Oxidation]K"))

    def test_mass_deltas(self):
        self.assertEqual(conv.to_proforma("AC[+57.02146]M[+15.99491]R")[1],
                         "AC[Carbamidomethyl]M[Oxidation]R")

    def test_nterm_acetyl_prefix(self):
        self.assertEqual(conv.to_proforma("A[UNIMOD:1]PEPK")[1], "[Acetyl]-APEPK")
        self.assertEqual(conv.to_proforma("S[+42.01057]PEPK")[1], "[Acetyl]-SPEPK")

    def test_summed_first_residue_delta_is_split(self):
        self.assertEqual(conv.to_proforma("C[+99.03203]PEPK")[1], "[Acetyl]-C[Carbamidomethyl]PEPK")
        self.assertEqual(conv.to_proforma("M[+58.00548]PEPK")[1], "[Acetyl]-M[Oxidation]PEPK")

    def test_unknown_modification_is_an_error(self):
        with self.assertRaises(ValueError):
            conv.to_proforma("PEPS[+79.96633]K")
        with self.assertRaises(ValueError):
            conv.to_proforma("PEP-K")

    def test_run_name(self):
        self.assertEqual(conv.run_name("/data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML"),
                         "LFQ_Astral_DDA_15min_50ng_Condition_A_REP1")


@unittest.skipIf(pa is None, "pyarrow not installed")
class EndToEnd(unittest.TestCase):
    def test_table(self):
        runs = conv.MODULE_RUNS["qe"]
        with tempfile.TemporaryDirectory() as d:
            qd = os.path.join(d, "q")
            os.mkdir(qd)
            pq.write_table(pa.table({
                "precursor_key": ["PEPC[UNIMOD:4]K/2", "PEPC[UNIMOD:4]K/2", "XXXK/2", "AAAK/3", "DDDR/2"],
                "peptidoform": ["PEPC[UNIMOD:4]K", "PEPC[UNIMOD:4]K", "XXXK", "AAAK", "DDDR"],
                "charge": pa.array([2, 2, 2, 3, 2], pa.uint8()),
                "proteins": ["sp|P1|A_HUMAN", "sp|P1|A_HUMAN", "XXX_sp|P9|Z_ECOLI",
                             "sp|P2|B_YEAST;sp|P3|C_YEAST", "sp|P4|D_ECOLI"],
                "is_decoy": [False, False, True, False, False],
            }), os.path.join(qd, f"{runs[0]}.quant.parquet"))
            mbr = os.path.join(d, "pb.mbr.parquet")
            rows = [  # key, run, evidence, quantity, quantity_normalized
                ("PEPC[UNIMOD:4]K/2", runs[0], "identified", 100.0, 90.0),
                ("PEPC[UNIMOD:4]K/2", runs[3], "transfer_ms1", 50.0, 55.0),
                ("PEPC[UNIMOD:4]K/2", runs[4], "below_threshold", 7.0, 7.0),      # not reported
                ("AAAK/3", runs[1], "identified_below_threshold", 10.0, 10.0),    # not reported
                ("AAAK/3", runs[2], "transfer_ms2", 20.0, None),                  # normalised missing
                ("DDDR/2", runs[5], "identified", 0.0, 0.0),                      # no signal
            ]
            pq.write_table(pa.table({
                "precursor_key": [r[0] for r in rows],
                "run": [f"/x/{r[1]}.mzML" for r in rows],
                "evidence_type": [r[2] for r in rows],
                "quantity": [r[3] for r in rows],
                "quantity_normalized": [r[4] for r in rows],
            }), mbr)
            out = os.path.join(d, "raw.tsv")
            self.assertEqual(conv.main(["--mbr", mbr, "--quant-dir", qd, "--module", "qe",
                                        "--column", "quantity", "-o", out]), 0)
            got = read_tsv(out)
            self.assertEqual([r["Modified sequence"] for r in got], ["AAAK", "PEPC[Carbamidomethyl]K"])
            pep = got[1]
            self.assertEqual((pep["Sequence"], pep["Charge"], pep["Proteins"]), ("PEPCK", "2", "sp|P1|A_HUMAN"))
            self.assertEqual([pep[r] for r in runs], ["100.0", "", "", "50.0", "", ""])
            self.assertEqual(got[0]["Proteins"], "sp|P2|B_YEAST;sp|P3|C_YEAST")
            self.assertEqual([got[0][r] for r in runs], ["", "", "20.0", "", "", ""])

            out = os.path.join(d, "norm.tsv")
            conv.main(["--mbr", mbr, "--quant-dir", qd, "--module", "qe", "-o", out])
            got = read_tsv(out)
            self.assertEqual(len(got), 1)  # AAAK has no normalised value
            self.assertEqual([got[0][r] for r in runs], ["90.0", "", "", "55.0", "", ""])

            out = os.path.join(d, "id.tsv")
            conv.main(["--mbr", mbr, "--quant-dir", qd, "--module", "qe", "--column", "quantity",
                       "--evidence", "identified", "-o", out])
            got = read_tsv(out)
            self.assertEqual([got[0][r] for r in runs], ["100.0", "", "", "", "", ""])

    def test_wrong_module_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            mbr = os.path.join(d, "pb.mbr.parquet")
            pq.write_table(pa.table({"precursor_key": ["AK/2"], "run": ["other_run"],
                                     "evidence_type": ["identified"], "quantity": [1.0],
                                     "quantity_normalized": [1.0]}), mbr)
            with self.assertRaises(SystemExit):
                conv.main(["--mbr", mbr, "--quant-dir", d, "--module", "astral",
                           "--column", "quantity", "-o", os.path.join(d, "o.tsv")])


if __name__ == "__main__":
    unittest.main()
