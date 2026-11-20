#!/usr/bin/env python3
"""Tests of andes_to_proteobench.py for the DIA Astral module.

Run: python3 -m unittest docs/benchmarks/proteobench/test_andes_to_proteobench.py
(needs pyarrow).
"""
import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import andes_to_proteobench as conv  # noqa: E402

RUNS = [f"LFQ_Astral_DIA_15min_50ng_Condition_{c}_REP{i}" for c in "AB" for i in (1, 2, 3)]


def write_dia_quant(path, rows):
    import pyarrow as pa
    import pyarrow.parquet as pq

    cols = ["filename", "peptide", "charge", "proteins", "evidence_type", "q_value",
            "quantity", "quantity_frag", "quantity_ms1", "quantity_normalized"]
    data = {c: [r.get(c) for r in rows] for c in cols}
    schema = pa.schema([
        ("filename", pa.string()), ("peptide", pa.string()), ("charge", pa.uint8()),
        ("proteins", pa.string()), ("evidence_type", pa.string()), ("q_value", pa.float64()),
        ("quantity", pa.float64()), ("quantity_frag", pa.float64()), ("quantity_ms1", pa.float64()),
        ("quantity_normalized", pa.float64()),
    ])
    t = pa.table(data, schema=schema).replace_schema_metadata({"andes.dia_quant.schema": "1.0"})
    pq.write_table(t, path)


def row(run, pep, z, v, ev="identified", prot="sp|P1|A_HUMAN"):
    return {"filename": run, "peptide": pep, "charge": z, "proteins": prot, "evidence_type": ev,
            "q_value": 0.001, "quantity": v, "quantity_frag": v,
            "quantity_ms1": None if v is None else 10 * v,
            "quantity_normalized": None if v is None else v / 2}


def read_tsv(path):
    with open(path) as fh:
        return list(csv.reader(fh, delimiter="\t"))


class DiaAstral(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.pq = os.path.join(self.dir.name, "study.dia_quant.parquet")
        self.out = os.path.join(self.dir.name, "custom_input.tsv")

    def tearDown(self):
        self.dir.cleanup()

    def run_conv(self, extra=()):
        return conv.main(["--dia-quant", self.pq, "--module", "astral_dia", "-o", self.out, *extra])

    def test_dia_rows_map_to_module_columns(self):
        rows = [row(r, "_.PEPTIDEK.-", 2, 100.0 * (i + 1)) for i, r in enumerate(RUNS)]
        rows.append(row(RUNS[0], "K.LESLIEK.-", 3, 7.0, ev="transfer",
                        prot="sp|P2|B_YEAST;sp|P3|C_YEAST"))
        rows.append(row(RUNS[1], "K.LESLIEK.-", 3, None))
        rows.append(row(RUNS[2], "K.NOQUANTK.-", 2, None))
        write_dia_quant(self.pq, rows)
        self.assertEqual(self.run_conv(), 0)
        t = read_tsv(self.out)
        self.assertEqual(t[0], ["Sequence", "Modified sequence", "Proteins", "Charge"] + RUNS)
        body = {(r[0], r[3]): r for r in t[1:]}
        self.assertEqual(len(body), 2, "a precursor without any quantity is not written")
        pep = body[("PEPTIDEK", "2")]
        self.assertEqual(pep[1], "PEPTIDEK")
        self.assertEqual(pep[4:], ["50.0", "100.0", "150.0", "200.0", "250.0", "300.0"])
        les = body[("LESLIEK", "3")]
        self.assertEqual(les[2], "sp|P2|B_YEAST;sp|P3|C_YEAST")
        self.assertEqual(les[4:], ["3.5", "", "", "", "", ""], "transfer kept, missing empty")
        self.assertEqual(self.run_conv(["--column", "quantity_ms1"]), 0)
        les = {(r[0], r[3]): r for r in read_tsv(self.out)[1:]}[("LESLIEK", "3")]
        self.assertEqual(les[4], "70.0")
        self.assertEqual(self.run_conv(["--evidence", "identified"]), 0)
        keys = {(r[0], r[3]) for r in read_tsv(self.out)[1:]}
        self.assertEqual(keys, {("PEPTIDEK", "2")})

    def test_modifications_spelled_as_module_expects(self):
        cases = {
            "_.C+57.02146PEPTIDEK.-": ("CPEPTIDEK", "C[Carbamidomethyl]PEPTIDEK"),
            "K.PEPM+15.99492TIDEK.-": ("PEPMTIDEK", "PEPM[Oxidation]TIDEK"),
            "-.M+42.01057PEPTIDEK.-": ("MPEPTIDEK", "[Acetyl]-MPEPTIDEK"),
            "-.C+99.03203PEPTM+15.99492IDEK.-": ("CPEPTMIDEK", "[Acetyl]-C[Carbamidomethyl]PEPTM[Oxidation]IDEK"),
            "-.M+58.00548PEPTIDEK.-": ("MPEPTIDEK", "[Acetyl]-M[Oxidation]PEPTIDEK"),
        }
        for pin, want in cases.items():
            self.assertEqual(conv.to_proforma(conv.pin_peptidoform(pin)), want, pin)
        with self.assertRaises(ValueError):
            conv.to_proforma(conv.pin_peptidoform("-.S+79.96633PEPTIDEK.-"))
        with self.assertRaises(ValueError):
            conv.pin_peptidoform("PEPTIDEK")

    def test_unknown_run_names_fail_with_expected_list(self):
        rows = [row(r, "_.PEPTIDEK.-", 2, 1.0) for r in RUNS[:5]]
        rows.append(row("LFQ_Astral_DIA_15min_50ng_Human_01", "_.PEPTIDEK.-", 2, 1.0))
        write_dia_quant(self.pq, rows)
        with self.assertRaises(SystemExit) as e:
            self.run_conv()
        msg = str(e.exception)
        self.assertIn("LFQ_Astral_DIA_15min_50ng_Human_01", msg)
        for r in RUNS:
            self.assertIn(r, msg)
        self.assertFalse(os.path.exists(self.out))


if __name__ == "__main__":
    unittest.main()
