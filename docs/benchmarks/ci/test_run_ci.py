"""Check of run_ci.quant on a synthetic design with known ratios: python3 test_run_ci.py"""
import math
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

import run_ci as R


def test_quant():
    d = Path(tempfile.mkdtemp())
    (d / "design.tsv").write_text("run\tcondition\nA1\tA\nA2\tA\nB1\tB\nB2\tB\n")
    (d / "amounts.tsv").write_text("species\tcondition\tamount\nHUMAN\t*\t1\nYEAST\tA\t2\nYEAST\tB\t1\n")
    (d / "q").mkdir()
    keys = {"PEPH/2": "sp|P1|X_HUMAN", "PEPY/2": "sp|P2|X_YEAST", "PEPM/2": "sp|P3|X_HUMAN;sp|P4|X_YEAST"}
    pq.write_table(pa.table({"precursor_key": list(keys), "proteins": list(keys.values()),
                             "is_decoy": [False] * 3}), d / "q" / "r.quant.parquet")
    rows = []
    for k, a in (("PEPH/2", 1), ("PEPY/2", 2), ("PEPM/2", 1)):
        for run, v in (("A1", a), ("A2", a), ("B1", 1), ("B2", 1)):
            rows.append((k, run + ".raw", "identified", 1e6 * v, None, 0.001))
    rows.append(("PEPH/2", "A1.raw", "transfer_ms1", 9e9, 0.5, 0.001))   # transfer_q too high: dropped
    pq.write_table(pa.table(dict(zip(["precursor_key", "run", "evidence_type", "quantity_normalized",
                                      "transfer_q", "library_q"], zip(*rows)))), d / "m.parquet")
    q = R.quant(d / "m.parquet", d / "q", d / "design.tsv", d / "amounts.tsv")
    by = {g["group"]: g for g in q["ratios"]}
    assert set(by) == {"HUMAN", "YEAST"}, by           # the two-species precursor is not a group
    assert by["HUMAN"]["comparison"] == "A/B" and abs(by["HUMAN"]["median"]) < 1e-9
    assert abs(by["YEAST"]["median"] - 1) < 1e-9 and abs(by["YEAST"]["expected"] - math.log2(2)) < 1e-9
    assert by["YEAST"]["mae"] < 1e-9 and q["quantified"] == 3
    assert all(c["median_cv"] < 1e-9 for c in q["cv"])


if __name__ == "__main__":
    test_quant()
    print("ok")
