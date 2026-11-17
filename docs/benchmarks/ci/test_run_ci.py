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
    keys = {"PEPH/2": "sp|P1|X_HUMAN", "PEPY/2": "sp|P2|X_YEAST", "PEPM/2": "sp|P3|X_HUMAN;sp|P4|X_YEAST",
            "PEPZ/2": "sp|P5|Z_YEAST"}
    pq.write_table(pa.table({"precursor_key": list(keys), "proteins": list(keys.values()),
                             "is_decoy": [False] * 4}), d / "q" / "r.quant.parquet")
    rows = []
    for k, a, lq in (("PEPH/2", 1, 0.001), ("PEPY/2", 2, 0.001), ("PEPM/2", 1, 0.001), ("PEPZ/2", 2, 0.0)):
        for run, v in (("A1", a), ("A2", a), ("B1", 1), ("B2", 1)):
            rows.append((k, run + ".raw", "identified", 1e6 * v, None, lq))
    rows.append(("PEPH/2", "A1.raw", "transfer_ms1", 9e9, 0.5, 0.001))   # transfer_q too high: dropped
    rows.append(("PEPY/2", "A1.raw", "identified", 9e9, None, None))    # no library_q: dropped
    pq.write_table(pa.table(dict(zip(["precursor_key", "run", "evidence_type", "quantity_normalized",
                                      "transfer_q", "library_q"], zip(*rows)))), d / "m.parquet")
    q = R.quant(d / "m.parquet", d / "q", d / "design.tsv", d / "amounts.tsv")
    by = {g["group"]: g for g in q["ratios"]}
    assert set(by) == {"HUMAN", "YEAST"}, by           # the two-species precursor is not a group
    assert by["HUMAN"]["comparison"] == "A/B" and abs(by["HUMAN"]["median"]) < 1e-9
    assert abs(by["YEAST"]["median"] - 1) < 1e-9 and abs(by["YEAST"]["expected"] - math.log2(2)) < 1e-9
    assert by["YEAST"]["mae"] < 1e-9 and by["YEAST"]["n"] == 2   # library_q 0.0 is kept
    assert q["quantified"] == 4
    assert all(c["median_cv"] < 1e-9 for c in q["cv"])


def test_bench_records():
    """A partial quant download is `missing`; a failing evaluation keeps IDs and cost as `failed`."""
    d = Path(tempfile.mkdtemp())
    (d / "databases").mkdir()
    (d / "databases" / "db.fasta").write_text(">sp|P1|X_HUMAN\nPEPTIDEK\n")
    (d / "s").mkdir()
    for run in ("A1", "A2", "B1"):                     # B2 of the design is absent
        (d / "s" / f"{run}.mzML").write_text("")
    (d / "design.tsv").write_text("run\tcondition\nA1\tA\nA2\tA\nB1\tB\nB2\tB\n")
    ds = {"fetch": "s", "database": "db.fasta", "quant": {"design": str(d / "design.tsv"), "amounts": "none"}}
    r = R.bench("x", ds, d, "false", 1, d / "w")
    assert r["status"] == "missing" and "B2" in r["error"], r

    (d / "s" / "B2.mzML").write_text("")
    fake = d / "fake-andes"                             # writes a statistics log, no mbr.parquet
    fake.write_text("#!/bin/sh\nwhile [ $# -gt 0 ]; do [ \"$1\" = --output-pin ] && p=$2; shift; done\n"
                    "printf '  Spectra searched : 10\\n  Identified PSMs : 4\\n' > \"${p%.pin}.statistics.log\"\n")
    fake.chmod(0o755)
    r = R.bench("x", ds, d, str(fake), 1, d / "w")
    assert r["status"] == "failed" and r["error"].startswith("evaluation:") and r["psms"] == 4, r


if __name__ == "__main__":
    test_quant()
    test_bench_records()
    print("ok")
