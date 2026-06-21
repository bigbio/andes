#!/usr/bin/env python3
"""Synthetic tests for the model-independence gate.

Uses tiny in-memory parquet stores (no MS-GF+ data) so CI can exercise the gate
logic without the seed ever being present in this repo. Run with `pytest` or
directly (`python test_check_models_independence.py`).
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_models_independence as gate  # noqa: E402

_KEYS = gate._KEY_COLS


def _write_store(path, models):
    """models: {model_id: [(values_list), ...]} -> one `table` record per list."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    cols = {c: [] for c in ["record_kind", "model_id", "values", *_KEYS]}
    for mid, tables in models.items():
        for ti, vals in enumerate(tables):
            cols["record_kind"].append("table")
            cols["model_id"].append(mid)
            cols["values"].append([float(x) for x in vals])
            # distinct key per table index; constant across stores so identical
            # `values` => identical fingerprint.
            cols["table_kind"].append(f"k{ti}")
            cols["ion_kind"].append("b")
            cols["ion_charge"].append(1)
            cols["ion_offset_bits"].append(0)
            cols["part_charge"].append(2)
            cols["part_mass_bits"].append(0)
            cols["part_seg"].append(0)
    pq.write_table(pa.table(cols), path)


def test_detects_seed_identical_and_passes_distinct():
    with tempfile.TemporaryDirectory() as d:
        seed = os.path.join(d, "seed.parquet")
        bundle = os.path.join(d, "bundle.parquet")
        # seed model A has tables [[1,2,3]]
        _write_store(seed, {"A": [[1, 2, 3]]})
        # bundle: X copies A's tables (leak), Y is distinct (clean),
        # Z is a relabel of A's tables under a different id (still a leak).
        _write_store(bundle, {"X": [[1, 2, 3]], "Y": [[4, 5, 6]], "Z": [[1, 2, 3]]})

        leaks = gate.find_leaks(bundle, seed)
        leaked_ids = {mid for mid, _sid, _h in leaks}
        assert leaked_ids == {"X", "Z"}, leaked_ids
        # every leak points back at seed model A
        assert all(sid == "A" for _mid, sid, _h in leaks)


def test_main_exit_codes_and_waiver():
    with tempfile.TemporaryDirectory() as d:
        seed = os.path.join(d, "seed.parquet")
        bundle = os.path.join(d, "bundle.parquet")
        _write_store(seed, {"A": [[1, 2, 3]]})
        _write_store(bundle, {"X": [[1, 2, 3]], "Y": [[4, 5, 6]]})

        # X leaks -> exit 1
        assert gate.main(["--models", bundle, "--seed", seed]) == 1
        # waiving X -> exit 0
        assert gate.main(["--models", bundle, "--seed", seed, "--allow", "X"]) == 0

        # a fully clean bundle -> exit 0
        clean = os.path.join(d, "clean.parquet")
        _write_store(clean, {"Y": [[4, 5, 6]]})
        assert gate.main(["--models", clean, "--seed", seed]) == 0


if __name__ == "__main__":
    test_detects_seed_identical_and_passes_distinct()
    test_main_exit_codes_and_waiver()
    print("ok: all gate tests passed")
