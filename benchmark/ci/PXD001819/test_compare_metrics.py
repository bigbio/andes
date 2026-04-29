#!/usr/bin/env python3
"""Unit tests for compare_metrics.py + extract_metrics.parse_pin.

Run with: python3 -m unittest benchmark.ci.PXD001819.test_compare_metrics
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
import textwrap
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("compare_metrics.py")


class CompareMetricsTest(unittest.TestCase):
    def _run(self, metrics_text: str, baseline_text: str) -> subprocess.CompletedProcess[str]:
        tmp = Path(self.id().replace(".", "_"))
        tmp.mkdir(exist_ok=True)
        metrics = tmp / "metrics.txt"
        baseline = tmp / "baseline.tsv"
        metrics.write_text(textwrap.dedent(metrics_text))
        baseline.write_text(textwrap.dedent(baseline_text))
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(metrics), str(baseline)],
            capture_output=True,
            text=True,
        )

    def tearDown(self) -> None:
        tmp = Path(self.id().replace(".", "_"))
        if tmp.exists():
            for p in tmp.iterdir():
                p.unlink()
            tmp.rmdir()

    def test_all_in_range_passes(self) -> None:
        r = self._run(
            "wall_time_sec=120\nnative_target_count=28000\n",
            "metric\tmin\tmax\toptional\nwall_time_sec\t60\t900\tno\nnative_target_count\t14000\t35000\tno\n",
        )
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("within baseline ranges", r.stdout)

    def test_out_of_range_fails(self) -> None:
        r = self._run(
            "wall_time_sec=2000\n",
            "metric\tmin\tmax\toptional\nwall_time_sec\t60\t900\tno\n",
        )
        self.assertEqual(r.returncode, 1)
        self.assertIn("outside", r.stderr)

    def test_missing_required_fails(self) -> None:
        r = self._run(
            "native_target_count=28000\n",
            "metric\tmin\tmax\toptional\nwall_time_sec\t60\t900\tno\n",
        )
        self.assertEqual(r.returncode, 1)
        self.assertIn("missing", r.stderr)

    def test_missing_optional_warns(self) -> None:
        r = self._run(
            "wall_time_sec=120\n",
            "metric\tmin\tmax\toptional\npeak_rss_kb\t0\t999999\tyes\n",
        )
        self.assertEqual(r.returncode, 0)
        self.assertIn("warning", r.stderr)

    def test_na_value_treated_as_missing(self) -> None:
        r = self._run(
            "peak_rss_kb=NA\n",
            "metric\tmin\tmax\toptional\npeak_rss_kb\t0\t999999\tyes\n",
        )
        self.assertEqual(r.returncode, 0)
        self.assertIn("warning", r.stderr)

    def test_non_numeric_fails(self) -> None:
        r = self._run(
            "wall_time_sec=abc\n",
            "metric\tmin\tmax\toptional\nwall_time_sec\t60\t900\tno\n",
        )
        self.assertEqual(r.returncode, 1)
        self.assertIn("not numeric", r.stderr)

    def test_empty_range_row_is_skipped(self) -> None:
        r = self._run(
            "wall_time_sec=120\n",
            "metric\tmin\tmax\toptional\ncpu_percent\t\t\tno\nwall_time_sec\t60\t900\tno\n",
        )
        self.assertEqual(r.returncode, 0, r.stderr)


def _load_extract_metrics():
    spec = importlib.util.spec_from_file_location(
        "extract_metrics", Path(__file__).with_name("extract_metrics.py")
    )
    em = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(em)
    return em


class ParsePinTest(unittest.TestCase):
    """Verify extract_metrics.parse_pin counts target / decoy rows correctly."""

    def setUp(self) -> None:
        self.em = _load_extract_metrics()
        self.tmp = Path(self.id().replace(".", "_"))
        self.tmp.mkdir(exist_ok=True)

    def tearDown(self) -> None:
        for p in self.tmp.iterdir():
            p.unlink()
        self.tmp.rmdir()

    def test_parse_pin_counts_labels(self) -> None:
        pin = self.tmp / "tiny.pin"
        pin.write_text(
            "SpecId\tLabel\tScanNr\tFeatures\n"
            "spec1\t1\t100\tx\n"
            "spec2\t-1\t101\tx\n"
            "spec3\t1\t102\tx\n"
            "spec4\t1\t103\tx\n"
            "spec5\t-1\t104\tx\n"
        )
        targets, decoys = self.em.parse_pin(pin)
        self.assertEqual(targets, 3)
        self.assertEqual(decoys, 2)

    def test_parse_pin_empty_returns_zeros(self) -> None:
        pin = self.tmp / "empty.pin"
        pin.write_text("SpecId\tLabel\tScanNr\tFeatures\n")
        targets, decoys = self.em.parse_pin(pin)
        self.assertEqual(targets, 0)
        self.assertEqual(decoys, 0)

    def test_parse_pin_skips_malformed_rows(self) -> None:
        pin = self.tmp / "malformed.pin"
        pin.write_text(
            "SpecId\tLabel\tScanNr\tFeatures\n"
            "spec1\t1\t100\tx\n"
            "incomplete\n"
            "spec2\t0\t102\tx\n"
            "spec3\t-1\t103\tx\n"
        )
        targets, decoys = self.em.parse_pin(pin)
        self.assertEqual(targets, 1)
        self.assertEqual(decoys, 1)


if __name__ == "__main__":
    unittest.main()
