#!/usr/bin/env python3
"""Unit tests for compare_metrics.py.

Run with: python3 -m pytest benchmark/ci/PXD001819/test_compare_metrics.py
or with stdlib: python3 -m unittest benchmark.ci.PXD001819.test_compare_metrics
"""
from __future__ import annotations

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
            "wall_time_sec=120\npsm_1pct_fdr=14000\n",
            "metric\tmin\tmax\toptional\nwall_time_sec\t60\t900\tno\npsm_1pct_fdr\t12000\t17000\tno\n",
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
            "psm_1pct_fdr=14000\n",
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
            "metric\tmin\tmax\toptional\ndistinct_peptides\t\t\tno\nwall_time_sec\t60\t900\tno\n",
        )
        self.assertEqual(r.returncode, 0, r.stderr)


if __name__ == "__main__":
    unittest.main()
