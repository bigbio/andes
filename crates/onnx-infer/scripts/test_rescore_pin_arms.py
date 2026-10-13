#!/usr/bin/env python3
"""rescore_pin_arms.py: an arm built only from MS2 blocks never needs RT calibration.

Requires: numpy, pandas. Run with pytest or as a script.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parent / "rescore_pin_arms.py"
MS2 = [
    "spec_pearson_norm",
    "ionb_pearson_norm",
    "iony_pearson_norm",
    "spec_mse_norm",
    "dotprod_norm",
    "cos_norm",
    "spectral_angle",
    "pred_explained",
]


def write_inputs(d: Path, n_rows: int, n_confident: int):
    rng = np.random.default_rng(0)
    ids = [f"s{i}" for i in range(n_rows)]
    labels = [1 if i % 2 == 0 else -1 for i in range(n_rows)]
    with open(d / "search.pin", "w") as f:
        f.write("SpecId\tLabel\tScanNr\tScore\tPeptide\tProteins\n")
        for i, (s, lab) in enumerate(zip(ids, labels)):
            f.write(f"{s}\t{lab}\t{i}\t{rng.random():.4f}\tK.PEPTIDEK.R\tP1\n")
    with open(d / "pdeep.tsv", "w") as f:
        f.write("\t".join(["SpecId", "observed_retention_time", "pdeep_irt", "pdeep_missing"] + MS2) + "\n")
        for i, s in enumerate(ids):
            ms2 = [f"{v:.4f}" for v in rng.random(len(MS2))]
            f.write("\t".join([s, f"{10 + i * 0.1:.4f}", f"{0.1 * i:.6f}", "0"] + ms2) + "\n")
    targets = [s for s, lab in zip(ids, labels) if lab == 1]
    with open(d / "base.psms", "w") as f:
        f.write("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n")
        for k, s in enumerate(targets):
            q = 0.001 if k < n_confident else 0.5
            f.write(f"{s}\t1.0\t{q}\t0.01\tK.PEPTIDEK.R\tP1\n")


def run_arms(d: Path, arm: str):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--pin", str(d / "search.pin"),
         "--features", f"pdeep={d / 'pdeep.tsv'}", "--baseline-psms", str(d / "base.psms"),
         "--out-prefix", str(d / "run"), "--arm", arm],
        capture_output=True, text=True,
    )


def check_ms2_only(n_confident: int):
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        write_inputs(d, 40, n_confident)
        r = run_arms(d, "ms2only=pdeep:ms2")
        assert r.returncode == 0, f"{n_confident} confident targets:\n{r.stderr}"
        header = (d / "run_ms2only.pin").read_text().splitlines()[0].split("\t")
        assert [c for c in header if c.startswith("pdeep_")] == [f"pdeep_{c}" for c in MS2] + ["pdeep_missing"]


def test_ms2_only_arm_without_confident_targets():
    check_ms2_only(0)


def test_ms2_only_arm_with_too_few_targets_for_the_rt_bins():
    check_ms2_only(10)


def test_rt_arm_still_calibrates():
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        write_inputs(d, 600, 300)
        r = run_arms(d, "rt=pdeep:rt")
        assert r.returncode == 0, r.stderr
        header = (d / "run_rt.pin").read_text().splitlines()[0].split("\t")
        assert {"pdeep_rt_diff", "pdeep_predicted_rt", "pdeep_missing"} <= set(header)
        assert "RT calibration on 300 confident targets" in r.stdout


if __name__ == "__main__":
    test_ms2_only_arm_without_confident_targets()
    test_ms2_only_arm_with_too_few_targets_for_the_rt_bins()
    test_rt_arm_still_calibrates()
    print("rescore_pin_arms tests passed")
