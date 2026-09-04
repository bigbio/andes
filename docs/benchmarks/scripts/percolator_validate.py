#!/usr/bin/env python3
"""Percolator-output sanity validator for MS-GF+ benchmark arms.

Purpose: after every Percolator run, flag pathological target/decoy distributions
or feature-weight artifacts that a raw ID count would miss. Runs on the standard
Percolator output files:

    <prefix>.percolator.target.psms.txt
    <prefix>.percolator.decoy.psms.txt
    <prefix>.percolator.weights.txt

Emits <label>_report.md + PNG plots. Exits non-zero on HARD check failure.

Usage:
    python percolator_validate.py <label> <percolator_output_dir> [<prefix>]

    <label>    Free-form tag (e.g. "baseline_pxd001819", "branch_pxd001819").
               Used as filename stem for the report and plots.
    <dir>      Directory containing Percolator .psms.txt / .weights.txt files.
    <prefix>   (Optional) Filename prefix before ".percolator.*.psms.txt".
               If omitted, auto-detected from files present.

HARD checks (exit 1 if any fails):
    - Mean target score > mean decoy score
    - ≥ 100 targets accepted at q <= 0.01
    - PEP monotonically non-increasing in score (up to float noise)
    - Decoy distribution is monomodal (no secondary peak > 50% of primary peak height)

WARN checks (reported but do not fail):
    - No feature weight with |w| > 10× median |w|
    - Target/decoy overlap at 1%-FDR score threshold < 30%
    - No q-value plateau > 5% of ranks (flat regions = degenerate rescoring)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


HARD_FAIL = "HARD FAIL"
WARN = "WARN"
OK = "OK"


def load_psms(target_path: Path, decoy_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    t = pd.read_csv(target_path, sep="\t")
    d = pd.read_csv(decoy_path, sep="\t")
    return t, d


def load_weights(weights_path: Path) -> pd.DataFrame:
    # Percolator weights file: 3 lines per iteration, columns = feature names,
    # rows = iteration weights. Last iteration's weights are the learned model.
    with weights_path.open() as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if len(lines) < 2:
        return pd.DataFrame()
    header = lines[0].split("\t")
    # Last non-header line is the final-iteration weights row.
    weights = [float(x) for x in lines[-1].split("\t")]
    return pd.DataFrame({"feature": header, "weight": weights})


def check_score_ordering(t: pd.DataFrame, d: pd.DataFrame) -> tuple[str, str]:
    mt = t["score"].mean()
    md = d["score"].mean()
    if mt <= md:
        return HARD_FAIL, f"mean(target score)={mt:.3f} <= mean(decoy score)={md:.3f} — catastrophic T/D inversion"
    return OK, f"mean(target)={mt:.3f} > mean(decoy)={md:.3f} (gap={mt - md:.3f})"


def check_target_count_1pct(t: pd.DataFrame) -> tuple[str, str]:
    n = int((t["q-value"] <= 0.01).sum())
    if n < 100:
        return HARD_FAIL, f"only {n} targets at 1% q-value (threshold 100)"
    return OK, f"{n} targets at 1% q-value"


def check_pep_monotone(t: pd.DataFrame) -> tuple[str, str]:
    if "posterior_error_prob" not in t.columns:
        return WARN, "no posterior_error_prob column — skipped"
    s = t.sort_values("score", ascending=False)["posterior_error_prob"].to_numpy()
    # PEP should be non-decreasing as we walk from best→worst score (i.e. non-increasing in score).
    diffs = np.diff(s)
    # Tolerance: small floating noise (1e-12). Anything bigger is a real monotonicity break.
    violations = int((diffs < -1e-12).sum())
    if violations > len(s) * 0.001:  # allow 0.1% numerical noise
        return HARD_FAIL, f"PEP monotonicity violated at {violations}/{len(s)} positions ({100*violations/len(s):.2f}%)"
    if violations > 0:
        return WARN, f"PEP monotonicity: {violations} minor violations (below threshold)"
    return OK, "PEP monotonically non-increasing in score"


def check_decoy_monomodal(d: pd.DataFrame, label: str, outdir: Path) -> tuple[str, str]:
    scores = d["score"].to_numpy()
    if len(scores) < 50:
        return WARN, f"only {len(scores)} decoys — histogram shape check skipped"
    hist, edges = np.histogram(scores, bins=50)
    # Smooth slightly so a jagged tail doesn't look like a peak.
    kernel = np.ones(3) / 3.0
    smooth = np.convolve(hist, kernel, mode="same")
    primary_idx = int(smooth.argmax())
    primary_h = smooth[primary_idx]
    # find_peaks returns indices of local maxima with a minimum prominence.
    peaks, _ = find_peaks(smooth, height=primary_h * 0.5)
    # Exclude the primary.
    secondaries = [p for p in peaks if p != primary_idx]
    if secondaries:
        return HARD_FAIL, (
            f"decoy distribution is multimodal — secondary peak at score "
            f"{edges[secondaries[0]]:.3f} is ≥ 50% of primary at {edges[primary_idx]:.3f}. "
            f"Likely causes: feature leakage, bad decoy generation, or score miscalibration."
        )
    return OK, "decoy distribution is monomodal"


def check_feature_weights(weights: pd.DataFrame) -> tuple[str, str]:
    if weights.empty:
        return WARN, "no weights file — skipped"
    w = weights["weight"].abs().to_numpy()
    if len(w) < 3:
        return OK, "fewer than 3 features — runaway check skipped"
    # Compute median excluding zero-weight features (they just mean "ignored by Percolator").
    nonzero = w[w > 1e-12]
    if len(nonzero) < 2:
        return WARN, "nearly all feature weights are zero — Percolator may not have converged"
    median = np.median(nonzero)
    runaway = weights[weights["weight"].abs() > 10 * median]
    if not runaway.empty:
        top = runaway.sort_values("weight", key=lambda s: s.abs(), ascending=False).iloc[0]
        return WARN, (
            f"feature '{top['feature']}' weight {top['weight']:.2f} is "
            f"{abs(top['weight']) / median:.1f}× median ({median:.3f}) — possible feature dominance"
        )
    return OK, f"no runaway weights; median |w|={median:.3f}, max |w|={w.max():.3f}"


def check_overlap_1pct(t: pd.DataFrame, d: pd.DataFrame) -> tuple[str, str]:
    if (t["q-value"] <= 0.01).sum() == 0:
        return WARN, "no targets at 1% q — overlap check skipped"
    threshold = t.loc[t["q-value"] <= 0.01, "score"].min()
    decoys_above = int((d["score"] >= threshold).sum())
    targets_above = int((t["score"] >= threshold).sum())
    overlap = decoys_above / max(targets_above, 1)
    if overlap > 0.30:
        return WARN, (
            f"T/D overlap at 1%-FDR score threshold ({threshold:.3f}) is {100*overlap:.1f}% — "
            f"higher than expected (decoys={decoys_above}, targets={targets_above})"
        )
    return OK, f"T/D overlap at 1%-FDR threshold = {100*overlap:.1f}% (decoys/targets at score ≥ {threshold:.3f})"


def check_qvalue_plateau(t: pd.DataFrame) -> tuple[str, str]:
    if len(t) < 100:
        return WARN, "too few targets for plateau check"
    s = t.sort_values("score", ascending=False)["q-value"].to_numpy()
    diffs = np.diff(s)
    flat_run = 0
    max_flat = 0
    for diff in diffs:
        if diff == 0:
            flat_run += 1
            max_flat = max(max_flat, flat_run)
        else:
            flat_run = 0
    fraction = max_flat / len(s)
    if fraction > 0.05:
        return WARN, (
            f"largest q-value plateau spans {max_flat} ranks ({100*fraction:.1f}% of targets) — "
            f"may indicate degenerate/tied features"
        )
    return OK, f"largest q-value plateau = {max_flat} ranks ({100*fraction:.2f}% of targets)"


def plot_td_histogram(t: pd.DataFrame, d: pd.DataFrame, label: str, outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(min(t["score"].min(), d["score"].min()),
                       max(t["score"].max(), d["score"].max()), 60)
    ax.hist(d["score"], bins=bins, alpha=0.6, label=f"decoy (n={len(d)})", color="#d62728")
    ax.hist(t["score"], bins=bins, alpha=0.6, label=f"target (n={len(t)})", color="#2ca02c")
    ax.set_xlabel("Percolator score")
    ax.set_ylabel("count")
    ax.set_title(f"T/D score distribution — {label}")
    ax.legend()
    ax.grid(alpha=0.3)
    path = outdir / f"{label}_td_hist.png"
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def plot_pep_curve(t: pd.DataFrame, label: str, outdir: Path) -> Path | None:
    if "posterior_error_prob" not in t.columns:
        return None
    s = t.sort_values("score", ascending=False)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(s["score"], s["posterior_error_prob"], linewidth=0.8)
    ax.set_xlabel("Percolator score")
    ax.set_ylabel("posterior error probability")
    ax.set_title(f"PEP vs score — {label}")
    ax.grid(alpha=0.3)
    ax.invert_xaxis()
    path = outdir / f"{label}_pep_curve.png"
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def plot_qvalue_curve(t: pd.DataFrame, label: str, outdir: Path) -> Path:
    s = t.sort_values("q-value").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(s.index + 1, s["q-value"], linewidth=1.0)
    ax.axhline(0.01, color="red", linestyle="--", linewidth=0.6, label="1% FDR")
    ax.set_xlabel("target rank (sorted by q-value)")
    ax.set_ylabel("q-value")
    ax.set_title(f"q-value vs rank — {label}")
    ax.set_yscale("log")
    ax.legend()
    ax.grid(alpha=0.3)
    path = outdir / f"{label}_qvalue_curve.png"
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def write_report(label: str, outdir: Path, results: list[tuple[str, str, str]],
                 plots: list[Path], t: pd.DataFrame, d: pd.DataFrame) -> Path:
    path = outdir / f"{label}_report.md"
    n_hard = sum(1 for _, s, _ in results if s == HARD_FAIL)
    n_warn = sum(1 for _, s, _ in results if s == WARN)
    verdict = "FAIL" if n_hard else ("PASS (with warnings)" if n_warn else "PASS")

    with path.open("w") as f:
        f.write(f"# Percolator validation — {label}\n\n")
        f.write(f"**Verdict:** {verdict} — {n_hard} hard failures, {n_warn} warnings\n\n")
        f.write(f"**Targets:** {len(t)}  |  **Decoys:** {len(d)}\n\n")
        f.write(f"**Targets at 1% q-value:** {int((t['q-value'] <= 0.01).sum())}\n")
        f.write(f"**Targets at 5% q-value:** {int((t['q-value'] <= 0.05).sum())}\n\n")
        f.write("## Checks\n\n")
        f.write("| Check | Status | Detail |\n|---|---|---|\n")
        for name, status, detail in results:
            f.write(f"| {name} | **{status}** | {detail} |\n")
        f.write("\n## Plots\n\n")
        for p in plots:
            f.write(f"![{p.stem}]({p.name})\n\n")
    return path


def autodetect_prefix(outdir: Path) -> str:
    """Derive the percolator file prefix from whatever files the user's pipeline produced."""
    target_files = list(outdir.glob("*.percolator.target.psms.txt"))
    if not target_files:
        sys.exit(f"no *.percolator.target.psms.txt found in {outdir}")
    name = target_files[0].name
    return name.split(".percolator.target.psms.txt")[0]


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    label = sys.argv[1]
    outdir = Path(sys.argv[2]).resolve()
    prefix = sys.argv[3] if len(sys.argv) > 3 else autodetect_prefix(outdir)

    target_path = outdir / f"{prefix}.percolator.target.psms.txt"
    decoy_path = outdir / f"{prefix}.percolator.decoy.psms.txt"
    weights_path = outdir / f"{prefix}.percolator.weights.txt"

    for p in (target_path, decoy_path):
        if not p.exists():
            sys.exit(f"missing required file: {p}")

    t, d = load_psms(target_path, decoy_path)
    weights = load_weights(weights_path) if weights_path.exists() else pd.DataFrame()

    checks = [
        ("T/D score ordering", *check_score_ordering(t, d)),
        ("Target count at 1% q", *check_target_count_1pct(t)),
        ("PEP monotonicity", *check_pep_monotone(t)),
        ("Decoy distribution monomodal", *check_decoy_monomodal(d, label, outdir)),
        ("Feature-weight sanity", *check_feature_weights(weights)),
        ("T/D overlap at 1% FDR", *check_overlap_1pct(t, d)),
        ("q-value plateau", *check_qvalue_plateau(t)),
    ]

    plots = [plot_td_histogram(t, d, label, outdir), plot_qvalue_curve(t, label, outdir)]
    pep_plot = plot_pep_curve(t, label, outdir)
    if pep_plot:
        plots.append(pep_plot)

    report_path = write_report(label, outdir, checks, plots, t, d)

    hard_fails = [(n, s, detail) for n, s, detail in checks if s == HARD_FAIL]
    warns = [(n, s, detail) for n, s, detail in checks if s == WARN]

    print(f"report: {report_path}")
    for name, _, detail in hard_fails:
        print(f"  FAIL  {name}: {detail}", file=sys.stderr)
    for name, _, detail in warns:
        print(f"  WARN  {name}: {detail}")

    return 1 if hard_fails else 0


if __name__ == "__main__":
    sys.exit(main())
