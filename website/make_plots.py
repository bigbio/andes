#!/usr/bin/env python3
"""Render the benchmark figures from benchmarks.json as SVG.

usage: make_plots.py <benchmarks.json> <out_dir>
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Colour-blind-safe (Okabe-Ito); one colour per engine on every figure.
COLOURS = {
    "andes": "#0072B2",
    "Comet 2025.01": "#E69F00",
    "Java MS-GF+ v20240326": "#999999",
}
INK = "#333333"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.edgecolor": "#888888",
    "axes.labelcolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
})


def grouped(ax, data, datasets, ylabel, fmt, log=False):
    engines = list(data)
    width = 0.8 / len(engines)
    for i, eng in enumerate(engines):
        xs = [j + (i - (len(engines) - 1) / 2) * width for j in range(len(datasets))]
        vals = [v if v is not None else 0 for v in data[eng]]
        bars = ax.bar(xs, vals, width * 0.92, label=eng, color=COLOURS[eng])
        for bar, v in zip(bars, data[eng]):
            if v is None:
                continue
            ax.annotate(fmt(v), (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                        ha="center", va="bottom", fontsize=7.5, color=INK,
                        xytext=(0, 2), textcoords="offset points")
    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels(datasets)
    ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale("log")
    ax.margins(y=0.12)
    ax.legend(frameon=False, fontsize=8)


def save(fig, out, name):
    fig.tight_layout()
    fig.savefig(out / name, format="svg", transparent=True)
    plt.close(fig)


def main(src, out_dir):
    d = json.loads(Path(src).read_text())
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ds = d["datasets"]

    fig, ax = plt.subplots(figsize=(8, 4))
    grouped(ax, d["psms"], ds, "PSMs at Percolator q ≤ 0.01", lambda v: f"{v:,.0f}")
    ax.set_title("Identifications", loc="left", fontsize=11, color=INK)
    save(fig, out, "psms.svg")

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    for ax, key, label, fmt, title in (
        (axes[0], "wall_s", "seconds (log scale)", lambda v: f"{v:,.0f}", "Wall time"),
        (axes[1], "cpu_s", "CPU seconds (log scale)", lambda v: f"{v:,.0f}", "CPU time"),
        (axes[2], "peak_gb", "GB", lambda v: f"{v:.1f}", "Peak memory"),
    ):
        grouped(ax, d[key], [x.split("\n")[0] for x in ds], label, fmt, log=key != "peak_gb")
        ax.set_title(title, loc="left", fontsize=11, color=INK)
        if ax is not axes[0]:
            ax.get_legend().remove()
        ax.tick_params(axis="x", labelrotation=30)
    save(fig, out, "cost.svg")

    rows = d["true_fdp"]["rows"]
    fig, ax = plt.subplots(figsize=(8, 4.4))
    for i, r in enumerate(rows):
        ax.plot([r["lo"], r["hi"]], [i, i], color=COLOURS["andes"], lw=5, solid_capstyle="round")
        ax.annotate(f'{r["lo"]:.2f}–{r["hi"]:.2f}%', (2.45, i), va="center", fontsize=8, color=INK)
    ax.axvline(1.0, color=INK, lw=1, ls="--")
    ax.annotate("nominal 1%", (1.0, -0.9), ha="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0.5, 2.95)
    ax.set_xticks([0.5, 1.0, 1.5, 2.0])
    ax.set_xlabel("entrapment-estimated true FDP at Percolator q ≤ 0.01 (%)")
    ax.set_title("Does the nominal 1% hold?", loc="left", fontsize=11, color=INK)
    save(fig, out, "true_fdp.svg")

    o = d["opt_in"]
    fig, ax = plt.subplots(figsize=(5, 3.6))
    mids = [(a + b) / 2 for a, b in zip(o["psms_lo"], o["psms_hi"])]
    errs = [[m - a for m, a in zip(mids, o["psms_lo"])], [b - m for m, b in zip(mids, o["psms_hi"])]]
    bars = ax.bar(o["labels"], mids, yerr=errs, color=[COLOURS["Java MS-GF+ v20240326"], COLOURS["andes"]],
                  width=0.55, capsize=4)
    for bar, a, b in zip(bars, o["psms_lo"], o["psms_hi"]):
        ax.annotate(f"{a:,}–{b:,}", (bar.get_x() + bar.get_width() / 2, b), ha="center",
                    xytext=(0, 4), textcoords="offset points", fontsize=8, color=INK)
    ax.set_ylabel("PSMs at a true FDP of 1%")
    ax.set_title("Chimeric mode (Astral)", loc="left", fontsize=11, color=INK)
    ax.margins(y=0.15)
    save(fig, out, "chimeric.svg")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
