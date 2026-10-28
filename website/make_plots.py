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


def range_plot(rows, xlabel, title, xmax, ref=None, ticks=None, figsize=(8, 3), label_x=None):
    fig, ax = plt.subplots(figsize=figsize)
    label_x = xmax * 0.84 if label_x is None else label_x
    for i, r in enumerate(rows):
        ax.plot([r["lo"], r["hi"]], [i, i], color=COLOURS["andes"], lw=5, solid_capstyle="round")
        ax.annotate(f'{r["lo"]:.2f}–{r["hi"]:.2f}%', (label_x, i), va="center", fontsize=8, color=INK)
    if ref is not None:
        ax.axvline(ref, color=INK, lw=1, ls="--")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, xmax)
    if ticks:
        ax.set_xticks(ticks)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    return fig


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

    k = d["knockoff"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    cols = ["#999999", "#999999", COLOURS["andes"], COLOURS["andes"]]
    for ax, key, title, fmt in (
        (axes[0], "true_fdp_pct", "True FDP of reported variant PSMs (%)", lambda v: f"{v:g}%"),
        (axes[1], "variant_psms", "Variant PSMs reported", lambda v: f"{v:,}"),
    ):
        bars = ax.bar(range(4), k[key], color=cols, width=0.65)
        for bar, v in zip(bars, k[key]):
            ax.annotate(fmt(v), (bar.get_x() + bar.get_width() / 2, v), ha="center",
                        xytext=(0, 2), textcoords="offset points", fontsize=8, color=INK)
        ax.set_xticks(range(4))
        ax.set_xticklabels(k["labels"], fontsize=7)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.margins(y=0.12)
    save(fig, out, "knockoff.svg")

    t = d["tiers"]
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2))
    for ax, key, title, fmt in (
        (axes[0], "carried_variants", "Carried variants found", lambda v: f"{v}"),
        (axes[1], "true_fdp_variants_pct", "True FDP, variants (%)", lambda v: f"{v:g}%"),
    ):
        bars = ax.bar(t["labels"], t[key], color=[COLOURS["Java MS-GF+ v20240326"], COLOURS["andes"]], width=0.55)
        for bar, v in zip(bars, t[key]):
            ax.annotate(fmt(v), (bar.get_x() + bar.get_width() / 2, v), ha="center",
                        xytext=(0, 2), textcoords="offset points", fontsize=8, color=INK)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.margins(y=0.15)
        ax.tick_params(axis="x", labelsize=8)
    save(fig, out, "tiers.svg")

    ti = d["tmt_interference"]
    fig, ax = plt.subplots(figsize=(8, 3.8))
    series = {"andes (quant gate)": COLOURS["andes"], "without the gate": COLOURS["Java MS-GF+ v20240326"],
              "stricter purity filter (≥ 0.9)": COLOURS["Comet 2025.01"]}
    w = 0.26
    for i, (name, col) in enumerate(series.items()):
        xs = [j + (i - 1) * w for j in range(3)]
        bars = ax.bar(xs, ti[name], w * 0.92, label=name, color=col)
        for bar, v in zip(bars, ti[name]):
            ax.annotate(f"{v:.3f}", (bar.get_x() + bar.get_width() / 2, v), ha="center",
                        xytext=(0, 2), textcoords="offset points", fontsize=7.5, color=INK)
    ax.set_xticks(range(3))
    ax.set_xticklabels(ti["acquisitions"])
    ax.set_ylabel("knockout-channel interference (lower is better)")
    ax.set_title("TMT reporter interference, yeast knockout (PXD017803)", loc="left", fontsize=11, color=INK)
    ax.legend(frameon=False, fontsize=8)
    ax.margins(y=0.15)
    save(fig, out, "tmt_interference.svg")

    r = d["retrieval"]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
    for ax, lab, b, a in zip(axes, r["labels"], r["before"], r["after"]):
        bars = ax.bar([r["before_label"].split(" / ")[axes.tolist().index(ax)], r["after_label"].split(" / ")[axes.tolist().index(ax)]],
                      [b, a], color=[COLOURS["Java MS-GF+ v20240326"], COLOURS["andes"]], width=0.55)
        for bar, v in zip(bars, (b, a)):
            ax.annotate(f"{v:g}", (bar.get_x() + bar.get_width() / 2, v), ha="center",
                        xytext=(0, 2), textcoords="offset points", fontsize=8, color=INK)
        ax.set_title(lab, loc="left", fontsize=10, color=INK)
        ax.tick_params(axis="x", labelsize=7.5)
        ax.margins(y=0.15)
    save(fig, out, "retrieval.svg")

    pf = d["protein_fdp"]["rows"]
    fig = range_plot(pf, "entrapment-estimated true FDP of protein groups (%)", "Protein groups, Astral",
                     8, ref=None, ticks=[0, 1, 2, 3, 4, 5, 6, 7], figsize=(7, 2.4), label_x=7.2)
    ax = fig.axes[0]
    for x in (1, 5):
        ax.axvline(x, color=INK, lw=1, ls="--")
    save(fig, out, "protein_fdp.svg")

    g = d["glyco_pxd065680"]
    fig, ax = plt.subplots(figsize=(8, 4))
    rows = g["ranked"]
    for i, r in enumerate(rows):
        ax.barh(i, r["v"], color=COLOURS["andes"] if r["andes"] else COLOURS["Java MS-GF+ v20240326"], height=0.7)
        ax.annotate(f'{r["v"]:,}', (r["v"], i), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["n"] for r in rows], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 1450)
    ax.set_xlabel("unique sequence + glycan composition identifications")
    ax.set_title("Glyco benchmark PXD065680: counts (not error rates)", loc="left", fontsize=11, color=INK)
    save(fig, out, "glyco_counts.svg")

    ov = g["overlap"]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for i, (t, tot, sh) in enumerate(zip(ov["tools"], ov["total"], ov["shared"])):
        ax.barh(i, sh, color=COLOURS["andes"], height=0.65)
        ax.barh(i, tot - sh, left=sh, color=COLOURS["Java MS-GF+ v20240326"], height=0.65)
        ax.annotate(f"{sh:,} of {tot:,} ({100 * sh / tot:.0f}%)", (tot, i), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(ov["tools"])))
    ax.set_yticklabels(ov["tools"], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 1550)
    ax.set_xlabel("glycopeptides in each tool's list (blue: also reported by andes arm A2)")
    ax.set_title("Overlap with each published tool", loc="left", fontsize=11, color=INK)
    save(fig, out, "glyco_overlap.svg")

    sup = g["support"]
    fig, ax = plt.subplots(figsize=(7, 3.4))
    labels = ["none\n(andes only)", "1 tool", "2 tools", "3 tools", "4 tools", "all 5"]
    bars = ax.bar(labels, sup, color=[COLOURS["Comet 2025.01"]] + [COLOURS["andes"]] * 5, width=0.65)
    for bar, v in zip(bars, sup):
        ax.annotate(f"{v}", (bar.get_x() + bar.get_width() / 2, v), ha="center", xytext=(0, 2),
                    textcoords="offset points", fontsize=8, color=INK)
    ax.set_ylabel("andes A2 identifications")
    ax.set_title("How many of the five tools also report each andes identification", loc="left", fontsize=10.5, color=INK)
    ax.margins(y=0.15)
    ax.tick_params(axis="x", labelsize=8)
    save(fig, out, "glyco_support.svg")

    tr = g["tiers"]
    fig, ax = plt.subplots(figsize=(8, 3.4))
    tier_cols = [COLOURS["andes"], "#56B4E9", COLOURS["Comet 2025.01"], "#D55E00"]
    left = [0.0] * len(tr["sets"])
    for key, name, col in zip(("t1", "t2", "t3", "t4"), ("Tier 1", "Tier 2", "Tier 3", "Tier 4"), tier_cols):
        ax.barh(range(len(tr["sets"])), tr[key], left=left, color=col, height=0.65, label=name)
        left = [a + b for a, b in zip(left, tr[key])]
    ax.set_yticks(range(len(tr["sets"])))
    ax.set_yticklabels(tr["sets"], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of PSMs")
    ax.set_title("Evidence tiers (Tier 1: oxonium, b/y and Y0/Y1 ions all present)", loc="left", fontsize=10.5, color=INK)
    ax.legend(frameon=False, fontsize=8, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.45))
    save(fig, out, "glyco_tiers.svg")

    q = d["quant_accuracy"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for ax, key, title in ((axes[0], "QExactive", "ProteoBench Q Exactive module"), (axes[1], "Astral", "ProteoBench Astral module")):
        rows = q[key]
        bars = ax.bar(range(len(rows)), [r["v"] for r in rows],
                      color=[COLOURS["andes"] if r["andes"] else COLOURS["Java MS-GF+ v20240326"] for r in rows], width=0.6)
        for bar, r in zip(bars, rows):
            extra = f'\n{r["k"]:,} precursors' if r["k"] else ""
            ax.annotate(f'{r["v"]:.3f}{extra}', (bar.get_x() + bar.get_width() / 2, r["v"]), ha="center",
                        xytext=(0, 2), textcoords="offset points", fontsize=7.5, color=INK)
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels([r["n"] for r in rows], fontsize=8)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.set_ylabel("median |ε| (lower is better)")
        ax.margins(y=0.3)
    save(fig, out, "quant_accuracy.svg")

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
