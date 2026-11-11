#!/usr/bin/env python3
"""Figures for the glyco LFQ + MBR benchmark.

Reads `glyco_quant_numbers.json` + the lfq table and writes:

  fig1_identification_gate.png   accepted-feature counts and entrapment by gate
  fig2_matrix_completeness.png   glycopeptides reaching >=k runs, identified vs +MBR
  fig3_cross_run.png             5x5 Pearson log2-area heatmaps, identified vs +MBR
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GREY = "#b9b7ae"
CMAP = LinearSegmentedColormap.from_list(
    "blue_seq", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
                 "#184f95", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "axes.spines.top": False,
    "axes.spines.right": False,
})


def fig1(d, out):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.4), dpi=200,
                                   sharey=True,
                                   gridspec_kw={"width_ratios": [1, 1]})
    fig.subplots_adjust(left=0.20, right=0.985, top=0.77, bottom=0.14,
                        wspace=0.10)
    n_id = d["gated_counts"]["identified"]
    n_tr = d["gated_counts"]["transferred"]
    f = d["features"]
    id_all = sum(x["acc_identified"] for x in f)
    tr_all = sum(x["acc_transferred"] for x in f)
    rows = [("identified · no ID gate", id_all, GREY, 0.55),
            ("transferred · no ID gate", tr_all, GREY, 0.55),
            ("identified · ID-gated", n_id, BLUE, 1.0),
            ("transferred · ID-gated", n_tr, ORANGE, 1.0)]
    y = np.arange(len(rows))[::-1]

    for yi, (_, v, col, al) in zip(y, rows):
        ax1.plot([v], [yi], "o", ms=9, color=col, alpha=al,
                 markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=4)
        ax1.text(v * 1.18, yi, f"{v:,}", va="center", ha="left",
                 fontsize=8.4, color=INK2)
    ax1.set_xscale("log")
    ax1.set_xlim(4e3, 4.5e5)
    ax1.set_yticks(y)
    ax1.set_yticklabels([r[0] for r in rows], fontsize=8.4)
    ax1.set_xlabel("accepted quant features (log scale)")
    ax1.set_title("What the tables report", fontsize=10.5, pad=8, loc="left")
    ax1.grid(axis="y", visible=False)

    # entrapment share, same row order
    ebe = d["entrap_by_evidence"]
    def share(ev, gated):
        for gl, info in ebe[ev].items():
            is_gated = "gate" not in gl
            if is_gated == gated:
                return 100 * info["raw"]
        return float("nan")
    vals = [share("Identified", False), share("Transferred", False),
            share("Identified", True), share("Transferred", True)]
    for yi, (row, v) in zip(y, zip(rows, vals)):
        ax2.plot([v], [yi], "o", ms=9, color=row[2], alpha=row[3],
                 markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=4)
        right = v >= 1.0
        ax2.text(v * (1.20 if right else 0.80), yi,
                 f"{v:.2f}%" if v < 10 else f"{v:.1f}%",
                 va="center", ha="left" if right else "right",
                 fontsize=8.4, color=INK2)
    factor = d["sequon_space"]["factor"]
    nominal = 100 * d["q_id"] / factor
    ax2.axvline(nominal, color=MUTED, lw=1.2, ls="--", zorder=3)
    ax2.text(nominal * 1.15, 3.35,
             f"1% FDP ⇒ {nominal:.2f}% entrapment",
             fontsize=7.4, color=INK2, va="top", ha="left")
    ax2.set_xscale("log")
    ax2.set_xlim(0.15, 100)
    ax2.set_xlabel("ENTRAP_ share of accepted features (%, log scale)")
    ax2.set_title("Entrapment: the gate recovers the nominal FDR",
                  fontsize=10.5, pad=8, loc="left")
    ax2.set_ylim(-0.6, 3.6)
    ax2.grid(axis="y", visible=False)

    fig.suptitle("Glycopeptide quantification: identification gate", y=0.975,
                 fontsize=12, x=0.20, ha="left", color=INK)
    fig.text(0.20, 0.885,
             "Mouse liver, 5 fractions pooled.  `Accepted` is a feature/XIC "
             "q-value; identification FDR comes from Percolator on the .glyco.pin",
             fontsize=8.2, color=INK2, ha="left")
    p = os.path.join(out, "fig1_identification_gate.png")
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)


def fig2(d, out):
    runs = d["runs"]
    bi = d["blocks"]["gated_identified"]
    bm = d["blocks"]["gated_with_mbr"]
    ks = [f">={k}" for k in range(1, len(runs) + 1)]
    x = np.arange(len(ks))
    w = 0.38
    fig, ax = plt.subplots(figsize=(9.6, 5.0), dpi=200)
    fig.subplots_adjust(left=0.10, right=0.985, top=0.80, bottom=0.16)
    vi = [bi["presence"][k] for k in ks]
    vm = [bm["presence"][k] for k in ks]
    ax.bar(x - w / 2 - 0.01, vi, w, color=BLUE, edgecolor=SURFACE, lw=1.5,
           zorder=4, label="identified")
    ax.bar(x + w / 2 + 0.01, vm, w, color=AQUA, edgecolor=SURFACE, lw=1.5,
           zorder=4, label="identified + MBR")
    for xi, (a, b) in enumerate(zip(vi, vm)):
        ax.text(xi - w / 2 - 0.01, a + 20, f"{a:,}", ha="center", va="bottom",
                fontsize=7.6, color=INK2, rotation=90)
        ax.text(xi + w / 2 + 0.01, b + 20, f"{b:,}", ha="center", va="bottom",
                fontsize=7.6, color=INK2, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([f"≥{k} of {len(runs)} runs" for k in range(1, len(runs) + 1)])
    ax.set_ylabel("glycopeptides (peptide, site, composition)")
    ax.set_ylim(0, max(vm) * 1.22)
    ax.grid(axis="x", visible=False)
    ax.legend(frameon=False, loc="upper right", fontsize=8.4)
    fig.suptitle("Quantification matrix completeness — matched between runs",
                 y=0.975, fontsize=11.5, x=0.10, ha="left", color=INK)
    fig.text(0.10, 0.905,
             f"ID-gated (Percolator q ≤ {d['q_id']}); matrix fill "
             f"{bi['grid_pct']:.1f}% → {bm['grid_pct']:.1f}% of the "
             f"{bm['n_glycopeptides']:,} × {len(runs)} grid",
             fontsize=8.4, color=INK2, ha="left")
    p = os.path.join(out, "fig2_matrix_completeness.png")
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)


def heat(r, runs, title, ax):
    n = len(runs)
    mat = np.full((n, n), np.nan)
    for x in r:
        i, j = runs.index(x["a"]), runs.index(x["b"])
        mat[i, j] = mat[j, i] = x["pearson"]
    np.fill_diagonal(mat, 1.0)
    im = ax.imshow(mat, cmap=CMAP, vmin=0.5, vmax=1.0)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([s.replace("MouseLiver-Z-T-", "T-") for s in runs],
                       fontsize=8)
    ax.set_yticklabels([s.replace("MouseLiver-Z-T-", "T-") for s in runs],
                       fontsize=8)
    ax.set_title(title, fontsize=10.5, pad=8)
    for i in range(n):
        for j in range(n):
            v = mat[i, j]
            ax.text(j, i, f"{v:.2f}" if np.isfinite(v) else "—",
                    ha="center", va="center", fontsize=7.6,
                    color=SURFACE if v > 0.82 else INK)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def fig3(d, out):
    runs = d["runs"]
    bi = d["blocks"]["gated_identified"]
    bm = d["blocks"]["gated_with_mbr"]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.9), dpi=200)
    fig.subplots_adjust(left=0.10, right=0.88, top=0.77, bottom=0.15,
                        wspace=0.30)
    im = heat(bi["pairs"], runs, "identified", axes[0])
    heat(bm["pairs"], runs, "identified + MBR", axes[1])
    cax = fig.add_axes([0.90, 0.15, 0.016, 0.62])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Pearson r, log2 area", fontsize=8, color=INK2)
    cb.ax.tick_params(labelsize=7, colors=MUTED)
    cb.outline.set_edgecolor(AXIS)
    fig.suptitle("Cross-run consistency (technical replicates, 5 fractions)",
                 y=0.965, fontsize=12, x=0.10, ha="left", color=INK)
    fig.text(0.10, 0.875,
             f"median pairwise r: {bi['median_pair_r']:.3f} identified → "
             f"{bm['median_pair_r']:.3f} with MBR.  Glycoform-share mean |Δ| "
             f"{bi['share_mae']:.3f} → {bm['share_mae']:.3f}",
             fontsize=8.4, color=INK2, ha="left")
    p = os.path.join(out, "fig3_cross_run.png")
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    with open(args.json) as fh:
        d = json.load(fh)
    fig1(d, args.out)
    fig2(d, args.out)
    fig3(d, args.out)


if __name__ == "__main__":
    main()
