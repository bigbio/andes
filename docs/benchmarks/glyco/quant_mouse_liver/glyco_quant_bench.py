#!/usr/bin/env python3
"""Glycopeptide LFQ + MBR benchmark — mouse liver, 5 fractions (deep tier).

Reads the tables written by `andes --glyco --lfq --mbr` and, when a
Percolator `out.psms` is supplied, gates every quantified feature on the
IDENTIFICATION q-value of the PSM it came from.  The `Accepted` flag in the
`.glyco.lfq.tsv` table is a feature/measurement q-value (XIC quality), NOT an
identification FDR; without the join the tables are dominated by random
rank-1 glyco matches (see the entrapment section).

Usage:
  python glyco_quant_bench.py --dir DIR --stem STEM --out OUTDIR \
      [--psms out.psms] [--fasta mouse_entrap.fasta] [--q-id 0.01]

Environment: python with pandas/numpy; pure stdlib otherwise.
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

Q_ACCEPT = 0.05          # andes' feature FDR default
Q_ID = 0.01              # Percolator identification q-value gate
ENTRAP_COMMON_DIR = os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))


def log2(x):
    return np.log2(np.asarray(x, float))


def stats(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return None
    return {"n": int(len(x)),
            "pearson": float(np.corrcoef(x, y)[0, 1]),
            "spearman": float(pd.Series(x).corr(pd.Series(y), method="spearman")),
            "median_abs": float(np.median(np.abs(x - y))),
            "p95_abs": float(np.percentile(np.abs(x - y), 95))}


def qstats(s):
    s = np.sort(np.asarray(pd.to_numeric(s, errors="coerce"), float))
    s = s[np.isfinite(s)]
    if not len(s):
        return {"n": 0}
    return {"n": int(len(s)), "median": float(np.median(s)),
            "p90": float(np.percentile(s, 90)),
            "p99": float(np.percentile(s, 99)),
            "max": float(s.max())}


def read_lfq(path):
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in ("Area", "MonoArea", "ApexRT", "FeatureQValue", "Cosine"):
        if c in d:
            d[c] = pd.to_numeric(d[c].replace("", np.nan), errors="coerce")
    d["Accepted"] = pd.to_numeric(d["Accepted"].replace("", "0"),
                                  errors="coerce").fillna(0).astype(int)
    # Transferred rows carry the SOURCE PSM's SpecId plus `_mbr_<target run>`.
    # Strip the suffix to recover the identification the feature inherits.
    d["SourceSpecId"] = d["SpecId"].str.replace(r"_mbr_.*$", "", regex=True)
    return d


def glycokey(d):
    return d["Peptide"] + "|" + d["Glycosite"] + "|" + d["Glycan"]


def matrix_of(d, runs):
    d = d.copy()
    d["key"] = glycokey(d)
    m = d.pivot_table(index="key", columns="Run", values="Area",
                      aggfunc="sum")
    return m.reindex(columns=runs)


def presence_counts(m, runs):
    pres = m.notna().sum(axis=1)
    return {f">={k}": int((pres >= k).sum()) for k in range(1, len(runs) + 1)}


def pairwise_corr_split(d, runs):
    """Cross-run r by evidence class: both identified vs at least one transfer."""
    d = d.copy()
    d["key"] = glycokey(d)
    area = d.pivot_table(index="key", columns="Run", values="Area",
                         aggfunc="sum").reindex(columns=runs)
    ev = d.pivot_table(index="key", columns="Run", values="Evidence",
                       aggfunc="first").reindex(columns=runs)
    out = []
    for i, a in enumerate(runs):
        for b in runs[i + 1:]:
            x, y = area[a], area[b]
            ok = x.notna() & y.notna() & (x > 0) & (y > 0)
            if not ok.any():
                continue
            ea = ev.loc[ok, a].to_numpy()
            eb = ev.loc[ok, b].to_numpy()
            for cat, mask in [
                    ("both identified", (ea == "Identified") & (eb == "Identified")),
                    ("one transferred", (ea == "Identified") != (eb == "Identified")),
                    ("both transferred", (ea == "Transferred") & (eb == "Transferred"))]:
                if mask.sum() < 3:
                    continue
                st = stats(log2(x[ok][mask]), log2(y[ok][mask]))
                if st:
                    out.append({"a": a, "b": b, "class": cat, **st})
    return out


def pairwise_corr(m, runs):
    rows = []
    for i, a in enumerate(runs):
        for b in runs[i + 1:]:
            x, y = m[a], m[b]
            ok = x.notna() & y.notna() & (x > 0) & (y > 0)
            st = stats(log2(x[ok]), log2(y[ok]))
            if st:
                rows.append({"a": a, "b": b, **st})
    return rows


def glycoform_consistency(d):
    """Pairwise |share_i - share_j| for site-glycan pairs seen in >= 2 runs."""
    d = d.copy()
    site = d.groupby(["Run", "Peptide", "Glycosite"])["Area"].transform("sum")
    d["share"] = np.where(site > 0, d["Area"] / site, np.nan)
    out = []
    for _, g in d.groupby(["Peptide", "Glycosite", "Glycan"]):
        s = g.groupby("Run")["share"].mean().dropna().values.astype(float)
        if len(s) < 2:
            continue
        iu = np.triu_indices(len(s), k=1)
        out.extend(np.abs(s[iu[0]] - s[iu[1]]))
    return out


def protein_rollup(d, runs):
    d = d.copy()
    d["key"] = glycokey(d)
    pep = d.groupby(["Run", "Proteins", "key"], as_index=False)["Area"].sum()
    pep = pep.sort_values(["Run", "Proteins", "Area"],
                          ascending=[True, True, False])
    pep["rn"] = pep.groupby(["Run", "Proteins"]).cumcount()
    top3 = (pep[pep["rn"] < 3].groupby(["Run", "Proteins"], as_index=False)
            ["Area"].sum().rename(columns={"Area": "top3"}))
    summed = (pep.groupby(["Run", "Proteins"], as_index=False)["Area"].sum()
              .rename(columns={"Area": "sum"}))
    j = top3.merge(summed, on=["Run", "Proteins"])
    per_run = {}
    for r in runs:
        g = j[j["Run"] == r]
        per_run[r] = stats(log2(g["top3"]), log2(g["sum"]))
    m = j.pivot_table(index="Proteins", columns="Run", values="top3")
    pairs = pairwise_corr(m.reindex(columns=runs), runs)
    return j, per_run, pairs


def sequon_factor(fasta):
    """1 + S/E by N-X-S/T sequon count, exactly as eval_entrap.py does."""
    if not fasta or not os.path.exists(fasta):
        return None
    pat = re.compile(r"N[^P][ST]")
    s = e = 0
    hdr, buf = None, []

    def flush(hdr, seq):
        nonlocal s, e
        if hdr is None or hdr.startswith(">XXX_") or hdr.startswith(">DECOY"):
            return
        n = len(pat.findall(seq))
        if "ENTRAP_" in hdr:
            e += n
        else:
            s += n

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(hdr, "".join(buf))
                hdr, buf = line.strip(), []
            else:
                buf.append(line.strip())
    flush(hdr, "".join(buf))
    if not e:
        return None
    return {"sample_sequons": s, "entrap_sequons": e, "factor": 1.0 + s / e}


def entrap_classifier(fasta):
    """(classify fn, SampleSequences or None) from the shared benchmark module."""
    if not fasta or not os.path.exists(fasta):
        def naive(accessions, peptide):
            return "entrap" if "ENTRAP_" in accessions else "sample"
        return naive, None
    sys.path.insert(0, ENTRAP_COMMON_DIR)
    from entrap_common import SampleSequences, classify  # noqa: E402
    samp = SampleSequences(fasta)

    def cls(accessions, peptide):
        return classify([accessions], peptide, samp)
    return cls, samp


def block(d, runs, classify, factor, L, add, label):
    """Matrix / consistency / protein / entrapment for one feature set."""
    add(f"### {label} — {len(d):,} features")
    add("")
    m = matrix_of(d, runs)
    p = presence_counts(m, runs)
    add("| presence | glycopeptides |")
    add("|---|---:|")
    for k, v in p.items():
        add(f"| {k} runs | {v:,} |")
    add(f"| **unique glycopeptides** | **{m.shape[0]:,}** |")
    add(f"| matrix cells filled | {int(m.notna().sum().sum()):,} "
        f"({100*m.notna().sum().sum()/max(1, m.shape[0]*len(runs)):.1f}% of grid) |")
    add("")
    pairs = pairwise_corr(m, runs)
    if pairs:
        add("| run pair | n | Pearson log2 area |")
        add("|---|---:|---:|")
        for x in pairs:
            add(f"| {x['a']} vs {x['b']} | {x['n']:,} | {x['pearson']:.3f} |")
        add("")
        add(f"Median pairwise r = "
            f"{np.median([x['pearson'] for x in pairs]):.3f}.")
        add("")
        if (d["Evidence"] == "Transferred").any():
            sp = pairwise_corr_split(d, runs)
            if sp:
                add("Same comparison split by whether a measurement was "
                    "identified or transferred in:")
                add("")
                add("| run pair | class | n | Pearson |")
                add("|---|---|---:|---:|")
                for x in sp:
                    add(f"| {x['a']} vs {x['b']} | {x['class']} | "
                        f"{x['n']:,} | {x['pearson']:.3f} |")
                add("")
    s = glycoform_consistency(d)
    if s:
        add(f"Glycoform share agreement at multi-run sites: mean |Δshare| "
            f"{np.mean(s):.4f} ({len(s):,} pairs).")
        add("")
    j, per_run, ppairs = protein_rollup(d, runs)
    if len(j):
        add(f"Proteins: {j['Proteins'].nunique():,}; median cross-run Top3 r = "
            f"{np.median([x['pearson'] for x in ppairs]):.3f}."
            if ppairs else f"Proteins: {j['Proteins'].nunique():,}.")
        add("")
    cls = d.apply(lambda r: classify(r["Proteins"], r["Peptide"]), axis=1)
    n_e = int((cls == "entrap").sum())
    n_s = int((cls != "entrap").sum())
    raw = n_e / max(1, n_e + n_s)
    add(f"Entrapment: {n_e:,} / {n_e + n_s:,} accepted features "
        f"({100*raw:.2f}% raw)" +
        (f"; sequon-corrected FDP ≈ {100*raw*factor:.2f}%"
         if factor else "; FDP not corrected (no FASTA)") + ".")
    add("")
    split = pairwise_corr_split(d, runs) if (d["Evidence"] == "Transferred").any() else []
    return {"n_features": int(len(d)), "presence": p, "split_pairs": split,
            "n_glycopeptides": int(m.shape[0]),
            "cells_filled": int(m.notna().sum().sum()),
            "grid_pct": float(100*m.notna().sum().sum()/max(1, m.shape[0]*len(runs))),
            "pairs": pairs,
            "median_pair_r": float(np.median([x["pearson"] for x in pairs]))
            if pairs else None,
            "share_mae": float(np.mean(s)) if s else None,
            "proteins": int(j["Proteins"].nunique()) if len(j) else 0,
            "protein_median_pair_r":
                float(np.median([x["pearson"] for x in ppairs])) if ppairs else None,
            "entrap": n_e, "sample": n_s, "entrap_raw": float(raw),
            "entrap_fdp": float(raw*factor) if factor else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--stem", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--lfq", default=None)
    ap.add_argument("--log", default=None)
    ap.add_argument("--psms", default=None)
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--q-id", type=float, default=Q_ID)
    args = ap.parse_args()
    base = os.path.join(args.dir, args.stem)
    lfq_path = args.lfq or base + ".glyco.lfq.tsv"
    log_path = args.log or os.path.join(args.dir, "andes5.log")
    psms_path = args.psms or os.path.join(args.dir, "out.psms")
    os.makedirs(args.out, exist_ok=True)

    d = read_lfq(lfq_path)
    runs = sorted(d["Run"].unique())
    sf = sequon_factor(args.fasta)
    factor = sf["factor"] if sf else None
    classify, sample = entrap_classifier(args.fasta)
    L = []
    add = L.append
    add("# Glycopeptide LFQ + MBR benchmark — mouse liver, 5 fractions (deep tier)")
    add("")
    add(f"Runs: {', '.join(runs)}.  Feature acceptance in the emitted tables: "
        f"feature q ≤ {Q_ACCEPT} (`Accepted`).  Identification gate applied here: "
        f"Percolator q ≤ {args.q_id} on the PSM each feature came from "
        "(for transferred rows, the source PSM).")
    add("")

    # ------------------------------------------------ identification layer
    qmap = None
    if os.path.exists(psms_path):
        p = pd.read_csv(psms_path, sep="\t")
        qmap = p.set_index("PSMId")["q-value"]
        d["QId"] = d["SourceSpecId"].map(qmap)
        add("## 0. Identification layer (pooled Percolator on the .glyco.pin)")
        add("")
        add(f"Percolator target PSMs: {len(p):,}; at q ≤ {args.q_id}: "
            f"{int((p['q-value'] <= args.q_id).sum()):,}.")
        if sample is not None:
            sys.path.insert(0, ENTRAP_COMMON_DIR)
            from entrap_common import count_percolator_psms  # noqa: E402
            add("")
            add("| Percolator gate | target PSMs | entrapment | corrected FDP |")
            add("|---|---:|---:|---:|")
            for q in (0.01, 0.05):
                n, ne, nr = count_percolator_psms(psms_path, q, sample)
                if n == 0:
                    continue
                raw = ne / n
                add(f"| q ≤ {q} | {n:,} | {ne:,} ({100*raw:.2f}%) | "
                    f"{100*raw*factor:.2f}% |" if factor else
                    f"| q ≤ {q} | {n:,} | {ne:,} ({100*raw:.2f}%) | — |")
        add("")

    # ------------------------------------------------ feature bookkeeping
    add("## 1. What the quant tables contain")
    add("")
    acc = d[d["Accepted"] == 1].copy()
    add("| run | identified rows | transferred rows | accepted identified "
        "| accepted transferred | accepted glycopeptides |")
    add("|---|---:|---:|---:|---:|---:|")
    feat_rows = []
    for r in runs:
        g = d[d["Run"] == r]
        idn = g[g["Evidence"] == "Identified"]
        trn = g[g["Evidence"] == "Transferred"]
        ia = int((idn["Accepted"] == 1).sum())
        ta = int((trn["Accepted"] == 1).sum())
        gp = g[g["Accepted"] == 1].pipe(glycokey).nunique()
        add(f"| {r} | {len(idn):,} | {len(trn):,} | {ia:,} | {ta:,} | {gp:,} |")
        feat_rows.append({"run": r, "identified": len(idn),
                          "transferred": len(trn), "acc_identified": ia,
                          "acc_transferred": ta, "glycopeptides": gp})
    add("")
    topo = acc.groupby(["Peptide", "Glycosite", "Glycan", "Topology"]).ngroups
    multi = (acc.pipe(glycokey).to_frame("k").assign(
        t=acc["Topology"].values).groupby("k")["t"].nunique() > 1).sum()
    add(f"Accepted: **{len(acc):,} features**, "
        f"**{acc.pipe(glycokey).nunique():,} glycopeptides** "
        f"(peptide, site, composition), **{topo:,}** at topology resolution; "
        f"{int(multi):,} compositions carry more than one retained topology.")
    add("")

    # ------------------------------------------------ the headline check
    add("## 2. Identification control: feature q is not identification q")
    add("")
    entrap_by_evidence = {}
    if qmap is not None:
        add("`Accepted` is a measurement/XIC q-value.  Joining each feature to "
            "the Percolator q-value of its PSM shows what it lets through:")
        add("")
        add("| feature set | features | entrapment raw | sequon-corrected FDP |")
        add("|---|---:|---:|---:|")
        for label, sel in [
                ("accepted, no ID gate", acc),
                (f"accepted + ID q ≤ {args.q_id}", acc[acc["QId"] <= args.q_id]),
                ("accepted + ID q > 0.05 / no q",
                 acc[~(acc["QId"] <= 0.05)])]:
            if not len(sel):
                continue
            cls = sel.apply(lambda r: classify(r["Proteins"], r["Peptide"]), axis=1)
            n_e = int((cls == "entrap").sum())
            raw = n_e / len(sel)
            add(f"| {label} | {len(sel):,} | {100*raw:.1f}% | "
                f"{100*raw*factor:.1f}% |" if factor else
                f"| {label} | {len(sel):,} | {100*raw:.1f}% | — |")
        add("")
        add("Entrapment by evidence class:")
        add("")
        add("| evidence | gate | features | entrapment |")
        add("|---|---|---:|---:|")
        entrap_by_evidence = {}
        for ev in ("Identified", "Transferred"):
            for gl, sel in [("no ID gate", acc[acc["Evidence"] == ev]),
                            (f"q ≤ {args.q_id}",
                             acc[(acc["Evidence"] == ev) &
                                 (acc["QId"] <= args.q_id)])]:
                cls = sel.apply(
                    lambda r: classify(r["Proteins"], r["Peptide"]), axis=1)
                n_e = int((cls == "entrap").sum())
                add(f"| {ev} | {gl} | {len(sel):,} | "
                    f"{n_e:,} ({100*n_e/max(1,len(sel)):.2f}%) |")
                entrap_by_evidence.setdefault(ev, {})[gl] = {
                    "n": int(len(sel)), "entrap": n_e,
                    "raw": float(n_e / max(1, len(sel)))}
        add("")
        drop = 100 * (1 - (acc["QId"] <= args.q_id).mean())
        add(f"Only **{int((acc['QId'] <= args.q_id).sum()):,}** of "
            f"{len(acc):,} accepted features ({100*(acc['QId'] <= args.q_id).mean():.1f}%) "
            f"pass the 1% identification gate; {drop:.0f}% of what the tables "
            "report as quantified glycopeptides come from PSMs that are mostly "
            "random rank-1 matches.  All sections after this one are computed "
            "**after** the identification gate.")
        add("")
    else:
        add("No Percolator `out.psms` supplied — sections below use the emitted "
            "acceptance and are NOT identification-FDR controlled.")
        add("")
    gated = acc[acc["QId"] <= args.q_id].copy() if qmap is not None else acc.copy()
    g_id = gated[gated["Evidence"] == "Identified"]
    g_tr = gated[gated["Evidence"] == "Transferred"]
    add(f"Gated features: {len(gated):,} "
        f"({len(g_id):,} identified + {len(g_tr):,} transferred).")
    add("")
    add("| run | identified | transferred | unique glycopeptides (composition) "
        "| at topology resolution |")
    add("|---|---:|---:|---:|---:|")
    for r in runs:
        gg = gated[gated["Run"] == r]
        add(f"| {r} | {int((gg['Evidence'] == 'Identified').sum()):,} | "
            f"{int((gg['Evidence'] == 'Transferred').sum()):,} | "
            f"{gg.pipe(glycokey).nunique():,} | "
            f"{gg.groupby(['Peptide', 'Glycosite', 'Glycan', 'Topology']).ngroups:,} |")
    add("")
    topo_g = gated.groupby(["Peptide", "Glycosite", "Glycan", "Topology"]).ngroups
    multi_g = int((gated.assign(k=glycokey(gated)).groupby("k")["Topology"]
                   .nunique() > 1).sum())
    add(f"Gated: **{gated.pipe(glycokey).nunique():,}** glycopeptides at composition "
        f"resolution, **{topo_g:,}** at topology resolution; "
        f"{multi_g:,} compositions carry more than one retained topology.")
    add("")

    # ------------------------------------------------ transfer stats
    add("## 3. Transfers (MBR)")
    add("")
    tr_all = d[d["Evidence"] == "Transferred"]
    if len(tr_all):
        qs = qstats(pd.to_numeric(tr_all.loc[tr_all["Accepted"] == 1,
                                             "FeatureQValue"], errors="coerce"))
        add(f"Accepted transfers (feature q ≤ {Q_ACCEPT}): {len(g_tr):,} after the "
            f"ID gate / {int((tr_all['Accepted'] == 1).sum()):,} before it.")
        if qs.get("n"):
            add(f"Transfer feature q: median {qs['median']:.4f}, "
                f"p90 {qs['p90']:.4f}, max {qs['max']:.4f}.")
        add("")

    # ------------------------------------------------ gated blocks
    add("## 4. Quantification matrix, cross-run consistency, proteins")
    add("")
    res = {}
    res["all_accepted"] = block(acc, runs, classify, factor, L, add,
                                "All accepted (as emitted, no ID gate)")
    res["gated_identified"] = block(g_id, runs, classify, factor, L, add,
                                    "Identified only, ID-gated")
    res["gated_with_mbr"] = block(gated, runs, classify, factor, L, add,
                                  "Identified + MBR, ID-gated")
    if qmap is not None:
        res["ungated_with_mbr"] = res["all_accepted"]

    # ------------------------------------------------ log
    logtxt = open(log_path).read() if os.path.exists(log_path) else ""
    add("## 5. Run log (quant / MBR lines)")
    add("")
    for ln in logtxt.splitlines():
        if re.search(r"glyco-quant|glyco-mbr|Quantification:|provenance", ln):
            add(f"- `{ln.strip()}`")
    add("")

    report = "\n".join(L)
    with open(os.path.join(args.out, "glyco_quant_report.md"), "w") as fh:
        fh.write(report + "\n")
    print(report)

    out = {"runs": runs, "q_id": args.q_id, "q_accept": Q_ACCEPT,
           "sequon_space": sf, "features": feat_rows,
           "gated_counts": {"identified": int(len(g_id)),
                            "transferred": int(len(g_tr))},
           "topology": {"compositions": int(acc.pipe(glycokey).nunique()),
                        "topologies": int(topo),
                        "multi_topology_compositions": int(multi)},
           "topology_gated": {"compositions": int(gated.pipe(glycokey).nunique()),
                              "topologies": int(topo_g),
                              "multi_topology_compositions": int(multi_g)},
           "entrap_by_evidence": entrap_by_evidence,
           "blocks": res}
    with open(os.path.join(args.out, "glyco_quant_numbers.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
