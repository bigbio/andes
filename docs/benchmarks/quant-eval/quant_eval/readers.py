"""Input readers, one per tool. Every reader returns a run-tagged table with the columns
key (andes key form), run (design run name), quantity (RAW, no normalisation of any kind), is_transfer, library_q, tq, thr.

  andes     finalize parquet: `quantity` (raw), `quant_q`, `transfer_q`, `library_q`, `evidence_type`
  MaxQuant  evidence.txt `Intensity`, summed over the evidence rows of a precursor x run; reverse and contaminants removed
  FragPipe  IonQuant combined_ion.tsv `<run> Intensity`, from a run with IonQuant normalisation off (`--normalization 0`)
"""
import csv
import math
import os
import sys

import numpy as np
import pandas as pd

from . import runguard as G
from .keys import fragpipe_key, maxquant_key, norm_key

IDQ = 0.01          # library_q limit applied to every andes row
TRANSFER_Q = 0.01   # transfer_q limit of the `default` rule
ANDES_RULES = ("default", "all", "hard-cut", "weights")


def _fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def _read_tsv(path):
    """Rows of a TSV as dicts (read fully, file closed) and its header."""
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        return list(r), r.fieldnames


def remap_runs(d, runs):
    """Map a tool's run names onto the design run names ('-' and '_' and letter case are ignored); unmatched rows are dropped."""
    norm = {r.replace("-", "_").lower(): r for r in runs}
    d = d.copy()
    d["run"] = d.run.map(lambda r: norm.get(str(r).replace("-", "_").lower()))
    n0 = len(d)
    d = d[d.run.notna()]
    print(f"remap {n0} -> {len(d)} rows on {d.run.nunique()} runs", file=sys.stderr)
    return d


def _finish(d, tag):
    d = d[(d.quantity > 0) & np.isfinite(d.quantity)].drop_duplicates(["key", "run"]).copy()
    d["library_q"] = 0.0
    # a tool that reports a per-transfer q-value (FlashLFQ: PIP q-value) passes it as the column pq
    d["tq"] = d["pq"].values if "pq" in d.columns else 0.0
    d["thr"] = d["tq"].values
    return G.tagged(d[["key", "run", "quantity", "is_transfer", "library_q", "tq", "thr"]], tag)


def read_maxquant(path, runs):
    r, cols = _read_tsv(path)
    for c in ("Modified sequence", "Charge", "Raw file", "Type", "Intensity"):
        if c not in cols:
            raise SystemExit(f"{path}: missing column {c!r}")
    # MULTI-MATCH evidence rows of some versions carry no 'Modified sequence': take it from 'Mod. peptide ID'
    modseq_of = {}
    if "Mod. peptide ID" in cols:
        for row in r:
            ms = (row.get("Modified sequence") or "").strip()
            if ms:
                modseq_of.setdefault(row["Mod. peptide ID"], ms)
    groups, drop = {}, {"reverse": 0, "contaminant": 0, "no_intensity": 0, "no_sequence": 0, "unmapped_mod": 0}
    for row in r:
        if not (row.get("Modified sequence") or "").strip():
            ms = modseq_of.get(row.get("Mod. peptide ID"))
            if not ms:
                drop["no_sequence"] += 1
                continue
            row["Modified sequence"] = ms
        if (row.get("Reverse") or "").strip() == "+":
            drop["reverse"] += 1
            continue
        if (row.get("Potential contaminant") or row.get("Contaminant") or "").strip() == "+":
            drop["contaminant"] += 1
            continue
        q = _fnum(row["Intensity"])
        if not (q > 0):
            drop["no_intensity"] += 1
            continue
        try:
            key = maxquant_key(row["Modified sequence"], int(float(row["Charge"])))
        except ValueError:
            drop["unmapped_mod"] += 1
            continue
        g = groups.setdefault((key, row["Raw file"]), [0.0, True])
        g[0] += q
        if not row["Type"].strip().endswith(("MULTI-MATCH", "MULTI-MATCH-MSMS")):
            g[1] = False
    print(f"maxquant {path}: precursor x run groups {len(groups)} dropped {drop}", file=sys.stderr)
    d = pd.DataFrame([(k, rn, g[0], g[1]) for (k, rn), g in groups.items()], columns=["key", "run", "quantity", "is_transfer"])
    return _finish(remap_runs(d, runs), G.run_tag("maxquant-evidence-Intensity", path))


def read_fragpipe(path, runs):
    r, cols = _read_tsv(path)
    found = sorted({c[: -len(" Intensity")] for c in cols if c.endswith(" Intensity") and not c.endswith("MaxLFQ Intensity")
                    and c[: -len(" Intensity")] + " Match Type" in cols})
    if not found:
        raise SystemExit(f"{path}: no '<run> Intensity' + '<run> Match Type' column pairs")
    rows, bad = {}, 0
    for row in r:
        try:
            key = fragpipe_key(row["Peptide Sequence"], row["Assigned Modifications"], int(float(row["Charge"])))
        except ValueError:
            bad += 1
            continue
        for x in found:
            q = _fnum(row[x + " Intensity"])
            if not (q > 0):
                continue
            o = rows.get((key, x))
            if o is None or o[0] < q:
                rows[(key, x)] = (q, (row[x + " Match Type"] or "").strip() == "MBR")
    print(f"fragpipe {path}: runs {len(found)} unmapped-modification ions dropped {bad}", file=sys.stderr)
    d = pd.DataFrame([(k, x, q, t) for (k, x), (q, t) in rows.items()], columns=["key", "run", "quantity", "is_transfer"])
    return _finish(remap_runs(d, runs), G.run_tag("fragpipe-ionquant-norm0", path))


def read_flashlfq(path, runs, pip_q=None):
    """FlashLFQ (PIP-ECHO MBR FDR) `QuantifiedPeaks.tsv` -> rows of one precursor key x run.
    Key = `Full Sequence` + `/` + `Precursor Charge` (the identification charge of the peak); the file's `Full Sequence` is the andes key form
    when the identifications were written from andes search tables. Decoy-peptide and random-RT peaks are dropped. `Peak Detection Type`
    MBR = transfer; its `PIP Q-Value` is kept in tq / thr (0 for MS/MS peaks); `pip_q` drops MBR peaks above that PIP q-value.
    FlashLFQ sums the charge states of a peptide into one peak intensity unless it is run with charge-specific quantification."""
    r, cols = _read_tsv(path)
    for c in ("File Name", "Full Sequence", "Precursor Charge", "Peak intensity", "Peak Detection Type"):
        if c not in cols:
            raise SystemExit(f"{path}: missing column {c!r}")
    rows, drop = {}, {"decoy": 0, "random_rt": 0, "no_intensity": 0, "above_pip_q": 0}
    for row in r:
        if (row.get("Decoy Peptide") or "").strip().lower() == "true":
            drop["decoy"] += 1
            continue
        if (row.get("Random RT") or "").strip().lower() == "true":
            drop["random_rt"] += 1
            continue
        q = _fnum(row["Peak intensity"])
        if not (q > 0):
            drop["no_intensity"] += 1
            continue
        is_tr = row["Peak Detection Type"].strip().upper() == "MBR"
        pq = _fnum(row.get("PIP Q-Value")) if is_tr else 0.0
        if is_tr and pip_q is not None and not (pq <= pip_q):
            drop["above_pip_q"] += 1
            continue
        key = norm_key(row["Full Sequence"] + "/" + str(int(float(row["Precursor Charge"]))))
        k = (key, row["File Name"])
        o = rows.get(k)
        if o is None or o[0] < q:
            rows[k] = (q, is_tr, pq if pq == pq else 0.0)
    print(f"flashlfq {path}: rows {len(rows)} dropped {drop}", file=sys.stderr)
    d = pd.DataFrame([(k, rn, v[0], v[1], v[2]) for (k, rn), v in rows.items()], columns=["key", "run", "quantity", "is_transfer", "pq"])
    return _finish(remap_runs(d, runs), G.run_tag("flashlfq-quantified-peaks", path))


def recover_p(tr):
    """Per-quantity error probability p recovered from quant_q = running mean of p over a run's accepted transfers sorted by p.
    tr: transfer rows with run and quant_q (unique index). p_i = i * q_i - (i - 1) * q_(i-1) within a run.
    The result is indexed like `tr` (input order), not in the sorted order used internally."""
    order = tr.index
    tr = tr.sort_values(["run", "quant_q"], kind="stable")
    i = tr.groupby("run").cumcount().values + 1
    qq = tr.quant_q.values
    prev = np.r_[0.0, qq[:-1]]
    prev[i == 1] = 0.0
    return pd.Series(np.clip(i * qq - (i - 1) * prev, 0.0, 1.0), index=tr.index).reindex(order)


def read_andes(path, runs, rule="default", hard_cut=0.01, library_q=IDQ):
    """andes finalize table -> rows of one view.
      default    identified rows + transfers with transfer_q <= 0.01
      all        identified rows + every transfer the finalize output carries
      hard-cut   identified rows + transfers with quant_q <= hard_cut
      weights    as `all`, plus the column p (recovered error probability of each transfer) for the precision-weight view
    Rows need library_q <= `library_q` (default 0.01; a run with a looser library entry gate passes its own limit) and a positive
    finite quantity; below-threshold rows are dropped."""
    if rule not in ANDES_RULES:
        raise SystemExit(f"andes rule {rule!r} not in {ANDES_RULES}")
    m = pd.read_parquet(path, columns=["precursor_key", "run", "evidence_type", "quantity", "quant_q", "transfer_q", "library_q"])
    tag = G.run_tag(f"andes-finalize:{rule}", path)
    m = m[np.isfinite(m.quantity) & (m.quantity > 0) & (m.evidence_type != "below_threshold")].copy()
    m["key"] = [norm_key(k) for k in m.precursor_key.values]
    m["is_transfer"] = ~m.evidence_type.str.startswith("identified")
    m = remap_runs(m.drop_duplicates(["key", "run"]), runs)
    m["thr"] = m.quant_q
    m["tq"] = np.where(m.is_transfer, m.transfer_q, 0.0)
    m = m[m.library_q <= library_q]
    if rule == "default":
        m = m[(~m.is_transfer) | (m.tq <= TRANSFER_Q)]
    elif rule == "hard-cut":
        m = m[(~m.is_transfer) | (m.quant_q <= hard_cut)]
    out = m[["key", "run", "quantity", "is_transfer", "library_q", "tq", "thr"]].copy()
    if rule == "weights":
        tr = m[m.is_transfer]
        ptab = G.tagged(pd.DataFrame({"key": tr.key.values, "run": tr.run.values, "p": recover_p(tr).values}), tag)
        out = G.join(G.tagged(out, tag), ptab, ["key", "run"], how="left")
    else:
        out["p"] = np.nan
        out = G.tagged(out, tag)
    return out


def read_andes_rt(path, runs):
    """Retention-time proxy of each precursor key: median apex RT (`apex_rt_s`) over the identified rows of an andes finalize table
    (library_q <= 0.01). One proxy for every arm of the RT-local normalisation; keys without an andes row have none."""
    m = pd.read_parquet(path, columns=["precursor_key", "run", "evidence_type", "apex_rt_s", "library_q"])
    m = m[m.evidence_type.str.startswith("identified") & (m.library_q <= IDQ) & np.isfinite(m.apex_rt_s)]
    m = remap_runs(m, runs)
    m["key"] = [norm_key(k) for k in m.precursor_key.values]
    return m.groupby("key").apex_rt_s.median()
