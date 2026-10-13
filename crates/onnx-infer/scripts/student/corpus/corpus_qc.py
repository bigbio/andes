#!/usr/bin/env python3
"""Quality control for training-corpus PSM sources.

Three input kinds, one set of metrics (written per PSM, per run and per dataset):

``flat``      MsNet-format flats (``seq, charge, prec_mz, res_mod_pos, res_mod_delta,
              nterm_delta, cterm_delta, mz[], intensity[]``), one parquet per run,
              file name ``<DATASET>__<run>.parquet`` (or pass ``--dataset``).
              Spectra are present, so all spectrum metrics are computed.
``psm_table`` qpx ``psm.parquet`` (ProForma peptidoform, charge, PEP, q-value,
              is_decoy, run_file_name, scan, rt). No spectra: score, decoy and
              composition metrics only.
``shard``     andes-ms2 shards (already annotated): annotation metrics re-aggregated.

Metrics
-------
Spectrum (flat):
  * ``n_peaks``, ``tic``                              peaks / total intensity per spectrum
  * ``ann_frac_ft`` / ``ann_frac_it``                 annotated-intensity fraction of plain
    b/y 1+/2+ at 20 ppm and at 0.5 Da (each peak counted once)
  * ``ann_frac``                                      the same at the tolerance of the
    *inferred* analyser (below) - the headline number
  * ``n_by``, ``by_coverage``                         matched b/y 1+ ions / (n - 1) x 2
  * ``frag_err_ppm`` / ``frag_err_da``                median signed error of the b/y 1+
    ions matched in a 0.5 Da window (the most intense peak per ion)
  * ``frag_abs_err_da``                               median |error| (analyser inference:
    a run whose median |error| is < 0.01 Da is FTMS-like, > 0.05 Da ion-trap-like)
Precursor (flat):
  * ``prec_err_ppm``, ``prec_iso``                    best isotope offset in -1..3 at the
    stated charge; ``prec_fit`` = |error| <= ``--prec-tol-ppm``
  * ``charge_alt_fit``                                the precursor fits *better* at z +/- 1
    (charge-state inconsistency)
Composition (all inputs):
  * length, charge, PTM inventory (Unimod via ``mod_composition``; unmapped deltas
    kept as ``mass:<delta>``), ``has_I``, tryptic C-terminus, missed cleavages
  * isoleucine handling: per dataset fraction of sequences with I; 0 % with L present
    means I was written as L (the qpx convention) and the I/L-merged key must be used
  * cross-source duplicates: ``peptidoforms`` subcommand writes the per-dataset set of
    I/L-merged peptidoform keys; ``overlap`` compares any number of them
Score (psm_table): q-value / PEP quantiles, decoy share, target/decoy score separation
  (AUC of the best available score column) when decoys are present.
Run-level outliers: robust z (median / 1.4826 MAD over the runs of a dataset) of PSM
  count, median ann_frac, peaks, precursor and fragment error; |z| > ``--outlier-z``
  is flagged.

Usage::

  python corpus_qc.py flat --out-dir QC/PXD009875 --dataset PXD009875 --workers 8 FLAT.parquet...
  python corpus_qc.py psm_table --out-dir QC/PXD012345 psm.parquet
  python corpus_qc.py peptidoforms --out keys.parquet QC/*/psm_qc.parquet
  python corpus_qc.py overlap keysA.parquet keysB.parquet ...
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "ms2"))
from mod_composition import MASSES, UNIMOD, match_mass  # noqa: E402
from ms2_common import H2O, PROTON, RESIDUE_MASS  # noqa: E402

C13 = 1.0033548378
CAM = 4
FT_TOL_PPM = 20.0
IT_TOL_DA = 0.5


# ----------------------------------------------------------------------------- helpers
def parse_flat_mods(seq, pos, delta, nterm, cterm):
    """-> list of (pos, label) with label = Unimod id or 'mass:+x.xxxx'; pos 0 = N-term."""
    mods = []
    for p, d in zip(pos or [], delta or []):
        d = float(d)
        r = seq[p - 1] if 1 <= p <= len(seq) else None
        uid = match_mass(d, r, False)
        if uid is None and p == 1:
            uid = match_mass(d, None, True)
            if uid is not None:
                mods.append((0, uid, d))
                continue
        mods.append((p, uid if uid is not None else f"mass:{d:+.4f}", d))
    for term, d in ((0, nterm), (len(seq) + 1, cterm)):
        if d and abs(d) > 1e-6:
            uid = match_mass(float(d), None, term == 0)
            mods.append((term, uid if uid is not None else f"mass:{float(d):+.4f}", float(d)))
    return mods


def residue_masses(seq, mods):
    m = np.array([RESIDUE_MASS.get(c, np.nan) for c in seq], np.float64)
    n = len(seq)
    for pos, _, d in mods:
        m[min(max(pos, 1), n) - 1] += d
    return m


def prec_fit(res, z, mz):
    """best (error ppm, isotope) at charge z."""
    theo = (res.sum() + H2O + z * PROTON) / z
    best = (np.inf, 0)
    for k in range(-1, 4):
        e = (mz - (theo + k * C13 / z)) / theo * 1e6
        if abs(e) < abs(best[0]):
            best = (e, k)
    return best


def by_mz(res):
    """(n-1, 4) b1 b2 y1 y2."""
    pre = np.cumsum(res)[:-1]
    tot = res.sum()
    b, y = pre, tot - pre + H2O
    return np.stack([(b + PROTON), (b + 2 * PROTON) / 2, (y + PROTON), (y + 2 * PROTON) / 2], axis=1)


def ann_fraction(theo, valid, mz, inten, tol, ppm):
    w = theo * tol * 1e-6 if ppm else np.full_like(theo, tol)
    lo = np.searchsorted(mz, theo - w, "left")
    hi = np.searchsorted(mz, theo + w, "right")
    used = np.zeros(len(mz), bool)
    hit = valid & (hi > lo)
    for a, b in zip(lo[hit], hi[hit]):
        used[a:b] = True
    tot = inten.sum()
    return (float(inten[used].sum() / tot) if tot > 0 else 0.0), hit


def is_tryptic(seq):
    return seq[-1] in "KR"


def missed_cleavages(seq):
    return sum(1 for i, c in enumerate(seq[:-1]) if c in "KR" and seq[i + 1] != "P")


def il_key(peptidoform: str) -> str:
    return peptidoform.replace("I", "L")


def mod_label(uid):
    return f"UNIMOD:{uid}" if isinstance(uid, int) else uid


def peptidoform_of(seq, mods):
    tag = collections.defaultdict(list)
    for pos, uid, _ in mods:
        tag[pos].append(f"[{mod_label(uid)}]")
    s = "".join(tag.get(0, []))
    s = s + "-" if s else ""
    for i, c in enumerate(seq, 1):
        s += c + "".join(tag.get(i, []))
    ct = "".join(tag.get(len(seq) + 1, []))
    return s + ("-" + ct if ct else "")


# ----------------------------------------------------------------------------- flats
def qc_flat_file(args):
    path, dataset, prec_tol, limit = args
    run = os.path.basename(path).rsplit(".parquet", 1)[0]
    if "__" in run:
        ds, run = run.split("__", 1)
        dataset = dataset or ds
    cols = ["seq", "charge", "prec_mz", "res_mod_pos", "res_mod_delta", "nterm_delta", "cterm_delta", "mz", "intensity"]
    out = collections.defaultdict(list)
    n = 0
    for batch in pq.ParquetFile(path).iter_batches(batch_size=2048, columns=cols):
        d = batch.to_pydict()
        for i in range(batch.num_rows):
            n += 1
            if limit and n > limit:
                break
            seq, z, pmz = d["seq"][i], int(d["charge"][i]), float(d["prec_mz"][i])
            mods = parse_flat_mods(seq, d["res_mod_pos"][i], d["res_mod_delta"][i], d["nterm_delta"][i], d["cterm_delta"][i])
            res = residue_masses(seq, mods)
            ok_aa = not np.isnan(res).any()
            added_cam = False
            if ok_aa:
                e, k = prec_fit(res, z, pmz)
                if abs(e) > prec_tol and "C" in seq:
                    cam = [(j + 1, CAM, MASSES[CAM]) for j, c in enumerate(seq)
                           if c == "C" and not any(p == j + 1 for p, _, _ in mods)]
                    res2 = residue_masses(seq, mods + cam)
                    e2, k2 = prec_fit(res2, z, pmz)
                    if abs(e2) <= prec_tol:
                        mods, res, e, k, added_cam = sorted(mods + cam, key=lambda x: x[0]), res2, e2, k2, True
                alt = min((abs(prec_fit(res, zz, pmz)[0]) for zz in (z - 1, z + 1) if zz >= 1), default=np.inf)
            else:
                e, k, alt = np.nan, 0, np.inf
            mz = np.asarray(d["mz"][i], np.float64)
            it = np.asarray(d["intensity"][i], np.float64)
            o = np.argsort(mz, kind="stable")
            mz, it = mz[o], it[o]
            L = len(seq)
            if ok_aa and L >= 2 and len(mz):
                theo = by_mz(res)
                valid = np.ones_like(theo, bool)
                valid[:, [1, 3]] = z >= 2
                f_ft, hit_ft = ann_fraction(theo, valid, mz, it, FT_TOL_PPM, True)
                f_it, hit_it = ann_fraction(theo, valid, mz, it, IT_TOL_DA, False)
                # fragment error of 1+ ions matched in 0.5 Da (most intense peak in window)
                t1 = theo[:, [0, 2]].ravel()
                lo = np.searchsorted(mz, t1 - IT_TOL_DA, "left")
                hi = np.searchsorted(mz, t1 + IT_TOL_DA, "right")
                errs, tm = [], []
                for t, a, b in zip(t1, lo, hi):
                    if b > a:
                        j = a + int(np.argmax(it[a:b]))
                        errs.append(mz[j] - t)
                        tm.append(t)
                errs, tm = np.asarray(errs), np.asarray(tm)
                n_by_ft = int(hit_ft[:, [0, 2]].sum())
                n_by_it = len(errs)
                fe_da = float(np.median(errs)) if len(errs) else np.nan
                fae_da = float(np.median(np.abs(errs))) if len(errs) else np.nan
                fe_ppm = float(np.median(errs / tm * 1e6)) if len(errs) else np.nan
            else:
                f_ft = f_it = fe_da = fae_da = fe_ppm = np.nan
                n_by_ft = n_by_it = 0
            unm = [mod_label(u) for _, u, _ in mods if not isinstance(u, int)]
            out["dataset"].append(dataset or "")
            out["run"].append(run)
            out["row"].append(n - 1)
            out["sequence"].append(seq)
            out["peptidoform"].append(peptidoform_of(seq, mods))
            out["mods"].append("|".join(sorted({mod_label(u) for _, u, _ in mods})))
            out["n_unmapped_mods"].append(len(unm))
            out["added_cam"].append(added_cam)
            out["length"].append(L)
            out["charge"].append(z)
            out["has_I"].append("I" in seq)
            out["has_L"].append("L" in seq)
            out["tryptic_cterm"].append(is_tryptic(seq))
            out["missed_cleavages"].append(missed_cleavages(seq))
            out["prec_err_ppm"].append(float(e))
            out["prec_iso"].append(int(k))
            out["prec_fit"].append(bool(abs(e) <= prec_tol))
            out["charge_alt_fit"].append(bool(alt < abs(e) and alt <= prec_tol))
            out["n_peaks"].append(len(mz))
            out["tic"].append(float(it.sum()))
            out["ann_frac_ft"].append(f_ft)
            out["ann_frac_it"].append(f_it)
            out["n_by_ft"].append(n_by_ft)
            out["n_by_it"].append(n_by_it)
            out["by_coverage_ft"].append(n_by_ft / max(1, 2 * (L - 1)))
            out["frag_err_da"].append(fe_da)
            out["frag_abs_err_da"].append(fae_da)
            out["frag_err_ppm"].append(fe_ppm)
        if limit and n > limit:
            break
    return pa.table(dict(out))


def infer_analyser(frag_abs_err_da_median, ann_ft, ann_it):
    """FT if fragments are accurate to < 0.01 Da, IT if > 0.05 Da or 20 ppm loses most of
    the 0.5 Da annotation; 'mixed?' otherwise."""
    if np.isnan(frag_abs_err_da_median):
        return "unknown"
    if frag_abs_err_da_median < 0.01 and ann_ft >= 0.8 * ann_it:
        return "FTMS"
    if frag_abs_err_da_median > 0.05 or ann_ft < 0.5 * ann_it:
        return "ITMS"
    return "mixed?"


def q(x, p):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    return float(np.quantile(x, p)) if len(x) else float("nan")


def robust_z(v):
    v = np.asarray(v, float)
    med = np.nanmedian(v)
    mad = np.nanmedian(np.abs(v - med)) * 1.4826
    return (v - med) / mad if mad > 0 else np.zeros_like(v)


def summarise(df, z_thr):
    """run-level table + dataset summary dict from the PSM-level frame (flats)."""
    runs = []
    for run, g in df.groupby("run", sort=True):
        fa = g.frag_abs_err_da.median()
        an = infer_analyser(fa, g.ann_frac_ft.median(), g.ann_frac_it.median())
        ann = g.ann_frac_ft if an == "FTMS" else g.ann_frac_it
        runs.append({
            "dataset": g.dataset.iloc[0], "run": run, "n_psms": len(g), "analyser_inferred": an,
            "ann_frac_median": float(ann.median()), "ann_frac_ft_median": float(g.ann_frac_ft.median()),
            "ann_frac_it_median": float(g.ann_frac_it.median()), "n_peaks_median": float(g.n_peaks.median()),
            "by_coverage_ft_median": float(g.by_coverage_ft.median()),
            "frag_err_da_median": float(g.frag_err_da.median()), "frag_abs_err_da_median": float(fa),
            "frag_err_ppm_median": float(g.frag_err_ppm.median()),
            "prec_err_ppm_median": float(g.prec_err_ppm.median()),
            "prec_err_ppm_iqr": q(g.prec_err_ppm, .75) - q(g.prec_err_ppm, .25),
            "prec_fit_frac": float(g.prec_fit.mean()), "prec_iso_nonzero_frac": float((g.prec_iso != 0).mean()),
            "charge_alt_fit_frac": float(g.charge_alt_fit.mean()),
            "added_cam_frac": float(g.added_cam.mean()), "unmapped_mod_frac": float((g.n_unmapped_mods > 0).mean()),
            "tryptic_frac": float(g.tryptic_cterm.mean()), "length_median": float(g.length.median()),
            "charge_mode": int(g.charge.mode().iloc[0]),
            "low_ann_frac": float((ann < 0.2).mean()), "unique_peptidoforms": int(g.peptidoform.nunique()),
        })
    import pandas as pd

    R = pd.DataFrame(runs)
    # robust z within the runs of the same inferred analyser: a mixed FTMS/ITMS dataset is
    # reported as mixed (dataset level), not as a set of outliers
    zcols = ("n_psms", "ann_frac_median", "n_peaks_median", "prec_err_ppm_median", "frag_err_da_median")
    for col in zcols:
        R[f"z_{col}"] = R.groupby("analyser_inferred")[col].transform(lambda v: pd.Series(robust_z(v), index=v.index))
    flags = []
    for _, r in R.iterrows():
        f = [c[2:] for c in R.columns if c.startswith("z_") and abs(r[c]) > z_thr]
        if r.analyser_inferred in ("mixed?", "unknown"):
            f.append("analyser_" + r.analyser_inferred.strip("?"))
        if r.prec_fit_frac < 0.95:
            f.append("precursor_join")
        flags.append(",".join(f))
    R["outlier_flags"] = flags
    return R


def composition_summary(df, seq_col="sequence", pf_col="peptidoform", mods_col="mods"):
    s = df[seq_col]
    has_I = s.str.contains("I").mean()
    has_L = s.str.contains("L").mean()
    mods = collections.Counter()
    for m in df[mods_col]:
        for x in str(m).split("|"):
            if x:
                mods[x] += 1
    lens = df[seq_col].str.len()
    return {
        "n_psms": int(len(df)),
        "unique_sequences": int(s.nunique()),
        "unique_peptidoforms": int(df[pf_col].nunique()),
        "unique_peptidoforms_IL_merged": int(df[pf_col].str.replace("I", "L").nunique()),
        "frac_seq_with_I": float(has_I), "frac_seq_with_L": float(has_L),
        "isoleucine_written_as_L": bool(has_I == 0 and has_L > 0),
        "length_quantiles": {p: q(lens, p) for p in (0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)},
        "length_hist": {int(k): int(v) for k, v in lens.clip(upper=40).value_counts().sort_index().items()},
        "charge_hist": {int(k): int(v) for k, v in df.charge.value_counts().sort_index().items()},
        "tryptic_cterm_frac": float(s.str[-1].isin(["K", "R"]).mean()),
        "mod_psm_counts": dict(mods.most_common()),
        "mod_names": {k: UNIMOD[int(k.split(":")[1])][0] for k in mods if k.startswith("UNIMOD:")
                      and int(k.split(":")[1]) in UNIMOD},
    }


def deamid_artifact(df):
    """Deamidation (+0.984) on a precursor picked one isotope low (-1.003) is the unmodified
    peptide at its monoisotopic peak: the classic mis-assignment. Report how many deamidated
    PSMs need isotope -1 (or +1 with two deamidations) vs the rate on unmodified PSMs."""
    dm = df.mods.str.contains("UNIMOD:7", regex=False)
    if not dm.any():
        return {"deamidated_psms": 0}
    return {
        "deamidated_psms": int(dm.sum()), "deamidated_frac": float(dm.mean()),
        "iso_minus1_frac_deamidated": float((df.prec_iso[dm] == -1).mean()),
        "iso_minus1_frac_other": float((df.prec_iso[~dm] == -1).mean()),
        "iso_nonzero_frac_deamidated": float((df.prec_iso[dm] != 0).mean()),
        "iso_nonzero_frac_other": float((df.prec_iso[~dm] != 0).mean()),
    }


def cmd_flat(a):
    import pandas as pd

    os.makedirs(a.out_dir, exist_ok=True)
    jobs = [(p, a.dataset, a.prec_tol_ppm, a.limit) for p in a.inputs]
    with ProcessPoolExecutor(a.workers) as ex:
        tables = list(ex.map(qc_flat_file, jobs))
    t = pa.concat_tables(tables)
    pq.write_table(t, os.path.join(a.out_dir, "psm_qc.parquet"), compression="zstd")
    df = t.to_pandas()
    R = summarise(df, a.outlier_z)
    R.to_csv(os.path.join(a.out_dir, "run_qc.tsv"), sep="\t", index=False, float_format="%.5g")
    S = composition_summary(df)
    per_an = {}
    for an, g in R.groupby("analyser_inferred"):
        sub = df[df.run.isin(g.run)]
        ann = sub.ann_frac_ft if an == "FTMS" else sub.ann_frac_it
        per_an[an] = {
            "runs": int(len(g)), "psms": int(len(sub)),
            "ann_frac_quantiles": {p: q(ann, p) for p in (0.1, 0.25, 0.5, 0.75, 0.9)},
            "n_peaks_quantiles": {p: q(sub.n_peaks, p) for p in (0.1, 0.5, 0.9)},
            "frag_err_da_quantiles": {p: q(sub.frag_err_da, p) for p in (0.05, 0.5, 0.95)},
            "frag_abs_err_da_median": q(sub.frag_abs_err_da, 0.5),
            "frag_err_ppm_quantiles": {p: q(sub.frag_err_ppm, p) for p in (0.05, 0.5, 0.95)},
            "low_ann_frac_lt_0p2": float((ann < 0.2).mean()),
        }
    S.update({
        "source": "flat", "dataset": df.dataset.iloc[0], "n_runs": int(df.run.nunique()),
        "by_inferred_analyser": per_an,
        "prec_err_ppm_quantiles": {p: q(df.prec_err_ppm, p) for p in (0.01, 0.05, 0.5, 0.95, 0.99)},
        "prec_fit_frac": float(df.prec_fit.mean()),
        "prec_iso_hist": {int(k): int(v) for k, v in df.prec_iso.value_counts().sort_index().items()},
        "charge_alt_fit_frac": float(df.charge_alt_fit.mean()),
        "deamidation_isotope_artifact": deamid_artifact(df),
        "added_cam_frac": float(df.added_cam.mean()),
        "unmapped_mod_psm_frac": float((df.n_unmapped_mods > 0).mean()),
        "missed_cleavage_hist": {int(k): int(v) for k, v in df.missed_cleavages.clip(upper=4).value_counts().sort_index().items()},
        "outlier_runs": R.loc[R.outlier_flags != "", ["run", "outlier_flags"]].values.tolist(),
        "decoy_separation": "n/a: flats carry targets only, no scores",
        "qvalue_pep": "n/a: flats are pre-filtered (MSFragger expect <= 0.01) and carry no score",
    })
    json.dump(S, open(os.path.join(a.out_dir, "dataset_qc.json"), "w"), indent=1, default=float)
    print(json.dumps({k: S[k] for k in ("dataset", "n_runs", "n_psms", "unique_peptidoforms", "isoleucine_written_as_L",
                                        "prec_fit_frac", "tryptic_cterm_frac")}, default=float))


# ----------------------------------------------------------------------------- psm_table
SCORE_CANDS = ["andromeda_score", "hyperscore", "search_engine_score", "score", "rank_score",
               "posterior_error_probability", "q_value", "global_qvalue"]


def auc(pos, neg):
    """P(target score > decoy score), Mann-Whitney."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if not len(pos) or not len(neg):
        return float("nan")
    allv = np.concatenate([pos, neg])
    ranks = allv.argsort().argsort().astype(float) + 1
    return float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def cmd_psm_table(a):
    import pandas as pd

    from ms2_common import parse_proforma  # noqa: E402

    os.makedirs(a.out_dir, exist_ok=True)
    frames = []
    for path in a.inputs:
        names = pq.ParquetFile(path).schema_arrow.names
        want = [c for c in ["peptidoform", "sequence", "precursor_charge", "charge", "posterior_error_probability",
                            "pep", "q_value", "global_qvalue", "is_decoy", "run_file_name", "rt", "scan",
                            "calculated_mz", "observed_mz"] + SCORE_CANDS if c in names]
        frames.append(pq.read_table(path, columns=sorted(set(want))).to_pandas())
    df = pd.concat(frames, ignore_index=True)
    if "charge" not in df and "precursor_charge" in df:
        df["charge"] = df.precursor_charge
    pep_col = next((c for c in ("posterior_error_probability", "pep") if c in df), None)
    qcol = next((c for c in ("q_value", "global_qvalue") if c in df), None)
    dec = df.is_decoy.fillna(False).astype(bool) if "is_decoy" in df else pd.Series(False, index=df.index)
    S = {"source": "psm_table", "dataset": a.dataset or "", "n_rows": int(len(df)), "decoy_frac": float(dec.mean())}
    for name, col in (("pep", pep_col), ("q", qcol)):
        if col:
            S[f"{name}_quantiles_targets"] = {p: q(df.loc[~dec, col], p) for p in (0.5, 0.9, 0.99, 1.0)}
    sep = {}
    if dec.any():
        for c in SCORE_CANDS:
            if c in df and df[c].notna().any():
                sgn = -1 if c in ("posterior_error_probability", "q_value", "global_qvalue") else 1
                sep[c] = auc(sgn * df.loc[~dec, c], sgn * df.loc[dec, c])
    S["target_decoy_auc"] = sep or "n/a: no decoys in table (pre-filtered)"
    if a.max_pep is not None and pep_col:
        tgt = df[~dec & (df[pep_col] <= a.max_pep)].copy()
    else:
        tgt = df[~dec].copy()
    parsed = tgt.peptidoform.map(parse_proforma)
    tgt["sequence"] = [p[0] if p else re.sub(r"\[[^\]]*\]|[^A-Z]", "", s) for p, s in zip(parsed, tgt.peptidoform)]
    tgt["mods"] = ["|".join(sorted({f"UNIMOD:{u}" for _, u in p[1]})) if p else "unparsed" for p in parsed]
    S.update(composition_summary(tgt))
    S["filter"] = f"targets, {pep_col} <= {a.max_pep}" if a.max_pep is not None and pep_col else "targets"
    if "run_file_name" in tgt:
        runs = tgt.groupby("run_file_name").agg(n_psms=("peptidoform", "size"),
                                                 unique_peptidoforms=("peptidoform", "nunique"),
                                                 rt_max=("rt", "max") if "rt" in tgt else ("peptidoform", "size"),
                                                 charge_mean=("charge", "mean"))
        runs["z_n_psms"] = robust_z(runs.n_psms)
        runs["outlier_flags"] = np.where(np.abs(runs.z_n_psms) > a.outlier_z, "n_psms", "")
        runs.to_csv(os.path.join(a.out_dir, "run_qc.tsv"), sep="\t", float_format="%.5g")
        S["n_runs"] = int(len(runs))
        S["outlier_runs"] = runs.index[runs.outlier_flags != ""].tolist()
    pq.write_table(pa.Table.from_pandas(tgt[["peptidoform", "sequence", "charge", "mods"]].assign(dataset=a.dataset or ""),
                                        preserve_index=False), os.path.join(a.out_dir, "psm_qc.parquet"), compression="zstd")
    json.dump(S, open(os.path.join(a.out_dir, "dataset_qc.json"), "w"), indent=1, default=float)
    print(json.dumps({k: S[k] for k in ("n_rows", "decoy_frac", "unique_peptidoforms", "isoleucine_written_as_L")}, default=float))


# ----------------------------------------------------------------------------- shard
def cmd_shard(a):
    import pandas as pd

    os.makedirs(a.out_dir, exist_ok=True)
    df = pd.concat([pq.read_table(p, columns=["dataset", "run", "peptidoform", "sequence", "charge", "activation",
                                              "n_matched", "explained", "pep"]).to_pandas() for p in a.inputs])
    df["mods"] = df.peptidoform.str.findall(r"\[(UNIMOD:\d+)\]").map(lambda x: "|".join(sorted(set(x))))
    S = composition_summary(df)
    S["source"] = "shard"
    S["explained_by_activation"] = {k: {p: q(g.explained, p) for p in (0.25, 0.5, 0.75)} for k, g in df.groupby("activation")}
    json.dump(S, open(os.path.join(a.out_dir, "dataset_qc.json"), "w"), indent=1, default=float)


# ----------------------------------------------------------------------------- overlap
def cmd_peptidoforms(a):
    import pandas as pd

    df = pd.concat([pq.read_table(p, columns=["dataset", "peptidoform"]).to_pandas() for p in a.inputs])
    df["key"] = df.peptidoform.map(il_key)
    df = df[["dataset", "key"]].drop_duplicates()
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False), a.out, compression="zstd")
    print(df.groupby("dataset").size().to_string())


def cmd_overlap(a):
    import pandas as pd

    df = pd.concat([pq.read_table(p).to_pandas() for p in a.inputs]).drop_duplicates()
    sets = {d: set(g.key) for d, g in df.groupby("dataset")}
    counts = collections.Counter(k for s in sets.values() for k in s)
    rows = []
    for d, s in sets.items():
        shared = sum(1 for k in s if counts[k] > 1)
        rows.append({"dataset": d, "peptidoforms": len(s), "shared_with_other_datasets": shared,
                     "unique_to_dataset": len(s) - shared})
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"union {len(counts)}; in >=2 datasets {sum(1 for v in counts.values() if v > 1)}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("flat", "psm_table", "shard"):
        s = sub.add_parser(name)
        s.add_argument("--out-dir", required=True)
        s.add_argument("--dataset", default=None)
        s.add_argument("--outlier-z", type=float, default=3.5)
        s.add_argument("inputs", nargs="+")
        if name == "flat":
            s.add_argument("--workers", type=int, default=os.cpu_count())
            s.add_argument("--prec-tol-ppm", type=float, default=20.0)
            s.add_argument("--limit", type=int, default=0, help="rows per file (0 = all)")
        if name == "psm_table":
            s.add_argument("--max-pep", type=float, default=0.01)
    s = sub.add_parser("peptidoforms")
    s.add_argument("--out", required=True)
    s.add_argument("inputs", nargs="+")
    s = sub.add_parser("overlap")
    s.add_argument("inputs", nargs="+")
    a = p.parse_args()
    {"flat": cmd_flat, "psm_table": cmd_psm_table, "shard": cmd_shard,
     "peptidoforms": cmd_peptidoforms, "overlap": cmd_overlap}[a.cmd](a)


if __name__ == "__main__":
    main()
