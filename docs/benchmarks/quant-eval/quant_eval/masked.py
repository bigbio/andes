"""Masked wrong-peak truth of the transfers, as a curve over the quant_q limit.

Inputs (all produced by andes, nothing is searched here):
  full   quant dir of the run on all identifications: library.parquet, alignment.parquet, <stem>.extract.parquet
  masked quant dir of the run in which 5 % of the identified (precursor, run) pairs were masked (their PSMs removed from the library), same files
  pairs  mask_pairs.parquet (columns precursor, run): the masked pairs
  fin    finalize parquet of the masked quant dir
A masked pair that finalize transfers again is a labelled fit when the full run also has a fit of it. The transfer picked the wrong peak
when its apex differs from the full run's identified apex by more than the identified fit's peak width (|t0 masked - t0 full| > w).
Level A = the transfer carries MS2 evidence (>= 4 matched fragments, angle >= 0.6), level B = MS1 only.

`wrong_curve` returns, for each quant_q limit tau of the grid, the wrong-peak counts of the labelled transfers with quant_q <= tau, pooled and
per level, and the per-run count above 1 % for runs with >= 200 pairs. The rate at the tau of interest and the tau of a target rate are read with
`rate_at` and `tau_for_rate`."""
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

GRID = [0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1]
EX_COLS = ["precursor", "run", "kind", "role", "t0", "w", "quantity", "ms2_n_matched", "ms2_angle"]


def _stems(qdir):
    a = pq.read_table(qdir + "/alignment.parquet", columns=["run", "stem"]).to_pandas().sort_values("run")
    return list(a.stem.values)


def _lib(qdir):
    return pq.read_table(qdir + "/library.parquet", columns=["key"]).to_pandas().key.values


def _extract(qdir, stems, role):
    parts = []
    for s in stems:
        have = pq.read_schema(f"{qdir}/{s}.extract.parquet").names
        d = pq.read_table(f"{qdir}/{s}.extract.parquet", columns=[c for c in EX_COLS if c in have]).to_pandas()
        parts.append(d[(d.kind == 0) & (d.role == role)])
    d = pd.concat(parts, ignore_index=True)
    d["code"] = d.precursor.values.astype(np.int64) * 65536 + d.run.values.astype(np.int64)
    return d.drop_duplicates("code").set_index("code")


def labelled_transfers(full_dir, masked_dir, pairs_path, fin_path):
    stems = _stems(full_dir)
    r2i = {s: i for i, s in enumerate(stems)}
    kf, km = _lib(full_dir), _lib(masked_dir)
    k2f = pd.Series(np.arange(len(kf)), index=kf)
    k2f = k2f[~k2f.index.duplicated()]
    k2m = pd.Series(np.arange(len(km)), index=km)
    k2m = k2m[~k2m.index.duplicated()]
    # masked finalize transfers
    F = pq.read_table(fin_path, columns=["precursor_key", "run", "evidence_type", "quant_q"]).to_pandas()
    F["pm"] = k2m.reindex(F.precursor_key.values).values
    F = F[F.pm.notna() & F.evidence_type.str.startswith("transfer")].copy()
    F["ri"] = F.run.map(r2i)
    F = F[F.ri.notna()]
    F["code"] = F.pm.astype(np.int64).values * 65536 + F.ri.astype(np.int64).values
    F = F.drop_duplicates("code").set_index("code")
    exm = _extract(masked_dir, stems, 1)
    exf = _extract(full_dir, stems, 0)
    mp = pd.read_parquet(pairs_path)
    mp["pf"] = k2f.reindex(mp.precursor.values).values
    mp["pm"] = k2m.reindex(mp.precursor.values).values
    mp["ri"] = mp.run.map(r2i)
    mp = mp[mp.pf.notna() & mp.pm.notna() & mp.ri.notna()].copy()
    mp["codef"] = mp.pf.astype(np.int64).values * 65536 + mp.ri.astype(np.int64).values
    mp["code"] = mp.pm.astype(np.int64).values * 65536 + mp.ri.astype(np.int64).values
    j = exf.index.get_indexer(mp.codef.values)
    for c, n in (("t0", "f_t0"), ("w", "f_w")):
        mp[n] = np.where(j >= 0, exf[c].values[np.maximum(j, 0)], np.nan)
    j = exm.index.get_indexer(mp.code.values)
    mp["m_t0"] = np.where(j >= 0, exm.t0.values[np.maximum(j, 0)], np.nan)
    nm = np.where(j >= 0, exm.ms2_n_matched.values[np.maximum(j, 0)], 0)
    an = np.where(j >= 0, exm.ms2_angle.values[np.maximum(j, 0)], 0)
    mp["levelA"] = (nm >= 4) & (an >= 0.6)
    j = F.index.get_indexer(mp.code.values)
    mp["quant_q"] = np.where(j >= 0, F.quant_q.values[np.maximum(j, 0)], np.nan)
    lab = np.isfinite(mp.quant_q) & np.isfinite(mp.f_t0) & np.isfinite(mp.m_t0) & np.isfinite(mp.f_w)
    mp = mp[lab].copy()
    mp["wrong"] = np.abs(mp.m_t0 - mp.f_t0) > mp.f_w
    mp["run_stem"] = [stems[int(i)] for i in mp.ri]
    return mp[["run_stem", "quant_q", "levelA", "wrong"]]


def wrong_curve(lt, grid=GRID):
    out = []
    for t in grid:
        k = lt[lt.quant_q <= t]
        o = {"tau": t, "n": int(len(k)), "wrong": int(k.wrong.sum())}
        for lab, m in (("A", k.levelA), ("B", ~k.levelA)):
            o[lab] = {"n": int(m.sum()), "wrong": int((k.wrong & m).sum())}
        g = k.groupby("run_stem").wrong.agg(["size", "sum"])
        g = g[g["size"] >= 200]
        o["runs_ge200"] = int(len(g))
        o["runs_over_1pct"] = int(((g["sum"] / g["size"]) > 0.01).sum())
        out.append(o)
    return out


def rate_at(curve, tau, level=None):
    p = min(curve, key=lambda x: abs(x["tau"] - tau))
    o = p if level is None else p[level]
    return o["wrong"] / o["n"] if o["n"] else float("nan")


def tau_for_rate(curve, rate, min_n=200):
    """Loosest grid tau whose pooled wrong-peak rate is <= rate (and that keeps >= min_n labelled transfers); None when none does."""
    best = None
    for p in curve:
        if p["n"] >= min_n and p["wrong"] / p["n"] <= rate and (best is None or p["tau"] > best):
            best = p["tau"]
    return best
