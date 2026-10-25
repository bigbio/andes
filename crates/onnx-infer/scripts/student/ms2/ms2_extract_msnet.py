#!/usr/bin/env python3
"""Annotate MsNet-format PSM tables (peak lists attached) into andes-ms2 shards.

Two MsNet layouts are read, by column name:

- layout A: ``charge, calculated_mz, observed_mz, run_file_name, rt`` and ``cv_params`` as a list of
  ``{cv_name, cv_value}``;
- layout B: ``precursor_charge, exp_mass_to_charge, reference_file_name, retention_time`` and
  ``cv_params`` as a map (``Instrument``, ``Fragmentation Method``, ``Collision Energy``, ``Mass Analyzer``).

Peptidoforms are written with named modifications (``N[Deamidated]``, ``[Acetyl]-PEP``); the names are
mapped to the built-in Unimod table. A PSM with an unknown name, a precursor m/z that does not fit the
parsed peptidoform (isotope offsets -1..3 within ``--prec-tol-ppm``), an unknown activation, or fewer than
``--min-matched`` matched plain b/y ions is dropped and counted. The annotation arrays of MsNet are not
used: ions are annotated with ``ms2_common.annotate`` from the peak list, as for every other source.

Activation and NCE come from each PSM's ``cv_params`` (dissociation method and mass analyzer give
``HCD-FTMS``, ``CID-ITMS``...; the collision energy is the NCE for orbitrap HCD and ion trap, NaN when
missing or for time-of-flight instruments). Fragment tolerance is 20 ppm for FTMS and 0.5 Da for ITMS.

Sampling: PSMs with PEP <= ``--pep-max``, targets only, at most ``--per-key`` PSMs per (peptidoform,
charge, NCE bin of 2) and ``--cap`` PSMs per dataset. Large files are sampled by row group (random, seeded
by the dataset name), with a per-group quota, so the sample spreads over the runs of the file.

Usage::

  python ms2_extract_msnet.py --dataset PXD004732 --out shards/PXD004732 --cap 60000 FILE.parquet...
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import re
import sys

import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import ACTIVATIONS, annotate, instrument_class, precursor_mz, shard_table, to_proforma  # noqa: E402
from mod_composition import UNIMOD  # noqa: E402

NAME_TO_ID = {v[0].lower(): k for k, v in UNIMOD.items()}
NAME_TO_ID.update({
    "gln->pyro-glu": 28, "glu->pyro-glu": 27, "carbamidomethyl": 4, "tmt6plex": 737, "tmt10plex": 737,
    "tmt11plex": 737, "tmtpro": 2016, "tmt": 739, "label:13c(6)15n(2)": 259, "label:13c(6)15n(4)": 267,
    "itraq4plex": 214, "itraq8plex": 730, "phosphorylation": 21, "oxidation": 35, "acetylation": 1,
})
_TOKEN = re.compile(r"\[([^\]]+)\]")


def parse_named(s: str):
    """``[Acetyl]-PEPM[Oxidation]K-[Amidated]`` -> (sequence, [(pos, unimod id)]) or None for an unknown
    modification name, mass tag or non-standard residue. Position 0 = N-terminus, n + 1 = C-terminus."""
    mods = []
    i = 0
    n_term = []
    if s.startswith("["):
        j = s.find("]-")
        if j < 0:
            return None
        n_term = _TOKEN.findall(s[: j + 1])
        i = j + 2
    seq = []
    while i < len(s):
        c = s[i]
        if c == "[":
            m = _TOKEN.match(s, i)
            if m is None or not seq:
                return None
            uid = NAME_TO_ID.get(m.group(1).strip().lower())
            if uid is None:
                return None
            mods.append((len(seq), uid))
            i = m.end()
        elif c == "-":
            tail = _TOKEN.findall(s[i + 1:])
            if not tail:
                return None
            for t in tail:
                uid = NAME_TO_ID.get(t.strip().lower())
                if uid is None:
                    return None
                mods.append((len(seq) + 1, uid))
            break
        elif "A" <= c <= "Z":
            seq.append(c)
            i += 1
        else:
            return None
    if not seq:
        return None
    for t in n_term:
        uid = NAME_TO_ID.get(t.strip().lower())
        if uid is None:
            return None
        mods.append((0, uid))
    return "".join(seq), sorted(mods)


def cv_dict(cv) -> dict:
    """Lower-case name -> value for either cv_params layout."""
    if not cv:
        return {}
    if isinstance(cv, dict):
        return {str(k).lower(): v for k, v in cv.items()}
    return {str(kv.get("cv_name", "")).lower(): kv.get("cv_value") for kv in cv if isinstance(kv, dict)}


def conditions(cv):
    """-> (activation name, analyser, NCE or nan, instrument model) from a cv_params value."""
    d = cv_dict(cv)
    analyser = str(d.get("mass analyzer") or "").upper()
    diss = str(d.get("dissociation method") or d.get("fragmentation method") or "").upper()
    name = f"{diss}-{analyser}"
    act = name if name in ACTIVATIONS else "unknown"
    inst = str(d.get("instrument") or "")
    nce = float("nan")
    raw = d.get("normalized collision energy", d.get("collision energy"))
    try:
        nce = float(raw)
    except (TypeError, ValueError):
        pass
    if "tims" in inst.lower() or analyser in ("TOF", "QTOF"):
        nce = float("nan")
    return act, analyser, nce, inst


def fits(seq, mods, z, mz, tol_ppm):
    theo = precursor_mz(seq, mods, z)
    return any(abs(mz - (theo + k * 1.0033548378 / z)) <= tol_ppm * 1e-6 * theo for k in range(-1, 4))


def first_present(names, have):
    for n in names:
        if n in have:
            return n
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cap", type=int, default=60000, help="PSMs kept per dataset")
    p.add_argument("--per-key", type=int, default=4, help="PSMs kept per (peptidoform, charge, NCE bin)")
    p.add_argument("--pep-max", type=float, default=0.01)
    p.add_argument("--prec-tol-ppm", type=float, default=20.0)
    p.add_argument("--min-matched", type=int, default=4)
    p.add_argument("--max-groups", type=int, default=12, help="row groups sampled per file")
    p.add_argument("--batch", type=int, default=20000)
    p.add_argument("files", nargs="+")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(int(hashlib.md5(a.dataset.encode()).hexdigest()[:8], 16))
    st = {"rows": 0, "decoy_or_pep": 0, "bad_charge": 0, "unmapped": 0, "precursor_mismatch": 0, "activation_unknown": 0,
          "too_few_ions": 0, "key_full": 0, "kept": 0, "files": 0}
    seen: dict = {}
    rows = []
    kept_total = 0
    for fi, path in enumerate(a.files):
        pf = pq.ParquetFile(path)
        names = set(pf.schema_arrow.names)
        c_z = first_present(["charge", "precursor_charge"], names)
        c_mz = first_present(["observed_mz", "exp_mass_to_charge"], names)
        c_run = first_present(["run_file_name", "reference_file_name"], names)
        need = [c_z, c_mz, c_run]
        if None in need or not {"peptidoform", "mz_array", "intensity_array", "cv_params", "posterior_error_probability",
                                "is_decoy", "scan"} <= names:
            print(f"{path}: unsupported layout {sorted(names)}", flush=True)
            continue
        st["files"] += 1
        cols = ["peptidoform", c_z, c_mz, c_run, "scan", "cv_params", "mz_array", "intensity_array",
                "posterior_error_probability", "is_decoy"]
        groups = list(range(pf.metadata.num_row_groups))
        if len(groups) > a.max_groups:
            groups = sorted(rng.sample(groups, a.max_groups))
        quota = max(1, (a.cap - kept_total) // max(1, len(groups))) if a.cap > kept_total else 0
        for g in groups:
            if kept_total >= a.cap:
                break
            got = 0
            for batch in pf.iter_batches(batch_size=a.batch, row_groups=[g], columns=cols):
                d = batch.to_pydict()
                for i in range(batch.num_rows):
                    st["rows"] += 1
                    pep = d["posterior_error_probability"][i]
                    if d["is_decoy"][i] or pep is None or not pep <= a.pep_max:
                        st["decoy_or_pep"] += 1
                        continue
                    z = d[c_z][i]
                    if z is None or not 1 <= int(z) <= 8:
                        st["bad_charge"] += 1
                        continue
                    z = int(z)
                    parsed = parse_named(d["peptidoform"][i])
                    if parsed is None or len(parsed[0]) < 5:
                        st["unmapped"] += 1
                        continue
                    seq, mods = parsed
                    if not fits(seq, mods, z, float(d[c_mz][i]), a.prec_tol_ppm):
                        st["precursor_mismatch"] += 1
                        continue
                    act, analyser, nce, inst = conditions(d["cv_params"][i])
                    if act == "unknown" and analyser not in ("TOF", "QTOF"):
                        st["activation_unknown"] += 1
                        continue
                    key = (to_proforma(seq, mods), z, -1 if math.isnan(nce) else int(nce // 2), act)
                    if seen.get(key, 0) >= a.per_key:
                        st["key_full"] += 1
                        continue
                    mz = np.asarray(d["mz_array"][i], np.float64)
                    inten = np.asarray(d["intensity_array"][i], np.float64)
                    if len(mz) < 10:
                        st["too_few_ions"] += 1
                        continue
                    order = np.argsort(mz, kind="stable")
                    ppm = analyser != "ITMS"
                    vals, n_base, explained = annotate(seq, mods, z, mz[order], inten[order], 20.0 if ppm else 0.5, ppm)
                    if vals is None or n_base < a.min_matched:
                        st["too_few_ions"] += 1
                        continue
                    seen[key] = seen.get(key, 0) + 1
                    scan = d["scan"][i]
                    scan = scan[0] if isinstance(scan, list) and scan else scan
                    rows.append({
                        "dataset": a.dataset, "run": str(d[c_run][i]), "scan": int(scan) if scan is not None else -1,
                        "peptidoform": key[0], "sequence": seq, "charge": z, "nce": nce, "activation": act,
                        "instrument": inst, "pep": float(pep), "n_matched": n_base, "explained": explained,
                        "intensities": vals,
                    })
                    st["kept"] += 1
                    kept_total += 1
                    got += 1
                    if got >= quota or kept_total >= a.cap:
                        break
                if got >= quota or kept_total >= a.cap:
                    break
            if len(rows) >= 50000:
                pq.write_table(shard_table(rows), os.path.join(a.out, f"part-msnet-{fi:03d}-{g:04d}.parquet"), compression="zstd")
                rows = []
        print(f"{path}: {json.dumps(st)}", flush=True)
    if rows:
        pq.write_table(shard_table(rows), os.path.join(a.out, "part-msnet-last.parquet"), compression="zstd")
    json.dump(st, open(os.path.join(a.out, "msnet.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
