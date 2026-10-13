#!/usr/bin/env python3
"""Annotate andes training flats (peak lists already attached) into andes-ms2 shards.

Flat schema: ``seq, charge, prec_mz, res_mod_pos (1-based), res_mod_delta,
nterm_delta, cterm_delta, mz[], intensity[]``. Mass deltas are mapped to the
built-in Unimod table (``mod_composition.match_mass``; an N-terminal delta to an
N-term-capable entry). Flats may omit fixed carbamidomethylation: when the
precursor m/z only fits with +57.0215 on every C, it is added.

The flats carry no run, scan, NCE or activation, so those come from arguments
(--activation, --instrument, --nce; NCE NaN = unknown). Rows whose precursor
does not fit the parsed peptidoform (isotope offsets -1..3 within
--prec-tol-ppm), carry an unmapped delta, or match fewer than --min-matched plain
b/y ions are dropped and counted.

Usage::

  python ms2_extract_flats.py --dataset PXD009875 --out shards/PXD009875 \\
      --activation HCD-FTMS --instrument "Q Exactive" --tol 20 --unit ppm FLAT.parquet...
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import ACTIVATIONS, annotate, precursor_mz, shard_table, to_proforma  # noqa: E402
from mod_composition import match_mass  # noqa: E402

CAM = 4


def parse_flat(seq, pos, delta, nterm, cterm):
    """-> (seq, [(pos, unimod)]) or None for an unmapped delta / residue."""
    mods = []
    for p, d in zip(pos or [], delta or []):
        if not 1 <= p <= len(seq):
            return None
        uid = match_mass(float(d), seq[p - 1], False)
        if uid is None:
            # an N-terminal mod folded onto residue 1
            uid = match_mass(float(d), None, True) if p == 1 else None
            if uid is None:
                return None
            mods.append((0, uid))
            continue
        mods.append((p, uid))
    if nterm and abs(nterm) > 1e-6:
        uid = match_mass(float(nterm), None, True)
        if uid is None:
            return None
        mods.append((0, uid))
    if cterm and abs(cterm) > 1e-6:
        uid = match_mass(float(cterm), None, False)
        if uid is None:
            return None
        mods.append((len(seq) + 1, uid))
    return seq, sorted(mods)


def fits(seq, mods, z, mz, tol_ppm):
    theo = precursor_mz(seq, mods, z)
    return any(abs(mz - (theo + k * 1.0033548378 / z)) <= tol_ppm * 1e-6 * theo for k in range(-1, 4))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", required=True)
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--activation", default="HCD-FTMS", choices=ACTIVATIONS)
    p.add_argument("--instrument", default="")
    p.add_argument("--nce", type=float, default=float("nan"))
    p.add_argument("--tol", type=float, default=20.0)
    p.add_argument("--unit", choices=["ppm", "da"], default="ppm")
    p.add_argument("--prec-tol-ppm", type=float, default=20.0)
    p.add_argument("--min-matched", type=int, default=4)
    p.add_argument("--limit", type=int, default=0, help="stop after this many input rows (0 = all)")
    p.add_argument("flats", nargs="+")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    st = {"rows": 0, "unmapped": 0, "added_cam": 0, "precursor_mismatch": 0, "too_few_ions": 0, "kept": 0}
    cols = ["seq", "charge", "prec_mz", "res_mod_pos", "res_mod_delta", "nterm_delta", "cterm_delta", "mz", "intensity"]
    for fi, path in enumerate(a.flats):
        rows = []
        row_in_file = -1
        for batch in pq.ParquetFile(path).iter_batches(batch_size=4096, columns=cols):
            d = batch.to_pydict()
            for i in range(batch.num_rows):
                st["rows"] += 1
                row_in_file += 1
                seq, z = d["seq"][i], int(d["charge"][i])
                parsed = parse_flat(seq, d["res_mod_pos"][i], d["res_mod_delta"][i], d["nterm_delta"][i], d["cterm_delta"][i])
                if parsed is None or z < 1:
                    st["unmapped"] += 1
                    continue
                seq, mods = parsed
                pmz = float(d["prec_mz"][i])
                if not fits(seq, mods, z, pmz, a.prec_tol_ppm):
                    cam_free = [k + 1 for k, c in enumerate(seq) if c == "C" and (k + 1, CAM) not in mods]
                    with_cam = sorted(mods + [(k, CAM) for k in cam_free])
                    if cam_free and fits(seq, with_cam, z, pmz, a.prec_tol_ppm):
                        mods = with_cam
                        st["added_cam"] += 1
                    else:
                        st["precursor_mismatch"] += 1
                        continue
                mz = np.asarray(d["mz"][i], np.float64)
                inten = np.asarray(d["intensity"][i], np.float64)
                order = np.argsort(mz, kind="stable")
                vals, n_base, explained = annotate(seq, mods, z, mz[order], inten[order], a.tol, a.unit == "ppm")
                if vals is None or n_base < a.min_matched:
                    st["too_few_ions"] += 1
                    continue
                st["kept"] += 1
                rows.append({
                    "dataset": a.dataset, "run": os.path.basename(path), "scan": row_in_file, "peptidoform": to_proforma(seq, mods),
                    "sequence": seq, "charge": z, "nce": a.nce, "activation": a.activation, "instrument": a.instrument,
                    "pep": float("nan"), "n_matched": n_base, "explained": explained, "intensities": vals,
                })
                if a.limit and st["rows"] >= a.limit:
                    break
            if a.limit and st["rows"] >= a.limit:
                break
        if rows:
            pq.write_table(shard_table(rows), os.path.join(a.out, f"part-flat-{fi:05d}.parquet"), compression="zstd")
        print(f"{path}: {json.dumps(st)}", flush=True)
        if a.limit and st["rows"] >= a.limit:
            break
    json.dump(st, open(os.path.join(a.out, "flats.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
