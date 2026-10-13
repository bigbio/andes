#!/usr/bin/env python3
"""Add AlphaPeptDeep (generic MS2 model) predictions to andes-ms2 shards, for distillation.

For every shard row the teacher predicts b/y 1+/2+ and the modloss channels at
the row's own charge, NCE (--default-nce when unknown) and instrument, and the
predictions are written as a ``teacher`` column in the shard layout
(``ms2_common.ION_TYPES``, -1 where the teacher has no value). Unique
(peptidoform, charge, NCE, instrument) requests are predicted once.

PeptDeep's ``modloss`` channels carry the loss of the modification on the
fragment: H3PO4 for pS/pT, CH4SO for oxidised M. They fill the H3PO4 channels
when the peptide has a pS/pT, otherwise the CH4SO channels when it has an
oxidised M; the other loss channels stay -1.

Requires ``peptdeep`` (and its ``alphabase``) in the environment, with the
pretrained models installed (``peptdeep install-models``). Peptidoforms
carrying a modification without an alphabase name below get no teacher (-1).

Usage::

  python ms2_teacher.py --out-dir SHARDS_T SHARDS/...parquet [--device cuda]
"""
from __future__ import annotations

import argparse
import glob
import math
import os
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import INSTRUMENTS, N_ION, instrument_class, parse_proforma, valid_mask  # noqa: E402

# Unimod id -> alphabase name template ({aa} = residue); N-terminal variants separately.
ALPHABASE = {
    1: "Acetyl@Protein_N-term", 4: "Carbamidomethyl@C", 7: "Deamidated@{aa}", 21: "Phospho@{aa}",
    35: "Oxidation@{aa}", 737: "TMT6plex@{aa}", 2016: "TMTpro@{aa}", 121: "GlyGly@K", 34: "Methyl@{aa}",
    36: "Dimethyl@{aa}", 37: "Trimethyl@K", 27: "Glu->pyro-Glu@E^Any_N-term", 28: "Gln->pyro-Glu@Q^Any_N-term",
    214: "iTRAQ4plex@{aa}", 730: "iTRAQ8plex@{aa}", 26: "Pyro-carbamidomethyl@C^Any_N-term",
}
NTERM = {737: "TMT6plex@Any_N-term", 2016: "TMTpro@Any_N-term", 214: "iTRAQ4plex@Any_N-term",
         730: "iTRAQ8plex@Any_N-term", 1: "Acetyl@Protein_N-term", 34: "Methyl@Any_N-term",
         36: "Dimethyl@Any_N-term"}
# andes-ms2 instrument family -> PeptDeep instrument name
PDEEP_INSTRUMENT = {"unknown": "Lumos", "QExactive": "QE", "Exploris": "QE", "Tribrid": "Lumos",
                    "VelosElite": "Lumos", "Astral": "Lumos", "LTQ": "Lumos"}
FRAG = ["b_z1", "b_z2", "y_z1", "y_z2", "b_modloss_z1", "b_modloss_z2", "y_modloss_z1", "y_modloss_z2"]


def alphabase(seq, mods):
    names, sites = [], []
    for pos, uid in mods:
        if pos == 0:
            name = NTERM.get(uid) or ALPHABASE.get(uid)
            site = 0
        elif pos == len(seq) + 1:
            return None
        else:
            name = ALPHABASE.get(uid)
            site = pos
        if name is None:
            return None
        if "^Any_N-term" in name:
            site = 1  # residue-specific N-terminal mods sit on residue 1 in alphabase
        names.append(name.format(aa=seq[max(pos, 1) - 1]))
        sites.append(str(site))
    return ";".join(names), ";".join(sites)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--default-nce", type=float, default=30.0)
    p.add_argument("--device", default="cpu")
    p.add_argument("--model-type", default="generic")
    p.add_argument("shards", nargs="+")
    a = p.parse_args()
    import pandas as pd
    from peptdeep.pretrained_models import ModelManager

    mm = ModelManager(mask_modloss=False, device=a.device)
    mm.load_installed_models(a.model_type)
    files = []
    for s in a.shards:
        files += sorted(glob.glob(os.path.join(s, "**", "*.parquet"), recursive=True)) if os.path.isdir(s) else [s]
    os.makedirs(a.out_dir, exist_ok=True)
    for f in files:
        t = pq.read_table(f)
        df = t.select(["peptidoform", "charge", "nce", "instrument"]).to_pandas()
        df["nce_t"] = df.nce.fillna(a.default_nce).round(1)
        df["instr_t"] = [PDEEP_INSTRUMENT[INSTRUMENTS[instrument_class(x)]] for x in df.instrument]
        keys = df[["peptidoform", "charge", "nce_t", "instr_t"]].drop_duplicates().reset_index(drop=True)
        rows = []
        for k, r in keys.iterrows():
            pm = parse_proforma(r.peptidoform)
            ab = alphabase(*pm) if pm else None
            if ab is None:
                continue
            rows.append({"key": k, "sequence": pm[0], "mods": ab[0], "mod_sites": ab[1], "charge": int(r.charge),
                         "nce": float(r.nce_t), "instrument": r.instr_t})
        pre = pd.DataFrame(rows)
        teach = {}
        if len(pre):
            res = mm.predict_all(pre, predict_items=["ms2"], frag_types=FRAG, multiprocessing=False)
            pdf, inten = res["precursor_df"], res["fragment_intensity_df"]
            stop = "frag_stop_idx" if "frag_stop_idx" in pdf.columns else "frag_end_idx"
            for key, lo, hi in zip(pdf.key, pdf.frag_start_idx, pdf[stop]):
                teach[key] = inten.iloc[lo:hi][FRAG].to_numpy(np.float32)
        out = []
        key_of = {tuple(x): i for i, x in enumerate(keys.itertuples(index=False, name=None))}
        for r in df.itertuples(index=False):
            k = key_of[(r.peptidoform, r.charge, r.nce_t, r.instr_t)]
            seq, mods = parse_proforma(r.peptidoform)
            n = len(seq)
            v = np.full((n - 1, N_ION), -1.0, np.float32)
            if k in teach:
                blk = teach[k]
                valid = valid_mask(seq, mods, int(r.charge))
                v[:, :4] = np.where(valid[:, :4], blk[:, :4], -1)
                has_ph = valid[:, 4:8].any()
                lo = 4 if has_ph else 8
                v[:, lo: lo + 4] = np.where(valid[:, lo: lo + 4], blk[:, 4:8], -1)
                top = v.max()
                if top > 0:
                    v = np.where(v >= 0, v / top, -1).astype(np.float32)
            out.append(v.ravel())
        t = t.append_column("teacher", pa.array(out, pa.list_(pa.float32())))
        dst = os.path.join(a.out_dir, os.path.relpath(f, os.path.commonpath(files)) if len(files) > 1 else os.path.basename(f))
        os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
        pq.write_table(t, dst, compression="zstd")
        print(f"{f}: {len(df)} rows, {len(keys)} requests, {len(teach)} predicted -> {dst}", flush=True)


if __name__ == "__main__":
    main()
