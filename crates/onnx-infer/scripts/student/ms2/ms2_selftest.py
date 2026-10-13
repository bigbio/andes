#!/usr/bin/env python3
"""End-to-end check of the andes-ms2 pipeline on synthetic data (CPU, about a minute).

Writes a qpx-like ``psm.parquet`` and two mzML runs whose MS2 spectra follow a
simple deterministic fragmentation rule plus noise peaks, then runs the real
command-line tools: ``ms2_manifest.py`` -> ``ms2_select.py`` ->
``ms2_extract_raw.py --mzml-dir`` -> ``ms2_train.py`` -> ``ms2_export.py``.
Checks that every selected PSM is annotated, that validation spectral angle
rises well above its starting value, and leaves the weights and the parity
reference in --work for ``ms2_native_check``.

Usage::

  python ms2_selftest.py --work DIR [--steps 600]
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import subprocess
import sys
import zlib

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ms2_common import fragment_mz, precursor_mz, to_proforma, valid_mask  # noqa: E402

AAS = "ACDEFGHIKLMNPQRSTVWY"
WEIGHT = {a: 0.2 + 0.8 * ((i * 7) % 20) / 19 for i, a in enumerate(AAS)}


def synth_peptide(rng):
    n = int(rng.integers(7, 21))
    seq = "".join(rng.choice(list(AAS.replace("K", "").replace("R", "")), n - 1)) + rng.choice(["K", "R"])
    mods = [(i + 1, 4) for i, c in enumerate(seq) if c == "C"]
    mods += [(i + 1, 35) for i, c in enumerate(seq) if c == "M" and rng.random() < 0.5]
    st = [i + 1 for i, c in enumerate(seq) if c in "ST"]
    if st and rng.random() < 0.3:
        mods.append((int(rng.choice(st)), 21))
    return seq, sorted(mods)


def synth_intensities(seq, mods, z):
    """The rule the model has to learn: y ions from the residue after the cleavage,
    b ions from the one before, P boosts the y ion N-terminal to it; 2+ fragments
    only at charge >= 3; neutral losses at a third of their parent."""
    n = len(seq)
    v = np.zeros((n - 1, 12))
    for i in range(n - 1):
        a, b = seq[i], seq[i + 1]
        y = WEIGHT[b] * (3.0 if b == "P" else 1.0) * (0.5 + i / n)
        bb = 0.6 * WEIGHT[a] * (1.0 - i / n)
        v[i, 0], v[i, 2] = bb, y
        if z >= 3:
            v[i, 1], v[i, 3] = 0.5 * bb * (i / n), 0.5 * y * (1 - i / n)
        v[i, 4:8] = v[i, 0:4] / 3
        v[i, 8:12] = v[i, 0:4] / 3
    return np.where(valid_mask(seq, mods, z), v, 0)


def b64(arr):
    return base64.b64encode(zlib.compress(np.asarray(arr, "<f8").tobytes())).decode()


def write_mzml(path, spectra):
    out = ['<?xml version="1.0" encoding="utf-8"?>',
           '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">',
           '<referenceableParamGroupList count="1"><referenceableParamGroup id="commonInstrumentParams">',
           '<cvParam cvRef="MS" accession="MS:1001911" name="Q Exactive" value=""/>',
           '<cvParam cvRef="MS" accession="MS:1000529" name="instrument serial number" value="X1"/>',
           '</referenceableParamGroup></referenceableParamGroupList>',
           f'<run id="r"><spectrumList count="{len(spectra)}">']
    for k, (scan, filt, pmz, z, mz, inten) in enumerate(spectra):
        out.append(
            f'<spectrum index="{k}" id="controllerType=0 controllerNumber=1 scan={scan}" defaultArrayLength="{len(mz)}">'
            '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>'
            f'<scanList count="1"><scan><cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="{filt}"/></scan></scanList>'
            '<precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>'
            f'<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{pmz:.6f}"/>'
            f'<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{z}"/>'
            '</selectedIon></selectedIonList><activation>'
            '<cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation" value=""/>'
            '<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="27.0"/>'
            '</activation></precursor></precursorList>'
            '<binaryDataArrayList count="2">'
            '<binaryDataArray><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>'
            '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>'
            f'<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/><binary>{b64(mz)}</binary></binaryDataArray>'
            '<binaryDataArray><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>'
            '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>'
            f'<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/><binary>{b64(inten)}</binary></binaryDataArray>'
            '</binaryDataArrayList></spectrum>')
    out.append("</spectrumList></run></mzML>")
    with open(path, "w") as fh:
        fh.write("\n".join(out))


def run(cmd):
    print("+", " ".join(cmd), flush=True)
    r = subprocess.run([sys.executable] + cmd, capture_output=True, text=True)
    sys.stdout.write(r.stdout[-3000:])
    if r.returncode != 0:
        sys.stderr.write(r.stderr[-5000:])
        sys.exit(f"failed: {cmd[0]}")
    return r.stdout


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", required=True)
    p.add_argument("--n-psms", type=int, default=4000)
    p.add_argument("--steps", type=int, default=600)
    p.add_argument("--threads", type=int, default=4)
    a = p.parse_args()
    rng = np.random.default_rng(1)
    w = a.work
    mz_dir = os.path.join(w, "mzml")
    os.makedirs(mz_dir, exist_ok=True)
    psm = {"peptidoform": [], "is_decoy": [], "run_file_name": [], "precursor_charge": [], "scan": [],
           "posterior_error_probability": []}
    runs = {"run_a": [], "run_b": []}
    for k in range(a.n_psms):
        seq, mods = synth_peptide(rng)
        z = int(rng.choice([2, 3], p=[0.6, 0.4]))
        run_name = "run_a" if k % 2 == 0 else "run_b"
        scan = 1000 + k
        theo = fragment_mz(seq, mods)
        inten = synth_intensities(seq, mods, z) * rng.lognormal(0, 0.15, theo.shape)
        sel = inten > 0
        noise_mz = rng.uniform(150, 1500, 40)
        mz = np.r_[theo[sel], noise_mz]
        it = np.r_[inten[sel] * 1e6, rng.uniform(0, 0.15, 40) * 1e6]
        order = np.argsort(mz)
        runs[run_name].append((scan, f"FTMS + p NSI d Full ms2 {precursor_mz(seq, mods, z):.2f}@hcd27.00 [100.00-2000.00]",
                               precursor_mz(seq, mods, z), z, mz[order], it[order]))
        psm["peptidoform"].append(to_proforma(seq, mods))
        psm["is_decoy"].append(False)
        psm["run_file_name"].append(run_name + ".raw")
        psm["precursor_charge"].append(z)
        psm["scan"].append(scan)
        psm["posterior_error_probability"].append(float(rng.uniform(0, 0.005)))
    # a decoy and a low-confidence row that selection must drop
    for dec, pep in ((True, 0.001), (False, 0.5)):
        for k_, v in (("peptidoform", "PEPTIDEK"), ("is_decoy", dec), ("run_file_name", "run_a.raw"),
                      ("precursor_charge", 2), ("scan", 1), ("posterior_error_probability", pep)):
            psm[k_].append(v)
    for name, spectra in runs.items():
        write_mzml(os.path.join(mz_dir, name + ".mzML"), spectra)
    psm_path = os.path.join(w, "PXD999999_psm.parquet")
    pq.write_table(pa.table(psm), psm_path)

    run([os.path.join(HERE, "ms2_manifest.py"), "--out", os.path.join(w, "datasets.tsv"),
         "--raw-template", os.path.join(w, "raw", "{dataset}"), f"PXD999999={psm_path}"])
    run([os.path.join(HERE, "ms2_select.py"), "--manifest", os.path.join(w, "datasets.tsv"),
         "--out-dir", os.path.join(w, "sel"), "--psms-per-task", "1500"])
    tasks = [l.split("\t")[0] for l in open(os.path.join(w, "sel", "tasks.tsv")).read().splitlines()[1:]]
    for t in tasks:
        run([os.path.join(HERE, "ms2_extract_raw.py"), "--tasks", os.path.join(w, "sel", "tasks.tsv"),
             "--task-index", t, "--out-dir", os.path.join(w, "shards"), "--work-dir", os.path.join(w, "scratch"),
             "--mzml-dir", mz_dir, "--overwrite"])
    kept = sum(json.load(open(os.path.join(w, "shards", "PXD999999", f"part-{int(t):05d}.json")))["total"]["kept"]
               for t in tasks)
    assert kept == a.n_psms, f"annotated {kept} of {a.n_psms} PSMs"
    out = run([os.path.join(HERE, "ms2_train.py"), "--shards", os.path.join(w, "shards"), "--out", os.path.join(w, "model"),
               "--holdout-datasets", "--steps", str(a.steps), "--batch", "64", "--eval-every", str(max(a.steps // 6, 1)),
               "--ch", "32", "--emb", "16", "--dilations", "1", "2", "4", "--n-heads", "2", "--lr", "3e-3",
               "--threads", str(a.threads)])
    hist = json.load(open(os.path.join(w, "model", "train.json")))["history"]
    first, best = hist[0]["val_median_sa_base"], max(h["val_median_sa_base"] for h in hist)
    print(f"val median SA: untrained {first:.3f}, best {best:.3f}")
    assert best > 0.8 and best > first + 0.2, "the model did not learn the synthetic rule"
    run([os.path.join(HERE, "ms2_export.py"), "--model-dir", os.path.join(w, "model"), "--shards", os.path.join(w, "shards"),
         "--n", "300", "--out-weights", os.path.join(w, "andes_ms2.ams2"),
         "--out-reference", os.path.join(w, "ms2_reference.tsv")])
    print("selftest ok;", out.strip().splitlines()[-1])


if __name__ == "__main__":
    main()
