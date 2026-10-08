#!/usr/bin/env python3
"""Generate the synthetic quantification fixtures in test-fixtures/.

Two small mzML files built from six BSA tryptic peptides with pyOpenMS
(`pip install pyopenms`):

* lfq_synthetic.mzML.gz   — label-free: MS1 scans at 1 Hz over 240 s with each
  peptide eluting as a Gaussian (sigma 4 s) of known apex intensity, plus one
  HCD MS2 per peptide 2 s before its apex. Used by the `--lfq` end-to-end test,
  which checks the integrated areas against the designed intensities.
* tmt10_synthetic.mzML.gz — TMT10: the same peptides carrying TMT on K and the
  N-terminus; every MS2 holds the ten reporter ions in the designed channel
  pattern (channels 6..10 at twice channels 1..5), an SPS-MS3 scan referencing
  the MS2 carries a different pattern (1..10 × scale), and the MS1 of the
  second peptide carries an interfering peak inside its isolation window so its
  precursor purity is 2/3.

The exact numbers the tests assert are the constants below.

Usage: scripts/make_quant_fixtures.py [--out test-fixtures]
"""
import argparse
import gzip
import math
import os
import random
import shutil

import pyopenms as p

PEPTIDES = [
    # sequence, apex RT (s), apex intensity (mono), charge
    ("LVNELTEFAK", 40.0, 1.0e6, 2),
    ("YLYEIAR", 70.0, 5.0e5, 2),
    ("AEFVEVTK", 100.0, 2.0e6, 2),
    ("HLVDEPQNLIK", 130.0, 8.0e5, 2),
    ("LGEYGFQNALIVR", 160.0, 3.0e5, 2),
    ("DAFLGSFLYEYSR", 190.0, 1.5e6, 2),
]
SIGMA_S = 4.0
GRADIENT_S = 240
MS2_BEFORE_APEX_S = 2.0
N_ISOTOPES = 5
TMT_REPORTERS = [126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
                 129.131471, 129.137790, 130.134825, 130.141145, 131.138180]
# MS2 reporter pattern: channels 1..5 at 1x, 6..10 at 2x (times the peptide scale)
TMT_MS2_PATTERN = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
# MS3 reporter pattern: 1..10 (times the peptide scale / 2)
TMT_MS3_PATTERN = list(range(1, 11))
TMT_SCALE = 2.0e4
# The second peptide gets an interfering MS1 peak at 0.3 Th above its mono
# with half the mono intensity: purity = envelope / (envelope + interferer).
INTERFERED_PEPTIDE = 1
INTERFERER_OFFSET = 0.3
INTERFERER_FRACTION = 0.5


def envelope(seq: p.AASequence, n: int):
    gen = p.CoarseIsotopePatternGenerator(n)
    dist = seq.getFormula().getIsotopeDistribution(gen)
    return [c.getIntensity() for c in dist.getContainer()]


def fragments(seq: p.AASequence, charge: int):
    """Theoretical b/y ions (charge 1 and 2) with a deterministic intensity shape."""
    tsg = p.TheoreticalSpectrumGenerator()
    params = tsg.getParameters()
    params.setValue("add_b_ions", "true")
    params.setValue("add_y_ions", "true")
    params.setValue("add_metainfo", "true")
    tsg.setParameters(params)
    spec = p.MSSpectrum()
    tsg.getSpectrum(spec, seq, 1, min(2, charge))
    names = spec.getStringDataArrays()[0]
    out = []
    n = seq.size()
    for i in range(spec.size()):
        name = names[i].decode() if isinstance(names[i], bytes) else str(names[i])
        mz = spec[i].getMZ()
        # y ions strong and rising with length, b ions weaker.
        try:
            idx = int(name[1:].split("+")[0].split("-")[0])
        except ValueError:
            idx = 1
        if name.startswith("y"):
            inten = 2000.0 * (0.4 + idx / n)
        else:
            inten = 600.0 * (0.4 + idx / n)
        if "++" in name:
            inten *= 0.3
        out.append((mz, inten))
    return out


def make_run(tmt: bool, path: str):
    rnd = random.Random(7 if tmt else 11)
    exp = p.MSExperiment()
    instrument = exp.getInstrument()
    analyzer = p.MassAnalyzer()
    analyzer.setType(p.MassAnalyzer.AnalyzerType.ORBITRAP)
    instrument.setMassAnalyzers([analyzer])
    instrument.setName("Q Exactive")
    exp.setInstrument(instrument)

    peps = []
    for seq_str, apex_rt, apex_int, z in PEPTIDES:
        if tmt:
            s = ".(TMT6plex)" + seq_str.replace("K", "K(TMT6plex)")
        else:
            s = seq_str
        seq = p.AASequence.fromString(s)
        mono = seq.getMonoWeight(p.Residue.ResidueType.Full, z) / z
        peps.append((seq, mono, apex_rt, apex_int, z, envelope(seq, N_ISOTOPES)))

    scan = 1
    spectra = []
    ms2_of = {}
    t = 0.0
    while t <= GRADIENT_S:
        ms1 = p.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.setRT(t)
        ms1.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
        ms1.setType(p.SpectrumSettings.SpectrumType.CENTROID)
        peaks = []
        # background: 200 random low-intensity peaks
        for _ in range(200):
            peaks.append((rnd.uniform(350.0, 1500.0), rnd.uniform(50.0, 400.0)))
        for k, (seq, mono, apex_rt, apex_int, z, env) in enumerate(peps):
            g = math.exp(-0.5 * ((t - apex_rt) / SIGMA_S) ** 2)
            if g < 1e-3:
                continue
            for i, e in enumerate(env):
                inten = apex_int * e / env[0] * g
                if inten > 10.0:
                    peaks.append((mono + i * 1.0033548 / z, inten))
            if k == INTERFERED_PEPTIDE:
                peaks.append((mono + INTERFERER_OFFSET, apex_int * INTERFERER_FRACTION * g))
        peaks.sort()
        ms1.set_peaks(([mz for mz, _ in peaks], [it for _, it in peaks]))
        spectra.append(ms1)
        scan += 1
        # MS2 (and MS3) for peptides whose trigger time falls in this cycle
        for k, (seq, mono, apex_rt, apex_int, z, env) in enumerate(peps):
            trigger = apex_rt - MS2_BEFORE_APEX_S
            if not (t <= trigger < t + 1.0):
                continue
            ms2 = p.MSSpectrum()
            ms2.setMSLevel(2)
            ms2.setRT(t + 0.4)
            ms2.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
            ms2.setType(p.SpectrumSettings.SpectrumType.CENTROID)
            pre = p.Precursor()
            pre.setMZ(mono)
            pre.setCharge(z)
            pre.setIntensity(apex_int)
            pre.setIsolationWindowLowerOffset(0.7)
            pre.setIsolationWindowUpperOffset(0.7)
            pre.setActivationMethods({p.Precursor.ActivationMethod.HCD})
            ms2.setPrecursors([pre])
            frag = fragments(seq, z)
            if tmt:
                scale = TMT_SCALE * (k + 1)
                frag += [(mz, scale * f) for mz, f in zip(TMT_REPORTERS, TMT_MS2_PATTERN)]
            frag.sort()
            ms2.set_peaks(([mz for mz, _ in frag], [it for _, it in frag]))
            spectra.append(ms2)
            ms2_of[k] = ms2.getNativeID()
            scan += 1
            if tmt:
                ms3 = p.MSSpectrum()
                ms3.setMSLevel(3)
                ms3.setRT(t + 0.6)
                ms3.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
                ms3.setType(p.SpectrumSettings.SpectrumType.CENTROID)
                notch = p.Precursor()
                notch.setMZ(frag[-3][0])
                notch.setCharge(1)
                notch.setMetaValue("spectrum_ref", ms2.getNativeID())
                notch.setActivationMethods({p.Precursor.ActivationMethod.HCD})
                ms3.setPrecursors([notch])
                scale = TMT_SCALE * (k + 1) / 2.0
                rep = [(mz, scale * f) for mz, f in zip(TMT_REPORTERS, TMT_MS3_PATTERN)]
                ms3.set_peaks(([mz for mz, _ in rep], [it for _, it in rep]))
                spectra.append(ms3)
                scan += 1
        t += 1.0

    exp.setSpectra(spectra)
    tmp = path[:-3] if path.endswith(".gz") else path
    p.MzMLFile().store(tmp, exp)
    if path.endswith(".gz"):
        with open(tmp, "rb") as src, gzip.open(path, "wb", compresslevel=9) as dst:
            shutil.copyfileobj(src, dst)
        os.remove(tmp)
    print(f"wrote {path}: {len(spectra)} spectra, MS2 scans {ms2_of}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="test-fixtures")
    args = ap.parse_args()
    make_run(False, os.path.join(args.out, "lfq_synthetic.mzML.gz"))
    make_run(True, os.path.join(args.out, "tmt10_synthetic.mzML.gz"))


if __name__ == "__main__":
    main()
