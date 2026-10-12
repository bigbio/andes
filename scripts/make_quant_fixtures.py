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

With `--mbr` it writes three label-free runs for the match-between-runs
end-to-end test instead (no pyOpenMS-specific options beyond the above):

* mbr_run{1,2,3}.mzML.gz — every BSA tryptic peptide of 7-25 residues with at
  most one missed cleavage and no methionine (cysteines carbamidomethylated),
  at charge 2, and the ones of 12 or more residues also at charge 3. Each
  precursor elutes as a Gaussian (sigma 4 s, in run time) at its consensus RT
  mapped through the run's warp: `t`, `1.05 t + 20` and `t + 15 sin(t / 600)`.
  The six peptides of the design above keep their roles: peptide 4
  (HLVDEPQNLIK) is never triggered for MS2 in run 3, and peptide 6
  (DAFLGSFLYEYSR) is absent from MS1 in every run and has an MS2 only in run 1.
  Every other precursor is triggered with probability 0.7 per run, so each run
  has untriggered precursors to transfer. The alignment needs 50 anchors per
  run and the transfer FDR needs 100 training rows per fold, which six
  peptides cannot give.

The exact numbers the tests assert are the constants below.

Usage: scripts/make_quant_fixtures.py [--out test-fixtures] [--mbr]
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


MBR_SIGMA_S = 4.0
MBR_FIRST_APEX_S = 60.0
MBR_LAST_APEX_S = 660.0
MBR_TAIL_S = 60.0
MBR_BACKGROUND_PEAKS = 10
MBR_TRIGGER_PROBABILITY = 0.7
MBR_CHARGE3_MIN_LENGTH = 12
MBR_CHARGE3_FRACTION = 0.4
# Design roles (0-based into PEPTIDES): never triggered in run 3 (index 2),
# and absent from MS1 everywhere with an MS2 in run 1 only.
MBR_UNTRIGGERED = (3, 2)
MBR_MS1_ABSENT = 5
MBR_WARPS = [
    lambda t: t,
    lambda t: 1.05 * t + 20.0,
    lambda t: t + 15.0 * math.sin(t / 600.0),
]


def hcd_activation():
    """The HCD activation enum value (its location differs between pyOpenMS releases)."""
    enum = getattr(p.Precursor, "ActivationMethod", None) or p.ActivationMethod
    return enum.HCD


def bsa_peptides(fasta: str):
    seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">"))
    pieces, cur = [], ""
    for i, aa in enumerate(seq):
        cur += aa
        nxt = seq[i + 1] if i + 1 < len(seq) else ""
        if aa in "KR" and nxt != "P":
            pieces.append(cur)
            cur = ""
    if cur:
        pieces.append(cur)
    out = set()
    for i in range(len(pieces)):
        for mc in (0, 1):
            if i + mc < len(pieces):
                pep = "".join(pieces[i : i + mc + 1])
                if 7 <= len(pep) <= 25 and "M" not in pep:
                    out.add(pep)
    return sorted(out)


def mbr_precursors(fasta: str):
    """(sequence, charge, consensus apex RT, apex intensity of the mono)."""
    rnd = random.Random(23)
    design = [s for s, *_ in PEPTIDES]
    fillers = [s for s in bsa_peptides(fasta) if s not in design]
    rnd.shuffle(fillers)
    peptides = design + fillers
    n = len(peptides)
    slots = [MBR_FIRST_APEX_S + k * (MBR_LAST_APEX_S - MBR_FIRST_APEX_S) / (n - 1) for k in range(n)]
    rnd.shuffle(slots)
    out = []
    for k, seq in enumerate(peptides):
        if k < len(PEPTIDES):
            apex_int = PEPTIDES[k][2]
        else:
            apex_int = 10 ** rnd.uniform(5.3, 6.7)
        out.append((seq, 2, slots[k], apex_int))
        if len(seq) >= MBR_CHARGE3_MIN_LENGTH:
            out.append((seq, 3, slots[k], apex_int * MBR_CHARGE3_FRACTION))
    return out


def make_mbr_run(run: int, precursors, path: str):
    rnd = random.Random(101 + run)
    warp = MBR_WARPS[run]
    exp = p.MSExperiment()
    instrument = exp.getInstrument()
    analyzer = p.MassAnalyzer()
    analyzer.setType(p.MassAnalyzer.AnalyzerType.ORBITRAP)
    instrument.setMassAnalyzers([analyzer])
    instrument.setName("Q Exactive")
    exp.setInstrument(instrument)

    peps = []
    design = [s for s, *_ in PEPTIDES]
    for seq_str, z, apex_c, apex_int in precursors:
        s = seq_str.replace("C", "C(Carbamidomethyl)")
        seq = p.AASequence.fromString(s)
        mono = seq.getMonoWeight(p.Residue.ResidueType.Full, z) / z
        k = design.index(seq_str) if seq_str in design else None
        in_ms1 = k != MBR_MS1_ABSENT
        if k == MBR_MS1_ABSENT:
            triggered = run == 0
        elif (k, run) == MBR_UNTRIGGERED:
            triggered = False
        elif k is not None:
            triggered = True
        else:
            triggered = rnd.random() < MBR_TRIGGER_PROBABILITY
        peps.append((seq, mono, warp(apex_c), apex_int, z, envelope(seq, N_ISOTOPES), in_ms1, triggered))

    end = max(warp(MBR_LAST_APEX_S), MBR_LAST_APEX_S) + MBR_TAIL_S
    scan = 1
    spectra = []
    t = 0.0
    while t <= end:
        ms1 = p.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.setRT(t)
        ms1.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
        ms1.setType(p.SpectrumSettings.SpectrumType.CENTROID)
        peaks = [(rnd.uniform(350.0, 1500.0), rnd.uniform(50.0, 400.0)) for _ in range(MBR_BACKGROUND_PEAKS)]
        for seq, mono, apex_rt, apex_int, z, env, in_ms1, _ in peps:
            if not in_ms1:
                continue
            g = math.exp(-0.5 * ((t - apex_rt) / MBR_SIGMA_S) ** 2)
            if g < 1e-3:
                continue
            for i, e in enumerate(env):
                inten = apex_int * e / env[0] * g
                if inten > 10.0:
                    peaks.append((mono + i * 1.0033548 / z, inten))
        peaks.sort()
        ms1.set_peaks(([mz for mz, _ in peaks], [it for _, it in peaks]))
        spectra.append(ms1)
        scan += 1
        j = 0
        for seq, mono, apex_rt, apex_int, z, env, in_ms1, triggered in peps:
            trigger = apex_rt - MS2_BEFORE_APEX_S
            if not triggered or not (t <= trigger < t + 1.0):
                continue
            j += 1
            ms2 = p.MSSpectrum()
            ms2.setMSLevel(2)
            ms2.setRT(t + 0.05 * j)
            ms2.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
            ms2.setType(p.SpectrumSettings.SpectrumType.CENTROID)
            pre = p.Precursor()
            pre.setMZ(mono)
            pre.setCharge(z)
            pre.setIntensity(apex_int)
            pre.setIsolationWindowLowerOffset(0.7)
            pre.setIsolationWindowUpperOffset(0.7)
            pre.setActivationMethods({hcd_activation()})
            ms2.setPrecursors([pre])
            frag = sorted(fragments(seq, z))
            ms2.set_peaks(([mz for mz, _ in frag], [it for _, it in frag]))
            spectra.append(ms2)
            scan += 1
        t += 1.0

    exp.setSpectra(spectra)
    tmp = path[:-3]
    opts = p.PeakFileOptions()
    opts.setIntensity32Bit(True)
    opts.setCompression(True)
    f = p.MzMLFile()
    f.setOptions(opts)
    f.store(tmp, exp)
    with open(tmp, "rb") as src, gzip.open(path, "wb", compresslevel=9) as dst:
        shutil.copyfileobj(src, dst)
    os.remove(tmp)
    n_ms2 = sum(1 for x in peps if x[7])
    print(f"wrote {path}: {len(spectra)} spectra, {len(peps)} precursors, {n_ms2} triggered")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="test-fixtures")
    ap.add_argument("--mbr", action="store_true", help="write the three MBR runs only")
    args = ap.parse_args()
    if args.mbr:
        precursors = mbr_precursors(os.path.join(args.out, "BSA.fasta"))
        for run in range(len(MBR_WARPS)):
            make_mbr_run(run, precursors, os.path.join(args.out, f"mbr_run{run + 1}.mzML.gz"))
        return
    make_run(False, os.path.join(args.out, "lfq_synthetic.mzML.gz"))
    make_run(True, os.path.join(args.out, "tmt10_synthetic.mzML.gz"))


if __name__ == "__main__":
    main()
