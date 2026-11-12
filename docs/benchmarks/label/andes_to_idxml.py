#!/usr/bin/env python3
"""Write andes PSMs as one idXML per run, for OpenMS MS1LabeledWorkflow.

    andes_to_idxml.py Q_TSV PIN RUN MZML OUT.idXML

Q_TSV is andes's `<stem>.q0p01.tsv` (PSMs at q <= 0.01; PSMId
`<run>/<native id>_<scan>_<rank>[_<row>]`, peptide `X.PEPTIDE.X` with mass
deltas after the residue), PIN its `--output-pin` file (the charge of each
PSM), RUN the run name in the PSMIds, MZML that run's spectra (retention time
and precursor m/z per native id). The score is the
andes posterior error probability, as MS1LabeledWorkflow needs; the mass deltas
become named Unimod modifications, so that the workflow sees the SILAC labels.
"""
import re
import sys

import pyopenms as oms

MODS = {
    ("K", "+8.01420"): "Label:13C(6)15N(2)",
    ("R", "+10.00827"): "Label:13C(6)15N(4)",
    ("K", "+6.02013"): "Label:13C(6)",
    ("R", "+6.02013"): "Label:13C(6)",
    ("K", "+4.02511"): "Label:2H(4)",
    ("C", "+57.02146"): "Carbamidomethyl",
    ("M", "+15.99491"): "Oxidation",
}
ACETYL = "+42.01056"
TOKEN = re.compile(r"([A-Z])((?:[+-][0-9]+\.[0-9]+)*)")


def to_openms(pep):
    # strip the flanking residues: andes writes `K.PEPTIDE.A`
    core = pep
    m = re.match(r"^[A-Z\-]\.(.*)\.[A-Z\-]$", pep)
    if m:
        core = m.group(1)
    out = []
    nterm = ""
    for k, (aa, deltas) in enumerate(TOKEN.findall(core)):
        names = []
        for d in re.findall(r"[+-][0-9]+\.[0-9]+", deltas):
            if d == ACETYL and k == 0:
                nterm = ".(Acetyl)"
            elif (aa, d) in MODS:
                names.append(MODS[(aa, d)])
            elif k == 0 and (aa, "%+.5f" % (float(d) - float(ACETYL))) in MODS:
                # protein N-terminal acetyl merged with the first residue's own modification
                nterm = ".(Acetyl)"
                names.append(MODS[(aa, "%+.5f" % (float(d) - float(ACETYL)))])
            else:
                raise ValueError(f"unknown modification {aa}{d} in {pep}")
        out.append(aa + "".join(f"({n})" for n in names))
    return nterm + "".join(out)


def main():
    q_tsv, pin, run, mzml, out = sys.argv[1:6]
    charges = {}
    with open(pin) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        zcols = [(i, int(c[6:])) for i, c in enumerate(header) if re.fullmatch(r"charge\d+", c)]
        for line in fh:
            r = line.split("\t", zcols[-1][0] + 1)
            if not r[0].startswith(run + "/"):
                continue
            z = [c for i, c in zcols if r[i] == "1"]
            charges[r[0]] = z[0] if z else 1
    exp = oms.MSExperiment()
    opts = oms.PeakFileOptions()
    opts.setFillData(False)
    opts.setMSLevels([2])
    f = oms.MzMLFile()
    f.setOptions(opts)
    f.load(mzml, exp)
    meta = {}
    for s in exp.getSpectra():
        p = s.getPrecursors()
        if p:
            meta[s.getNativeID()] = (s.getRT(), p[0].getMZ())

    peps = []
    skipped = 0
    with open(q_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            sid = r[col["PSMId"]]
            if not sid.startswith(run + "/"):
                continue
            native = re.match(r"^(.*?)_\d+_\d+(?:_\d+)?$", sid[len(run) + 1:]).group(1)
            charge = charges[sid]
            rt, mz = meta[native]
            pid = oms.PeptideIdentification()
            pid.setRT(rt)
            pid.setMZ(mz)
            pid.setScoreType("Posterior Error Probability")
            pid.setHigherScoreBetter(False)
            pid.setMetaValue("spectrum_reference", native)
            try:
                seq = to_openms(r[col["peptide"]])
            except ValueError as e:
                skipped += 1
                print("skipped:", e, file=sys.stderr)
                continue
            h = oms.PeptideHit()
            h.setSequence(oms.AASequence.fromString(seq))
            h.setCharge(int(charge))
            h.setScore(float(r[col["posterior_error_probability"]]))
            h.setMetaValue("target_decoy", "target")
            h.setMetaValue("q-value", float(r[col["q-value"]]))
            evs = []
            for acc in r[col["proteinIds"]].split(","):
                ev = oms.PeptideEvidence()
                ev.setProteinAccession(acc)
                evs.append(ev)
            h.setPeptideEvidences(evs)
            pid.setHits([h])
            pid.setIdentifier("andes")
            peps.append(pid)

    prot = oms.ProteinIdentification()
    prot.setIdentifier("andes")
    prot.setSearchEngine("andes")
    prot.setPrimaryMSRunPath([mzml.encode()])
    accs = sorted({e.getProteinAccession() for p in peps for e in p.getHits()[0].getPeptideEvidences()})
    hits = []
    for a in accs:
        ph = oms.ProteinHit()
        ph.setAccession(a)
        ph.setMetaValue("target_decoy", "target")
        hits.append(ph)
    prot.setHits(hits)
    sp = prot.getSearchParameters()
    sp.fixed_modifications = [b"Carbamidomethyl (C)"]
    sp.variable_modifications = [b"Oxidation (M)", b"Acetyl (Protein N-term)",
                                 b"Label:13C(6)15N(2) (K)", b"Label:13C(6)15N(4) (R)"]
    sp.digestion_enzyme = oms.ProteaseDB().getEnzyme("Trypsin/P")
    sp.missed_cleavages = 2
    prot.setSearchParameters(sp)
    oms.IdXMLFile().store(out, [prot], peps)
    print(run, len(peps), "PSMs,", skipped, "skipped (unknown modification)", file=sys.stderr)


if __name__ == "__main__":
    main()
