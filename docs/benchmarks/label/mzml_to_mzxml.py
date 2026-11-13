#!/usr/bin/env python3
"""Stream an mzML file into the indexed mzXML layout MaxQuant's mzXML reader expects.

MaxQuant 2.0.3's mzML reader (PluginRawMzMl) keeps only the MS2 isolation
windows that repeat once per cycle, i.e. a DIA/SWATH layout. On DDA mzML every
MS2 spectrum is dropped and the .apl peak lists come out empty. Its mzXML
reader has no such filter. For Thermo files it takes the analyzer and the
profile/centroid state from `filterLine`, so this writes one filterLine per
scan.

    mzml_to_mzxml.py IN.mzML OUT.mzXML [--vendor "Thermo Scientific"] [--analyzer FTMS]
                     [--ms2-from CENTROIDED.mzML] [--ms1-profile-resolution R]

Reads uncompressed or zlib arrays, 32- or 64-bit. Writes 32-bit network-order
m/z-int pairs (MaxQuant 2.0.3 misreads 64-bit mzXML peaks) without
compression, plus the scan index. Empty spectra are left
out (MaxQuant's MS1 peak detection and top-N MS2 filter crash on them) and scans are renumbered 1..n.
"""

import argparse
import base64
import struct
import sys
import xml.etree.ElementTree as ET
import zlib
from xml.sax.saxutils import quoteattr

NS = "{http://psi.hupo.org/ms/mzml}"


def cv(elem, groups=None):
    """The cvParams of `elem` and its descendants, with those of every
    referenced `referenceableParamGroup` (inline ones win)."""
    if elem is None:
        return {}
    out = {}
    for ref in elem.iter(NS + "referenceableParamGroupRef"):
        out.update((groups or {}).get(ref.get("ref"), {}))
    out.update({p.get("accession"): p.get("value") for p in elem.iter(NS + "cvParam")})
    return out


def decode(bda, groups):
    params = cv(bda, groups)
    raw = base64.b64decode(bda.find(NS + "binary").text or "")
    if "MS:1000574" in params:
        raw = zlib.decompress(raw)
    fmt = "d" if "MS:1000523" in params else "f"
    n = len(raw) // struct.calcsize(fmt)
    vals = struct.unpack("<%d%s" % (n, fmt), raw)
    kind = "mz" if "MS:1000514" in params else "int" if "MS:1000515" in params else None
    return kind, vals


def peaks(el, groups):
    arrays = {}
    for bda in el.iter(NS + "binaryDataArray"):
        k, v = decode(bda, groups)
        if k:
            arrays[k] = v
    return arrays


def param_group(el, groups):
    """Record a `referenceableParamGroup` (they precede the spectra)."""
    groups[el.get("id")] = {p.get("accession"): p.get("value") for p in el.iter(NS + "cvParam")}


def spectra(path):
    """(profile, arrays) of every spectrum of an mzML, in file order."""
    groups = {}
    for _, el in ET.iterparse(path, events=("end",)):
        if el.tag == NS + "referenceableParamGroup":
            param_group(el, groups)
        elif el.tag == NS + "spectrum":
            yield "MS:1000128" in cv(el, groups), peaks(el, groups)
            el.clear()


def pseudo_profile(mz, it, resolution):
    """Each centroid as a 5-point Gaussian-like peak of FWHM mz/resolution, zeros
    on its flanks unless a neighbour is that close."""
    pts = []
    for k, (m, i) in enumerate(zip(mz, it)):
        d = m / resolution / 2.0
        if k == 0 or mz[k - 1] + 4 * d < m - 2 * d:
            pts.append((m - 2 * d, 0.0))
        pts += [(m - d, 0.5 * i), (m, i), (m + d, 0.5 * i)]
        if k + 1 == len(mz) or mz[k + 1] - 4 * d > m + 2 * d:
            pts.append((m + 2 * d, 0.0))
    pts.sort()
    return [p[0] for p in pts], [p[1] for p in pts]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--vendor", default="Thermo Scientific")
    ap.add_argument("--analyzer", default="FTMS")
    ap.add_argument("--ms2-from", help="take the MS2 peaks from this mzML (same spectra, same order), "
                    "e.g. the input centroided by OpenMS PeakPickerHiRes")
    ap.add_argument("--ms1-profile-resolution", type=float,
                    help="write centroided MS1 as synthetic profile peaks of this resolving power "
                    "(MaxQuant computes NaN masses from centroided MS1 and crashes)")
    a = ap.parse_args()

    out = open(a.out, "wb")
    pos = 0

    def w(s):
        nonlocal pos
        b = s.encode()
        out.write(b)
        pos += len(b)

    w('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    w('<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">\n')
    hdr = pos
    w('<msRun scanCount="%010d" startTime="PT%014.4fS" endTime="PT%014.4fS">\n' % (0, 0.0, 0.0))
    w('<parentFile fileName=%s fileType="RAWData" fileSha1="0000000000000000000000000000000000000000"/>\n'
      % quoteattr(a.inp))
    w('<msInstrument><msManufacturer category="msManufacturer" value=%s/>'
      '<msModel category="msModel" value="unknown"/><msIonisation category="msIonisation" value="ESI"/>'
      '<msMassAnalyzer category="msMassAnalyzer" value=%s/><msDetector category="msDetector" value="unknown"/></msInstrument>\n'
      % (quoteattr(a.vendor), quoteattr(a.analyzer)))
    w('<dataProcessing centroided="0"><software type="conversion" name="mzml_to_mzxml.py" version="1"/></dataProcessing>\n')

    offsets = []
    num = 0
    t0 = t1 = 0.0
    groups = {}
    other = spectra(a.ms2_from) if a.ms2_from else None
    for ev, el in ET.iterparse(a.inp, events=("end",)):
        if el.tag == NS + "referenceableParamGroup":
            param_group(el, groups)
            continue
        if el.tag != NS + "spectrum":
            continue
        p = cv(el, groups)
        level = int(p.get("MS:1000511", "1"))
        profile = "MS:1000128" in p
        scan = el.find(NS + "scanList/" + NS + "scan")
        sp = cv(scan, groups)
        rt = float(sp.get("MS:1000016") or 0.0)
        # Minutes by the unit accession (unitName is optional).
        if scan is not None:
            for c in scan.iter(NS + "cvParam"):
                if c.get("accession") == "MS:1000016" and (
                        c.get("unitAccession") == "UO:0000031" or c.get("unitName") == "minute"):
                    rt *= 60.0
        arrays = peaks(el, groups)
        if other is not None:
            o_profile, o_arrays = next(other)
            if level == 2:
                profile, arrays = o_profile, o_arrays
        mz = arrays.get("mz", ())
        it = arrays.get("int", ())
        n = min(len(mz), len(it))
        if n == 0:
            # MaxQuant indexes past the end of an empty spectrum (MS1 peak detection,
            # MS2 top-N filter).
            el.clear()
            continue
        if level == 1 and not profile and a.ms1_profile_resolution:
            mz, it = pseudo_profile(mz[:n], it[:n], a.ms1_profile_resolution)
            n = len(mz)
            profile = True
        num += 1
        if n:
            bp = max(range(n), key=lambda i: it[i])
            bpmz, bpi, tic, lo, hi = mz[bp], it[bp], sum(it[:n]), mz[0], mz[n - 1]
        else:
            bpmz = bpi = tic = lo = hi = 0.0
        pc = "p" if profile else "c"
        prec = ""
        ce = ""
        if level == 2:
            pre = el.find(NS + "precursorList/" + NS + "precursor")
            if pre is None:
                pre = ET.Element(NS + "precursor")
            # Every part is optional in mzML: a missing one has no parameters.
            iso = cv(pre.find(NS + "isolationWindow"), groups)
            si = cv(pre.find(NS + "selectedIonList/" + NS + "selectedIon"), groups)
            act = cv(pre.find(NS + "activation"), groups)
            pmz = float(si.get("MS:1000744") or iso.get("MS:1000827"))
            z = si.get("MS:1000041")
            pint = si.get("MS:1000042") or "0"
            method = "HCD" if "MS:1000422" in act else "ETD" if "MS:1000598" in act else "CID"
            energy = float(act.get("MS:1000045") or 0.0)
            ce = ' collisionEnergy="%g"' % energy
            filt = "%s + %s NSI d Full ms2 %.4f@%s%.2f [%.2f-%.2f]" % (
                a.analyzer, pc, pmz, method.lower(), energy, lo, hi)
            zattr = ' precursorCharge="%s"' % z if z else ""
            win = ""
            if "MS:1000828" in iso and "MS:1000829" in iso:
                win = ' windowWideness="%g"' % (float(iso["MS:1000828"]) + float(iso["MS:1000829"]))
            prec = '<precursorMz precursorIntensity="%s"%s activationMethod="%s"%s>%.6f</precursorMz>\n' % (
                pint, zattr, method, win, pmz)
        else:
            filt = "%s + %s NSI Full ms [%.2f-%.2f]" % (a.analyzer, pc, lo, hi)
        pairs = bytearray()
        for i in range(n):
            pairs += struct.pack(">ff", mz[i], it[i])
        offsets.append((num, pos))
        t0 = rt if num == 1 else t0
        t1 = rt
        w('<scan num="%d" scanType="Full" centroided="%d" msLevel="%d" peaksCount="%d" polarity="+" '
          'retentionTime="PT%.4fS"%s lowMz="%.4f" highMz="%.4f" basePeakMz="%.4f" basePeakIntensity="%.4g" '
          'totIonCurrent="%.4g" filterLine=%s>\n'
          % (num, 0 if profile else 1, level, n, rt, ce, lo, hi, bpmz, bpi, tic, quoteattr(filt)))
        w(prec)
        w('<peaks compressionType="none" compressedLen="0" precision="32" byteOrder="network" '
          'contentType="m/z-int">%s</peaks>\n' % base64.b64encode(bytes(pairs)).decode())
        w('</scan>\n')
        el.clear()
        if num % 2000 == 0:
            print(num, file=sys.stderr)
    if num == 0:
        out.close()
        sys.exit(f"{a.inp}: no spectrum with peaks; nothing written")
    w('</msRun>\n')
    idx = pos
    w('<index name="scan">\n')
    for s, o in offsets:
        w('<offset id="%d">%d</offset>\n' % (s, o))
    w('</index>\n<indexOffset>%d</indexOffset>\n<sha1>0</sha1>\n</mzXML>\n' % idx)
    out.seek(hdr)
    out.write(('<msRun scanCount="%010d" startTime="PT%014.4fS" endTime="PT%014.4fS">\n' % (num, t0, t1)).encode())
    out.close()
    print("scans", num, file=sys.stderr)


if __name__ == "__main__":
    main()
