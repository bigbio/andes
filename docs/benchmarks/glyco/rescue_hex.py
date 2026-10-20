#!/usr/bin/env python3
"""The reference Hex diagnostic-ion rescue, applied to andes's Hex-containing GPSMs.

`eval_hex_entrap.py` reports the raw Hex-containing rate — an UPPER bound on the
error, because it counts every Hex-containing O-GPSM as wrong without the
diagnostic-ion rescue reference engine P applies. This script applies that rescue so the two
numbers are directly comparable.

Engine P counts a Hex-containing O-glycopeptide as an error ONLY if its Hex
diagnostic ions fail validation. The exact criterion: the summed intensity of the two Hex diagnostic ions
  m/z 163.060  (Hex oxonium, [Hex+H-H2O]+)
  m/z 366.139  (Hex-HexNAc oxonium, [HexHexNAc+H]+)
must reach >= 10% of the base-peak intensity.

The SpectId is parsed for `frac` (the prefix before `_controllerType`, added by
`pool_pins.py`) and `scan`, and the raw `.mgf` for that fraction is streamed to
fetch the spectrum. `mgfdir` must contain `<frac>.mgf` for each fraction present.

Usage:
    rescue_hex.py <pooled.glyco.pin> <percolator.psms> <mgfdir>
"""
import os
import re
import sys

if len(sys.argv) < 4:
    sys.exit(__doc__)
pin, psms, mgfdir = sys.argv[1], sys.argv[2], sys.argv[3]

COMP = re.compile(r"HexNAc(\d+)Hex(\d+)Fuc(\d+)NeuAc(\d+)NeuGc(\d+)")
HEX_IONS = (163.0604, 366.1395)
TOL_PPM = 20.0
THRESH = 0.10  # engine P: sum(Hex ion intensities) >= 10% of base peak


def comp_of(pep):
    m = COMP.search(pep)
    return tuple(int(x) for x in m.groups()) if m else None


# SpecId -> glycan composition (Label=1 rows only).
rows = {}
hdr = None
with open(pin) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            ip, il = hdr.index("Peptide"), hdr.index("Label")
            continue
        try:
            if int(f[il]) != 1:
                continue
        except (ValueError, IndexError):
            continue
        c = comp_of(f[ip])
        if c:
            rows[f[0]] = c

accepted = []
total_accepted = 0
need = {}
with open(psms) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if f[0] == "PSMId" or len(f) < 3:
            continue
        try:
            if float(f[2]) > 0.01:
                continue
        except ValueError:
            continue
        c = rows.get(f[0])
        if not c:
            continue
        total_accepted += 1
        if c[1] == 0:  # Hex == 0: not a Hex-containing candidate
            continue
        frac = f[0].split("_controllerType")[0]
        scan = int(re.search(r"scan=(\d+)", f[0]).group(1))
        accepted.append((frac, scan, c))
        need.setdefault(frac, set()).add(scan)

# Stream each .mgf once, harvesting only the needed spectra.
spectra = {}
for frac, scans in need.items():
    path = os.path.join(mgfdir, frac + ".mgf")
    cur_scan = None
    cur_peaks = []
    grab = False
    found = set()
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith("BEGIN IONS"):
                cur_scan = None
                cur_peaks = []
                grab = False
            elif ln.startswith(("TITLE=", "SCANS=")):
                m = re.search(r"scan=(\d+)", ln)
                if m:
                    cur_scan = int(m.group(1))
                    grab = cur_scan in scans and cur_scan not in found
            elif grab and ln and not ln.startswith(("END", "CHARGE", "PEPMASS", "RTINSECONDS")):
                p = ln.split()
                if len(p) >= 2:
                    try:
                        cur_peaks.append((float(p[0]), float(p[1])))
                    except ValueError:
                        pass
            elif ln.startswith("END IONS") and grab and cur_scan is not None:
                spectra[(frac, cur_scan)] = cur_peaks
                found.add(cur_scan)
                grab = False


def check(peaks):
    """Return (verdict, summed_frac, base, matched_ions) for one spectrum."""
    if not peaks:
        return ("ERROR", 0.0, 0.0, None)
    base = max(i for _, i in peaks)
    got = []
    for ion in HEX_IONS:
        best = min(peaks, key=lambda p: abs(p[0] - ion))
        if abs(best[0] - ion) <= ion * TOL_PPM * 1e-6:
            got.append((ion, best[0], best[1]))
    if not got:
        return ("ERROR", 0.0, base, None)
    frac = sum(g[2] for g in got) / base
    verdict = "RESCUED" if frac >= THRESH else "ERROR"
    return (verdict, frac, base, got)


rescued = presence = 0
detail = []
for frac, scan, c in accepted:
    peaks = spectra.get((frac, scan))
    if peaks is None:
        detail.append((frac, scan, c, "NO-SPECTRUM", 0.0))
        continue
    verdict, frac_i, _base, got = check(peaks)
    if got is not None:
        presence += 1
    if verdict == "RESCUED":
        rescued += 1
    detail.append((frac, scan, c, verdict, frac_i))

total = len(accepted)
print(f"accepted O-GPSMs @1%          : {total_accepted}")
print(f"Hex-containing accepted GPSMs : {total}")
print(f"  ion present (any of 163/366) : {presence}")
print(f"  rescued (>=10% of base peak) : {rescued}")
print(f"  unresolved (Hex error)       : {total - rescued}")
print(
    f"  andes Hex-error (engine P rule): {total - rescued}/{total_accepted} "
    f"= {100 * (total - rescued) / total_accepted:.2f}%"
)
print("  engine P: 2/484 = 0.41%; other published engines 2.9%, 10%, 20%")
print()
for frac, scan, c, v, fr in detail:
    print(f"  {frac} scan={scan} {c} -> {v}  (Hex-ion sum/base = {100 * fr:.1f}%)")
