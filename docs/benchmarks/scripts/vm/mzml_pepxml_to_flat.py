#!/usr/bin/env python3
"""MSFragger pepXML + mzML → flat training parquet (Lumos CID self-labels)."""
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

UNIMOD = {
    "Carbamidomethyl": 57.021464,
    "Oxidation": 15.994915,
    "TMT6plex": 229.162932,
    "TMTpro": 304.207146,
}
PROTON = 1.007276

# Monoisotopic residue masses (Da) — used to convert MSFragger's absolute
# `mod_aminoacid_mass` (residue + mod) into the mod *delta* the flat schema uses.
AA_MONO = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
    "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
NTERM_H = 1.0078250319    # standard peptide N-terminus (H), folded into mod_nterm_mass
CTERM_OH = 17.0027396542  # standard peptide C-terminus (OH), folded into mod_cterm_mass


def local(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


def decode_binary(b64: str, is64: bool, zlibc: bool):
    """Decode an mzML <binary> array. Handles zlib compression and 32/64-bit
    floats (ThermoRawFileParser emits 64-bit float + zlib by default, which the
    old fixed `<n>f` unpack could not read)."""
    import base64
    import struct
    import zlib

    raw = base64.b64decode(b64)
    if zlibc:
        raw = zlib.decompress(raw)
    if is64:
        return struct.unpack(f"{len(raw) // 8}d", raw)
    return struct.unpack(f"{len(raw) // 4}f", raw)


def parse_mods(modinfo, seq: str):
    """Return (res_pos, res_delta, nterm_delta, cterm_delta) from a pepXML
    `<modification_info>` element.

    Residue deltas are the absolute `mod_aminoacid_mass` minus the residue's
    monoisotopic mass, so fixed (C+57) and variable mods both appear as deltas —
    matching the MSnet flat convention. Terminal deltas strip the standard H/OH
    that MSFragger folds into `mod_{n,c}term_mass`. Positions are 1-based.
    """
    res_pos: list[int] = []
    res_delta: list[float] = []
    nterm = cterm = 0.0
    if modinfo is None:
        return res_pos, res_delta, nterm, cterm
    nt = modinfo.get("mod_nterm_mass")
    if nt is not None:
        nterm = float(nt) - NTERM_H
    ct = modinfo.get("mod_cterm_mass")
    if ct is not None:
        cterm = float(ct) - CTERM_OH
    for child in modinfo:
        if local(child.tag) != "mod_aminoacid_mass":
            continue
        pos = int(child.get("position", "0"))
        if pos < 1 or pos > len(seq):
            continue
        base = AA_MONO.get(seq[pos - 1])
        if base is None:
            continue
        delta = float(child.get("mass", "0")) - base
        if abs(delta) < 1e-6:
            continue
        res_pos.append(pos)
        res_delta.append(delta)
    return res_pos, res_delta, nterm, cterm


def load_ms2_by_scan(mzml: Path, need_scans: set[int]) -> dict[int, tuple[float, int, list, list]]:
    """Iterparse mzML; return scan -> (prec_mz, charge, mz[], intensity[])."""
    out: dict[int, tuple[float, int, list, list]] = {}
    if not need_scans:
        return out
    ctx = ET.iterparse(mzml, events=("end",))
    for _event, elem in ctx:
        if local(elem.tag) != "spectrum":
            continue
        ms_level = None
        # Scan number: ThermoRawFileParser/Thermo mzML carries it in the
        # spectrum `id` ("controllerType=0 controllerNumber=1 scan=N"), NOT a
        # "scan number" cvParam. Key on the id (matching the Rust reader and the
        # pepXML `start_scan`); the cvParam below stays as a fallback for other
        # converters.
        scan_num = None
        sid_match = re.search(r"scan=(\d+)", elem.get("id", ""))
        if sid_match:
            scan_num = int(sid_match.group(1))
        prec_mz = 0.0
        charge = 0
        mzs: list[float] = []
        ints: list[float] = []
        for child in elem:
            ln = local(child.tag)
            if ln == "cvParam" and child.get("name") == "ms level":
                ms_level = int(child.get("value", "0"))
            if ln == "cvParam" and child.get("name") == "scan number":
                scan_num = int(child.get("value", "0"))
            if ln == "precursorList":
                for prec in child.iter():
                    if local(prec.tag) == "cvParam" and prec.get("name") == "charge state":
                        charge = int(float(prec.get("value", "0")))
                    if local(prec.tag) == "cvParam" and prec.get("accession") == "MS:1000744":
                        prec_mz = float(prec.get("value", "0"))
            if ln == "binaryDataArrayList":
                arr_mz = arr_i = None
                for arr in child:
                    if local(arr.tag) != "binaryDataArray":
                        continue
                    kind = None
                    is64 = False
                    zlibc = False
                    b64 = None
                    for a in arr:
                        al = local(a.tag)
                        if al == "cvParam":
                            nm = a.get("name")
                            if nm == "m/z array":
                                kind = "mz"
                            elif nm == "intensity array":
                                kind = "int"
                            elif nm == "64-bit float":
                                is64 = True
                            elif nm == "zlib compression":
                                zlibc = True
                        if al == "binary":
                            b64 = a.text
                    if kind and b64:
                        vals = decode_binary(b64, is64, zlibc)
                        if kind == "mz":
                            arr_mz = vals
                        else:
                            arr_i = vals
                if arr_mz and arr_i:
                    mzs = list(arr_mz)
                    ints = list(arr_i)
        if ms_level == 2 and scan_num in need_scans and len(mzs) >= 10 and charge >= 2:
            out[scan_num] = (prec_mz, charge, mzs, ints)
            if len(out) >= len(need_scans):
                break
        elem.clear()
    return out


def parse_pepxml_hits(pepxml: Path, expect_max: float, max_psms: int):
    hits = []
    ctx = ET.iterparse(pepxml, events=("end",))
    for _event, elem in ctx:
        if local(elem.tag) != "search_hit":
            continue
        expect = float(elem.get("expect", "999"))
        if expect > expect_max:
            elem.clear()
            continue
        pep = elem.get("peptide", "")
        charge = int(elem.get("assumed_charge", "0"))
        mod_pep = pep
        for child in elem:
            if local(child.tag) == "modification_info":
                mod_pep = child.get("modified_peptide", pep)
        # parent spectrum_query
        parent = elem
        # walk up isn't easy with iterparse; spectrum attr on search_hit missing
        # MSFragger puts spectrum on spectrum_query; capture via stack in full parse
        elem.clear()
        if pep and charge >= 2:
            hits.append((pep, charge, mod_pep, expect))
        if len(hits) >= max_psms:
            break
    return hits


def parse_pepxml_full(pepxml: Path, expect_max: float, max_psms: int):
    tree = ET.parse(pepxml)
    root = tree.getroot()
    hits = []
    for sq in root.iter():
        if local(sq.tag) != "spectrum_query":
            continue
        # Scan key: MSFragger's `spectrum` attr is `base.SCAN.SCAN.charge`
        # (no `scan=`), so the old `scan=(\d+)` regex matched nothing and the
        # flat came out empty. Prefer the integer `start_scan` attribute; fall
        # back to `scan=` in `spectrumNativeID`, then the `.SCAN.SCAN.charge`
        # spectrum string. Must yield the same integer the mzML side keys on.
        scan = None
        ss = sq.get("start_scan")
        if ss and ss.isdigit():
            scan = int(ss)
        else:
            nid = sq.get("spectrumNativeID", "")
            m = re.search(r"scan=(\d+)", nid)
            if m is None:
                m = re.search(r"\.(\d+)\.\d+\.\w+$", sq.get("spectrum", ""))
            if m:
                scan = int(m.group(1))
        if scan is None:
            continue
        # `assumed_charge` is on the spectrum_query, not the search_hit.
        charge = int(sq.get("assumed_charge", "0"))
        # search_hit is nested under <search_result>, so descend with iter()
        # rather than iterating direct children of spectrum_query.
        for sh in sq.iter():
            if local(sh.tag) != "search_hit":
                continue
            if sh.get("hit_rank", "1") != "1":
                continue  # keep only the best hit per scan as the self-label
            pep = sh.get("peptide", "")
            # `expect` is a `<search_score name="expect">` child, NOT an
            # attribute; mods come from the structured `<modification_info>`.
            expect = 999.0
            modinfo = None
            for child in sh:
                cl = local(child.tag)
                if cl == "search_score" and child.get("name") == "expect":
                    expect = float(child.get("value", "999"))
                elif cl == "modification_info":
                    modinfo = child
            if expect > expect_max:
                continue
            seq = re.sub(r"[^A-Z]", "", pep)
            if seq and charge >= 2:
                hits.append((scan, seq, charge, parse_mods(modinfo, seq), expect))
            if len(hits) >= max_psms:
                return hits
    return hits


def main() -> None:
    mzml = Path(sys.argv[1])
    pepxml = Path(sys.argv[2])
    out = Path(sys.argv[3])
    expect_max = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
    max_psms = int(sys.argv[5]) if len(sys.argv) > 5 else 50000

    hits = parse_pepxml_full(pepxml, expect_max, max_psms)
    scans = {h[0] for h in hits}
    print(f"pepXML hits={len(hits)} unique_scans={len(scans)} expect<={expect_max}")
    spectra = load_ms2_by_scan(mzml, scans)
    print(f"mzML matched scans={len(spectra)}")

    schema = pa.schema(
        [
            ("seq", pa.string()),
            ("charge", pa.int32()),
            ("prec_mz", pa.float64()),
            ("res_mod_pos", pa.list_(pa.int32())),
            ("res_mod_delta", pa.list_(pa.float64())),
            ("nterm_delta", pa.float64()),
            ("cterm_delta", pa.float64()),
            ("mz", pa.list_(pa.float32())),
            ("intensity", pa.list_(pa.float32())),
        ]
    )
    c = {k: [] for k in schema.names}
    kept = 0
    for scan, seq, charge, mods, _exp in hits:
        if scan not in spectra:
            continue
        prec_mz, ch, mzs, ints = spectra[scan]
        if ch != charge:
            charge = ch
        if len(seq) < 6:
            continue
        rp, rd, nt, ct = mods
        c["seq"].append(seq)
        c["charge"].append(int(charge))
        c["prec_mz"].append(float(prec_mz))
        c["res_mod_pos"].append(rp)
        c["res_mod_delta"].append(rd)
        c["nterm_delta"].append(nt)
        c["cterm_delta"].append(ct)
        c["mz"].append([float(x) for x in mzs])
        c["intensity"].append([float(x) for x in ints])
        kept += 1

    pq.write_table(pa.table(c, schema=schema), out)
    print(f"kept={kept} -> {out}")


if __name__ == "__main__":
    main()
