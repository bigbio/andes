"""Streaming mzML reader for selected MS2 scans (standard library + numpy only).

Yields one dict per wanted scan: ``scan``, ``ms_level``, ``filter``,
``collision_energy``, ``precursor_mz``, ``mz`` (sorted float64), ``intensity``
(float64). Binary arrays are decoded only for the requested scans. Handles
plain and gzip-compressed (``.mzML.gz``) files, 32/64-bit floats, zlib or no
compression (numpress is rejected).
"""
from __future__ import annotations

import base64
import gzip
import re
import zlib
import xml.etree.ElementTree as ET

import numpy as np

_SCAN = re.compile(r"scan=(\d+)")
_ACC = {
    "MS:1000511": "ms_level",
    "MS:1000512": "filter",
    "MS:1000045": "collision_energy",
    "MS:1000744": "precursor_mz",
    "MS:1000041": "precursor_charge",
}


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _decode(bda) -> tuple[str | None, np.ndarray]:
    dtype, zlibbed, kind, text = np.float64, False, None, None
    for ch in bda:
        t = _local(ch.tag)
        if t == "cvParam":
            acc = ch.get("accession")
            if acc == "MS:1000521":
                dtype = np.float32
            elif acc == "MS:1000523":
                dtype = np.float64
            elif acc == "MS:1000574":
                zlibbed = True
            elif acc == "MS:1000576":
                zlibbed = False
            elif acc in ("MS:1002312", "MS:1002313", "MS:1002314"):
                raise ValueError("numpress-compressed mzML is not supported")
            elif acc == "MS:1000514":
                kind = "mz"
            elif acc == "MS:1000515":
                kind = "intensity"
        elif t == "binary":
            text = ch.text or ""
    raw = base64.b64decode(text) if text else b""
    if zlibbed and raw:
        raw = zlib.decompress(raw)
    return kind, np.frombuffer(raw, dtype=dtype).astype(np.float64)


def _open(path: str):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def instrument_name(path: str) -> str | None:
    """First instrument-model cvParam name in the file header, if any."""
    from ms2_common import instrument_class

    names = []
    with _open(path) as fh:
        for ev, el in ET.iterparse(fh, events=("start",)):
            t = _local(el.tag)
            if t == "cvParam" and el.get("accession", "").startswith("MS:"):
                names.append(el.get("name", ""))
            if t == "run":
                break
    for n in names:
        if instrument_class(n) > 0:
            return n
    return None


def read_scans(path: str, wanted: set[int] | None = None):
    """Yield wanted MS2 scans in file order (all MS2 scans when ``wanted`` is None)."""
    with _open(path) as fh:
        it = ET.iterparse(fh, events=("start", "end"))
        slist = None
        for ev, el in it:
            t = _local(el.tag)
            if ev == "start":
                if t == "spectrumList":
                    slist = el
                continue
            if t != "spectrum":
                continue
            m = _SCAN.search(el.get("id", ""))
            scan = int(m.group(1)) if m else int(el.get("index", "-1")) + 1
            if wanted is None or scan in wanted:
                rec = {"scan": scan, "ms_level": None, "filter": None, "collision_energy": None,
                       "precursor_mz": None, "precursor_charge": None}
                arrays = {}
                for sub in el.iter():
                    st = _local(sub.tag)
                    if st == "cvParam":
                        key = _ACC.get(sub.get("accession"))
                        if key and rec[key] is None:
                            rec[key] = sub.get("value")
                    elif st == "binaryDataArray":
                        kind, arr = _decode(sub)
                        if kind:
                            arrays[kind] = arr
                lvl = int(rec["ms_level"]) if rec["ms_level"] else None
                if lvl == 2 and "mz" in arrays and "intensity" in arrays:
                    mz, inten = arrays["mz"], arrays["intensity"]
                    order = np.argsort(mz, kind="stable")
                    rec["ms_level"] = lvl
                    rec["mz"], rec["intensity"] = mz[order], inten[order]
                    for k in ("collision_energy", "precursor_mz"):
                        rec[k] = float(rec[k]) if rec[k] is not None else float("nan")
                    rec["precursor_charge"] = int(rec["precursor_charge"]) if rec["precursor_charge"] else 0
                    yield rec
            el.clear()
            if slist is not None:
                slist.clear()
