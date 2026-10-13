"""Shared definitions for andes-ms2: ion channels, fragment masses, b/y annotation,
acquisition metadata and the parquet shard schema.

Channel layout per cleavage site (``N_ION`` = 12), site ``i`` = cleavage after
residue ``i + 1`` (0-based), b = N-terminal fragment of length ``i + 1``,
y = the complementary C-terminal fragment::

    0 b1   1 b2   2 y1   3 y2                      plain b/y, fragment charge 1/2
    4..7   the same four, minus H3PO4 (97.9769)    only if the fragment holds pS/pT
    8..11  the same four, minus CH4SO (63.9983)    only if the fragment holds oxidised M

The first four channels are the PeptDeep ``b_z1, b_z2, y_z1, y_z2`` layout used by
``examples/rescore_pin.rs``. A channel is *invalid* (stored as -1) when its
fragment charge exceeds the precursor charge or, for a loss channel, when the
fragment carries no residue that can lose that group. Valid channels hold the
observed intensity divided by the most intense valid annotated ion, so values lie
in [0, 1].

``onnx-infer::ms2_native`` mirrors ``ION_TYPES``, ``ACTIVATIONS``,
``INSTRUMENTS``, ``instrument_class`` and ``valid_mask``.
"""
from __future__ import annotations

import hashlib
import os
import re
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from mod_composition import MASSES, UNIMOD  # noqa: E402

PROTON = 1.007276466812
H2O = 18.0105646837
H3PO4 = 97.976895573
CH4SO = 63.998285
LOSS_MASS = {None: 0.0, "H3PO4": H3PO4, "CH4SO": CH4SO}

# (ion, fragment charge, neutral loss)
ION_TYPES = [
    ("b", 1, None), ("b", 2, None), ("y", 1, None), ("y", 2, None),
    ("b", 1, "H3PO4"), ("b", 2, "H3PO4"), ("y", 1, "H3PO4"), ("y", 2, "H3PO4"),
    ("b", 1, "CH4SO"), ("b", 2, "CH4SO"), ("y", 1, "CH4SO"), ("y", 2, "CH4SO"),
]
N_ION = len(ION_TYPES)
N_BASE_ION = 4

# Activation = dissociation + MS2 analyser; index 0 = unknown.
ACTIVATIONS = ["unknown", "HCD-FTMS", "HCD-ITMS", "CID-ITMS", "CID-FTMS", "HCD-ASTMS"]
# Instrument family; index 0 = unknown.
INSTRUMENTS = ["unknown", "QExactive", "Exploris", "Tribrid", "VelosElite", "Astral", "LTQ"]
# Precursor charges 1..MAX_CHARGE are one-hot encoded; higher charges use the last slot.
MAX_CHARGE = 6

RESIDUE_MASS = {
    "G": 57.02146372, "A": 71.03711379, "S": 87.03202841, "P": 97.05276385,
    "V": 99.06841391, "T": 101.0476785, "C": 103.0091845, "L": 113.084064,
    "I": 113.084064, "N": 114.0429274, "D": 115.026943, "Q": 128.0585775,
    "K": 128.094963, "E": 129.0425931, "M": 131.0404846, "H": 137.0589119,
    "F": 147.0684139, "R": 156.101111, "Y": 163.0633286, "W": 186.0793129,
}

PHOSPHO = 21
OXIDATION = 35


def instrument_class(name: str | None) -> int:
    """Index in ``INSTRUMENTS`` for a free-text instrument model name."""
    if not name:
        return 0
    s = name.lower().replace(" ", "").replace("-", "")
    if "astral" in s:
        return INSTRUMENTS.index("Astral")
    if "exploris" in s:
        return INSTRUMENTS.index("Exploris")
    if "exactive" in s:
        return INSTRUMENTS.index("QExactive")
    if any(k in s for k in ("fusion", "lumos", "eclipse", "ascend", "tribrid")):
        return INSTRUMENTS.index("Tribrid")
    if "velos" in s or "elite" in s:
        return INSTRUMENTS.index("VelosElite")
    if "ltq" in s:
        return INSTRUMENTS.index("LTQ")
    return 0


def activation_index(name: str | None) -> int:
    return ACTIVATIONS.index(name) if name in ACTIVATIONS else 0


def split_of(seq: str) -> str:
    """train / val / test by an md5 of the stripped sequence (same rule as andes-rt; I as L)."""
    h = int(hashlib.md5(seq.replace("I", "L").encode()).hexdigest()[:8], 16) % 100
    return "test" if h < 5 else "val" if h < 10 else "train"


_PF = re.compile(r"\[UNIMOD:(\d+)\]")


def parse_proforma(s: str):
    """``[UNIMOD:737]-PEPS[UNIMOD:21]K`` -> (sequence, [(pos, unimod id)]); pos 0 = N-term,
    n + 1 = C-term. Returns None for anything else (named mods, mass tags, ambiguity)."""
    mods = []
    i = 0
    if s.startswith("["):
        j = s.find("]-")
        if j < 0:
            return None
        head = s[: j + 1]
        ids = _PF.findall(head)
        if "".join(f"[UNIMOD:{m}]" for m in ids) != head:
            return None
        mods += [(0, int(m)) for m in ids]
        i = j + 2
    seq = []
    while i < len(s):
        c = s[i]
        if c == "-":
            tail = s[i + 1:]
            ids = _PF.findall(tail)
            if not ids or "".join(f"[UNIMOD:{m}]" for m in ids) != tail:
                return None
            mods += [(len(seq) + 1, int(m)) for m in ids]
            break
        if c == "[":
            m = _PF.match(s, i)
            if m is None or not seq:
                return None
            mods.append((len(seq), int(m.group(1))))
            i = m.end()
            continue
        if not ("A" <= c <= "Z"):
            return None
        seq.append(c)
        i += 1
    if not seq:
        return None
    return "".join(seq), sorted(mods)


def to_proforma(seq: str, mods) -> str:
    tag = {}
    for pos, uid in mods:
        tag.setdefault(pos, []).append(f"[UNIMOD:{uid}]")
    out = "".join(tag.get(0, []))
    if out:
        out += "-"
    for i, c in enumerate(seq, 1):
        out += c + "".join(tag.get(i, []))
    ct = "".join(tag.get(len(seq) + 1, []))
    return out + ("-" + ct if ct else "")


def supported(seq: str, mods) -> bool:
    return all(c in RESIDUE_MASS for c in seq) and all(uid in UNIMOD for _, uid in mods)


def residue_masses(seq: str, mods):
    """Per-residue masses with mods; N-term mods on residue 1, C-term mods on residue n."""
    n = len(seq)
    m = np.array([RESIDUE_MASS[c] for c in seq], dtype=np.float64)
    for pos, uid in mods:
        m[min(max(pos, 1), n) - 1] += MASSES[uid]
    return m


def precursor_mz(seq: str, mods, charge: int) -> float:
    return (residue_masses(seq, mods).sum() + H2O + charge * PROTON) / charge


def loss_carriers(seq: str, mods):
    """Per residue: (carries pS/pT, carries oxidised M)."""
    n = len(seq)
    ph = np.zeros(n, bool)
    ox = np.zeros(n, bool)
    for pos, uid in mods:
        if 1 <= pos <= n:
            if uid == PHOSPHO and seq[pos - 1] in "ST":
                ph[pos - 1] = True
            if uid == OXIDATION and seq[pos - 1] == "M":
                ox[pos - 1] = True
    return ph, ox


def valid_mask(seq: str, mods, charge: int) -> np.ndarray:
    """(n - 1, N_ION) bool: which channels exist for this precursor."""
    n = len(seq)
    ph, ox = loss_carriers(seq, mods)
    ph_pre, ox_pre = np.cumsum(ph), np.cumsum(ox)
    v = np.zeros((n - 1, N_ION), bool)
    for i in range(n - 1):
        has_b = {"H3PO4": ph_pre[i] > 0, "CH4SO": ox_pre[i] > 0}
        has_y = {"H3PO4": ph_pre[-1] - ph_pre[i] > 0, "CH4SO": ox_pre[-1] - ox_pre[i] > 0}
        for k, (ion, z, loss) in enumerate(ION_TYPES):
            if z > charge:
                continue
            if loss is None:
                v[i, k] = True
            else:
                v[i, k] = has_b[loss] if ion == "b" else has_y[loss]
    return v


def fragment_mz(seq: str, mods) -> np.ndarray:
    """(n - 1, N_ION) theoretical m/z for every channel (valid or not)."""
    res = residue_masses(seq, mods)
    total = res.sum()
    prefix = np.cumsum(res)[:-1]
    out = np.zeros((len(prefix), N_ION))
    for k, (ion, z, loss) in enumerate(ION_TYPES):
        neutral = prefix if ion == "b" else total - prefix + H2O
        out[:, k] = (neutral - LOSS_MASS[loss] + z * PROTON) / z
    return out


def annotate(seq: str, mods, charge: int, mz: np.ndarray, intensity: np.ndarray,
             tol: float, ppm: bool):
    """Observed intensities per channel, normalised to the most intense valid ion.

    ``mz`` must be sorted. The observed value of an ion is the most intense peak
    within the tolerance (as in ``rescore_pin.rs``). Returns
    ``(values (n-1, N_ION) float32 with -1 for invalid channels, n_matched_base,
    explained)``; ``explained`` is the intensity of distinct matched peaks over the
    total spectrum intensity. ``values`` is None when nothing matched.
    """
    theo = fragment_mz(seq, mods)
    valid = valid_mask(seq, mods, charge)
    w = theo * tol * 1e-6 if ppm else np.full_like(theo, tol)
    lo = np.searchsorted(mz, theo - w, side="left")
    hi = np.searchsorted(mz, theo + w, side="right")
    obs = np.zeros(theo.shape, np.float64)
    used = np.zeros(len(mz), bool)
    for idx in zip(*np.nonzero(valid & (hi > lo))):
        a, b = lo[idx], hi[idx]
        j = a + int(np.argmax(intensity[a:b]))
        obs[idx] = intensity[j]
        used[j] = True
    top = obs[valid].max() if valid.any() else 0.0
    n_base = int(np.count_nonzero(obs[:, :N_BASE_ION] > 0))
    total = float(intensity.sum())
    explained = float(intensity[used].sum() / total) if total > 0 else 0.0
    if top <= 0:
        return None, n_base, explained
    vals = np.where(valid, obs / top, -1.0).astype(np.float32)
    return vals, n_base, explained


def parse_filter(filt: str | None):
    """Thermo filter string -> (activation name, NCE or nan).

    ``FTMS + p NSI d Full ms2 500.00@hcd28.00 [...]`` -> ("HCD-FTMS", 28.0).
    EThcD / ETD / multi-stage activations return ("unknown", nan)."""
    if not filt:
        return "unknown", float("nan")
    f = filt.strip()
    analyser = f.split(" ", 1)[0].upper()
    acts = re.findall(r"@([a-z]+)(\d+(?:\.\d+)?)", f)
    if len(acts) != 1 or acts[0][0] not in ("hcd", "cid"):
        return "unknown", float("nan")
    name = f"{acts[0][0].upper()}-{analyser}"
    return (name if name in ACTIVATIONS else "unknown"), float(acts[0][1])


# Shard schema (one row per annotated PSM).
SHARD_COLUMNS = [
    "dataset", "run", "scan", "peptidoform", "sequence", "charge", "nce", "activation",
    "instrument", "pep", "n_matched", "explained", "intensities",
]


def shard_table(rows: list[dict]):
    import pyarrow as pa

    return pa.table({
        "dataset": pa.array([r["dataset"] for r in rows], pa.string()),
        "run": pa.array([r["run"] for r in rows], pa.string()),
        "scan": pa.array([r["scan"] for r in rows], pa.int64()),
        "peptidoform": pa.array([r["peptidoform"] for r in rows], pa.string()),
        "sequence": pa.array([r["sequence"] for r in rows], pa.string()),
        "charge": pa.array([r["charge"] for r in rows], pa.int8()),
        "nce": pa.array([r["nce"] for r in rows], pa.float32()),
        "activation": pa.array([r["activation"] for r in rows], pa.string()),
        "instrument": pa.array([r["instrument"] for r in rows], pa.string()),
        "pep": pa.array([r["pep"] for r in rows], pa.float32()),
        "n_matched": pa.array([r["n_matched"] for r in rows], pa.int16()),
        "explained": pa.array([r["explained"] for r in rows], pa.float32()),
        "intensities": pa.array([r["intensities"].ravel() for r in rows], pa.list_(pa.float32())),
    })


if __name__ == "__main__":
    assert parse_proforma("[UNIMOD:737]-PEPS[UNIMOD:21]K[UNIMOD:737]") == ("PEPSK", [(0, 737), (4, 21), (5, 737)])
    assert parse_proforma("AC[UNIMOD:4]K-[UNIMOD:2]") == ("ACK", [(2, 4), (4, 2)])
    assert parse_proforma("PEP[Oxidation]K") is None
    for s in ["[UNIMOD:737]-PEPS[UNIMOD:21]K[UNIMOD:737]", "AC[UNIMOD:4]K-[UNIMOD:2]", "PEPTIDE"]:
        assert to_proforma(*parse_proforma(s)) == s
    # PEPTIDE: b2 = P+E+proton, y1 = E + H2O + proton
    t = fragment_mz("PEPTIDE", [])
    assert abs(t[1, 0] - (97.05276385 + 129.0425931 + PROTON)) < 1e-6
    assert abs(t[-1, 2] - (129.0425931 + H2O + PROTON)) < 1e-6
    assert abs(precursor_mz("PEPTIDE", [], 2) - 400.687258) < 1e-4
    v = valid_mask("MPSK", [(1, 35), (3, 21)], 1)
    assert v[:, 1].sum() == 0 and v[:, 3].sum() == 0
    assert list(v[:, 4]) == [False, False, True] and list(v[:, 6]) == [True, True, False]
    assert list(v[:, 8]) == [True, True, True] and list(v[:, 10]) == [False, False, False]
    # synthetic spectrum: every b1/y1 of PEPTIDEK at intensity 1, y1 at 4
    th = fragment_mz("PEPTIDEK", [])
    mz = np.sort(np.r_[th[:, 0], th[:, 2]])
    inten = np.ones_like(mz)
    inten[np.argmin(np.abs(mz - th[-1, 2]))] = 4.0
    vals, nb, ex = annotate("PEPTIDEK", [], 2, mz, inten, 20, True)
    assert nb == 14 and abs(ex - 1.0) < 1e-9 and vals[-1, 2] == 1.0 and abs(vals[0, 0] - 0.25) < 1e-6
    assert (vals[:, 4:] == -1).all() and (vals[:, 1] == 0).all()
    assert parse_filter("FTMS + p NSI d Full ms2 500.00@hcd28.00 [100.00-1000.00]") == ("HCD-FTMS", 28.0)
    assert parse_filter("ITMS + c NSI r d Full ms2 445.12@cid35.00 [110.00-905.00]") == ("CID-ITMS", 35.0)
    assert parse_filter("FTMS + p NSI d Full ms2 445.12@etd50.00@hcd20.00 [110.00-905.00]")[0] == "unknown"
    assert instrument_class("Q Exactive HF-X") == 1 and instrument_class("Orbitrap Fusion Lumos") == 3
    assert instrument_class("Orbitrap Astral") == 5 and instrument_class("LTQ Orbitrap Velos") == 4
    print("ms2_common ok")
