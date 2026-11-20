"""Shared definitions for andes-ms2: ion channels, fragment masses, b/y annotation,
acquisition metadata and the parquet shard schema.

Channel layout per cleavage site (``N_ION`` = 12), site ``i`` = cleavage after
residue ``i + 1`` (0-based), b = N-terminal fragment of length ``i + 1``,
y = the complementary C-terminal fragment::

    0 b1   1 b2   2 y1   3 y2                      plain b/y, fragment charge 1/2
    4..7   the same four, minus H3PO4 (97.9769)    only if the fragment holds pS/pT
    8..11  the same four, minus CH4SO (63.9983)    only if the fragment holds oxidised M

The first four channels are the PeptDeep ``b_z1, b_z2, y_z1, y_z2`` layout used by
``examples/rescore_pin.rs``.

The v2 layout (``ION_TYPES_V2``, ``N_ION_V2`` = 19) keeps those 12 channels first, in
the same order, and appends::

    12 b3  13 y3                                   plain b/y, fragment charge 3
    14 a1                                          b - CO
    15 c1  16 c2                                   b + NH3
    17 z1  18 z2                                   z-dot: y - NH3 + H

c and z-dot channels are valid only for electron-based activation (ETD, EThcD) and
never at a cleavage N-terminal to proline (the proline ring keeps the fragments
together). A channel is *invalid* (stored as -1) when its
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
from mod_composition import UNIMOD  # noqa: E402

# Element and group masses: the same f64 expressions as ``crates/model/src/mass.rs`` (bit-identical;
# ``test_v2_masses_equal_model_mass_rs`` parses that file).
EL_H = 1.00782503223
EL_C = 12.0
EL_N = 14.00307400443
EL_O = 15.99491461957
EL_S = 31.9720711744
EL_P = 30.97376199842
EL_ISOTOPE = {"2H": 2.01410177812, "13C": 13.00335483507, "15N": 15.00010889888, "18O": 17.99915961286}
PROTON = 1.0072764665789
H2O = EL_H * 2.0 + EL_O
NH3 = EL_N + EL_H * 3.0
CO = EL_C + EL_O
H3PO4 = EL_H * 3.0 + EL_P + EL_O * 4.0
CH4SO = EL_C + EL_H * 4.0 + EL_S + EL_O
H_ATOM = EL_H
Z_DOT_OFFSET = -(EL_N + EL_H * 2.0)
LOSS_MASS = {None: 0.0, "H3PO4": H3PO4, "CH4SO": CH4SO}

# (ion, fragment charge, neutral loss)
ION_TYPES = [
    ("b", 1, None), ("b", 2, None), ("y", 1, None), ("y", 2, None),
    ("b", 1, "H3PO4"), ("b", 2, "H3PO4"), ("y", 1, "H3PO4"), ("y", 2, "H3PO4"),
    ("b", 1, "CH4SO"), ("b", 2, "CH4SO"), ("y", 1, "CH4SO"), ("y", 2, "CH4SO"),
]
N_ION = len(ION_TYPES)
N_BASE_ION = 4
ION_TYPES_V2 = ION_TYPES + [
    ("b", 3, None), ("y", 3, None),
    ("a", 1, None),
    ("c", 1, None), ("c", 2, None),
    ("z", 1, None), ("z", 2, None),
]
N_ION_V2 = len(ION_TYPES_V2)
PREFIX_SERIES = ("a", "b", "c")
ELECTRON_SERIES = ("c", "z")

# Activation = dissociation + MS2 analyser; index 0 = unknown. v2 entries only after the v1 ones.
ACTIVATIONS = ["unknown", "HCD-FTMS", "HCD-ITMS", "CID-ITMS", "CID-FTMS", "HCD-ASTMS", "ETD", "EThcD", "HCD-TOF"]
ELECTRON_ACTIVATIONS = ("ETD", "EThcD")
# Instrument family; index 0 = unknown. v2 entries only after the v1 ones.
INSTRUMENTS = ["unknown", "QExactive", "Exploris", "Tribrid", "VelosElite", "Astral", "LTQ", "timsTOF", "TOF"]
# Rows of the activation / instrument embeddings in a v1 (AMS20001) model.
V1_N_ACTIVATIONS = 6
V1_N_INSTRUMENTS = 7
# Isobaric label condition (v2 only); index 0 = none.
LABELS = ["none", "TMT", "iTRAQ"]
TMT_MODS = frozenset({737, 738, 739, 2016})
ITRAQ_MODS = frozenset({214, 730})
CARBAMIDOMETHYL = 4
# Precursor charges 1..MAX_CHARGE are one-hot encoded; higher charges use the last slot.
MAX_CHARGE = 6

# (C, H, N, O, S) per residue, as ``model::amino_acid::standard_composition``.
RESIDUE_COMPOSITION = {
    "G": (2, 3, 1, 1, 0), "A": (3, 5, 1, 1, 0), "S": (3, 5, 1, 2, 0), "P": (5, 7, 1, 1, 0),
    "V": (5, 9, 1, 1, 0), "T": (4, 7, 1, 2, 0), "C": (3, 5, 1, 1, 1), "L": (6, 11, 1, 1, 0),
    "I": (6, 11, 1, 1, 0), "N": (4, 6, 2, 2, 0), "D": (4, 5, 1, 3, 0), "Q": (5, 8, 2, 2, 0),
    "K": (6, 12, 2, 1, 0), "E": (5, 7, 1, 3, 0), "M": (5, 9, 1, 1, 1), "H": (6, 7, 3, 1, 0),
    "F": (9, 9, 1, 1, 0), "R": (6, 12, 4, 1, 0), "Y": (9, 9, 1, 2, 0), "W": (11, 10, 2, 1, 0),
}
RESIDUE_MASS = {
    aa: float(c) * EL_C + float(h) * EL_H + float(n) * EL_N + float(o) * EL_O + float(sx) * EL_S
    for aa, (c, h, n, o, sx) in RESIDUE_COMPOSITION.items()
}


def _composition_mass(comp: dict) -> float:
    """Monoisotopic mass of a Unimod composition with the ``mass.rs`` element masses (others from
    ``mod_composition.MONO``)."""
    from mod_composition import MONO

    el = {"H": EL_H, "C": EL_C, "N": EL_N, "O": EL_O, "S": EL_S, "P": EL_P, **EL_ISOTOPE}
    return sum(el.get(e, MONO[e]) * n for e, n in comp.items())


def _mod_masses() -> dict:
    from mod_composition import parse_composition

    return {k: _composition_mass(parse_composition(c)) for k, (_, c, _) in UNIMOD.items()}


MOD_MASS = _mod_masses()

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
    if "timstof" in s:
        return INSTRUMENTS.index("timsTOF")
    if "tof" in s:
        return INSTRUMENTS.index("TOF")
    return 0


def activation_index(name: str | None) -> int:
    return ACTIVATIONS.index(name) if name in ACTIVATIONS else 0


def label_index(name: str | None) -> int:
    """Index in ``LABELS`` for a label name (``none``, ``TMT``, ``TMTpro``, ``iTRAQ8plex`` ...)."""
    s = (name or "").lower()
    if "itraq" in s:
        return LABELS.index("iTRAQ")
    if "tmt" in s:
        return LABELS.index("TMT")
    return 0


def ion_types(version: int):
    """Channel layout of an AMS2 model of the given format version."""
    return ION_TYPES_V2 if version >= 2 else ION_TYPES


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
        m[min(max(pos, 1), n) - 1] += MOD_MASS[uid]
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


def valid_mask(seq: str, mods, charge: int, ions=ION_TYPES, activation: int = 0) -> np.ndarray:
    """(n - 1, len(ions)) bool: which channels exist for this precursor.

    ``activation`` is an ``ACTIVATIONS`` index; it only matters for c / z-dot channels."""
    n = len(seq)
    ph, ox = loss_carriers(seq, mods)
    ph_pre, ox_pre = np.cumsum(ph), np.cumsum(ox)
    electron = 0 <= activation < len(ACTIVATIONS) and ACTIVATIONS[activation] in ELECTRON_ACTIVATIONS
    v = np.zeros((n - 1, len(ions)), bool)
    for i in range(n - 1):
        has_b = {"H3PO4": ph_pre[i] > 0, "CH4SO": ox_pre[i] > 0}
        has_y = {"H3PO4": ph_pre[-1] - ph_pre[i] > 0, "CH4SO": ox_pre[-1] - ox_pre[i] > 0}
        before_proline = seq[i + 1] == "P"
        for k, (ion, z, loss) in enumerate(ions):
            if z > charge:
                continue
            if ion in ELECTRON_SERIES and (not electron or before_proline):
                continue
            if loss is None:
                v[i, k] = True
            else:
                v[i, k] = has_b[loss] if ion in PREFIX_SERIES else has_y[loss]
    return v


def fragment_mz(seq: str, mods, ions=ION_TYPES) -> np.ndarray:
    """(n - 1, len(ions)) theoretical m/z for every channel (valid or not)."""
    res = residue_masses(seq, mods)
    total = res.sum()
    prefix = np.cumsum(res)[:-1]
    suffix = total - prefix + H2O
    neutral_of = {"b": prefix, "a": prefix - CO, "c": prefix + NH3, "y": suffix, "z": suffix + Z_DOT_OFFSET}
    out = np.zeros((len(prefix), len(ions)))
    for k, (ion, z, loss) in enumerate(ions):
        out[:, k] = (neutral_of[ion] - LOSS_MASS[loss] + z * PROTON) / z
    return out


def annotate(seq: str, mods, charge: int, mz: np.ndarray, intensity: np.ndarray,
             tol: float, ppm: bool, ions=ION_TYPES, activation: int = 0):
    """Observed intensities per channel, normalised to the most intense valid ion.

    ``mz`` must be sorted. The observed value of an ion is the most intense peak
    within the tolerance (as in ``rescore_pin.rs``). Returns
    ``(values (n-1, len(ions)) float32 with -1 for invalid channels, n_matched_base,
    explained)``; ``explained`` is the intensity of distinct matched peaks over the
    total spectrum intensity. ``values`` is None when nothing matched.
    """
    theo = fragment_mz(seq, mods, ions)
    valid = valid_mask(seq, mods, charge, ions, activation)
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


def parse_filter_v2(filt: str | None):
    """Thermo filter string -> (v2 activation name, NCE or nan).

    As ``parse_filter`` for one HCD/CID step; ``@etd..@hcd..`` (or ``@cid``) is ``EThcD`` with the
    supplemental collision energy as NCE, ``@etd..`` alone is ``ETD`` with NCE nan."""
    if not filt:
        return "unknown", float("nan")
    acts = re.findall(r"@([a-z]+)(\d+(?:\.\d+)?)", filt.strip())
    kinds = [k for k, _ in acts]
    if kinds and kinds[0] == "etd":
        if len(acts) == 1:
            return "ETD", float("nan")
        if len(acts) == 2 and kinds[1] in ("hcd", "cid"):
            return "EThcD", float(acts[1][1])
        return "unknown", float("nan")
    return parse_filter(filt)


def v2_activation(dissociation, analyser, instrument_cls: str | None) -> str:
    """v2 activation name from a dissociation method, an MS2 analyser and an ``INSTRUMENTS`` name.

    Collisional activation on a time-of-flight analyser is HCD-TOF; an Orbitrap Astral MS2 is HCD-ASTMS;
    Q Exactive and Exploris have only HCD with Orbitrap detection (HCD-FTMS) and timsTOF only beam-type
    collisional activation (HCD-TOF), whatever the label says. 'unknown' when it cannot be decided."""
    d = str(dissociation or "").upper().replace("-", "").replace(" ", "")
    an = str(analyser or "").upper().replace("-", "")
    if "ETHCD" in d or "ETCID" in d:
        return "EThcD"
    if "ETD" in d or "ELECTRONTRANSFER" in d:
        return "ETD"
    tof = an in ("TOF", "QTOF", "TIMEOFFLIGHT") or instrument_cls in ("timsTOF", "TOF")
    collisional = d in ("HCD", "CID", "CAD", "BEAMTYPECID", "COLLISIONINDUCEDDISSOCIATION",
                        "HIGHERENERGYBEAMTYPECOLLISIONINDUCEDDISSOCIATION", "BEAMTYPECOLLISIONINDUCEDDISSOCIATION")
    if instrument_cls in ("QExactive", "Exploris") and (not d or collisional):
        return "HCD-FTMS"
    if instrument_cls == "timsTOF" and (not d or collisional):
        return "HCD-TOF"
    if not collisional:
        return "unknown"
    if tof:
        return "HCD-TOF"
    hcd = d in ("HCD", "BEAMTYPECID", "HIGHERENERGYBEAMTYPECOLLISIONINDUCEDDISSOCIATION",
                "BEAMTYPECOLLISIONINDUCEDDISSOCIATION")
    if an == "ASTMS" or (instrument_cls == "Astral" and an not in ("FTMS", "ITMS")):
        return "HCD-ASTMS"
    if an in ("FTMS", "ITMS"):
        return f"{'HCD' if hcd else 'CID'}-{an}"
    return "unknown"


def label_of_mods(mods) -> str:
    """``LABELS`` name of the isobaric tag a peptidoform carries ("none" without one)."""
    ids = {uid for _, uid in mods}
    if ids & TMT_MODS:
        return "TMT"
    if ids & ITRAQ_MODS:
        return "iTRAQ"
    return "none"


SHARD_SCHEMA_VERSION = "2.0"
# Shard schema v2 (one row per annotated PSM); ``intensities`` holds ``N_ION_V2`` channels per site.
SHARD_COLUMNS = [
    "dataset", "run", "scan", "peptidoform", "sequence", "charge", "nce", "activation",
    "instrument", "pep", "n_matched", "explained", "intensities", "label", "hard", "agree",
]


def shard_table(rows: list[dict]):
    """Rows -> a v2 shard table with the key metadata ``andes.ms2_shard.schema`` = ``2.0``."""
    import pyarrow as pa

    t = pa.table({
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
        "label": pa.array([r["label"] for r in rows], pa.string()),
        "hard": pa.array([bool(r["hard"]) for r in rows], pa.bool_()),
        "agree": pa.array([r["agree"] for r in rows], pa.bool_()),
    })
    return t.replace_schema_metadata({"andes.ms2_shard.schema": SHARD_SCHEMA_VERSION})


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
    assert instrument_class("timsTOF Pro 2") == 7 and instrument_class("TripleTOF 6600") == 8
    assert activation_index("EThcD") == 7 and label_index("TMTpro") == 1 and label_index("iTRAQ8plex") == 2
    assert ION_TYPES_V2[:N_ION] == ION_TYPES and N_ION_V2 == 19
    # PEPTIDE, hand-computed: c1 = P + NH3 + proton, a2 = P + E - CO + proton, z-dot1 = E + H2O - NH3 + H + proton
    t2 = fragment_mz("PEPTIDE", [], ION_TYPES_V2)
    assert abs(t2[0, 15] - (97.05276385 + 17.02654910 + PROTON)) < 1e-6
    assert abs(t2[1, 14] - (97.05276385 + 129.0425931 - 27.99491462 + PROTON)) < 1e-6
    assert abs(t2[-1, 17] - (129.0425931 + H2O - 17.02654910 + 1.00782503 + PROTON)) < 1e-6
    assert abs(t2[0, 12] * 3 - (97.05276385 + 3 * PROTON)) < 1e-6
    v2 = valid_mask("PEPTIDE", [], 3, ION_TYPES_V2, activation_index("ETD"))
    assert v2[:, 12].all() and list(v2[:, 15]) == [True, False, True, True, True, True]
    assert not valid_mask("PEPTIDE", [], 3, ION_TYPES_V2, activation_index("HCD-FTMS"))[:, 15:].any()
    assert not valid_mask("PEPTIDE", [], 2, ION_TYPES_V2)[:, 12:14].any()
    print("ms2_common ok")
