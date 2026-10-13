"""Built-in Unimod elemental compositions for the andes-rt / andes-ms2 encoders.

No database lookup: a silently missing unimod sqlite is how other predictors end
up encoding unknown mods as zeros. Each entry is (name, composition delta,
allowed sites); heavy isotopes are folded into their element for the feature
vector, while the monoisotopic mass uses the isotope masses.

The per-position feature vector is ``[C, H, N, O, S, P, other, mass / 100]``.
The Rust module ``rt_native::mods`` carries the same table; the parity test
(``rt_export.py`` reference file) fails if the two disagree.
"""
from __future__ import annotations

import re

MONO = {
    "C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956,
    "S": 31.97207100, "P": 30.97376163, "Na": 22.9897692809, "K": 38.9637066,
    "Se": 79.9165213, "13C": 13.0033548378, "15N": 15.0001088982, "18O": 17.9991610,
    "2H": 2.0141017778,
}
ELEMENTS = ["C", "H", "N", "O", "S", "P"]
FOLD = {"13C": "C", "15N": "N", "18O": "O", "2H": "H"}
N_COMP = 8  # C H N O S P other mass/100

# id: (name, composition, sites); "^" = peptide N-term, "$" = C-term
UNIMOD = {
    1: ("Acetyl", "H2 C2 O", "^K"),
    2: ("Amidated", "H N O-1", "$"),
    4: ("Carbamidomethyl", "H3 C2 N O", "C"),
    5: ("Carbamyl", "H C N O", "^K"),
    7: ("Deamidated", "H-1 N-1 O", "NQR"),
    21: ("Phospho", "H O3 P", "STYH"),
    23: ("Dehydrated", "H-2 O-1", "STDY"),
    24: ("Propionamide", "H5 C3 N O", "C"),
    26: ("Pyro-carbamidomethyl", "C2 O", "^"),
    27: ("Glu->pyro-Glu", "H-2 O-1", "^"),
    28: ("Gln->pyro-Glu", "H-3 N-1", "^"),
    30: ("Cation:Na", "H-1 Na", "DE$"),
    34: ("Methyl", "H2 C", "KRHDE^"),
    35: ("Oxidation", "O", "MWHC"),
    36: ("Dimethyl", "H4 C2", "KR^"),
    37: ("Trimethyl", "H6 C3", "K"),
    40: ("Sulfo", "O3 S", "STY"),
    58: ("Propionyl", "H4 C3 O", "K^"),
    64: ("Succinyl", "H4 C4 O3", "K^"),
    121: ("GlyGly", "H6 C4 N2 O2", "K"),
    122: ("Formyl", "C O", "K^"),
    214: ("iTRAQ4plex", "H12 C4 13C3 N 15N O", "K^Y"),
    259: ("Label:13C(6)15N(2)", "13C6 C-6 15N2 N-2", "K"),
    267: ("Label:13C(6)15N(4)", "13C6 C-6 15N4 N-4", "R"),
    312: ("Cysteinyl", "H5 C3 N O2 S", "C"),
    354: ("Nitro", "H-1 N O2", "YW"),
    385: ("Ammonia-loss", "H-3 N-1", "^"),
    425: ("Dioxidation", "O2", "MWC"),
    730: ("iTRAQ8plex", "H24 C7 13C7 N3 15N O3", "K^Y"),
    737: ("TMT6plex", "H20 C8 13C4 N 15N O2", "K^"),
    738: ("TMT2plex", "H20 C11 13C N2 O2", "K^"),
    739: ("TMT", "H20 C12 N2 O2", "K^"),
    747: ("Malonyl", "H2 C3 O3", "K"),
    1363: ("Crotonyl", "H4 C4 O", "K"),
    2016: ("TMTpro", "H25 C8 13C7 N 15N2 O3", "K^"),
}

_TERM = re.compile(r"^(\d*[A-Z][a-z]?)(-?\d+)?$")


def parse_composition(s: str) -> dict:
    out: dict = {}
    for tok in s.split():
        m = re.match(r"^(\d+[A-Z][a-z]?|[A-Z][a-z]?)(-?\d+)?$", tok)
        if not m:
            raise ValueError(f"bad composition token {tok!r}")
        el, n = m.group(1), int(m.group(2) or 1)
        out[el] = out.get(el, 0) + n
    return out


def mass_of(comp: dict) -> float:
    return sum(MONO[e] * n for e, n in comp.items())


def comp_vector(comp: dict) -> list:
    v = [0.0] * N_COMP
    for e, n in comp.items():
        e2 = FOLD.get(e, e)
        if e2 in ELEMENTS:
            v[ELEMENTS.index(e2)] += n
        else:
            v[6] += n
    v[7] = mass_of(comp) / 100.0
    return v


VECTORS = {k: comp_vector(parse_composition(c)) for k, (_, c, _) in UNIMOD.items()}
MASSES = {k: mass_of(parse_composition(c)) for k, (_, c, _) in UNIMOD.items()}


def match_mass(delta: float, residue: str | None, nterm: bool, tol: float = 2e-3, cterm: bool = False):
    """Unimod id for a mass delta at a residue (or the N- or C-terminus), else None.

    Among site-compatible entries within ``tol``, the smallest mass error wins
    and the lowest id breaks exact ties (e.g. Oxidation over other O-adding
    entries). TMTpro and iTRAQ8plex are 1.8 mDa apart, so the error, not the
    id, has to separate them. Mirrors ``rt_native::mods::match_mass``."""
    ok = [k for k, m in MASSES.items()
          if abs(m - delta) <= tol
          and ((nterm and "^" in UNIMOD[k][2]) or (cterm and "$" in UNIMOD[k][2])
               or (residue and residue in UNIMOD[k][2]))]
    return min(ok, key=lambda k: (abs(MASSES[k] - delta), k)) if ok else None


def split_nterm_residue(delta: float, residue: str, tol: float = 2e-3):
    """(n-term id, residue id) whose masses sum to a delta printed on residue 1.

    andes prints stacked N-terminal + residue mods as one delta. The pair
    closest in mass wins; exact ties go to the lowest N-terminal id, then the
    lowest residue id. Mirrors ``rt_native::mods::split_nterm_residue``."""
    best = None
    for n in sorted(k for k in UNIMOD if "^" in UNIMOD[k][2]):
        for r in sorted(k for k in UNIMOD if residue in UNIMOD[k][2]):
            err = abs(MASSES[n] + MASSES[r] - delta)
            if err <= tol and (best is None or err < best[0]):
                best = (err, n, r)
    return None if best is None else (best[1], best[2])


if __name__ == "__main__":
    ref = {1: 42.010565, 4: 57.021464, 7: 0.984016, 21: 79.966331, 27: -18.010565, 28: -17.026549,
           35: 15.994915, 121: 114.042927, 214: 144.102063, 354: 44.985078, 730: 304.205360,
           737: 229.162932, 2016: 304.207146, 259: 8.014199, 739: 224.152478, 30: 21.981943, 738: 225.155833}
    for k, m in ref.items():
        assert abs(MASSES[k] - m) < 2e-5, (k, MASSES[k], m)
    assert match_mass(79.96633, "S", False) == 21
    assert match_mass(229.1629, "K", False) == 737 and match_mass(229.1629, "P", True) == 737
    assert match_mass(42.01057, "M", True) == 1
    assert match_mass(-0.98402, None, False, cterm=True) == 2 and match_mass(-0.98402, None, False) is None
    assert match_mass(304.20715, "K", False) == 2016 and match_mass(304.20715, None, True) == 2016
    assert match_mass(304.20536, "K", False) == 730
    print("composition table ok:", len(UNIMOD), "entries")
