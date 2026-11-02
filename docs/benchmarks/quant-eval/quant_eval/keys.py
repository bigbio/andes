"""Precursor keys: every tool's modified sequence + charge is mapped to the andes key form.

Key = residues in order, each modified residue followed by `[+m.mmmmm]` (the summed mass delta of its modifications; a
peptide or protein N-terminal modification is folded onto residue 1), then `/charge`. Only Carbamidomethyl, Oxidation and
Acetyl are mapped; any other modification raises ValueError and the reader drops (and counts) that precursor."""
import re

MASS = {"Carbamidomethyl": 57.02146, "Oxidation": 15.99491, "Acetyl": 42.010565}
UNIMOD = {"UNIMOD:4": "Carbamidomethyl", "UNIMOD:35": "Oxidation", "UNIMOD:1": "Acetyl"}
_ORDER = ("Carbamidomethyl", "Oxidation", "Acetyl")
_COMBOS = [(sum(MASS[n] for n in c), c) for c in
           [(a,) for a in _ORDER] + [(a, b) for i, a in enumerate(_ORDER) for b in _ORDER[i + 1:]] + [_ORDER]]


def names_of_delta(m):
    """Numeric residue tag -> the modification names it is made of (tolerance 0.002 Da)."""
    hits = [c for s, c in _COMBOS if abs(s - m) < 0.002]
    if len(hits) != 1:
        raise ValueError(f"cannot read modification delta {m}")
    return list(hits[0])


def fmt_key(seq, deltas, charge):
    out = []
    for i, aa in enumerate(seq):
        out.append(aa)
        if i in deltas:
            names = sorted(deltas[i], key=_ORDER.index)
            out.append("[{:+.5f}]".format(sum(MASS[n] for n in names)))
    return "".join(out) + "/" + str(int(charge))


def andes_key(key):
    """Normalise an andes precursor key: UNIMOD accession tags or numeric tags -> one numeric form."""
    pep, z = key.rsplit("/", 1)
    seq, deltas, i = [], {}, -1
    for tok in re.findall(r"[A-Z]|\[[^\]]*\]", pep):
        if tok.startswith("["):
            t = tok[1:-1]
            names = [UNIMOD[t]] if t in UNIMOD else names_of_delta(float(t))
            deltas.setdefault(max(i, 0), []).extend(names)
        else:
            seq.append(tok)
            i += 1
    return fmt_key("".join(seq), deltas, z)


def norm_key(k):
    try:
        return andes_key(k)
    except Exception:
        return k


def _mod_name_from_mass(m):
    for name, ref in MASS.items():
        if abs(m - ref) < 0.01:
            return name
    raise ValueError(f"unmapped modification mass {m}")


# FragPipe 'Assigned Modifications': "N-term(42.0106), 5M(15.9949), 7C(57.0215)" (1-based positions)
_FP_MOD = re.compile(r"^\s*(N-term|n|C-term|c|(\d+)([A-Z]))\((-?[\d.]+)\)\s*$")


def fragpipe_key(seq, assigned, charge):
    deltas = {}
    for part in [p for p in (assigned or "").split(",") if p.strip()]:
        m = _FP_MOD.match(part)
        if not m:
            raise ValueError(f"cannot parse FragPipe modification {part!r} of {seq}")
        site, pos, aa, mass = m.group(1), m.group(2), m.group(3), float(m.group(4))
        name = _mod_name_from_mass(mass)
        if site in ("N-term", "n"):
            idx = 0
        elif site in ("C-term", "c"):
            idx = len(seq) - 1
        else:
            idx = int(pos) - 1
            if seq[idx] != aa:
                raise ValueError(f"FragPipe modification {part!r} does not match {seq}")
        deltas.setdefault(idx, []).append(name)
    return fmt_key(seq, deltas, charge)


# MaxQuant 'Modified sequence': "_(Acetyl (Protein N-term))M(Oxidation (M))PEPTIDEK_" (variable modifications only; the fixed
# Carbamidomethyl (C) is not printed and is added on every C). Older notations "(ox)", "(ac)" are accepted.
_MQ_NAMES = {"Oxidation (M)": "Oxidation", "ox": "Oxidation", "Acetyl (Protein N-term)": "Acetyl", "ac": "Acetyl",
             "Acetyl (N-term)": "Acetyl", "Carbamidomethyl (C)": "Carbamidomethyl", "cam": "Carbamidomethyl"}


def maxquant_key(modseq, charge, fixed_cam=True):
    s = modseq.strip("_")
    if not s:
        raise ValueError(f"empty MaxQuant modified sequence (charge {charge})")
    seq, deltas, i, pos = [], {}, -1, 0
    while pos < len(s):
        c = s[pos]
        if c == "(":
            depth, j = 0, pos
            while True:  # balanced parentheses: names contain "(M)"
                if s[j] == "(":
                    depth += 1
                elif s[j] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            name = s[pos + 1:j]
            if name not in _MQ_NAMES:
                raise ValueError(f"unmapped MaxQuant modification {name!r} in {modseq}")
            deltas.setdefault(max(i, 0), []).append(_MQ_NAMES[name])
            pos = j + 1
        elif c == "_":
            pos += 1
        else:
            seq.append(c)
            i += 1
            pos += 1
    seq = "".join(seq)
    if fixed_cam:
        for k, aa in enumerate(seq):
            if aa == "C" and "Carbamidomethyl" not in deltas.get(k, []):
                deltas.setdefault(k, []).append("Carbamidomethyl")
    return fmt_key(seq, deltas, charge)


# FlashLFQ QuantifiedPeaks 'Full Sequence' in mzLib (MetaMorpheus) notation: each modification in brackets after its residue,
# "[Type:Name on Target]"; an N-terminal modification stands before residue 1, a C-terminal one after "-" at the end, e.g.
# "[Common Biological:Acetylation on X]M[Common Variable:Oxidation on M]PEC[Common Fixed:Carbamidomethyl on C]K".
# Fixed modifications are printed (no implicit Carbamidomethyl). Numeric tags ("[+15.9949]"), UNIMOD accessions and bare names
# are accepted as well.
_MZLIB_NAMES = {"carbamidomethyl": "Carbamidomethyl", "carbamidomethylation": "Carbamidomethyl", "oxidation": "Oxidation",
                "acetyl": "Acetyl", "acetylation": "Acetyl"}


def _mzlib_mod(tag):
    """One bracket tag -> modification names. UniProt protein N-terminal acetylations (N-acetylserine, ...) count as Acetyl."""
    t = tag.strip()
    if t in UNIMOD:
        return [UNIMOD[t]]
    if t in MASS:
        return [t]
    if ":" in t:
        n = t.split(":", 1)[1].rsplit(" on ", 1)[0].strip().lower()
        if n in _MZLIB_NAMES:
            return [_MZLIB_NAMES[n]]
        if n.startswith("n-acetyl"):
            return ["Acetyl"]
        raise ValueError(f"unmapped mzLib modification {tag!r}")
    try:
        return names_of_delta(float(t))
    except ValueError:
        raise ValueError(f"unmapped mzLib modification {tag!r}")


def flashlfq_key(fullseq, charge):
    s = fullseq.strip()
    seq, deltas, i, pos = [], {}, -1, 0
    while pos < len(s):
        c = s[pos]
        if c == "[":
            depth, j = 0, pos
            while j < len(s):  # balanced brackets: UniProt names may contain "[...]"
                depth += {"[": 1, "]": -1}.get(s[j], 0)
                if depth == 0:
                    break
                j += 1
            if j == len(s):
                raise ValueError(f"unbalanced bracket in {fullseq!r}")
            deltas.setdefault(max(i, 0), []).extend(_mzlib_mod(s[pos + 1:j]))
            pos = j + 1
        elif c == "-":  # separator before a C-terminal modification
            pos += 1
        elif "A" <= c <= "Z":
            seq.append(c)
            i += 1
            pos += 1
        else:
            raise ValueError(f"cannot read FlashLFQ sequence {fullseq!r}")
    if not seq:
        raise ValueError(f"empty FlashLFQ sequence (charge {charge})")
    return fmt_key("".join(seq), deltas, charge)


def bare(k):
    """Bare sequence of an andes key (modification tags and charge removed)."""
    return re.sub(r"\[[^\]]*\]", "", k.rsplit("/", 1)[0]).replace("-", "")
