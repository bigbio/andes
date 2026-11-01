"""FASTA: species and protein of each precursor.

A precursor maps to a protein only when its peptide (I = L) is contained in exactly one target entry; its species is set when
every entry containing it is of one species. Species come from the entry name: `_HUMAN`, `_YEAST`, `_ECOLI` suffixes, and
UPS for names containing `UPS` that do not end in `_YEAST`. `CONT_` entries are skipped."""
import re

import pandas as pd

from .keys import bare


def read_fasta(path):
    out, name, cur = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(cur)))
                name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line.strip())
    if name:
        out.append((name, "".join(cur)))
    return out


def species_of_entry(name, classes):
    u = name.upper()
    if "CONT_" in u:
        return None
    if "UPS" in classes and "UPS" in u and not u.endswith("_YEAST"):
        return "UPS"
    m = re.search(r"_(HUMAN|YEAST|ECOLI)$", u)
    return m.group(1) if (m and m.group(1) in classes) else None


def map_precursors(keys, fasta, classes, k=6):
    """-> (species per key, protein per key, species per entry name). Species: a class name, MIXED, or NONE."""
    ent = [(n, s.replace("I", "L"), species_of_entry(n, classes)) for n, s in read_fasta(fasta)]
    ent = [e for e in ent if e[2] is not None]
    idx = {}
    for j, (_, s, _) in enumerate(ent):
        for i in range(len(s) - k + 1):
            idx.setdefault(s[i:i + k], set()).add(j)
    sp, prot, cache = {}, {}, {}
    for key in keys:
        b = bare(key).replace("I", "L")
        if b not in cache:
            hits = [j for j in idx.get(b[:k], ()) if b in ent[j][1]]
            sps = {ent[j][2] for j in hits}
            cache[b] = ((sps.pop() if len(sps) == 1 else ("MIXED" if sps else "NONE")), (ent[hits[0]][0] if len(hits) == 1 else None))
        sp[key], prot[key] = cache[b]
    return pd.Series(sp, dtype=object), pd.Series(prot, dtype=object), {n: s for n, _, s in ent}
