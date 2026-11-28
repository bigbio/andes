"""Residue-pair keys in the andes form and the protein positions of peptides.

andes writes a residue pair as `ACC:pos -- ACC:pos`: each side lists every accession and 1-based
protein position of the linked residue, joined by `;`, and the two sides are sorted
(crates/andes/src/bin/andes/xl/fdr.rs). Keys of every tool are brought to this form, with
UniProt-style identifiers (`sp|P02768|ALBU_HUMAN`) reduced to their accession, so that the
tables of different tools join on the same string.
"""

SEP = " -- "


def acc_of(token):
    """Accession of a FASTA identifier: the middle field of `db|ACC|NAME`, else the token."""
    parts = token.strip().split("|")
    return parts[1] if len(parts) >= 3 else token.strip()


def read_fasta(path):
    """Accession -> sequence of a FASTA file."""
    proteins, acc, seq = {}, None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if acc is not None:
                    proteins[acc] = "".join(seq)
                acc, seq = acc_of(line[1:].split()[0]), []
            elif line:
                seq.append(line)
    if acc is not None:
        proteins[acc] = "".join(seq)
    return proteins


def side(proteins, peptide, site):
    """Residue-pair side of a peptide linked at 0-based `site`: every `ACC:pos` it maps to."""
    hits = set()
    for acc, seq in proteins.items():
        at = seq.find(peptide)
        while at >= 0:
            hits.add(f"{acc}:{at + site + 1}")
            at = seq.find(peptide, at + 1)
    return ";".join(sorted(hits))


def pair_key(a, b):
    """The andes residue-pair (or protein-pair) key of two sides."""
    return f"{a}{SEP}{b}" if a <= b else f"{b}{SEP}{a}"


def normalise_key(key):
    """A residue-pair key of any tool in the canonical andes form, accessions reduced."""
    sides = []
    for s in str(key).split(SEP):
        items = []
        for item in s.split(";"):
            acc, _, pos = item.strip().rpartition(":")
            items.append(f"{acc_of(acc)}:{pos}" if acc else item.strip())
        sides.append(";".join(sorted(set(items))))
    if len(sides) != 2:
        raise ValueError(f"not a residue-pair key: {key!r}")
    return pair_key(*sides)


def accessions(side_text):
    """The accessions of one key side (`ACC:pos;ACC:pos` or `ACC;ACC`)."""
    return {item.rpartition(":")[0] or item for item in side_text.split(";") if item}
