"""Residue-pair keys in the andes form and the protein positions of peptides.

andes writes a residue pair as `ACC:pos -- ACC:pos`: each side lists every accession and 1-based
protein position of the linked residue, joined by `;`, and the two sides are sorted
(crates/andes/src/bin/andes/xl/fdr.rs). Keys of every tool are brought to this form, with
UniProt-style identifiers (`sp|P02768|ALBU_HUMAN`) reduced to their accession, so that the
tables of different tools join on the same string.
"""

SEP = " -- "
MATCH_EXAMPLES = 5


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


def _item_pairs(key):
    a, b = key.split(SEP)
    return {tuple(sorted((x, y))) for x in a.split(";") for y in b.split(";")}


def comparison_ids(keys_by_arm):
    """One comparison id per residue-pair key, for choosing the residue pairs the arms share.

    Tools can list different protein positions for one peptide: andes lists the positions of the
    peptide's digestion occurrences, `side` every substring occurrence. Keys that mean the same
    residue pair can then differ as strings. Keys that share one pair of `ACC:pos` items, one from
    each side, are collected into a set, and the keys of a set get one id (its smallest key) when

    - no arm has more than one key in the set, and
    - one pair of items is shared by every key of the set. Overlap that only chains (A with B, B
      with C, but A and C share nothing) is not a match.

    Otherwise the set is ambiguous and each of its keys keeps itself as its id.

    The ids are for the `common` selection only. The keys themselves are never changed: an arm's
    own metrics and its truth read the positions and proteins the arm reported.

    Returns (id of every key, report for the result file).
    """
    parent = {}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    owner, pairs = {}, {}
    for key in sorted(set().union(*keys_by_arm.values())) if keys_by_arm else []:
        parent[key] = key
        pairs[key] = _item_pairs(key)
        for pair in sorted(pairs[key]):
            if pair in owner:
                parent[find(key)] = find(owner[pair])
            else:
                owner[pair] = key
    sets = {}
    for key in parent:
        sets.setdefault(find(key), []).append(key)
    ids, matched, ambiguous, examples = {key: key for key in parent}, 0, 0, []
    for members in sets.values():
        if len(members) == 1:
            continue
        members.sort()
        one_per_arm = all(sum(k in keys for k in members) <= 1 for keys in keys_by_arm.values())
        shared = set.intersection(*(pairs[k] for k in members))
        if not (one_per_arm and shared):
            ambiguous += 1
            continue
        matched += 1
        ids.update({k: members[0] for k in members})
        if len(examples) < MATCH_EXAMPLES:
            examples.append(members)
    return ids, {"matched_sets": matched, "ambiguous_sets": ambiguous, "examples": examples}
