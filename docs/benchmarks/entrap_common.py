#!/usr/bin/env python3
"""Shared entrapment rules for the benchmark evaluators.

One definition of "entrapment hit", used by `reproduce/run.sh`, `reproduce/build_databases.sh`,
`glyco/eval_entrap.py` and `glyco/compare_engines.py`:

* An identification is entrapment only if EVERY protein it maps to is an `ENTRAP_` protein.
  A peptide shared with a sample protein is not evidence of an error.
* ...and only if its bare peptide, with I folded to L, is NOT a substring of any I/L-folded
  sample protein (non-`ENTRAP_`, non-decoy). Isoleucine and leucine are isobaric, so such a
  peptide explains the spectrum exactly as well as the sample sequence does. The rule is a
  plain substring test and does not depend on the enzyme or on the search's digest settings.
  Such PSMs are counted as target ("rescued").

The T/E factor (`te_ratio`) sizes the entrapment correction in peptide space and therefore
uses a digest that mirrors the andes defaults. Its entrapment half applies the same rule as
the counting: an entrapment peptide that is a substring of a folded sample protein is moved
to the target half, because a hit on it would be counted as target.

Command line (used by the shell scripts):

    entrap_common.py count <percolator.psms> [search.fasta] [q-threshold]
        prints "<accepted> <entrapment hits> <rescued>" for Percolator target PSMs
    entrap_common.py te <search.fasta>
        prints the target and entrapment peptide counts and T/E
"""
import re
import sys

ENTRAP_TAG = "ENTRAP_"
DECOY_HEADER_PREFIXES = (">XXX_", ">DECOY")


def fold_il(seq):
    """Fold I to L: the two are isobaric (same precursor and fragment masses)."""
    return seq.replace("I", "L")


# `pre.CORE.post`, optionally followed by ONE bracketed tag after the flank.
_FLANKED = re.compile(r"^(.)\.(.*)\.(.)(\[[^\[\]]*\])?$")
_BRACKETED = re.compile(r"\[[^\[\]]*\]")
_SIGNED_MASS = re.compile(r"[+-]\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")


def bare_peptide(cell):
    """Bare amino-acid sequence from an andes PIN / Percolator `Peptide` cell.

    Accepted format: `pre.SEQ.post`, where `pre`/`post` are single flanking residues (`-` at
    a protein terminus) and modifications follow their residue either as a signed mass
    (`M+15.99491`, `Q-17.02655`) or as a bracketed accession (`M[UNIMOD:35]`). A glyco PIN
    appends one tag after the C-terminal flank (`K.TAANGTR.K[HexNAc2Hex5Fuc0NeuAc0NeuGc0@N4]`);
    only that trailing tag is dropped before the flanks are removed, so a bracket inside the
    sequence never truncates it. A cell without flanks is treated as the sequence itself.
    """
    s = cell.strip()
    m = _FLANKED.match(s)
    core = m.group(2) if m else s
    core = _BRACKETED.sub("", core)
    core = _SIGNED_MASS.sub("", core)
    return "".join(c for c in core if "A" <= c <= "Z")


def is_decoy_header(hdr):
    return hdr.startswith(DECOY_HEADER_PREFIXES)


def is_entrap_accession(acc):
    """Substring, not prefix: the tag may sit inside a UniProt-style accession."""
    return ENTRAP_TAG in acc


def all_entrapment(accessions):
    """True only if there is at least one accession and every one is entrapment."""
    accs = [a for a in accessions if a]
    return bool(accs) and all(is_entrap_accession(a) for a in accs)


def read_fasta(path):
    """Yield (header line, sequence) pairs."""
    hdr, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(buf)
                hdr, buf = line.strip(), []
            else:
                buf.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(buf)


def protein_class(hdr):
    """'decoy', 'entrap' or 'sample'. Decoys are tested first: a decoy built from an
    entrapment protein keeps the tag inside its accession (`>XXX_sp|ENTRAP_...`)."""
    if is_decoy_header(hdr):
        return "decoy"
    if ENTRAP_TAG in hdr:
        return "entrap"
    return "sample"


class SampleSequences:
    """I/L-folded sample protein sequences, searchable by substring.

    The sequences are joined with a newline, which is not an amino-acid letter, so a match
    can never span two proteins. Loading is deferred to the first query: a database without
    entrapment hits never pays for it. Answers are cached per peptide.
    """

    def __init__(self, fasta):
        self._fasta = fasta
        self._text = None
        self._cache = {}
        self._kmers = None

    def _load(self):
        if self._text is None:
            self._text = "\n".join(
                fold_il(seq) for hdr, seq in read_fasta(self._fasta)
                if protein_class(hdr) == "sample"
            )
        return self._text

    def contains(self, peptide):
        """Is the I/L-folded peptide a substring of a folded sample protein?"""
        p = fold_il(peptide)
        hit = self._cache.get(p)
        if hit is None:
            hit = bool(p) and p in self._load()
            self._cache[p] = hit
        return hit

    def contains_many(self, peptides, k=6):
        """Subset of `peptides` (already folded, length >= k) found in the sample text.

        For large inputs (the T/E digest) a set of sample k-mers screens out most
        candidates before the exact substring test.
        """
        text = self._load()
        if self._kmers is None:
            self._kmers = {text[i:i + k] for i in range(len(text) - k + 1)}
        km = self._kmers
        found = set()
        for p in peptides:
            if len(p) < k:
                if p in text:
                    found.add(p)
                continue
            mid = (len(p) - k) // 2
            if p[:k] in km and p[-k:] in km and p[mid:mid + k] in km and p in text:
                found.add(p)
        return found


def classify(accessions, peptide_cell, sample):
    """'sample', 'entrap' or 'rescued' for one accepted target identification.

    `sample` is a SampleSequences, or None to disable the I/L rescue.
    """
    if not all_entrapment(accessions):
        return "sample"
    if sample is not None and sample.contains(bare_peptide(peptide_cell)):
        return "rescued"
    return "entrap"


def count_percolator_psms(psms_path, qmax=0.01, sample=None):
    """(accepted, entrapment hits, rescued) from a Percolator target PSM file.

    Proteins occupy the `proteinIds` column and every column after it.
    """
    n = n_ent = n_res = 0
    with open(psms_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_q = header.index("q-value")
        i_pep = header.index("peptide")
        i_prot = header.index("proteinIds")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_prot:
                continue
            try:
                if float(f[i_q]) > qmax:
                    continue
            except ValueError:
                continue
            n += 1
            c = classify(f[i_prot:], f[i_pep], sample)
            if c == "entrap":
                n_ent += 1
            elif c == "rescued":
                n_res += 1
    return n, n_ent, n_res


def tryptic_peptides(seq):
    """I/L-folded tryptic peptides at the andes defaults: <= 1 missed cleavage, 6..50
    residues, cleavage after every K/R (no proline rule)."""
    s = fold_il(seq)
    sites = sorted(set([0] + [i + 1 for i in range(len(s) - 1) if s[i] in "KR"] + [len(s)]))
    out = set()
    for i in range(len(sites) - 1):
        for j in range(i + 1, min(i + 3, len(sites))):
            p = s[sites[i]:sites[j]]
            if 6 <= len(p) <= 50:
                out.add(p)
    return out


def te_ratio(fasta):
    """(T, E, moved): target and entrapment-only peptide counts for the FDP correction.

    T: folded tryptic peptides of the sample proteins. E: folded tryptic peptides of the
    entrapment proteins that are not in T and are not a substring of a folded sample protein
    (`moved` counts the latter; they join T, as a hit on one counts as target).
    """
    t, e = set(), set()
    for hdr, seq in read_fasta(fasta):
        cls = protein_class(hdr)
        if cls == "sample":
            t.update(tryptic_peptides(seq))
        elif cls == "entrap":
            e.update(tryptic_peptides(seq))
    e -= t
    moved = SampleSequences(fasta).contains_many(e)
    e -= moved
    t |= moved
    return len(t), len(e), len(moved)


def _main(argv):
    if len(argv) >= 2 and argv[0] == "count":
        psms = argv[1]
        fasta = argv[2] if len(argv) > 2 and argv[2] else None
        qmax = float(argv[3]) if len(argv) > 3 else 0.01
        sample = SampleSequences(fasta) if fasta else None
        print(*count_percolator_psms(psms, qmax, sample))
        return 0
    if len(argv) == 2 and argv[0] == "te":
        t, e, moved = te_ratio(argv[1])
        r = t / max(e, 1)
        print(f"    target peptides {t}  entrapment-only peptides {e}"
              f"  ({moved} entrapment peptides found in a sample protein, counted as target)")
        print(f"    T/E = {r:.3f}  ->  FDP = (hits/total) x {1 + r:.3f}")
        print("    NOTE: this factor is a property of YOUR build. Never assume 1:1 (factor 2).")
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
