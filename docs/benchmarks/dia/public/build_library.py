"""Paired-entrapment one-peptide FASTA and library table from the UniProt human proteome.

Deterministic for a given proteome file, --n and --seed:

- targets: N distinct fully tryptic human peptides (cleavage after K/R, not before P, no missed
  cleavage, C-terminal K or R, length 7-35, only the 20 standard residues, first residue not M),
  distinct with I and L taken as equal, sampled uniformly with random.Random(seed) from the
  sorted set of all such peptides;
- entrapment: for each target one partner, its residues except the C-terminal one shuffled with
  the same generator, accepted when (I/L folded) it differs from its target, does not start
  with M, is not a fully tryptic human peptide of any length, does not occur anywhere in a
  human protein sequence, and is not already in the database; a target whose partner fails
  100 draws is dropped.
- No entry starts with M: andes removes an initiator Met from a FASTA entry, which would add
  clipped forms that are not in the library. Partner ids are <target id>_p_target, so every
  entrapment precursor has a target precursor of identical composition and r = 1;
- library table (eval_dia.py): every peptidoform of every entry with Carbamidomethyl (C) fixed
  and Oxidation (M) variable, at most one, at charges 1..5 whose m/z lies in [--mz-lo, --mz-hi].

Usage: build_library.py <uniprot.fasta[.gz]> <out.fasta> <out_lib.tsv> [--n N] [--seed S]
"""
import argparse, gzip, random, sys

MONO = {
    'G': 57.021463720, 'A': 71.037113785, 'S': 87.032028405, 'P': 97.052763850,
    'V': 99.068413945, 'T': 101.047678505, 'C': 103.009184505, 'L': 113.084064015,
    'I': 113.084064015, 'N': 114.042927470, 'D': 115.026943065, 'Q': 128.058577540,
    'K': 128.094963050, 'E': 129.042593135, 'M': 131.040484645, 'H': 137.058911875,
    'F': 147.068413945, 'R': 156.101111050, 'Y': 163.063328575, 'W': 186.079312980,
}
WATER = 18.010564684
PROTON = 1.007276466812
CAM = 57.021464
OX = 15.994915
K = 7


def read_fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    seq = []
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if seq:
                    yield ''.join(seq)
                seq = []
            else:
                seq.append(line.strip())
    if seq:
        yield ''.join(seq)


def tryptic(seq):
    start = 0
    for i, aa in enumerate(seq):
        if aa in 'KR' and (i + 1 == len(seq) or seq[i + 1] != 'P'):
            yield seq[start:i + 1]
            start = i + 1
    if start < len(seq):
        yield seq[start:]


def fold(p):
    return p.replace('I', 'L')


def forms(pep):
    base = [('C(UniMod:4)' if a == 'C' else a) for a in pep]
    yield ''.join(base), 0
    for i, a in enumerate(pep):
        if a == 'M':
            yield ''.join(base[:i] + ['M(UniMod:35)'] + base[i + 1:]), 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('proteome')
    ap.add_argument('fasta')
    ap.add_argument('lib')
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--seed', type=int, default=20261119)
    ap.add_argument('--min-length', type=int, default=7)
    ap.add_argument('--max-length', type=int, default=35)
    ap.add_argument('--mz-lo', type=float, default=370.0)
    ap.add_argument('--mz-hi', type=float, default=990.0)
    a = ap.parse_args()

    all_tryptic = set()
    candidates = {}
    proteins = list(read_fasta(a.proteome))
    folded = '|'.join(proteins).replace('I', 'L')
    kmers = {folded[i:i + K] for i in range(len(folded) - K + 1)}

    def in_proteome(fp):
        return (all(fp[i:i + K] in kmers for i in range(len(fp) - K + 1))
                and fp in folded)

    for prot in proteins:
        for p in tryptic(prot):
            all_tryptic.add(fold(p))
            if (a.min_length <= len(p) <= a.max_length and p[-1] in 'KR' and p[0] != 'M'
                    and all(c in MONO for c in p)):
                candidates.setdefault(fold(p), p)
    pool = sorted(candidates.values())
    rng = random.Random(a.seed)
    picked = rng.sample(pool, a.n)
    used = {fold(p) for p in picked}
    pairs = []
    dropped = 0
    for i, t in enumerate(picked):
        body = list(t[:-1])
        for _ in range(100):
            rng.shuffle(body)
            e = ''.join(body) + t[-1]
            fe = fold(e)
            if (e[0] != 'M' and fe != fold(t) and fe not in all_tryptic and fe not in used
                    and not in_proteome(fe)):
                used.add(fe)
                pairs.append((f'T{i:06d}', t, e))
                break
        else:
            dropped += 1
    with open(a.fasta, 'w') as fh:
        for pid, t, e in pairs:
            fh.write(f'>{pid}\n{t}\n>{pid}_p_target\n{e}\n')
    n_prec = [0, 0]
    with open(a.lib, 'w') as fh:
        fh.write('ModifiedPeptideSequence\tPrecursorCharge\tProteinId\tPrecursorMz\n')
        for pid, t, e in pairs:
            for k, (prot, pep) in enumerate(((pid, t), (pid + '_p_target', e))):
                m0 = sum(MONO[c] for c in pep) + WATER + CAM * pep.count('C')
                for mod, nox in forms(pep):
                    m = m0 + OX * nox
                    for z in range(1, 6):
                        mz = (m + z * PROTON) / z
                        if a.mz_lo <= mz <= a.mz_hi:
                            fh.write(f'{mod}\t{z}\t{prot}\t{mz:.6f}\n')
                            n_prec[k] += 1
    print(f'tryptic peptides {len(all_tryptic)}; candidates {len(pool)}; targets {len(pairs)} '
          f'(dropped {dropped}); target precursors {n_prec[0]}; entrapment precursors '
          f'{n_prec[1]}; r = {n_prec[1] / n_prec[0]:.6f}', file=sys.stderr)


main()
