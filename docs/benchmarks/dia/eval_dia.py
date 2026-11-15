"""Precursor-level evaluation of an andes --dia run on the ProteoBench Astral entrapment PAIRLIB targets,
matching OpenDIA's diann-decoys-20261114/evaluate.py: entrapment = every protein of the precursor carries
_p_target; FDP_c = E (1 + 1/r) / N with r = entrapment / human target precursors of the PAIRLIB library;
N@FDP_c 1 %/2 % ranks targets by score. Only precursors of the PAIRLIB library count (the universe
DIA-NN searched).

A precursor is (modified sequence in UniMod notation, charge); its score is the Percolator score of its
best PSM. The q-value is precursor-level target-decoy competition over those best scores (decoy
precursors from the decoy PSM file), (D + 1) / T, monotone.

Usage: eval_dia.py <pin> <percolator.target.psms.txt> <percolator.decoy.psms.txt>
"""
import sys, re, numpy as np, pandas as pd

B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
lib = pd.read_csv(B + '/diann-pairlib-20261113/lib_targets.tsv', sep='\t',
                  usecols=['ModifiedPeptideSequence', 'PrecursorCharge', 'ProteinId'], dtype=str)
lib = lib.drop_duplicates(['ModifiedPeptideSequence', 'PrecursorCharge'])
lib['ent'] = lib.ProteinId.str.split(';').map(lambda ps: all('_p_target' in p for p in ps))
lib['key'] = lib.ModifiedPeptideSequence + '/' + lib.PrecursorCharge
ent = lib.set_index('key').ent
r = ent.sum() / (~ent).sum()
k = 1 + 1 / r


def unimod(pep):
    core = pep[2:-2] if len(pep) > 4 and pep[1] == '.' and pep[-2] == '.' else pep
    return core.replace('+57.02146', '(UniMod:4)').replace('+15.99491', '(UniMod:35)')


def n_at(e, t):
    f = np.cumsum(e) * k / np.arange(1, len(e) + 1)
    ok = np.where(f <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0


def load(path, decoy):
    d = pd.read_csv(path, sep='\t', usecols=['PSMId', 'score', 'peptide'])
    d['decoy'] = decoy
    return d


pin, tpath, dpath = sys.argv[1:4]
charge = {}
with open(pin) as fh:
    h = fh.readline().rstrip('\n').split('\t')
    cz = [(i, int(c[6:])) for i, c in enumerate(h) if re.fullmatch(r'charge\d+', c)]
    last = max(i for i, _ in cz)
    for line in fh:
        x = line.split('\t', last + 1)
        charge[x[0]] = next(z for i, z in cz if x[i] == '1')
d = pd.concat([load(tpath, False), load(dpath, True)])
d['key'] = d.peptide.map(unimod) + '/' + d.PSMId.map(charge).astype(str)
d = d.sort_values('score', ascending=False).drop_duplicates(['key', 'decoy'])
dec = d.decoy.values
q = np.minimum.accumulate(((np.cumsum(dec) + 1) / np.maximum(np.cumsum(~dec), 1))[::-1])[::-1]
d['q'] = q
t = d[~d.decoy].copy()
t['ent'] = t.key.map(ent)
unmapped = int(t.ent.isna().sum())
t = t.dropna(subset=['ent'])
e = t.ent.values.astype(bool)
out = [f'precursors {len(t)} (not in library {unmapped})']
for qq in (0.01, 0.05):
    s = t[t.q <= qq]
    out.append(f'q<={qq}: N {len(s)} E {int(s.ent.sum())} FDP_c {s.ent.sum() * k / max(len(s), 1):.4f}')
out.append(f'N@FDP_c 1%/2%: {n_at(e, 0.01)}/{n_at(e, 0.02)}')
print(' | '.join(out))
