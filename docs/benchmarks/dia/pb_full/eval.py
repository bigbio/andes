"""Precursor-level entrapment evaluation of andes --dia (or DIA-NN) runs, with the quantities of
../eval_dia.py, for any paired-entrapment FASTA.

A precursor is (peptidoform in UniMod notation, charge), scored by its best PSM. The q-value is
precursor-level target-decoy competition (D + 1) / T, made monotone. Entrapment FDP_c = E (1 + 1/r) / N.
N@FDP_c 1 % / 2 % ranks the target precursors by score and counts how many are kept before the
combined FDP passes 1 % / 2 %.

Which precursors are entrapment, and r:
  - default: a precursor is entrapment when every protein of its PSM carries the entrapment tag
    (--tag, default `_p_target`); r is --r (default 1, the ProteoBench paired FASTA);
  - --library <tsv> (the PAIRLIB layout, columns ModifiedPeptideSequence, PrecursorCharge, ProteinId):
    only precursors of the library count, entrapment comes from the library's proteins and r from
    its entrapment / target precursor counts. This reproduces ../eval_dia.py.

Usage:
  eval.py andes <pin> <percolator target psms> <percolator decoy psms> [--library TSV] [--r R] [--tag T]
  eval.py diann <report.parquet> <run name substring> [--tag T] [--r R]
Prints one line: precursors, N / E / FDP_c at q <= 0.01 and 0.05, N@FDP_c 1 % / 2 %.
"""
import argparse
import re

import numpy as np
import pandas as pd


def unimod(pep):
    core = pep[2:-2] if len(pep) > 4 and pep[1] == '.' and pep[-2] == '.' else pep
    core = re.sub(r'\[?\+15\.99\d*\]?', '(UniMod:35)', core)
    return re.sub(r'\[?\+57\.02\d*\]?', '(UniMod:4)', core)


def n_at(ent, k, t):
    f = np.cumsum(ent) * k / np.arange(1, len(ent) + 1)
    ok = np.where(f <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0


def summary(name, q, ent, k, extra=''):
    out = [name, f'precursors {len(q)}{extra}']
    for qq in (0.01, 0.05):
        m = q <= qq
        e = int(ent[m].sum())
        out.append(f'q<={qq}: N {int(m.sum())} E {e} FDP_c {e * k / max(int(m.sum()), 1):.4f}')
    out.append(f'N@FDP_c 1%/2%: {n_at(ent, k, 0.01)}/{n_at(ent, k, 0.02)}')
    print(' | '.join(out))


def tagged(proteins, tag):
    ps = [p for p in re.split(r'[;,\t]', str(proteins)) if p]
    return bool(ps) and all(tag in p for p in ps)


def read_psms(path):
    """Percolator --only-psms output: the protein ids fill every column after the fifth, so rows have
    as many fields as the PSM has proteins and cannot be read as a fixed-width table."""
    rows = []
    with open(path) as fh:
        h = fh.readline().rstrip('\n').split('\t')
        n = h.index('proteinIds')
        for line in fh:
            x = line.rstrip('\n').split('\t')
            rows.append(x[:n] + ['\t'.join(x[n:])])
    d = pd.DataFrame(rows, columns=h[:n + 1])
    d['score'] = d.score.astype(float)
    d['q-value'] = d['q-value'].astype(float)
    return d


def pin_charges(pin):
    charge = {}
    with open(pin) as fh:
        h = fh.readline().rstrip('\n').split('\t')
        cz = [(i, int(c[6:])) for i, c in enumerate(h) if re.fullmatch(r'charge\d+', c)]
        last = max(i for i, _ in cz)
        for line in fh:
            x = line.split('\t', last + 1)
            charge[x[0]] = next((z for i, z in cz if x[i] == '1'), 0)
    return charge


def andes(a):
    charge = pin_charges(a.pin)
    parts = []
    for path, decoy in ((a.targets, False), (a.decoys, True)):
        d = read_psms(path)[['PSMId', 'score', 'peptide', 'proteinIds']]
        d['decoy'] = decoy
        parts.append(d)
    d = pd.concat(parts)
    d['key'] = d.peptide.map(unimod) + '/' + d.PSMId.map(charge).fillna(0).astype(int).astype(str)
    d = d.sort_values('score', ascending=False, kind='mergesort').drop_duplicates(['key', 'decoy'])
    dec = d.decoy.values
    d['q'] = np.minimum.accumulate(((np.cumsum(dec) + 1) / np.maximum(np.cumsum(~dec), 1))[::-1])[::-1]
    t = d[~d.decoy].copy()
    extra = ''
    if a.library:
        lib = pd.read_csv(a.library, sep='\t', dtype=str,
                          usecols=['ModifiedPeptideSequence', 'PrecursorCharge', 'ProteinId'])
        lib = lib.drop_duplicates(['ModifiedPeptideSequence', 'PrecursorCharge'])
        lib['ent'] = lib.ProteinId.str.split(';').map(lambda ps: all(a.tag in p for p in ps))
        ent = lib.set_index(lib.ModifiedPeptideSequence + '/' + lib.PrecursorCharge).ent
        r = ent.sum() / (~ent).sum()
        t['ent'] = t.key.map(ent)
        extra = f' (not in library {int(t.ent.isna().sum())})'
        t = t.dropna(subset=['ent'])
    else:
        r = a.r
        t['ent'] = t.proteinIds.map(lambda p: tagged(p, a.tag))
    summary('andes', t.q.values, t.ent.values.astype(bool), 1 + 1 / r, extra)


def diann(a):
    r = pd.read_parquet(a.report, columns=['Run', 'Precursor.Id', 'Protein.Ids', 'Q.Value', 'PEP'])
    r = r[r.Run.str.contains(a.run, regex=False)]
    r = r.sort_values(['Q.Value', 'PEP']).drop_duplicates('Precursor.Id')
    ent = r['Protein.Ids'].map(lambda p: tagged(p, a.tag)).values.astype(bool)
    summary(f'diann {a.run}', r['Q.Value'].values, ent, 1 + 1 / a.r)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='tool', required=True)
    pa = sub.add_parser('andes')
    pa.add_argument('pin')
    pa.add_argument('targets')
    pa.add_argument('decoys')
    pa.add_argument('--library')
    pd_ = sub.add_parser('diann')
    pd_.add_argument('report')
    pd_.add_argument('run')
    for s in (pa, pd_):
        s.add_argument('--r', type=float, default=1.0)
        s.add_argument('--tag', default='_p_target')
    a = p.parse_args()
    andes(a) if a.tool == 'andes' else diann(a)
