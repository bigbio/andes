"""Where andes loses DIA-NN's precursors in one run, by DIA-NN intensity quartile.

gap.py <pin> <percolator target psms> <percolator decoy psms> <DIA-NN report.parquet> <run substring>
       <scans.parquet from scans.py> [--tag T]

DIA-NN's precursors at run-level q <= 0.01 whose proteins are not all entrapment (--tag) fall into:
  found              andes precursor q <= 0.01 (best PSM, precursor-level target-decoy competition);
  ranked first at RT a PIN row within 0.1 min of DIA-NN's RT that is its scan's best row (RawScore);
  in PIN at RT       a PIN row within 0.1 min, never its scan's best;
  away from RT       PIN rows only further from DIA-NN's RT;
  not in PIN         no PIN row.
A precursor is (stripped sequence, charge); the first-pass report of a DIA-NN run with --reanalyse is
the per-run (no match-between-runs) comparison.
"""
import argparse
import re

import numpy as np
import pandas as pd

a = argparse.ArgumentParser()
for n in ('pin', 'targets', 'decoys', 'report', 'run', 'scans'):
    a.add_argument(n)
a.add_argument('--tag', default='_p_target')
a = a.parse_args()


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


def core(p):
    c = p[2:-2] if len(p) > 4 and p[1] == '.' and p[-2] == '.' else p
    return re.sub(r'[^A-Z]', '', c)


def tagged(proteins):
    ps = [x for x in re.split(r'[;,\t]', str(proteins)) if x]
    return bool(ps) and all(a.tag in x for x in ps)


h = open(a.pin).readline().rstrip('\n').split('\t')
cz = [c for c in h if re.fullmatch(r'charge\d+', c)]
p = pd.read_csv(a.pin, sep='\t', usecols=['SpecId', 'Label', 'ScanNr', 'RawScore', 'Peptide'] + cz)
p['z'] = [int(cz[i][6:]) for i in p[cz].values.argmax(1)]
p['key'] = p.Peptide.map(core) + '/' + p.z.astype(str)
sc = pd.read_parquet(a.scans)
p['rt'] = p.ScanNr.map(dict(zip(sc.scan, sc.rt)))
p['first'] = p.RawScore >= p.groupby('ScanNr').RawScore.transform('max')

parts = []
for path, decoy in ((a.targets, False), (a.decoys, True)):
    d = read_psms(path)[['PSMId', 'score']]
    d['decoy'] = decoy
    parts.append(d)
s = pd.concat(parts).merge(p[['SpecId', 'key']], left_on='PSMId', right_on='SpecId')
s = s.sort_values('score', ascending=False).drop_duplicates(['key', 'decoy'])
dec = s.decoy.values
s['q'] = np.minimum.accumulate(((np.cumsum(dec) + 1) / np.maximum(np.cumsum(~dec), 1))[::-1])[::-1]
found = set(s[~s.decoy & (s.q <= 0.01)].key)

r = pd.read_parquet(a.report, columns=['Run', 'Stripped.Sequence', 'Precursor.Charge', 'RT', 'Q.Value',
                                       'Precursor.Quantity', 'Protein.Ids'])
r = r[r.Run.str.contains(a.run, regex=False) & (r['Q.Value'] <= 0.01) & ~r['Protein.Ids'].map(tagged)]
r['key'] = r['Stripped.Sequence'] + '/' + r['Precursor.Charge'].astype(str)
r = r.sort_values('Q.Value').drop_duplicates('key')
r['quartile'] = pd.qcut(r['Precursor.Quantity'].rank(method='first'), 4, labels=['Q1 low', 'Q2', 'Q3', 'Q4 high'])
tp = p[p.Label == 1].merge(r[['key', 'RT']], on='key', how='left')
at = tp[(tp.rt - tp.RT).abs() <= 0.1]
r['category'] = np.select(
    [r.key.isin(found), r.key.isin(set(at[at['first']].key)), r.key.isin(set(at.key)), r.key.isin(set(tp.key))],
    ['found', 'ranked first at RT', 'in PIN at RT', 'away from RT'], 'not in PIN')
print(f'DIA-NN precursors at q <= 0.01 (non-entrapment), run {a.run}: {len(r)}')
print(pd.crosstab(r.category, r.quartile, margins=True).to_string())
