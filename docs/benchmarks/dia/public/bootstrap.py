"""bootstrap.py <runA_dir> <runB_dir> <raw> [B] [seed]: bootstrap over precursors of the paired
difference B - A in N@FDP_c 1 % and 2 % (mean over the Percolator seeds both runs share).

A precursor is scored as in eval_dia.py (best target PSM per (UniMod sequence, charge), library
precursors only, entrapment = every protein _p_target, FDP_c = E (1 + 1/r) / N). Each bootstrap
replicate draws a Poisson(1) weight for every precursor of the union of both runs and all seeds,
so a precursor keeps its scores in both runs and in every seed together; N@FDP_c is computed on
the weighted ranked lists, per seed, and the seed mean of B - A is the replicate's value.
Prints: point estimate (weights 1), bootstrap mean, sd, 2.5 / 97.5 % quantiles.
With --dump <npz>, also writes the replicate differences, for pooling over slices.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.environ.get('BENCH_ROOT', '/home/user/dia-bench')
LIB = ROOT + '/data/frozen/pairlib_public_lib.tsv'

args = [a for a in sys.argv[1:] if not a.startswith('--')]
dump = None
if '--dump' in sys.argv:
    dump = sys.argv[sys.argv.index('--dump') + 1]
    args.remove(dump)
run_a, run_b, raw = args[:3]
n_boot = int(args[3]) if len(args) > 3 else 1000
rng = np.random.default_rng(int(args[4]) if len(args) > 4 else 20261120)

lib = pd.read_csv(LIB, sep='\t', usecols=['ModifiedPeptideSequence', 'PrecursorCharge', 'ProteinId'],
                  dtype=str).drop_duplicates(['ModifiedPeptideSequence', 'PrecursorCharge'])
lib['ent'] = lib.ProteinId.str.split(';').map(lambda ps: all('_p_target' in p for p in ps))
lib['key'] = lib.ModifiedPeptideSequence + '/' + lib.PrecursorCharge
ent = lib.set_index('key').ent
k = 1 + 1 / (ent.sum() / (~ent).sum())


def unimod(pep):
    core = pep[2:-2] if len(pep) > 4 and pep[1] == '.' and pep[-2] == '.' else pep
    return re.sub(r'\+15\.9949[12]', '(UniMod:35)', core).replace('+57.02146', '(UniMod:4)')


def charges(pin):
    out = {}
    with open(pin) as fh:
        h = fh.readline().rstrip('\n').split('\t')
        cz = [(i, int(c[6:])) for i, c in enumerate(h) if re.fullmatch(r'charge\d+', c)]
        last = max(i for i, _ in cz)
        for line in fh:
            x = line.split('\t', last + 1)
            out[x[0]] = next(z for i, z in cz if x[i] == '1')
    return out


def targets(run, seed, charge):
    d = pd.read_csv(f'{run}/{raw}/seed{seed}/t.psms', sep='\t', usecols=['PSMId', 'score', 'peptide'])
    d['key'] = d.peptide.map(unimod) + '/' + d.PSMId.map(charge).astype(str)
    d = d.sort_values('score', ascending=False, kind='mergesort').drop_duplicates('key')
    d['ent'] = d.key.map(ent)
    return d.dropna(subset=['ent'])[['key', 'ent']]


seeds = ['42', '1', '2', '3', '7']
lists = {}
for run in (run_a, run_b):
    ch = charges(f'{run}/{raw}/{raw}.pin')
    lists[run] = [targets(run, s, ch) for s in seeds]
keys = pd.Index(sorted(set().union(*(set(t.key) for r in lists.values() for t in r))))
ranked = {run: [(keys.get_indexer(t.key), t.ent.values.astype(bool)) for t in ts]
          for run, ts in lists.items()}


def n_at(idx, e, w, t):
    ww = w[idx]
    cw = np.cumsum(ww)
    ce = np.cumsum(ww * e)
    ok = (ce * k <= t * cw) & (cw > 0)
    return cw[ok].max() if ok.any() else 0.0


def diffs(w):
    out = []
    for t in (0.01, 0.02):
        d = [n_at(*ranked[run_b][s], w, t) - n_at(*ranked[run_a][s], w, t) for s in range(len(seeds))]
        out.append(np.mean(d))
    return out


point = diffs(np.ones(len(keys)))
reps = np.array([diffs(rng.poisson(1.0, len(keys)).astype(float)) for _ in range(n_boot)])
if dump:
    np.savez(dump, point=np.array(point), reps=reps)
for j, name in enumerate(('N@FDP_c 1%', 'N@FDP_c 2%')):
    lo, hi = np.quantile(reps[:, j], [0.025, 0.975])
    print(f'{raw} {name}: diff {point[j]:+.0f}, boot sd {reps[:, j].std():.0f}, '
          f'95% [{lo:+.0f}, {hi:+.0f}]')
