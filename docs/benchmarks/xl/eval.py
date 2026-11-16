"""Score an andes --xl result (<base>.xl.tsv) against a dataset's ground truth.

Truth kinds
  lib:<name>      synthetic peptide library (IMP-X-FDR xlsx): a TT CSM is true when both
                  peptides map (substring) onto library peptides of a common group with the
                  link on each library peptide's designated residue.
  clasen:<batchN> recombinant protein standard (Clasen et al. 2024): a TT CSM is true when
                  both proteins are in one interaction group of that batch (intra: the
                  protein is in that batch).
Reports per link class (intra, inter, all) at andes q <= 1 %: accepted TT, false TT, false
share, true unique residue pairs, and the true CSM count at an actual false share of 1 %
ranked by andes score (discrimination independent of the FDR estimate).
"""
import collections, json, re, sys
import pandas as pd

import os
TRUTH = os.path.join(os.environ.get('XLBENCH_DATA', '.'), 'truth') + '/'


def load_lib(name):
    df = pd.read_excel(f'{TRUTH}{name}.xlsx', header=None)
    groups, cur = collections.defaultdict(list), None
    for a, b in df.itertuples(index=False):
        a = str(a).strip()
        if a.lower().startswith('group'):
            cur = a
            continue
        if a != 'nan':
            groups[cur].append((a, int(b)))
    return groups


def load_clasen(batch):
    df = pd.read_excel(f'{TRUTH}protein-construct-sequencesandinteraction-groups-corrected.xlsx')
    g = collections.defaultdict(set)
    for acc, grp in zip(df['Uniprot ID'], df['group']):
        if str(grp).endswith('_' + batch):
            g[grp].add(str(acc).strip())
    return g


def bare(proforma):
    s = re.sub(r'\[[^\]]*\]', '', str(proforma))
    parts = s.split('.')
    if len(parts) == 3:
        s = parts[1]
    s = s.replace('-', '')
    return ''.join(c for c in s if c.isalpha() and c.isupper())


def acc_of(token):
    p = token.split('|')
    return p[1] if len(p) >= 3 else token


class LibTruth:
    def __init__(self, name):
        self.groups = load_lib(name)

    def groups_of(self, pep, pos, kind):
        if kind not in ('K',):
            return set()
        out = set()
        for g, members in self.groups.items():
            for lib, idx in members:
                o = lib.find(pep)
                while o >= 0:
                    if o + pos == idx:
                        out.add(g)
                    o = lib.find(pep, o + 1)
        return out

    def is_true(self, r):
        a = self.groups_of(bare(r.alpha_peptide), int(r.alpha_site) - 1, r.alpha_site_kind)
        b = self.groups_of(bare(r.beta_peptide), int(r.beta_site) - 1, r.beta_site_kind)
        return bool(a & b)


class ClasenTruth:
    def __init__(self, batch):
        self.groups = load_clasen(batch)
        self.batch_prots = set().union(*self.groups.values())

    def is_true(self, r):
        pa = {acc_of(t) for t in str(r.alpha_proteins).split(';')}
        pb = {acc_of(t) for t in str(r.beta_proteins).split(';')}
        for a in pa:
            for b in pb:
                if a == b and a in self.batch_prots:
                    return True
                if any(a in g and b in g for g in self.groups.values()):
                    return True
        return False


def truth_for(spec):
    kind, arg = spec.split(':', 1)
    return LibTruth(arg) if kind == 'lib' else ClasenTruth(arg)


def at_actual(df, fdp=0.01):
    d = df.sort_values('score', ascending=False)
    t = f = best = 0
    for v in d['true']:
        t += v
        f += not v
        if f <= fdp * (t + f):
            best = max(best, t)
    return best


def evaluate(tsv, truth_spec, q=0.01):
    df = pd.read_csv(tsv, sep='\t', low_memory=False)
    tr = truth_for(truth_spec)
    tt = df[df['class'] == 'TT'].copy()
    tt['true'] = [tr.is_true(r) for r in tt.itertuples(index=False)]
    tt['rp'] = [tuple(sorted([(bare(r.alpha_peptide), int(r.alpha_site)), (bare(r.beta_peptide), int(r.beta_site))]))
                for r in tt.itertuples(index=False)]
    out = {}
    for cls in ('intra', 'homotypic', 'inter', 'all'):
        sub = tt if cls == 'all' else tt[tt['link'] == cls]
        acc = sub[sub['q_value'].astype(float) <= q]
        n, nf = len(acc), int((~acc['true']).sum())
        out[cls] = dict(
            csm=n, false=nf, fdp=round(nf / n, 4) if n else 0.0,
            true_rp=int(acc[acc['true']]['rp'].nunique()),
            false_rp=int(acc[~acc['true']]['rp'].nunique()),
            true_at_1pct_actual=at_actual(sub) if len(sub) else 0)
    out['all']['csm_rows'] = len(df)
    return out


if __name__ == '__main__':
    res = evaluate(sys.argv[1], sys.argv[2])
    print(json.dumps(res))
