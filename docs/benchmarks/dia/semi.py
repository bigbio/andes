# Offline semi-supervised rescoring of a Phase B PIN, PyProphet style: start from a hand-built score,
# then rounds of (positives = best peak group of targets at q <= 1 %, negatives = best peak group of decoys),
# 3-fold CV by precursor, linear discriminant (LDA) or GBDT, score every peak group, re-pick best per precursor.
# Usage: semi.py <dia.pin> [lda|gbdt]
import sys, pandas as pd, numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import HistGradientBoostingClassifier
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
lib = pd.read_csv(B + '/diann-pairlib-20261113/lib_targets.tsv', sep='\t', usecols=['ModifiedPeptideSequence', 'PrecursorCharge', 'ProteinId'], dtype=str).drop_duplicates(['ModifiedPeptideSequence', 'PrecursorCharge'])
lib['ent'] = lib.ProteinId.str.split(';').map(lambda ps: all('_p_target' in p for p in ps))
ent = dict(zip(lib.ModifiedPeptideSequence + '/' + lib.PrecursorCharge, lib.ent)); k = 2.0
kind = sys.argv[2] if len(sys.argv) > 2 else 'lda'
d = pd.read_csv(sys.argv[1], sep='\t')
cz = [c for c in d.columns if c.startswith('charge')]; d['z'] = [int(cz[i][6:]) for i in d[cz].values.argmax(1)]
feats = [c for c in d.columns if c.startswith(('pc_', 'pa_'))] + cz
X = d[feats].values.astype(float)
d['key'] = d.Peptide.str[2:-2].str.replace('+57.02146', '(UniMod:4)', regex=False).str.replace('+15.99491', '(UniMod:35)', regex=False) + '/' + d.z.astype(str)
d['pid'] = d.key + '|' + d.Label.astype(str)
pid_codes, pid_idx = np.unique(d.pid.values, return_inverse=True)
rng = np.random.default_rng(0); fold_of = rng.integers(0, 3, len(pid_codes)); fold = fold_of[pid_idx]
dec = (d.Label == -1).values


def best_rows(score):
    o = np.lexsort((-score, pid_idx))
    first = np.ones(len(o), bool); first[1:] = pid_idx[o][1:] != pid_idx[o][:-1]
    return o[first]


def qvals(rows, score):
    o = rows[np.argsort(-score[rows])]
    dd = dec[o]
    q = np.minimum.accumulate(((np.cumsum(dd) + 1) / np.maximum(np.cumsum(~dd), 1))[::-1])[::-1]
    return o, q


def report(score, tag):
    rows = best_rows(score); o, q = qvals(rows, score)
    t = d.iloc[o][~dec[o]].copy(); t['q'] = q[~dec[o]]
    t['ent'] = t.key.map(ent); t = t.dropna(subset=['ent']); e = t.ent.values.astype(bool)
    fd = np.cumsum(e) * k / np.arange(1, len(e) + 1); ok = np.where(fd <= 0.01)[0]
    s = t[t.q <= 0.01]
    print(f'{tag}: q<=1% N {len(s)} FDP_c {s.ent.sum() * k / max(len(s), 1):.4f} | N@FDP_c1% {ok.max() + 1 if len(ok) else 0}', flush=True)
    return o, q


s0 = (d.pc_corr_mean.clip(0) * d.pc_area_angle * d.pc_apex_matched_frac * (0.5 + d.pc_ms1_found) * np.exp(-0.5 * (d.pc_rt_err / 2) ** 2)).values
score = s0.copy()
o, q = report(score, 'start')
for it in range(4):
    new = np.zeros(len(d))
    for f in range(3):
        tr_rows = o[fold[o] != f]
        tq = q[fold[o] != f]
        pos = tr_rows[(~dec[tr_rows]) & (tq <= 0.01)]
        neg = tr_rows[dec[tr_rows]]
        if len(pos) < 20:
            new[fold == f] = score[fold == f]
            continue
        idx = np.concatenate([pos, neg]); y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
        if kind == 'gbdt':
            m = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.1, class_weight='balanced').fit(X[idx], y)
            new[fold == f] = m.decision_function(X[fold == f])
        else:
            m = LinearDiscriminantAnalysis().fit(X[idx], y)
            new[fold == f] = m.decision_function(X[fold == f])
        # fold-wise scale: decoy best-row median -> 0, decoy best-row sd -> 1
        sub = best_rows(np.where(fold == f, new, -np.inf))
        dd = sub[(fold[sub] == f) & dec[sub]]
        mu, sd = np.median(new[dd]), np.std(new[dd]) + 1e-9
        new[fold == f] = (new[fold == f] - mu) / sd
    score = new
    o, q = report(score, f'{kind} it{it + 1}')
