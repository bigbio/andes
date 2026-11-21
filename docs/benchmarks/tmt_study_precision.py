"""tmt_study_precision.py <name> <study A> <study B> [null pairs a/b,...] : protein-level precision of two TMT study tables.

Rows used: quant_gate_passed, quant_weight > 0, unique peptide, in a protein group. Protein quantity per
(plex, protein, channel) = sum of quant_weight x reporter intensity. Null pairs default to the SDRF
reference channels (two references per plex)."""
import sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

name, pa, pb = sys.argv[1], sys.argv[2], sys.argv[3]
pairs_arg = sys.argv[4] if len(sys.argv) > 4 else None
rng = np.random.default_rng(7)

def load(p):
    t = pq.read_table(p).to_pandas()
    ch = [c[4:] for c in t.columns if c.startswith("tmt_")]
    u = t[(t.quant_weight > 0) & (t["unique"] == True) & t.anchor_protein.notna()]
    w = u.quant_weight.to_numpy()[:, None]
    q = pd.DataFrame(u[[f"tmt_{c}" for c in ch]].to_numpy() * w, columns=ch)
    q["plex"] = u.plex.to_numpy(); q["prot"] = u.anchor_protein.to_numpy()
    prot = q.groupby(["plex", "prot"]).sum()
    design = t.groupby("plex").first()[[f"sample_tmt_{c}" for c in ch] + [f"role_tmt_{c}" for c in ch]]
    return t, u, prot, ch, design

A, B = load(pa), load(pb)
ch = A[3]

def null_pairs(design):
    if pairs_arg:
        return {plex: [tuple(x.split("/")) for x in pairs_arg.split(",")] for plex in design.index}
    out = {}
    for plex, r in design.iterrows():
        refs = [c for c in ch if r[f"role_tmt_{c}"] == "reference"]
        out[plex] = [(refs[0], refs[1])] if len(refs) >= 2 else []
    return out

def nulls(prot, pairs, keep=None):
    vals = {}
    for plex, ps in pairs.items():
        if plex not in prot.index.get_level_values(0):
            continue
        sub = prot.loc[plex]
        for a, b in ps:
            ok = (sub[a] > 0) & (sub[b] > 0)
            x = np.log2(sub[a][ok] / sub[b][ok])
            x = x - np.median(x)
            for k, v in x.items():
                vals[(plex, a, b, k)] = v
    s = pd.Series(vals)
    if keep is not None:
        s = s[s.index.isin(keep)]
    return s

pairs = null_pairs(A[4])
na, nb = nulls(A[2], pairs), nulls(B[2], pairs)
common = na.index.intersection(nb.index)
def stats(x):
    return f"SD {x.std():.4f}  robust SD {1.4826*np.median(np.abs(x-np.median(x))):.4f}  |x|>0.5 {(x.abs()>0.5).mean()*100:.2f} %  n {len(x)}"
print(f"== {name}: null pairs {pairs}")
print(f"rows used A {len(A[1])} B {len(B[1])}; proteins quantified A {A[2].index.get_level_values(1).nunique()} B {B[2].index.get_level_values(1).nunique()}")
print(f"null log2 ratio, all proteins   A: {stats(na)}")
print(f"null log2 ratio, all proteins   B: {stats(nb)}")
ca, cb = na[common], nb[common]
print(f"null log2 ratio, common proteins A: {stats(ca)}")
print(f"null log2 ratio, common proteins B: {stats(cb)}")
d = []
idx = np.arange(len(common))
for _ in range(500):
    i = rng.choice(idx, len(idx))
    d.append(cb.iloc[i].std() - ca.iloc[i].std())
lo, hi = np.percentile(d, [2.5, 97.5])
print(f"SD(B) - SD(A) on common proteins: {cb.std()-ca.std():+.4f}  95 % bootstrap CI [{lo:+.4f}, {hi:+.4f}]")

def cross_instrument(prot, design):
    plexes = list(design.index)
    if len(plexes) < 2:
        return None
    res = []
    samp = {p: {design.loc[p, f"sample_tmt_{c}"]: c for c in ch if design.loc[p, f"role_tmt_{c}"] == "sample"} for p in plexes}
    shared = set(samp[plexes[0]]) & set(samp[plexes[1]])
    for s in sorted(shared):
        v = []
        for p in plexes[:2]:
            sub = prot.loc[p]
            refs = [c for c in ch if design.loc[p, f"role_tmt_{c}"] == "reference"]
            ref = sub[refs].mean(axis=1)
            x = sub[samp[p][s]]
            ok = (x > 0) & (ref > 0)
            v.append(np.log2(x[ok] / ref[ok]))
        j = v[0].index.intersection(v[1].index)
        res.append((s, len(j), np.corrcoef(v[0][j], v[1][j])[0, 1]))
    return res

for arm, X in (("A", A), ("B", B)):
    r = cross_instrument(X[2], X[4])
    if r:
        print(f"same sample on both instruments, Pearson of log2(sample/reference) {arm}: " + "; ".join(f"{s} n {n} r {c:.4f}" for s, n, c in r) + f"; mean r {np.mean([c for _,_,c in r]):.4f}")
    for plex, row in X[4].iterrows():
        blanks = [c for c in ch if row[f"role_tmt_{c}"] == "blank"]
        if not blanks or plex not in X[2].index.get_level_values(0):
            continue
        sub = X[2].loc[plex]
        tot = sub[ch].sum(axis=1)
        share = sub[blanks].sum(axis=1)[tot > 0] / tot[tot > 0]
        print(f"blank share {arm} plex {plex} ({blanks}): median per protein {share.median()*100:.3f} %, intensity-weighted {sub[blanks].to_numpy().sum()/tot.sum()*100:.3f} %, proteins {len(share)}")

def ratios(prot, design, s):
    out = []
    for p in list(design.index)[:2]:
        sub = prot.loc[p]
        samp = {design.loc[p, f"sample_tmt_{c}"]: c for c in ch if design.loc[p, f"role_tmt_{c}"] == "sample"}
        refs = [c for c in ch if design.loc[p, f"role_tmt_{c}"] == "reference"]
        ref = sub[refs].mean(axis=1); x = sub[samp[s]]
        ok = (x > 0) & (ref > 0)
        out.append(np.log2(x[ok] / ref[ok]))
    return out

if len(A[4]) >= 2:
    plexes = list(A[4].index)
    shared = [s for s in A[4].loc[plexes[0]].values if s in set(A[4].loc[plexes[1]].values) and str(s).startswith("PXD") and "GIS" not in str(s) and "blank" not in str(s)]
    for s in sorted(set(shared)):
        ra, rb = ratios(A[2], A[4], s), ratios(B[2], B[4], s)
        j = ra[0].index.intersection(ra[1].index).intersection(rb[0].index).intersection(rb[1].index)
        va = (ra[0][j].to_numpy(), ra[1][j].to_numpy()); vb = (rb[0][j].to_numpy(), rb[1][j].to_numpy())
        r = lambda v, i: np.corrcoef(v[0][i], v[1][i])[0, 1]
        full = np.arange(len(j)); d = []
        for _ in range(500):
            i = rng.choice(full, len(full)); d.append(r(vb, i) - r(va, i))
        lo, hi = np.percentile(d, [2.5, 97.5])
        print(f"common proteins {s}: n {len(j)} r A {r(va, full):.4f} B {r(vb, full):.4f} diff {r(vb, full)-r(va, full):+.4f} CI [{lo:+.4f}, {hi:+.4f}]")
