import sys, numpy as np, pandas as pd, pyarrow.parquet as pq
lab, run, ffm, design = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
t = pq.read_table(lab).to_pandas(); t = t[t.filename == run]
ref = t[t.channel == 1][['peptide','charge','mz','rt_apex']]
rows=[]
for _, r in pd.read_csv(ffm, sep='\t').iterrows():
    ch = {int(a):(float(b),float(c)) for a,b,c in (x.split(':') for x in r.channels.split(';'))}
    if 0 in ch and 1 in ch and ch[0][1]>0 and ch[1][1]>0:
        rows.append((r.rt, r.charge, ch[0][0], np.log2(ch[1][1]/ch[0][1])))
f = pd.DataFrame(rows, columns=['rt','charge','mz','o']).sort_values('mz').reset_index(drop=True)
iqr = lambda x: x.quantile(.75)-x.quantile(.25)
print(f'{run}: OpenMS complete multiplets {len(f)} median {f.o.median():.3f} IQR {iqr(f.o):.3f} |error| median {np.median(np.abs(f.o-design)):.3f}')
a = t[(t.channel==2) & t.log2_ratio.notna()]
print(f'{run}: andes channel ratios {len(a)} median {a.log2_ratio.median():.3f} IQR {iqr(a.log2_ratio):.3f} |error| median {np.median(np.abs(a.log2_ratio-design)):.3f}')
m=[]; mz=f.mz.values
for _, g in ref.iterrows():
    lo, hi = np.searchsorted(mz, [g.mz*(1-20e-6), g.mz*(1+20e-6)])
    c = f.iloc[lo:hi]; c = c[(c.charge==g.charge)&((c.rt-g.rt_apex).abs()<=30)]
    if len(c): m.append((g.peptide, g.charge, c.iloc[(c.rt-g.rt_apex).abs().argmin()].o))
m = pd.DataFrame(m, columns=['peptide','charge','o']).merge(a[['peptide','charge','log2_ratio']], on=['peptide','charge']).dropna()
print(f'{run}: same peptides {len(m)}: OpenMS median {m.o.median():.3f} IQR {iqr(m.o):.3f} |error| {np.median(np.abs(m.o-design)):.3f} | andes median {m.log2_ratio.median():.3f} IQR {iqr(m.log2_ratio):.3f} |error| {np.median(np.abs(m.log2_ratio-design)):.3f}')
