import os, pandas as pd, re, sys, numpy as np
tag=sys.argv[1]
B=os.environ['B']
r=pd.read_parquet(B+'/diann-pairlib-20261113/report_raw01.parquet',columns=['Modified.Sequence','Precursor.Charge','Q.Value','RT'])
r=r[(r['Q.Value']<=0.01)&(r.RT>=10.1)&(r.RT<=12.9)]
dk=set(r['Modified.Sequence']+'/'+r['Precursor.Charge'].astype(str))
cols=open(tag+'.pin').readline().rstrip('\n').split('\t'); k=cols.index('Peptide')
cz=[c for c in cols if re.fullmatch(r'charge\d',c)]
p=pd.read_csv(tag+'.pin',sep='\t',usecols=['SpecId','Label','ScanNr','RawScore','Peptide']+cz)
p['z']=[int(cz[i][6:]) for i in p[cz].values.argmax(1)]
p['rank']=p.SpecId.str.extract(r'_(\d+)(?:_\d+)?$')[0].astype(int)
p['key']=p.Peptide.str[2:-2].str.replace('+57.02146','(UniMod:4)',regex=False).str.replace('+15.99491','(UniMod:35)',regex=False)+'/'+p.z.astype(str)
t=pd.concat([pd.read_csv(tag+'.t.psms',sep='\t',usecols=['PSMId','q-value'])])
q=dict(zip(t.PSMId,t['q-value']))
p['q']=p.SpecId.map(q)
inpin=p[p.key.isin(dk)]
best=inpin.groupby('key').agg(n=('SpecId','size'),minrank=('rank','min'),q=('q','min'))
print(f'DIA-NN in-slice {len(dk)}; in PIN {len(best)}; rank1 somewhere {(best.minrank==1).sum()}; q<=1% {(best.q<=0.01).sum()}; q<=5% {(best.q<=0.05).sum()}')
print('PSM rows per DIA-NN precursor', best.n.describe()[['25%','50%','75%']].tolist())
