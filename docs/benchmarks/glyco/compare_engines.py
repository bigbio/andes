#!/usr/bin/env python3
"""Cross-engine glycopeptide set comparison (PXD031032 / PXD005553, MouseLiver-Z-T-1).

Normalises five engines to a canonical glycan composition `(HexNAc, Hex, Fuc,
NeuAc, NeuGc)` and a bare peptide sequence (flanking residues, modification
masses and the glycan tag stripped), then reports (a) the distinct
(peptide, composition) count per engine and (b) the pairwise overlap / Jaccard.

Inputs (see CONFIG below) are the published result files, none of which is
produced by andes except the andes PIN + native-rescore PSMs:
  * Glyco-Decipher  *_Glyco_Decipher_GPSM.txt   (PeptideFDR <= 0.01)
  * StrucGP         *_result_StrucGP.xlsx        (all rows)
  * Byonic          *.raw_*_Byonic.xlsx          (Score >= 300)
  * pGlyco 2.0      *-FDR.txt                    (TotalFDR <= 0.01, target only)
  * andes           <stem>.glyco.pin + native-rescore .psms (q <= 0.01)
"""
import csv, re, json, collections
from xml.etree import ElementTree as ET
M='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

# ---------- input files (edit for a different run) ----------
CONFIG = {
    'Glyco-Decipher': '/tmp/pxd031032_ref/MouseLiver-Z-T-1_Glyco_Decipher_GPSM.txt',
    'StrucGP':        '/tmp/pxd031032_ref/MouseLiver-Z-T-1_result_StrucGP.xlsx',
    'Byonic':         '/tmp/pxd031032_ref/MouseLiver-Z-T-1.raw_20210216_Byonic.xlsx',
    'pGlyco2':        '/tmp/pxd031032_ref/MouseLiver-Z-T-1-FDR.txt',
    'andes_pin':      '/tmp/andes_mouse/z1_full_gdb.glyco.pin',
    'andes_psms':     '/tmp/andes_mouse/z1_full_gdb.psms',
    'out':            '/tmp/final_comparison.json',
}

# ---------- glycan comp normalization: all -> (HexNAc,Hex,Fuc,NeuAc,NeuGc) ----------
def parse_comp_simple(s):
    d={'Hex':0,'HexNAc':0,'Fuc':0,'NeuAc':0,'NeuGc':0,'Sia':0,'Neu5Ac':0,'Neu5Gc':0}
    for m in re.finditer(r'([A-Za-z0-9]+)\((\d+)\)', s):
        k,v=m.group(1),int(m.group(2))
        if k in d: d[k]+=v
    d['NeuAc']+=d.get('Neu5Ac',0)+d.get('Sia',0); d['NeuGc']+=d.get('Neu5Gc',0)
    return (d['HexNAc'],d['Hex'],d['Fuc'],d['NeuAc'],d['NeuGc'])

def parse_comp_strucgp(s):
    d={'N':0,'H':0,'F':0,'A':0,'G':0,'S':0}
    for m in re.finditer(r'([NHFAGS])(\d+)', s): d[m.group(1)]+=int(m.group(2))
    return (d['N'],d['H'],d['F'],d['A']+d['S'],d['G'])

def parse_comp_andes(s):
    # "HexNAc2Hex8Fuc0NeuAc0NeuGc0"
    d={'HexNAc':0,'Hex':0,'Fuc':0,'NeuAc':0,'NeuGc':0}
    for m in re.finditer(r'(HexNAc|Hex|Fuc|NeuAc|NeuGc)(\d+)', s):
        d[m.group(1)]+=int(m.group(2))
    return (d['HexNAc'],d['Hex'],d['Fuc'],d['NeuAc'],d['NeuGc'])

def parse_comp_pglyco(s):
    # pGlyco2 "Hex HexNAc Fuc NeuAc NeuGc" e.g. "6 2 0 0 0" (order confirmed by GlyMass)
    t=s.split()
    t=t+[0]*5
    hex_,hn,f,a,g=[int(float(x)) for x in t[:5]]
    return (hn,hex_,f,a,g)

def acc(prot):
    m=re.search(r'\|([A-Za-z0-9_]+)\|', prot)
    return m.group(1) if m else prot

def acc_strucgp(prot): return prot  # already accession

def strip_byonic_pep(p):
    # "K.HLLEN[+1864.634]ATASVSEAER.K" -> strip the single-AA flanks by first/last
    # dot, NOT a full split (the glycan/mod masses contain decimal points).
    p=p.strip()
    if '.' in p:
        p=p.split('.', 1)[1]    # drop leading "X."
        p=p.rsplit('.', 1)[0]   # drop trailing ".X"
    return re.sub(r'\[[^\[\]]*\]','',p)

def load_sheet(path, sheetfile):
    z=__import__('zipfile').ZipFile(path)
    sst=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall(M+'si'):
            sst.append(''.join(t.text or '' for t in si.iter(M+'t')))
    root=ET.fromstring(z.read(sheetfile)); rows=[]
    for row in root.findall('.//'+M+'row'):
        cells={}
        for c in row.findall(M+'c'):
            ref=c.get('r') or ''; t=c.get('t'); val=''
            v=c.find(M+'v')
            if v is not None and v.text is not None: val=sst[int(v.text)] if t=='s' else v.text
            else:
                isn=c.find(M+'is')
                if isn is not None: val=''.join(x.text or '' for x in isn.iter(M+'t'))
            n=0
            for ch in ''.join(x for x in ref if x.isalpha()): n=n*26+(ord(ch.upper())-64)
            cells[n]=val
        rows.append(cells)
    return rows

ref={}

# 1. Glyco-Decipher (PeptideFDR<=0.01)
gd=[]
with open(CONFIG['Glyco-Decipher']) as f:
    r=csv.reader(f,delimiter='\t'); next(r)
    for line in r:
        if len(line)<19 or float(line[11])>0.01: continue
        gd.append((line[5], parse_comp_simple(line[17]), acc(line[8])))
ref['Glyco-Decipher']={'gps':set((p,g) for p,g,_ in gd),'prots':set(x[2] for x in gd),'n':len(gd)}

# 2. StrucGP
sg=load_sheet(CONFIG['StrucGP'],'xl/worksheets/sheet1.xml')
sgd=[]
for r in sg[1:]:
    p=r.get(12,'').strip(); pr=r.get(15,'').strip(); gc=r.get(20,'').strip()
    if not p: continue
    sgd.append((p, parse_comp_strucgp(gc), pr))
ref['StrucGP']={'gps':set((p,g) for p,g,_ in sgd),'prots':set(x[2] for x in sgd),'n':len(sgd)}

# 3. Byonic (Score>=300)
by=load_sheet(CONFIG['Byonic'],'xl/worksheets/sheet2.xml')
byd=[]
for r in by[1:]:
    gly=r.get(4,'').strip(); pep=r.get(3,'').strip(); pr=r.get(19,'').strip()
    if not gly: continue
    if float(r.get(14,0) or 0)<300: continue
    byd.append((strip_byonic_pep(pep), parse_comp_simple(gly), acc(pr)))
ref['Byonic']={'gps':set((p,g) for p,g,_ in byd),'prots':set(x[2] for x in byd),'n':len(byd)}

# 3b. pGlyco 2.0 (PXD005553, same raw files; TotalFDR<=0.01, target only)
pg=[]
with open(CONFIG['pGlyco2']) as f:
    r=csv.reader(f,delimiter='\t'); next(r)
    for line in r:
        if len(line)<25: continue
        if line[19].strip()!='0' or line[20].strip()!='0': continue  # GlyDecoy/PepDecoy
        if float(line[23])>0.01: continue  # TotalFDR
        pep=line[5].replace('J','N')          # J = glycosylated Asn
        g=parse_comp_pglyco(line[8])
        pg.append((pep, g, acc(line[24])))
ref['pGlyco2']={'gps':set((p,g) for p,g,_ in pg),'prots':set(x[2] for x in pg),'n':len(pg)}

# 4. andes: read pin + psms q-values
def load_andes(pin_path, psms_path, qmax=0.01):
    # q-values by SpecId
    q={}
    with open(psms_path) as f:
        r=csv.reader(f,delimiter='\t'); next(r)
        for line in r:
            if len(line)<3: continue
            if float(line[2])<=qmax: q[line[0]]=float(line[2])
    # pin rows
    andes=[]
    n_entrap=0
    with open(pin_path) as f:
        r=csv.reader(f,delimiter='\t'); hdr=next(r)
        ci={h:i for i,h in enumerate(hdr)}
        for line in r:
            if len(line)<76: continue
            spec=line[ci['SpecId']]; label=line[ci['Label']]
            if label!='1': continue
            if spec not in q: continue
            pep=line[ci['Peptide']]; prot=line[ci['Proteins']]
            # parse "K.TAANGTR.K[HexNAc2Hex8Fuc0NeuAc0NeuGc0@N4]"
            m=re.match(r'^(.*?)\[([^\]]*)\]$', pep)
            if not m: continue
            pseq, gann = m.group(1), m.group(2)
            # strip flanks "K.X.K" -> "X" (mods like +15.99491 also contain dots)
            if '.' in pseq:
                pseq=pseq.split('.', 1)[1]
                pseq=pseq.rsplit('.', 1)[0]
            pseq=re.sub(r'\+[\d.]+','',pseq)  # strip mods like +15.99491
            # glycan comp
            gm=re.match(r'(HexNAc\d+Hex\d+Fuc\d+NeuAc\d+NeuGc\d+)', gann)
            if not gm: continue
            gc=parse_comp_andes(gm.group(1))
            a=acc(prot)
            is_entrap = a.startswith('ENTRAP_') or 'ENTRAP' in a
            if is_entrap: n_entrap+=1
            andes.append((pseq, gc, a))
    return andes, n_entrap

andes_rows, n_entrap = load_andes(CONFIG['andes_pin'],CONFIG['andes_psms'],0.01)
andes_gps=set((p,g) for p,g,_ in andes_rows)
andes_prots=set(a for _,_,a in andes_rows)
andes_prots_noentrap=set(a for _,_,a in andes_rows if not a.startswith('ENTRAP_'))

print("="*72)
print("MouseLiver-Z-T-1  |  unique intact glycopeptides (peptide + glycan composition)")
print("="*72)
allnames=['Glyco-Decipher','StrucGP','Byonic','pGlyco2','andes']
for nm in ['Glyco-Decipher','StrucGP','Byonic','pGlyco2']:
    print(f"  {nm:16s}  {len(ref[nm]['gps']):5d}  glycopeptides   (PSM rows: {ref[nm]['n']})   proteins: {len(ref[nm]['prots'])}")
print(f"  {'andes (1% FDR)':16s}  {len(andes_gps):5d}  glycopeptides   (PSM rows: {len(andes_rows)})   proteins: {len(andes_prots)} (excl ENTRAP: {len(andes_prots_noentrap)})")
print(f"\n  [andes entrapment-mapped rows at 1% FDR: {n_entrap}]")

# overlap (peptide,glycan)
print("\n=== pairwise overlap on (peptide, glycan comp) ===")
sets={'Glyco-Decipher':ref['Glyco-Decipher']['gps'],'StrucGP':ref['StrucGP']['gps'],'Byonic':ref['Byonic']['gps'],'pGlyco2':ref['pGlyco2']['gps'],'andes':andes_gps}
for i in range(len(allnames)):
    for j in range(i+1,len(allnames)):
        a=allnames[i]; b=allnames[j]
        inter=len(sets[a]&sets[b]); union=len(sets[a]|sets[b])
        print(f"  {a:16s} vs {b:8s}: {inter:5d} shared  (Jaccard {inter/union:.3f})")

# protein overlap
print("\n=== protein overlap ===")
psets={'Glyco-Decipher':ref['Glyco-Decipher']['prots'],'StrucGP':ref['StrucGP']['prots'],'Byonic':ref['Byonic']['prots'],'pGlyco2':ref['pGlyco2']['prots'],'andes':andes_prots_noentrap}
for i in range(len(allnames)):
    for j in range(i+1,len(allnames)):
        a=allnames[i]; b=allnames[j]
        print(f"  {a:16s} vs {b:8s}: {len(psets[a]&psets[b]):5d} shared proteins")

json.dump({'andes_gps':sorted(andes_gps),'andes_prots':sorted(andes_prots_noentrap),
           'ref':{k:{'gps':sorted(v['gps']),'prots':sorted(v['prots']),'n':v['n']} for k,v in ref.items()}},
          open(CONFIG['out'],'w'))
print(f"\n[saved {CONFIG['out']}]")
