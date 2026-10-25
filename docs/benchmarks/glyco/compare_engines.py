#!/usr/bin/env python3
"""Cross-engine glycopeptide set comparison (PXD031032 / PXD005553, MouseLiver-Z-T-1).

Normalises six engines to a canonical glycan composition `(HexNAc, Hex, Fuc,
NeuAc, NeuGc)` and a bare peptide sequence (flanking residues, modification
masses and the glycan tag stripped), then reports (a) the distinct
(peptide, composition) count per engine and (b) the pairwise overlap / Jaccard.

Inputs (command-line paths) are the published result files, none of which is
produced by andes except the andes PIN + native-rescore PSMs:
  * Glyco-Decipher   GPSM text export                 (PeptideFDR <= 0.01)
  * StrucGP          result workbook (.xlsx)          (all rows)
  * Byonic           result workbook (.xlsx)          (Score >= 300)
  * pGlyco2          *-FDR.txt, PXD005553             (TotalFDR <= 0.01, target only)
  * MSFragger-Glyco  psm.tsv, PXD031032               (Expectation <= 0.01, Delta
                     Mass -> gdb composition within 0.02 Da)
  * andes            <stem>.glyco.pin + native-rescore .psms (q <= 0.01)

An andes row is entrapment only if every protein it maps to is an ENTRAP_ protein. Pass
`--fasta` (the search database) to also rescue entrapment-only rows whose I/L-folded bare
peptide occurs in a sample protein (isobaric, therefore not really entrapment). Both rules
are shared with eval_entrap.py and reproduce/run.sh via ../entrap_common.py.
"""
import argparse, csv, re, json, collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from entrap_common import SampleSequences, bare_peptide, classify  # noqa: E402
from xml.etree import ElementTree as ET
M='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

# ---------- input files ----------
_ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
_ap.add_argument('--glyco-decipher', required=True, help='Glyco-Decipher GPSM text export')
_ap.add_argument('--strucgp', required=True, help='StrucGP result workbook (.xlsx)')
_ap.add_argument('--byonic', required=True, help='Byonic result workbook (.xlsx)')
_ap.add_argument('--pglyco2', required=True, help='pGlyco2 *-FDR.txt')
_ap.add_argument('--msfragger-glyco', required=True, help='MSFragger-Glyco psm.tsv')
_ap.add_argument('--gdb', required=True, help='glycan .gdb used to map MSFragger-Glyco delta masses')
_ap.add_argument('--andes-pin', required=True)
_ap.add_argument('--andes-psms', required=True)
_ap.add_argument('--fasta', default=None, help='search FASTA; enables I/L-equivalent entrapment rescue')
_ap.add_argument('--out', default='final_comparison.json')
_args = _ap.parse_args()
CONFIG = {
    'Glyco-Decipher':  _args.glyco_decipher,
    'StrucGP':         _args.strucgp,
    'Byonic':          _args.byonic,
    'pGlyco2':         _args.pglyco2,
    'MSFragger-Glyco': _args.msfragger_glyco,
    'gdb':             _args.gdb,
    'andes_pin':       _args.andes_pin,
    'andes_psms':      _args.andes_psms,
    'fasta':           _args.fasta,
    'out':             _args.out,
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

def parse_comp_pglyco2(s):
    # pGlyco2 "Hex HexNAc Fuc NeuAc NeuGc" e.g. "6 2 0 0 0" (order confirmed by GlyMass)
    t=s.split()
    t=t+[0]*5
    hex_,hn,f,a,g=[int(float(x)) for x in t[:5]]
    return (hn,hex_,f,a,g)

def acc(prot):
    m=re.search(r'\|([^|]+)\|', prot)
    a=m.group(1) if m else prot
    return re.sub(r'-\d+$', '', a)  # UniProt isoform -> canonical accession

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

# 3b. pGlyco2 (PXD005553, same raw files; TotalFDR<=0.01, target only)
pg=[]
with open(CONFIG['pGlyco2']) as f:
    r=csv.reader(f,delimiter='\t'); next(r)
    for line in r:
        if len(line)<25: continue
        if line[19].strip()!='0' or line[20].strip()!='0': continue  # GlyDecoy/PepDecoy
        if float(line[23])>0.01: continue  # TotalFDR
        pep=line[5].replace('J','N')          # J = glycosylated Asn
        g=parse_comp_pglyco2(line[8])
        pg.append((pep, g, acc(line[24])))
ref['pGlyco2']={'gps':set((p,g) for p,g,_ in pg),'prots':set(x[2] for x in pg),'n':len(pg)}

# 4. andes: read pin + psms q-values
def load_andes(pin_path, psms_path, qmax=0.01, sample=None):
    """Accepted andes glyco rows as (peptide, composition, accession, is_entrapment).

    A row is entrapment by the shared rule in ../entrap_common.py: every protein column is an
    ENTRAP_ protein and, when `sample` (a SampleSequences) is given, its I/L-folded bare
    peptide does not occur in a sample protein.
    """
    # q-values by SpecId
    q={}
    with open(psms_path) as f:
        r=csv.reader(f,delimiter='\t'); qi=next(r).index('q-value')  # a `filename` column may precede it
        for line in r:
            if len(line)<=qi: continue
            if float(line[qi])<=qmax: q[line[0]]=float(line[qi])
    # pin rows
    andes=[]
    n_entrap=0
    with open(pin_path) as f:
        r=csv.reader(f,delimiter='\t'); hdr=next(r)
        ci={h:i for i,h in enumerate(hdr)}
        i_pep=ci['Peptide']
        for line in r:
            if len(line)<=i_pep+1: continue  # no protein column
            spec=line[ci['SpecId']]; label=line[ci['Label']]
            if label!='1': continue
            if spec not in q: continue
            pep=line[i_pep]; prots=line[i_pep+1:]
            # "K.TAANGTR.K[HexNAc2Hex8Fuc0NeuAc0NeuGc0@N4]": the glycan tag is the trailing
            # bracket after the C-terminal flank.
            m=re.search(r'\.[^.]\[([^\[\]]*)\]$', pep)
            if not m: continue
            gm=re.match(r'(HexNAc\d+Hex\d+Fuc\d+NeuAc\d+NeuGc\d+)', m.group(1))
            if not gm: continue
            gc=parse_comp_andes(gm.group(1))
            pseq=bare_peptide(pep)
            a=acc(prots[0])
            is_entrap = classify(prots, pep, sample) == 'entrap'
            if is_entrap: n_entrap+=1
            andes.append((pseq, gc, a, is_entrap))
    return andes, n_entrap

# 5. MSFragger-Glyco: the export carries the glycan ONLY as a "Delta Mass"
#    (residue-sum mass, no water subtracted), so we map it back to a canonical
#    composition via the --gdb composition space. Expectation is
#    the FDR proxy (the file has no glycan q-value / no decoys).
def _gdb_comp_masses(gdb):
    M={'H':162.052824,'N':203.079373,'F':146.057909,'A':291.095417,'G':307.090607}
    comps=set()
    with open(gdb) as f:
        f.readline()
        for line in f:
            c=collections.Counter(ch for ch in line.strip() if ch in 'HNFAG')
            if c: comps.add((c['N'],c['H'],c['F'],c['A'],c['G']))
    items=sorted((round(M['N']*t[0]+M['H']*t[1]+M['F']*t[2]+M['A']*t[3]+M['G']*t[4],5),t) for t in comps)
    return [m for m,_ in items], items

def load_msfragger_glyco(tsv, gdb, tissue='MouseLiver-Z-T-1', emax=1e-2, tol=0.02):
    import bisect
    masses, items = _gdb_comp_masses(gdb)
    C13=1.003355  # one ^13C: MSFragger-Glyco reports M+1/M+2 precursors as their own rows
    def lookup(dm):
        lo=bisect.bisect_left(masses, dm-tol); hi=bisect.bisect_right(masses, dm+tol)
        cand=items[lo:hi]
        if not cand: return None
        cand.sort(key=lambda x: abs(x[0]-dm))
        return cand[0][1]
    rows=[]; prots=set()
    with open(tsv, encoding='utf-8', errors='replace') as f:
        r=csv.reader(f, delimiter='\t'); hdr=next(r)
        ci={h:i for i,h in enumerate(hdr)}
        si=ci['Spectrum']; pi=ci['Peptide']; di=ci['Delta Mass']; ei=ci['Expectation']; pri=ci['Protein']
        for row in r:
            if len(row)<=max(si,pi,di,ei,pri): continue
            if not row[si].startswith(tissue): continue
            try:
                if float(row[ei])>emax: continue
                dm=float(row[di])
            except (ValueError,TypeError): continue
            if abs(dm)<400: continue          # non-glycosylated peptide
            comp=None
            for shift in (0.0,-C13,-2*C13,C13):
                comp=lookup(dm+shift)
                if comp is not None: break
            if comp is None: continue
            rows.append((row[pi].strip(), comp, acc(row[pri]))); prots.add(acc(row[pri]))
    return rows, prots

_sample = SampleSequences(CONFIG['fasta']) if CONFIG['fasta'] else None
andes_rows, n_entrap = load_andes(CONFIG['andes_pin'],CONFIG['andes_psms'],0.01,_sample)
ms_rows, ms_prots = load_msfragger_glyco(CONFIG['MSFragger-Glyco'], CONFIG['gdb'])
ref['MSFragger-Glyco']={'gps':set((p,g) for p,g,_ in ms_rows),'prots':ms_prots,'n':len(ms_rows)}
# Entrapment hits are known-false by construction and no other engine searched an
# entrapment database, so they must not inflate the headline count.
andes_gps=set((p,g) for p,g,a,e in andes_rows if not e)
andes_gps_with_entrap=set((p,g) for p,g,a,e in andes_rows)
andes_prots=set(a for _,_,a,_ in andes_rows)
andes_prots_noentrap=set(a for _,_,a,e in andes_rows if not e)

print("="*72)
print("MouseLiver-Z-T-1  |  unique intact glycopeptides (peptide + glycan composition)")
print("="*72)
allnames=['Glyco-Decipher','StrucGP','Byonic','pGlyco2','MSFragger-Glyco','andes']
for nm in ['Glyco-Decipher','StrucGP','Byonic','pGlyco2','MSFragger-Glyco']:
    print(f"  {nm:16s}  {len(ref[nm]['gps']):5d}  glycopeptides   (PSM rows: {ref[nm]['n']})   proteins: {len(ref[nm]['prots'])}")
print(f"  {'andes (1% FDR)':16s}  {len(andes_gps):5d}  glycopeptides   (PSM rows: {len(andes_rows)-n_entrap})   proteins: {len(andes_prots_noentrap)}")
print(f"  {'':16s}  {len(andes_gps_with_entrap):5d}  including entrapment   (PSM rows: {len(andes_rows)})   proteins: {len(andes_prots)})")
print(f"\n  [andes entrapment-mapped rows at 1% FDR: {n_entrap}]")

# overlap (peptide,glycan)
print("\n=== pairwise overlap on (peptide, glycan comp) ===")
sets={'Glyco-Decipher':ref['Glyco-Decipher']['gps'],'StrucGP':ref['StrucGP']['gps'],'Byonic':ref['Byonic']['gps'],'pGlyco2':ref['pGlyco2']['gps'],'MSFragger-Glyco':ref['MSFragger-Glyco']['gps'],'andes':andes_gps}
for i in range(len(allnames)):
    for j in range(i+1,len(allnames)):
        a=allnames[i]; b=allnames[j]
        inter=len(sets[a]&sets[b]); union=len(sets[a]|sets[b])
        print(f"  {a:16s} vs {b:8s}: {inter:5d} shared  (Jaccard {inter/union:.3f})")

# protein overlap
print("\n=== protein overlap ===")
psets={'Glyco-Decipher':ref['Glyco-Decipher']['prots'],'StrucGP':ref['StrucGP']['prots'],'Byonic':ref['Byonic']['prots'],'pGlyco2':ref['pGlyco2']['prots'],'MSFragger-Glyco':ref['MSFragger-Glyco']['prots'],'andes':andes_prots_noentrap}
for i in range(len(allnames)):
    for j in range(i+1,len(allnames)):
        a=allnames[i]; b=allnames[j]
        print(f"  {a:16s} vs {b:8s}: {len(psets[a]&psets[b]):5d} shared proteins")

json.dump({'andes_gps':sorted(andes_gps),'andes_prots':sorted(andes_prots_noentrap),
           'ref':{k:{'gps':sorted(v['gps']),'prots':sorted(v['prots']),'n':v['n']} for k,v in ref.items()}},
          open(CONFIG['out'],'w'))
print(f"\n[saved {CONFIG['out']}]")
