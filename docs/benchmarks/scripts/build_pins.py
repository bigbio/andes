import sys, math, re
import xml.etree.ElementTree as ET

def lnsafe(x):
    try: x=float(x)
    except: return 0.0
    return math.log(x) if x>0 else math.log(1e-300)

def build_java(tsv, out):
    rows=[]
    with open(tsv) as f:
        f.readline()
        for ln in f:
            c=ln.rstrip("\n").split("\t")
            if len(c)<14: continue
            specid=c[1]; scan=c[2]; iso=c[5]; perr=c[6]; charge=int(c[7]); pep=c[8]; prot=c[9]
            denovo=float(c[10]); raw=float(c[11]); spece=c[12]; evalue=c[13]
            isdec = 1 if prot and all(p.startswith("XXX_") for p in prot.split(";") if p) else 0
            label = -1 if isdec else 1
            seq=re.sub(r'[^A-Z]','',pep); plen=len(seq)
            try: aperr=abs(float(perr))
            except: aperr=0.0
            feats=[raw, denovo, raw-denovo, lnsafe(spece), lnsafe(evalue), float(iso), aperr, plen,
                   1 if charge==2 else 0, 1 if charge==3 else 0, 1 if charge>=4 else 0]
            rows.append((specid,label,scan,feats,seq,prot))
    fn=["RawScore","DeNovoScore","ScoreDiff","lnSpecEValue","lnEValue","IsotopeError","absPrecErrPpm","PepLen","Charge2","Charge3","Charge4"]
    with open(out,"w") as o:
        o.write("SpecId\tLabel\tScanNr\t"+"\t".join(fn)+"\tPeptide\tProteins\n")
        for specid,label,scan,feats,seq,prot in rows:
            o.write(f"{specid}\t{label}\t{scan}\t"+"\t".join(f"{x:.6g}" for x in feats)+f"\t-.{seq}.-\t"+"\t".join(prot.split(";"))+"\n")
    print(f"java pin: {len(rows)} rows")

def build_prose(idxml, out):
    NUM=["delta_score","fragment_mz_error_median_ppm","precursor_mz_error_ppm","matched_ion_current",
         "matched_prefix_ions_fraction","matched_suffix_ions_fraction","num_matched_peaks",
         "matched_prefix_ions","matched_suffix_ions","longest_peptide_ion_sequence","isotope_error"]
    rows=[]; n=0; prot_map={}
    for ev,el in ET.iterparse(idxml, events=("end",)):
        tag=el.tag.split('}')[-1]
        if tag=="ProteinHit":
            pid=el.get("id"); acc=el.get("accession")
            if pid and acc: prot_map[pid]=acc
        elif tag=="PeptideIdentification":
            spec=None
            for up in el.findall("./UserParam"):
                if up.get("name")=="spectrum_reference": spec=up.get("value")
            best=None; bestsc=None
            for h in el.findall("./PeptideHit"):
                try: sc=float(h.get("score"))
                except: continue
                if bestsc is None or sc>bestsc: bestsc=sc; best=h
            if best is not None:
                h=best
                seq=re.sub(r'[^A-Z]','',h.get("sequence","")); charge=int(h.get("charge","2"))
                td={up.get("name"):up.get("value") for up in h.findall("./UserParam")}
                accs=[prot_map.get(r,r) for r in h.get("protein_refs","").split() if r] or ["unknown"]
                isdec=1 if td.get("target_decoy","")=="decoy" or all(a.startswith("DECOY_") for a in accs) else 0
                label=-1 if isdec else 1
                feats=[float(h.get("score"))]
                for k in NUM:
                    try: fv=float(td.get(k,"0"))
                    except: fv=0.0
                    if "error" in k: fv=abs(fv)
                    feats.append(fv)
                feats += [len(seq), 1 if charge==2 else 0, 1 if charge==3 else 0, 1 if charge>=4 else 0]
                n+=1
                rows.append((f"prose_{spec}_{n}",label,n,feats,seq,"\t".join(accs)))
            el.clear()
    fn=["hyperscore"]+NUM+["PepLen","Charge2","Charge3","Charge4"]
    with open(out,"w") as o:
        o.write("SpecId\tLabel\tScanNr\t"+"\t".join(fn)+"\tPeptide\tProteins\n")
        for specid,label,scan,feats,seq,prot in rows:
            o.write(f"{specid}\t{label}\t{scan}\t"+"\t".join(f"{x:.6g}" for x in feats)+f"\t-.{seq}.-\t{prot}\n")
    print(f"prose pin: {len(rows)} rows")

if sys.argv[1]=="java": build_java(sys.argv[2], sys.argv[3])
else: build_prose(sys.argv[2], sys.argv[3])
