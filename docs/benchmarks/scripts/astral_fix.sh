#!/bin/bash
set -uo pipefail
AD=/srv/data/msgf-bench/astral-data
RAW=$AD/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw
MZML=$AD/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
FASTA=$AD/ProteoBenchFASTA_MixedSpecies_HYE.fasta
FRAGGER=/srv/data/msgf-bench/engines/msfragger-env/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
SAGE=/srv/data/msgf-bench/engines/sage-v0.14.7-x86_64-unknown-linux-gnu/sage
SAGECFG=/srv/data/msgf-bench/repo/msgf-rust/docs/benchmarks/configs/sage-astral.json
FRGCFG=/srv/data/msgf-bench/engines/msfragger-cfg/astral_fragger.params
RES=/srv/data/msgf-bench/repo/bench-astral; mkdir -p $RES
WORK=/tmp/astral-fix; mkdir -p $WORK
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
perc(){ docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y --results-psms /r/$1.t.psms --decoy-results-psms /r/$1.d.psms --only-psms=false /r/$1.pin > $RES/$1.perc.log 2>&1; echo "    perc $1 exit=$?"; }
count(){ tp=$RES/$1.t.psms
  q=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")print i}' "$tp")
  rc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="proteinIds")print i}' "$tp")
  pc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="peptide")print i}' "$tp")
  ps=$(awk -F"\t" -v q="$q" 'NR>1&&$q<=0.01{c++}END{print c+0}' "$tp")
  pep=$(awk -F"\t" -v q="$q" -v p="$pc" 'NR>1&&$q<=0.01{s=$p;gsub(/^[A-Z-]\./,"",s);gsub(/\.[A-Z-]$/,"",s);gsub(/\[[^]]*\]/,"",s);print s}' "$tp"|sort -u|wc -l)
  pr=$(awk -F"\t" -v q="$q" -v r="$rc" 'NR>1&&$q<=0.01{print $r}' "$tp"|tr "\t" "\n"|grep -v "^XXX_\|^rev_\|^DECOY\|^$"|sort -u|wc -l)
  printf "  RESULT %-22s PSMs@1%%=%s peptides@1%%=%s proteins@1%%=%s\n" "$1" "$ps" "$pep" "$pr"; }
wallrss(){ grep -E "Elapsed \(wall|Maximum resident" $1 | sed 's/^[[:space:]]*/    /'; }

echo "################ ASTRAL FIX (Sage+MSFragger+ProSE) $(date -Is) ################"

echo "=== Sage v0.14.7 chimera=true (native .raw) $(date -Is) ==="
sed -E "s#\"fasta\": *\"[^\"]*\"#\"fasta\": \"$FASTA\"#; s#\"output_directory\": *\"[^\"]*\"#\"output_directory\": \"$WORK/sage_astral\"#" $SAGECFG > $WORK/sage_astral.json
/usr/bin/time -v $SAGE $WORK/sage_astral.json --write-pin -o $WORK/sage_astral $RAW > $RES/sage_astral.log 2>&1
echo "  exit=$?"; wallrss $RES/sage_astral.log
[ -f $WORK/sage_astral/results.sage.pin ] && cp -f $WORK/sage_astral/results.sage.pin $RES/sage_astral.pin
echo "  rows=$(($(wc -l < $RES/sage_astral.pin 2>/dev/null || echo 1)-1))"; perc sage_astral; count sage_astral

echo "=== MSFragger 4.2 DDA+ (mzML) $(date -Is) ==="
/usr/bin/time -v java -Xmx14g -jar $FRAGGER $FRGCFG $MZML > $RES/fragger_astral.log 2>&1
echo "  exit=$?"; wallrss $RES/fragger_astral.log
PINOUT="${MZML%.mzML}.pin"
[ -f "$PINOUT" ] && cp -f "$PINOUT" $RES/fragger_astral.pin
echo "  rows=$(($(wc -l < $RES/fragger_astral.pin 2>/dev/null || echo 1)-1))"; perc fragger_astral; count fragger_astral
rm -f "${MZML%.mzML}.pin" "${MZML%.mzML}.pepXML" "${MZML%.mzML}.tsv"

echo "=== ProSE native-TDC (from existing search) $(date -Is) ==="
T=$RES/prose-out/prose_raw.tsv
python3 - "$T" <<'PY'
import sys
rows=[]
for ln in open(sys.argv[1]):
    if not ln.startswith("PEPTIDE\t"): continue
    f=ln.rstrip("\n").split("\t")
    score=float(f[3]); rank=f[4]; seq=f[5]; acc=f[11]
    if rank!="0": continue
    accs=[a for a in acc.split(";") if a]
    isdec=1 if accs and all(a.startswith("DECOY_") for a in accs) else 0
    rows.append((score,isdec,seq))
rows.sort(key=lambda r:-r[0])
d=t=0; fdrs=[]
for s,isdec,seq in rows:
    if isdec: d+=1
    else: t+=1
    fdrs.append(d/max(t,1))
# cumulative min from bottom -> q-value
q=[0]*len(rows); m=1.0
for i in range(len(rows)-1,-1,-1):
    m=min(m,fdrs[i]); q[i]=m
acc=0; peps=set()
for (s,isdec,seq),qi in zip(rows,q):
    if isdec==0 and qi<=0.01:
        acc+=1; peps.add(seq)
print(f"  RESULT prose_astral           PSMs@1%={acc} peptides@1%={len(peps)} proteins@1%=NA (native TDC, DECOY_)")
PY
echo "################ ASTRAL_FIX_DONE $(date -Is) ################"
