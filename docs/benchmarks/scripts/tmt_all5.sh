#!/bin/bash
set -uo pipefail
TD=$BENCH/tmt-data
MZML=$TD/a05058.mzML
FASTA=$TD/PXD007683_UP000005640_UP000002311_reviewed.fasta
MODS=$TD/mods-numeric.txt
REPO=$BENCH/repo/msgf-rust
# Binary renamed to andes; fall back to the old cimas/simas names on pre-rename checkouts.
ANDES="${ANDES:-$REPO/target/release/andes}"; [ -x "$ANDES" ] || ANDES="$REPO/target/release/cimas"; [ -x "$ANDES" ] || ANDES="$REPO/target/release/simas"
SIMAS=$ANDES
SMODEL=$BENCH/repo/msgf-rust/resources/ionstat/models.parquet
JAR=$BENCH/MSGFPlus_v20240326.jar
FRAGGER=$BENCH/engines/msfragger-env/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
SAGE=$BENCH/engines/sage-v0.14.7-x86_64-unknown-linux-gnu/sage
SAGECFG=$TD/sage-a05058.json
FRGCFG=$TD/fragger_a05058.params
RES=$BENCH/repo/bench-tmt; mkdir -p $RES
WORK=/tmp/tmt-bench; mkdir -p $WORK
OIMG=ghcr.io/openms/openms-tools-thirdparty:latest
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
export DOTNET_ROOT=/opt/dotnet8; export PATH=/opt/dotnet8:$PATH
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

echo "################ a05058 TMT 5-ENGINE BENCHMARK $(date -Is) ################"
echo "low-res CID TMT: trypsin 2MC len6-40 z2-4 fixC+57 fixK+229.16 fixNterm+229.16 varM+ox prec20ppm frag0.4Da iso-1..2 FDR@1% Percolator(seed42)"

echo "=== [1] simas top-1 (mzML) $(date -Is) ==="
/usr/bin/time -v $SIMAS --spectrum $MZML --database $FASTA --mods $MODS --model-store $SMODEL --fragmentation CID --instrument low-res --protocol TMT --precursor-tol-ppm 20 --isotope-error-min -1 --isotope-error-max 2 --enzyme-specificity fully --max-missed-cleavages 2 --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 --top-n 1 --min-peaks 10 --threads 8 --output-pin $RES/simas_top1.pin > $RES/simas_top1.log 2>&1
echo "  exit=$?"; wallrss $RES/simas_top1.log; perc simas_top1; count simas_top1

echo "=== [2] Java MS-GF+ v20240326 (mzML) $(date -Is) ==="
/usr/bin/time -v java -Xmx14g -jar $JAR -s $MZML -d $FASTA -mod $MODS -o $WORK/java_tmt.mzid -tda 1 -t 20ppm -ti -1,2 -m 1 -inst 0 -e 1 -protocol 4 -ntt 2 -minLength 6 -maxLength 40 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -thread 8 > $RES/java_tmt.log 2>&1
echo "  exit=$?"; wallrss $RES/java_tmt.log
docker run --rm --platform linux/amd64 -v "$WORK":/w -v "$RES":/r $PIMG msgf2pin -P XXX_ -e trypsin -o /r/java_tmt.pin /w/java_tmt.mzid > $RES/java_msgf2pin.log 2>&1
echo "  msgf2pin exit=$? rows=$(($(wc -l < $RES/java_tmt.pin 2>/dev/null || echo 1)-1))"; perc java_tmt; count java_tmt

echo "=== [3] Sage v0.14.7 (mzML) $(date -Is) ==="
sed "s#SAGE_OUT_PLACEHOLDER#$WORK/sage_tmt#g" $SAGECFG > $WORK/sage_tmt.json
/usr/bin/time -v $SAGE $WORK/sage_tmt.json --write-pin -o $WORK/sage_tmt $MZML > $RES/sage_tmt.log 2>&1
echo "  exit=$?"; wallrss $RES/sage_tmt.log
[ -f $WORK/sage_tmt/results.sage.pin ] && cp -f $WORK/sage_tmt/results.sage.pin $RES/sage_tmt.pin
echo "  rows=$(($(wc -l < $RES/sage_tmt.pin 2>/dev/null || echo 1)-1))"; perc sage_tmt; count sage_tmt

echo "=== [4] MSFragger 4.2 (mzML) $(date -Is) ==="
cp -f $MZML $WORK/a05058.mzML
/usr/bin/time -v java -Xmx14g -jar $FRAGGER $FRGCFG $WORK/a05058.mzML > $RES/fragger_tmt.log 2>&1
echo "  exit=$?"; wallrss $RES/fragger_tmt.log
[ -f $WORK/a05058.pin ] && cp -f $WORK/a05058.pin $RES/fragger_tmt.pin
echo "  rows=$(($(wc -l < $RES/fragger_tmt.pin 2>/dev/null || echo 1)-1))"; perc fragger_tmt; count fragger_tmt
rm -f $WORK/a05058.mzML $WORK/a05058.pin $WORK/a05058.pepXML $WORK/a05058.tsv

echo "=== [5] ProSE (mzML, OpenMS) $(date -Is) ==="
mkdir -p $RES/prose-out
/usr/bin/time -v docker run --rm -v $TD:/data:ro -v $RES/prose-out:/out $OIMG bash -c '
/opt/OpenMS/bin/ProSE -in /data/a05058.mzML -database /data/PXD007683_UP000005640_UP000002311_reviewed.fasta -out_idxml /out/prose.idXML -Search:decoys -Search:enzyme Trypsin -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 6 -Search:peptide:max_size 40 -Search:peptide:enzyme_specificity full -Search:precursor:mass_tolerance_lower 20 -Search:precursor:mass_tolerance_upper 20 -Search:precursor:mass_tolerance_unit ppm -Search:precursor:min_charge 2 -Search:precursor:max_charge 4 -Search:precursor:isotope_error_min -1 -Search:precursor:isotope_error_max 2 -Search:fragment:mass_tolerance 0.4 -Search:fragment:mass_tolerance_unit Da -Search:modifications:fixed "Carbamidomethyl (C)" "TMT6plex (K)" "TMT6plex (N-term)" -Search:modifications:variable "Oxidation (M)" -Search:modifications:variable_max_per_peptide 3 -threads 8 2>&1 | tail -4
echo ---TextExporter---; /opt/OpenMS/bin/TextExporter -in /out/prose.idXML -out /out/prose_raw.tsv 2>&1 | tail -2
' > $RES/prose.log 2>&1
echo "  ProSE search exit=$?"; wallrss $RES/prose.log
python3 - "$RES/prose-out/prose_raw.tsv" <<'PY'
import sys
rows=[]
for ln in open(sys.argv[1]):
    if not ln.startswith("PEPTIDE\t"): continue
    f=ln.rstrip("\n").split("\t")
    if f[4]!="0": continue
    accs=[a for a in f[11].split(";") if a]
    isdec=1 if accs and all(a.startswith("DECOY_") for a in accs) else 0
    rows.append((float(f[3]),isdec,f[5],accs))
rows.sort(key=lambda r:-r[0]); d=t=0; fdr=[]
for s,isdec,seq,accs in rows:
    if isdec: d+=1
    else: t+=1
    fdr.append(d/max(t,1))
q=[0]*len(rows); m=1.0
for i in range(len(rows)-1,-1,-1): m=min(m,fdr[i]); q[i]=m
acc=0; peps=set(); prots=set()
for (s,isdec,seq,accs),qi in zip(rows,q):
    if isdec==0 and qi<=0.01:
        acc+=1; peps.add(seq)
        for a in accs:
            if not a.startswith("DECOY_"): prots.add(a)
print(f"  RESULT prose_tmt             PSMs@1%={acc} peptides@1%={len(peps)} proteins@1%={len(prots)} (native TDC, DECOY_)")
PY
echo "################ TMT_5ENGINE_DONE $(date -Is) ################"
