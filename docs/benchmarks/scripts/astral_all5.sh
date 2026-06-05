#!/bin/bash
set -uo pipefail
AD=/srv/data/msgf-bench/astral-data
RAW=$AD/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw
MZML=$AD/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
FASTA=$AD/ProteoBenchFASTA_MixedSpecies_HYE.fasta
MODS=/srv/data/msgf-bench/astral_mods_rust.txt
REPO=/srv/data/msgf-bench/repo/msgf-rust
# Binary renamed simas -> cimas; fall back to the old name on pre-rename checkouts.
SIMAS=$REPO/target/release/cimas; [ -x "$SIMAS" ] || SIMAS=$REPO/target/release/simas
SMODEL=/srv/data/msgf-bench/repo/msgf-rust/resources/ionstat/models.parquet
JAR=/srv/data/msgf-bench/MSGFPlus_v20240326.jar
FRAGGER=/srv/data/msgf-bench/engines/msfragger-env/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
SAGE=/srv/data/msgf-bench/engines/sage-v0.14.7-x86_64-unknown-linux-gnu/sage
SAGECFG=/srv/data/msgf-bench/repo/msgf-rust/docs/benchmarks/configs/sage-astral.json
FRGCFG=/srv/data/msgf-bench/engines/msfragger-cfg/astral_fragger.params
RES=/srv/data/msgf-bench/repo/bench-astral; mkdir -p $RES
WORK=/tmp/astral-bench; mkdir -p $WORK
OIMG=ghcr.io/openms/openms-tools-thirdparty:latest
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
export DOTNET_ROOT=/opt/dotnet8; export PATH=/opt/dotnet8:$PATH
perc(){ docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG \
  percolator --seed 42 -Y --results-psms /r/$1.t.psms --decoy-results-psms /r/$1.d.psms --only-psms=false /r/$1.pin > $RES/$1.perc.log 2>&1; echo "    perc $1 exit=$?"; }
count(){ tp=$RES/$1.t.psms
  q=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")print i}' "$tp")
  rc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="proteinIds")print i}' "$tp")
  pc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="peptide")print i}' "$tp")
  ps=$(awk -F"\t" -v q="$q" 'NR>1&&$q<=0.01{c++}END{print c+0}' "$tp")
  pep=$(awk -F"\t" -v q="$q" -v p="$pc" 'NR>1&&$q<=0.01{s=$p;gsub(/^[A-Z-]\./,"",s);gsub(/\.[A-Z-]$/,"",s);gsub(/\[[^]]*\]/,"",s);print s}' "$tp"|sort -u|wc -l)
  pr=$(awk -F"\t" -v q="$q" -v r="$rc" 'NR>1&&$q<=0.01{print $r}' "$tp"|tr "\t" "\n"|grep -v "^XXX_\|^rev_\|^DECOY\|^$"|sort -u|wc -l)
  printf "  RESULT %-22s PSMs@1%%=%s peptides@1%%=%s proteins@1%%=%s\n" "$1" "$ps" "$pep" "$pr"; }
wallrss(){ grep -E "Elapsed \(wall|Maximum resident" $1 | sed 's/^[[:space:]]*/    /'; }

echo "################ ASTRAL 5-ENGINE BENCHMARK $(date -Is) ################"
echo "trypsin 2MC len7-40 z2-4 fixC+57 varM+ox varProtNterm+acetyl prec10ppm iso-1..2 FDR@1% Percolator(seed42)"

echo "=== [1] simas top-1 (native .raw) $(date -Is) ==="
/usr/bin/time -v $SIMAS --spectrum $RAW --database $FASTA --mods $MODS --model-store $SMODEL --fragmentation auto --precursor-tol-ppm 10 --isotope-error-min -1 --isotope-error-max 2 --enzyme-specificity fully --max-missed-cleavages 2 --min-length 7 --max-length 40 --charge-min 2 --charge-max 4 --top-n 1 --min-peaks 10 --threads 8 --output-pin $RES/simas_top1.pin > $RES/simas_top1.log 2>&1
echo "  exit=$?"; wallrss $RES/simas_top1.log; perc simas_top1; count simas_top1

echo "=== [2] simas --chimeric (native .raw) $(date -Is) ==="
/usr/bin/time -v $SIMAS --spectrum $RAW --database $FASTA --mods $MODS --model-store $SMODEL --fragmentation auto --precursor-tol-ppm 10 --isotope-error-min -1 --isotope-error-max 2 --enzyme-specificity fully --max-missed-cleavages 2 --min-length 7 --max-length 40 --charge-min 2 --charge-max 4 --top-n 1 --min-peaks 10 --threads 8 --chimeric --output-pin $RES/simas_chim.pin > $RES/simas_chim.log 2>&1
echo "  exit=$?"; wallrss $RES/simas_chim.log; perc simas_chim; count simas_chim

echo "=== [3] Java MS-GF+ v20240326 (mzML) $(date -Is) ==="
/usr/bin/time -v java -Xmx14g -jar $JAR -s $MZML -d $FASTA -mod $MODS -o $WORK/java_astral.mzid -tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 7 -maxLength 40 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -thread 8 > $RES/java_astral.log 2>&1
echo "  exit=$?"; wallrss $RES/java_astral.log
docker run --rm --platform linux/amd64 -v "$WORK":/w -v "$RES":/r $PIMG msgf2pin -P XXX_ -e trypsin -o /r/java_astral.pin /w/java_astral.mzid > $RES/java_msgf2pin.log 2>&1
echo "  msgf2pin exit=$? rows=$(($(wc -l < $RES/java_astral.pin 2>/dev/null || echo 1)-1))"; perc java_astral; count java_astral

echo "=== [4] Sage v0.14.7 chimera=true (mzML) $(date -Is) ==="
# Sage's generic release lacks the proprietary Thermo reader, so it runs on the mzML.
sed "s#\"fasta\":.*#\"fasta\": \"$FASTA\",#; s#\"output_directory\":.*#\"output_directory\": \"$WORK/sage_astral\"#" $SAGECFG > $WORK/sage_astral.json
/usr/bin/time -v $SAGE $WORK/sage_astral.json --write-pin -o $WORK/sage_astral $MZML > $RES/sage_astral.log 2>&1
echo "  exit=$?"; wallrss $RES/sage_astral.log
[ -f $WORK/sage_astral/results.sage.pin ] && cp -f $WORK/sage_astral/results.sage.pin $RES/sage_astral.pin
echo "  rows=$(($(wc -l < $RES/sage_astral.pin 2>/dev/null || echo 1)-1))"; perc sage_astral; count sage_astral

echo "=== [5] MSFragger 4.2 DDA+ data_type=3 (mzML) $(date -Is) ==="
# MSFragger runs on the mzML here (the ext/thermo Batmass-IO reader is not installed).
cp -f $MZML $WORK/astral.mzML
/usr/bin/time -v java -Xmx14g -jar $FRAGGER $FRGCFG $WORK/astral.mzML > $RES/fragger_astral.log 2>&1
echo "  exit=$?"; wallrss $RES/fragger_astral.log
[ -f $WORK/astral.pin ] && cp -f $WORK/astral.pin $RES/fragger_astral.pin
echo "  rows=$(($(wc -l < $RES/fragger_astral.pin 2>/dev/null || echo 1)-1))"; perc fragger_astral; count fragger_astral
rm -f $WORK/astral.mzML $WORK/astral.pin $WORK/astral.pepXML $WORK/astral.tsv

echo "=== [6] ProSE (native .raw, OpenMS) $(date -Is) ==="
mkdir -p $RES/prose-out
/usr/bin/time -v docker run --rm -v $AD:/data:ro -v $RES/prose-out:/out $OIMG bash -c '
/opt/OpenMS/bin/ProSE -in /data/'"$(basename $RAW)"' -database /data/'"$(basename $FASTA)"' -out_idxml /out/prose.idXML -Search:decoys -Search:enzyme Trypsin -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 7 -Search:peptide:max_size 40 -Search:peptide:enzyme_specificity full -Search:precursor:mass_tolerance_lower 10 -Search:precursor:mass_tolerance_upper 10 -Search:precursor:mass_tolerance_unit ppm -Search:precursor:min_charge 2 -Search:precursor:max_charge 4 -Search:precursor:isotope_error_min -1 -Search:precursor:isotope_error_max 2 -Search:fragment:mass_tolerance 20 -Search:fragment:mass_tolerance_unit ppm -Search:modifications:fixed "Carbamidomethyl (C)" -Search:modifications:variable "Oxidation (M)" "Acetyl (Protein N-term)" -Search:modifications:variable_max_per_peptide 3 -threads 8 2>&1 | tail -4
echo ---PSMFeatureExtractor---; /opt/OpenMS/bin/PSMFeatureExtractor -in /out/prose.idXML -out /out/prose_feat.idXML 2>&1 | tail -2
echo ---PercolatorAdapter---; /opt/OpenMS/bin/PercolatorAdapter -in /out/prose_feat.idXML -out /out/prose_perc.idXML -seed 42 -post_processing_tdc -percolator_executable /opt/OpenMS/thirdparty/Percolator/percolator 2>&1 | tail -4
echo ---TextExporter---; /opt/OpenMS/bin/TextExporter -in /out/prose_perc.idXML -out /out/prose.tsv 2>&1 | tail -2
echo ---FDRfallback---; /opt/OpenMS/bin/FalseDiscoveryRate -in /out/prose.idXML -out /out/prose_fdr.idXML 2>&1 | tail -2; /opt/OpenMS/bin/TextExporter -in /out/prose_fdr.idXML -out /out/prose_fdr.tsv 2>&1 | tail -1
' > $RES/prose.log 2>&1
echo "  ProSE chain exit=$?"; wallrss $RES/prose.log
echo "  prose.idXML PeptideHits: $(grep -c PeptideHit $RES/prose-out/prose.idXML 2>/dev/null)"
ls -la $RES/prose-out/ 2>&1 | grep -E "tsv|idXML"
echo "################ ASTRAL_5ENGINE_DONE $(date -Is) ################"
