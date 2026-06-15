#!/bin/bash
# Fresh competitor benchmark — Astral (mzML, all engines uniform). Engines: andes + Java MS-GF+, Sage, Comet, ProSE.
# Engines: Comet, Sage, ProSE, Java MS-GF+ (slowest last). Each -> .pin -> uniform
# Percolator(seed42) -> PSMs/peptides/proteins@1%. andes already measured separately.
set -uo pipefail
AD=/srv/data/msgf-bench/astral-data
MZML=$AD/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
FASTA=$AD/ProteoBenchFASTA_MixedSpecies_HYE.fasta
MODS=/srv/data/msgf-bench/astral_mods_rust.txt
JAR=/srv/data/msgf-bench/MSGFPlus_v20240326.jar
SAGE=/srv/data/msgf-bench/engines/sage-v0.14.7-x86_64-unknown-linux-gnu/sage
SAGECFG=/srv/data/msgf-bench/repo/msgf-rust/docs/benchmarks/configs/sage-astral.json
BP=/srv/data/msgf-bench/repo/build_pins.py
RES=/srv/data/msgf-bench/bench-pub-astral; mkdir -p $RES
WORK=/tmp/pub-astral; mkdir -p $WORK
OIMG=ghcr.io/openms/openms-tools-thirdparty:latest
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
export DOTNET_ROOT=/opt/dotnet8; export PATH=/opt/dotnet8:$PATH
perc(){ docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y --results-psms /r/$1.t.psms --decoy-results-psms /r/$1.d.psms --only-psms=false /r/$1.pin > $RES/$1.perc.log 2>&1; }
count(){ tp=$RES/$1.t.psms; [ -s "$tp" ] || { echo "  RESULT $1 NO_PSMS"; return; }
  q=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")print i}' "$tp")
  pc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="peptide")print i}' "$tp")
  rc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="proteinIds")print i}' "$tp")
  ps=$(awk -F"\t" -v q="$q" 'NR>1&&$q<=0.01{c++}END{print c+0}' "$tp")
  pep=$(awk -F"\t" -v q="$q" -v p="$pc" 'NR>1&&$q<=0.01{s=$p;gsub(/^[A-Z-]\./,"",s);gsub(/\.[A-Z-]$/,"",s);gsub(/\[[^]]*\]/,"",s);print s}' "$tp"|sort -u|wc -l)
  pr=$(awk -F"\t" -v q="$q" -v r="$rc" 'NR>1&&$q<=0.01{print $r}' "$tp"|tr "\t" "\n"|grep -v "^XXX_\|^rev_\|^DECOY\|^$"|sort -u|wc -l)
  printf "  RESULT %-16s PSMs@1%%=%s peptides@1%%=%s proteins@1%%=%s\n" "$1" "$ps" "$pep" "$pr"; }
wall(){ grep -E "Elapsed \(wall" $1 2>/dev/null | sed 's/^[[:space:]]*/    /'; }

echo "################ ASTRAL PUBLIC BENCHMARK $(date -Is) ################"; df -h /srv/data|tail -1

echo "=== [Comet 2025.01] $(date -Is) ==="
/usr/bin/time -v docker run --rm -v $AD:/data -v $RES:/out $OIMG bash -c '
cd /out && C=/opt/OpenMS/thirdparty/Comet/comet.exe; $C -p >/dev/null 2>&1; P=comet.params.new
sed -i -E "s#^database_name = .*#database_name = /data/'"$(basename $FASTA)"'#; s#^decoy_search = .*#decoy_search = 1#; s#^num_threads = .*#num_threads = 8#; s#^peptide_mass_tolerance_upper = .*#peptide_mass_tolerance_upper = 10.0#; s#^peptide_mass_tolerance_lower = .*#peptide_mass_tolerance_lower = -10.0#; s#^peptide_mass_units = .*#peptide_mass_units = 2#; s#^isotope_error = .*#isotope_error = 4#; s#^search_enzyme_number = .*#search_enzyme_number = 1#; s#^allowed_missed_cleavage = .*#allowed_missed_cleavage = 2#; s#^num_enzyme_termini = .*#num_enzyme_termini = 2#; s#^fragment_bin_tol = .*#fragment_bin_tol = 0.02#; s#^fragment_bin_offset = .*#fragment_bin_offset = 0.0#; s#^theoretical_fragment_ions = .*#theoretical_fragment_ions = 0#; s#^output_percolatorfile = .*#output_percolatorfile = 1#; s#^output_txtfile = .*#output_txtfile = 0#; s#^output_pepxmlfile = .*#output_pepxmlfile = 0#; s#^add_C_cysteine = .*#add_C_cysteine = 57.021464#; s#^precursor_charge = .*#precursor_charge = 2 4#; s#^variable_mod01 = .*#variable_mod01 = 15.994915 M 0 3 -1 0 0 0.0#; s#^variable_mod02 = .*#variable_mod02 = 42.010565 n 0 1 0 0 0 0.0#;" $P
grep -q "^peptide_length_range" $P && sed -i -E "s#^peptide_length_range = .*#peptide_length_range = 7 40#" $P
cp /data/'"$(basename $MZML)"' /out/comet_in.mzML; $C -P$P /out/comet_in.mzML 2>&1 | tail -4; rm -f /out/comet_in.mzML
' > $RES/comet.log 2>&1; echo "  comet exit=$?"; wall $RES/comet.log
[ -f $RES/comet_in.pin ] && mv -f $RES/comet_in.pin $RES/comet.pin
echo "  rows=$(($(wc -l < $RES/comet.pin 2>/dev/null||echo 1)-1))"; perc comet; count comet; rm -f $RES/comet.pin

echo "=== [Sage 0.14.7] $(date -Is) ==="
sed "s#\"fasta\":.*#\"fasta\": \"$FASTA\",#; s#\"output_directory\":.*#\"output_directory\": \"$WORK/sage\"#" $SAGECFG > $WORK/sage.json
/usr/bin/time -v $SAGE $WORK/sage.json --write-pin -o $WORK/sage $MZML > $RES/sage.log 2>&1; echo "  sage exit=$?"; wall $RES/sage.log
[ -f $WORK/sage/results.sage.pin ] && cp -f $WORK/sage/results.sage.pin $RES/sage.pin
echo "  rows=$(($(wc -l < $RES/sage.pin 2>/dev/null||echo 1)-1))"; perc sage; count sage; rm -f $RES/sage.pin

echo "=== [ProSE / OpenMS] (mzML; native .raw not on disk) $(date -Is) ==="
mkdir -p $RES/prose-out
/usr/bin/time -v docker run --rm -v $AD:/data:ro -v $RES/prose-out:/out $OIMG bash -c '
/opt/OpenMS/bin/ProSE -in /data/'"$(basename $MZML)"' -database /data/'"$(basename $FASTA)"' -out_idxml /out/prose.idXML -Search:decoys -Search:enzyme Trypsin -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 7 -Search:peptide:max_size 40 -Search:peptide:enzyme_specificity full -Search:precursor:mass_tolerance_lower 10 -Search:precursor:mass_tolerance_upper 10 -Search:precursor:mass_tolerance_unit ppm -Search:precursor:min_charge 2 -Search:precursor:max_charge 4 -Search:precursor:isotope_error_min -1 -Search:precursor:isotope_error_max 2 -Search:fragment:mass_tolerance 20 -Search:fragment:mass_tolerance_unit ppm -Search:modifications:fixed "Carbamidomethyl (C)" -Search:modifications:variable "Oxidation (M)" "Acetyl (Protein N-term)" -Search:modifications:variable_max_per_peptide 3 -threads 8 2>&1 | tail -3
' > $RES/prose.log 2>&1; echo "  ProSE exit=$?"; wall $RES/prose.log
python3 $BP prose $RES/prose-out/prose.idXML $RES/prose.pin > $RES/prose_build.log 2>&1
echo "  prose pin rows=$(($(wc -l < $RES/prose.pin 2>/dev/null||echo 1)-1))"; perc prose; count prose; rm -f $RES/prose.pin

echo "=== [Java MS-GF+ v20240326] (slow ~2h on Astral) $(date -Is) ==="
/usr/bin/time -v java -Xmx14g -jar $JAR -s $MZML -d $FASTA -mod $MODS -o $WORK/java.mzid -tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 7 -maxLength 40 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -thread 8 > $RES/java.log 2>&1; echo "  java exit=$?"; wall $RES/java.log
docker run --rm --platform linux/amd64 -v "$WORK":/w -v "$RES":/r $PIMG msgf2pin -P XXX_ -e trypsin -o /r/java.pin /w/java.mzid > $RES/java_msgf2pin.log 2>&1
echo "  msgf2pin exit=$? rows=$(($(wc -l < $RES/java.pin 2>/dev/null||echo 1)-1))"; perc java; count java; rm -f $RES/java.pin $WORK/java.mzid

echo "################ ASTRAL_PUB_DONE $(date -Is) ################"; df -h /srv/data|tail -1