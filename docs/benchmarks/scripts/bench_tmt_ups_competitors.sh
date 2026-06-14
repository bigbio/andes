#!/bin/bash
# Fresh competitor benchmark — TMT a05058 + UPS1 (mzML). Engines:
# Comet, Java MS-GF+, Sage, ProSE -> .pin -> uniform Percolator(seed42) ->
# PSMs/peptides/proteins@1%. andes measured separately. Waits for the Astral
# competitor job to finish first (clean wall-times, no VM contention).
set -uo pipefail
B=/srv/data/msgf-bench
JAR=$B/MSGFPlus_v20240326.jar
SAGE=$B/engines/sage-v0.14.7-x86_64-unknown-linux-gnu/sage
BP=$B/repo/msgf-rust/docs/benchmarks/scripts/build_pins.py
OIMG=ghcr.io/openms/openms-tools-thirdparty:latest
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
export DOTNET_ROOT=/opt/dotnet8; export PATH=/opt/dotnet8:$PATH

echo "### waiting for Astral competitor job to finish (avoid contention) ###"
for i in $(seq 1 360); do grep -q "ASTRAL_PUB_DONE" $B/bench-pub-astral/run.log 2>/dev/null && break; sleep 30; done
echo "### Astral done (or 3h timeout); starting TMT+UPS $(date -Is) ###"

run_perc_count(){ # $1=RES dir  $2=label
  local RES=$1 lab=$2
  docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y --results-psms /r/$lab.t.psms --decoy-results-psms /r/$lab.d.psms --only-psms=false /r/$lab.pin > $RES/$lab.perc.log 2>&1
  local tp=$RES/$lab.t.psms; [ -s "$tp" ] || { echo "  RESULT $lab NO_PSMS"; return; }
  local q pc rc ps pep pr
  q=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")print i}' "$tp")
  pc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="peptide")print i}' "$tp")
  rc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="proteinIds")print i}' "$tp")
  ps=$(awk -F"\t" -v q="$q" 'NR>1&&$q<=0.01{c++}END{print c+0}' "$tp")
  pep=$(awk -F"\t" -v q="$q" -v p="$pc" 'NR>1&&$q<=0.01{s=$p;gsub(/^[A-Z-]\./,"",s);gsub(/\.[A-Z-]$/,"",s);gsub(/\[[^]]*\]/,"",s);print s}' "$tp"|sort -u|wc -l)
  pr=$(awk -F"\t" -v q="$q" -v r="$rc" 'NR>1&&$q<=0.01{print $r}' "$tp"|tr "\t" "\n"|grep -v "^XXX_\|^rev_\|^DECOY\|^$"|sort -u|wc -l)
  printf "  RESULT %-16s PSMs@1%%=%s peptides@1%%=%s proteins@1%%=%s\n" "$lab" "$ps" "$pep" "$pr"; }
wall(){ grep -E "Elapsed \(wall" $1 2>/dev/null | sed 's/^[[:space:]]*/    /'; }

comet_run(){ # $1=datadir $2=fasta $3=mzml $4=RES $5=extra-sed (TMT mods etc)  $6=tol  $7=binTol $8=binOff $9=theo
  local DD=$1 FA=$2 MZ=$3 RES=$4 EXTRA=$5 TOL=$6 BT=$7 BO=$8 TH=$9
  /usr/bin/time -v docker run --rm -v $DD:/data -v $RES:/out $OIMG bash -c "
cd /out && C=/opt/OpenMS/thirdparty/Comet/comet.exe; \$C -p >/dev/null 2>&1; P=comet.params.new
sed -i -E \"s#^database_name = .*#database_name = /data/$FA#; s#^decoy_search = .*#decoy_search = 1#; s#^num_threads = .*#num_threads = 8#; s#^peptide_mass_tolerance_upper = .*#peptide_mass_tolerance_upper = $TOL#; s#^peptide_mass_tolerance_lower = .*#peptide_mass_tolerance_lower = -$TOL#; s#^peptide_mass_units = .*#peptide_mass_units = 2#; s#^isotope_error = .*#isotope_error = 4#; s#^search_enzyme_number = .*#search_enzyme_number = 1#; s#^allowed_missed_cleavage = .*#allowed_missed_cleavage = 2#; s#^num_enzyme_termini = .*#num_enzyme_termini = 2#; s#^fragment_bin_tol = .*#fragment_bin_tol = $BT#; s#^fragment_bin_offset = .*#fragment_bin_offset = $BO#; s#^theoretical_fragment_ions = .*#theoretical_fragment_ions = $TH#; s#^output_percolatorfile = .*#output_percolatorfile = 1#; s#^output_txtfile = .*#output_txtfile = 0#; s#^output_pepxmlfile = .*#output_pepxmlfile = 0#; s#^add_C_cysteine = .*#add_C_cysteine = 57.021464#; s#^precursor_charge = .*#precursor_charge = 2 4#; s#^variable_mod01 = .*#variable_mod01 = 15.994915 M 0 3 -1 0 0 0.0#; $EXTRA\" \$P
grep -q '^peptide_length_range' \$P && sed -i -E 's#^peptide_length_range = .*#peptide_length_range = 6 40#' \$P
cp /data/$MZ /out/comet_in.mzML; \$C -P\$P /out/comet_in.mzML 2>&1 | tail -4; rm -f /out/comet_in.mzML
" > $RES/comet.log 2>&1; echo "  comet exit=$?"; wall $RES/comet.log
  [ -f $RES/comet_in.pin ] && mv -f $RES/comet_in.pin $RES/comet.pin
  echo "  rows=$(($(wc -l < $RES/comet.pin 2>/dev/null||echo 1)-1))"; run_perc_count $RES comet; rm -f $RES/comet.pin; }

############################## TMT a05058 ##############################
TD=$B/tmt-data; RES=$B/bench-pub-tmt; mkdir -p $RES; W=/tmp/pub-tmt; mkdir -p $W
MZ=a05058.mzML; FA=PXD007683_UP000005640_UP000002311_reviewed.fasta; MODS=$TD/mods-numeric.txt
echo "################ TMT a05058 PUBLIC $(date -Is) ################"
echo "=== [Comet] ==="
comet_run $TD $FA $MZ $RES 's#^add_K_lysine = .*#add_K_lysine = 229.162932#; s#^add_Nterm_peptide = .*#add_Nterm_peptide = 229.162932#;' 20.0 1.0005 0.4 1
echo "=== [Java MS-GF+] (-inst 0 -m 1 -protocol 4) ==="
/usr/bin/time -v java -Xmx14g -jar $JAR -s $TD/$MZ -d $TD/$FA -mod $MODS -o $W/java.mzid -tda 1 -t 20ppm -ti -1,2 -m 1 -inst 0 -e 1 -protocol 4 -ntt 2 -minLength 6 -maxLength 40 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -thread 8 > $RES/java.log 2>&1; echo "  java exit=$?"; wall $RES/java.log
if docker run --rm --platform linux/amd64 -v "$W":/w -v "$RES":/r $PIMG msgf2pin -P XXX_ -e trypsin -o /r/java.pin /w/java.mzid > $RES/java_msgf2pin.log 2>&1 && [ -s $RES/java.pin ]; then echo "  msgf2pin OK"; else echo "  msgf2pin failed -> MzIDToTsv+build_pins"; java -Xmx14g -cp $JAR edu.ucsd.msjava.ui.MzIDToTsv -i $W/java.mzid -o $W/java.tsv -showDecoy 1 > $RES/java_mzidtotsv.log 2>&1; python3 $BP java $W/java.tsv $RES/java.pin > $RES/java_buildpin.log 2>&1; fi
echo "  rows=$(($(wc -l < $RES/java.pin 2>/dev/null||echo 1)-1))"; run_perc_count $RES java; rm -f $RES/java.pin $W/java.mzid $W/java.tsv
echo "=== [Sage] ==="
sed "s#SAGE_OUT_PLACEHOLDER#$W/sage#g" $TD/sage-a05058.json > $W/sage.json
/usr/bin/time -v $SAGE $W/sage.json --write-pin -o $W/sage $TD/$MZ > $RES/sage.log 2>&1; echo "  sage exit=$?"; wall $RES/sage.log
[ -f $W/sage/results.sage.pin ] && cp -f $W/sage/results.sage.pin $RES/sage.pin
echo "  rows=$(($(wc -l < $RES/sage.pin 2>/dev/null||echo 1)-1))"; run_perc_count $RES sage; rm -f $RES/sage.pin
echo "=== [ProSE] (0.4 Da, TMT mods) ==="
mkdir -p $RES/prose-out
/usr/bin/time -v docker run --rm -v $TD:/data:ro -v $RES/prose-out:/out $OIMG bash -c '
/opt/OpenMS/bin/ProSE -in /data/a05058.mzML -database /data/PXD007683_UP000005640_UP000002311_reviewed.fasta -out_idxml /out/prose.idXML -Search:decoys -Search:enzyme Trypsin -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 6 -Search:peptide:max_size 40 -Search:peptide:enzyme_specificity full -Search:precursor:mass_tolerance_lower 20 -Search:precursor:mass_tolerance_upper 20 -Search:precursor:mass_tolerance_unit ppm -Search:precursor:min_charge 2 -Search:precursor:max_charge 4 -Search:precursor:isotope_error_min -1 -Search:precursor:isotope_error_max 2 -Search:fragment:mass_tolerance 0.4 -Search:fragment:mass_tolerance_unit Da -Search:modifications:fixed "Carbamidomethyl (C)" "TMT6plex (K)" "TMT6plex (N-term)" -Search:modifications:variable "Oxidation (M)" -Search:modifications:variable_max_per_peptide 3 -threads 8 2>&1 | tail -3
' > $RES/prose.log 2>&1; echo "  prose exit=$?"; wall $RES/prose.log
python3 $BP prose $RES/prose-out/prose.idXML $RES/prose.pin > $RES/prose_build.log 2>&1
echo "  prose rows=$(($(wc -l < $RES/prose.pin 2>/dev/null||echo 1)-1))"; run_perc_count $RES prose; rm -f $RES/prose.pin
echo "################ TMT_PUB_DONE $(date -Is) ################"; df -h /srv/data|tail -1

############################## UPS1 PXD001819 ##############################
DD=$B/data; RES=$B/bench-pub-ups1; mkdir -p $RES; W=/tmp/pub-ups1; mkdir -p $W
MZ=UPS1_5000amol_R1.mzML; FA=PXD001819_uniprot_yeast_ups.fasta; MODS=$DD/mods-ups1-clean.txt
echo "################ UPS1 PXD001819 PUBLIC $(date -Is) ################"
echo "=== [Comet] (20ppm to match other engines) ==="
comet_run $DD $FA $MZ $RES 's#^variable_mod02 = .*#variable_mod02 = 42.010565 n 0 1 0 0 0 0.0#;' 20.0 1.0005 0.4 1
echo "=== [Java MS-GF+] (-inst 0 -m 1 -protocol 0) ==="
/usr/bin/time -v java -Xmx14g -jar $JAR -s $DD/$MZ -d $DD/$FA -mod $MODS -o $W/java.mzid -tda 1 -t 20ppm -ti -1,2 -m 1 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -thread 8 > $RES/java.log 2>&1; echo "  java exit=$?"; wall $RES/java.log
if docker run --rm --platform linux/amd64 -v "$W":/w -v "$RES":/r $PIMG msgf2pin -P XXX_ -e trypsin -o /r/java.pin /w/java.mzid > $RES/java_msgf2pin.log 2>&1 && [ -s $RES/java.pin ]; then echo "  msgf2pin OK"; else echo "  msgf2pin failed -> MzIDToTsv+build_pins"; java -Xmx14g -cp $JAR edu.ucsd.msjava.ui.MzIDToTsv -i $W/java.mzid -o $W/java.tsv -showDecoy 1 > $RES/java_mzidtotsv.log 2>&1; python3 $BP java $W/java.tsv $RES/java.pin > $RES/java_buildpin.log 2>&1; fi
echo "  rows=$(($(wc -l < $RES/java.pin 2>/dev/null||echo 1)-1))"; run_perc_count $RES java; rm -f $RES/java.pin $W/java.mzid $W/java.tsv
echo "=== [Sage] ==="
sed "s#SAGE_OUT_PLACEHOLDER#$W/sage#g" $DD/sage-ups1.json > $W/sage.json
/usr/bin/time -v $SAGE $W/sage.json --write-pin -o $W/sage $DD/$MZ > $RES/sage.log 2>&1; echo "  sage exit=$?"; wall $RES/sage.log
[ -f $W/sage/results.sage.pin ] && cp -f $W/sage/results.sage.pin $RES/sage.pin
echo "  rows=$(($(wc -l < $RES/sage.pin 2>/dev/null||echo 1)-1))"; run_perc_count $RES sage; rm -f $RES/sage.pin
echo "=== [ProSE] (0.4 Da) ==="
mkdir -p $RES/prose-out
/usr/bin/time -v docker run --rm -v $DD:/data:ro -v $RES/prose-out:/out $OIMG bash -c '
/opt/OpenMS/bin/ProSE -in /data/UPS1_5000amol_R1.mzML -database /data/PXD001819_uniprot_yeast_ups.fasta -out_idxml /out/prose.idXML -Search:decoys -Search:enzyme Trypsin -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 6 -Search:peptide:max_size 40 -Search:peptide:enzyme_specificity full -Search:precursor:mass_tolerance_lower 20 -Search:precursor:mass_tolerance_upper 20 -Search:precursor:mass_tolerance_unit ppm -Search:precursor:min_charge 2 -Search:precursor:max_charge 4 -Search:precursor:isotope_error_min -1 -Search:precursor:isotope_error_max 2 -Search:fragment:mass_tolerance 0.4 -Search:fragment:mass_tolerance_unit Da -Search:modifications:fixed "Carbamidomethyl (C)" -Search:modifications:variable "Oxidation (M)" -Search:modifications:variable_max_per_peptide 3 -threads 8 2>&1 | tail -3
' > $RES/prose.log 2>&1; echo "  prose exit=$?"; wall $RES/prose.log
python3 $BP prose $RES/prose-out/prose.idXML $RES/prose.pin > $RES/prose_build.log 2>&1
echo "  prose rows=$(($(wc -l < $RES/prose.pin 2>/dev/null||echo 1)-1))"; run_perc_count $RES prose; rm -f $RES/prose.pin
echo "################ UPS1_PUB_DONE $(date -Is) ################"; df -h /srv/data|tail -1
echo "################ ALL_PUB_COMPETITORS_DONE $(date -Is) ################"