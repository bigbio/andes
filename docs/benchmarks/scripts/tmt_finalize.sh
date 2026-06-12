#!/bin/bash
set -uo pipefail
TD=$BENCH/tmt-data
RES=$BENCH/repo/bench-tmt
JAR=$BENCH/MSGFPlus_v20240326.jar
OIMG=ghcr.io/openms/openms-tools-thirdparty:latest
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
# build_pins.py lives next to this script, not at the repo root.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd $BENCH/repo
echo "################ TMT FINALIZE (uniform Percolator) $(date -Is) ################"

echo "=== ProSE re-run at 0.1 Da (its max; high-res engine on low-res data) ==="
mkdir -p $RES/prose-out
docker run --rm -v $TD:/data:ro -v $RES/prose-out:/out $OIMG bash -c '
/opt/OpenMS/bin/ProSE -in /data/a05058.mzML -database /data/PXD007683_UP000005640_UP000002311_reviewed.fasta -out_idxml /out/prose.idXML -Search:decoys -Search:enzyme Trypsin -Search:peptide:missed_cleavages 2 -Search:peptide:min_size 6 -Search:peptide:max_size 40 -Search:peptide:enzyme_specificity full -Search:precursor:mass_tolerance_lower 20 -Search:precursor:mass_tolerance_upper 20 -Search:precursor:mass_tolerance_unit ppm -Search:precursor:min_charge 2 -Search:precursor:max_charge 4 -Search:precursor:isotope_error_min -1 -Search:precursor:isotope_error_max 2 -Search:fragment:mass_tolerance 0.1 -Search:fragment:mass_tolerance_unit Da -Search:modifications:fixed "Carbamidomethyl (C)" "TMT6plex (K)" "TMT6plex (N-term)" -Search:modifications:variable "Oxidation (M)" -Search:modifications:variable_max_per_peptide 3 -threads 8 2>&1 | tail -3
' > $RES/prose_rerun.log 2>&1
echo "  ProSE idXML PeptideHits: $(grep -c PeptideHit $RES/prose-out/prose.idXML 2>/dev/null || echo FAILED)"

echo "=== Java MzIDToTsv ==="
java -Xmx8g -cp $JAR edu.ucsd.msjava.ui.MzIDToTsv -i /tmp/tmt-bench/java_tmt.mzid -o $RES/java_tmt.tsv -showQValue 1 -showDecoy 1 -unroll 0 > $RES/java_mzidtotsv.log 2>&1
echo "  java tsv rows=$(($(wc -l < $RES/java_tmt.tsv 2>/dev/null || echo 1)-1))"

echo "=== build PINs ==="
python3 "$SCRIPT_DIR/build_pins.py" java $RES/java_tmt.tsv $RES/java_tmt.pin
python3 "$SCRIPT_DIR/build_pins.py" prose $RES/prose-out/prose.idXML $RES/prose_tmt.pin

perc(){ docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y --results-psms /r/$1.t.psms --decoy-results-psms /r/$1.d.psms --only-psms=false /r/$1.pin > $RES/$1.perc.log 2>&1; echo "  perc $1 exit=$?"; }
count(){ tp=$RES/$1.t.psms; [ -f "$tp" ] || { printf "  %-16s (no percolator output)\n" "$1"; return; }
  q=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")print i}' "$tp")
  rc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="proteinIds")print i}' "$tp")
  pc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="peptide")print i}' "$tp")
  ps=$(awk -F"\t" -v q="$q" 'NR>1&&$q<=0.01{c++}END{print c+0}' "$tp")
  pep=$(awk -F"\t" -v q="$q" -v p="$pc" 'NR>1&&$q<=0.01{s=$p;gsub(/^[A-Z-]\./,"",s);gsub(/\.[A-Z-]$/,"",s);gsub(/\[[^]]*\]/,"",s);gsub(/[^A-Z]/,"",s);print s}' "$tp"|sort -u|wc -l)
  pr=$(awk -F"\t" -v q="$q" -v r="$rc" 'NR>1&&$q<=0.01{for(i=r;i<=NF;i++)print $i}' "$tp"|sed -E "s/\(pre=.*//"|grep -vE "^XXX_|^rev_|^DECOY_|^$|^unknown$"|sort -u|wc -l)
  printf "  %-16s PSMs@1%%=%-7s peptides@1%%=%-7s proteins@1%%=%s\n" "$1" "$ps" "$pep" "$pr"; }

perc java_tmt; perc prose_tmt
echo "================ TMT UNIFORM PERCOLATOR RESULTS (Percolator 3.7.1, seed 42) ================"
for e in simas_top1 fragger_tmt sage_tmt comet_tmt java_tmt prose_tmt; do count $e; done
echo "################ TMT_FINALIZE_DONE $(date -Is) ################"
