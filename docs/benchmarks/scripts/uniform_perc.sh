#!/bin/bash
set -uo pipefail
RES=/srv/data/msgf-bench/repo/bench-astral
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
# build_pins.py lives next to this script, not at the repo root.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
echo "################ UNIFORM PERCOLATOR (Astral, all 7) $(date -Is) ################"
echo "=== building Java PIN from MzIDToTsv output ==="
python3 "$SCRIPT_DIR/build_pins.py" java "$RES/java_astral.tsv" "$RES/java_astral.pin"
echo "=== building ProSE PIN from idXML annotations ==="
python3 "$SCRIPT_DIR/build_pins.py" prose "$RES/prose-out/prose.idXML" "$RES/prose_astral.pin"

perc(){ docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y --results-psms /r/$1.t.psms --decoy-results-psms /r/$1.d.psms --only-psms=false /r/$1.pin > $RES/$1.perc.log 2>&1; echo "  perc $1 exit=$?"; }
count(){ tp=$RES/$1.t.psms
  q=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")print i}' "$tp")
  rc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="proteinIds")print i}' "$tp")
  pc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="peptide")print i}' "$tp")
  ps=$(awk -F"\t" -v q="$q" 'NR>1&&$q<=0.01{c++}END{print c+0}' "$tp")
  pep=$(awk -F"\t" -v q="$q" -v p="$pc" 'NR>1&&$q<=0.01{s=$p;gsub(/^[A-Z-]\./,"",s);gsub(/\.[A-Z-]$/,"",s);gsub(/\[[^]]*\]/,"",s);gsub(/[^A-Z]/,"",s);print s}' "$tp"|sort -u|wc -l)
  pr=$(awk -F"\t" -v q="$q" -v r="$rc" 'NR>1&&$q<=0.01{print $r}' "$tp"|tr "\t" "\n"|grep -vE "^XXX_|^rev_|^DECOY_|^$|^unknown$"|sort -u|wc -l)
  printf "  %-16s PSMs@1%%=%-7s peptides@1%%=%-7s proteins@1%%=%s\n" "$1" "$ps" "$pep" "$pr"; }

echo "=== percolating Java + ProSE (new PINs) ==="
perc java_astral; perc prose_astral
echo "================ UNIFORM PERCOLATOR RESULTS (all via Percolator 3.7.1, seed 42) ================"
for e in simas_top1 simas_chim fragger_astral sage_astral comet_astral java_astral prose_astral; do count $e; done
echo "################ UNIFORM_PERC_DONE $(date -Is) ################"
