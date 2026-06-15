#!/bin/bash
# Astral entrapment-FDP experiment — definitively test whether andes's Percolator
# q-values on Astral are honest. Runs andes top-1 + chimeric (params identical to
# the benchmark) against ASTRAL_entrapment.fasta (1:1 real:ENT_), percolates, and
# computes true FDP = 2*ENT/total at a q-vs-FDP calibration curve.
# WAITS for the full competitor benchmark to finish first (no CPU contention ->
# protects the published wall-times).
set -uo pipefail
B=/srv/data/msgf-bench
A=$B/repo/msgf-rust/target/release/andes
SMODEL=$B/repo/msgf-rust/resources/ionstat/models.parquet
MODS=$B/astral_mods_rust.txt
MZ=$B/astral-data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
ENTDB=$B/astral-data/ASTRAL_entrapment.fasta
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
RES=$B/astral-entrapment; mkdir -p $RES
export DOTNET_ROOT=/opt/dotnet8; export PATH=/opt/dotnet8:$PATH

echo "### waiting for competitor benchmark (ALL_PUB_COMPETITORS_DONE) before running (clean wall-times) ###"
for i in $(seq 1 480); do
  grep -q "ALL_PUB_COMPETITORS_DONE" $B/bench-pub-tmt-ups.run.log 2>/dev/null && break
  sleep 30
done
echo "############ ASTRAL ENTRAPMENT-FDP $(date -Is) ############"; df -h /srv/data|tail -1
echo "DB: $(grep -c '^>' $ENTDB) proteins (1:1 real:ENT_); andes adds XXX_ decoys"; echo "binary: $A"; $A --version 2>/dev/null | head -1

run_one(){  # $1=label  $2=extra-flags
  local lab=$1 extra=$2
  echo "=== andes $lab (entrapment DB) $(date -Is) ==="
  /usr/bin/time -v $A --spectrum $MZ --database $ENTDB --mods $MODS --model-store $SMODEL \
     --fragmentation auto --precursor-tol-ppm 10 --isotope-error-min -1 --isotope-error-max 2 \
     --enzyme-specificity fully --max-missed-cleavages 2 --min-length 7 --max-length 40 \
     --charge-min 2 --charge-max 4 --top-n 1 --min-peaks 10 --threads 8 $extra \
     --output-pin $RES/$lab.pin > $RES/$lab.log 2>&1
  echo "  andes $lab exit=$? rows=$(($(wc -l < $RES/$lab.pin 2>/dev/null||echo 1)-1))"; grep -i "Param model" $RES/$lab.log | head -1
  docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y \
    --results-psms /r/$lab.t.psms --decoy-results-psms /r/$lab.d.psms --only-psms=false /r/$lab.pin > $RES/$lab.perc.log 2>&1
  echo "  percolator $lab exit=$?"
  rm -f $RES/$lab.pin   # free disk immediately (entrapment PINs are large)
  df -h /srv/data|tail -1
}
run_one top1 ""
run_one chim "--chimeric"

echo "==================== TRUE ENTRAPMENT-FDP (q-vs-FDP calibration) ===================="
python3 $B/entrap_fdp.py
echo "############ ASTRAL_ENTRAPMENT_DONE $(date -Is) ############"; df -h /srv/data|tail -1
