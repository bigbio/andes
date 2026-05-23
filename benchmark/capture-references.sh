#!/usr/bin/env bash
# Capture Java reference .pin outputs for the three sign-off datasets at both
# precursorCal modes. Run via `mvn -Pcapture-references package`. Output lands
# in references/ (gitignored).
#
# This script assumes the bench-machine layout already in use during the
# Phase B / msnet trainer work:
#   /srv/data/msgf-bench/astral-data/    Astral mzML + FASTA + mods
#   /srv/data/msgf-bench/tmt-data/       TMT mzML + FASTA + mods
#   /srv/data/msgf-bench/data/           PXD001819 mzML + FASTA
#
# Locally (macOS/Linux dev), pass DATA_ROOT explicitly to override.
set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/srv/data/msgf-bench}"
OUT_DIR="${OUT_DIR:-references}"
JAR="${JAR:-target/MSGFPlus.jar}"

mkdir -p "$OUT_DIR"

run_one() {
  local label="$1"
  local mzml="$2"
  local fasta="$3"
  local mods="$4"
  local args="$5"
  local cal="$6"
  local out="$OUT_DIR/${label}_cal-${cal}.pin"
  echo "[$label cal=$cal] -> $out"
  java -Xmx8192m -jar "$JAR" \
    -s "$mzml" -d "$fasta" -mod "$mods" -o "$out" $args -precursorCal "$cal"
}

# Astral
run_one astral \
  "$DATA_ROOT/astral-data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML" \
  "$DATA_ROOT/astral-data/ProteoBenchFASTA_MixedSpecies_HYE.fasta" \
  "$DATA_ROOT/astral-data/mods.txt" \
  "-tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread 8" \
  off
run_one astral \
  "$DATA_ROOT/astral-data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML" \
  "$DATA_ROOT/astral-data/ProteoBenchFASTA_MixedSpecies_HYE.fasta" \
  "$DATA_ROOT/astral-data/mods.txt" \
  "-tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread 8" \
  auto

# TMT
run_one tmt \
  "$DATA_ROOT/tmt-data/a05058.mzML" \
  "$DATA_ROOT/tmt-data/PXD007683_UP000005640_UP000002311_reviewed.fasta" \
  "$DATA_ROOT/tmt-data/mods.txt" \
  "-tda 1 -t 20ppm -ti -1,2 -m 1 -inst 1 -e 1 -protocol 4 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread 8" \
  off
run_one tmt \
  "$DATA_ROOT/tmt-data/a05058.mzML" \
  "$DATA_ROOT/tmt-data/PXD007683_UP000005640_UP000002311_reviewed.fasta" \
  "$DATA_ROOT/tmt-data/mods.txt" \
  "-tda 1 -t 20ppm -ti -1,2 -m 1 -inst 1 -e 1 -protocol 4 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread 8" \
  auto

# PXD001819
run_one pxd001819 \
  "$DATA_ROOT/data/UPS1_5000amol_R1.mzML" \
  "$DATA_ROOT/data/PXD001819_uniprot_yeast_ups.fasta" \
  "$DATA_ROOT/mods.txt" \
  "-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread 8" \
  off
run_one pxd001819 \
  "$DATA_ROOT/data/UPS1_5000amol_R1.mzML" \
  "$DATA_ROOT/data/PXD001819_uniprot_yeast_ups.fasta" \
  "$DATA_ROOT/mods.txt" \
  "-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread 8" \
  auto

echo "All reference captures done. Outputs in $OUT_DIR/."
