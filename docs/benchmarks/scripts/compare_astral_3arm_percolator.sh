#!/bin/bash
# Run Percolator on Astral 3-arm pin files and print comparison.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

RESULTS_DIR="$SCRIPT_DIR/results/Astral"
PERCOLATOR_DIR="$RESULTS_DIR/percolator_3arm"
mkdir -p "$PERCOLATOR_DIR"

for arm in astralA_baseline astralB_AonlyOff astralC_AplusB; do
    PIN="$RESULTS_DIR/${arm}_output.pin"
    [ -f "$PIN" ] || { echo "skip $arm"; continue; }
    bash "$SCRIPT_DIR/run_percolator_docker.sh" "$PIN" "$PERCOLATOR_DIR" "$arm"
done

echo ""
echo "=== Astral 3-arm comparison @ 1% FDR ==="
printf "%-22s %12s %12s %12s %12s\n" arm targets_1pct targets_5pct decoys_total targets_total
printf "%-22s %12s %12s %12s %12s\n" ---------- ------------ ------------ ------------ ------------
for arm in astralA_baseline astralB_AonlyOff astralC_AplusB; do
    T="$PERCOLATOR_DIR/${arm}.target.psms.txt"
    D="$PERCOLATOR_DIR/${arm}.decoy.psms.txt"
    [ -f "$T" ] || { printf "%-22s %12s\n" "$arm" "FAILED"; continue; }
    TTOT=$(($(wc -l < "$T") - 1))
    T1=$(benchmark_percolator_targets_at_threshold "$T" 0.01)
    T5=$(benchmark_percolator_targets_at_threshold "$T" 0.05)
    DTOT=$(($(wc -l < "$D") - 1))
    printf "%-22s %12d %12d %12d %12d\n" "$arm" "$T1" "$T5" "$DTOT" "$TTOT"
done
