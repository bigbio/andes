#!/bin/bash
# After run_pxd001819_3arm.sh produces the 3 pin files, this script runs
# Percolator on each and prints a comparison matrix.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

RESULTS_DIR="$SCRIPT_DIR/results/PXD001819"
PERCOLATOR_DIR="$RESULTS_DIR/percolator_3arm"
mkdir -p "$PERCOLATOR_DIR"

for arm in armA_baseline armB_AonlyOff armC_AplusB; do
    PIN="$RESULTS_DIR/${arm}_output.pin"
    [ -f "$PIN" ] || { echo "skip $arm (no pin)"; continue; }
    bash "$SCRIPT_DIR/run_percolator_docker.sh" "$PIN" "$PERCOLATOR_DIR" "$arm"
done

echo ""
echo "=== 3-arm comparison @ 1% FDR ==="
printf "%-18s %12s %12s %12s %12s\n" arm targets_1pct targets_5pct decoys_total targets_total
printf "%-18s %12s %12s %12s %12s\n" ------------ ------------ ------------ ------------ ------------
for arm in armA_baseline armB_AonlyOff armC_AplusB; do
    T="$PERCOLATOR_DIR/${arm}.target.psms.txt"
    D="$PERCOLATOR_DIR/${arm}.decoy.psms.txt"
    [ -f "$T" ] || { printf "%-18s %12s\n" "$arm" "FAILED"; continue; }
    TTOT=$(($(wc -l < "$T") - 1))
    T1=$(benchmark_percolator_targets_at_threshold "$T" 0.01)
    T5=$(benchmark_percolator_targets_at_threshold "$T" 0.05)
    DTOT=$(($(wc -l < "$D") - 1))
    printf "%-18s %12d %12d %12d %12d\n" "$arm" "$T1" "$T5" "$DTOT" "$TTOT"
done
