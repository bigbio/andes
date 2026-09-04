#!/bin/bash
# Runs baseline + new JARs on PXD001819 with -outputFormat 3 (pin) to validate
# DirectPinWriter schema changes (Achievement A: lnDeltaSpecEValue + matchedIonRatio
# + enzN/enzC/enzInt/mass + OpenMS header renames).
#
# Output: benchmark/results/PXD001819/{baseline,branch}_output.pin + logs
#
# Usage: bash run_pxd001819_pin.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data/PXD001819"
RESULTS_DIR="$SCRIPT_DIR/results/PXD001819"
BASELINE_JAR="$SCRIPT_DIR/baseline-fork-pre-pr22/MSGFPlus.jar"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"

MZML="$DATA_DIR/UPS1_5000amol_R1.mzML"
FASTA="$DATA_DIR/PXD001819_uniprot_yeast_ups.fasta"
MODS="$DATA_DIR/mods.txt"

THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-4096m}"

# Same SDRF-derived params as the .mzid benchmark + -outputFormat 3 for pin.
SEARCH_ARGS="-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -outputFormat 3 -thread $THREAD_COUNT"

mkdir -p "$RESULTS_DIR"

for f in "$MZML" "$FASTA" "$MODS"; do
    [ -f "$f" ] || { echo "ERROR: $f missing (run download_data.sh first)"; exit 1; }
done

run_one() {
    local JAR="$1"
    local LABEL="$2"
    # MS-GF+ requires -o to end in .mzid; it internally rewrites to .pin when outputFormat=pin.
    local OUTPUT_MZID="$RESULTS_DIR/${LABEL}_output.mzid"
    local OUTPUT="$RESULTS_DIR/${LABEL}_output.pin"
    local LOG="$RESULTS_DIR/${LABEL}_pin.log"

    [ -f "$JAR" ] || { echo "ERROR: $JAR missing"; return 1; }

    echo "=== $LABEL ==="
    # Clear cached suffix array for fair timing.
    rm -f "$DATA_DIR"/*.canno "$DATA_DIR"/*.cnlcp "$DATA_DIR"/*.csarr "$DATA_DIR"/*.cseq 2>/dev/null || true

    local START=$(python3 -c "import time; print(time.time())")

    /usr/bin/time -l java "-Xmx${JVM_MEMORY}" -jar "$JAR" \
        -s "$MZML" -d "$FASTA" -mod "$MODS" -o "$OUTPUT_MZID" \
        $SEARCH_ARGS > "$LOG" 2>&1 || true

    local END=$(python3 -c "import time; print(time.time())")
    local ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")
    local RSS=$(awk '/maximum resident set size/ {print $1; exit}' "$LOG" 2>/dev/null || echo 0)
    local RSS_MB=$(python3 -c "print(f'{int($RSS or 0)/1048576:.0f}')")

    if [ -f "$OUTPUT" ]; then
        local LINES=$(wc -l < "$OUTPUT")
        local TARGETS=$(awk -F'\t' 'NR>1 && $2==1 {c++} END {print c+0}' "$OUTPUT")
        local DECOYS=$(awk -F'\t' 'NR>1 && $2==-1 {c++} END {print c+0}' "$OUTPUT")
        echo "[$LABEL] wall=${ELAPSED}s peak_rss=${RSS_MB}MB rows=$LINES targets=$TARGETS decoys=$DECOYS"
    else
        echo "[$LABEL] FAILED: no .pin produced. See $LOG"
        tail -20 "$LOG" || true
        return 1
    fi
}

run_one "$BASELINE_JAR" "baseline" || true
run_one "$NEW_JAR" "branch" || true

echo ""
echo "Done. Artifacts in $RESULTS_DIR/"
ls -la "$RESULTS_DIR"/*.pin "$RESULTS_DIR"/*_pin.log 2>/dev/null | awk '{print $NF, $5}' || true
