#!/bin/bash
# Baseline-only PXD001819 pin run (rebuild baseline against origin/dev 1d481aa).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data/PXD001819"
RESULTS_DIR="$SCRIPT_DIR/results/PXD001819"
BASELINE_JAR="$SCRIPT_DIR/baseline-fork-pre-pr22/MSGFPlus.jar"
MZML="$DATA_DIR/UPS1_5000amol_R1.mzML"
FASTA="$DATA_DIR/PXD001819_uniprot_yeast_ups.fasta"
MODS="$DATA_DIR/mods.txt"
THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-4096m}"
SEARCH_ARGS="-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -outputFormat 3 -thread $THREAD_COUNT"
mkdir -p "$RESULTS_DIR"
OUTPUT_MZID="$RESULTS_DIR/baseline_output.mzid"
OUTPUT="$RESULTS_DIR/baseline_output.pin"
LOG="$RESULTS_DIR/baseline_pin.log"
rm -f "$DATA_DIR"/*.canno "$DATA_DIR"/*.cnlcp "$DATA_DIR"/*.csarr "$DATA_DIR"/*.cseq 2>/dev/null || true
START=$(python3 -c "import time; print(time.time())")
/usr/bin/time -l java "-Xmx${JVM_MEMORY}" -jar "$BASELINE_JAR" \
    -s "$MZML" -d "$FASTA" -mod "$MODS" -o "$OUTPUT_MZID" \
    $SEARCH_ARGS > "$LOG" 2>&1 || true
END=$(python3 -c "import time; print(time.time())")
ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")
RSS=$(awk '/maximum resident set size/ {print $1; exit}' "$LOG" 2>/dev/null || echo 0)
RSS_MB=$(python3 -c "print(f'{int($RSS or 0)/1048576:.0f}')")
if [ -f "$OUTPUT" ]; then
    TARGETS=$(awk -F'\t' 'NR>1 && $2==1 {c++} END {print c+0}' "$OUTPUT")
    DECOYS=$(awk -F'\t' 'NR>1 && $2==-1 {c++} END {print c+0}' "$OUTPUT")
    echo "[baseline] wall=${ELAPSED}s peak_rss=${RSS_MB}MB targets=$TARGETS decoys=$DECOYS"
else
    echo "[baseline] FAILED"; tail -15 "$LOG"
fi
