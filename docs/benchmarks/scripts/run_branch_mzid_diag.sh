#!/bin/bash
# Diagnostic: run the branch (52cab74) JAR with -outputFormat 0 (mzid) on PXD001819
# to isolate whether the 6x slowdown is in the pin writer or the scoring path.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data/PXD001819"
RESULTS_DIR="$SCRIPT_DIR/results/PXD001819"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"
MZML="$DATA_DIR/UPS1_5000amol_R1.mzML"
FASTA="$DATA_DIR/PXD001819_uniprot_yeast_ups.fasta"
MODS="$DATA_DIR/mods.txt"
THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-4096m}"
SEARCH_ARGS="-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -outputFormat 0 -thread $THREAD_COUNT"
mkdir -p "$RESULTS_DIR"
OUTPUT="$RESULTS_DIR/branch_mzid_diag_output.mzid"
LOG="$RESULTS_DIR/branch_mzid_diag.log"
rm -f "$DATA_DIR"/*.canno "$DATA_DIR"/*.cnlcp "$DATA_DIR"/*.csarr "$DATA_DIR"/*.cseq 2>/dev/null || true
START=$(python3 -c "import time; print(time.time())")
/usr/bin/time -l java "-Xmx${JVM_MEMORY}" -jar "$NEW_JAR" \
    -s "$MZML" -d "$FASTA" -mod "$MODS" -o "$OUTPUT" \
    $SEARCH_ARGS > "$LOG" 2>&1 || true
END=$(python3 -c "import time; print(time.time())")
ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")
RSS=$(awk '/maximum resident set size/ {print $1; exit}' "$LOG" 2>/dev/null || echo 0)
RSS_MB=$(python3 -c "print(f'{int($RSS or 0)/1048576:.0f}')")
echo "[branch-mzid] wall=${ELAPSED}s peak_rss=${RSS_MB}MB"
[ -f "$OUTPUT" ] && ls -la "$OUTPUT"
