#!/bin/bash
# 3-engine benchmark on PXD007683 TMT (Lumos, TMT6plex):
#   arm A — MS-GF+ baseline JAR (origin/dev @ 1d481aa, pre PR #22)
#   arm B — MS-GF+ current dev JAR (origin/dev @ 78af285, post PR #22)
#   arm C — Sage 0.14.7 via biocontainers docker
# All three emit .pin; Percolator 3.7.1 rescores each with the same seed.
# This dataset is the closest match to the Sage paper's 500x wall-time claim
# (PXD016766 was TMT Lumos; this is a sibling TMT run).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data/PXD007683"
RESULTS_DIR="$SCRIPT_DIR/results/PXD007683/3engine"
BASELINE_JAR="$SCRIPT_DIR/baseline-fork-pre-pr22/MSGFPlus.jar"
CURRENT_JAR="$SCRIPT_DIR/current-dev/MSGFPlus.jar"
SAGE_BIN="$SCRIPT_DIR/sage-bin/sage"

MZML="$DATA_DIR/a05058.mzML"
FASTA="$DATA_DIR/PXD007683_UP000005640_UP000002311_reviewed.fasta"
REVCAT_FASTA="$DATA_DIR/PXD007683_UP000005640_UP000002311_reviewed.revCat.fasta"
MODS="$DATA_DIR/mods.txt"
THREADS=4

# Matches run_tmt_3arm.sh: -t 20ppm -ti -1,2 -m 1 (CID) -inst 1 (Orbitrap) -protocol 4 (TMT)
MSGFPLUS_ARGS_COMMON="-tda 1 -t 20ppm -ti -1,2 -m 1 -inst 1 -e 1 -protocol 4 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread $THREADS"

mkdir -p "$RESULTS_DIR"

run_msgf() {
    local JAR="$1"; local LABEL="$2"; local OUT_FMT="$3"; local OUT_EXT="$4"
    local OUT_FILE="$RESULTS_DIR/${LABEL}.${OUT_EXT}"
    local PIN="$RESULTS_DIR/${LABEL}.pin"
    local LOG="$RESULTS_DIR/${LABEL}.log"
    echo "=== $LABEL ==="
    rm -f "$DATA_DIR"/*.canno "$DATA_DIR"/*.cnlcp "$DATA_DIR"/*.csarr "$DATA_DIR"/*.cseq 2>/dev/null || true
    local START=$(python3 -c "import time; print(time.time())")
    /usr/bin/time -l java -Xmx8192m -jar "$JAR" \
        -s "$MZML" -d "$FASTA" -mod "$MODS" -o "$OUT_FILE" \
        $([[ "$OUT_FMT" == "upstream" ]] && echo "" || echo "-outputFormat $OUT_FMT") $MSGFPLUS_ARGS_COMMON > "$LOG" 2>&1 || true
    local END=$(python3 -c "import time; print(time.time())")
    local ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")
    local RSS=$(awk '/maximum resident set size/ {print $1; exit}' "$LOG" 2>/dev/null || echo 0)
    local RSS_MB=$(python3 -c "print(f'{int($RSS or 0)/1048576:.0f}')")
    if [ -f "$PIN" ]; then
        local T=$(awk -F'\t' 'NR>1 && $2==1 {c++} END {print c+0}' "$PIN")
        local D=$(awk -F'\t' 'NR>1 && $2==-1 {c++} END {print c+0}' "$PIN")
        echo "[$LABEL] wall=${ELAPSED}s peak_rss=${RSS_MB}MB targets=$T decoys=$D pin=$(wc -c < "$PIN")B"
    else
        echo "[$LABEL] FAILED"; tail -15 "$LOG"
    fi
}

run_sage() {
    local LABEL="sage"
    local OUT_DIR="$RESULTS_DIR/${LABEL}"
    local LOG="$RESULTS_DIR/${LABEL}.log"
    mkdir -p "$OUT_DIR"
    echo "=== $LABEL ==="

    local CFG_SRC="$SCRIPT_DIR/sage/tmt_config.json"
    local CFG_RUN=$(mktemp -t sage_cfg.XXXXXX).json
    sed "s|/data/|$DATA_DIR/|g; s|/out|$OUT_DIR|g" "$CFG_SRC" > "$CFG_RUN"

    local START=$(python3 -c "import time; print(time.time())")
    /usr/bin/time -l "$SAGE_BIN" "$CFG_RUN" --write-pin -o "$OUT_DIR" "$MZML" > "$LOG" 2>&1 || true
    local END=$(python3 -c "import time; print(time.time())")
    local ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")
    local RSS=$(awk '/maximum resident set size/ {print $1; exit}' "$LOG" 2>/dev/null)
    local RSS_MB=$(python3 -c "v=${RSS:-0}; print(f'{int(v)/1048576:.0f}')")
    local PEAK=$(awk '/peak memory footprint/ {print $1; exit}' "$LOG" 2>/dev/null)
    local PEAK_MB=$(python3 -c "v=${PEAK:-0}; print(f'{int(v)/1048576:.0f}')")
    local PIN=$(ls "$OUT_DIR"/*.pin 2>/dev/null | head -1)
    if [ -n "$PIN" ] && [ -f "$PIN" ]; then
        local T=$(awk -F'\t' 'NR>1 && $2==1 {c++} END {print c+0}' "$PIN")
        local D=$(awk -F'\t' 'NR>1 && $2==-1 {c++} END {print c+0}' "$PIN")
        echo "[$LABEL] wall=${ELAPSED}s rss=${RSS_MB}MB peak=${PEAK_MB}MB targets=$T decoys=$D pin=$PIN"
        cp "$PIN" "$RESULTS_DIR/${LABEL}.pin"
    else
        echo "[$LABEL] FAILED"; tail -20 "$LOG"
    fi
    rm -f "$CFG_RUN"
}

run_msgf "$BASELINE_JAR"  "A_msgfplus_baseline" upstream mzid
run_msgf "$CURRENT_JAR"   "B_msgfplus_current"  0 pin
run_sage
echo ""
ls -la "$RESULTS_DIR"/*.pin 2>/dev/null
