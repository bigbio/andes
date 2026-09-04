#!/bin/bash
# Local Benchmark: Baseline (dev) vs New (feature/streaming-mzml-parser)
# Uses the test.mgf + human-uniprot-contaminants.fasta from src/test/resources
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
RESULTS_DIR="$SCRIPT_DIR/results"
BASELINE_JAR="$SCRIPT_DIR/baseline/MSGFPlus.jar"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"

SPECTRUM="$REPO_DIR/src/test/resources/test.mgf"
FASTA="$REPO_DIR/src/test/resources/human-uniprot-contaminants.fasta"
MODS="$REPO_DIR/src/test/resources/Mods.txt"
THREAD_COUNT=4
JVM_MEMORY="4096m"
SEARCH_ARGS="-tda 1 -t 20ppm -ti -1,2 -e 1 -ntt 2 -minLength 6 -maxLength 40 -minCharge 2 -maxCharge 4 -n 1 -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"

mkdir -p "$RESULTS_DIR"

run_one() {
    local JAR="$1"
    local LABEL="$2"
    local OUTPUT="$RESULTS_DIR/${LABEL}_output.mzid"
    local LOG="$RESULTS_DIR/${LABEL}.log"

    echo ""
    echo "========================================"
    echo "  $LABEL"
    echo "  JAR: $JAR"
    echo "========================================"

    # Clean suffix array index for fair comparison
    rm -f "$REPO_DIR/src/test/resources"/*.canno \
          "$REPO_DIR/src/test/resources"/*.cnlcp \
          "$REPO_DIR/src/test/resources"/*.csarr \
          "$REPO_DIR/src/test/resources"/*.cseq 2>/dev/null || true

    local START
    START=$(python3 -c "import time; print(time.time())")

    java "-Xmx${JVM_MEMORY}" -jar "$JAR" \
        -s "$SPECTRUM" \
        -d "$FASTA" \
        -mod "$MODS" \
        -o "$OUTPUT" \
        $SEARCH_ARGS \
        2>&1 | tee "$LOG"

    local END
    END=$(python3 -c "import time; print(time.time())")
    local ELAPSED
    ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")

    echo ""
    echo "[$LABEL] Wall time: ${ELAPSED}s"

    if [ ! -f "$OUTPUT" ]; then
        echo "[$LABEL] ERROR: Output not created"
        echo "label=$LABEL" > "$RESULTS_DIR/${LABEL}_metrics.txt"
        echo "wall_time=$ELAPSED" >> "$RESULTS_DIR/${LABEL}_metrics.txt"
        echo "error=output_not_created" >> "$RESULTS_DIR/${LABEL}_metrics.txt"
        return 1
    fi

    local SII_COUNT
    SII_COUNT=$(grep -c 'SpectrumIdentificationItem' "$OUTPUT" 2>/dev/null || echo 0)
    local FILE_SIZE
    FILE_SIZE=$(du -sh "$OUTPUT" | cut -f1)
    local QVALUE_COUNT
    QVALUE_COUNT=$(grep -c 'MS:1002054' "$OUTPUT" 2>/dev/null || echo 0)

    local PSM_1PCT
    PSM_1PCT=$(python3 -c "
import re
count = 0
with open('$OUTPUT') as f:
    for line in f:
        m = re.search(r'accession=\"MS:1002054\".*?value=\"([^\"]+)\"', line)
        if m and float(m.group(1)) <= 0.01:
            count += 1
print(count)
" 2>/dev/null || echo "N/A")

    echo "[$LABEL] Results:"
    echo "  Output size: $FILE_SIZE"
    echo "  SII count: $SII_COUNT"
    echo "  QValues: $QVALUE_COUNT"
    echo "  PSMs at 1% FDR: $PSM_1PCT"

    cat > "$RESULTS_DIR/${LABEL}_metrics.txt" <<EOF
label=$LABEL
wall_time=$ELAPSED
file_size=$FILE_SIZE
sii_count=$SII_COUNT
qvalue_count=$QVALUE_COUNT
psm_1pct_fdr=$PSM_1PCT
EOF
}

echo "========== LOCAL MS-GF+ BENCHMARK =========="
echo "Date: $(date)"
echo "Java: $(java -version 2>&1 | head -1)"
echo "Spectrum: $(basename "$SPECTRUM") ($(du -sh "$SPECTRUM" | cut -f1))"
echo "Database: $(basename "$FASTA") ($(du -sh "$FASTA" | cut -f1))"
echo ""

# Run baseline
run_one "$BASELINE_JAR" "baseline"

# Run new
run_one "$NEW_JAR" "new"

# Summary
echo ""
echo "========================================"
echo "  COMPARISON"
echo "========================================"

BT=$(grep '^wall_time=' "$RESULTS_DIR/baseline_metrics.txt" | cut -d= -f2)
NT=$(grep '^wall_time=' "$RESULTS_DIR/new_metrics.txt" | cut -d= -f2)
BP=$(grep '^psm_1pct_fdr=' "$RESULTS_DIR/baseline_metrics.txt" | cut -d= -f2)
NP=$(grep '^psm_1pct_fdr=' "$RESULTS_DIR/new_metrics.txt" | cut -d= -f2)
BS=$(grep '^sii_count=' "$RESULTS_DIR/baseline_metrics.txt" | cut -d= -f2)
NS=$(grep '^sii_count=' "$RESULTS_DIR/new_metrics.txt" | cut -d= -f2)

SPEEDUP=$(python3 -c "print(f'{(1-$NT/$BT)*100:.1f}%')" 2>/dev/null || echo "N/A")
PSM_CHANGE=$(python3 -c "
b, n = $BP, $NP
if b > 0: print(f'{(n/b-1)*100:+.1f}%')
else: print('N/A')
" 2>/dev/null || echo "N/A")

printf "\n%-12s %10s %10s %10s\n" "" "Time(s)" "PSMs@1%" "SII"
printf "%-12s %10s %10s %10s\n" "----------" "----------" "----------" "----------"
printf "%-12s %10s %10s %10s\n" "baseline" "${BT}s" "$BP" "$BS"
printf "%-12s %10s %10s %10s\n" "new" "${NT}s" "$NP" "$NS"
printf "%-12s %10s %10s\n" "delta" "$SPEEDUP" "$PSM_CHANGE"
echo ""
