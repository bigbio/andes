#!/bin/bash
# Benchmark: Baseline vs New MS-GF+ on Orbitrap Astral (ProteoBench Module 8)
#
# Dataset: ProteoBench Module 8 — LFQ, HCD, Orbitrap Astral, DDA 15min
#   Spectra: https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/lfq/OrbitrapAstral/ProteoBench_Module_8/
#   FASTA:   https://github.com/bigbio/quantms-test-datasets (databases/ProteoBenchFASTA_MixedSpecies_HYE.fasta)
#   SDRF:    https://github.com/bigbio/quantms-test-datasets (benchmarks/lfq/OrbitrapAstral/noPXD_ProteoBench_Module_8/)
#
# Search parameters extracted from SDRF:
#   Enzyme: Trypsin, Precursor: 10ppm, Fragment: 20ppm, HCD
#   Fixed: Carbamidomethyl (C), Variable: Oxidation (M), Acetyl (Protein N-term)
#   Database: Mixed-species HYE (Human 20537 + Yeast 6722 + E.coli 4401 + Contaminants 381 = 31889 proteins)
#
# Usage:
#   bash run_astral_benchmark.sh              # Run both baseline + new
#   bash run_astral_benchmark.sh new-only      # Run only the new JAR
#   bash run_astral_benchmark.sh baseline-only  # Run only the baseline JAR
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data/Astral_ProteoBench_Module_8"
RESULTS_DIR="$SCRIPT_DIR/results/Astral"
BASELINE_JAR="$SCRIPT_DIR/baseline/MSGFPlus.jar"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"

MZML="$DATA_DIR/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML"
FASTA="$DATA_DIR/ProteoBenchFASTA_MixedSpecies_HYE.fasta"
MODS="$DATA_DIR/mods.txt"

THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-8192m}"

# Search parameters from SDRF:
#   Trypsin (-e 1), HCD (-m 3), Orbitrap Astral high-res MS2 (-inst 3 = high-res)
#   Precursor: 10ppm, Fragment: 20ppm
SEARCH_ARGS="-tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"

mkdir -p "$RESULTS_DIR"

for f in "$MZML" "$FASTA" "$MODS"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        echo "Run download_data.sh first, or decompress the .mzML.gz file."
        exit 1
    fi
done

run_one() {
    local JAR="$1"
    local LABEL="$2"
    local OUTPUT="$RESULTS_DIR/${LABEL}_output.mzid"
    local LOG="$RESULTS_DIR/${LABEL}.log"

    if [ ! -f "$JAR" ]; then
        echo "ERROR: JAR not found: $JAR"
        return 1
    fi

    echo ""
    echo "========================================"
    echo "  $LABEL"
    echo "  JAR: $JAR"
    echo "  Spectrum: $(basename "$MZML") ($(du -sh "$MZML" | cut -f1))"
    echo "  Database: $(basename "$FASTA") ($(du -sh "$FASTA" | cut -f1))"
    echo "========================================"

    # Clean suffix array index for fair comparison
    rm -f "$DATA_DIR"/*.canno "$DATA_DIR"/*.cnlcp "$DATA_DIR"/*.csarr "$DATA_DIR"/*.cseq 2>/dev/null || true

    local START
    START=$(python3 -c "import time; print(time.time())")

    /usr/bin/time -l java "-Xmx${JVM_MEMORY}" -jar "$JAR" \
        -s "$MZML" \
        -d "$FASTA" \
        -mod "$MODS" \
        -o "$OUTPUT" \
        $SEARCH_ARGS \
        > "$LOG" 2>&1 || true

    local END
    END=$(python3 -c "import time; print(time.time())")
    local ELAPSED
    ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")

    local PEAK_RSS_BYTES
    PEAK_RSS_BYTES=$(awk '/maximum resident set size/ {print $1; exit}' "$LOG" 2>/dev/null)
    PEAK_RSS_BYTES=${PEAK_RSS_BYTES:-0}
    local PEAK_RSS_MB
    PEAK_RSS_MB=$(python3 -c "print(f'{int($PEAK_RSS_BYTES)/1048576:.0f}')" 2>/dev/null || echo "N/A")

    echo ""
    echo "[$LABEL] Wall time: ${ELAPSED}s"
    echo "[$LABEL] Peak RSS: ${PEAK_RSS_MB} MB"

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
    local RAWSCORE_COUNT
    RAWSCORE_COUNT=$(grep -c 'MS:1002049' "$OUTPUT" 2>/dev/null || echo 0)

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
    echo "  RawScores: $RAWSCORE_COUNT"
    echo "  QValues: $QVALUE_COUNT"
    echo "  PSMs at 1% FDR: $PSM_1PCT"

    cat > "$RESULTS_DIR/${LABEL}_metrics.txt" <<EOF
label=$LABEL
wall_time=$ELAPSED
peak_rss_mb=$PEAK_RSS_MB
file_size=$FILE_SIZE
sii_count=$SII_COUNT
rawscore_count=$RAWSCORE_COUNT
qvalue_count=$QVALUE_COUNT
psm_1pct_fdr=$PSM_1PCT
EOF
}

# ============================================================
# Main
# ============================================================

FILTER="${1:-all}"

echo "========== ASTRAL BENCHMARK (ProteoBench Module 8) =========="
echo "Date: $(date)"
echo "Java: $(java -version 2>&1 | head -1)"
echo "CPUs: $(sysctl -n hw.ncpu 2>/dev/null || nproc) (using $THREAD_COUNT threads)"
echo "Memory: $(sysctl -n hw.memsize 2>/dev/null | awk '{printf "%.0fG", $1/1073741824}' || free -h | awk '/Mem:/{print $2}')"
echo "mzML: $(basename "$MZML") ($(du -sh "$MZML" 2>/dev/null | cut -f1 || echo 'not found'))"
echo "FASTA: $(basename "$FASTA") ($(du -sh "$FASTA" | cut -f1))"
echo "Database: Mixed-species HYE — H(20537) + Y(6722) + E(4401) + Cont(381) = 31889 proteins"
echo ""

if [ "$FILTER" != "new-only" ]; then
    run_one "$BASELINE_JAR" "baseline"
fi

if [ "$FILTER" != "baseline-only" ]; then
    run_one "$NEW_JAR" "new"
fi

# Comparison
if [ -f "$RESULTS_DIR/baseline_metrics.txt" ] && [ -f "$RESULTS_DIR/new_metrics.txt" ]; then
    echo ""
    echo "========================================"
    echo "  COMPARISON"
    echo "========================================"

    BT=$(grep '^wall_time=' "$RESULTS_DIR/baseline_metrics.txt" | cut -d= -f2)
    NT=$(grep '^wall_time=' "$RESULTS_DIR/new_metrics.txt" | cut -d= -f2)
    BR=$(grep '^peak_rss_mb=' "$RESULTS_DIR/baseline_metrics.txt" | cut -d= -f2)
    NR=$(grep '^peak_rss_mb=' "$RESULTS_DIR/new_metrics.txt" | cut -d= -f2)
    BP=$(grep '^psm_1pct_fdr=' "$RESULTS_DIR/baseline_metrics.txt" | cut -d= -f2)
    NP=$(grep '^psm_1pct_fdr=' "$RESULTS_DIR/new_metrics.txt" | cut -d= -f2)
    BS=$(grep '^sii_count=' "$RESULTS_DIR/baseline_metrics.txt" | cut -d= -f2)
    NS=$(grep '^sii_count=' "$RESULTS_DIR/new_metrics.txt" | cut -d= -f2)

    SPEEDUP=$(python3 -c "print(f'{(1-$NT/$BT)*100:.1f}%')" 2>/dev/null || echo "N/A")
    RSS_CHANGE=$(python3 -c "
b, n = $BR, $NR
if b > 0: print(f'{(n/b-1)*100:+.1f}%')
else: print('N/A')
" 2>/dev/null || echo "N/A")
    PSM_CHANGE=$(python3 -c "
b, n = $BP, $NP
if b > 0: print(f'{(n/b-1)*100:+.1f}%')
else: print('N/A')
" 2>/dev/null || echo "N/A")

    printf "\n%-12s %10s %10s %10s %10s\n" "" "Time(s)" "RSS(MB)" "PSMs@1%" "SII"
    printf "%-12s %10s %10s %10s %10s\n" "----------" "----------" "----------" "----------" "----------"
    printf "%-12s %10s %10s %10s %10s\n" "baseline" "${BT}s" "${BR}" "$BP" "$BS"
    printf "%-12s %10s %10s %10s %10s\n" "new" "${NT}s" "${NR}" "$NP" "$NS"
    printf "%-12s %10s %10s %10s\n" "delta" "$SPEEDUP" "$RSS_CHANGE" "$PSM_CHANGE"
    echo ""
fi
