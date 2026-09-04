#!/bin/bash
# Multi-dataset 3-way benchmark: baseline vs parser-only vs primitives-gf
#
# Datasets:
#   1. yeast  — PXD001819, CID, LTQ Orbitrap Velos, small yeast+UPS DB
#   2. tmt    — PXD007683, TMT-labeled, HCD
#   3. astral — ProteoBench Module 8, Astral DDA, mixed-species
#
# Usage:
#   bash run_local_all.sh           # Run all 3 datasets
#   bash run_local_all.sh yeast     # Run yeast only
#   bash run_local_all.sh tmt       # Run TMT only
#   bash run_local_all.sh astral    # Run Astral only
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

RESULTS_DIR="$SCRIPT_DIR/results/primitives_gf_multi"
BASELINE_JAR="$SCRIPT_DIR/baseline/MSGFPlus.jar"
PARSER_JAR="$SCRIPT_DIR/parser_only/MSGFPlus.jar"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"

THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-8192m}"
MATCHES_PER_SPEC="${MSGFPLUS_MATCHES_PER_SPEC:-1}"
BENCHMARK_MODE="${MSGFPLUS_BENCHMARK_MODE:-warm}"

benchmark_validate_mode "$BENCHMARK_MODE"

mkdir -p "$RESULTS_DIR"

# Dataset definitions
YEAST_MZML="$SCRIPT_DIR/data/PXD001819/UPS1_5000amol_R1.mzML"
YEAST_FASTA="$SCRIPT_DIR/data/PXD001819/PXD001819_uniprot_yeast_ups.fasta"
YEAST_MODS="$SCRIPT_DIR/data/PXD001819/mods.txt"
YEAST_ARGS="-tda 1 -t 0.8Da -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -tasks 0 -msLevel 2 -thread $THREAD_COUNT"

TMT_MZML="$SCRIPT_DIR/data/PXD007683/a05058.mzML"
TMT_FASTA="$SCRIPT_DIR/data/PXD007683/PXD007683_UP000005640_UP000002311_reviewed.fasta"
TMT_MODS="$SCRIPT_DIR/data/PXD007683/mods.txt"
TMT_ARGS="-tda 1 -t 10ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol 4 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -tasks 0 -msLevel 2 -thread $THREAD_COUNT"

ASTRAL_MZML="$SCRIPT_DIR/data/Astral_ProteoBench_Module_8/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML"
ASTRAL_FASTA="$SCRIPT_DIR/data/Astral_ProteoBench_Module_8/ProteoBenchFASTA_MixedSpecies_HYE.fasta"
ASTRAL_MODS="$SCRIPT_DIR/data/Astral_ProteoBench_Module_8/mods.txt"
ASTRAL_ARGS="-tda 1 -t 10ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -tasks 0 -msLevel 2 -thread $THREAD_COUNT"

run_one() {
    local JAR="$1"
    local LABEL="$2"
    local MZML="$3"
    local FASTA="$4"
    local MODS="$5"
    local EXTRA_ARGS="$6"
    local OUTPUT="$RESULTS_DIR/${LABEL}_output.mzid"
    local LOG="$RESULTS_DIR/${LABEL}.log"
    local PREBUILD_LOG="$RESULTS_DIR/${LABEL}_prebuild.log"
    local TIME_FLAG

    echo ""
    echo "========================================"
    echo "  $LABEL — $(basename "$JAR") ($BENCHMARK_MODE)"
    echo "========================================"

    if [ ! -f "$JAR" ]; then
        echo "[$LABEL] SKIP: JAR not found: $JAR"
        return 0
    fi
    if [ ! -f "$MZML" ]; then
        echo "[$LABEL] SKIP: spectrum file not found: $MZML"
        return 0
    fi

    local DATA_DIR
    DATA_DIR=$(dirname "$FASTA")
    benchmark_prepare_search_cache "$BENCHMARK_MODE" "$DATA_DIR" "$FASTA" "$JAR" "$JVM_MEMORY" "$PREBUILD_LOG"

    local START
    START=$(python3 -c "import time; print(time.time())")

    local TIME_LOG="$RESULTS_DIR/${LABEL}_time.txt"
    TIME_FLAG="$(benchmark_time_flag)"
    /usr/bin/time "$TIME_FLAG" java "-Xmx${JVM_MEMORY}" \
        "-Xlog:gc*:file=$RESULTS_DIR/${LABEL}_gc.log:time,uptime,level,tags" \
        -jar "$JAR" \
        -s "$MZML" \
        -d "$FASTA" \
        -mod "$MODS" \
        -o "$OUTPUT" \
        $EXTRA_ARGS \
        2> >(tee "$TIME_LOG" >&2) | tee "$LOG"
    local EXIT_CODE=${PIPESTATUS[0]:-0}

    local END
    END=$(python3 -c "import time; print(time.time())")
    local ELAPSED
    ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")

    local PEAK_RSS_MB=0
    if [ -f "$TIME_LOG" ]; then
        PEAK_RSS_MB="$(benchmark_peak_rss_mb_from_log "$TIME_LOG")"
    fi

    local GC_PAUSE="N/A"
    if [ -f "$RESULTS_DIR/${LABEL}_gc.log" ]; then
        GC_PAUSE=$(python3 -c "
import re
total = 0
with open('$RESULTS_DIR/${LABEL}_gc.log') as f:
    for line in f:
        m = re.search(r'Pause\s+\w+.*?([0-9.]+)ms', line)
        if m:
            total += float(m.group(1)) / 1000.0
            continue
        m = re.search(r'real=([0-9.]+) secs', line)
        if m:
            total += float(m.group(1))
print(f'{total:.2f}')
" 2>/dev/null || echo "N/A")
    fi

    echo ""
    echo "[$LABEL] Wall time: ${ELAPSED}s | Peak RSS: ${PEAK_RSS_MB}MB | GC: ${GC_PAUSE}s"

    if [ -f "$OUTPUT" ]; then
        local TOTAL_SII
        TOTAL_SII=$(grep -c 'SpectrumIdentificationItem' "$OUTPUT" 2>/dev/null || echo 0)
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
        echo "[$LABEL] SII: $TOTAL_SII, PSMs@1%FDR: $PSM_1PCT"

        cat > "$RESULTS_DIR/${LABEL}_metrics.txt" <<EOF
label=$LABEL
wall_time=$ELAPSED
peak_rss_mb=$PEAK_RSS_MB
gc_pause_s=$GC_PAUSE
sii_count=$TOTAL_SII
psm_1pct_fdr=$PSM_1PCT
exit_code=$EXIT_CODE
EOF
    else
        echo "[$LABEL] ERROR: no output!"
        cat > "$RESULTS_DIR/${LABEL}_metrics.txt" <<EOF
label=$LABEL
wall_time=$ELAPSED
peak_rss_mb=$PEAK_RSS_MB
gc_pause_s=$GC_PAUSE
error=no_output
exit_code=$EXIT_CODE
EOF
    fi
}

run_dataset() {
    local DATASET=$1
    local MZML=$2
    local FASTA=$3
    local MODS=$4
    local ARGS=$5

    echo ""
    echo "########################################"
    echo "  DATASET: $DATASET"
    echo "########################################"

    run_one "$BASELINE_JAR" "${DATASET}_baseline"    "$MZML" "$FASTA" "$MODS" "$ARGS"
    run_one "$PARSER_JAR"   "${DATASET}_parser_only" "$MZML" "$FASTA" "$MODS" "$ARGS"
    run_one "$NEW_JAR"      "${DATASET}_primitives"  "$MZML" "$FASTA" "$MODS" "$ARGS"
}

print_summary() {
    echo ""
    echo "================================================================"
    echo "  MULTI-DATASET 3-WAY SUMMARY"
    echo "================================================================"
    echo ""
    printf "%-30s %10s %10s %8s %12s %12s\n" "Variant" "Time(s)" "RSS(MB)" "GC(s)" "PSMs@1%" "SII"
    printf "%-30s %10s %10s %8s %12s %12s\n" "------------------------------" "----------" "----------" "--------" "------------" "------------"

    for f in "$RESULTS_DIR"/*_metrics.txt; do
        [ -f "$f" ] || continue
        label=$(grep '^label=' "$f" | cut -d= -f2)
        wt=$(grep '^wall_time=' "$f" | cut -d= -f2)
        rss=$(grep '^peak_rss_mb=' "$f" | cut -d= -f2)
        gc=$(grep '^gc_pause_s=' "$f" | cut -d= -f2)
        psm=$(grep '^psm_1pct_fdr=' "$f" | cut -d= -f2)
        sii=$(grep '^sii_count=' "$f" | cut -d= -f2)
        printf "%-30s %10s %10s %8s %12s %12s\n" "$label" "${wt}s" "${rss}MB" "${gc}s" "$psm" "$sii"
    done

    echo ""
    echo "--- Per-dataset speedup breakdown ---"
    echo ""

    for dataset in yeast tmt astral; do
        bf="$RESULTS_DIR/${dataset}_baseline_metrics.txt"
        pf="$RESULTS_DIR/${dataset}_parser_only_metrics.txt"
        nf="$RESULTS_DIR/${dataset}_primitives_metrics.txt"
        if [ -f "$bf" ] && [ -f "$pf" ] && [ -f "$nf" ]; then
            bt=$(grep '^wall_time=' "$bf" | cut -d= -f2)
            pt=$(grep '^wall_time=' "$pf" | cut -d= -f2)
            nt=$(grep '^wall_time=' "$nf" | cut -d= -f2)
            br=$(grep '^peak_rss_mb=' "$bf" | cut -d= -f2)
            pr=$(grep '^peak_rss_mb=' "$pf" | cut -d= -f2)
            nr=$(grep '^peak_rss_mb=' "$nf" | cut -d= -f2)
            bp=$(grep '^psm_1pct_fdr=' "$bf" | cut -d= -f2)
            np=$(grep '^psm_1pct_fdr=' "$nf" | cut -d= -f2)

            python3 -c "
bt, pt, nt = $bt, $pt, $nt
br, pr, nr = int('${br:-0}' or 0), int('${pr:-0}' or 0), int('${nr:-0}' or 0)
bp, np = '${bp}', '${np}'
parser_sp = (1 - pt/bt) * 100 if bt > 0 else 0
gf_sp = (1 - nt/pt) * 100 if pt > 0 else 0
total_sp = (1 - nt/bt) * 100 if bt > 0 else 0
mem_change = (nr - pr) / pr * 100 if pr > 0 else 0
print(f'  $dataset:')
print(f'    Parser:     {parser_sp:+.1f}%  ({bt}s -> {pt}s)')
print(f'    GF prims:   {gf_sp:+.1f}%  ({pt}s -> {nt}s)')
print(f'    Combined:   {total_sp:+.1f}%  ({bt}s -> {nt}s)')
print(f'    Memory:     {pr}MB -> {nr}MB ({mem_change:+.1f}%)')
print(f'    PSMs@1%FDR: {bp} -> {np} (same={bp==np})')
print()
"
        fi
    done
}

# Main
FILTER="${1:-all}"

echo "========== MS-GF+ MULTI-DATASET 3-WAY BENCHMARK =========="
echo "Date: $(date)"
echo "Java: $(java -version 2>&1 | head -1)"
echo "CPUs: $(sysctl -n hw.ncpu 2>/dev/null || nproc) (using $THREAD_COUNT threads)"
echo "JVM max heap: $JVM_MEMORY"
echo "Benchmark mode: $BENCHMARK_MODE"

if [ "$FILTER" = "all" ] || [ "$FILTER" = "yeast" ]; then
    run_dataset "yeast" "$YEAST_MZML" "$YEAST_FASTA" "$YEAST_MODS" "$YEAST_ARGS"
fi

if [ "$FILTER" = "all" ] || [ "$FILTER" = "tmt" ]; then
    run_dataset "tmt" "$TMT_MZML" "$TMT_FASTA" "$TMT_MODS" "$TMT_ARGS"
fi

if [ "$FILTER" = "all" ] || [ "$FILTER" = "astral" ]; then
    run_dataset "astral" "$ASTRAL_MZML" "$ASTRAL_FASTA" "$ASTRAL_MODS" "$ASTRAL_ARGS"
fi

print_summary
