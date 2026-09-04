#!/bin/bash
# Phase 3 validation benchmark on PXD001819.
# Arms:
#   A  = baseline JAR
#   B  = branch JAR + -useFragmentIndex off (correctness gate — must bit-match A)
#   D  = branch JAR + -useFragmentIndex on -precursorCal auto (full new pipeline)
#
# Skips arm C (P2-cal only) since Achievement B was already validated separately.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

DATA_DIR="$SCRIPT_DIR/data/PXD001819"
RESULTS_DIR="$SCRIPT_DIR/results/PXD001819"
BASELINE_JAR="$SCRIPT_DIR/baseline-fork-pre-pr22/MSGFPlus.jar"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"
MZML="$DATA_DIR/UPS1_5000amol_R1.mzML"
FASTA="$DATA_DIR/PXD001819_uniprot_yeast_ups.fasta"
MODS="$DATA_DIR/mods.txt"

THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-4096m}"
MATCHES_PER_SPEC="${MSGFPLUS_MATCHES_PER_SPEC:-1}"
BENCHMARK_MODE="${1:-${MSGFPLUS_BENCHMARK_MODE:-warm}}"
COMPUTE_TOPK_RECALL="${MSGFPLUS_COMPUTE_TOPK_RECALL:-0}"

benchmark_validate_mode "$BENCHMARK_MODE"

COMMON_ARGS_BASELINE="-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -msLevel 2 -outputFormat 3 -thread $THREAD_COUNT"
COMMON_ARGS_NEW="-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"

mkdir -p "$RESULTS_DIR"

run_one() {
    local JAR="$1"
    local KIND="$2"
    local LABEL="$3"
    local EXTRA_ARGS="$4"
    local OUTPUT_FILE
    local OUTPUT
    local LOG="$RESULTS_DIR/${LABEL}_pin.log"
    local PREBUILD_LOG="$RESULTS_DIR/${LABEL}_prebuild.log"
    local TIME_FLAG
    local SEARCH_ARGS

    [ -f "$JAR" ] || { echo "ERROR: $JAR missing"; return 1; }

    case "$KIND" in
        baseline)
            OUTPUT_FILE="$RESULTS_DIR/${LABEL}_output.mzid"
            OUTPUT="$RESULTS_DIR/${LABEL}_output.pin"
            SEARCH_ARGS="$COMMON_ARGS_BASELINE"
            ;;
        current)
            OUTPUT_FILE="$RESULTS_DIR/${LABEL}_output.pin"
            OUTPUT="$OUTPUT_FILE"
            SEARCH_ARGS="$COMMON_ARGS_NEW"
            ;;
        *)
            echo "ERROR: unknown MSGF+ runner kind: $KIND" >&2
            return 1
            ;;
    esac

    echo "=== $LABEL ($BENCHMARK_MODE) ==="
    rm -f "$OUTPUT_FILE" "$OUTPUT" "$LOG" "$PREBUILD_LOG"
    if [ "$BENCHMARK_MODE" = "warm" ] && [ "$KIND" = "current" ] && [[ "$EXTRA_ARGS" == *"-useFragmentIndex on"* ]]; then
        benchmark_clean_dataset_cache "$DATA_DIR" "$FASTA"
        benchmark_build_frag_index_cache "$JAR" "$FASTA" "$JVM_MEMORY" "$PREBUILD_LOG" \
            -mod "$MODS" -inst 0 -e 1 -minLength 6 -maxLength 40 -maxMissedCleavages 2 -ignoreMetCleavage 0 -isoforms 128
    else
        benchmark_prepare_search_cache "$BENCHMARK_MODE" "$DATA_DIR" "$FASTA" "$JAR" "$JVM_MEMORY" "$PREBUILD_LOG"
    fi
    TIME_FLAG="$(benchmark_time_flag)"

    local START=$(python3 -c "import time; print(time.time())")
    /usr/bin/time "$TIME_FLAG" java "-Xmx${JVM_MEMORY}" -jar "$JAR" \
        -s "$MZML" -d "$FASTA" -mod "$MODS" -o "$OUTPUT_FILE" \
        $SEARCH_ARGS $EXTRA_ARGS > "$LOG" 2>&1 || true
    local END=$(python3 -c "import time; print(time.time())")
    local ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")
    local RSS_MB
    RSS_MB="$(benchmark_peak_rss_mb_from_log "$LOG")"
    if [ -f "$OUTPUT" ]; then
        local TARGETS=$(awk -F'\t' 'NR>1 && $2==1 {c++} END {print c+0}' "$OUTPUT")
        local DECOYS=$(awk -F'\t' 'NR>1 && $2==-1 {c++} END {print c+0}' "$OUTPUT")
        echo "[$LABEL] wall=${ELAPSED}s peak_rss=${RSS_MB}MB targets=$TARGETS decoys=$DECOYS"
    else
        echo "[$LABEL] FAILED"; tail -20 "$LOG"
        return 1
    fi
}

run_one "$BASELINE_JAR" "baseline" "ph3A_baseline" ""
run_one "$NEW_JAR" "current" "ph3B_allOff" "-precursorCal off -useFragmentIndex off"
run_one "$NEW_JAR" "current" "ph3D_AplusBplusC" "-precursorCal auto -useFragmentIndex on"
echo ""
echo "Benchmark mode: $BENCHMARK_MODE"
ls -la "$RESULTS_DIR"/ph3*.pin 2>/dev/null

if [ "$COMPUTE_TOPK_RECALL" = "1" ] && \
   [ -f "$RESULTS_DIR/ph3A_baseline_output.pin" ] && \
   [ -f "$RESULTS_DIR/ph3D_AplusBplusC_output.pin" ]; then
    echo ""
    echo "=== Baseline top-1 recovered in ON-mode top-$MATCHES_PER_SPEC ==="
    python3 "$SCRIPT_DIR/compare_pin_topk_recall.py" \
        "$RESULTS_DIR/ph3A_baseline_output.pin" \
        "$RESULTS_DIR/ph3D_AplusBplusC_output.pin" \
        --baseline-rank 1 \
        --experimental-rank "$MATCHES_PER_SPEC" \
        --only-label 1
fi
