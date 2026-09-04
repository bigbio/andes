#!/bin/bash
# Three-arm Orbitrap Astral benchmark — the real test of Achievement B.
# Astral often has 5-10 ppm systematic precursor bias → P2-cal should deliver
# a larger lift than the already-well-calibrated PXD001819 case.
#
# Arms:
#   A  = baseline JAR
#   B  = branch JAR + -precursorCal off (Achievement A only)
#   C  = branch JAR + -precursorCal auto (Achievement A + B)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

DATA_DIR="$SCRIPT_DIR/data/Astral_ProteoBench_Module_8"
RESULTS_DIR="$SCRIPT_DIR/results/Astral"
BASELINE_JAR="$SCRIPT_DIR/baseline-fork-pre-pr22/MSGFPlus.jar"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"
MZML="$DATA_DIR/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML"
FASTA="$DATA_DIR/ProteoBenchFASTA_MixedSpecies_HYE.fasta"
MODS="$DATA_DIR/mods.txt"

THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-8192m}"
MATCHES_PER_SPEC="${MSGFPLUS_MATCHES_PER_SPEC:-1}"
BENCHMARK_MODE="${1:-${MSGFPLUS_BENCHMARK_MODE:-cold}}"

benchmark_validate_mode "$BENCHMARK_MODE"

COMMON_ARGS_BASELINE="-tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -msLevel 2 -outputFormat 3 -thread $THREAD_COUNT"
COMMON_ARGS_NEW="-tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"

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
    benchmark_prepare_search_cache "$BENCHMARK_MODE" "$DATA_DIR" "$FASTA" "$JAR" "$JVM_MEMORY" "$PREBUILD_LOG"
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
        echo "[$LABEL] FAILED"; tail -15 "$LOG"
        return 1
    fi
}

run_one "$BASELINE_JAR" "baseline" "astralA_baseline" ""
run_one "$NEW_JAR" "current" "astralB_AonlyOff" "-precursorCal off"
run_one "$NEW_JAR" "current" "astralC_AplusB" "-precursorCal auto"

echo ""
echo "Benchmark mode: $BENCHMARK_MODE"
echo "=== Astral 3-arm comparison ==="
ls -la "$RESULTS_DIR"/astral[ABC]*.pin 2>/dev/null
