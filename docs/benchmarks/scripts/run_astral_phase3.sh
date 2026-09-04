#!/bin/bash
# Astral Phase 3 validation: does -useFragmentIndex on net a speedup
# on a larger dataset where the per-spectrum savings should scale?
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

DATA_DIR="$SCRIPT_DIR/data/Astral_ProteoBench_Module_8"
RESULTS_DIR="$SCRIPT_DIR/results/Astral"
BASELINE_JAR="$SCRIPT_DIR/baseline/MSGFPlus.jar"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"
MZML="$DATA_DIR/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML"
FASTA="$DATA_DIR/ProteoBenchFASTA_MixedSpecies_HYE.fasta"
MODS="$DATA_DIR/mods.txt"
THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-8192m}"
MATCHES_PER_SPEC="${MSGFPLUS_MATCHES_PER_SPEC:-1}"
BENCHMARK_MODE="${1:-${MSGFPLUS_BENCHMARK_MODE:-warm}}"

benchmark_validate_mode "$BENCHMARK_MODE"

COMMON_ARGS_NEW="-tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"

mkdir -p "$RESULTS_DIR"

run_one() {
    local JAR="$1"; local LABEL="$2"; local EXTRA_ARGS="$3"
    local OUTPUT="$RESULTS_DIR/${LABEL}_output.pin"
    local LOG="$RESULTS_DIR/${LABEL}_pin.log"
    local PREBUILD_LOG="$RESULTS_DIR/${LABEL}_prebuild.log"
    local TIME_FLAG
    [ -f "$JAR" ] || { echo "ERROR: $JAR missing"; return 1; }
    echo "=== $LABEL ($BENCHMARK_MODE) ==="
    rm -f "$OUTPUT" "$LOG" "$PREBUILD_LOG"
    if [ "$BENCHMARK_MODE" = "warm" ] && [[ "$EXTRA_ARGS" == *"-useFragmentIndex on"* ]]; then
        benchmark_clean_dataset_cache "$DATA_DIR" "$FASTA"
        benchmark_build_frag_index_cache "$JAR" "$FASTA" "$JVM_MEMORY" "$PREBUILD_LOG" \
            -mod "$MODS" -inst 3 -e 1 -minLength 6 -maxLength 40 -maxMissedCleavages 2 -ignoreMetCleavage 0 -isoforms 128
    else
        benchmark_prepare_search_cache "$BENCHMARK_MODE" "$DATA_DIR" "$FASTA" "$JAR" "$JVM_MEMORY" "$PREBUILD_LOG"
    fi
    TIME_FLAG="$(benchmark_time_flag)"
    local START=$(python3 -c "import time; print(time.time())")
    /usr/bin/time "$TIME_FLAG" java "-Xmx${JVM_MEMORY}" -jar "$JAR" \
        -s "$MZML" -d "$FASTA" -mod "$MODS" -o "$OUTPUT" \
        $COMMON_ARGS_NEW $EXTRA_ARGS > "$LOG" 2>&1 || true
    local END=$(python3 -c "import time; print(time.time())")
    local ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")
    local RSS_MB
    RSS_MB="$(benchmark_peak_rss_mb_from_log "$LOG")"
    if [ -f "$OUTPUT" ]; then
        local TARGETS=$(awk -F'\t' 'NR>1 && $2==1 {c++} END {print c+0}' "$OUTPUT")
        local DECOYS=$(awk -F'\t' 'NR>1 && $2==-1 {c++} END {print c+0}' "$OUTPUT")
        echo "[$LABEL] wall=${ELAPSED}s peak_rss=${RSS_MB}MB targets=$TARGETS decoys=$DECOYS"
    else
        echo "[$LABEL] FAILED"; tail -15 "$LOG"; return 1
    fi
}

# Skip baseline & all-off — already validated earlier
run_one "$NEW_JAR" "astralPh3_on" "-precursorCal off -useFragmentIndex on"
echo "Benchmark mode: $BENCHMARK_MODE"
ls -la "$RESULTS_DIR"/astralPh3_on*.pin 2>/dev/null
