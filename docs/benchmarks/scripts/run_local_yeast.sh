#!/bin/bash
# Local single-file benchmark: baseline vs parser-only vs primitives-gf on yeast (PXD001819)
#
# Three-way comparison to isolate the GF primitives speedup from the mzML parser speedup:
#   1. baseline     = original master (no parser, no GF changes)
#   2. parser_only  = streaming mzML parser only
#   3. primitives   = parser + primitive GF arrays
#
# Also tracks peak memory (RSS) via a background sampler.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data/PXD001819"
RESULTS_DIR="$SCRIPT_DIR/results/primitives_gf"
BASELINE_JAR="$SCRIPT_DIR/baseline/MSGFPlus.jar"
PARSER_JAR="$SCRIPT_DIR/parser_only/MSGFPlus.jar"
NEW_JAR="$SCRIPT_DIR/new/MSGFPlus.jar"

MZML="$DATA_DIR/UPS1_5000amol_R1.mzML"
FASTA="$DATA_DIR/PXD001819_uniprot_yeast_ups.fasta"
MODS="$DATA_DIR/mods.txt"

THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-8192m}"
COMMON_ARGS="-tda 1 -t 0.8Da -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -tasks 0 -msLevel 2 -thread $THREAD_COUNT"

mkdir -p "$RESULTS_DIR"

run_one() {
    local JAR="$1"
    local LABEL="$2"
    local OUTPUT="$RESULTS_DIR/${LABEL}_output.mzid"
    local LOG="$RESULTS_DIR/${LABEL}.log"

    echo ""
    echo "========================================"
    echo "  $LABEL — $(basename "$JAR")"
    echo "========================================"

    if [ ! -f "$JAR" ]; then
        echo "[$LABEL] SKIP: JAR not found: $JAR"
        return 0
    fi

    # Clean SA index for fair comparison
    rm -f "$DATA_DIR"/*.canno "$DATA_DIR"/*.cnlcp "$DATA_DIR"/*.csarr "$DATA_DIR"/*.cseq 2>/dev/null || true

    local START
    START=$(python3 -c "import time; print(time.time())")

    # Use /usr/bin/time for reliable peak RSS measurement
    local TIME_LOG="$RESULTS_DIR/${LABEL}_time.txt"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        /usr/bin/time -l java "-Xmx${JVM_MEMORY}" \
            "-Xlog:gc*:file=$RESULTS_DIR/${LABEL}_gc.log:time,uptime,level,tags" \
            -jar "$JAR" \
            -s "$MZML" \
            -d "$FASTA" \
            -mod "$MODS" \
            -o "$OUTPUT" \
            $COMMON_ARGS \
            2> >(tee "$TIME_LOG" >&2) | tee "$LOG"
    else
        /usr/bin/time -v java "-Xmx${JVM_MEMORY}" \
            "-Xlog:gc*:file=$RESULTS_DIR/${LABEL}_gc.log:time,uptime,level,tags" \
            -jar "$JAR" \
            -s "$MZML" \
            -d "$FASTA" \
            -mod "$MODS" \
            -o "$OUTPUT" \
            $COMMON_ARGS \
            2> >(tee "$TIME_LOG" >&2) | tee "$LOG"
    fi
    local EXIT_CODE=${PIPESTATUS[0]:-0}

    local END
    END=$(python3 -c "import time; print(time.time())")
    local ELAPSED
    ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")

    # Read peak RSS from /usr/bin/time output
    local PEAK_RSS_MB=0
    local TIME_LOG="$RESULTS_DIR/${LABEL}_time.txt"
    if [ -f "$TIME_LOG" ]; then
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS: "maximum resident set size" in bytes
            local RSS_BYTES
            RSS_BYTES=$(grep 'maximum resident set size' "$TIME_LOG" | awk '{print $1}' | head -1)
            PEAK_RSS_MB=$(python3 -c "print(f'{int(${RSS_BYTES:-0})/1048576:.0f}')" 2>/dev/null || echo 0)
        else
            # Linux: "Maximum resident set size (kbytes)"
            local RSS_KB
            RSS_KB=$(grep 'Maximum resident set size' "$TIME_LOG" | awk '{print $NF}' | head -1)
            PEAK_RSS_MB=$(python3 -c "print(f'{int(${RSS_KB:-0})/1024:.0f}')" 2>/dev/null || echo 0)
        fi
    fi

    # Parse GC stats
    local GC_PAUSE="N/A"
    if [ -f "$RESULTS_DIR/${LABEL}_gc.log" ]; then
        GC_PAUSE=$(python3 -c "
import re
total = 0
with open('$RESULTS_DIR/${LABEL}_gc.log') as f:
    for line in f:
        # JDK17 unified logging: 'Pause Young ... 12.345ms' or 'real=X.XX secs' (JDK8)
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
    echo "[$LABEL] Wall time: ${ELAPSED}s | Peak RSS: ${PEAK_RSS_MB}MB | GC pause: ${GC_PAUSE}s"

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

echo "========== PRIMITIVES-GF BENCHMARK (yeast PXD001819) =========="
echo "Date: $(date)"
echo "Java: $(java -version 2>&1 | head -1)"
echo "CPUs: $(sysctl -n hw.ncpu 2>/dev/null || nproc) (using $THREAD_COUNT threads)"
echo "JVM max heap: $JVM_MEMORY"
echo ""

# Run all three JARs
run_one "$BASELINE_JAR" "baseline"
run_one "$PARSER_JAR"   "parser_only"
run_one "$NEW_JAR"      "primitives"

# Summary
echo ""
echo "================================================================"
echo "  3-WAY COMPARISON (yeast PXD001819)"
echo "================================================================"
echo ""
printf "%-15s %10s %10s %10s %12s %12s\n" "Variant" "Time(s)" "RSS(MB)" "GC(s)" "PSMs@1%" "SII"
printf "%-15s %10s %10s %10s %12s %12s\n" "---------------" "----------" "----------" "----------" "------------" "------------"

for label in baseline parser_only primitives; do
    f="$RESULTS_DIR/${label}_metrics.txt"
    if [ -f "$f" ]; then
        wt=$(grep '^wall_time=' "$f" | cut -d= -f2)
        rss=$(grep '^peak_rss_mb=' "$f" | cut -d= -f2)
        gc=$(grep '^gc_pause_s=' "$f" | cut -d= -f2)
        psm=$(grep '^psm_1pct_fdr=' "$f" | cut -d= -f2)
        sii=$(grep '^sii_count=' "$f" | cut -d= -f2)
        printf "%-15s %10s %10s %10s %12s %12s\n" "$label" "${wt}s" "${rss}MB" "${gc}s" "$psm" "$sii"
    fi
done

echo ""
echo "--- Speedup breakdown ---"

BT=$(grep '^wall_time=' "$RESULTS_DIR/baseline_metrics.txt" 2>/dev/null | cut -d= -f2 || echo "0")
PT=$(grep '^wall_time=' "$RESULTS_DIR/parser_only_metrics.txt" 2>/dev/null | cut -d= -f2 || echo "0")
NT=$(grep '^wall_time=' "$RESULTS_DIR/primitives_metrics.txt" 2>/dev/null | cut -d= -f2 || echo "0")

python3 -c "
bt, pt, nt = $BT, $PT, $NT
if bt > 0 and pt > 0 and nt > 0:
    parser_speedup = (1 - pt/bt) * 100
    gf_speedup = (1 - nt/pt) * 100
    total_speedup = (1 - nt/bt) * 100
    print(f'  Parser improvement (baseline -> parser_only): {parser_speedup:+.1f}%  ({bt}s -> {pt}s)')
    print(f'  GF primitives     (parser_only -> primitives): {gf_speedup:+.1f}%  ({pt}s -> {nt}s)')
    print(f'  Combined          (baseline -> primitives):    {total_speedup:+.1f}%  ({bt}s -> {nt}s)')
    print()
    parser_delta = bt - pt
    gf_delta = pt - nt
    total_delta = bt - nt
    print(f'  Parser saved:  {parser_delta:.1f}s')
    print(f'  GF saved:      {gf_delta:.1f}s')
    print(f'  Total saved:   {total_delta:.1f}s')
else:
    print('  (insufficient data for breakdown)')
"

# Memory comparison
echo ""
echo "--- Memory comparison ---"
BR=$(grep '^peak_rss_mb=' "$RESULTS_DIR/baseline_metrics.txt" 2>/dev/null | cut -d= -f2 || echo "0")
PR=$(grep '^peak_rss_mb=' "$RESULTS_DIR/parser_only_metrics.txt" 2>/dev/null | cut -d= -f2 || echo "0")
NR=$(grep '^peak_rss_mb=' "$RESULTS_DIR/primitives_metrics.txt" 2>/dev/null | cut -d= -f2 || echo "0")

python3 -c "
br, pr, nr = int('${BR:-0}' or 0), int('${PR:-0}' or 0), int('${NR:-0}' or 0)
if br > 0 and pr > 0 and nr > 0:
    print(f'  Baseline peak RSS:     {br}MB')
    print(f'  Parser-only peak RSS:  {pr}MB  ({(pr-br)/br*100:+.1f}% vs baseline)')
    print(f'  Primitives peak RSS:   {nr}MB  ({(nr-br)/br*100:+.1f}% vs baseline, {(nr-pr)/pr*100:+.1f}% vs parser)')
else:
    print('  (insufficient memory data)')
"
echo ""
