#!/bin/bash
# Multi-Dataset Benchmark: Baseline vs New MS-GF+
#
# Tests three experiment types to cover different instruments, fragmentation,
# databases, and labeling strategies:
#
#   1. TMT (Human)         — TMT-labeled, HCD, large DB (41K proteins)
#   2. LFQ (Yeast)         — PXD001819, CID, low-res LTQ, small yeast DB (6.8K proteins)
#   3. LFQ (Mixed-species) — PXD028735, HCD, Q Exactive HF, medium DB (31.9K proteins)
#
# Data sources:
#   Spectra: https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/
#   FASTA/SDRF: https://github.com/bigbio/quantms-test-datasets
#
# Usage:
#   bash run_benchmark_multi.sh              # Run all 3 datasets
#   bash run_benchmark_multi.sh tmt          # Run TMT only
#   bash run_benchmark_multi.sh yeast        # Run yeast only
#   bash run_benchmark_multi.sh lfq          # Run LFQ mixed-species only
#   bash run_benchmark_multi.sh new-only     # Run only the new JAR (skip baseline)
#
set -euo pipefail

# ============================================================
# Configuration — adjust paths to match your environment
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Override BENCHMARK_DIR via env var, or default to /srv path for remote server
BENCHMARK_DIR="${MSGFPLUS_BENCHMARK_DIR:-/srv/data/msgf_improvements/benchmark}"
RESULTS_DIR="$BENCHMARK_DIR/results/multi"
BASELINE_JAR="${MSGFPLUS_BASELINE_JAR:-$SCRIPT_DIR/baseline/MSGFPlus.jar}"
NEW_JAR="${MSGFPLUS_NEW_JAR:-$SCRIPT_DIR/new/MSGFPlus.jar}"
THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEMORY="${MSGFPLUS_MEMORY:-8192m}"

mkdir -p "$RESULTS_DIR"

# ============================================================
# Dataset definitions
# ============================================================

# Dataset 1: TMT Human (current benchmark)
TMT_MZML="$BENCHMARK_DIR/data/f07531_Prot_02_R2_F08.mzML"
TMT_FASTA="$BENCHMARK_DIR/data/Homo-sapiens-target-only.fasta"
TMT_MODS="$BENCHMARK_DIR/data/mods.txt"
TMT_ARGS="-tda 1 -t 30ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 4 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -tasks 0 -msLevel 2 -thread $THREAD_COUNT"

# Dataset 2: Yeast + UPS1 spike-in (PXD001819) — CID, low-res
# Spectra: https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/lfq/LTQOrbitrapVelos/PXD001819/
# FASTA: https://github.com/bigbio/quantms-test-datasets (benchmarks/lfq/PXD001819/)
YEAST_MZML="${BENCHMARK_DIR}/data/PXD001819/UPS1_5000amol_R1.mzML"
YEAST_FASTA="${BENCHMARK_DIR}/data/PXD001819/yeast_ups1_reviewed.fasta"
YEAST_MODS="$BENCHMARK_DIR/data/mods_lfq.txt"
YEAST_ARGS="-tda 1 -t 0.8Da -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -tasks 0 -msLevel 2 -thread $THREAD_COUNT"

# Dataset 3: LFQ mixed-species ProteoBench (PXD028735) — HCD, Q Exactive HF
# Spectra: https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/
# FASTA: https://github.com/bigbio/quantms-test-datasets (benchmarks/lfq/PXD028735/)
LFQ_MZML="${BENCHMARK_DIR}/data/PXD028735/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML"
LFQ_FASTA="${BENCHMARK_DIR}/data/PXD028735/ProteoBenchFASTA_MixedSpecies_HYE.fasta"
LFQ_MODS="$BENCHMARK_DIR/data/mods_lfq.txt"
LFQ_ARGS="-tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -tasks 0 -msLevel 2 -thread $THREAD_COUNT"

# ============================================================
# Helper functions
# ============================================================

run_msgfplus() {
    local JAR="$1"
    local LABEL="$2"
    local MZML="$3"
    local FASTA="$4"
    local MODS="$5"
    local EXTRA_ARGS="$6"
    local OUTPUT="$RESULTS_DIR/${LABEL}_output.mzid"
    local LOG="$RESULTS_DIR/${LABEL}.log"

    echo ""
    echo "========================================"
    echo "  $LABEL"
    echo "  JAR: $(basename "$JAR")"
    echo "  Spectrum: $(basename "$MZML")"
    echo "  Database: $(basename "$FASTA")"
    echo "========================================"

    if [ ! -f "$JAR" ]; then
        echo "[$LABEL] ERROR: JAR not found: $JAR"
        return 1
    fi
    if [ ! -f "$MZML" ]; then
        echo "[$LABEL] ERROR: Spectrum file not found: $MZML"
        echo "  Run download_data.sh or check MSGFPLUS_BENCHMARK_DIR"
        return 1
    fi
    if [ ! -f "$FASTA" ]; then
        echo "[$LABEL] ERROR: FASTA file not found: $FASTA"
        return 1
    fi

    # Clean up suffix array index files for fair comparison
    local DATA_DIR
    DATA_DIR=$(dirname "$FASTA")
    rm -f "$DATA_DIR"/*.canno "$DATA_DIR"/*.cnlcp "$DATA_DIR"/*.csarr "$DATA_DIR"/*.cseq 2>/dev/null || true

    local START
    START=$(date +%s)

    java "-Xmx${JVM_MEMORY}" -jar "$JAR" \
        -s "$MZML" \
        -d "$FASTA" \
        -mod "$MODS" \
        -o "$OUTPUT" \
        $EXTRA_ARGS \
        2>&1 | tee "$LOG"

    local END
    END=$(date +%s)
    local ELAPSED=$(( END - START ))

    echo ""
    echo "[$LABEL] Wall time: ${ELAPSED}s"

    if [ ! -f "$OUTPUT" ]; then
        echo "[$LABEL] ERROR: Output file not created!"
        cat > "$RESULTS_DIR/${LABEL}_metrics.txt" <<EOF
label=$LABEL
wall_time=${ELAPSED}
error=output_not_created
EOF
        return 1
    fi

    local TOTAL_SII
    TOTAL_SII=$(grep -c 'SpectrumIdentificationItem' "$OUTPUT" 2>/dev/null || echo 0)
    local TOTAL_RAWSCORE
    TOTAL_RAWSCORE=$(grep -c 'MS:1002049' "$OUTPUT" 2>/dev/null || echo 0)
    local TOTAL_DENOVO
    TOTAL_DENOVO=$(grep -c 'MS:1002050' "$OUTPUT" 2>/dev/null || echo 0)
    local TOTAL_SPECEVALUE
    TOTAL_SPECEVALUE=$(grep -c 'MS:1002052' "$OUTPUT" 2>/dev/null || echo 0)
    local TOTAL_EVALUE
    TOTAL_EVALUE=$(grep -c 'MS:1002053' "$OUTPUT" 2>/dev/null || echo 0)
    local TOTAL_QVALUE
    TOTAL_QVALUE=$(grep -c 'MS:1002054' "$OUTPUT" 2>/dev/null || echo 0)
    local FILE_SIZE
    FILE_SIZE=$(du -sh "$OUTPUT" | cut -f1)

    local PSM_1PCT_FDR
    PSM_1PCT_FDR=$(python3 -c "
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
    echo "  Output file size: $FILE_SIZE"
    echo "  SpectrumIdentificationItem count: $TOTAL_SII"
    echo "  RawScore (MS:1002049): $TOTAL_RAWSCORE"
    echo "  DeNovoScore (MS:1002050): $TOTAL_DENOVO"
    echo "  SpecEValue (MS:1002052): $TOTAL_SPECEVALUE"
    echo "  EValue (MS:1002053): $TOTAL_EVALUE"
    echo "  QValue (MS:1002054): $TOTAL_QVALUE"
    echo "  PSMs at 1% FDR: $PSM_1PCT_FDR"

    cat > "$RESULTS_DIR/${LABEL}_metrics.txt" <<EOF
label=$LABEL
wall_time=${ELAPSED}
file_size=$FILE_SIZE
sii_count=$TOTAL_SII
raw_score_count=$TOTAL_RAWSCORE
denovo_score_count=$TOTAL_DENOVO
spec_evalue_count=$TOTAL_SPECEVALUE
evalue_count=$TOTAL_EVALUE
qvalue_count=$TOTAL_QVALUE
psm_1pct_fdr=$PSM_1PCT_FDR
EOF
}

run_dataset() {
    local DATASET_NAME="$1"
    local MZML="$2"
    local FASTA="$3"
    local MODS="$4"
    local EXTRA_ARGS="$5"
    local SKIP_BASELINE="${6:-false}"

    echo ""
    echo "########################################"
    echo "  DATASET: $DATASET_NAME"
    echo "########################################"

    if [ "$SKIP_BASELINE" != "true" ]; then
        run_msgfplus "$BASELINE_JAR" "${DATASET_NAME}_baseline" "$MZML" "$FASTA" "$MODS" "$EXTRA_ARGS"
    fi
    run_msgfplus "$NEW_JAR" "${DATASET_NAME}_new" "$MZML" "$FASTA" "$MODS" "$EXTRA_ARGS"
}

print_summary() {
    echo ""
    echo "========================================================"
    echo "  MULTI-DATASET BENCHMARK SUMMARY"
    echo "========================================================"
    echo ""
    printf "%-25s %10s %10s %12s %12s\n" "Dataset" "Time(s)" "PSMs@1%" "SII Count" "Scores OK"
    printf "%-25s %10s %10s %12s %12s\n" "-------------------------" "----------" "----------" "------------" "------------"

    for f in "$RESULTS_DIR"/*_metrics.txt; do
        [ -f "$f" ] || continue
        label=$(grep '^label=' "$f" | cut -d= -f2)
        wtime=$(grep '^wall_time=' "$f" | cut -d= -f2)
        psm=$(grep '^psm_1pct_fdr=' "$f" | cut -d= -f2)
        sii=$(grep '^sii_count=' "$f" | cut -d= -f2)
        raw=$(grep '^raw_score_count=' "$f" | cut -d= -f2)
        sii_half=$(( sii / 2 ))
        scores_ok="YES"
        if [ "$raw" != "$sii_half" ] 2>/dev/null; then
            scores_ok="MISMATCH"
        fi
        printf "%-25s %10s %10s %12s %12s\n" "$label" "${wtime}s" "$psm" "$sii" "$scores_ok"
    done

    echo ""
    echo "--- Per-dataset comparison ---"
    echo ""

    for dataset in tmt yeast lfq; do
        baseline_f="$RESULTS_DIR/${dataset}_baseline_metrics.txt"
        new_f="$RESULTS_DIR/${dataset}_new_metrics.txt"
        if [ -f "$baseline_f" ] && [ -f "$new_f" ]; then
            bt=$(grep '^wall_time=' "$baseline_f" | cut -d= -f2)
            nt=$(grep '^wall_time=' "$new_f" | cut -d= -f2)
            bp=$(grep '^psm_1pct_fdr=' "$baseline_f" | cut -d= -f2)
            np=$(grep '^psm_1pct_fdr=' "$new_f" | cut -d= -f2)
            speedup=$(python3 -c "print(f'{(1-$nt/$bt)*100:.1f}%')" 2>/dev/null || echo "N/A")
            psm_change=$(python3 -c "print(f'{($np/$bp-1)*100:+.1f}%')" 2>/dev/null || echo "N/A")
            echo "  $dataset: ${bt}s -> ${nt}s (${speedup} speedup), PSMs: $bp -> $np (${psm_change})"
        elif [ -f "$new_f" ]; then
            nt=$(grep '^wall_time=' "$new_f" | cut -d= -f2)
            np=$(grep '^psm_1pct_fdr=' "$new_f" | cut -d= -f2)
            echo "  $dataset: ${nt}s, PSMs@1%FDR: $np (new only, no baseline)"
        fi
    done
    echo ""
}

# ============================================================
# Main
# ============================================================

FILTER="${1:-all}"
SKIP_BASELINE=false
if [ "$FILTER" = "new-only" ]; then
    SKIP_BASELINE=true
    FILTER="all"
fi

echo "========== MS-GF+ MULTI-DATASET BENCHMARK =========="
echo "Date: $(date)"
echo "Java: $(java -version 2>&1 | head -1)"
echo "CPUs: $(nproc 2>/dev/null || sysctl -n hw.ncpu) (using $THREAD_COUNT threads)"
echo "Memory: $(free -h 2>/dev/null | awk '/Mem:/{print $2}' || sysctl -n hw.memsize 2>/dev/null | awk '{printf "%.0fG", $1/1073741824}')"
echo "Benchmark dir: $BENCHMARK_DIR"
echo "Filter: $FILTER"
echo "Baseline JAR: $BASELINE_JAR"
echo "New JAR: $NEW_JAR"

if [ "$FILTER" = "all" ] || [ "$FILTER" = "tmt" ]; then
    run_dataset "tmt" "$TMT_MZML" "$TMT_FASTA" "$TMT_MODS" "$TMT_ARGS" "$SKIP_BASELINE"
fi

if [ "$FILTER" = "all" ] || [ "$FILTER" = "yeast" ]; then
    run_dataset "yeast" "$YEAST_MZML" "$YEAST_FASTA" "$YEAST_MODS" "$YEAST_ARGS" "$SKIP_BASELINE"
fi

if [ "$FILTER" = "all" ] || [ "$FILTER" = "lfq" ]; then
    run_dataset "lfq" "$LFQ_MZML" "$LFQ_FASTA" "$LFQ_MODS" "$LFQ_ARGS" "$SKIP_BASELINE"
fi

print_summary
