#!/bin/bash
# 2-arm validation: master vs branch.
#
# For each dataset, runs three phases per arm:
#   1) BuildSA  — DB index construction (cold; FASTA copied to a per-arm work dir)
#   2) Search   — MSGF+ search on the prebuilt index, emits .pin
#   3) Percolator — biocontainers Docker image; emits target/decoy PSM tables
#
# All phases are wrapped in /usr/bin/time -l (macOS) or -v (Linux) to capture
# wall, max RSS, user/sys CPU. Per-arm metrics written to a TSV.
#
# Usage:
#   bash run_2arm_validation.sh <dataset>             # one of: pxd001819 | tmt | astral | all
#   bash run_2arm_validation.sh <dataset> --skip-percolator
#   MSGFPLUS_THREADS=8 MSGFPLUS_MEMORY=16384m bash run_2arm_validation.sh astral
#
# Layout produced under benchmark/results/validation/<dataset>/:
#   metrics.tsv         all phases x both arms
#   psm_summary.tsv     pin PSM counts + percolator 1% FDR PSMs / peptides
#   <arm>/build.log
#   <arm>/search.log
#   <arm>/output.pin
#   <arm>/percolator/{target,decoy}.psms.txt
#
# Both arms use a private working FASTA copy so BuildSA caches don't collide.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
source "$SCRIPT_DIR/lib/measure.sh"

DATASET="${1:-}"
SKIP_PERC=0
[ "${2:-}" = "--skip-percolator" ] && SKIP_PERC=1

if [ -z "$DATASET" ]; then
    echo "usage: $0 <pxd001819|tmt|astral|all> [--skip-percolator]" >&2
    exit 2
fi

if [ "$DATASET" = "all" ]; then
    for ds in pxd001819 tmt astral; do
        bash "$0" "$ds" "${2:-}"
    done
    bash "$SCRIPT_DIR/render_validation_report.sh"
    exit 0
fi

THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
MATCHES_PER_SPEC="${MSGFPLUS_MATCHES_PER_SPEC:-1}"

case "$DATASET" in
    pxd001819)
        DATA_DIR="$SCRIPT_DIR/data/PXD001819"
        MZML="$DATA_DIR/UPS1_5000amol_R1.mzML"
        FASTA_SRC="$DATA_DIR/PXD001819_uniprot_yeast_ups.fasta"
        MODS="$DATA_DIR/mods.txt"
        SEARCH_ARGS="-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"
        JVM_MEMORY="${MSGFPLUS_MEMORY:-4096m}"
        ;;
    tmt)
        DATA_DIR="$SCRIPT_DIR/data/PXD007683"
        MZML="$DATA_DIR/a05058.mzML"
        FASTA_SRC="$DATA_DIR/PXD007683_UP000005640_UP000002311_reviewed.fasta"
        MODS="$DATA_DIR/mods.txt"
        SEARCH_ARGS="-tda 1 -t 20ppm -ti -1,2 -m 1 -inst 1 -e 1 -protocol 4 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"
        JVM_MEMORY="${MSGFPLUS_MEMORY:-8192m}"
        ;;
    astral)
        DATA_DIR="$SCRIPT_DIR/data/Astral_ProteoBench_Module_8"
        MZML="$DATA_DIR/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML"
        FASTA_SRC="$DATA_DIR/ProteoBenchFASTA_MixedSpecies_HYE.fasta"
        MODS="$DATA_DIR/mods.txt"
        SEARCH_ARGS="-tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n $MATCHES_PER_SPEC -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"
        JVM_MEMORY="${MSGFPLUS_MEMORY:-8192m}"
        ;;
    *)
        echo "ERROR: unknown dataset '$DATASET'" >&2; exit 2 ;;
esac

# Baseline is origin/dev (PR #24 fork point) by default. The fork's master predates
# DirectPinWriter (PRs #20/22) so it can't emit .pin and is not a viable apples-to-apples
# baseline. Override BASELINE_JAR if you want a different reference build.
BASELINE_JAR="${BASELINE_JAR:-$SCRIPT_DIR/baseline-dev/MSGFPlus.jar}"
BRANCH_JAR="${BRANCH_JAR:-$SCRIPT_DIR/branch-mzid/MSGFPlus.jar}"
BASELINE_LABEL="${BASELINE_LABEL:-baseline-dev}"
[ -f "$BASELINE_JAR" ] || { echo "missing $BASELINE_JAR"; exit 1; }
[ -f "$BRANCH_JAR" ]   || { echo "missing $BRANCH_JAR (build improve-mzid-suffix and copy to $BRANCH_JAR)"; exit 1; }
[ -f "$MZML" ]       || { echo "missing $MZML"; exit 1; }
[ -f "$FASTA_SRC" ]  || { echo "missing $FASTA_SRC"; exit 1; }
[ -f "$MODS" ]       || { echo "missing $MODS"; exit 1; }

OUT_ROOT="$SCRIPT_DIR/results/validation/$DATASET"
METRICS="$OUT_ROOT/metrics.tsv"
SUMMARY="$OUT_ROOT/psm_summary.tsv"

mkdir -p "$OUT_ROOT"
measure_emit_header "$METRICS"
printf "dataset\tarm\tmsgf_targets\tmsgf_decoys\tperc_psms_1pct\tperc_peptides_1pct\n" > "$SUMMARY"

run_arm() {
    local ARM="$1"        # master | branch
    local JAR="$2"
    local ARM_DIR="$OUT_ROOT/$ARM"
    mkdir -p "$ARM_DIR"

    # Each arm gets a private FASTA so BuildSA caches don't collide.
    local FASTA="$ARM_DIR/$(basename "$FASTA_SRC")"
    rm -f "$ARM_DIR"/*.canno "$ARM_DIR"/*.cnlcp "$ARM_DIR"/*.csarr "$ARM_DIR"/*.cseq
    rm -f "$ARM_DIR"/*.revCat.fasta
    cp -f "$FASTA_SRC" "$FASTA"

    echo
    echo "=== [$DATASET / $ARM] BuildSA (DB construction) ==="
    local BUILD_LOG="$ARM_DIR/build.log"
    local BUILD_RC=0
    measure_run "$BUILD_LOG" BUILD_RC \
        java "-Xmx${JVM_MEMORY}" -cp "$JAR" edu.ucsd.msjava.msdbsearch.BuildSA \
        -d "$FASTA" -tda 1
    measure_record "$METRICS" "$DATASET" "$ARM" "build" "$BUILD_LOG" "$BUILD_RC" ""
    if [ "$BUILD_RC" -ne 0 ]; then
        echo "[$ARM] BuildSA FAILED — last lines of $BUILD_LOG:"
        tail -20 "$BUILD_LOG" >&2
        return 1
    fi
    echo "[$ARM] BuildSA done. wall=$(awk -F'\t' -v a="$ARM" '$2==a && $3=="build"{print $4}' "$METRICS")s"

    echo
    echo "=== [$DATASET / $ARM] Search ==="
    local SEARCH_LOG="$ARM_DIR/search.log"
    local PIN="$ARM_DIR/output.pin"
    local SEARCH_RC=0
    rm -f "$PIN"
    measure_run "$SEARCH_LOG" SEARCH_RC \
        java "-Xmx${JVM_MEMORY}" -jar "$JAR" \
        -s "$MZML" -d "$FASTA" -mod "$MODS" -o "$PIN" \
        $SEARCH_ARGS
    if [ ! -f "$PIN" ]; then
        # Master arm: pin is auto-emitted next to the .mzid output by DirectPinWriter.
        # Find any .pin under ARM_DIR that the search produced.
        local found
        found=$(find "$ARM_DIR" -maxdepth 2 -name '*.pin' -print -quit)
        if [ -n "$found" ]; then
            cp -f "$found" "$PIN"
        fi
    fi
    local PIN_COUNTS
    PIN_COUNTS=$(measure_count_pin_psms "$PIN")
    local TGT="${PIN_COUNTS%%	*}"
    local DCY="${PIN_COUNTS#*	}"
    measure_record "$METRICS" "$DATASET" "$ARM" "search" "$SEARCH_LOG" "$SEARCH_RC" "targets=$TGT decoys=$DCY"
    if [ "$SEARCH_RC" -ne 0 ] || [ ! -f "$PIN" ]; then
        echo "[$ARM] Search FAILED — last lines:"
        tail -20 "$SEARCH_LOG" >&2
        return 1
    fi
    echo "[$ARM] Search done. wall=$(awk -F'\t' -v a="$ARM" '$2==a && $3=="search"{print $4}' "$METRICS")s, MSGF+ PSMs T=$TGT D=$DCY"

    local PERC_PSMS_1=0 PERC_PEP_1=0
    if [ "$SKIP_PERC" -eq 0 ]; then
        echo
        echo "=== [$DATASET / $ARM] Percolator ==="
        local PERC_DIR="$ARM_DIR/percolator"
        mkdir -p "$PERC_DIR"
        local PERC_LOG="$PERC_DIR/percolator.log"
        local PERC_TIME_LOG="$PERC_DIR/percolator_time.log"
        local PERC_RC=0
        local TGT_FILE="$PERC_DIR/target.psms.txt"
        local DCY_FILE="$PERC_DIR/decoy.psms.txt"
        local WTS_FILE="$PERC_DIR/weights.txt"
        rm -f "$TGT_FILE" "$DCY_FILE" "$WTS_FILE"

        measure_run "$PERC_TIME_LOG" PERC_RC \
            docker run --rm --platform linux/amd64 \
                -v "$ARM_DIR":/pin:ro \
                -v "$PERC_DIR":/out \
                quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 \
                percolator \
                    --seed 42 \
                    --results-psms /out/target.psms.txt \
                    --decoy-results-psms /out/decoy.psms.txt \
                    --weights /out/weights.txt \
                    --only-psms \
                    /pin/output.pin
        cp -f "$PERC_TIME_LOG" "$PERC_LOG" 2>/dev/null || true

        measure_record "$METRICS" "$DATASET" "$ARM" "percolator" "$PERC_TIME_LOG" "$PERC_RC" ""

        if [ "$PERC_RC" -eq 0 ] && [ -f "$TGT_FILE" ]; then
            PERC_PSMS_1=$(measure_perc_targets_at "$TGT_FILE" 0.01)
            PERC_PEP_1=$(measure_perc_unique_peptides_at "$TGT_FILE" 0.01)
            echo "[$ARM] Percolator OK. 1% FDR PSMs=$PERC_PSMS_1 peptides=$PERC_PEP_1"
        else
            echo "[$ARM] Percolator FAILED — last lines:"
            tail -20 "$PERC_TIME_LOG" >&2
        fi
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$DATASET" "$ARM" "$TGT" "$DCY" "$PERC_PSMS_1" "$PERC_PEP_1" >> "$SUMMARY"
}

run_arm "$BASELINE_LABEL" "$BASELINE_JAR"
run_arm branch "$BRANCH_JAR"

echo
echo "=== [$DATASET] DONE ==="
echo "Metrics: $METRICS"
echo "Summary: $SUMMARY"
echo
column -t -s $'\t' "$SUMMARY"
echo
column -t -s $'\t' "$METRICS"
