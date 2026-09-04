#!/bin/bash
# Run Percolator on the .pin files produced by run_2arm_validation.sh.
# Useful when you finished searches with --skip-percolator and want to
# run Percolator separately (e.g. once Docker is available).
#
# For each dataset arm under results/validation/<dataset>/<arm>/output.pin:
#   - Runs biocontainers Percolator 3.7.1
#   - Records wall/RSS/CPU into the dataset's metrics.tsv (phase=percolator)
#   - Updates psm_summary.tsv with 1% FDR PSMs and unique peptides
#
# Usage:
#   bash run_2arm_percolator.sh <dataset>             # pxd001819 | tmt | astral | all
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib/measure.sh"

DATASET="${1:-}"
[ -n "$DATASET" ] || { echo "usage: $0 <pxd001819|tmt|astral|all>" >&2; exit 2; }

if [ "$DATASET" = "all" ]; then
    for ds in pxd001819 tmt astral; do
        bash "$0" "$ds" || true
    done
    exit 0
fi

if ! docker info >/dev/null 2>&1; then
    echo "ERROR: docker daemon is not reachable. Start Rancher Desktop / Docker first." >&2
    exit 1
fi

OUT_ROOT="$SCRIPT_DIR/results/validation/$DATASET"
METRICS="$OUT_ROOT/metrics.tsv"
SUMMARY="$OUT_ROOT/psm_summary.tsv"
[ -f "$METRICS" ] || { echo "ERROR: $METRICS not found — run searches first" >&2; exit 1; }

# Drop any prior percolator rows (so reruns don't pile up).
TMP_M="$(mktemp)"
awk -F'\t' 'NR==1 || $3!="percolator"' "$METRICS" > "$TMP_M" && mv "$TMP_M" "$METRICS"

# Rebuild summary in place.
TMP_S="$(mktemp)"
awk -F'\t' 'NR==1' "$SUMMARY" > "$TMP_S"
mv "$TMP_S" "$SUMMARY.new"

for ARM_DIR in "$OUT_ROOT"/*/; do
    ARM="$(basename "$ARM_DIR")"
    PIN="$ARM_DIR/output.pin"
    [ -f "$PIN" ] || { echo "skipping $ARM — no output.pin"; continue; }

    echo "=== [$DATASET / $ARM] Percolator ==="
    PERC_DIR="$ARM_DIR/percolator"
    mkdir -p "$PERC_DIR"
    PERC_TIME_LOG="$PERC_DIR/percolator_time.log"
    PERC_RC=0
    rm -f "$PERC_DIR"/target.psms.txt "$PERC_DIR"/decoy.psms.txt "$PERC_DIR"/weights.txt

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

    measure_record "$METRICS" "$DATASET" "$ARM" "percolator" "$PERC_TIME_LOG" "$PERC_RC" ""

    PSMS_1=0; PEP_1=0
    if [ "$PERC_RC" -eq 0 ] && [ -f "$PERC_DIR/target.psms.txt" ]; then
        PSMS_1=$(measure_perc_targets_at "$PERC_DIR/target.psms.txt" 0.01)
        PEP_1=$(measure_perc_unique_peptides_at "$PERC_DIR/target.psms.txt" 0.01)
        echo "[$ARM] Percolator OK. 1% FDR PSMs=$PSMS_1 peptides=$PEP_1"
    else
        echo "[$ARM] Percolator FAILED — $PERC_TIME_LOG"
        tail -10 "$PERC_TIME_LOG" >&2 || true
    fi

    # Look up MSGF+ PSM counts from the prior summary (or compute fresh).
    TGT=$(awk -F'\t' -v a="$ARM" '$1=="'"$DATASET"'" && $2==a {print $3}' "$SUMMARY" | tail -1)
    DCY=$(awk -F'\t' -v a="$ARM" '$1=="'"$DATASET"'" && $2==a {print $4}' "$SUMMARY" | tail -1)
    if [ -z "$TGT" ] || [ -z "$DCY" ]; then
        PIN_COUNTS=$(measure_count_pin_psms "$PIN")
        TGT="${PIN_COUNTS%%	*}"
        DCY="${PIN_COUNTS#*	}"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$DATASET" "$ARM" "$TGT" "$DCY" "$PSMS_1" "$PEP_1" >> "$SUMMARY.new"
done

mv "$SUMMARY.new" "$SUMMARY"
echo
echo "=== [$DATASET] DONE ==="
column -t -s $'\t' "$SUMMARY"
