#!/bin/bash
# Run Percolator (via biocontainers docker image) on a .pin file.
# Usage: bash run_percolator_docker.sh <pin-file> <output-dir> <label>
#
# Output:
#   <output-dir>/<label>.target.psms.txt
#   <output-dir>/<label>.decoy.psms.txt
#   <output-dir>/<label>.weights.txt
#   <output-dir>/<label>.log
set -euo pipefail

PIN_FILE="$(realpath "$1")"
OUTPUT_DIR="$(realpath "$2")"
LABEL="$3"

mkdir -p "$OUTPUT_DIR"
PIN_DIR="$(dirname "$PIN_FILE")"
PIN_NAME="$(basename "$PIN_FILE")"

# Mount the pin's directory as /pin and the output dir as /out so percolator
# can read/write in-place without requiring a copy.
docker run --rm --platform linux/amd64 \
    -v "$PIN_DIR":/pin:ro \
    -v "$OUTPUT_DIR":/out \
    quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 \
    percolator \
        --seed 42 \
        --results-psms "/out/${LABEL}.target.psms.txt" \
        --decoy-results-psms "/out/${LABEL}.decoy.psms.txt" \
        --weights "/out/${LABEL}.weights.txt" \
        --only-psms \
        "/pin/${PIN_NAME}" \
    > "$OUTPUT_DIR/${LABEL}.log" 2>&1

if [ ! -f "$OUTPUT_DIR/${LABEL}.target.psms.txt" ]; then
    echo "[$LABEL] PERCOLATOR FAILED — last log lines:"
    tail -20 "$OUTPUT_DIR/${LABEL}.log"
    exit 1
fi

# Summary — locate the q-value column by header name, because Sage output
# has an extra FileName column (q-value is col 4) while MS-GF+ has it at col 3.
TGT_FILE="$OUTPUT_DIR/${LABEL}.target.psms.txt"
DCY_FILE="$OUTPUT_DIR/${LABEL}.decoy.psms.txt"
T_TOTAL=$(($(wc -l < "$TGT_FILE") - 1))
D_TOTAL=$(($(wc -l < "$DCY_FILE") - 1))
Q_COL=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="q-value") {print i; exit}}' "$TGT_FILE")
if [ -z "$Q_COL" ]; then
    echo "[$LABEL] WARNING: could not find q-value column; header was:" >&2
    head -1 "$TGT_FILE" >&2
    Q_COL=3
fi
T_1PCT=$(awk -F'\t' -v q="$Q_COL" 'NR>1 && $q<=0.01 {c++} END {print c+0}' "$TGT_FILE")
T_5PCT=$(awk -F'\t' -v q="$Q_COL" 'NR>1 && $q<=0.05 {c++} END {print c+0}' "$TGT_FILE")

echo "[$LABEL] targets_total=$T_TOTAL targets_1pct=$T_1PCT targets_5pct=$T_5PCT decoys_total=$D_TOTAL q_col=$Q_COL"
