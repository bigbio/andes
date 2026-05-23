#!/usr/bin/env bash
# Bisect oracle for the score_psm under-scoring regression.
#
# - Builds msgf-rust at the current commit
# - Runs it on PXD001819 single-threaded with --max-spectra 30000
# - Greps scan=28787's RawScore from the pin (column 7)
# - Appends <sha>,<rawscore> to /tmp/bisect-trace.csv (cumulative log)
# - Exits 0 (good) if RawScore >= 290
# - Exits 1 (bad)  if RawScore <  200
# - Exits 125 (skip) on build failure or missing scan in pin
#
# Determinism: --threads 1 eliminates rayon nondeterminism. The same
# commit produces the same RawScore across runs.

set -uo pipefail

REPO_ROOT="/Users/yperez/work/msgfplus-workspace/astral-speed-score-fix"
PXD_MZML="/Users/yperez/work/msgfplus-workspace/benchmark/data/PXD001819/UPS1_5000amol_R1.mzML"
PXD_FASTA="/Users/yperez/work/msgfplus-workspace/benchmark/data/PXD001819/PXD001819_uniprot_yeast_ups.fasta"
TRACE_CSV="/tmp/bisect-trace.csv"
PIN_OUT="/tmp/bisect.pin"

cd "$REPO_ROOT/rust"
SHA=$(git rev-parse --short HEAD)

# Skip non-existent inputs (would lead to false bad).
if [ ! -f "$PXD_MZML" ] || [ ! -f "$PXD_FASTA" ]; then
    echo "[$SHA] missing PXD001819 fixture — skip"
    exit 125
fi

# Build. Use full build (not --quiet) so cargo errors are visible in
# `git bisect run` logs.
if ! cargo build --release --bin msgf-rust 2>&1 | tail -5; then
    echo "[$SHA] build failed — skip"
    echo "$SHA,BUILD_FAIL" >> "$TRACE_CSV"
    exit 125
fi

BIN="$REPO_ROOT/rust/target/release/msgf-rust"
rm -f "$PIN_OUT"

if ! "$BIN" \
        --spectrum "$PXD_MZML" \
        --database "$PXD_FASTA" \
        --output-pin "$PIN_OUT" \
        --precursor-tol-ppm 5 \
        --isotope-error-min=0 \
        --isotope-error-max=1 \
        --top-n 1 \
        --threads 1 \
        --max-spectra 30000 \
        > /tmp/bisect.log 2>&1; then
    echo "[$SHA] msgf-rust run failed — skip"
    echo "$SHA,RUN_FAIL" >> "$TRACE_CSV"
    exit 125
fi

# Column 7 of the pin is RawScore.
RAW=$(awk -F'\t' 'NR>1 && $3 == 28787 {print $7; exit}' "$PIN_OUT")

if [ -z "$RAW" ]; then
    echo "[$SHA] scan=28787 not in pin output — skip"
    echo "$SHA,MISSING_SCAN" >> "$TRACE_CSV"
    exit 125
fi

echo "$SHA,$RAW" >> "$TRACE_CSV"
echo "[$SHA] scan=28787 RawScore=$RAW"

if [ "$RAW" -ge 290 ] 2>/dev/null; then
    exit 0  # good
fi
if [ "$RAW" -lt 200 ] 2>/dev/null; then
    exit 1  # bad
fi

# In the dead-band 200..290: skip to avoid mis-bisecting on intermediate.
echo "[$SHA] RawScore=$RAW in dead band 200..290 — skip"
exit 125
