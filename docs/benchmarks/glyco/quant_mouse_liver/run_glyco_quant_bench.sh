#!/usr/bin/env bash
# Glyco LFQ + MBR benchmark: 5 mouse liver fractions in one invocation.
set -uo pipefail
if [ "$#" -lt 5 ] || [ "$#" -gt 6 ]; then
  echo "usage: $0 <andes-binary> <mouse_entrap.fasta> <glycan.gdb> <dotnet-root> <raw-dir> [out-dir]" >&2
  exit 2
fi
BIN=$1
DB=$2
GDB=$3
export DOTNET_ROOT=$4

RAWDIR=$5
OUT=${6:-/tmp/glyco_quant_bench}
PIN=$OUT/mouseliver5.glyco.pin
LOG=$OUT/andes5.log

mkdir -p "$OUT"
echo "provenance: andes $($BIN --version 2>&1 | head -1) | $(uname -srm) | $(nproc) threads | $(date -u +%FT%TZ)"

"$BIN" \
  --spectrum "$RAWDIR/MouseLiver-Z-T-1.raw" \
  --spectrum "$RAWDIR/MouseLiver-Z-T-2.raw" \
  --spectrum "$RAWDIR/MouseLiver-Z-T-3.raw" \
  --spectrum "$RAWDIR/MouseLiver-Z-T-4.raw" \
  --spectrum "$RAWDIR/MouseLiver-Z-T-5.raw" \
  --database "$DB" \
  --decoy-strategy sequon \
  --glyco --glyco-glycan-gdb "$GDB" \
  --fragmentation HCD --fragment-tol-ppm 20 \
  --lfq --mbr \
  --threads 112 \
  --output-pin "$PIN" \
  > "$LOG" 2>&1
rc=$?
echo "andes exit $rc"
touch "$OUT/.andes_done"
exit $rc
