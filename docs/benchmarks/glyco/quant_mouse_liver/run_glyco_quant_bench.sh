#!/usr/bin/env bash
# Glyco LFQ + MBR benchmark: 5 mouse liver fractions in one invocation.
set -uo pipefail
export DOTNET_ROOT=/home/daichengxin/.dotnet

BIN=/home/daichengxin/gitrepo/andes/target/release/andes
DB=/home/daichengxin/andes-bench/databases/mouse_entrap.fasta
GDB=/home/daichengxin/gitrepo/pro/pGlyco-N-Mouse.gdb
RAWDIR=/Data/daicx/glyco-mouse
OUT=/tmp/glyco_quant_bench
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
  --lfq --glyco-mbr \
  --threads 112 \
  --output-pin "$PIN" \
  > "$LOG" 2>&1
rc=$?
echo "andes exit $rc"
touch "$OUT/.andes_done"
exit $rc
