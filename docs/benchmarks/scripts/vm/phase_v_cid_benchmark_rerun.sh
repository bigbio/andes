#!/usr/bin/env bash
# Re-run andes (rank/strong) on existing CID benchmark output.
# Keeps MSFragger PINs from a prior run; skips Comet (SKIP_COMET=1).
#
# Usage:
#   OUT=/srv/data/msgf-bench/cid-bench-20260606-1907 \
#   DATASETS=PXD016999,PXD007683 \
#   bash /srv/data/msgf-bench/phase_v_cid_benchmark_rerun.sh
set -euo pipefail

export OUT="${OUT:-/srv/data/msgf-bench/cid-bench-20260606-1907}"
export SKIP_MSFRAGGER=1
export SKIP_COMET=1
export DATASETS="${DATASETS:-PXD016999,PXD007683}"

echo "=== CID BENCHMARK RERUN $(date -Is) OUT=$OUT ==="
exec bash /srv/data/msgf-bench/phase_v_cid_benchmark.sh
