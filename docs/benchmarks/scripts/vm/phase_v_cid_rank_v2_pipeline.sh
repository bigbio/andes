#!/usr/bin/env bash
# Retrain rank (MSNet + Lumos self-label) then re-benchmark vs MSFragger.
#
# Usage:
#   nohup bash /srv/data/msgf-bench/phase_v_cid_rank_v2_pipeline.sh \
#     > /srv/data/msgf-bench/phase_v_cid_rank_v2_pipeline.log 2>&1 &
set -euo pipefail

LOG="${LOG:-/srv/data/msgf-bench/phase_v_cid_rank_v2_pipeline.log}"
exec > >(tee -a "$LOG") 2>&1

echo "=== RANK_V2_PIPELINE_START $(date -Is) ==="

export REFLAT=1
export REFLAT_IDS="${REFLAT_IDS:-PXD016986,PXD000865-CID}"
export PREPARE_LUMOS=1
export SAMPLE_865="${SAMPLE_865:-800000}"

bash /srv/data/msgf-bench/phase_v_cid_train_models.sh

python3 /srv/data/msgf-bench/diagnose_rank_dilution.py || true

export OUT="${OUT:-/srv/data/msgf-bench/cid-bench-rank-v2}"
export SKIP_MSFRAGGER=1
export SKIP_COMET=1
export DATASETS="${DATASETS:-PXD007683,PXD016999}"
bash /srv/data/msgf-bench/phase_v_cid_benchmark_rerun.sh

python3 /srv/data/msgf-bench/summarize_phase_v_cid_benchmark.py "$OUT" || true

echo "=== RANK_V2_PIPELINE_DONE $(date -Is) OUT=$OUT ==="
