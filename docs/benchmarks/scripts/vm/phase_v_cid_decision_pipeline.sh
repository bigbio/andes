#!/usr/bin/env bash
# Full Phase V CID decision pipeline:
#   1. Train CID models (MSNet → intensity + scoring store)
#   2. Benchmark MSFragger/Comet vs andes rank/strong (3×4 raw files)
#
# Usage:
#   nohup bash /srv/data/msgf-bench/phase_v_cid_decision_pipeline.sh \
#     > /srv/data/msgf-bench/phase_v_cid_decision.log 2>&1 &
set -euo pipefail

LOG="${LOG:-/srv/data/msgf-bench/phase_v_cid_decision.log}"
echo "=== PHASE_V_CID_DECISION_START $(date -Is) ===" | tee -a "$LOG"

bash /srv/data/msgf-bench/phase_v_cid_train_models.sh 2>&1 | tee -a "$LOG"

export INTENSITY_MODEL=/srv/data/msnet/intensity_model_cid.parquet
export SCORING_STORE=/srv/data/msnet/models_cid_msnet.parquet
bash /srv/data/msgf-bench/phase_v_cid_benchmark.sh 2>&1 | tee -a "$LOG"

echo "=== PHASE_V_CID_DECISION_DONE $(date -Is) ===" | tee -a "$LOG"
