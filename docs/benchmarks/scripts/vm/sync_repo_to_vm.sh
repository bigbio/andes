#!/usr/bin/env bash
# Rsync local msgf-rust tree to the bench VM repo path.
#
# Usage (from msgf-rust/):
#   bash benchmark/vm/sync_repo_to_vm.sh
#
# Override target:
#   VM_HOST=user@your-bench-host VM_REPO=/path/to/msgf-rust bash docs/benchmarks/scripts/vm/sync_repo_to_vm.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
VM_HOST="${VM_HOST:?set VM_HOST=user@host — no default, this is not our infrastructure to publish}"
VM_REPO="${VM_REPO:-/srv/data/msgf-bench/repo/msgf-rust}"
SSH_CTRL="${SSH_CTRL:-/tmp/msgfplus-bench.sock}"
RSYNC_SSH="ssh -S ${SSH_CTRL}"

echo "Syncing ${REPO_ROOT} -> ${VM_HOST}:${VM_REPO}"
rsync -az --delete \
  --exclude target/ \
  --exclude .git/ \
  --exclude 'benchmark/results/' \
  --exclude 'test-fixtures/' \
  -e "$RSYNC_SSH" \
  "${REPO_ROOT}/" "${VM_HOST}:${VM_REPO}/"

# Deploy VM gate scripts alongside existing bench harness.
rsync -az \
  -e "$RSYNC_SSH" \
  "${SCRIPT_DIR}/phase_v_train_intensity.sh" \
  "${SCRIPT_DIR}/phase_v_strong_score_gate.sh" \
  "${SCRIPT_DIR}/phase_v_cid_train_models.sh" \
  "${SCRIPT_DIR}/phase_v_cid_benchmark.sh" \
  "${SCRIPT_DIR}/phase_v_cid_benchmark_rerun.sh" \
  "${SCRIPT_DIR}/comet_cid_run.sh" \
  "${SCRIPT_DIR}/phase_v_cid_decision_pipeline.sh" \
  "${SCRIPT_DIR}/summarize_phase_v_gate.py" \
  "${SCRIPT_DIR}/summarize_phase_v_cid_benchmark.py" \
  "${VM_HOST}:/srv/data/msgf-bench/"

# Deploy intensity aggregation tooling.
rsync -az \
  -e "$RSYNC_SSH" \
  "${REPO_ROOT}/../internal-docs/msnet-model-catalog/msnet_intensity_agg.py" \
  "${REPO_ROOT}/../internal-docs/msnet-model-catalog/msnet_to_flat.py" \
  "${VM_HOST}:/srv/data/msnet/"

echo "Sync done. Build on VM:"
echo "  ssh -S ${SSH_CTRL} ${VM_HOST} 'cd ${VM_REPO} && cargo +1.95.0 build --release --features thermo -p andes'"
