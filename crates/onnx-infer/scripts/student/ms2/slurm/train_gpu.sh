#!/usr/bin/env bash
# andes-ms2 training job (one GPU), then the parity reference for the Rust model.
#
#   sbatch --gres=gpu:1 --cpus-per-task=8 --mem=128G --time=24:00:00 [site flags] \
#       train_gpu.sh PYTHON SHARD_DIR RUN_DIR [extra ms2_train.py args]
#
# Writes RUN_DIR/{model.pt, andes_ms2.ams2, train.json, ms2_reference.tsv, train.log}.
set -euo pipefail
PY=$1; SHARDS=$2; RUN=$3; shift 3
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p "$RUN"
"$PY" -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())" | tee "$RUN/train.log"
"$PY" "$HERE/ms2_train.py" --shards "$SHARDS" --out "$RUN" --threads "${SLURM_CPUS_PER_TASK:-8}" "$@" 2>&1 | tee -a "$RUN/train.log"
"$PY" "$HERE/ms2_export.py" --model-dir "$RUN" --shards "$SHARDS" --n 2000 \
    --out-weights "$RUN/andes_ms2.ams2" --out-reference "$RUN/ms2_reference.tsv" 2>&1 | tee -a "$RUN/train.log"
