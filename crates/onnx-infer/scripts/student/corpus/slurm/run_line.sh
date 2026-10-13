#!/usr/bin/env bash
# Run line $SLURM_ARRAY_TASK_ID of a commands file (one shell command per line), e.g. the
# qc_commands.sh written by msnet_inventory.py or a list of corpus_qc.py psm_table calls.
# Usage: sbatch --array=1-N%20 --cpus-per-task=8 --mem=32G --time=4:00:00 -o LOGDIR/%A_%a.log \
#            run_line.sh COMMANDS_FILE [WORKDIR]
set -euo pipefail
cmds=$1
[ $# -ge 2 ] && cd "$2"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$cmds")
[ -z "$line" ] && { echo "no line ${SLURM_ARRAY_TASK_ID} in $cmds"; exit 1; }
echo "[$(date -Is)] task ${SLURM_ARRAY_TASK_ID}: ${line:0:300}"
bash -c "$line"
echo "[$(date -Is)] done"
