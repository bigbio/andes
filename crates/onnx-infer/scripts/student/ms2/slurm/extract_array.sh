#!/usr/bin/env bash
# One SLURM array element of the andes-ms2 spectrum extraction.
#
#   sbatch --array=1-N%P --cpus-per-task=2 --mem=8G --time=12:00:00 [site flags] \
#       extract_array.sh PYTHON TASKS_TSV OUT_DIR SCRATCH_DIR "CONVERTER_COMMAND_PREFIX" [extra ms2_extract_raw.py args]
#
# PYTHON              interpreter with numpy + pyarrow + pandas
# TASKS_TSV           tasks.tsv from ms2_select.py (row N = array index N)
# OUT_DIR             shard root; the task writes OUT_DIR/<dataset>/part-<N>.parquet (+ .json)
# SCRATCH_DIR         node-local scratch for the converted mzML (deleted per run)
# CONVERTER_COMMAND_PREFIX how to run the raw-to-mzML converter, e.g.
#                     "singularity exec --bind DATA_ROOT IMAGE.sif CONVERTER"
# A task whose runs all succeeded is skipped on resubmission; a task with a failed run writes
# no part and exits 1, so re-sending the same array retries it.
set -euo pipefail
PY=$1; TASKS=$2; OUT=$3; SCRATCH=$4; CONVERTER=$5; shift 5
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${SLURM_ARRAY_TASK_ID:?run as a SLURM array job}"
exec "$PY" "$HERE/ms2_extract_raw.py" --tasks "$TASKS" --task-index "$SLURM_ARRAY_TASK_ID" \
    --out-dir "$OUT" --work-dir "$SCRATCH" --converter "$CONVERTER" "$@"
