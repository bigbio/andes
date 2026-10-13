#!/usr/bin/env bash
# One SLURM array element of the andes-ms2 spectrum extraction.
#
#   sbatch --array=1-N%P --cpus-per-task=2 --mem=8G --time=12:00:00 [site flags] \
#       extract_array.sh PYTHON TASKS_TSV OUT_DIR SCRATCH_DIR "TRFP_COMMAND_PREFIX" [extra ms2_extract_raw.py args]
#
# PYTHON              interpreter with numpy + pyarrow + pandas
# TASKS_TSV           tasks.tsv from ms2_select.py (row N = array index N)
# OUT_DIR             shard root; the task writes OUT_DIR/<dataset>/part-<N>.parquet (+ .json)
# SCRATCH_DIR         node-local scratch for the converted mzML (deleted per run)
# TRFP_COMMAND_PREFIX how to run ThermoRawFileParser, e.g.
#                     "singularity exec --bind DATA_ROOT IMAGE.sif ThermoRawFileParser.sh"
# A finished task is skipped on resubmission, so the same array can be re-sent after failures.
set -euo pipefail
PY=$1; TASKS=$2; OUT=$3; SCRATCH=$4; TRFP=$5; shift 5
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${SLURM_ARRAY_TASK_ID:?run as a SLURM array job}"
exec "$PY" "$HERE/ms2_extract_raw.py" --tasks "$TASKS" --task-index "$SLURM_ARRAY_TASK_ID" \
    --out-dir "$OUT" --work-dir "$SCRATCH" --trfp "$TRFP" "$@"
