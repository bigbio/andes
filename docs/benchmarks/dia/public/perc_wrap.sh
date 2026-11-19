#!/bin/bash
# perc_wrap.sh: Percolator stand-in for `andes --percolator-bin`. Runs $BENCH_PERC (a command,
# word-split, default Percolator 3.7.1 from bioconda) with the arguments andes passes, under
# /usr/bin/time -v, writing one numbered time file per call into $BENCH_PERC_LOG.
# andes already runs inside the benchmark lock, so this does not lock again.
set -u
PERC=${BENCH_PERC:-${BENCH_ROOT:-/home/user/dia-bench}/tools/env-perc/bin/percolator}
dir=${BENCH_PERC_LOG:?BENCH_PERC_LOG not set}
mkdir -p "$dir"
n=$(ls "$dir" | grep -c "^call.*\.args$" || true)
n=$((n + 1))
printf '%s\n' "$PERC $*" > "$dir/call$n.args"
exec /usr/bin/time -v -o "$dir/call$n.time" $PERC "$@"
