#!/bin/bash
# Public DIA benchmark harness (andes --dia on fixed RT slices of ProteoBench Astral DIA raws 01/02,
# paired-entrapment library from UniProt human, Percolator 3.7.1, 5 seeds).
#
#   bench.sh [options] <tag> <andes-binary> [quick|full] [extra andes args...]
#   bench.sh [options] --pin-from <tag-or-dir> <tag> - [quick|full]
#
# Options (before <tag>):
#   --perc "CMD"        Percolator command for every Percolator run (in-andes passes and the
#                       seed re-runs); default Percolator 3.7.1 (bioconda). The native engine of
#                       main is $BENCH_ROOT/tools/percnative-main (percnative/).
#   --seeds "S ..."     override the tier's seeds.
#   --raws "01 02"      override the tier's raws.
#   --pin-from X        skip the search: re-score the final PINs of runs/X/<raw>/<raw>.pin (or
#                       X/<raw>/<raw>.pin when X is a directory). The binary argument is ignored.
#   --no-results        do not append to runs/RESULTS.tsv.
#
# Tiers: quick = raw 01, seeds 42 1 2; full = raws 01 02, seeds 42 1 2 3 7.
# Per raw: one `andes --dia --rescore` (two Percolator passes with the per-run RT refinement
# between them, seed 42, Percolator through perc_wrap.sh), then Percolator re-run on the final
# (RT-refined) PIN for each seed with the arguments andes uses (-Y --trainFDR 0.05
# --num-threads 4 --only-psms), each evaluated with eval_dia.py. Seed 42 reuses andes's own
# second pass (its last non-static Percolator call), which is that same command on that same
# PIN (checked byte-identical). Every search and every
# Percolator run holds $BENCH_ROOT/bench.lock (flock), so measurements never overlap.
# BENCH_ROOT (default /home/user/dia-bench) holds data/, tools/, bin/ and runs/.
set -euo pipefail
ROOT=${BENCH_ROOT:-/home/user/dia-bench}
S=$(cd "$(dirname "$0")" && pwd)
LOCK=$ROOT/bench.lock
EVAL=$S/eval_dia.py; [ -f "$EVAL" ] || EVAL=$S/../eval_dia.py
source "${BENCH_FROZEN:-$S/frozen.env}"
PERC_DEFAULT=$ROOT/tools/env-perc/bin/percolator
PERC=$PERC_DEFAULT; SEEDS=; RAWS=; PIN_FROM=; RESULTS=$ROOT/runs/RESULTS.tsv
while [ $# -gt 0 ]; do
  case $1 in
    --perc) PERC=$2; shift 2;;
    --seeds) SEEDS=$2; shift 2;;
    --raws) RAWS=$2; shift 2;;
    --pin-from) PIN_FROM=$2; shift 2;;
    --no-results) RESULTS=; shift;;
    -h|--help) sed -n 2,25p "$0"; exit 0;;
    --) shift; break;;
    -*) echo "unknown option $1" >&2; exit 2;;
    *) break;;
  esac
done
[ $# -ge 2 ] || { sed -n 2,25p "$0"; exit 2; }
tag=$1; bin=$2; tier=${3:-full}; shift 2; [ $# -gt 0 ] && shift
extra=("$@")
case $tier in
  quick) : "${SEEDS:=42 1 2}"; : "${RAWS:=01}";;
  full) : "${SEEDS:=42 1 2 3 7}"; : "${RAWS:=01 02}";;
  *) echo "tier must be quick or full" >&2; exit 2;;
esac
THREADS=4
out=$ROOT/runs/$tag
mkdir -p "$out"
if [ -n "$PIN_FROM" ]; then
  src=$PIN_FROM; [ -d "$src" ] || src=$ROOT/runs/$PIN_FROM
  [ -d "$src" ] || { echo "no run $PIN_FROM" >&2; exit 2; }
  commit=$(cat "$src/commit" 2>/dev/null || echo unknown)
  bin=$(cat "$src/binary" 2>/dev/null || echo -)
else
  bin=$(readlink -f "$bin")
  [ -x "$bin" ] || { echo "not executable: $bin" >&2; exit 2; }
  commit=$(cat "$bin.commit" 2>/dev/null || echo unknown)
  commit=${commit:0:7}
fi
echo "$commit" > "$out/commit"; echo "$bin" > "$out/binary"
{
  echo "date $(date -Is)"; echo "tag $tag"; echo "tier $tier"; echo "binary $bin"
  echo "commit $commit"; [ -f "$bin" ] && echo "binary_sha256 $(sha256sum "$bin" | cut -d' ' -f1)"
  echo "pin_from ${PIN_FROM:--}"; echo "perc $PERC"; echo "seeds $SEEDS"; echo "raws $RAWS"
  echo "extra ${extra[*]:-}"; echo "host $(nproc) cores, $(free -g | awk '/Mem/{print $2}') GB"
  echo "cpu $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
} > "$out/meta.txt"

for raw in $RAWS; do
  d=$out/$raw; mkdir -p "$d"
  slice_var=SLICE_$raw; slice=${!slice_var}
  pin=$d/$raw.pin
  if [ -n "$PIN_FROM" ]; then
    ln -sf "$src/$raw/$raw.pin" "$pin"
    cp -f "$src/$raw/search.time" "$d/search.time" 2>/dev/null || true
    rm -rf "$d/perc_calls"; cp -r "$src/$raw/perc_calls" "$d/perc_calls" 2>/dev/null || true
  else
    rm -rf "$d/perc_calls"
    echo "[$tag/$raw] search $(date +%T)" >&2
    ( cd "$d" && BENCH_PERC="$PERC" BENCH_PERC_LOG="$d/perc_calls" \
      flock "$LOCK" /usr/bin/time -v -o "$d/search.time" "$bin" --spectrum "$slice" --database "$FASTA" \
        --fixed-mods "Carbamidomethyl (C)" --variable-mods "Oxidation (M)" --max-mods 1 \
        --dia --enzyme nocleavage --model hcd_astral_tryp --charge 1..5 --min-length 7 --max-length 35 \
        --predict rt,ms2 --threads $THREADS --rescore --percolator-bin "$S/perc_wrap.sh" \
        --output-pin "$pin" "${extra[@]}" > "$d/search.log" 2>&1 ) \
      || { echo "[$tag/$raw] search failed, see $d/search.log" >&2; tail -5 "$d/search.log" >&2; exit 1; }
    rm -f "$d"/andes-candidx-*.bin
  fi
  for s in $SEEDS; do
    sd=$d/seed$s; mkdir -p "$sd"
    last=$(grep -L -- "--static" $(ls "$d"/perc_calls/call*.args 2>/dev/null | sort -V) /dev/null 2>/dev/null \
      | grep -v '^/dev/null$' | tail -1 || true)
    if [ -z "$PIN_FROM" ] && [ "$s" = 42 ] && [ -n "$last" ] && grep -q -- "--seed 42 " "$last" \
       && grep -q -- "--results-psms $d/$raw.percolator.target.psms.txt " "$last" \
       && [ -f "$d/$raw.percolator.target.psms.txt" ]; then
      echo "[$tag/$raw] seed 42: reusing andes's second Percolator pass (same command)" >&2
      cp "$d/$raw.percolator.target.psms.txt" "$sd/t.psms"; cp "$d/$raw.percolator.decoy.psms.txt" "$sd/d.psms"
      cp "$d/$raw.percolator.weights.txt" "$sd/weights"; echo "$last" > "$sd/perc.reused"; rm -f "$sd/perc.time"
    else
    echo "[$tag/$raw] percolator seed $s $(date +%T)" >&2
    flock "$LOCK" /usr/bin/time -v -o "$sd/perc.time" $PERC --seed $s -Y --trainFDR 0.05 \
      --num-threads $THREADS --only-psms --results-psms "$sd/t.psms" --decoy-results-psms "$sd/d.psms" \
      --weights "$sd/weights" "$pin" > "$sd/perc.log" 2>&1 \
      || { echo "[$tag/$raw] percolator seed $s failed, see $sd/perc.log" >&2; exit 1; }
    fi
    python3 -I "$EVAL" "$pin" "$sd/t.psms" "$sd/d.psms" "$LIB" > "$sd/eval.txt"
    echo "[$tag/$raw] seed $s: $(cat "$sd/eval.txt")" >&2
  done
  python3 -I "$S/summarize.py" "$d" "$tag" "$commit" "$raw" "$tier" "${extra[*]:-}" "$PERC" $RESULTS \
    | tee "$d/summary.tsv"
done
