#!/bin/bash
# PB-Astral-full: andes --dia --rescore on one raw, Percolator 3.7.1 re-run per seed, entrapment eval.
#
#   run.sh search <tag> <spectra> [extra andes args...]   andes --dia --rescore -> <tag>.pin
#   run.sh perc   <tag> <seed>                            Percolator on <tag>.pin, then eval.py
#   run.sh eval   <tag> <seed>                            eval.py on an existing Percolator output
#
# Environment:
#   FASTA       the paired entrapment FASTA (ProteoBench entrapment FASTA, or PAIRLIB's pairlib.fasta)
#   ANDES       the andes binary (default: andes)
#   PERCOLATOR  a Percolator 3.7.1 executable, or a wrapper script that runs one in a container
#   PY          python with pandas and pyarrow (default: python3)
#   THREADS     threads for andes and Percolator (default 14)
#   EVAL_ARGS   extra eval.py arguments (e.g. "--library lib_targets.tsv" for PAIRLIB)
# Results: one line per seed appended to RESULTS.txt; run time and peak memory in <tag>.time.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
: "${FASTA:?set FASTA}" "${PERCOLATOR:?set PERCOLATOR}"
ANDES=${ANDES:-andes}; PY=${PY:-python3}; THREADS=${THREADS:-14}
search=(--enzyme nocleavage --fixed-mods "Carbamidomethyl (C)" --variable-mods "" --charge 1..4
        --min-length 7 --max-length 30)
stage=$1; tag=$2; shift 2
case $stage in
search)
  mz=$1; shift
  /usr/bin/time -v -o "$tag.time" "$ANDES" --spectrum "$mz" --database "$FASTA" "${search[@]}" \
    --dia --model hcd_astral_tryp --predict rt,ms2 --rescore --percolator-bin "$PERCOLATOR" \
    --threads "$THREADS" --output-pin "$tag.pin" "$@" > "$tag.log" 2>&1
  grep -E "Elapsed|Maximum resident" "$tag.time" ;;
perc)
  seed=$1
  "$PERCOLATOR" --seed "$seed" -Y --trainFDR 0.05 --num-threads "$THREADS" \
    -m "$tag.s$seed.t.psms" -M "$tag.s$seed.d.psms" -w "$tag.s$seed.weights" --only-psms "$tag.pin" \
    > "$tag.s$seed.perc.log" 2>&1
  exec bash "$0" eval "$tag" "$seed" ;;
eval)
  seed=$1
  # shellcheck disable=SC2086
  res=$("$PY" "$HERE/eval.py" andes "$tag.pin" "$tag.s$seed.t.psms" "$tag.s$seed.d.psms" ${EVAL_ARGS:-})
  echo "$tag seed $seed: $res" | tee -a RESULTS.txt ;;
*) echo "unknown stage $stage" >&2; exit 2 ;;
esac
