#!/bin/bash
# run.sh <tag> <mzML> [extra andes args...]: andes --dia, Percolator 3.7.1 (--trainFDR 0.05, as
# `andes --dia --rescore` passes it; SEED= overrides 42), entrapment eval.
# B= the OpenDIA benchmark root (PAIRLIB FASTA, lib_targets.tsv); ANDES= the binary; PY= python with pandas.
set -e
tag=$1; mz=$2; shift 2
B=${B:?set B to the OpenDIA benchmark root}
HERE=$(cd "$(dirname "$0")" && pwd)
${ANDES:-andes} --spectrum $mz --database $B/pairlib-20261113/lib/pairlib.fasta \
  --fixed-mods "Carbamidomethyl (C)" --variable-mods "Oxidation (M)" --max-mods 1 \
  --output-pin $tag.pin --dia --charge 1..5 --min-length 7 --max-length 35 --threads 14 "$@" > $tag.log 2>&1
docker run --rm -u $(id -u):$(id -g) -v $PWD:/w -w /w quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 \
  percolator --seed ${SEED:-42} -Y --trainFDR 0.05 --num-threads 14 -m $tag.t.psms -M $tag.d.psms -w $tag.weights --only-psms $tag.pin > $tag.perc.log 2>&1
res=$(B=$B ${PY:-python3} $HERE/eval_dia.py $tag.pin $tag.t.psms $tag.d.psms)
echo "$tag: $res" | tee -a RESULTS.txt
