#!/bin/bash
# run.sh <tag> <mzML> [extra andes args...]: andes --dia, Percolator 3.7.1, entrapment eval
set -e
tag=$1; mz=$2; shift 2
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
/home/singjc/Documents/github/andes/target/release/andes --spectrum $mz --database $B/pairlib-20261113/lib/pairlib.fasta \
  --mods mods.txt --output-pin $tag.pin --dia --charge 1..5 --min-length 7 --max-length 35 --threads 14 "$@" > $tag.log 2>&1
docker run --rm -u $(id -u):$(id -g) -v $PWD:/w -w /w quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 \
  percolator --seed 42 -Y --num-threads 14 -m $tag.t.psms -M $tag.d.psms -w $tag.weights --only-psms $tag.pin > $tag.perc.log 2>&1
echo "$tag: $(~/.venvs/redeem-external-dev-uv-20261022/bin/python /home/singjc/Documents/github/andes/docs/benchmarks/dia/eval_dia.py $tag.pin $tag.t.psms $tag.d.psms)" | tee -a RESULTS.txt
