#!/bin/bash
# runb.sh <tag> <mzML> [extra andes args]: andes --dia (Phases A + B), Percolator on both PINs, entrapment eval.
set -e
tag=$1; mz=$2; shift 2
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
EV=/home/singjc/Documents/github/andes/docs/benchmarks/dia/eval_dia.py
PY=~/.venvs/redeem-external-dev-uv-20261022/bin/python
${ANDES_BIN:-/home/singjc/Documents/github/andes/target/release/andes} --spectrum $mz --database $B/pairlib-20261113/lib/pairlib.fasta \
  --mods mods.txt --output-pin $tag.pin --dia --charge ${CHARGE:-1..5} --min-length 7 --max-length 35 --threads 14 \
  --enzyme nocleavage --model hcd_astral_tryp --predict rt,ms2 "$@" > $tag.log 2>&1
perc() { docker run --rm -u $(id -u):$(id -g) -v $PWD:/w -w /w quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2 \
  percolator --seed 42 -Y --num-threads 14 -m $1.t.psms -M $1.d.psms -w $1.weights --only-psms $1.pin > $1.perc.log 2>&1; }
[ -f $tag.t.psms ] || perc $tag
perc $tag.dia
echo "$tag A: $($PY $EV $tag.pin $tag.t.psms $tag.d.psms)" | tee -a RESULTS.txt
echo "$tag B: $($PY $EV $tag.dia.pin $tag.dia.t.psms $tag.dia.d.psms)" | tee -a RESULTS.txt
