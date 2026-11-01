#!/bin/sh
# Absence check only (no amounts.tsv): rows in the blank run. Replace the placeholder paths as in the other examples.
set -eu
mkdir -p results
python3 -m quant_eval.downstream --name SC_blank \
  --design examples/sc/design.tsv --fasta data/sc/human.fasta --classes HUMAN \
  --andes andes:default:data/sc/andes/finalize.parquet \
  --maxquant MQ:data/sc/maxquant/evidence.txt \
  --fragpipe FP:data/sc/fragpipe/combined_ion.tsv \
  --out results/SC_blank.json
