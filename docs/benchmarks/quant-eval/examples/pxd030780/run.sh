#!/bin/sh
# PXD030780, 25 runs, 5 HYE groups A-E. No longer held out: context only.
set -eu
python3 -m quant_eval.downstream --name PXD030780 --design examples/pxd030780/design.tsv --amounts examples/pxd030780/amounts.tsv \
  --fasta data/PXD030780/target.fasta --classes HUMAN,YEAST,ECOLI --null HUMAN \
  --andes andes:default:data/PXD030780/andes/finalize.parquet --fragpipe FP:data/PXD030780/fragpipe/combined_ion.tsv --out results/PXD030780.json
