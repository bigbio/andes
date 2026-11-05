#!/bin/sh
# PXD001819, 27 runs, 9 UPS1 levels (level is confounded with run order). Use the suite manifest (quant_eval.suite) for the full comparison.
set -eu
python3 -m quant_eval.downstream --name PXD001819 --design examples/ups27/design.tsv --amounts examples/ups27/amounts.tsv \
  --fasta data/PXD001819/yeast_ups1.fasta --classes UPS,YEAST --null YEAST \
  --andes andes:default:data/PXD001819/andes/finalize.parquet --maxquant MQ:data/PXD001819/maxquant/evidence.txt \
  --fragpipe FP:data/PXD001819/fragpipe/combined_ion.tsv --out results/PXD001819.json
