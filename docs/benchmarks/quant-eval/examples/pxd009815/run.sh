#!/bin/sh
# Run from docs/benchmarks/quant-eval. Replace the placeholder paths below (relative to the working directory) with the
# andes finalize output, the MaxQuant evidence.txt and the FragPipe combined_ion.tsv (IonQuant normalisation off) of this dataset,
# and the run names in design.tsv with the andes run stems if they differ.
set -eu
mkdir -p results
python3 -m quant_eval.downstream --name PXD009815 \
  --design examples/pxd009815/design.tsv --amounts examples/pxd009815/amounts.tsv \
  --fasta data/PXD009815/yeast_ups1.fasta --classes UPS,YEAST --null YEAST \
  --andes andes:default:data/PXD009815/andes/finalize.parquet \
  --andes andes_V:weights:data/PXD009815/andes/finalize.parquet \
  --andes andes_H:hard-cut:data/PXD009815/andes/finalize.parquet \
  --maxquant MQ:data/PXD009815/maxquant/evidence.txt \
  --fragpipe FP:data/PXD009815/fragpipe/combined_ion.tsv \
  --out results/PXD009815.json
