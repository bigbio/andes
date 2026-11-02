#!/bin/sh
# Absence check only (no amounts.tsv): yeast rows in the 20 human-only runs (absent_species in design.tsv). Replace the placeholder paths.
set -eu
mkdir -p results
python3 -m quant_eval.downstream --name PXD014415 \
  --design examples/p14/design.tsv --fasta data/PXD014415/human_yeast.fasta --classes HUMAN,YEAST \
  --andes andes:default:data/PXD014415/andes/finalize.parquet \
  --fragpipe FP:data/PXD014415/fragpipe/combined_ion.tsv \
  --out results/PXD014415.json
