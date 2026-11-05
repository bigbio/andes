#!/bin/sh
# Run from docs/benchmarks/quant-eval. Replace the placeholder paths below (relative to the working directory) with the
# andes finalize output, the MaxQuant evidence.txt and the FragPipe combined_ion.tsv (IonQuant normalisation off) of this dataset,
# and the run names in design.tsv with the andes run stems if they differ.
set -eu
mkdir -p results
python3 -m quant_eval.downstream --name ProteoBench_Astral \
  --design examples/proteobench_astral/design.tsv --amounts examples/proteobench_astral/amounts.tsv \
  --fasta data/ProteoBench_Astral/hye.fasta --classes HUMAN,YEAST,ECOLI --null HUMAN \
  --andes andes:default:data/ProteoBench_Astral/andes/finalize.parquet \
  --andes andes_V:weights:data/ProteoBench_Astral/andes/finalize.parquet \
  --andes andes_H:hard-cut:data/ProteoBench_Astral/andes/finalize.parquet \
  --maxquant MQ:data/ProteoBench_Astral/maxquant/evidence.txt \
  --fragpipe FP:data/ProteoBench_Astral/fragpipe/combined_ion.tsv \
  --out results/ProteoBench_Astral.json
