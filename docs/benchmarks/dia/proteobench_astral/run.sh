#!/bin/bash
# ProteoBench "Quant LFQ DIA ion" Astral module (Module 7), six runs as one andes study, both
# quantity arms scored with the proteobench package.
#
#   run.sh <andes> <raw_dir> <fasta> <out_dir> <python> <proteobench_dir> [threads]
#
# <raw_dir>: the six LFQ_Astral_DIA_15min_50ng_Condition_{A,B}_REP{1,2,3}.raw (ProteoBench server);
# <fasta>: ProteoBenchFASTA_MixedSpecies_HYE.fasta; <python>: Python >= 3.11 with pyarrow;
# <proteobench_dir>: `pip install --target <dir> proteobench==0.17.6 "plotly>=6,<7" "pandas>=2,<3"`.
# Reading .raw needs the .NET 8 runtime (DOTNET_ROOT). The six searches and extracts are
# independent; run them in parallel where possible.
set -euo pipefail
ANDES=$1 RAW=$2 FASTA=$3 OUT=$4 PY=$5 PB=$6 THREADS=${7:-14}
HERE=$(cd "$(dirname "$0")" && pwd)
RUNS=(LFQ_Astral_DIA_15min_50ng_Condition_A_REP1 LFQ_Astral_DIA_15min_50ng_Condition_A_REP2
      LFQ_Astral_DIA_15min_50ng_Condition_A_REP3 LFQ_Astral_DIA_15min_50ng_Condition_B_REP1
      LFQ_Astral_DIA_15min_50ng_Condition_B_REP2 LFQ_Astral_DIA_15min_50ng_Condition_B_REP3)
mkdir -p "$OUT"
for s in "${RUNS[@]}"; do
  "$ANDES" --spectrum "$RAW/$s.raw" --database "$FASTA" \
    --enzyme trypsin --max-missed-cleavages 1 --min-length 6 --max-length 30 --charge 1..5 \
    --fixed-mods "Carbamidomethyl (C)" --variable-mods "Oxidation (M)" --max-mods 1 \
    --dia --model hcd_astral_tryp --predict rt,ms2 --rescore --threads "$THREADS" \
    --quant-dir "$OUT/q" --output-pin "$OUT/$s.pin" > "$OUT/$s.search.log" 2>&1
done
for s in "${RUNS[@]}"; do
  "$ANDES" quant extract --quant-dir "$OUT/q" --run "$s" --threads "$THREADS" > "$OUT/$s.extract.log" 2>&1
done
for seed in 42 1 2; do
  F=$OUT/fin_s$seed; mkdir -p "$F"
  "$ANDES" quant finalize --quant-dir "$OUT/q" -o "$F/pb7" --database "$FASTA" --seed $seed \
    --threads "$THREADS" > "$F/fin.log" 2>&1
  PYTHONPATH=$PB "$PY" "$HERE/arms.py" "$F/pb7.dia_quant.parquet" "$F/arms"
  for arm in frag ms1; do
    PYTHONPATH=$PB "$PY" "$HERE/../../proteobench/andes_to_proteobench.py" \
      --dia-quant "$F/arms/$arm.dia_quant.parquet" --module astral_dia --column quantity_normalized \
      -o "$F/arms/custom_$arm.tsv"
    PYTHONPATH=$PB "$PY" "$HERE/score.py" "$F/arms/custom_$arm.tsv" "seed $seed $arm"
  done
done
