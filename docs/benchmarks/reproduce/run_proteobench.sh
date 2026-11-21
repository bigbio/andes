#!/usr/bin/env bash
set -euo pipefail

usage="usage: run_proteobench.sh <qe|astral> <database.fasta> <run.mzML ...>"
MODULE="${1:?$usage}"
shift
case "$MODULE" in
  qe | astral) ;;
  *) echo "$usage" >&2; exit 1 ;;
esac
FASTA="${1:?missing database fasta}"
shift
[ $# -ge 1 ] || { echo "$usage" >&2; exit 1; }
ANDES="${ANDES:-andes}"
THREADS="${THREADS:-8}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
command -v "$ANDES" >/dev/null || { echo "andes not found (set \$ANDES)" >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found" >&2; exit 1; }
echo "provenance: andes $("$ANDES" --version 2>&1 | head -1) | $(uname -srm) | ${THREADS} threads | $(date -u +%FT%TZ)"
mkdir -p q
PSMS=()
for spec in "$@"; do
  R="$(basename "$spec")"
  R="${R%.*}"
  "$ANDES" --spectrum "$spec" --database "$FASTA" --rescore --percolator-docker --output-pin "$R.pin" --quant-dir q/ --threads "$THREADS"
  PSMS+=(--psms "$R=$R.percolator.target.psms.txt" --decoy-psms "$R=$R.percolator.decoy.psms.txt")
done
"$ANDES" quant library --quant-dir q/ "${PSMS[@]}"
for spec in "$@"; do
  R="$(basename "$spec")"
  R="${R%.*}"
  "$ANDES" quant extract --quant-dir q/ --run "$R"
done
"$ANDES" quant finalize --quant-dir q/ -o pb
"$ANDES" quant finalize --quant-dir q/ -o pb05 --quant-fdr 0.05
python3 "$HERE/../proteobench/andes_to_proteobench.py" --mbr pb.mbr.parquet --quant-dir q/ \
      --module "$MODULE" --column quantity_normalized -o custom_input.tsv
