#!/usr/bin/env bash
# Simulated-substitution benchmark, following Yu et al., Nat. Commun. 2020
# ("Identification of modified peptides using localization-aware open search").
#
#   ./run_substitution_bench.sh <data-dir> [mzML]
#
# The paper builds a "modified" database by introducing random amino-acid
# substitutions and searches the UNMODIFIED spectra against it with open search:
# the mass offset is then the substitution's mass difference, and the
# substituted residue is the ground truth for localization.
#
# The PXD001468 download in the benchmark data directory already carries such a
# database (`UP000005640_9606_modified.fasta`, 58,015 substitutions over 8,832
# proteins), so this script only reproduces the search and the evaluation:
#
#   reference  narrow search of the spectra against the plain human database
#              (the paper's variable mods: M oxidation, protein N-term acetyl;
#              C carbamidomethyl fixed). Its PSMs are the truth peptides.
#   open       open search of the same spectra against the substituted database
#              (fixed mods only), offset window -150..150 Da.
#
# Needs: release `andes` on PATH (or $ANDES), Docker for Percolator, python3.
set -euo pipefail

DATA="${1:?usage: run_substitution_bench.sh <data-dir> [mzML]}"
FILE="${2:-$DATA/PXD001468/b1906_293T_proteinID_01A_QE3_122212.mzML}"
PLAIN="$DATA/PXD001468/UP000005640_9606.fasta"
MUT="$DATA/PXD001468/UP000005640_9606_modified.fasta"
OUT="${OUT:-$DATA/results-substitution}"
ANDES="${ANDES:-andes}"
THREADS="${THREADS:-112}"
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUT"
printf '57.021464,C,fix,any,Carbamidomethyl\n' > "$OUT/fixed.mods"

search() {  # $1 = name, $2 = database, rest = extra flags
  local name="$1" db="$2"; shift 2
  echo "search: $name ($db)"
  "$ANDES" --spectrum "$FILE" --database "$db" --threads "$THREADS" \
    --output-pin "$OUT/$name.pin" "$@" > "$OUT/$name.log" 2>&1
  grep -E "Search complete|open search:" "$OUT/$name.log" | sed 's/^/  /'
}

search reference "$PLAIN"
search open "$MUT" --mods "$OUT/fixed.mods" --open-search --open-search-delta=-150..150

for name in reference open; do
  echo "percolator: $name"
  docker run --rm -v "$OUT":/r "$PIMG" percolator --seed 42 -Y --only-psms=false \
    --results-psms "/r/$name.psms" --decoy-results-psms "/r/$name.dpsms" \
    "/r/$name.pin" > "$OUT/$name.percolator.log" 2>&1
done

python3 -I "$HERE/analyze_substitution.py" \
  --dir "$OUT" --plain "$PLAIN" --mutated "$MUT" --out "$OUT/report.json"
