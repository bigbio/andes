#!/usr/bin/env bash
# Open-search benchmark: oxidised-M recovery and entrapment error on one
# mouse-liver fraction.
#
#   ./run_open_bench.sh <data-dir> [raw-file]
#
# Needs: a release `andes` on PATH (or $ANDES) built with `--features thermo`
# for native `.raw`, Docker for Percolator, python3.
#
# Two searches of the same spectra:
#   truth   defaults (Carbamidomethyl-C fixed; Oxidation-M and N-term-acetyl
#           variable). Its oxidised PSMs are the ground truth for the open run.
#   open    fixed mods only + `--open-search`. Oxidation is then an unknown
#           offset, and the run reports it (plus its residue) instead of
#           enumerating it.
set -euo pipefail

DATA="${1:?usage: run_open_bench.sh <data-dir> [raw-file]}"
RAW="${2:-$DATA/glyco-mouse/MouseLiver-Z-T-1.raw}"
DB="${DB:-$DATA/databases/mouse_entrap.fasta}"
OUT="${OUT:-$DATA/results-open}"
ANDES="${ANDES:-andes}"
THREADS="${THREADS:-112}"
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUT"
printf '57.021464,C,fix,any,Carbamidomethyl\n' > "$OUT/fixed.mods"

run_andes() {  # $1 = name, rest = extra flags
  local name="$1"; shift
  echo "search: $name"
  "$ANDES" --spectrum "$RAW" --database "$DB" --threads "$THREADS" \
    --output-pin "$OUT/$name.pin" "$@" > "$OUT/$name.log" 2>&1
  echo "  $(grep -E 'Search complete' "$OUT/$name.log" | tail -1)"
}

run_andes truth
run_andes open --mods "$OUT/fixed.mods" --open-search

for name in truth open; do
  echo "percolator: $name"
  docker run --rm -v "$OUT":/r "$PIMG" percolator --seed 42 -Y --only-psms=false \
    --results-psms "/r/$name.psms" --decoy-results-psms "/r/$name.dpsms" \
    "/r/$name.pin" > "$OUT/$name.percolator.log" 2>&1
  # accepted PSMs, entrapment hits, I/L-rescued PSMs at q <= 0.01
  printf '  entrapment %s: ' "$name"
  python3 -I "$HERE/../entrap_common.py" count "$OUT/$name.psms" "$DB" 0.01
done

# T/E is a property of the database build; computing it scans every tryptic
# peptide of the FASTA (minutes on a mammalian database). Pass it in with
# TE=<ratio> when known, or analyze.py computes it once.
if [ -n "${TE:-}" ]; then
  python3 -I "$HERE/analyze.py" --dir "$OUT" --fasta "$DB" --te "$TE" --out "$OUT/report.json"
else
  python3 -I "$HERE/analyze.py" --dir "$OUT" --fasta "$DB" --out "$OUT/report.json"
fi
echo "report: $OUT/report.json"
