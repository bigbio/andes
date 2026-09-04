#!/usr/bin/env bash
# Run the benchmarks and print the results table.
#
#   ./run.sh <data-dir> [dataset ...]      # astral tmt ups1  (default: all three)
#
# Needs: an `andes` binary on PATH (or $ANDES), Docker for Percolator. No absolute paths,
# nothing specific to our hosts. Every number it prints carries its own provenance line.
set -euo pipefail

DATA_DIR="${1:?usage: run.sh <data-dir> [dataset ...]}"; shift || true
SETS=("$@"); [ ${#SETS[@]} -eq 0 ] && SETS=(astral tmt ups1)
ANDES="${ANDES:-andes}"; THREADS="${THREADS:-8}"
DB="$DATA_DIR/databases"; OUT="$DATA_DIR/results"; mkdir -p "$OUT"
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v "$ANDES" >/dev/null || { echo "andes not found (set \$ANDES)" >&2; exit 1; }
command -v docker  >/dev/null || { echo "docker not found (needed for Percolator)" >&2; exit 1; }

echo "provenance: andes $("$ANDES" --version 2>&1 | head -1) | $(uname -srm) | ${THREADS} threads | $(date -u +%FT%TZ)"
printf '\n%-8s %10s %12s %14s\n' dataset wall "PSMs@q0.01" "entrap hits"

for ds in "${SETS[@]}"; do
  case "$ds" in
    astral) spec=$(ls "$DATA_DIR"/astral/*.mzML "$DATA_DIR"/astral/*.raw 2>/dev/null | head -1)
            fa="$DB/hye.fasta";          extra=(--mods "$HERE/../configs/astral_mods.txt" --precursor-tol 10ppm --enzyme trypsin) ;;
    tmt)    spec=$(ls "$DATA_DIR"/tmt/*.mzML "$DATA_DIR"/tmt/*.raw 2>/dev/null | head -1)
            fa="$DB/tmt_db.fasta";       extra=(--mods "$HERE/../configs/mods-tmt.txt") ;;
    ups1)   spec=$(ls "$DATA_DIR"/ups1/*.mzML "$DATA_DIR"/ups1/*.raw 2>/dev/null | head -1)
            fa="$DB/yeast_entrap.fasta"; extra=() ;;
    *) echo "unknown dataset '$ds'" >&2; exit 1 ;;
  esac
  [ -n "${spec:-}" ] && [ -s "$spec" ] || { echo "$ds: no spectrum file - run fetch_spectra.sh" >&2; continue; }
  [ -s "$fa" ] || { echo "$ds: no database - run build_databases.sh" >&2; continue; }

  s=$(date +%s)
  "$ANDES" --spectrum "$spec" --database "$fa" "${extra[@]}" \
           --threads "$THREADS" --output-pin "$OUT/$ds.pin" > "$OUT/$ds.log" 2>&1
  e=$(date +%s)
  docker run --rm -v "$OUT":/r $PIMG percolator --seed 42 -Y --only-psms=false \
    --results-psms "/r/$ds.t.psms" --decoy-results-psms "/r/$ds.d.psms" "/r/$ds.pin" \
    > "$OUT/$ds.perc.log" 2>&1
  read -r n ent < <(awk -F'\t' '
    NR==1{for(i=1;i<=NF;i++) if($i=="q-value") q=i; p=NF; next}
    $q<=0.01{c++; for(i=p;i<=NF;i++) if($i ~ /ENTRAP_/){e++; break}}
    END{print c+0, e+0}' "$OUT/$ds.t.psms")
  printf '%-8s %9ss %12s %14s\n' "$ds" "$((e-s))" "$n" "$ent"
done

cat <<'NOTE'

Reading these numbers
  * PSMs@q0.01 is Percolator's claim. It is comparable across engines run through this
    same protocol; it is NOT a measured error rate.
  * A true FDP needs an entrapment database and its own scaling factor:
        FDP = (entrap hits / total) x (1 + T/E)
    build_databases.sh prints T/E for YOUR build. Never assume 1:1 (factor 2) -- doing so
    has understated the true error in this project twice.
  * Only ups1 has an entrapment component here, so only it yields a true FDP.
  * Wall times compare only within one host and one session. Counts are portable.
NOTE
