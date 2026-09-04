#!/usr/bin/env bash
# Phase V gate: A/B --score rank (baseline) vs --score strong on Astral, TMT, PTM.
#
# For each dataset × score mode:
#   - NORMAL db  → Percolator PSMs@1%
#   - ENTRAP db  → entrapment FDP @1%
#   - search wall time (Pass-1 gate: strong ≤ 110% of rank)
#
# Prereq: intensity model at INTENSITY_MODEL (run phase_v_train_intensity.sh).
#
# Usage:
#   bash /srv/data/msgf-bench/phase_v_strong_score_gate.sh
#   SKIP_PTM=1 bash ...   # Astral + TMT only
set -uo pipefail

cd /srv/data/msgf-bench

REPO="${REPO:-/srv/data/msgf-bench/repo/msgf-rust}"
BIN="${BIN:-$REPO/target/release/andes}"
BUNDLED="${BUNDLED:-$REPO/resources/ionstat/models.parquet}"
INTENSITY_MODEL="${INTENSITY_MODEL:-/srv/data/msnet/intensity_model.parquet}"
OUT="${OUT:-/srv/data/msgf-bench/phase-v-$(date +%Y%m%d-%H%M)}"
SKIP_PTM="${SKIP_PTM:-0}"

export DOTNET_ROOT="${DOTNET_ROOT:-/opt/dotnet8}"
export PATH="${DOTNET_ROOT}:${PATH}"

mkdir -p "$OUT"
SUMMARY="$OUT/summary.tsv"
printf "dataset\tscore_mode\tdb\tpsms_1pct\tentrap_fdp\tentrap_combined_fdp\twall_s\tmax_rss_kb\tpin_rows\n" > "$SUMMARY"

psm1() {
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")q=i;next} q&&$q<=0.01{c++} END{print c+0}' "$1" 2>/dev/null
}

wall_s() {
  grep -E 'match_spectra wall:|Elapsed \(wall clock time\)' "$1" 2>/dev/null | tail -1 | awk -F': ' '{gsub(/s/,"",$NF); print $NF}'
}

max_rss_kb() {
  grep 'Maximum resident set size' "$1" 2>/dev/null | awk -F': ' '{print $NF}'
}

echo "=== BUILD $(date -Is) ==="
( cd "$REPO" && cargo +1.95.0 build --release --features thermo -p andes 2>&1 | tail -5 )
[ -x "$BIN" ] || { echo "missing $BIN"; exit 1; }
[ -f "$INTENSITY_MODEL" ] || { echo "missing intensity model: $INTENSITY_MODEL (run phase_v_train_intensity.sh)"; exit 1; }

COMMON="--model-store $BUNDLED --model hcd_qexactive_tryp --intensity-model $INTENSITY_MODEL \
  --enzyme-specificity fully --max-missed-cleavages 2 --min-length 7 --max-length 40 \
  --charge-min 2 --charge-max 4 --top-n 1 --min-peaks 10 --threads 8"

run_one() {
  local ds="$1" mode="$2" dbkind="$3" tag="$4"
  local spectra="$5" fasta="$6" mods="$7"
  shift 7
  local extra=("$@")

  local score_flag=(--score rank)
  [ "$mode" = "strong" ] && score_flag=(--score strong)

  local pin="$OUT/${tag}.pin"
  local log="$OUT/${tag}.log"
  echo "=== $tag $(date -Is) ==="
  /usr/bin/time -v "$BIN" \
    --spectrum "$spectra" \
    --database "$fasta" \
    --output-pin "$pin" \
    --mods "$mods" \
    "${score_flag[@]}" \
    $COMMON \
    "${extra[@]}" \
    > "$log" 2>&1
  local rc=$?
  local rows=0
  [ -f "$pin" ] && rows=$(($(wc -l < "$pin") - 1))
  echo "  exit=$rc rows=$rows wall=$(wall_s "$log")s"

  local psms=0 fdp="" cfdp=""
  if [ "$dbkind" = "normal" ]; then
    bash run_percolator_docker.sh "$pin" "$OUT" "$tag" > "$OUT/${tag}.perc.log" 2>&1
    psms=$(psm1 "$OUT/${tag}.target.psms.txt")
    echo "  PSMs@1%=$psms"
  else
    bash run_percolator_docker.sh "$pin" "$OUT" "$tag" > "$OUT/${tag}.perc.log" 2>&1
    local fdp_line
    fdp_line=$(python3 compute_entrapment_fdp.py "$OUT/${tag}.target.psms.txt" 0.01 "$tag" 2>/dev/null | head -1)
    echo "  $fdp_line"
    fdp=$(echo "$fdp_line" | sed -n 's/.*entrapment_fraction=\([0-9.]*\).*/\1/p')
    cfdp=$(echo "$fdp_line" | sed -n 's/.*combined_FDP=\([0-9.]*\).*/\1/p')
    psms=$(psm1 "$OUT/${tag}.target.psms.txt")
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$ds" "$mode" "$dbkind" "$psms" "$fdp" "$cfdp" "$(wall_s "$log")" "$(max_rss_kb "$log")" "$rows" \
    >> "$SUMMARY"
}

# ── Astral (high-res, thermo RAW) ───────────────────────────────────────────
AD=astral-data
ASTRAL_EXTRA="--precursor-tol-ppm 10 --isotope-error-min -1 --isotope-error-max 2"
for mode in rank strong; do
  run_one astral "$mode" normal "astral-${mode}-normal" \
    "$AD/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw" \
    "$AD/ProteoBenchFASTA_MixedSpecies_HYE.fasta" \
    astral_mods_rust.txt $ASTRAL_EXTRA
  run_one astral "$mode" entrap "astral-${mode}-entrap" \
    "$AD/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw" \
    "$AD/ASTRAL_entrapment.fasta" \
    astral_mods_rust.txt $ASTRAL_EXTRA
done

# ── TMT (low-res CID) ───────────────────────────────────────────────────────
TD=tmt-data
TMT_EXTRA="--precursor-tol-ppm 20 --isotope-error-min -1 --isotope-error-max 2 \
  --fragmentation CID --instrument low-res --protocol TMT"
for mode in rank strong; do
  run_one tmt "$mode" normal "tmt-${mode}-normal" \
    "$TD/a05058.mzML" \
    "$TD/PXD007683_UP000005640_UP000002311_reviewed.fasta" \
    "$TD/mods-numeric.txt" $TMT_EXTRA
  run_one tmt "$mode" entrap "tmt-${mode}-entrap" \
    "$TD/a05058.mzML" \
    "$TD/TMT_entrapment.fasta" \
    "$TD/mods-numeric.txt" $TMT_EXTRA
done

# ── PTM (phospho) — optional until fixtures are staged ─────────────────────
if [ "$SKIP_PTM" != "1" ] && [ -n "${PTM_MZML:-}" ] && [ -n "${PTM_FASTA:-}" ] && [ -n "${PTM_MODS:-}" ]; then
  PTM_EXTRA="${PTM_EXTRA:---precursor-tol-ppm 10 --isotope-error-min -1 --isotope-error-max 2 \
    --fragmentation HCD --instrument QExactive --protocol phospho}"
  PTM_ENTRAP_DB="${PTM_ENTRAP:-${PTM_FASTA%.fasta}_entrapment.fasta}"
  for mode in rank strong; do
    run_one ptm "$mode" normal "ptm-${mode}-normal" \
      "$PTM_MZML" "$PTM_FASTA" "$PTM_MODS" $PTM_EXTRA
    if [ -f "$PTM_ENTRAP_DB" ]; then
      run_one ptm "$mode" entrap "ptm-${mode}-entrap" \
        "$PTM_MZML" "$PTM_ENTRAP_DB" "$PTM_MODS" $PTM_EXTRA
    else
      echo "WARN: PTM entrap FASTA missing ($PTM_ENTRAP_DB); skipping entrap arm"
    fi
  done
else
  echo "WARN: PTM skipped (set PTM_MZML/PTM_FASTA/PTM_MODS or SKIP_PTM=0 after staging)"
fi

echo "=== SUMMARY $(date -Is) ==="
column -t -s $'\t' "$SUMMARY"
python3 /srv/data/msgf-bench/summarize_phase_v_gate.py "$OUT"
echo "OUT=$OUT"
echo "=== PHASE_V_DONE $(date -Is) ==="
