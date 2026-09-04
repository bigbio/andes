#!/usr/bin/env bash
# Phase V CID decision benchmark: MSFragger vs andes (rank/strong) on
# 3 CID datasets × 4 raw files each. Pools per-dataset PINs, Percolator @1%.
# Comet is optional (SKIP_COMET=1 by default).
#
# Prereq: phase_v_cid_train_models.sh (intensity + scoring stores).
#
# Usage:
#   bash /srv/data/msgf-bench/phase_v_cid_benchmark.sh
#   DATASETS=PXD007683 bash ...   # run one dataset only
set -uo pipefail

cd /srv/data/msgf-bench

REPO="${REPO:-/srv/data/msgf-bench/repo/msgf-rust}"
BIN="${BIN:-$REPO/target/release/andes}"
BUNDLED="${BUNDLED:-$REPO/resources/ionstat/models.parquet}"
INTENSITY_MODEL="${INTENSITY_MODEL:-/srv/data/msnet/intensity_model_cid.parquet}"
SCORING_STORE="${SCORING_STORE:-/srv/data/msnet/models_cid_msnet.parquet}"
OUT="${OUT:-/srv/data/msgf-bench/cid-bench-$(date +%Y%m%d-%H%M)}"
DATASETS="${DATASETS:-PXD016999,PXD007683,PXD001819}"

FRAGGER="${FRAGGER:-/srv/data/msgf-bench/engines/msfragger-env/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar}"
TRFP="${TRFP:-/srv/data/.conda/envs/nextflow/bin/ThermoRawFileParser}"
COMET_RUN="${COMET_RUN:-/srv/data/msgf-bench/comet_cid_run.sh}"
SKIP_MSFRAGGER="${SKIP_MSFRAGGER:-0}"
SKIP_COMET="${SKIP_COMET:-1}"
OIMG="${OIMG:-ghcr.io/openms/openms-tools-thirdparty:latest}"
PIMG="${PIMG:-quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2}"
export DOTNET_ROOT="${DOTNET_ROOT:-/opt/dotnet8}"
export PATH="${DOTNET_ROOT}:${PATH}"

WORK=/tmp/cid-bench-work
mkdir -p "$WORK" "$OUT"
SUMMARY="$OUT/summary.tsv"
printf "dataset\tengine\tpsms_1pct\tpeptides_1pct\tfiles\n" > "$SUMMARY"

echo "=== BUILD $(date -Is) ==="
( cd "$REPO" && cargo +1.95.0 build --release --features thermo -p andes 2>&1 | tail -5 )
[ -x "$BIN" ] || { echo "missing $BIN"; exit 1; }
[ -f "$INTENSITY_MODEL" ] || { echo "missing $INTENSITY_MODEL — run phase_v_cid_train_models.sh"; exit 1; }
[ -f "$SCORING_STORE" ] || { echo "missing $SCORING_STORE — run phase_v_cid_train_models.sh"; exit 1; }
if [ "$SKIP_COMET" != 1 ]; then
  [ -x "$COMET_RUN" ] || { echo "missing $COMET_RUN"; exit 1; }
fi

count_pin() {
  local tp="$1"
  local qcol pcol
  qcol=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value"){print i;exit}}' "$tp")
  pcol=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="peptide"){print i;exit}}' "$tp")
  local psms peps
  psms=$(awk -F"\t" -v q="$qcol" 'NR>1&&$q<=0.01{c++}END{print c+0}' "$tp")
  peps=$(awk -F"\t" -v q="$qcol" -v p="$pcol" \
    'NR>1&&$q<=0.01{s=$p;gsub(/^[A-Z-]\./,"",s);gsub(/\.[A-Z-]$/,"",s);gsub(/\[[^]]*\]/,"",s);gsub(/[^A-Z]/,"",s);print s}' \
    "$tp" | sort -u | wc -l)
  echo "$psms $peps"
}

percolate_pool() {
  local ds="$1" eng="$2" pooled="$3"
  local tps="$OUT/${ds}.${eng}.t.psms" dps="$OUT/${ds}.${eng}.d.psms"
  docker run --rm --platform linux/amd64 -v "$OUT":/r "$PIMG" \
    percolator --seed 42 -Y \
      --results-psms "/r/$(basename "$tps")" \
      --decoy-results-psms "/r/$(basename "$dps")" \
      --only-psms=false "/r/$(basename "$pooled")" \
    > "$OUT/${ds}.${eng}.perc.log" 2>&1
}

pool_engine() {
  local ds="$1" eng="$2" nfiles="${3:-0}"
  local dir="$OUT/pins/$ds"
  local pooled="$OUT/${ds}.${eng}.pooled.pin"
  local batched="$dir/${eng}.pooled.pin"
  if [ -f "$batched" ]; then
    cp -f "$batched" "$pooled"
    echo "  $ds/$eng batched PIN $(($(wc -l < "$pooled")-1)) rows"
  else
    local files=("$dir"/"$eng".*.pin)
    [ -e "${files[0]}" ] || return 1
    head -1 "${files[0]}" > "$pooled"
    for f in "${files[@]}"; do tail -n +2 "$f" >> "$pooled"; done
    nfiles=${#files[@]}
    echo "  $ds/$eng pooled $(($(wc -l < "$pooled")-1)) rows from ${nfiles} files"
  fi
  percolate_pool "$ds" "$eng" "$pooled"
  local counts
  counts=$(count_pin "$OUT/${ds}.${eng}.t.psms")
  [ "$nfiles" -eq 0 ] && nfiles=1
  printf "%s\t%s\t%s\t%s\t%d\n" "$ds" "$eng" ${counts// /$'\t'} "$nfiles" >> "$SUMMARY"
  echo "  RESULT $ds/$eng PSMs@1%=${counts%% *} peptides@1%=${counts##* }"
}

run_comet_tmt() { bash "$COMET_RUN" tmt "$1" "$3" "$4"; }
run_comet_std() { bash "$COMET_RUN" std "$1" "$3" "$4"; }

download_raw() {
  local url="$1" dest="$2"
  echo "  [download] $url"
  curl -sS -L -C - -o "$dest" "$url"
  local sz
  sz=$(stat -c%s "$dest" 2>/dev/null || echo 0)
  if [ "$sz" -lt 1000000 ]; then
    echo "  ERROR: download too small (${sz} bytes): $dest" >&2
    rm -f "$dest"
    return 1
  fi
  echo "  [download] OK $(du -h "$dest" | cut -f1)"
}

ensure_mzml() {
  local raw="$1" mzml="$2"
  if [ -f "$mzml" ]; then return 0; fi
  [ -f "$raw" ] || return 1
  "$TRFP" -i "$raw" -b "$mzml" -f 2 -l 2
  [ -s "$mzml" ]
}

run_andes_dataset() {
  local ds="$1" fasta="$2" mods="$3"
  shift 3
  local extra_andes=("$@")
  local pin_dir="$OUT/pins/$ds"
  mkdir -p "$pin_dir"

  local spec_args=()
  local n=0
  for raw in "${ANDES_RAWS[@]}"; do
    [ -f "$raw" ] || continue
    spec_args+=(--spectrum "$raw")
    n=$((n + 1))
  done
  [ "$n" -gt 0 ] || { echo "  [andes] no raw inputs for $ds"; return 1; }

  echo "  [andes rank] $n file(s) in one run"
  "$BIN" "${spec_args[@]}" --database "$fasta" --mods "$mods" \
    --output-pin "$pin_dir/andes_rank.pooled.pin" --score rank \
    --fragmentation CID --instrument low-res \
    --model-store "$SCORING_STORE" --model cid_msnet_tryp \
    --intensity-model "$INTENSITY_MODEL" \
    --enzyme-specificity fully --max-missed-cleavages 2 --min-length 6 --max-length 40 \
    --charge-min 2 --charge-max 4 --top-n 1 --min-peaks 10 --threads 8 \
    "${extra_andes[@]}" \
    > "$OUT/${ds}.andes_rank.log" 2>&1 || true

  echo "  [andes strong] $n file(s) in one run"
  "$BIN" "${spec_args[@]}" --database "$fasta" --mods "$mods" \
    --output-pin "$pin_dir/andes_strong.pooled.pin" --score strong \
    --fragmentation CID --instrument low-res \
    --model-store "$SCORING_STORE" --model cid_msnet_tryp \
    --intensity-model "$INTENSITY_MODEL" \
    --enzyme-specificity fully --max-missed-cleavages 2 --min-length 6 --max-length 40 \
    --charge-min 2 --charge-max 4 --top-n 1 --min-peaks 10 --threads 8 \
    "${extra_andes[@]}" \
    > "$OUT/${ds}.andes_strong.log" 2>&1 || true
}

process_file_baseline() {
  local ds="$1" tag="$2" raw="$3" mzml="$4" fasta="$5" _td_fasta="$6" mods="$7"
  local frg_params="$8" comet_fn="$9"
  local pin_dir="$OUT/pins/$ds"
  mkdir -p "$pin_dir"

  echo "=== $ds / $tag $(date -Is) ==="
  local need_mzml=0
  [ "$SKIP_MSFRAGGER" != 1 ] && need_mzml=1
  [ "$SKIP_COMET" != 1 ] && need_mzml=1
  if [ "$need_mzml" = 1 ]; then
    ensure_mzml "$raw" "$mzml" || { echo "  mzML missing for $tag"; return 1; }
  fi

  if [ "$SKIP_MSFRAGGER" != 1 ]; then
    echo "  [msfragger]"
    local frg_mzml="$WORK/${tag}_frg.mzML"
    cp -f "$mzml" "$frg_mzml"
    java -Xmx8g -jar "$FRAGGER" "$frg_params" "$frg_mzml" > "$OUT/${ds}.${tag}.fragger.log" 2>&1 || true
    [ -f "$WORK/${tag}_frg.pin" ] && cp -f "$WORK/${tag}_frg.pin" "$pin_dir/msfragger.${tag}.pin"
    [ -f "$WORK/${tag}.pin" ] && cp -f "$WORK/${tag}.pin" "$pin_dir/msfragger.${tag}.pin"
    rm -f "$frg_mzml" "$WORK/${tag}_frg.pin" "$WORK/${tag}.pin" "$WORK/${tag}_frg.pepXML" "$WORK/${tag}_frg.tsv"
  else
    echo "  [msfragger] skip (SKIP_MSFRAGGER=1)"
  fi

  if [ "$SKIP_COMET" != 1 ]; then
    echo "  [comet]"
    "$comet_fn" "$mzml" "$tag" "$fasta" "$pin_dir/comet.${tag}.pin" \
      > "$OUT/${ds}.${tag}.comet.log" 2>&1 || true
  else
    echo "  [comet] skip (SKIP_COMET=1)"
  fi

  rm -f "$mzml"
  echo "  done $tag free=$(df -h /srv/data | tail -1 | awk '{print $4}')"
}

run_pxd016999() {
  local ds=PXD016999
  local base=/srv/data/msgf-bench/pxd016999
  local urlbase=https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/11/PXD016999
  local files=(
    "Instrument1_sample01_121115_Fr01|I1Fr01"
    "Instrument1_sample01_121115_Fr02|I1Fr02"
    "SecondInstrument_Sample10_S2R6_061516_Fr01|S2Fr01"
    "SecondInstrument_Sample10_S2R6_061516_Fr03|S2Fr03"
  )
  ANDES_RAWS=()
  local entries=()
  for entry in "${files[@]}"; do
    IFS='|' read -r name tag <<< "$entry"
    local raw="$WORK/$name.raw"
    download_raw "$urlbase/$name.raw" "$raw"
    ANDES_RAWS+=("$raw")
    entries+=("$tag|$raw|$WORK/$tag.mzML")
  done
  echo "=== $ds andes batch $(date -Is) ==="
  run_andes_dataset "$ds" "$base/human_sp.fasta" "$base/mods_msgf.txt" \
    --protocol TMT --precursor-tol-ppm 20 --isotope-error-min -1 --isotope-error-max 2
  for entry in "${entries[@]}"; do
    IFS='|' read -r tag raw mzml <<< "$entry"
    process_file_baseline "$ds" "$tag" "$raw" "$mzml" \
      "$base/human_sp.fasta" "$base/human_sp_td.fasta" "$base/mods_msgf.txt" \
      "$base/fragger_tmt.params" run_comet_tmt
    rm -f "$raw"
  done
  pool_all_engines "$ds" "${#files[@]}"
}

pool_all_engines() {
  local ds="$1"
  local nfiles="${2:-4}"
  local engines=(msfragger andes_rank andes_strong)
  [ "$SKIP_COMET" != 1 ] && engines=(msfragger comet andes_rank andes_strong)
  for eng in "${engines[@]}"; do
    if [[ "$eng" == andes_* ]]; then
      pool_engine "$ds" "$eng" 1 || true
    else
      pool_engine "$ds" "$eng" "$nfiles" || true
    fi
  done
}

run_pxd007683() {
  local ds=PXD007683
  local base=/srv/data/msgf-bench/tmt-data
  local urlbase=https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD007683
  local files=(a05058 a05059 a05060 a05061)
  local frg="$base/fragger_a05058.params"
  [ -f "$frg" ] || frg=/srv/data/msgf-bench/repo/msgf-rust/docs/benchmarks/configs/msfragger-a05058.params
  ANDES_RAWS=()
  local entries=()
  for name in "${files[@]}"; do
    local tag="$name" raw="$WORK/$name.raw" mzml="$WORK/$name.mzML"
    download_raw "$urlbase/$name.raw" "$raw"
    ANDES_RAWS+=("$raw")
    entries+=("$tag|$raw|$mzml")
  done
  echo "=== $ds andes batch $(date -Is) ==="
  run_andes_dataset "$ds" \
    "$base/PXD007683_UP000005640_UP000002311_reviewed.fasta" \
    "$base/mods-numeric.txt" \
    --protocol TMT --precursor-tol-ppm 20 --isotope-error-min -1 --isotope-error-max 2
  for entry in "${entries[@]}"; do
    IFS='|' read -r tag raw mzml <<< "$entry"
    if [ "$tag" = "a05058" ] && [ -f "$base/a05058.mzML" ]; then
      cp -f "$base/a05058.mzML" "$mzml"
    fi
    process_file_baseline "$ds" "$tag" "$raw" "$mzml" \
      "$base/PXD007683_UP000005640_UP000002311_reviewed.fasta" \
      "$base/PXD007683_reviewed.td.fasta" \
      "$base/mods-numeric.txt" \
      "$frg" run_comet_tmt
    rm -f "$raw"
  done
  pool_all_engines "$ds" "${#files[@]}"
}

run_pxd001819() {
  local ds=PXD001819
  local base=/srv/data/msgf-bench/data
  local urlbase=https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/12/PXD001819
  local mods=/srv/data/msgf-bench/mods.txt
  [ -f "$mods" ] || mods="$base/../mods.txt"
  local frg="$OUT/fragger_ups1.params"
  cat > "$frg" <<EOF
database_name = $base/PXD001819_uniprot_yeast_ups.revCat.fasta
num_threads = 8
precursor_mass_lower = -5
precursor_mass_upper = 5
precursor_mass_units = 1
fragment_mass_tolerance = 0.4
fragment_mass_units = 0
isotope_error = 0/1
search_enzyme_name_1 = stricttrypsin
allowed_missed_cleavage_1 = 2
num_enzyme_termini = 2
variable_mod_01 = 15.9949 M 3
add_C_cysteine = 57.02146
output_format = pin
output_report_topN = 2
minimum_peaks = 10
digest_min_length = 6
digest_max_length = 40
precursor_charge = 2 4
EOF
  local files=(UPS1_5000amol_R1 UPS1_5000amol_R2 UPS1_5000amol_R3)
  ANDES_RAWS=()
  local entries=()
  for name in "${files[@]}"; do
    local tag="$name" raw="$WORK/$name.raw" mzml="$WORK/$name.mzML"
    download_raw "$urlbase/$name.raw" "$raw"
    ANDES_RAWS+=("$raw")
    entries+=("$tag|$raw|$mzml")
  done
  echo "=== $ds andes batch $(date -Is) ==="
  run_andes_dataset "$ds" \
    "$base/PXD001819_uniprot_yeast_ups.fasta" \
    "$mods" \
    --precursor-tol-ppm 5 --isotope-error-min 0 --isotope-error-max 1
  for entry in "${entries[@]}"; do
    IFS='|' read -r tag raw mzml <<< "$entry"
    if [ "$tag" = "UPS1_5000amol_R1" ] && [ -f "$base/UPS1_5000amol_R1.mzML" ]; then
      cp -f "$base/UPS1_5000amol_R1.mzML" "$mzml"
    fi
    process_file_baseline "$ds" "$tag" "$raw" "$mzml" \
      "$base/PXD001819_uniprot_yeast_ups.fasta" \
      "$base/PXD001819_uniprot_yeast_ups.revCat.fasta" \
      "$mods" "$frg" run_comet_std
    rm -f "$raw"
  done
  pool_all_engines "$ds" "${#files[@]}"
}

IFS=',' read -ra WANT <<< "$DATASETS"
for ds in "${WANT[@]}"; do
  case "$ds" in
    PXD016999) run_pxd016999 ;;
    PXD007683) run_pxd007683 ;;
    PXD001819) run_pxd001819 ;;
    *) echo "unknown dataset $ds"; exit 1 ;;
  esac
done

echo "=== SUMMARY $(date -Is) ==="
column -t "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
python3 /srv/data/msgf-bench/summarize_phase_v_cid_benchmark.py "$OUT" || true
echo "OUT=$OUT"
echo "=== PHASE_V_CID_BENCHMARK_DONE $(date -Is) ==="
