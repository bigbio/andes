#!/usr/bin/env bash
# Phase V CID model training: add 3 CID MSNet datasets to intensity + scoring stores.
#
# Intensity: merge existing HCD partials (if present) with 3 new CID partials.
# Scoring: train-from-msnet → cid_msnet_tryp in an extended model store.
#
# Usage:
#   bash /srv/data/msgf-bench/phase_v_cid_train_models.sh
#
# Output:
#   /srv/data/msnet/intensity_model_cid.parquet
#   /srv/data/msnet/models_cid_msnet.parquet
set -euo pipefail

cd /srv/data/msnet
mkdir -p raw partial flat

SCRIPT_DIR="${SCRIPT_DIR:-/srv/data/msgf-bench}"
if [ "${PREPARE_LUMOS:-1}" = 1 ] && [ -f "$SCRIPT_DIR/prepare_lumos_cid_training.sh" ]; then
  echo "=== LUMOS SELF-LABEL PREP $(date -Is) ==="
  bash "$SCRIPT_DIR/prepare_lumos_cid_training.sh"
fi

REPO="${REPO:-/srv/data/msgf-bench/repo/msgf-rust}"
BIN="${BIN:-$REPO/target/release/andes}"
AGG="${AGG:-/srv/data/msnet/msnet_intensity_agg.py}"
FLAT="${FLAT:-/srv/data/msnet/msnet_to_flat.py}"
INTENSITY_OUT="${INTENSITY_OUT:-/srv/data/msnet/intensity_model_cid.parquet}"
SCORING_OUT="${SCORING_OUT:-/srv/data/msnet/models_cid_msnet.parquet}"
BUNDLED="${BUNDLED:-$REPO/resources/ionstat/models.parquet}"
PEP_MAX="${PEP_MAX:-0.001}"
SAMPLE="${SAMPLE:-300000}"
THREADS="${THREADS:-8}"
TRAIN_MIN_MATCHED="${TRAIN_MIN_MATCHED:-6}"
TRAIN_MIN_EXPLAINED="${TRAIN_MIN_EXPLAINED:-0.50}"
# Sharper estimator than defaults (see internal-docs/2026-06-07-cid-rank-lumos-plan.md).
TRAIN_PSEUDO="${TRAIN_PSEUDO:-0.5}"
TRAIN_BACKOFF="${TRAIN_BACKOFF:-12}"
TRAIN_MIN_COUNT="${TRAIN_MIN_COUNT:-50}"
# Extra sample rows for large mixed HCD+CID parquets (CID filter keeps ~1–2% of rows).
SAMPLE_865="${SAMPLE_865:-800000}"
export THREADS FRAG_FILTER=cid

# dataset_id|url|sample_override (optional third field)
CID_DATASETS=(
  "PXD000865-CID|https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/msnet/PXD000865/PXD000865-MSNet.parquet|${SAMPLE_865}"
  "PXD016986|https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/msnet/PXD016986/PXD016986-MSNet.parquet|${SAMPLE}"
  "PXD009449-CID|https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/msnet/PXD009449-CID/PXD009449-CID-MSNet.parquet|${SAMPLE}"
)

echo "=== BUILD andes $(date -Is) ==="
( cd "$REPO" && cargo +1.95.0 build --release -p andes 2>&1 | tail -5 )
[ -x "$BIN" ] || { echo "missing $BIN"; exit 1; }

download_with_progress() {
  local id="$1" url="$2" dest="$3"
  echo "[$id] download start $(date -Is)"
  curl -sS -L -o "$dest" "$url" &
  local cpid=$!
  local last_mb=-1
  while kill -0 "$cpid" 2>/dev/null; do
    if [ -f "$dest" ]; then
      local mb
      mb=$(du -m "$dest" 2>/dev/null | cut -f1)
      if [ -n "$mb" ] && [ "$mb" != "$last_mb" ]; then
        echo "[$id] download ${mb}MB..."
        last_mb=$mb
      fi
    fi
    sleep 3
  done
  wait "$cpid"
  echo "[$id] download done $(du -h "$dest" | cut -f1) $(date -Is)"
}

stage_cid_one() {
  local idx="$1" total="$2" id="$3" url="$4" sample="$5"
  local raw="raw/${id}.parquet"
  local partial="partial/${id}.parquet"
  local flat_out="flat/${id}.parquet"
  echo "[dataset ${idx}/${total}] [$id] sample=$sample $(date -Is)"

  if [ "${REFLAT:-0}" = 1 ]; then
    rm -f "$partial" "$flat_out" "$raw"
  fi

  if [ ! -f "$partial" ]; then
    download_with_progress "$id" "$url" "$raw"
    echo "[$id] intensity aggregate pep_max=$PEP_MAX sample=$sample frag=cid"
    python3 "$AGG" "$raw" "$partial" "$PEP_MAX" "$sample" cid
    rm -f "$raw"
    echo "[$id] partial keys=$(python3 -c "import pyarrow.parquet as pq; print(pq.read_table('$partial').num_rows)")"
  else
    echo "[$id] partial exists, skip aggregate"
  fi

  local mm="$TRAIN_MIN_MATCHED" me="$TRAIN_MIN_EXPLAINED"
  case "$id" in
    PXD016986) me=0.45 ;;
  esac
  if [ ! -f "$flat_out" ] || [[ "${REFLAT_IDS:-}" == *"$id"* ]] || [ "${REFLAT:-0}" = 1 ]; then
    [ -f "$raw" ] || download_with_progress "$id" "$url" "$raw"
    rm -f "$flat_out"
    echo "[$id] flat convert sample=$sample frag=cid min_matched=$mm min_explained=$me"
    python3 "$FLAT" "$raw" "$flat_out" "$PEP_MAX" "$sample" 0 "$mm" "$me" 0 cid
    rm -f "$raw"
    echo "[$id] flat rows=$(python3 -c "import pyarrow.parquet as pq; print(pq.read_table('$flat_out').num_rows)")"
  else
    echo "[$id] flat exists, skip convert"
  fi
  echo "[$id] free=$(df -h /srv/data | tail -1 | awk '{print $4}')"
}

TOTAL=${#CID_DATASETS[@]}
IDX=0
for entry in "${CID_DATASETS[@]}"; do
  IDX=$((IDX + 1))
  IFS='|' read -r id url sample <<< "$entry"
  sample="${sample:-$SAMPLE}"
  stage_cid_one "$IDX" "$TOTAL" "$id" "$url" "$sample"
done

# Merge HCD partials from Phase T (if present) with new CID partials.
PARTIALS=(partial/*.parquet)
echo "=== MERGE ${#PARTIALS[@]} intensity partials $(date -Is) ==="
python3 "$AGG" --merge "${PARTIALS[@]}" -o "$INTENSITY_OUT.partial_merged.parquet"

echo "=== train-intensity $(date -Is) ==="
"$BIN" train-intensity \
  --in "$INTENSITY_OUT.partial_merged.parquet" \
  --out "$INTENSITY_OUT"

# Flat parquets from CID_DATASETS ids only (skip legacy HCD flats).
IN_ARGS=()
for entry in "${CID_DATASETS[@]}"; do
  IFS='|' read -r id _url <<< "$entry"
  f="flat/${id}.parquet"
  if [ -f "$f" ]; then
    rows=$(python3 -c "import pyarrow.parquet as pq; print(pq.read_table('$f').num_rows)")
    if [ "$rows" -gt 100 ]; then
      IN_ARGS+=(--in "$f")
      echo "  scoring input $f rows=$rows"
    else
      echo "  WARN skip $f (rows=$rows)"
    fi
  else
    echo "  WARN missing $f"
  fi
done
# Self-labeled Lumos CID flats (from prepare_lumos_cid_training.sh).
for f in flat/PXD007683-Lumos-*.parquet; do
  [ -f "$f" ] || continue
  rows=$(python3 -c "import pyarrow.parquet as pq; print(pq.read_table('$f').num_rows)")
  if [ "$rows" -gt 100 ]; then
    IN_ARGS+=(--in "$f")
    echo "  scoring input $f rows=$rows (Lumos self-label)"
  else
    echo "  WARN skip $f (rows=$rows)"
  fi
done
if [ ${#IN_ARGS[@]} -eq 0 ]; then
  echo "ERROR: no CID flat training files with rows>100"
  exit 1
fi
echo "=== train-from-msnet (${#IN_ARGS[@]} args) $(date -Is) ==="
cp -f "$BUNDLED" "$SCORING_OUT"
"$BIN" train-from-msnet \
  "${IN_ARGS[@]}" \
  --out-store "$SCORING_OUT" \
  --model-id cid_msnet_tryp \
  --seed-model cid_lowres_tryp \
  --fragment-tol-da 0.4 \
  --train-pseudo "$TRAIN_PSEUDO" \
  --train-backoff-weight "$TRAIN_BACKOFF" \
  --train-min-count "$TRAIN_MIN_COUNT" \
  --threads "$THREADS"

echo "=== SANITY $(date -Is) ==="
python3 - <<'PY' "$INTENSITY_OUT" "$SCORING_OUT"
import sys
import pyarrow.parquet as pq
int_path, store_path = sys.argv[1], sys.argv[2]
t = pq.read_table(int_path)
print(f"intensity keys={t.num_rows:,}")
store = pq.read_table(store_path, columns=["model_id"])
print(f"scoring models={store.column('model_id').to_pylist()}")
PY

echo "INTENSITY_MODEL=$INTENSITY_OUT"
echo "SCORING_STORE=$SCORING_OUT"
echo "=== PHASE_V_CID_TRAIN_DONE $(date -Is) ==="
