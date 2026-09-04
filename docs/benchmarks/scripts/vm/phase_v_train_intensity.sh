#!/usr/bin/env bash
# Phase T: train intensity_model.parquet from MSNet RAW parquets.
#
# Stages 3 HCD tryptic datasets, aggregates annotated b/y log-relative intensities,
# merges partials, and finalizes via `andes train-intensity`.
#
# Usage:
#   bash /srv/data/msgf-bench/phase_v_train_intensity.sh
#
# Output:
#   /srv/data/msnet/intensity_model.parquet
set -euo pipefail

cd /srv/data/msnet
mkdir -p raw partial

REPO="${REPO:-/srv/data/msgf-bench/repo/msgf-rust}"
BIN="${BIN:-$REPO/target/release/andes}"
AGG="${AGG:-/srv/data/msnet/msnet_intensity_agg.py}"
OUT="${OUT:-/srv/data/msnet/intensity_model.parquet}"
PEP_MAX="${PEP_MAX:-0.001}"
SAMPLE="${SAMPLE:-300000}"
THREADS="${THREADS:-8}"
export THREADS

# dataset_id|url
DATASETS=(
  "PXD000865|https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/msnet/PXD000865/PXD000865-MSNet.parquet"
  "IPX0002031000|https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/msnet/IPX0002031000/IPX0002031000-MSNet.parquet"
  "PXD024364-Trypsin|https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/msnet/PXD024364-Trypsin/PXD024364-Trypsin-MSNet.parquet"
)

echo "=== BUILD andes $(date -Is) ==="
( cd "$REPO" && cargo +1.95.0 build --release -p andes 2>&1 | tail -5 )

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

stage_one() {
  local idx="$1" total="$2" id="$3" url="$4"
  local raw="raw/${id}.parquet"
  local partial="partial/${id}.parquet"
  echo "[dataset ${idx}/${total}] [$id] $(date -Is)"
  if [ -f "$partial" ]; then
    echo "[$id] partial exists, skip"
    return 0
  fi
  download_with_progress "$id" "$url" "$raw"
  echo "[$id] aggregate pep_max=$PEP_MAX sample=$SAMPLE threads=$THREADS"
  python3 "$AGG" "$raw" "$partial" "$PEP_MAX" "$SAMPLE"
  rm -f "$raw"
  echo "[$id] partial keys=$(python3 -c "import pyarrow.parquet as pq; print(pq.read_table('$partial').num_rows)")"
  echo "[$id] free=$(df -h /srv/data | tail -1 | awk '{print $4}')"
}

TOTAL=${#DATASETS[@]}
IDX=0
for entry in "${DATASETS[@]}"; do
  IDX=$((IDX + 1))
  IFS='|' read -r id url <<< "$entry"
  stage_one "$IDX" "$TOTAL" "$id" "$url"
done

PARTIALS=(partial/*.parquet)
echo "=== MERGE ${#PARTIALS[@]} partials $(date -Is) ==="
python3 "$AGG" --merge "${PARTIALS[@]}" -o "$OUT.partial_merged.parquet"

echo "=== train-intensity $(date -Is) ==="
"$BIN" train-intensity \
  --in "$OUT.partial_merged.parquet" \
  --out "$OUT"

echo "=== SANITY $(date -Is) ==="
python3 - <<'PY' "$OUT"
import sys
import pyarrow.parquet as pq
t = pq.read_table(sys.argv[1])
print(f"keys={t.num_rows:,} cols={t.column_names}")
# y after K/R should dominate b1 on average (train-intensity also checks this).
rows = t.to_pylist()
from statistics import mean
y_kr = [r["mean_log_rel"] for r in rows if r["ion_type"] == "y" and r["flank_n"] in "KR" and r["count"] >= 100]
b1 = [r["mean_log_rel"] for r in rows if r["ion_type"] == "b" and r["pos_bin"] <= 1 and r["count"] >= 100]
if y_kr and b1:
    print(f"y(K|R) mean={mean(y_kr):.3f}  b1 mean={mean(b1):.3f}")
PY

echo "INTENSITY_MODEL=$OUT"
echo "=== PHASE_T_DONE $(date -Is) ==="
