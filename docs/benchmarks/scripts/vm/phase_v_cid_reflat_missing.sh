#!/usr/bin/env bash
# Re-build flat parquets that failed strict gates (LysC, PXD016986).
set -euo pipefail

cd /srv/data/msnet
FLAT="${FLAT:-/srv/data/msnet/msnet_to_flat.py}"
PEP_MAX="${PEP_MAX:-0.001}"
SAMPLE="${SAMPLE:-300000}"
export FRAG_FILTER=cid

# dataset_id|url|min_matched|min_explained
REFLAT=(
  "PXD000865-LysC|https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/msnet/PXD000865-LysC/PXD000865-LysC-MSNet.parquet|6|0.50"
  "PXD016986|https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-collections/msnet/PXD016986/PXD016986-MSNet.parquet|6|0.45"
)

for entry in "${REFLAT[@]}"; do
  IFS='|' read -r id url mm me <<< "$entry"
  raw="raw/${id}.reflat.parquet"
  flat_out="flat/${id}.parquet"
  echo "=== REFLAT $id min_matched=$mm min_explained=$me $(date -Is) ==="
  curl -sS -L -o "$raw" "$url"
  rm -f "$flat_out"
  python3 "$FLAT" "$raw" "$flat_out" "$PEP_MAX" "$SAMPLE" 0 "$mm" "$me" 0 cid
  rm -f "$raw"
  rows=$(python3 -c "import pyarrow.parquet as pq; print(pq.read_table('$flat_out').num_rows)")
  echo "[$id] flat rows=$rows"
done
echo "=== REFLAT_DONE $(date -Is) ==="
