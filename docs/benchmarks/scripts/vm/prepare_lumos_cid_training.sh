#!/usr/bin/env bash
# Build Lumos TMT CID flat training parquets from MSFragger labels + mzML peaks.
# No MSNet parquet exists for PXD007683; this is self-labeled rank training data.
set -euo pipefail

OUT_DIR="${OUT_DIR:-/srv/data/msnet/flat}"
WORK="${WORK:-/srv/data/msgf-bench/lumos-train-work}"
BASE="${BASE:-/srv/data/msgf-bench/tmt-data}"
FRAGGER="${FRAGGER:-/srv/data/msgf-bench/engines/msfragger-env/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar}"
TRFP="${TRFP:-/srv/data/.conda/envs/nextflow/bin/ThermoRawFileParser}"
FLAT_PY="${FLAT_PY:-/srv/data/msgf-bench/mzml_pepxml_to_flat.py}"
URLBASE="${URLBASE:-https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD007683}"
# a05058 has cached mzML; others download raw → mzML
FILES="${FILES:-a05058 a05059 a05060 a05061}"
EXPECT_MAX="${EXPECT_MAX:-0.01}"
MAX_PSMS="${MAX_PSMS:-40000}"

mkdir -p "$OUT_DIR" "$WORK"

run_one() {
  local tag="$1"
  local mzml="$WORK/${tag}.mzML"
  local pepxml="$WORK/${tag}.pepXML"
  local flat="$OUT_DIR/PXD007683-Lumos-${tag}.parquet"

  echo "=== Lumos prep $tag $(date -Is) ==="
  if [ "$tag" = "a05058" ] && [ -f "$BASE/a05058.mzML" ]; then
    cp -f "$BASE/a05058.mzML" "$mzml"
  elif [ ! -f "$mzml" ]; then
    local raw="$WORK/${tag}.raw"
    echo "  download $tag.raw"
    curl -sS -L -o "$raw" "$URLBASE/${tag}.raw"
    "$TRFP" -i "$raw" -b "$mzml" -f 2 -l 2
    rm -f "$raw"
  fi

  if [ ! -f "$pepxml" ]; then
    echo "  [msfragger] $tag → pepXML"
    local params="$WORK/fragger_${tag}.params"
    cp -f "$BASE/fragger_a05058.params" "$params"
    sed -i 's/^output_format = .*/output_format = pepXML/' "$params"
    java -Xmx8g -jar "$FRAGGER" "$params" "$mzml" > "$WORK/${tag}.fragger.log" 2>&1
    [ -f "$WORK/${tag}.pepXML" ] && mv -f "$WORK/${tag}.pepXML" "$pepxml"
    [ -f "${mzml%.mzML}.pepXML" ] && mv -f "${mzml%.mzML}.pepXML" "$pepxml"
  fi
  [ -f "$pepxml" ] || { echo "  ERROR: no pepXML for $tag"; return 1; }

  python3 "$FLAT_PY" "$mzml" "$pepxml" "$flat" "$EXPECT_MAX" "$MAX_PSMS"
  rows=$(python3 -c "import pyarrow.parquet as pq; print(pq.read_table('$flat').num_rows)")
  echo "  flat $flat rows=$rows"
  rm -f "$pepxml" "$WORK/${tag}.pin" "$WORK/${tag}.tsv"
}

for tag in $FILES; do
  run_one "$tag" || true
done

echo "=== LUMOS_CID_PREP_DONE $(date -Is) ==="
ls -lh "$OUT_DIR"/PXD007683-Lumos-*.parquet 2>/dev/null || true
