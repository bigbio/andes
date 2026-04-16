#!/usr/bin/env bash
# CI/local: download PXD001819 public files, run MS-GF+, write dataset-scoped ci_metrics.txt
# Prereq: GNU time (/usr/bin/time -v), curl, gzip; JAR at target/MSGFPlus.jar (override with JAR=).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
DATA_DIR="${DATA_DIR:-$REPO_ROOT/benchmark/data/PXD001819}"
OUT_DIR="${OUT_DIR:-$REPO_ROOT/benchmark/results/PXD001819/ci}"
JAR="${JAR:-$REPO_ROOT/target/MSGFPlus.jar}"
MODS="${MODS:-$REPO_ROOT/src/test/resources/benchmark/PXD001819/mods.txt}"
THREAD_COUNT="${MSGFPLUS_THREADS:-4}"
JVM_MEM="${MSGFPLUS_MEMORY:-4096m}"

PRIDE_MZML_GZ="https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks/lfq/LTQOrbitrapVelos/PXD001819/UPS1_5000amol_R1.mzML.gz"
FASTA_URL="https://raw.githubusercontent.com/bigbio/quantms-test-datasets/quantms/databases/PXD001819_uniprot_yeast_ups.fasta"

MZML_GZ="$DATA_DIR/UPS1_5000amol_R1.mzML.gz"
MZML="$DATA_DIR/UPS1_5000amol_R1.mzML"
FASTA="$DATA_DIR/PXD001819_uniprot_yeast_ups.fasta"
MZID="$OUT_DIR/ci_output.mzid"
TIME_TXT="$OUT_DIR/gnu_time.txt"
METRICS="$OUT_DIR/ci_metrics.txt"

SEARCH_ARGS="-tda 1 -t 5ppm -ti 0,1 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1 -msLevel 2 -thread $THREAD_COUNT"

mkdir -p "$DATA_DIR" "$OUT_DIR"

if [[ ! -f "$JAR" ]]; then
  echo "ERROR: JAR not found: $JAR (run mvn package first)" >&2
  exit 1
fi
if [[ ! -f "$MODS" ]]; then
  echo "ERROR: Mods file not found: $MODS" >&2
  exit 1
fi

download_if_missing() {
  local url="$1" dest="$2"
  if [[ -f "$dest" ]]; then
    echo "OK (cached): $dest"
    return 0
  fi
  echo "Downloading $(basename "$dest") ..."
  curl -fL --retry 3 --connect-timeout 30 -o "$dest" "$url"
}

download_if_missing "$PRIDE_MZML_GZ" "$MZML_GZ"
if [[ ! -f "$MZML" ]]; then
  echo "Decompressing mzML ..."
  gunzip -c "$MZML_GZ" >"$MZML"
fi
download_if_missing "$FASTA_URL" "$FASTA"

rm -f "$DATA_DIR"/*.canno "$DATA_DIR"/*.cnlcp "$DATA_DIR"/*.csarr "$DATA_DIR"/*.cseq 2>/dev/null || true

echo "Running MS-GF+ (wall clock via \$SECONDS) ..."
START_SECONDS=$SECONDS
set +e
/usr/bin/time -v -o "$TIME_TXT" \
  java "-Xmx${JVM_MEM}" -jar "$JAR" \
    -s "$MZML" \
    -d "$FASTA" \
    -mod "$MODS" \
    -o "$MZID" \
    $SEARCH_ARGS \
    >"$OUT_DIR/run.stdout.log" 2>"$OUT_DIR/run.stderr.log"
JAVA_RC=$?
set -e
WALL=$((SECONDS - START_SECONDS))

if [[ ! -f "$MZID" ]]; then
  echo "ERROR: mzIdentML not created (java exit $JAVA_RC)" >&2
  {
    echo "dataset=PXD001819"
    echo "error=missing_mzid"
    echo "java_exit=$JAVA_RC"
    echo "wall_time_sec=$WALL"
  } >"$METRICS"
  exit 1
fi

if [[ "$JAVA_RC" -ne 0 ]]; then
  echo "ERROR: MS-GF+ exited with code $JAVA_RC" >&2
  exit "$JAVA_RC"
fi

python3 - "$TIME_TXT" "$MZID" "$METRICS" "$WALL" "$JAVA_RC" <<'PY'
import re, sys
from pathlib import Path

time_path, mzid_path, metrics_path, wall, java_rc = sys.argv[1:6]
wall_sec = float(wall)
text = Path(time_path).read_text(errors="replace")
m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
rss_kb = m.group(1) if m else ""
m2 = re.search(r"Percent of CPU this job got: (\d+)", text)
cpu_pct = m2.group(1) if m2 else ""

mzid = Path(mzid_path)
sii = 0
for line in mzid.open(errors="replace"):
    if "SpectrumIdentificationItem" in line:
        sii += 1

psm_1pct = 0
with mzid.open(errors="replace") as f:
    for line in f:
        m = re.search(r'accession="MS:1002054".*?value="([^"]+)"', line)
        if m:
            try:
                if float(m.group(1)) <= 0.01:
                    psm_1pct += 1
            except ValueError:
                pass

lines = [
    "dataset=PXD001819",
    f"wall_time_sec={wall_sec}",
    f"java_exit={java_rc}",
    f"sii_count={sii}",
    f"psm_1pct_fdr={psm_1pct}",
]
if rss_kb:
    lines.append(f"peak_rss_kb={rss_kb}")
if cpu_pct:
    lines.append(f"cpu_percent={cpu_pct}")
Path(metrics_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
print(Path(metrics_path).read_text())
PY

echo "Wrote $METRICS"
