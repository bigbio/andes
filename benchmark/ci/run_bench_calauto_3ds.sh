#!/usr/bin/env bash
#
# Three-dataset precursor-cal bench harness (G1 ship gate).
#
# Runs simas with `--precursor-cal auto` against LFQ (PXD001819), Astral,
# and TMT fixtures. Pair with `compare_*_3arm_percolator.sh` to compute
# Percolator @1% FDR vs the Java baseline.
#
# G1 ship gate: Rust @1% FDR within ±1% of Java on all three datasets.
# See `docs/parity-analysis/notes/2026-05-25-precursor-cal-ship-gates.md`
# for the current gate status (G1 NOT yet closed in PR-A).
#
# Override dataset paths via env vars; defaults match the bigbio bench VM
# layout.

set -euo pipefail

LFQ_MGF="${LFQ_MGF:-/srv/data/msgf-bench/PXD001819/sample.mgf}"
LFQ_FASTA="${LFQ_FASTA:-/srv/data/msgf-bench/PXD001819/human-uniprot-contaminants.fasta}"
LFQ_PARAM="${LFQ_PARAM:-HCD_QExactive_Tryp.param}"

ASTRAL_MZML="${ASTRAL_MZML:-/srv/data/msgf-bench/astral/sample.mzML}"
ASTRAL_FASTA="${ASTRAL_FASTA:-/srv/data/msgf-bench/astral/human.fasta}"
ASTRAL_PARAM="${ASTRAL_PARAM:-HCD_HighRes_Tryp.param}"

TMT_MGF="${TMT_MGF:-/srv/data/msgf-bench/tmt/sample.mgf}"
TMT_FASTA="${TMT_FASTA:-/srv/data/msgf-bench/tmt/uniprot.fasta}"
TMT_PARAM="${TMT_PARAM:-HCD_HighRes_Tryp_TMT.param}"

OUT_DIR="${OUT_DIR:-./bench-results/calauto-$(date +%Y%m%d-%H%M)}"
SIMAS="${SIMAS:-./target/release/simas}"
MODE="${MODE:-auto}"

mkdir -p "${OUT_DIR}"

if [ ! -x "${SIMAS}" ]; then
  echo "ERROR: simas binary not found at ${SIMAS}" >&2
  echo "Run: cargo build --release -p simas" >&2
  exit 1
fi

run_one() {
  local label="$1" spectra="$2" fasta="$3" param="$4" extra="${5:-}"
  if [ ! -f "${spectra}" ]; then
    echo "WARN: skipping ${label} (spectra missing: ${spectra})" >&2
    return 0
  fi
  if [ ! -f "${fasta}" ]; then
    echo "WARN: skipping ${label} (fasta missing: ${fasta})" >&2
    return 0
  fi
  echo "=== ${label} (--precursor-cal ${MODE}) ==="
  local pin_path="${OUT_DIR}/${label}.pin"
  local log_path="${OUT_DIR}/${label}.log"
  /usr/bin/time -v "${SIMAS}" \
    --spectrum "${spectra}" \
    --database "${fasta}" \
    --param-file "${param}" \
    --precursor-cal "${MODE}" \
    --output-pin "${pin_path}" \
    ${extra} \
    > "${log_path}" 2>&1
  echo "wrote ${pin_path} (log: ${log_path})"
}

run_one "lfq"    "${LFQ_MGF}"    "${LFQ_FASTA}"    "${LFQ_PARAM}"
run_one "astral" "${ASTRAL_MZML}" "${ASTRAL_FASTA}" "${ASTRAL_PARAM}"
run_one "tmt"    "${TMT_MGF}"    "${TMT_FASTA}"    "${TMT_PARAM}"

echo
echo "Bench complete. PINs in ${OUT_DIR}/"
echo "Next: feed each PIN to Percolator and compare 1% FDR target counts vs Java."
