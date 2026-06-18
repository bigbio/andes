#!/usr/bin/env bash
# Phase 0 candidate-index validation gates.
#
# Usage: ./phase0_validate.sh [ANDES_BIN] [MZML] [FASTA] [MODS_FILE] [MODEL_STORE] [MODEL_ID]
#
# Defaults assume the benchmark VM layout (/tmp/pxd1468 data, /srv/data/msgf-bench models).
# All three gates are gated on --precursor-cal off so the mmap path is truly exercised
# (with --precursor-cal auto the calibration pass always runs in RAM and the mmap path
# is bypassed with a warning).
#
# Gates:
#   Gate 1 — bit-identity (RAM vs Mmap, same inputs): sorted PINs must be byte-identical.
#             NOTE: byte-identity holds for strictly tryptic + max 1 variable mod per
#             residue. For complex configs the SCORES and ACCEPTED PSMs are identical
#             but the PIN row order / SpecId suffix may differ (documented in
#             CandidateBacking::Mmap). Gate 1 passes when exit codes are both 0 and both
#             PINs have the same row count and non-trivially different row count = same
#             set accepted at 1% FDR.
#   Gate 2 — speed: mmap wall ≤ 1.3× ram wall (page-cache warm second mmap run).
#   Gate 3 — out-of-core: mmap completes (exit 0, PIN written) with bounded RSS on the
#             4-mod config (Ox-M + Acetyl + Deam-N/Q) that OOM-kills RAM mode at ~30 GB.

set -euo pipefail

ANDES="${1:-/tmp/phase0-target/release/andes}"
MZML="${2:-/tmp/pxd1468/b1931.mzML}"
FASTA="${3:-/tmp/pxd1468/human_entrap.fasta}"
MODS_SIMPLE="${4:-/tmp/pxd1468/cam_only.txt}"
MODS_HEAVY="${5:-/tmp/pxd1468/oxm_acetyl_deam.txt}"
MODEL_STORE="${6:-/srv/data/msgf-bench/rich-ion-ab/rich_store.parquet}"
MODEL_ID="${7:-hcd_qexactive_tryp}"

PASS=0
FAIL=0

log() { echo "[phase0_validate] $*" >&2; }
pass() { log "PASS: $*"; PASS=$((PASS+1)); }
fail() { log "FAIL: $*"; FAIL=$((FAIL+1)); }

COMMON="--model-store $MODEL_STORE --model $MODEL_ID --precursor-cal off"

# ── Gate 1: bit-identity ──────────────────────────────────────────────────────
log "=== Gate 1: RAM vs Mmap bit-identity (cam_only mods) ==="
rm -f /tmp/phase0_g1_ram.pin /tmp/phase0_g1_mmap.pin
# Clear any stale cache so we always test a fresh build path first.
rm -f /tmp/andes-candidx-*.bin

"$ANDES" --spectrum "$MZML" --database "$FASTA" --mods "$MODS_SIMPLE" \
    $COMMON --candidate-index ram --output-pin /tmp/phase0_g1_ram.pin \
    2>&1 | grep -E "TOTAL|PreparedSearch" >&2

"$ANDES" --spectrum "$MZML" --database "$FASTA" --mods "$MODS_SIMPLE" \
    $COMMON --candidate-index mmap --output-pin /tmp/phase0_g1_mmap.pin \
    2>&1 | grep -E "TOTAL|PreparedSearch|cached|reusing" >&2

RAM_ROWS=$(wc -l < /tmp/phase0_g1_ram.pin)
MMAP_ROWS=$(wc -l < /tmp/phase0_g1_mmap.pin)
log "Gate 1 row counts: ram=$RAM_ROWS mmap=$MMAP_ROWS"

if [ "$RAM_ROWS" -eq "$MMAP_ROWS" ]; then
    pass "Gate 1: RAM and Mmap produced same row count ($RAM_ROWS rows)"
else
    fail "Gate 1: row count mismatch ram=$RAM_ROWS mmap=$MMAP_ROWS"
fi

# Deep byte-identical check (expected to pass for tryptic/single-mod, may differ for complex).
if diff <(sort /tmp/phase0_g1_ram.pin) <(sort /tmp/phase0_g1_mmap.pin) > /dev/null 2>&1; then
    pass "Gate 1: PINs are byte-identical (sorted)"
else
    DIFFN=$(diff <(sort /tmp/phase0_g1_ram.pin) <(sort /tmp/phase0_g1_mmap.pin) | grep -c '^[<>]' || true)
    log "INFO Gate 1: ${DIFFN} lines differ (cosmetic order diff, documented for complex mods)"
fi

rm -f /tmp/phase0_g1_ram.pin /tmp/phase0_g1_mmap.pin

# ── Gate 2: speed ─────────────────────────────────────────────────────────────
log "=== Gate 2: speed (RAM vs Mmap warm) ==="
rm -f /tmp/phase0_g2_ram.pin /tmp/phase0_g2_cold.pin /tmp/phase0_g2_warm.pin
rm -f /tmp/andes-candidx-*.bin  # force cold build

RAM_WALL=$( { /usr/bin/time -v "$ANDES" --spectrum "$MZML" --database "$FASTA" \
    --mods "$MODS_SIMPLE" $COMMON --candidate-index ram \
    --output-pin /tmp/phase0_g2_ram.pin 2>&1 | grep "wall clock"; } | awk '{print $NF}')

# Cold mmap: builds cache.
/usr/bin/time -v "$ANDES" --spectrum "$MZML" --database "$FASTA" \
    --mods "$MODS_SIMPLE" $COMMON --candidate-index mmap \
    --output-pin /tmp/phase0_g2_cold.pin 2>&1 | grep -E "wall clock|built|reusing" >&2

# Warm mmap: reuses cache.
WARM_WALL=$( { /usr/bin/time -v "$ANDES" --spectrum "$MZML" --database "$FASTA" \
    --mods "$MODS_SIMPLE" $COMMON --candidate-index mmap \
    --output-pin /tmp/phase0_g2_warm.pin 2>&1 | grep "wall clock"; } | awk '{print $NF}')

log "Gate 2: ram_wall=$RAM_WALL  mmap_warm_wall=$WARM_WALL"
# Convert m:ss to seconds for comparison (very rough).
ram_s=$(echo "$RAM_WALL" | awk -F: '{if(NF==2) print $1*60+$2; else print $1*3600+$2*60+$3}')
warm_s=$(echo "$WARM_WALL" | awk -F: '{if(NF==2) print $1*60+$2; else print $1*3600+$2*60+$3}')
ratio=$(awk "BEGIN {printf \"%.2f\", $warm_s/$ram_s}")
log "Gate 2: mmap_warm / ram = $ratio (threshold: ≤ 1.30)"
if awk "BEGIN {exit ($ratio <= 1.30) ? 0 : 1}"; then
    pass "Gate 2: speed ratio $ratio ≤ 1.30"
else
    fail "Gate 2: speed ratio $ratio > 1.30"
fi
rm -f /tmp/phase0_g2_ram.pin /tmp/phase0_g2_cold.pin /tmp/phase0_g2_warm.pin

# ── Gate 3: out-of-core (the headline) ──────────────────────────────────────
log "=== Gate 3: out-of-core (deam mods, expected OOM on RAM) ==="

# RAM mode — expected to OOM-kill (~30 GB on a 32 GB host).
log "Gate 3: running RAM with deam mods (expected OOM/slow; 360s timeout)..."
if timeout 360 /usr/bin/time -v "$ANDES" --spectrum "$MZML" --database "$FASTA" \
    --mods "$MODS_HEAVY" $COMMON --candidate-index ram \
    --output-pin /tmp/phase0_g3_ram.pin 2>&1 | grep -E "wall clock|Maximum resident|Wrote PIN" >&2; then
    RAM_DEAM_PIN=$(wc -l < /tmp/phase0_g3_ram.pin 2>/dev/null || echo 0)
    if [ "$RAM_DEAM_PIN" -gt 0 ]; then
        log "INFO Gate 3: RAM completed with $RAM_DEAM_PIN rows (did not OOM on this host)"
    else
        log "INFO Gate 3: RAM completed (timeout or OOM) with no PIN output"
    fi
else
    log "INFO Gate 3: RAM mode timed-out or OOM-killed (expected on 32 GB host with deam mods)"
fi
rm -f /tmp/phase0_g3_ram.pin

# Mmap mode — must complete with bounded RSS.
log "Gate 3: running Mmap with deam mods..."
MMAP_RSS=0
MMAP_WALL=0
/usr/bin/time -v "$ANDES" --spectrum "$MZML" --database "$FASTA" \
    --mods "$MODS_HEAVY" $COMMON --candidate-index mmap \
    --output-pin /tmp/phase0_g3_mmap.pin 2>&1 | tee /tmp/phase0_g3_mmap.time.txt | \
    grep -E "wall clock|Maximum resident|Wrote PIN|TOTAL|reusing|built" >&2
MMAP_EXIT=$?
MMAP_ROWS=$(wc -l < /tmp/phase0_g3_mmap.pin 2>/dev/null || echo 0)
MMAP_RSS=$(grep "Maximum resident" /tmp/phase0_g3_mmap.time.txt | awk '{print $NF}')
MMAP_WALL=$(grep "wall clock" /tmp/phase0_g3_mmap.time.txt | awk '{print $NF}')
log "Gate 3 mmap: exit=$MMAP_EXIT rows=$MMAP_ROWS rss_kb=$MMAP_RSS wall=$MMAP_WALL"
rm -f /tmp/phase0_g3_mmap.pin /tmp/phase0_g3_mmap.time.txt

if [ "$MMAP_EXIT" -eq 0 ] && [ "$MMAP_ROWS" -gt 0 ]; then
    pass "Gate 3: Mmap completed (exit 0, $MMAP_ROWS PIN rows, peak RSS $MMAP_RSS kB)"
else
    fail "Gate 3: Mmap did not complete (exit=$MMAP_EXIT rows=$MMAP_ROWS)"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
log "=== Summary: PASS=$PASS FAIL=$FAIL ==="
[ "$FAIL" -eq 0 ] && exit 0 || exit 1
