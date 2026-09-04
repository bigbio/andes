#!/bin/bash
# Renders a single Markdown report from the per-dataset metrics.tsv + psm_summary.tsv files.
# Usage: bash render_validation_report.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$SCRIPT_DIR/results/validation"
REPORT="$ROOT/REPORT.md"

mkdir -p "$ROOT"

{
    echo "# Master vs branch validation"
    echo
    echo "Branch: feature/improve-mzid-suffix-big-fasta"
    echo
    echo "Generated: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo
    for ds in pxd001819 tmt astral; do
        local_summary="$ROOT/$ds/psm_summary.tsv"
        local_metrics="$ROOT/$ds/metrics.tsv"
        [ -f "$local_summary" ] || continue
        echo "## $ds"
        echo
        echo "### PSM / peptide counts"
        echo
        echo '```'
        column -t -s $'\t' "$local_summary"
        echo '```'
        echo
        echo "### Per-phase metrics (wall_s, max_rss_mb, user_s, sys_s, cpu_pct)"
        echo
        echo '```'
        column -t -s $'\t' "$local_metrics"
        echo '```'
        echo
        # Speedup / memory deltas
        python3 - "$local_metrics" <<'PY' || true
import csv, sys, collections
path = sys.argv[1]
rows = list(csv.DictReader(open(path), delimiter='\t'))
by = collections.defaultdict(dict)
for r in rows:
    by[r['phase']][r['arm']] = r
print('| phase | wall master | wall branch | speedup | rss master | rss branch | rss delta |')
print('|---|---|---|---|---|---|---|')
for phase in ['build', 'search', 'percolator']:
    m, b = by.get(phase, {}).get('master'), by.get(phase, {}).get('branch')
    if not m or not b: continue
    try:
        wm, wb = float(m['wall_s']), float(b['wall_s'])
        rm, rb = float(m['max_rss_mb']), float(b['max_rss_mb'])
        speedup = (wm/wb) if wb>0 else float('nan')
        rss_delta = (rb-rm)
        print(f"| {phase} | {wm:.1f}s | {wb:.1f}s | {speedup:.2f}x | {rm:.0f} MB | {rb:.0f} MB | {rss_delta:+.0f} MB |")
    except Exception as e:
        print(f"| {phase} | err: {e} | | | | | |")
PY
        echo
    done
} > "$REPORT"

echo "wrote $REPORT"
cat "$REPORT"
