# CI benchmark (PXD001819)

- **Workflow:** `.github/workflows/benchmark-pxd001819.yml` (`workflow_dispatch`)
- **Run locally:** `mvn -B package -DskipTests && bash benchmark/ci/PXD001819/run_ci.sh`
- **Compare:** `python3 benchmark/ci/PXD001819/compare_metrics.py benchmark/results/PXD001819/ci/ci_metrics.txt benchmark/ci/PXD001819/baseline.tsv`
- **Self-test:** `python3 -m unittest benchmark.ci.PXD001819.test_compare_metrics`

GitHub Actions: workflow **Benchmark PXD001819** (`workflow_dispatch`) on `self-hosted,linux,msgf-benchmark`. Python 3.11 is pinned via `actions/setup-python`.

## What gets measured

Each CI run produces a `ci_metrics.txt` with key=value pairs:

| Metric | Source | Notes |
|---|---|---|
| `wall_time_sec` | `$SECONDS` around the `java -jar` invocation | End-to-end wall-time, integer seconds |
| `peak_rss_kb` | `/usr/bin/time -v` (Linux) | Optional: not all runners expose this |
| `cpu_percent` | `/usr/bin/time -v` | Optional: parsed but not gated yet |
| `native_target_count` | Count of `Label==1` rows in the `.pin` | Deterministic across runs given same inputs |
| `native_decoy_count` | Count of `Label==-1` rows in the `.pin` | Deterministic across runs |

`baseline.tsv` declares acceptable `[min, max]` ranges per metric. `compare_metrics.py` exits non-zero if any required metric is outside its range.

**The CI gate is search-correctness, not 1 % FDR sensitivity.** Native target/decoy counts are deterministic — same inputs → identical numbers across runs — so they cleanly catch search-code regressions. For 1 % FDR PSM counts you need Percolator on the `.pin`; that's stochastic (seed 42 stabilises it) and is a separate downstream gate, not in this CI.

**Why PIN, not mzIdentML.** PR #23 removed mzIdentML reader/writer entirely; `.pin` is the only modern output format. The CI script outputs `.pin` and parses it directly (one stream-pass, two integer counts) — no XML, no Percolator, no flakiness.

## Scripts

| Script | Purpose |
|--------|---------|
| `run_ci.sh` | Downloads public inputs (mzML.gz from PRIDE, FASTA from `quantms-test-datasets`), runs MS-GF+ with fixed search args, invokes `extract_metrics.py` |
| `extract_metrics.py` | Counts target / decoy rows from the `.pin` (streaming, line-at-a-time); pulls RSS / CPU% from `/usr/bin/time -v` output |
| `compare_metrics.py` | Compares key=value metrics to the baseline TSV; required metrics out of range → exit 1; optional metrics missing → warning |
| `test_compare_metrics.py` | Unit tests: 7 for the comparator, 3 for `parse_pin`. Run with `python3 -m unittest benchmark.ci.PXD001819.test_compare_metrics` |

## Tightening the baseline after a green run

The current `baseline.tsv` ranges are intentionally wide (e.g. wall 60–900 s) to land a first green workflow on whatever runner you provision. After 3–5 successful runs with consistent numbers, narrow each `[min, max]` to roughly ±10 % of the observed median. This is what lets the CI catch real regressions.

## Future iterations

The retrospective `.claude/plans/astral-phase-a-retrospective.md` documents that this single-run CI scaffold is *insufficient* for measuring per-spectrum or thread-pool optimizations on Astral, where wall-time variance from machine state is ~30 %. For those iterations, future agents should:

1. Build a multi-run wrapper that runs N≥5 measurements back-to-back, reports median + IQR, and only flags a regression if the new median is outside the historical IQR.
2. Add CI scaffolds for TMT (PXD007683) and Astral (ProteoBench Module 8) following the same shape.
3. Use a reserved runner with thermal headroom; benchmark output is meaningless on a machine that's been running benchmarks for hours.
