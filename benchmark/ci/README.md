# CI benchmark (PXD001819)

- **Workflow:** `.github/workflows/benchmark-pxd001819.yml` (`workflow_dispatch`)
- **Run locally:** `mvn -B package -DskipTests && bash benchmark/ci/PXD001819/run_ci.sh`
- **Compare:** `python3 benchmark/ci/PXD001819/compare_metrics.py benchmark/results/PXD001819/ci/ci_metrics.txt benchmark/ci/PXD001819/baseline.tsv`
- **Self-test:** `python3 -m unittest benchmark.ci.PXD001819.test_compare_metrics`

GitHub Actions: workflow **Benchmark PXD001819** (`workflow_dispatch`) on `self-hosted,linux,msgf-benchmark`. Python 3.11 is pinned via `actions/setup-python`.

## Scripts

| Script | Purpose |
|--------|---------|
| `run_ci.sh` | Downloads public inputs, runs MS-GF+, invokes `extract_metrics.py` |
| `extract_metrics.py` | Streams the mzIdentML (ElementTree `iterparse`) to count SII and PSMs at 1% FDR; also extracts RSS/CPU from `time -v` |
| `compare_metrics.py` | Compares key=value metrics to the baseline TSV |
| `test_compare_metrics.py` | Unit tests for the comparator |
