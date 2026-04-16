# CI benchmark (PXD001819)

- **Workflow:** `.github/workflows/benchmark-pxd001819.yml` (`workflow_dispatch`)
- **Run locally:** `mvn -B package -DskipTests && bash benchmark/ci/PXD001819/run_ci.sh`
- **Compare:** `python3 benchmark/ci/PXD001819/compare_metrics.py benchmark/results/PXD001819/ci/ci_metrics.txt benchmark/ci/PXD001819/baseline.tsv`

GitHub Actions: workflow **Benchmark PXD001819** (`workflow_dispatch`) on `self-hosted,linux,msgf-benchmark`.
