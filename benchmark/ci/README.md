# CI benchmark (PXD001819)

- **Plan:** [`.claude/plans/pxd001819-ci-benchmark.md`](../../.claude/plans/pxd001819-ci-benchmark.md)
- **Run locally:** `mvn -B package -DskipTests && bash benchmark/ci/PXD001819/run_ci.sh`
- **Compare:** `python3 benchmark/ci/PXD001819/compare_metrics.py benchmark/results/PXD001819/ci/ci_metrics.txt benchmark/ci/PXD001819/baseline.tsv`

GitHub Actions: workflow **Benchmark PXD001819** (`workflow_dispatch`) on `self-hosted,linux,msgf-benchmark`.
