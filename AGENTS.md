# Andes agent instructions

## Code comments

- Do not add short single-line `//` comments to the code.
- Doc comments (`///`, `//!`) are fine, for example to describe a function's parameters.
- Existing comments stay; do not remove them only to follow these rules.
- Put design decisions and their rationale in commit messages and pull request descriptions, not in code comments.

## Rules by task

Read the page that matches the task before starting it.

| Page | Use when |
|---|---|
| [`docs/dev/bug-prevention.md`](docs/dev/bug-prevention.md) | editing any Rust code, fixing a bug, adding a flag |
| [`docs/dev/clean-logs.md`](docs/dev/clean-logs.md) | adding or changing stderr output, warnings, errors |
| [`docs/dev/cli-reports.md`](docs/dev/cli-reports.md) | changing the run summary, statistics.log or `--help` |
| [`docs/dev/performance.md`](docs/dev/performance.md) | touching a hot path or a structure that scales with the data |
| [`docs/dev/open-data-models.md`](docs/dev/open-data-models.md) | training, porting or shipping any model or weights |
| [`docs/dev/scientific-citations.md`](docs/dev/scientific-citations.md) | writing docs, NOTICE or doc comments that credit prior work |
| [`docs/dev/fdr-validation.md`](docs/dev/fdr-validation.md) | changing scores, features, decoys, rescoring or quantification |
| [`docs/dev/benchmark-reproducibility.md`](docs/dev/benchmark-reproducibility.md) | running or publishing benchmarks |
| [`docs/dev/parameters-and-constants.md`](docs/dev/parameters-and-constants.md) | proposing an algorithm, adding a flag, or a number, threshold, default or mass |
| [`docs/dev/cli-flags-not-env-vars.md`](docs/dev/cli-flags-not-env-vars.md) | adding a setting or switch, or code that reads `std::env` |
| [`docs/dev/mass-spec-correctness.md`](docs/dev/mass-spec-correctness.md) | masses, tolerances, modifications, digestion, decoys, readers |

## Enforcement

These rules are checked by tests, not only by review: `crates/andes/tests/structure.rs`
checks that every page above exists and is listed here, and further structure tests in the
same crate and in `crates/search/tests` fail the build when a rule is broken. See
[`CONTRIBUTING.md`](CONTRIBUTING.md) for how to work with them.
