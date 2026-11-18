# Andes

@AGENTS.md

## Skills

Task-specific guidance lives in `.claude/skills/` and loads when a task needs it:

| Skill | Use when |
|---|---|
| `bug-prevention` | editing any Rust code, fixing a bug, adding a flag |
| `clean-logs` | adding or changing stderr output, warnings, errors |
| `cli-reports` | changing the run summary, statistics.log or `--help` |
| `performance` | touching a hot path or a structure that scales with the data |
| `open-data-models` | training, porting or shipping any model or weights |
| `scientific-citations` | writing docs, NOTICE or doc comments that credit prior work |
| `fdr-validation` | changing scores, features, decoys, rescoring or quantification |
| `benchmark-reproducibility` | running or publishing benchmarks |
| `cli-flags-not-env-vars` | adding a setting or switch, or code that reads `std::env` |
| `mass-spec-correctness` | masses, tolerances, modifications, digestion, decoys, readers |
