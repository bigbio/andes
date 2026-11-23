# CLI flags, not environment variables

## Policy

**No environment variable changes a search result.** Everything that affects a search or a
model-training run is a CLI flag with a documented default, so a result can be reproduced
from the command line that produced it. The shipped binary reads no `ANDES_*` variable.

Environment variables are invisible in `--help`, untyped, unvalidated, and easy to tune
around and then forget. This codebase once re-tuned selector weights around a Y-ladder size
bias whose correction sat in the tree behind an unset variable, and shipped a c/z charge
guard that was never wired to a call site.

## Rules

- A switch that changes results is either **deleted** once its A/B is settled (the winning
  behaviour becomes the only behaviour) or becomes **a typed CLI flag**, documented in
  DOCS.md and `--help` in the same commit.
- Values that reach hot inner functions are installed once at startup from validated CLI
  input (`ScoringSettings`, `CzSettings`, `init_y_max_charge`, `init_ethcd_as_etd`,
  `init_dense_noise`), not read from the environment on each call.
- Flags that only affect presentation (`--quiet`, `--verbose`) are still flags, never
  environment variables.
- The only `ANDES_*` variables allowed are **test-harness** variables that select optional
  fixtures at `cargo test` time. They are listed in `TEST_HARNESS_VARS` in
  `crates/search/tests/env_registry.rs`, which fails when code reads a variable missing
  from that list or the list names one the code no longer reads. Add a new test-harness
  variable there with the file that reads it; never add a variable read by shipped code.
- Variables of other tools (`DOTNET_ROOT`, `PATH`) are fine where andes must find that tool.
