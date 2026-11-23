# Clean logs

andes writes its log to stderr and its results to files. A user should be able to read a
whole default run log in one screen and know what happened.

## System

All stderr output goes through the `log` facade. The binary installs one logger,
`crates/andes/src/bin/andes/logging.rs`; library crates call `log::info!`/`warn!`/`debug!`
and never print. `eprintln!`/`eprint!` are not allowed outside tests and examples.

- A disabled call costs one relaxed atomic load and formats nothing, so `log::debug!` is
  safe near hot code. Never compute an expensive value just to log it; guard with
  `log::log_enabled!(log::Level::Debug)` when you must.
- Each line is formatted once and written in a single `write_all` on locked stderr, so
  threads never interleave lines.
- `--quiet` shows warnings and errors only, `--verbose` adds debug lines.

## Format

One style for every line, produced by the logger:

```
[   12.3s] INFO  index: done in 3.42 s
[   12.9s] WARN  calibration: only 41 confident PSMs in run1.mzML; ...
[   13.0s] ERROR andes: cannot open run2.mzML: No such file or directory
```

- Time since the run started, level, stage, message. Never put a level word, a stage
  prefix, `[PHASE ...]` or `error:` in the message text itself.
- The stage is the call's `target:`, or the last segment of the module path. Set
  `target:` when the line belongs to another stage than its module: `input`, `database`,
  `index`, `search`, `calibration`, `model`, `rescore`, `quant`, `mbr`, `train`, `output`,
  `run` (start-of-run banner), `summary` (end-of-run report).
- Stage timings: `log::info!(target: "index", "done in {:.2} s", secs)`.
- Numbers: thousands separators for counts (`52,169 PSMs`), units always (`s`, `GB`, `ppm`,
  `Da`), at most 3 significant decimals.
- Name flags as the user types them (`--precursor-cal`), files by path, spectra by
  `file:scan`.

## Volume

- Default level: stage start/finish, what was decided automatically (model picked,
  instrument detected, tolerance calibrated, memory budget), warnings, the final summary.
- Never log inside a per-spectrum, per-peptide or per-fragment loop. Count instead and
  print one aggregated line at the end of the stage:
  `WARN  input: 312 spectra skipped (no peaks: 290, charge 0: 22)`.
- A warning is printed once per run, not once per file or per thread. Deduplicate with a
  counter or a `std::sync::Once`.
- Diagnostics for developers are `log::debug!`, shown only with `--verbose`.
- No progress bars or carriage-return animation when stderr is not a terminal
  (`std::io::IsTerminal`); CI and cluster logs must stay line-based.

## Warnings and errors

A warning says what happened, its consequence, and what to do:

```
log::warn!(target: "calibration", "only {n} confident PSMs in {file}; precursor calibration skipped (needs 100). The search used the uncalibrated tolerance");
```

An error names the input that caused it and exits non-zero. Never print a Rust `Debug`
dump, a panic backtrace or `called Option::unwrap() on a None value` to a user.

## Check

Run a small search (test-fixtures/) and read the full stderr. Every line should be
necessary; if two lines say the same thing, merge them.
