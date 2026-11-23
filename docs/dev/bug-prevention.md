# Bug prevention

## Before writing code

1. Find the existing behaviour first: grep for the function, its callers and its tests.
   Most andes bugs came from a second code path (DIA, `--refine`, `--chimeric`, out-of-core,
   per-file model selection) that the change did not reach.
2. For a bug fix, write the failing test first and watch it fail. The test name states the
   rule in a sentence (`a_terminal_site_is_reported_only_where_the_peptide_is_terminal`).
3. List every mode the change touches: DDA, `--dia`, `--chimeric`, `--refine`, `--glyco`,
   `--open`, xl, `--lfq/--tmt`, `train-from-search`. Check each one, or say in the PR which
   were not checked.

## While writing code

- No `unwrap()`/`expect()` on data from a file, the CLI or the network. Return an error
  that names the file, the scan or the flag. `unwrap()` is acceptable on invariants the
  code itself establishes, with a `///` doc on why it holds where it is not obvious. A new
  `unwrap`/`expect` in `crates/*/src` still fails the unwrap ratchet in
  `docs/dev/structure-checks.md`; a documented invariant is not exempt yet.
- No silent fallbacks. If input is malformed, either reject it with a clear message or
  count it and report the count once (see `docs/dev/clean-logs.md`). Never drop spectra,
  PSMs or proteins without a counter.
- Floats: compare with `total_cmp` for sorting, never `partial_cmp().unwrap()`. Watch for
  NaN from empty spectra, zero intensities and `ln(0)`.
- Integer widths: scan numbers, charges and indices from files are untrusted. Use checked
  or saturating arithmetic at the boundary.
- Every new result-changing switch is a typed CLI flag with a documented default
  (`docs/dev/cli-flags-not-env-vars.md`). Never read an `ANDES_*` environment variable in shipped code; only locating external tools (`PATH`, `DOTNET_ROOT`) may use the environment.
- Determinism: results must not depend on thread count, `HashMap` iteration order or
  rayon scheduling. Sort before writing; use `rustc_hash`/`BTreeMap` where order leaks out.

## Edge cases every input path must survive

Empty file, one spectrum, spectrum with zero peaks, charge 0 or missing, precursor m/z 0,
duplicate scan numbers, gzip input, CRLF line endings, very long FASTA headers, `X`/`B`/`Z`/`U`
residues, protein shorter than the minimum peptide length. `crates/andes/tests/empty_input.rs`
and `input_guards.rs` are the templates.

## Before committing

```bash
cargo clippy --locked --workspace --all-targets -- -D warnings
```

```bash
cargo test --release --workspace
```

If the change can touch the PIN, also run the bit-identity goldens:

```bash
cargo test -p andes --release --test pin_goldens -- --ignored
```

A change that is meant to be result-neutral (refactor, speed-up) must keep the PIN
bit-identical. A change that is meant to alter results must state the new PSM counts and
the entrapment FDP in the commit message (see `fdr-validation`).

## Review your own diff

Read `git diff` once as a reviewer before committing: dead code, leftover debug prints,
new `//` comments (not allowed, see AGENTS.md), flags missing from `--help` or DOCS.md.
