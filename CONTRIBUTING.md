# Contributing to andes

## Build and test

```bash
cargo build --release --workspace
cargo test --release --workspace
cargo clippy --workspace --all-targets -- -D warnings
```

The PIN goldens are ignored by default because they are slow; CI runs them with

```bash
cargo test -p andes --release --test out_of_core_pin_goldens -- --ignored
```

Format only the files you change (`rustfmt <file>`); never run `cargo fmt --all`, because
`main` is not fully formatted.

## Rules

The rules for contributors, human or agent, are in [`AGENTS.md`](AGENTS.md), which lists one
page per kind of task under [`docs/dev/`](docs/dev/). Read the page that matches your change
before starting. Tests in `crates/andes/tests/structure.rs` and the other structure tests
keep these pages and the code in step.

## Ratchets

Some structure tests hold a baseline: a list or a count of existing violations that may
shrink but never grow. A test fails when code adds a new violation, or when a baseline entry
no longer applies, so the baseline cannot go stale.

- To lower a baseline, fix the violation and delete its entry in the same commit.
- To move a line that a baseline names, update the entry to the new location in the same
  commit; the move must not add a violation.
- Never raise a baseline to make a test pass. Fix the code, or discuss the exception in the
  pull request.

## Pull requests

- One variable per pull request, so a result change has one cause.
- Split related work into stacked pull requests, each based on the previous one, and put a
  retarget note in every body that says which pull request it is stacked on and that it
  must be retargeted to `main` when its base merges.
- Any change in results (identifications, FDR, quantities, speed or memory) needs benchmark
  evidence in the description, produced as described in
  [`docs/dev/benchmark-reproducibility.md`](docs/dev/benchmark-reproducibility.md).
- Design decisions and their rationale go in the commit message and the description, not
  in code comments.

## Citations

Credit prior work as described in
[`docs/dev/scientific-citations.md`](docs/dev/scientific-citations.md).
