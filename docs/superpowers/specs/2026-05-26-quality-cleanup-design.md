# Design — Quality cleanup sweep (PR-Q1)

**Date:** 2026-05-26
**Branch:** `feat/quality-perf-id-rate`
**Status:** Spec for review
**First sub-project of three:** Q1 (this) → S1 (speed) → I1 (ID rate)

## Problem

Post-cutover the codebase carries stale historical references and lint debt accumulated across the Java→Rust port iterations. Specifically:

- **42 dangling `Xxx.java:LINE` pointers** in source comments cite Java code that no longer exists in this repo (removed in cutover commit `b4565b8e chore: remove Java tool sources`). They read as broken hyperlinks.
- **File-header `port of Java MS-GF+ X` framing** introduces modules as ports of files that no longer live in-tree; misleading for new contributors.
- **`MSGFRUST_RSS_PROBE` env var** and any remaining `java_*` / `msgf_*` symbol names carry iter-era naming that doesn't reflect the current binary identity (`msgf-rust`).
- **26 clippy warnings** across the workspace, plus a known dead `mut` and undiscovered `unused_*` items. CI lint job runs `continue-on-error: true` because the codebase isn't yet clippy-clean.

These don't affect runtime behavior, but they:
1. Confuse new contributors trying to read context-laden comments that point at non-existent code.
2. Block the CI lint job from being a real gate.
3. Make refactoring noisier than necessary (every modification trips a re-formatter or a stylistic warning).

## Goal

Single low-risk PR that lands a post-cutover quality sweep. Logic-preserving; bit-identical PIN/TSV output for `--precursor-cal off`. Lifts CI lint from advisory to required.

## Non-goals

- Speed or performance work (PR-S1, separate brainstorm).
- ID-rate work (PR-I1, separate multi-PR project).
- Parity test files (`tests/*_java_parity.rs`, `tests/gf_bsa_parity.rs`) — their identity IS Java parity; refs stay.
- `docs/parity-analysis/notes/` iter notes — historical; not edited.
- Renaming production public APIs across crate boundaries.
- Rust edition / toolchain bumps.

## Scope — 7 logical groups (6 in-PR + 1 out-of-repo)

### Group 1 — Dangling Java source pointers (42 refs)

Replace `Xxx.java:LINE` citations with intent-only comments. The semantic intent stays; the broken pointer goes.

Before:
```rust
// per-SpecKey raw-score retention (DBScanner.java:534).
```
After:
```rust
// per-SpecKey raw-score retention (Java parity).
```

Files (counts from initial scan):
- `crates/search/src/match_engine.rs` — 12 refs
- `crates/output/src/pin.rs` — 2 refs
- `crates/input/src/mzml.rs` — 2 refs
- `crates/search/src/psm.rs` — 1 ref
- `crates/search/src/mass_calibrator.rs` — 1 ref
- Others — smaller counts

**Excluded:** `crates/search/tests/gf_java_parity.rs`, `crates/search/tests/match_engine_java_parity.rs`. These tests' purpose is documenting Java parity; their citations are load-bearing.

### Group 2 — Stale "port of MS-GF+" framing

File-header `//!` intros and CLI `--help` strings that introduce modules/flags by reference to Java code go neutral.

Targets:
- `crates/search/src/lib.rs`, `crates/scoring/src/lib.rs`, `crates/output/src/lib.rs` headers
- `crates/msgf-rust/src/bin/msgf-rust.rs` CLI help strings
- A few in-source `//!` modules across `model/`, `search/`

Keep:
- `README.md` provenance section ("evolved from the Java MS-GF+ tradition" or equivalent)
- `DOCS.md` benchmarking-comparison table (explicitly cites Java numbers)
- All `docs/parity-analysis/` content

### Group 3 — Stale identifier renames

- `MSGFRUST_RSS_PROBE` env var → `MSGF_RSS_PROBE` (or just `RSS_PROBE`)
  - Accept BOTH the old and new name during one release; emit a one-line deprecation eprintln if the old name is set, then drop in the next quality cleanup
- Audit for any remaining `java_*` or `msgf_*` named items in source (excluding test fixtures)
- The binary name (`msgf-rust`) and crate name (`msgf-rust`) stay — those are the product identity

### Group 4 — Clippy 26 warnings + `unused_*` sweep

| Warning class | Count | Fix approach |
|---|---:|---|
| `too_many_arguments` (8/7 or 11/7) | 5 | Wrap shared args in a small struct; one cohesive grouping per call site |
| Complex type → `type` alias | 6 | 2-3 reusable type aliases (`SegmentPartitionCache`, etc.) |
| `map_or` simplification | 6 | Mechanical rewrite |
| `doc_list_item_without_indentation` | 4 | Reformat bullet indents |
| `unused_mut` (real dead) | 1 | Drop `mut` |
| Manual `?` rewrite | 1 | Apply |
| Manual `split_once` | 1 | Apply |
| Loop-index borrow | 1 | `iter().enumerate()` |
| Crate summaries | 4 | Mostly auto-fixable via `cargo clippy --fix --lib` |

Additionally:
- Run `cargo +nightly -W unused_variables -W dead_code -W unused_imports --workspace` and clean any findings the stable compiler missed.
- Where a finding is intentional, add `#[allow(...)]` with a one-line justification.

### Group 5 — Lift CI lint to required

`.github/workflows/ci.yml` currently runs the `lint` job with `continue-on-error: true`. After Groups 1-4, the workspace is clippy-clean. Drop the `continue-on-error` so lint becomes a real gate.

### Group 6 — Remove outdated in-repo docs

Tracked docs under `docs/superpowers/` exist for SHIPPED features that no longer need a public spec/plan to reference. Remove:

- `docs/superpowers/specs/2026-05-23-iter39-docs-rewrite-design.md` — iter39 shipped 2026-05-23 in PR #30; design no longer in-flight.
- `docs/superpowers/plans/2026-05-23-iter39-docs-rewrite.md` — same.

Keep:
- `docs/superpowers/specs/2026-05-26-quality-cleanup-design.md` — THIS spec; in-flight.
- All `docs/parity-analysis/notes/2026-05-25-*.md` — referenced by the in-flight ship-gates discussion (precursor-cal G1 still deferred).
- `docs/parity-analysis/snapshots/cal-shifts-2026-05-25.json` — current bench artifact.

Future protocol (documented in this spec so reviewers can apply it): when a `docs/superpowers/{specs,plans}/*.md` file references a feature that has fully shipped + been benched + closed any deferred gate, remove it in the next quality cleanup.

### Group 7 — Update project auto-memory (out-of-repo)

Auto-memory lives at `~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/`. Out of the PR's diff but in the cleanup sweep. To be done by the controller alongside PR-Q1:

- Update `MEMORY.md` index: PR #29 (rust-implement → dev) MERGED, not OPEN; PR #33 (precursor-cal-pr-a → dev) MERGED 2026-05-26; PR #32 (review/bug-hunt → dev) MERGED 2026-05-26.
- Add new entry referencing the 2026-05-25/26 bench numbers (LFQ 14,721 / Astral 36,771 / TMT 9,565 with `--precursor-cal auto`).
- Mark iter32-38 entries as historical / shipped.
- Note the new PR-Q1 / PR-S1 / PR-I1 sequencing.

## File-by-file inventory (estimate)

| File | Change kind | Risk |
|---|---|---|
| `crates/search/src/match_engine.rs` | Java-ref scrub + 1 `too_many_arguments` fix | Medium (hot path) |
| `crates/search/src/mass_calibrator.rs` | Java-ref scrub | Low |
| `crates/search/src/psm.rs` | Java-ref scrub | Low |
| `crates/scoring/src/scoring/scored_spectrum.rs` | 1 `too_many_arguments` fix + complex-type alias | Medium (hot path) |
| `crates/scoring/src/gf/*` | Clippy stylistic | Low |
| `crates/output/src/pin.rs` | Java-ref scrub + 1 `too_many_arguments` fix | Low |
| `crates/output/src/tsv.rs` | Clippy stylistic | Low |
| `crates/input/src/mzml.rs` | Java-ref scrub | Low |
| `crates/msgf-rust/src/bin/msgf-rust.rs` | CLI-help neutral + env var rename | Low |
| `crates/search/src/lib.rs`, `crates/scoring/src/lib.rs`, `crates/output/src/lib.rs` | Header neutral | Low |
| `crates/model/src/*` | Stylistic + 1 `unused_mut` | Low |
| `.github/workflows/ci.yml` | `continue-on-error` removed | Low |

Estimated total: ~30 files modified + 2 files deleted (Group 6), ~200 LOC of comment/identifier/structural change, 0 functional behavior change.

## Verification / ship criteria

| Gate | Threshold | How |
|---|---|---|
| Clippy clean on stable | 0 warnings on `cargo clippy --workspace --release` | CI lint job (now required) |
| Nightly unused-lints clean | 0 (or `#[allow]` justified) | `cargo +nightly -W unused_variables -W dead_code -W unused_imports --workspace` locally |
| Workspace tests | 0 failures under existing skip list | `cargo test --release --workspace -- --skip ...` |
| Off-path bit-identical | `precursor_cal_bit_identical` passes | Already in tree |
| Sanity bench | LFQ / Astral / TMT PSM count within ±5 of pre-cleanup on `--precursor-cal off` | Optional VM run; deferred to reviewer if rayon noise alone explains drift |

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Java-ref scrub accidentally rewords a load-bearing semantic note | Replace IN PLACE preserving comment lines around it; reviewer (or CodeRabbit) flags semantic drift |
| `too_many_arguments` refactor introduces a parameter-ordering bug | The 5 refactors each touch ≤ 1 hot-path function; bench gate catches PSM drift |
| `MSGFRUST_RSS_PROBE` rename breaks an external bench script | Accept BOTH old + new env var name for one release with deprecation eprintln |
| Lifting CI lint surfaces platform-specific warnings (macOS / Windows) | Run `cargo clippy --workspace` locally with `--target x86_64-pc-windows-gnu` and `--target x86_64-apple-darwin` before PR open |

## Sequencing (Q1 only)

```
feat/quality-perf-id-rate (current HEAD: a8ad6ddd)
  ↓
Group 1: Java-ref scrub          (commit 1)
  ↓
Group 2: Header / CLI framing    (commit 2)
  ↓
Group 3: Identifier renames      (commit 3)
  ↓
Group 4: Clippy + unused sweep   (commit 4)
  ↓
Group 5: CI lint required        (commit 5)
  ↓
Group 6: Remove shipped specs    (commit 6)
  ↓
Group 7: Memory update           (out-of-repo, separate)
  ↓
Verification: tests + bit-identical gate + local clippy on 3 platforms
  ↓
Push + open PR-Q1 → dev
```

6 in-PR commits (Groups 1-6) + 1 out-of-repo memory update (Group 7) — keeps reverts easy per-group if any one surfaces an issue.

## Open questions

None — all design points resolved in brainstorming.

## Related documents

- `docs/superpowers/specs/2026-05-25-precursor-cal-ship-design.md` — PR-A spec (the precursor calibrator port)
- `docs/parity-analysis/notes/2026-05-25-precursor-cal-ship-gates.md` — current bench numbers + G1 gate status
- `.github/workflows/ci.yml` — CI policy, including the existing test skip list and the lint job's `continue-on-error`
- `DOCS.md` — primary user-facing reference (touched by Group 2 only where stale framing appears)
