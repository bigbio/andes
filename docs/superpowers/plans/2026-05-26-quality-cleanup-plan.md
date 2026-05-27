# Quality cleanup (PR-Q1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land a single low-risk cleanup PR on `feat/quality-perf-id-rate` → `dev` that removes 32 dangling `Xxx.java:LINE` references in non-test source, neutralizes stale "port of MS-GF+" framing in module headers + CLI help, renames the `MSGFRUST_RSS_PROBE` env var (legacy-compatible), fixes all 37+ stable clippy warnings, lifts CI lint from advisory to required, and deletes 2 shipped design specs.

**Architecture:** Six in-PR commits (Groups 1-6 from the design spec) plus one out-of-repo memory update (Group 7, already completed by the controller during the brainstorm phase). Logic-preserving — `precursor_cal_bit_identical` regression gate is the safety net. Parity test files (`tests/*_java_parity.rs`, `tests/gf_bsa_parity.rs`, `tests/*_match_java.rs`) are NOT touched — their identity IS Java parity.

**Tech Stack:** Rust 2024 edition pinned to 1.87.0 (`rust-toolchain.toml`), cargo workspace, clippy (stable), `cargo test --release --workspace`, GitHub Actions CI.

**Spec:** `docs/superpowers/specs/2026-05-26-quality-cleanup-design.md`

---

## File map

**Group 1 — dangling Java refs (8 non-test files, 32 refs):**
- `crates/input/src/mzml.rs:63, 351` (2 refs)
- `crates/output/src/pin.rs:354, 417` (2 refs)
- `crates/search/src/mass_calibrator.rs:176` (1 ref)
- `crates/search/src/psm.rs:77, 92, 232, 247, 248, 445` (6 refs)
- `crates/search/src/match_engine.rs:346, 466, 479, 515, 691, 692, 789, 823, 825, 901, 975, 1324` (11 refs)
- `crates/scoring/src/scoring/scored_spectrum.rs:196, 223, 245, 901, 1239` (5 refs)
- `crates/scoring/src/scoring/psm_score.rs:45` (1 ref)
- `crates/msgf-rust/src/bin/msgf-rust.rs:990, 1008, 1118, 1331` (4 refs)

**Group 2 — stale framing:**
- `crates/search/src/lib.rs`, `crates/scoring/src/lib.rs`, `crates/output/src/lib.rs`, `crates/input/src/lib.rs`, `crates/model/src/lib.rs` — top-of-file `//!` headers
- `crates/msgf-rust/src/bin/msgf-rust.rs` — CLI `--help` strings (specifically `#[command(about = ...)]` and any `#[arg(help = ...)]` that compares behavior to Java)

**Group 3 — identifier renames:**
- `crates/msgf-rust/src/bin/msgf-rust.rs` — `MSGFRUST_RSS_PROBE` env var → support `MSGF_RSS_PROBE` AS WELL (accept both for one release)

**Group 4 — clippy 37+ warnings (per crate):**
- `crates/model/src/aa_set.rs:269` (1 warning: manual `split_once`)
- `crates/scoring/src/param_model.rs:365` (1 `map_or`)
- `crates/scoring/src/scoring/scored_spectrum.rs` (12 warnings: 6 complex types, 4 `map_or`, 1 too-many-args, 1 loop index)
- `crates/scoring/src/scoring/scored_spectrum.rs:133-134` (doc list items)
- `crates/search/src/precursor_cal.rs:95` (1 dead `mut`)
- `crates/search/src/match_engine.rs:297, 415` (1 too-many-args, 1 `map_or`)
- `crates/search/src/sa_walk.rs:165` (1 `?` rewrite)
- `crates/output/src/tsv.rs:45, 64, 125` (3 too-many-args)
- `crates/msgf-rust/src/bin/msgf-rust.rs` (13 warnings: 11 doc-indentation, 1 loop counter, 1 misc)

**Group 5 — CI lint required:**
- `.github/workflows/ci.yml` — drop `continue-on-error: true` from the `lint` job

**Group 6 — delete shipped specs:**
- `docs/superpowers/specs/2026-05-23-iter39-docs-rewrite-design.md` — DELETE
- `docs/superpowers/plans/2026-05-23-iter39-docs-rewrite.md` — DELETE

**Group 7 — out-of-repo (DONE during brainstorm):**
- `~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/MEMORY.md` — already updated
- `~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/project_pr_a_precursor_cal_shipped.md` — created
- `~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/project_quality_cleanup_pr_q1_active.md` — created
- `~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/project_next_sub_projects_sequencing.md` — created

Verification only at task 7.

---

## Pre-flight (verify before Task 1)

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git branch --show-current  # must be feat/quality-perf-id-rate
git log origin/dev..HEAD --oneline | wc -l  # expect 2 (a8ad6ddd + 55cff3fa)
git status --short  # expect clean tree
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | grep -E "^test result|error" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -10
# Expect: all `test result:` lines show `0 failed`.
```

If any non-skipped test fails, STOP — pre-flight failed.

---

## Task 1: Group 1 — Scrub dangling `.java:LINE` references

**Files:**
- Modify: `crates/input/src/mzml.rs`
- Modify: `crates/output/src/pin.rs`
- Modify: `crates/search/src/mass_calibrator.rs`
- Modify: `crates/search/src/psm.rs`
- Modify: `crates/search/src/match_engine.rs`
- Modify: `crates/scoring/src/scoring/scored_spectrum.rs`
- Modify: `crates/scoring/src/scoring/psm_score.rs`
- Modify: `crates/msgf-rust/src/bin/msgf-rust.rs`

**Rule:** Replace each `Xxx.java:LINE` or `Xxx.java` citation with intent-only text. Preserve the surrounding sentence's semantic meaning. Pattern:
- `// foo (DBScanner.java:534)` → `// foo (Java parity)`
- `// Java's NewScoredSpectrum.java:253 …` → `// Java parity: …`
- `/// MSGFPlus.java post-cal block` → `/// matching Java's post-cal block`

DO NOT touch:
- `crates/search/tests/gf_java_parity.rs`
- `crates/search/tests/match_engine_java_parity.rs`
- `crates/search/tests/gf_bsa_parity.rs`
- `crates/model/tests/*_match_java.rs`
- `docs/parity-analysis/**`

- [ ] **Step 1: Inventory and confirm exact ref count**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
grep -rEn "\.java:[0-9]+|\.java\b" crates/ --include='*.rs' 2>/dev/null \
  | grep -v "tests/.*java_parity\|tests/gf_bsa_parity\|tests/.*_match_java" \
  | tee /tmp/q1-task1-refs.txt | wc -l
```

Expected: 32 lines (matches the design spec).

- [ ] **Step 2: Scrub `crates/input/src/mzml.rs`**

Open the file. Find line 63:
```rust
// `msutil/ActivationMethod.java` — we map each to one of our five
```
Replace with:
```rust
// Java parity for activation method names — we map each to one of our five
```

Find line 351:
```rust
            // Selection rule (mirrors `StaxMzMLParser.java:595-605`):
```
Replace with:
```rust
            // Selection rule (Java parity):
```

- [ ] **Step 3: Scrub `crates/output/src/pin.rs`**

Find line 354:
```rust
    // enzN, enzC, enzInt — C-4 (2026-05-19): Java DirectPinWriter.java:199-203
```
Replace with:
```rust
    // enzN, enzC, enzInt — C-4 (2026-05-19): Java parity
```

Find line 417:
```rust
    // emits one accession per index — matching Java DirectPinWriter.java:237.
```
Replace with:
```rust
    // emits one accession per index — Java parity.
```

- [ ] **Step 4: Scrub `crates/search/src/mass_calibrator.rs`**

Find line 176:
```rust
/// `MSGFPlus.java` post-cal block). No-op when stats are unreliable or
```
Replace with:
```rust
/// matching Java's post-cal block). No-op when stats are unreliable or
```

- [ ] **Step 5: Scrub `crates/search/src/psm.rs`**

Find line 77:
```rust
    /// `DirectPinWriter.java:237`.
```
Replace with:
```rust
    /// (Java parity for PIN protein-list emission.)
```

Find line 92:
```rust
    /// `DBScanScorer.getScore` returns `node + edge` and `DBScanner.java:533`
```
Replace with:
```rust
    /// Java's score returns `node + edge` (Java parity)
```

Find line 232:
```rust
    /// Java's `DBScanner.java:540` (`size < n OR score == worst → add`).
```
Replace with:
```rust
    /// Java parity (`size < n OR score == worst → add`).
```

Find lines 247-248:
```rust
                    // R-1 (2026-05-18): Java's DBScanner.java:540 keeps tied
                    // PSMs at capacity (and DBScanner.java:745 keeps SpecE
```
Replace with:
```rust
                    // R-1 (2026-05-18): Java parity — keeps tied
                    // PSMs at capacity (and keeps SpecE
```

Find line 445:
```rust
        // (DBScanner.java:540 raw-score retention; DBScanner.java:745 SpecE
```
Replace with:
```rust
        // (Java parity — raw-score retention; SpecE
```

- [ ] **Step 6: Scrub `crates/search/src/match_engine.rs`**

This file has 11 refs. Use the inventory from Step 1 to locate each line. For each:
1. Use `grep -n "\.java:" crates/search/src/match_engine.rs` to confirm current text.
2. Replace `Xxx.java:LINE` patterns with `Java parity` or `Java's behavior` depending on grammar fit.
3. Preserve surrounding comment context — only the citation itself goes.

Example transformations (apply to each of the 11 refs):

```rust
// per-SpecKey raw-score retention (DBScanner.java:534).
```
→
```rust
// per-SpecKey raw-score retention (Java parity).
```

```rust
// Java's `DBScanner.java:619-621` reads
```
→
```rust
// Java parity reads
```

```rust
// `DirectPinWriter.java:165` does
```
→
```rust
// Java parity does
```

```rust
// Java parity (PSMFeatureFinder.java:51-54): feature-counting uses a
```
→
```rust
// Java parity: feature-counting uses a
```

After all 11 replacements, verify:
```bash
grep -c "\.java:" crates/search/src/match_engine.rs
# Expect: 0
```

- [ ] **Step 7: Scrub `crates/scoring/src/scoring/scored_spectrum.rs`**

5 refs at lines 196, 223, 245, 901, 1239. Apply same replacement pattern. Special case for line 901:
```rust
/// `astral-speed/src/main/java/edu/ucsd/msjava/msutil/Spectrum.java`.
```
→
```rust
/// (Java parity for spectrum filtering semantics.)
```

After:
```bash
grep -c "\.java" crates/scoring/src/scoring/scored_spectrum.rs
# Expect: 0
```

- [ ] **Step 8: Scrub `crates/scoring/src/scoring/psm_score.rs`**

Find line 45:
```rust
/// Mirrors Java's `DBScanner.java:513` call: fromIndex=1, toIndex=n+1 →
```
Replace with:
```rust
/// Java parity call: fromIndex=1, toIndex=n+1 →
```

- [ ] **Step 9: Scrub `crates/msgf-rust/src/bin/msgf-rust.rs`**

4 refs at lines 990, 1008, 1118, 1331. Same pattern. For line 990:
```rust
    // (NewScorerFactory.java line ~120). For (CID, HighRes, Tryp, TMT) this
```
→
```rust
    // (Java parity for scorer factory routing). For (CID, HighRes, Tryp, TMT) this
```

After all 4, verify:
```bash
grep -c "\.java" crates/msgf-rust/src/bin/msgf-rust.rs
# Expect: 0
```

- [ ] **Step 10: Final verification — zero dangling java refs in non-test code**

```bash
grep -rEn "\.java:[0-9]+|\.java\b" crates/ --include='*.rs' 2>/dev/null \
  | grep -v "tests/.*java_parity\|tests/gf_bsa_parity\|tests/.*_match_java"
```

Expected output: empty. If anything appears, fix it before committing.

Also verify parity tests untouched:
```bash
git diff -- crates/search/tests/gf_java_parity.rs crates/search/tests/match_engine_java_parity.rs crates/search/tests/gf_bsa_parity.rs crates/model/tests/chemistry_constants_match_java.rs crates/model/tests/standard_aa_masses_match_java.rs crates/model/tests/common_mod_masses_match_java.rs
# Expect: empty (no diffs)
```

- [ ] **Step 11: Run workspace tests**

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | grep -E "^test result|error" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -10
```

Expected: every `test result:` shows `0 failed`. Comment-only changes do not affect test outcomes.

- [ ] **Step 12: Commit**

```bash
git add crates/
git commit -m "$(cat <<'COMMIT_EOF'
chore: scrub 32 dangling .java:LINE references in non-test source

The Java source tree was removed in commit b4565b8e during the
Rust-cutover; the inline citations to specific Java line numbers now
point at code that does not exist in this repo. Replace each citation
with intent-only "Java parity" comments. Preserves semantic meaning;
removes the broken hyperlinks.

Parity-test files (tests/*_java_parity.rs, tests/gf_bsa_parity.rs,
tests/*_match_java.rs) untouched — their identity is Java parity and
the citations are load-bearing documentation.

8 non-test files touched, 32 refs replaced, 0 functional changes.
COMMIT_EOF
)"
```

Expected: commit created.

---

## Task 2: Group 2 — Neutralize "port of MS-GF+" framing

**Files:**
- Modify: `crates/search/src/lib.rs`
- Modify: `crates/scoring/src/lib.rs`
- Modify: `crates/output/src/lib.rs`
- Modify: `crates/input/src/lib.rs`
- Modify: `crates/model/src/lib.rs`
- Modify: `crates/msgf-rust/src/bin/msgf-rust.rs` (CLI help strings only)

**Rule:** Module headers (`//!`) and CLI `--help` strings that introduce a module/flag by reference to Java code should switch to neutral framing. The codebase is post-cutover; we ship `msgf-rust`, not a "port".

**Keep:** `README.md` and `DOCS.md` provenance sections that explain the project's lineage in user-facing context. Those stay.

- [ ] **Step 1: Inventory headers + help strings with stale framing**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
# crate-lib headers
head -10 crates/search/src/lib.rs crates/scoring/src/lib.rs crates/output/src/lib.rs crates/input/src/lib.rs crates/model/src/lib.rs

# CLI help strings
grep -nE "(MS-GF\+|MSGFPlus|port of.*MS-GF|Java MS-GF|mirrors? Java)" crates/msgf-rust/src/bin/msgf-rust.rs
```

Capture the output for Step 2.

- [ ] **Step 2: Edit each module header**

For each of the five `crates/*/src/lib.rs` files, if the top `//!` doc block opens with phrases like "Port of Java MS-GF+ X" or "Rust reimplementation of MSGFPlus", replace the opening sentence with a neutral description of what the crate does. The rest of the doc block stays.

Example (`crates/search/src/lib.rs`):

Current style (if present):
```rust
//! Port of Java MS-GF+ database search engine.
//!
//! Re-exports the public search surface.
```

Neutral:
```rust
//! Peptide database search engine: candidate enumeration,
//! precursor matching, scoring, and PSM aggregation.
//!
//! Re-exports the public search surface.
```

Apply analogous neutral framing to:
- `crates/scoring/src/lib.rs` ("Scoring model, ion prediction, and generating-function DP")
- `crates/output/src/lib.rs` ("Output writers: Percolator PIN, TSV")
- `crates/input/src/lib.rs` ("Input readers: MGF, mzML, FASTA")
- `crates/model/src/lib.rs` ("Core domain types: spectra, peptides, modifications, amino-acid sets, masses")

If a file does NOT have a stale "port of" opener, leave it alone.

- [ ] **Step 3: Edit CLI `--help` strings**

In `crates/msgf-rust/src/bin/msgf-rust.rs`, find `#[command(about = ...)]` near the `Cli` struct. If it mentions Java behavior comparison, replace with a behavior-only description.

Example:
```rust
about = "Rust port of MS-GF+: database search of MGF/mzML spectra against FASTA",
```
→
```rust
about = "msgf-rust: database search of MGF/mzML spectra against FASTA",
```

Then walk through the `#[arg(...)]` attributes. Any `help = "..."` string that explicitly says "matches Java -X behavior" or "Java MS-GF+ default" gets reworded to describe what the flag does without the comparison. Mention of Java numeric legacy values (`-protocol 0`, etc.) **stays** because that's user-facing migration info.

- [ ] **Step 4: Verify CLI still parses + tests pass**

```bash
cargo build --release -p msgf-rust 2>&1 | tail -3
./target/release/msgf-rust --help 2>&1 | head -5
# Expect: builds clean; --help opens with neutral about line.

cargo test --release -p msgf-rust 2>&1 | grep -E "^test result" | tail -5
# Expect: all PASS.
```

- [ ] **Step 5: Commit**

```bash
git add crates/
git commit -m "$(cat <<'COMMIT_EOF'
chore: neutralize "port of MS-GF+" framing in headers and CLI help

The codebase is post-cutover; new contributors should read crate-lib
top-of-file doc comments as descriptions of what each crate does, not
as port-bookkeeping. CLI --help strings that compared behavior to
Java's command-line options now describe behavior directly.

README.md and DOCS.md provenance sections kept (those are intentional
user-facing project lineage). docs/parity-analysis/** kept.

5 crate-lib headers + msgf-rust CLI help touched.
COMMIT_EOF
)"
```

---

## Task 3: Group 3 — Identifier renames + legacy compat

**Files:**
- Modify: `crates/msgf-rust/src/bin/msgf-rust.rs`

- [ ] **Step 1: Locate the `MSGFRUST_RSS_PROBE` env var**

```bash
grep -n "MSGFRUST_RSS_PROBE" crates/msgf-rust/src/bin/msgf-rust.rs
```

Expected: 1-3 sites (var read + maybe doc).

- [ ] **Step 2: Add legacy compat support**

Find the `log_rss` function (or equivalent that reads the env var). Replace the env-var read with both names:

```rust
fn log_rss(label: &str) {
    let new_name = std::env::var_os("MSGF_RSS_PROBE");
    let legacy = std::env::var_os("MSGFRUST_RSS_PROBE");
    if legacy.is_some() && new_name.is_none() {
        eprintln!(
            "WARN: MSGFRUST_RSS_PROBE is deprecated; use MSGF_RSS_PROBE \
             (legacy name accepted in this release, will be removed next)"
        );
    }
    if new_name.is_none() && legacy.is_none() {
        return;
    }
    // ... existing RSS-reading logic unchanged ...
}
```

If the original function used a different control-flow (e.g., early return when the var is unset), preserve that flow — only the env-var name reading changes.

- [ ] **Step 3: Update any in-source doc references to use the new name**

```bash
grep -n "MSGFRUST_RSS_PROBE" crates/msgf-rust/src/bin/msgf-rust.rs
```

For each remaining reference, if it's a doc comment, update to mention the new name with the legacy note. Example:
```rust
/// Memory probe (set MSGF_RSS_PROBE=1; legacy MSGFRUST_RSS_PROBE accepted).
```

- [ ] **Step 4: Verify**

```bash
cargo build --release -p msgf-rust 2>&1 | tail -3
# Sanity check both env-var names:
MSGF_RSS_PROBE=1 ./target/release/msgf-rust --help 2>&1 | grep -E "^startup\s|RSS" | head -3
# (header should print)
MSGFRUST_RSS_PROBE=1 ./target/release/msgf-rust --help 2>&1 | grep -E "WARN.*deprecated|^startup" | head -3
# (should print deprecation warning AND the rss-probe header)
```

- [ ] **Step 5: Commit**

```bash
git add crates/msgf-rust/src/bin/msgf-rust.rs
git commit -m "$(cat <<'COMMIT_EOF'
chore: rename MSGFRUST_RSS_PROBE -> MSGF_RSS_PROBE (legacy accepted)

The "MSGFRUST_" prefix dates from an early iter-era naming and doesn't
match the binary's identity (msgf-rust). Switch to MSGF_RSS_PROBE and
keep the legacy name accepted for this release with a deprecation
warning on stderr. The legacy name will be removed in the next quality
cleanup.

Side-effect-only env var; no functional change.
COMMIT_EOF
)"
```

---

## Task 4: Group 4 — Clippy + unused-lints sweep

This task is the largest. Sub-divided into Tasks 4a-4d by warning class. After each sub-task, run the relevant `cargo clippy` and verify counts drop.

### Task 4a: Auto-fixable simplifications (`map_or`, `?`, `split_once`, indentation)

**Files (per the clippy inventory):**
- `crates/model/src/aa_set.rs` (1 split_once)
- `crates/scoring/src/param_model.rs` (1 map_or)
- `crates/scoring/src/scoring/scored_spectrum.rs` (4 map_or, 2 doc indentation)
- `crates/search/src/match_engine.rs` (1 map_or)
- `crates/search/src/sa_walk.rs` (1 ? rewrite)
- `crates/msgf-rust/src/bin/msgf-rust.rs` (11 doc indentation)

- [ ] **Step 1: Apply per-crate `clippy --fix`**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
for c in model scoring search output msgf-rust; do
  cargo clippy --fix --lib -p "$c" --allow-dirty --allow-staged 2>&1 | tail -3
done
```

cargo-clippy will auto-apply the fixable lints (`map_or`, manual `split_once`, `?` rewrite, some doc-indent cases). Manual lints that don't have a machine-applicable fix remain.

- [ ] **Step 2: Verify fixes look correct**

```bash
git diff --stat | head -10
# Expect: ~5-10 files changed with small line counts.

# Sanity-check one of the rewrites:
grep -nE "manual.*split_once|map_or" crates/model/src/aa_set.rs crates/scoring/src/param_model.rs
```

If any `clippy --fix` result looks semantically wrong, revert that hunk with `git checkout <file>` and apply the fix manually instead.

- [ ] **Step 3: Workspace tests**

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | grep -E "^test result|error" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -10
```

Expected: 0 failed.

- [ ] **Step 4: Stage but don't commit yet** (commit at end of Task 4)

```bash
git add crates/
```

### Task 4b: Complex-type aliases in scored_spectrum.rs

**Files:**
- Modify: `crates/scoring/src/scoring/scored_spectrum.rs`

Six warnings at lines 108, 233, 272, 367, 390, 672 about "very complex type used". Introduce 1-2 type aliases near the top of the file that name the recurring complex type.

- [ ] **Step 1: Identify the recurring shape**

```bash
grep -B 1 "very complex type" /tmp/clippy-output.log 2>/dev/null \
  || cargo clippy --lib -p scoring 2>&1 | grep -A 8 "complex type" | head -40
```

Pattern (typical): `Vec<(Partition, Vec<(IonType, Vec<f32>)>)>` — the segment-partition cache. May also be a `&[(K, V)]` slice variant.

- [ ] **Step 2: Add a `type SegmentPartitionCache = ...;` near the top**

Open `crates/scoring/src/scoring/scored_spectrum.rs`. Find the existing `use ...;` block (lines 1-50 area). After the imports, before the first item, add:

```rust
/// Per-segment partition entries: `(Partition, Vec<(IonType, log-probs)>)`.
pub(crate) type SegmentPartitionCache = Vec<(Partition, Vec<(IonType, Vec<f32>)>)>;
```

If a slice-borrow shape is also complained-about, also add:
```rust
pub(crate) type SegmentPartitionSlice<'a> = &'a [(Partition, Vec<(IonType, Vec<f32>)>)];
```

- [ ] **Step 3: Substitute the alias at each warning site**

For each of the 6 lines flagged by clippy, replace the inline complex type with the alias. Example:

Before:
```rust
fn compute(...
    segment_partition_cache: &Vec<(Partition, Vec<(IonType, Vec<f32>)>)>,
) -> ... {
```

After:
```rust
fn compute(...
    segment_partition_cache: SegmentPartitionSlice<'_>,
) -> ... {
```

(Or `&SegmentPartitionCache` if the lifetime form doesn't fit.)

- [ ] **Step 4: Verify clippy is happy**

```bash
cargo clippy --lib -p scoring 2>&1 | grep "complex type" | wc -l
# Expect: 0
```

- [ ] **Step 5: Tests**

```bash
cargo test --release -p scoring 2>&1 | grep -E "^test result" | tail -3
# Expect: 0 failed.
```

- [ ] **Step 6: Stage**

```bash
git add crates/scoring/src/scoring/scored_spectrum.rs
```

### Task 4c: `too_many_arguments` refactors (5 sites)

**Files:**
- Modify: `crates/scoring/src/scoring/scored_spectrum.rs` (2 sites: line 381 has 11/7, line 669 has 8/7)
- Modify: `crates/search/src/match_engine.rs` (1 site: line 297 has 8/7)
- Modify: `crates/output/src/tsv.rs` (3 sites: lines 45, 64, 125)

**Pattern:** Group the shared args into a small struct passed by `&` reference; keep the caller side ergonomic.

- [ ] **Step 1: Refactor `scored_spectrum.rs:381` (11-arg fn)**

Locate the function (likely `Self::new` or `Self::compute_caches`). Identify which 3-5 args are passed together everywhere it's called. Common groupings:

```rust
struct ScoredSpectrumBuildContext<'a> {
    spec: &'a Spectrum,
    scorer: &'a RankScorer,
    charge: u8,
    fragment_tolerance_da: f64,
    deconv_peaks: Option<&'a [(f64, f32)]>,
}
```

Then change the function signature from 11 args to ~6 (the new ctx struct + the remaining standalone args).

Update all callers (use `cargo build` errors to find them):
```bash
cargo build -p scoring 2>&1 | grep -E "error\[E" | head
```

- [ ] **Step 2: Refactor `scored_spectrum.rs:669` (8-arg fn)**

Similar approach. If the function is `directional_node_score_inner`, the args fall into:
- Spectrum data: `peaks`, `ranks`, `precursor_filtered`
- Scoring context: `segment_partition_cache`, `scorer`, `nominal_mass`, `parent_mass`, etc.

Group whichever feels cohesive. Don't force one cohesive grouping if the args are genuinely independent — `#[allow(clippy::too_many_arguments)]` with a one-line justification is acceptable for hot-path functions where wrapping in a struct hurts readability.

- [ ] **Step 3: Refactor `match_engine.rs:297` (8-arg fn)**

This is in `PreparedSearch::run_chunk_inner`. The args are inherent to the search loop; `#[allow(clippy::too_many_arguments)]` with a comment is probably the right call here since the function is private and not called from many places.

```rust
#[allow(
    clippy::too_many_arguments,
    reason = "private inner driver; args reflect the search-loop state"
)]
fn run_chunk_inner(
    ...
) -> Vec<TopNQueue> { ... }
```

- [ ] **Step 4: Refactor `tsv.rs:45, 64, 125` (3 writer fns)**

Likely `write_tsv`, `write_psm_row`, etc. Args fall into:
- Output target: `writer`
- Data: `spectra`, `queues`, `candidates`, `params`, `idx`
- Format: `spec_file_name`, `use_mgf_specid`

Group into:
```rust
struct TsvWriteContext<'a> {
    spectra: &'a [Spectrum],
    queues: &'a [TopNQueue],
    candidates: &'a [Candidate],
    params: &'a SearchParams,
    idx: &'a SearchIndex,
}
```

Or alternatively, since this is public API across crate boundaries, use `#[allow(clippy::too_many_arguments)]` with a justification: "Writer API mirrors PIN writer; grouping into a context struct would diverge."

Pick whichever produces fewer touched call sites.

- [ ] **Step 5: Workspace tests**

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | grep -E "^test result|error" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -10
```

Expected: 0 failed.

- [ ] **Step 6: Stage**

```bash
git add crates/
```

### Task 4d: Dead `mut`, loop counter, doc indentation, remaining warnings

**Files:**
- Modify: `crates/search/src/precursor_cal.rs` (line 95: dead `mut`)
- Modify: `crates/scoring/src/scoring/scored_spectrum.rs` (line 693: loop index)
- Modify: `crates/msgf-rust/src/bin/msgf-rust.rs` (lines 179-183, 923, 1059, 1129-1135: doc indentation + loop counter)

- [ ] **Step 1: Fix dead `mut` in `precursor_cal.rs`**

Open `crates/search/src/precursor_cal.rs` at line 95. Find the `let mut ...` that isn't actually mutated. Remove the `mut`:

```rust
// before
let mut deviations: Vec<f64> = values.iter().map(|v| (v - center).abs()).collect();
// after
let deviations: Vec<f64> = values.iter().map(|v| (v - center).abs()).collect();
```

- [ ] **Step 2: Fix loop-index warning in `scored_spectrum.rs:693`**

This says "the loop variable `seg` is used to index `segment_partition_cache`". Replace `for seg in 0..cache.len() { let entry = &cache[seg]; ... }` with `for entry in &cache { ... }` (using `iter().enumerate()` if the index is also needed).

- [ ] **Step 3: Fix the 11 doc-indentation warnings in `msgf-rust.rs`**

Lines 179-183 and 1129-1135 are in doc-comment blocks (probably bullet lists). Reformat the bullets so the second line aligns with the first character after `* ` or `- `:

Before:
```rust
    /// * **First item:** description
    ///  description continues
```
After:
```rust
    /// * **First item:** description
    ///   description continues
```

(Note: 3 spaces after `///` for second line to align with the text after `* `.)

Apply to all flagged lines.

- [ ] **Step 4: Fix loop-counter warning at `msgf-rust.rs:1059`**

The warning says "the variable `seen` is used as a loop counter". Replace with the recommended pattern (e.g., `.enumerate()` or a separate counter outside the loop).

- [ ] **Step 5: Confirm clippy is clean**

```bash
cargo clippy --workspace --release 2>&1 | grep -cE "^warning:"
# Expect: 0 (or VERY close to 0 — any residual would be in transitive dep build script noise, which we can't fix)
```

- [ ] **Step 6: Workspace tests**

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | grep -E "^test result|error" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -10
```

Expected: 0 failed.

- [ ] **Step 7: Commit Task 4 (all sub-tasks)**

```bash
git add crates/
git commit -m "$(cat <<'COMMIT_EOF'
chore: fix all clippy warnings (workspace)

Brings the workspace to clippy-clean on stable 1.87.0 so the CI lint
job can be lifted from advisory to required.

Changes by class:
- map_or simplifications (6 sites): mechanical rewrite
- complex-type aliases (6 sites): SegmentPartitionCache/Slice
- too_many_arguments (5 sites): context structs OR justified allow
- doc-list indentation (15 sites): align bullet continuations
- unused_mut (1 site): drop unused mut
- ? rewrite, manual split_once, loop-counter, loop-index: per clippy hint

No functional behavior change; PIN/TSV bit-identical regression gate
in tree (precursor_cal_bit_identical) is the verification.
COMMIT_EOF
)"
```

---

## Task 5: Group 5 — Lift CI lint to required

**Files:**
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Locate the lint job's `continue-on-error`**

```bash
grep -n "continue-on-error\|lint:" .github/workflows/ci.yml | head -10
```

Should show the `lint:` job near line 75-80 with a `continue-on-error: true` immediately under it.

- [ ] **Step 2: Remove the line**

Open `.github/workflows/ci.yml`. Find:

```yaml
  lint:
    name: Lint (clippy + rustfmt)
    runs-on: ubuntu-latest
    # Advisory only — the iter1-38 codebase isn't fmt-clean / clippy-clean
    # yet (~11k lines of fmt churn pending). Surfaces the warnings without
    # blocking PRs while that cleanup is sequenced separately.
    continue-on-error: true
```

Replace with:

```yaml
  lint:
    name: Lint (clippy + rustfmt)
    runs-on: ubuntu-latest
```

(Both the `continue-on-error` line and the trailing comment block become obsolete.)

- [ ] **Step 3: Confirm the lint job still passes the test locally**

The CI lint job typically runs `cargo clippy --workspace --release -- -D warnings`. Simulate:

```bash
cargo clippy --workspace --release -- -D warnings 2>&1 | tail -10
```

Expected: `Finished` with no errors. If clippy fails, return to Task 4 — something was missed.

- [ ] **Step 4: Also verify rustfmt is clean (if the job runs it)**

```bash
grep "rustfmt\|cargo fmt" .github/workflows/ci.yml
```

If `cargo fmt --check` is part of the job, run it locally:

```bash
cargo fmt --check 2>&1 | head -20
```

If it fails, run `cargo fmt --all` and stage the formatting changes. Fmt changes can be folded into THIS commit since they're part of "make lint required".

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/ci.yml
git diff --cached --stat | head
git commit -m "$(cat <<'COMMIT_EOF'
ci: lift lint job from advisory to required

After the workspace clippy clean-up landed in the preceding commits,
the lint job can become a real PR gate. Drop continue-on-error: true
and the explanatory comment block.

Going forward, new clippy warnings or rustfmt drift will block PRs.
COMMIT_EOF
)"
```

---

## Task 6: Group 6 — Delete shipped design specs

**Files:**
- Delete: `docs/superpowers/specs/2026-05-23-iter39-docs-rewrite-design.md`
- Delete: `docs/superpowers/plans/2026-05-23-iter39-docs-rewrite.md`

- [ ] **Step 1: Verify the files exist and the iter39 work shipped**

```bash
ls docs/superpowers/specs/2026-05-23-*.md docs/superpowers/plans/2026-05-23-*.md
git log --oneline | grep -iE "iter39|docs.rewrite" | head -5
```

Expected: both files present; git log shows the iter39 merge (PR #30).

- [ ] **Step 2: Delete both files**

```bash
git rm docs/superpowers/specs/2026-05-23-iter39-docs-rewrite-design.md \
       docs/superpowers/plans/2026-05-23-iter39-docs-rewrite.md
```

Note: this uses `git rm` so the deletion is staged automatically.

- [ ] **Step 3: Confirm nothing references the deleted files**

```bash
grep -rEn "2026-05-23-iter39-docs-rewrite" docs/ crates/ README.md DOCS.md .github/ 2>/dev/null
```

Expected: empty. (If anything points at the deleted files, update the reference.)

- [ ] **Step 4: Commit**

```bash
git diff --cached --stat
git commit -m "$(cat <<'COMMIT_EOF'
docs: remove shipped iter39 design+plan specs

The iter39 docs-rewrite spec and plan shipped via PR #30 in 2026-05-23.
Now that the feature is in dev and being relied on, the design docs
no longer need to be discoverable in the repo. Their lineage is in
git history.

Future protocol: when a docs/superpowers/{specs,plans}/*.md file
references a feature that has fully shipped and closed any deferred
gate, remove it in the next quality cleanup.
COMMIT_EOF
)"
```

---

## Task 7: Final verification + push + open PR

- [ ] **Step 1: Confirm commit count**

```bash
git log origin/dev..HEAD --oneline
# Expect 8 commits:
#   1. a8ad6ddd docs: remove BUG_REVIEW.md; move CLI_MIGRATION.md to docs/  (pre-existing)
#   2. 55cff3fa docs(spec): PR-Q1 quality cleanup design + finalize CLI_MIGRATION refs  (pre-existing)
#   3. Group 1: java refs scrub
#   4. Group 2: framing neutralized
#   5. Group 3: env var rename
#   6. Group 4: clippy clean
#   7. Group 5: CI lint required
#   8. Group 6: shipped specs removed
```

- [ ] **Step 2: Full workspace test sweep**

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts 2>&1 | tee /tmp/q1-final-tests.log | grep -E "^test result|error" | grep -vE "0 passed.*0 failed.*0 ignored" | tail -15
```

Expected: every `test result:` shows `0 failed`. No errors.

- [ ] **Step 3: Bit-identical regression gate**

```bash
cargo test --release -p msgf-rust --test precursor_cal_bit_identical 2>&1 | tail -5
```

Expected: `test result: ok. 1 passed`.

- [ ] **Step 4: Confirm CI lint will pass under -D warnings**

```bash
cargo clippy --workspace --release -- -D warnings 2>&1 | tail -5
```

Expected: `Finished` with no errors.

- [ ] **Step 5: Confirm auto-memory still consistent**

```bash
ls ~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/project_pr_a_precursor_cal_shipped.md \
   ~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/project_quality_cleanup_pr_q1_active.md \
   ~/.claude/projects/-Users-yperez-work-msgfplus-workspace/memory/project_next_sub_projects_sequencing.md
```

Expected: all 3 present. (Group 7 was done during brainstorm; verification only.)

- [ ] **Step 6: Push the branch**

```bash
git push -u origin feat/quality-perf-id-rate 2>&1 | tail -5
```

Expected: branch pushed; URL printed.

- [ ] **Step 7: Open the PR**

```bash
gh pr create --base dev --head feat/quality-perf-id-rate \
  --title "chore: quality cleanup (Q1) — dangling Java refs, clippy clean, lint required" \
  --body "$(cat <<'PR_BODY'
## Summary

Post-cutover code-quality sweep. First of three sequential sub-projects
(Q1 quality → S1 speed → I1 ID-rate +5%/dataset).

Logic-preserving: PIN/TSV output for `--precursor-cal off` is identical
to dev (sorted-row regression gate in tree).

## What changed (6 commits)

- **Group 1 (java refs scrub):** 32 dangling `Xxx.java:LINE` citations
  in non-test source replaced with intent-only "Java parity" comments.
  Parity-test files (`tests/*_java_parity.rs`, `tests/gf_bsa_parity.rs`,
  `tests/*_match_java.rs`) untouched.
- **Group 2 (framing):** 5 crate-lib `//!` headers + CLI `--help`
  strings reworded to describe behavior directly (not as a port).
- **Group 3 (env var):** `MSGFRUST_RSS_PROBE` → `MSGF_RSS_PROBE`,
  legacy name accepted with deprecation warning for one release.
- **Group 4 (clippy):** All workspace warnings cleaned. New type
  aliases (`SegmentPartitionCache`, etc.), 5 `too_many_arguments`
  refactors / justified `#[allow]`, dead `mut`, doc indentation, etc.
- **Group 5 (CI):** Lint job lifted from `continue-on-error: true` to
  required.
- **Group 6 (docs):** Removed 2 shipped design specs from
  `docs/superpowers/`.

## What's NOT in scope

- Speed work (PR-S1, separate brainstorm)
- ID-rate work (PR-I1, multi-PR research project)
- Parity test files (deliberately preserved)
- `docs/parity-analysis/notes/` (current iter notes)

## Verification

- `cargo test --release --workspace` green under existing CI skip list
- `cargo clippy --workspace --release -- -D warnings` clean
- `precursor_cal_bit_identical` regression gate green
- Auto-memory updated (out-of-repo) with PR-A merged status + Q1/S1/I1 sequencing

Spec: `docs/superpowers/specs/2026-05-26-quality-cleanup-design.md`
Plan: `docs/superpowers/plans/2026-05-26-quality-cleanup-plan.md`
PR_BODY
)"
```

Expected: PR URL printed. Record the PR number.

- [ ] **Step 8: Verify CI starts**

```bash
sleep 30
gh pr view --json number,statusCheckRollup --jq '{number, checks: [.statusCheckRollup[]? | {name, status, conclusion}]}'
```

Expected: PR open; CI checks `IN_PROGRESS` or starting. Watch for `Lint (clippy + rustfmt)` to now be a hard gate (not skipped).

---

## Self-review

I checked the plan against the spec section-by-section:

**1. Spec coverage:**
- Group 1 (dangling Java refs) → Task 1 ✓
- Group 2 (stale framing) → Task 2 ✓
- Group 3 (identifier renames) → Task 3 ✓
- Group 4 (clippy + unused sweep) → Task 4 (4a-4d) ✓
- Group 5 (CI lint required) → Task 5 ✓
- Group 6 (remove shipped specs) → Task 6 ✓
- Group 7 (auto-memory) → Pre-done during brainstorm; verified at Task 7 Step 5 ✓
- All ship criteria → Task 7 Steps 2-4 ✓

**2. Placeholder scan:** Scanned for "TBD", "TODO", "fill in", "implement later". None present. Every Task 4 sub-task references a specific file/line from the clippy inventory.

**3. Type consistency:** `SegmentPartitionCache` introduced in Task 4b is used by name in subsequent steps. CI lint job name consistent (`Lint (clippy + rustfmt)`). Commit messages refer to the same commit SHAs (`a8ad6ddd`, `55cff3fa`) used in pre-flight expectations.

**Known soft spots:**
- The exact `cargo clippy --fix` output in Task 4a may vary slightly across clippy versions. If a `--fix` rewrite produces semantically suspect code, Step 2 of Task 4a documents the manual-revert procedure.
- The CLI `--help` strings in Task 2 are inspected by `head` and `grep` rather than enumerated up-front — the implementer reads the actual current content. The plan doesn't pre-script the exact replacements because the strings can drift between plan-writing and execution; the rule is "replace any Java-comparison phrasing with behavior-only".
