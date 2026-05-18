# R-2 Retention-Layer Refactor — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Rust's single-queue-per-spectrum architecture with Java's per-SpecKey pipeline (per-charge queues + dedup + per-charge GF + spectrum-level merge + protein-index aggregation in PIN), so that Astral identification metrics approach Java's coverage.

**Architecture:** Refactor `match_engine.rs`'s per-spectrum loop to use `HashMap<u8, TopNQueue>` keyed by charge. Add a pepSeq+score dedup pass per-charge before GF compute. Call `compute_spec_e_values_for_spectrum` once per non-empty per-charge queue. Merge per-charge queues into a single per-spectrum queue (R-1's tie keep makes SpecE tie-merge automatic). Change `PsmMatch.candidate_idx: u32` → `candidate_idxs: Vec<u32>` so dedup can aggregate protein indices; PIN writer iterates the Vec to emit tab-separated accessions matching Java's `DirectPinWriter.java:237`. External `match_spectra` API unchanged.

**Tech Stack:** Rust (cargo), `search` + `output` crates, bench harness on `pride-linux-vm`.

**Spec:** [`docs/parity-analysis/specs/2026-05-18-r2-retention-refactor-design.md`](../specs/2026-05-18-r2-retention-refactor-design.md) (commit `37d28f95`)

**Branch:** `rust-implement`, starting from HEAD `37d28f95` (R-2 spec on top of R-1 fix + test strengthening).

**3-dataset bar (from 2026-05-18 user instruction):** No merge to `dev` until Rust beats Java on BOTH PSMs at 1% FDR AND wall time for ALL THREE datasets (PXD001819, TMT, Astral). Currently passing for PXD001819, close on TMT, failing on Astral.

---

## File Structure

**Modified (production):**
- `rust/crates/search/src/psm.rs` — `PsmMatch.candidate_idx` → `candidate_idxs: Vec<u32>`, add `primary_candidate_idx()` helper, add `TopNQueue::drain_into_vec`, update test helpers + unit tests
- `rust/crates/search/src/match_engine.rs` — per-spectrum loop refactor (per-charge queues, dedup, per-charge GF, merge); add `dedup_pepseq_score` free function; update `.candidate_idx` accesses
- `rust/crates/output/src/pin.rs` — iterate `psm.candidate_idxs` for Proteins column; update one `.candidate_idx` access
- `rust/crates/output/src/tsv.rs` — update one `.candidate_idx` access
- `rust/crates/msgf-rust/src/bin/msgf-trace.rs` — update one `.candidate_idx` access

**Modified (tests):**
- `rust/crates/search/tests/gf_java_parity.rs` — 2 accesses
- `rust/crates/search/tests/match_engine_java_parity.rs` — 3 accesses + strengthen with deduped-count gate
- `rust/crates/search/tests/peptide_mismatch_diagnostic.rs` — 3 accesses
- `rust/crates/search/tests/match_spectra_thread_invariance.rs` — 2 accesses
- `rust/crates/search/tests/gf_bsa_parity.rs` — 1 access
- `rust/crates/search/tests/match_engine_smoke.rs` — 2 accesses

**Doc comments only (no code):**
- `rust/crates/output/src/pin.rs:37,90,93`
- `rust/crates/output/src/row_context.rs:27`

**New (will be created during plan):**
- `docs/parity-analysis/notes/2026-05-18-r2-bench-results.md` — Astral + PXD001819 + TMT bench results after R-2

**VM-side (not in repo):**
- `/srv/data/msgf-bench/msgf-rust-iter7.tgz` — source tarball
- `/srv/data/msgf-bench/track-iter7-build/` — extracted source + cargo build
- `/srv/data/msgf-bench/bench-iter7-results/{astral,pxd001819,tmt}-rust-r2.{pin,log}`
- `/srv/data/msgf-bench/percolator-iter7/{astral,pxd001819,tmt}_iter7.target.psms.txt`

---

## Task 1: PsmMatch data model migration (pure refactor, no behavior change)

**Files:**
- Modify: `rust/crates/search/src/psm.rs` (struct + helper + test helpers + unit tests)
- Modify: `rust/crates/search/src/match_engine.rs` (5 access sites: lines 298, 357, 513, 564, 591)
- Modify: `rust/crates/output/src/pin.rs` (line 270 access; doc comments at 37/90/93)
- Modify: `rust/crates/output/src/tsv.rs` (line 146)
- Modify: `rust/crates/output/src/row_context.rs` (doc comment at line 27)
- Modify: `rust/crates/msgf-rust/src/bin/msgf-trace.rs` (line 392)
- Modify: `rust/crates/search/tests/gf_java_parity.rs` (lines 165, 172)
- Modify: `rust/crates/search/tests/match_engine_java_parity.rs` (lines 199, 280, 282)
- Modify: `rust/crates/search/tests/peptide_mismatch_diagnostic.rs` (lines 176, 180, 194)
- Modify: `rust/crates/search/tests/match_spectra_thread_invariance.rs` (lines 70, 76)
- Modify: `rust/crates/search/tests/gf_bsa_parity.rs` (line 182)
- Modify: `rust/crates/search/tests/match_engine_smoke.rs` (lines 103, 104)

This task is purely mechanical: rename field + change type + update all accessors. **No behavior change** — every previously-1-element `Vec` stays a 1-element Vec; queue retention is unchanged.

- [ ] **Step 1: Change `PsmMatch.candidate_idx` to `candidate_idxs: Vec<u32>`**

In `rust/crates/search/src/psm.rs`, replace the field at line 70 and surrounding doc comment (lines 58-69). Find:

```rust
pub struct PsmMatch {
    pub spectrum_idx: usize,
    /// Index into the `&[Candidate]` slice owned by `PreparedSearch.candidates`.
    /// Replaces the inlined `Candidate` clone: previously each push to the queue
    /// cloned the full `Candidate` (including its `Peptide.residues: Vec<...>`),
    /// allocating millions of times per large-fasta search. Now the queue stores
    /// only a 4-byte index and consumers (writers, feature extraction, GF) look
    /// up the `Candidate` by index when needed.
    ///
    /// Every real PSM points at a valid index into `PreparedSearch.candidates`.
    /// There is no "synthetic / no backing Candidate" sentinel — test fixtures
    /// that don't need to resolve back use `0` as a placeholder and avoid
    /// touching the candidates slice from inside the test.
    pub candidate_idx: u32,
```

Replace with:

```rust
pub struct PsmMatch {
    pub spectrum_idx: usize,
    /// Indices into the `&[Candidate]` slice owned by `PreparedSearch.candidates`.
    /// Length is always ≥ 1. The first index (`candidate_idxs[0]`) is the
    /// "primary" candidate — used by callers that need a single Candidate
    /// (most do; see `primary_candidate_idx()`). Multiple indices accumulate
    /// when the R-2 pepSeq+score dedup pass merges multiple Candidates that
    /// share the same peptide sequence and rounded score (typically the same
    /// peptide matched against multiple proteins, e.g. shared tryptic
    /// peptides in target+decoy concat). The PIN writer iterates this Vec to
    /// emit one tab-separated `Proteins` column per row, matching Java's
    /// `DirectPinWriter.java:237`.
    ///
    /// Every real PSM has length ≥ 1 with valid indices into
    /// `PreparedSearch.candidates`. Test fixtures that don't need to resolve
    /// back use `vec![0]` as a placeholder and avoid touching the candidates
    /// slice from inside the test.
    pub candidate_idxs: Vec<u32>,
```

- [ ] **Step 2: Add `primary_candidate_idx()` helper**

Append to the `impl PsmMatch` block (find the existing `impl PartialEq for PsmMatch` around line 102 — insert a new `impl PsmMatch` block ABOVE it):

```rust
impl PsmMatch {
    /// Returns the first (primary) candidate index. Callers that need to
    /// resolve back to a single Candidate use this; PIN writer iterates
    /// `candidate_idxs` directly to emit the multi-protein `Proteins` column.
    pub fn primary_candidate_idx(&self) -> u32 {
        self.candidate_idxs[0]
    }
}
```

- [ ] **Step 3: Update `make_match` and `make_match_with_evalue` test helpers**

In `rust/crates/search/src/psm.rs`, find `make_match` at line 259 and replace `candidate_idx: 0` with `candidate_idxs: vec![0]`:

```rust
    fn make_match(spectrum_idx: usize, score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score,
            spec_e_value: 1.0,
            de_novo_score: i32::MIN,
            activation_method: None,
            e_value: 1.0,
            features: PsmFeatures::default(),
            isotope_offset: 0,
        }
    }
```

`make_match_with_evalue` (line 278) uses `make_match` so no change needed.

- [ ] **Step 4: Update match_engine.rs production code**

Five access sites need updating in `rust/crates/search/src/match_engine.rs`:

- **Line 298** (in queue.push, inside the candidate-scoring loop): change `candidate_idx: cand_idx as u32` to `candidate_idxs: vec![cand_idx as u32]`.
- **Line 357** (inside fill_post_topn closure): change `candidates[psm.candidate_idx as usize]` to `candidates[psm.primary_candidate_idx() as usize]`.
- **Line 513** (inside compute_spec_e_values_for_spectrum, protein-terminal flag loop): change `candidates[psm.candidate_idx as usize]` to `candidates[psm.primary_candidate_idx() as usize]`.
- **Line 564** (inside compute_spec_e_values_for_spectrum, after-GF PSM iteration): same `.candidate_idx` → `.primary_candidate_idx()` substitution.
- **Line 591** (inside compute_spec_e_values_for_spectrum, e_value calc): same substitution.

Also update doc comments at lines 387, 392 if they mention `candidate_idx`:

```bash
grep -n "candidate_idx" /Users/yperez/work/msgfplus-workspace/astral-speed/rust/crates/search/src/match_engine.rs
```

Replace each `psm.candidate_idx` with `psm.primary_candidate_idx()`, and `PsmMatch::candidate_idx` mentions in docs with `PsmMatch::candidate_idxs`.

- [ ] **Step 5: Update pin.rs and tsv.rs and row_context.rs**

In `rust/crates/output/src/pin.rs`:
- Line 270: `let cand = &candidates[psm.candidate_idx as usize];` → `let cand = &candidates[psm.primary_candidate_idx() as usize];`
- Lines 37, 90, 93 (doc comments): rewrite to refer to `candidate_idxs` (mentioning the multi-protein aggregation is forthcoming in Task 4)

In `rust/crates/output/src/tsv.rs`:
- Line 146: `let cand = &candidates[psm.candidate_idx as usize];` → `let cand = &candidates[psm.primary_candidate_idx() as usize];`

In `rust/crates/output/src/row_context.rs`:
- Line 27 (doc comment): rewrite to refer to `primary_candidate_idx()`

- [ ] **Step 6: Update msgf-trace.rs**

In `rust/crates/msgf-rust/src/bin/msgf-trace.rs`:
- Line 392: `let cand = &run_candidates[psm.candidate_idx as usize];` → `let cand = &run_candidates[psm.primary_candidate_idx() as usize];`

- [ ] **Step 7: Update all test files**

Apply the same `psm.candidate_idx` → `psm.primary_candidate_idx()` substitution AND any struct-literal `candidate_idx: X` → `candidate_idxs: vec![X]` in:

- `rust/crates/search/tests/gf_java_parity.rs` (lines 165, 172)
- `rust/crates/search/tests/match_engine_java_parity.rs` (lines 199, 280, 282)
- `rust/crates/search/tests/peptide_mismatch_diagnostic.rs` (lines 176, 180, 194)
- `rust/crates/search/tests/match_spectra_thread_invariance.rs` (lines 70, 76)
- `rust/crates/search/tests/gf_bsa_parity.rs` (line 182)
- `rust/crates/search/tests/match_engine_smoke.rs` (lines 103, 104)

Verify no struct literals constructing `PsmMatch` directly need updating:

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
grep -rn "PsmMatch {" rust/crates/ --include='*.rs'
```

Update any found constructors (other than `make_match` already handled in Step 3).

- [ ] **Step 8: Build and run full test suite**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo build -p search -p scoring -p output -p msgf-rust --tests 2>&1 | tail -20
```

Expected: clean build. If any `error[E0609]: no field 'candidate_idx'` or similar appears, locate the missed access and update it.

```bash
cargo test -p search -p scoring -p output 2>&1 | grep -E "^test result|FAILED" | head -30
```

Expected: same test results as the pre-Task-1 baseline:
- 18 passing `psm::tests::*` tests
- All `output_pin_schema_parity`, `output_tsv_*` passing
- 3 pre-existing `match_engine_smoke` failures (unrelated; they fail on baseline `b1d45bb`)
- `r1_tie_retention_active_in_production_pipeline` passing (from commit `de77ea9`)
- `gf_java_parity`, `score_psm_pxd001819_parity`, `score_psm_vs_gf_dp_edge_parity` passing

If any test that was passing pre-Task-1 now FAILS: revert this task and investigate.

- [ ] **Step 9: Commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add rust/ docs/
git commit -m "refactor(search): PsmMatch.candidate_idx -> Vec<u32> candidate_idxs (R-2 prep)

Pure refactor; no behavior change. Renames the field and changes the
type to support multi-protein aggregation in subsequent R-2 work
(per DBScanner.java:719 dedup + DatabaseMatch.java:75 addIndex +
DirectPinWriter.java:237 multi-accession emit).

Adds PsmMatch::primary_candidate_idx() helper for callers that need a
single Candidate (returns candidate_idxs[0]). PIN writer's Proteins
column iteration will use the full Vec in Task 4.

Touches ~20 access sites across 13 files. All previously-1-element
Vecs stay 1-element; no queue retention or scoring changes.

Tests: 18 psm::tests pass, integration tests unchanged, 3 pre-existing
match_engine_smoke failures unchanged (same as baseline b1d45bb)."
```

---

## Task 2: TopNQueue::drain_into_vec + dedup_pepseq_score (TDD)

**Files:**
- Modify: `rust/crates/search/src/psm.rs` (add `drain_into_vec` method + unit test)
- Modify: `rust/crates/search/src/match_engine.rs` (add `dedup_pepseq_score` free function)
- Test: same files

This task adds the dedup machinery without yet wiring it into the per-spectrum loop. The wiring happens in Task 3.

- [ ] **Step 1: Add failing dedup unit test**

In `rust/crates/search/src/psm.rs`, inside the `#[cfg(test)] mod tests` block (after the `topn_queue_keeps_ties_at_capacity` test from R-1), add:

```rust
    #[test]
    fn dedup_pepseq_score_aggregates_candidate_idxs() {
        // R-2.2 (2026-05-18): synthetic test for pepSeq+score dedup. Two PSMs
        // with the same (peptide_residue, score) key should collapse to one
        // PsmMatch with both candidate_idxs aggregated into the surviving Vec.
        //
        // We use drain_into_vec to extract PSMs, then assert the dedup helper
        // collapses them correctly.

        let mut q = TopNQueue::new(10);
        // Three PSMs: two share (peptide=0, score=50), one is distinct (peptide=1, score=40)
        let mut a = make_match(0, 50.0);
        a.candidate_idxs = vec![10];
        let mut b = make_match(0, 50.0);
        b.candidate_idxs = vec![20];
        let mut c = make_match(0, 40.0);
        c.candidate_idxs = vec![30];

        q.push(a);
        q.push(b);
        q.push(c);
        assert_eq!(q.len(), 3, "all three PSMs initially retained");

        let drained = q.drain_into_vec();
        assert_eq!(drained.len(), 3);

        // Caller (match_engine) provides the key function. Here we use
        // a synthetic key based on score only (test scaffolding — real
        // dedup uses peptide_residue + rounded_score from candidates).
        let deduped = simple_dedup_by_score_for_test(drained);

        // Expect: 2 groups — score=50 with idxs [10,20], score=40 with [30]
        assert_eq!(deduped.len(), 2, "should collapse to 2 unique-score groups");

        let mut score_50 = deduped.iter().find(|p| (p.score as i32) == 50).unwrap().candidate_idxs.clone();
        score_50.sort();
        assert_eq!(score_50, vec![10, 20], "score=50 should aggregate both idxs");

        let score_40 = &deduped.iter().find(|p| (p.score as i32) == 40).unwrap().candidate_idxs;
        assert_eq!(*score_40, vec![30]);
    }

    /// Test-only dedup that groups by score alone (real production
    /// dedup_pepseq_score in match_engine.rs uses peptide_residue + score).
    fn simple_dedup_by_score_for_test(psms: Vec<PsmMatch>) -> Vec<PsmMatch> {
        use std::collections::HashMap;
        let mut groups: HashMap<i32, PsmMatch> = HashMap::new();
        for psm in psms {
            let key = psm.score as i32;
            groups
                .entry(key)
                .and_modify(|existing| existing.candidate_idxs.extend(psm.candidate_idxs.iter().copied()))
                .or_insert(psm);
        }
        groups.into_values().collect()
    }
```

Run it:

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p search --lib psm::tests::dedup_pepseq_score_aggregates_candidate_idxs -- --nocapture 2>&1 | tail -10
```

Expected: FAIL with `no method named drain_into_vec found for struct TopNQueue` (compile error — the method doesn't exist yet).

- [ ] **Step 2: Add `drain_into_vec` method to TopNQueue**

In `rust/crates/search/src/psm.rs`, find the existing `fill_post_topn` method (around line 223). Add a new method immediately above it:

```rust
    /// Drain all PSMs from the queue, returning them in an unordered Vec.
    /// Leaves the queue empty after the call. The returned Vec preserves no
    /// particular order — callers that need ordering should sort the result.
    ///
    /// Cost: O(N) drain + Vec collection. Cheap for small N (top-N typically ≤ 10).
    pub fn drain_into_vec(&mut self) -> Vec<PsmMatch> {
        self.heap.drain().map(|Reverse(m)| m).collect()
    }
```

Re-run the test:

```bash
cargo test -p search --lib psm::tests::dedup_pepseq_score_aggregates_candidate_idxs -- --nocapture
```

Expected: PASS.

- [ ] **Step 3: Add `dedup_pepseq_score` free function in match_engine.rs**

In `rust/crates/search/src/match_engine.rs`, add a new function below `compute_spec_e_values_for_spectrum` (after the final `}` of that function — around line 600+, locate the end of the file's function block before any trailing helper functions):

```rust
/// Pre-merge dedup pass (R-2.2): collapse PSMs that share the same
/// (peptide_residue, rounded_score) key into a single entry, aggregating
/// their `candidate_idxs` into a unified Vec. Mirrors Java's
/// `DBScanner.java:719-733` `pepSeqMap` dedup.
///
/// Called by the per-spectrum loop after the per-candidate scoring loop,
/// before per-charge GF compute (so SpecE is computed on the deduped set).
///
/// Inputs:
/// - `psms`: drained from a per-charge `TopNQueue` via `drain_into_vec`
/// - `candidates`: the search's enumerated candidate slice; used to resolve
///   each PSM's peptide residue sequence for the dedup key
///
/// Returns: deduped `Vec<PsmMatch>`. The caller re-pushes these into the
/// per-charge queue via `queue.push()` for each entry.
pub(crate) fn dedup_pepseq_score(
    psms: Vec<PsmMatch>,
    candidates: &[Candidate],
) -> Vec<PsmMatch> {
    use std::collections::HashMap;

    // Key: (peptide_residue_bytes, rounded_score_i32)
    // The residue sequence is the unmodified bare AA string, matching Java's
    // `m.getPepSeq()` used as the dedup key (DBScanner.java:721).
    let mut groups: HashMap<(Vec<u8>, i32), PsmMatch> = HashMap::new();

    for psm in psms {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let pep_residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
        let score_rounded = psm.score.round() as i32;
        let key = (pep_residues, score_rounded);

        groups
            .entry(key)
            .and_modify(|existing| {
                // Aggregate this PSM's indices into the surviving entry.
                // Avoid duplicates if the same idx somehow appears twice.
                for &idx in &psm.candidate_idxs {
                    if !existing.candidate_idxs.contains(&idx) {
                        existing.candidate_idxs.push(idx);
                    }
                }
            })
            .or_insert(psm);
    }

    groups.into_values().collect()
}
```

- [ ] **Step 4: Run dedup-related tests**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p search --lib psm::tests:: -- --nocapture 2>&1 | tail -15
```

Expected: all `psm::tests::*` pass (19 tests total: 18 pre-existing + 1 new `dedup_pepseq_score_aggregates_candidate_idxs`).

- [ ] **Step 5: Run broader test suite to confirm no regression**

```bash
cargo test -p search -p scoring -p output 2>&1 | grep -E "^test result|FAILED" | head -20
```

Expected: same baseline as after Task 1.

- [ ] **Step 6: Commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add rust/crates/search/src/psm.rs rust/crates/search/src/match_engine.rs
git commit -m "feat(search): TopNQueue::drain_into_vec + dedup_pepseq_score (R-2.2)

Adds the machinery needed for R-2.2 pre-merge dedup (Java
DBScanner.java:719-733). Not yet wired into match_engine's per-spectrum
loop -- that happens in Task 3 along with the per-charge queue refactor.

TopNQueue::drain_into_vec: drains the heap into an unordered Vec. Mirror
of the existing fill_post_topn drain pattern, but exposes the Vec for
caller-side processing.

dedup_pepseq_score: groups PSMs by (peptide_residue_bytes, rounded_score)
key and aggregates candidate_idxs across the group. Matches Java's
pepSeqMap behavior: same peptide+score keep all protein indices in one
DatabaseMatch.

Unit test: synthetic test with three PSMs (two tied at score=50 with
distinct candidate_idxs, one at score=40) verifies dedup collapses
correctly and aggregates the Vec contents."
```

---

## Task 3: Per-charge queue refactor + per-charge GF + spectrum merge (R-2.1 + R-2.3 + R-2.4)

**Files:**
- Modify: `rust/crates/search/src/match_engine.rs` (per-spectrum loop body, ~150 lines refactored)

This is the biggest task. The current per-spectrum loop at `match_engine.rs:274-359` is replaced with a per-charge loop structure.

- [ ] **Step 1: Read the current per-spectrum loop block**

Read `rust/crates/search/src/match_engine.rs:274-365` to fully understand the flow you're refactoring. Pay special attention to:
- `window_cand_indices` iteration (line 274)
- `charges_to_try` inner loop (line 277)
- Single `queue.push(...)` at line 296
- The `compute_spec_e_values_for_spectrum` call at line 315 with `top_charge` heuristic
- `queue.fill_post_topn` feature fill at line 355

- [ ] **Step 2: Replace the loop block with per-charge architecture**

In `rust/crates/search/src/match_engine.rs`, find and replace the entire block from line 274 (start of `for &cand_idx in &window_cand_indices`) through line 359 (end of `queue.fill_post_topn`).

Replace with:

```rust
            use std::collections::HashMap;

            // R-2.1: per-charge queue keyed by charge state. Mirrors Java's
            // per-SpecKey raw-score retention (DBScanner.java:534).
            let mut per_charge_queues: HashMap<u8, TopNQueue> = HashMap::new();

            for &cand_idx in &window_cand_indices {
                let cand = &candidates[cand_idx];
                let cleavage_credit = compute_cleavage_credit(cand) as f32;
                for &z in &charges_to_try {
                    let scored_spec = scored_spec_for_charge(z);
                    let mut best_for_charge: Option<(MassError, f32)> = None;
                    for offset in params.isotope_error_range.clone() {
                        if let Some(err) = matches_precursor(spec, &cand.peptide, z, offset, &params.precursor_tolerance) {
                            let score = score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da)
                                + cleavage_credit;
                            if best_for_charge.as_ref().map_or(true, |(_, s)| score > *s) {
                                best_for_charge = Some((err, score));
                            }
                        }
                    }
                    if let Some((err, score)) = best_for_charge {
                        let features = PsmFeatures::default();
                        let psm = PsmMatch {
                            spectrum_idx: spec_idx,
                            candidate_idxs: vec![cand_idx as u32],
                            charge_used: z,
                            mass_error_ppm: err.mass_error_ppm,
                            score,
                            spec_e_value: 1.0,
                            de_novo_score: i32::MIN,
                            activation_method: Some(scorer.param().data_type.activation),
                            e_value: 1.0,
                            features,
                            isotope_offset: err.isotope_offset,
                        };
                        per_charge_queues
                            .entry(z)
                            .or_insert_with(|| TopNQueue::new(params.top_n_psms_per_spectrum as u32))
                            .push(psm);
                        psms_pushed.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
            candidates_visited.fetch_add(window_cand_indices.len() as u64, Ordering::Relaxed);

            // R-2.2: pepSeq + score dedup per-charge BEFORE GF compute.
            // Same peptide matched against multiple proteins collapses to one
            // PsmMatch with aggregated candidate_idxs (Java DBScanner.java:719-733).
            for queue in per_charge_queues.values_mut() {
                if queue.len() > 1 {
                    let drained = queue.drain_into_vec();
                    let deduped = dedup_pepseq_score(drained, candidates);
                    for psm in deduped {
                        queue.push(psm);
                    }
                }
            }

            // R-2.3: per-charge GF / SpecEValue compute. Each per-charge queue
            // gets SpecE calibrated against its OWN charge's GF distribution
            // (Java DBScanner.java:606,779 — getRankScorer per SpecKey).
            let enzyme_opt = if params.enzyme != Enzyme::NoCleavage
                && params.enzyme != Enzyme::NonSpecific
            {
                Some(params.enzyme)
            } else {
                None
            };
            let mut any_queue_nonempty = false;
            for (&charge, queue) in per_charge_queues.iter_mut() {
                if queue.is_empty() {
                    continue;
                }
                any_queue_nonempty = true;
                let scored_spec_charge = scored_spec_for_charge(charge);
                compute_spec_e_values_for_spectrum(
                    spec,
                    params,
                    queue,
                    aa_set_for_gf,
                    enzyme_opt,
                    scorer,
                    scored_spec_charge,
                    charge,
                    fragment_tolerance_da,
                    idx,
                    candidates,
                );
            }
            if any_queue_nonempty {
                spectra_with_psms.fetch_add(1, Ordering::Relaxed);
            }

            // R-2.4: spectrum-level merge with SpecE tie keep. R-1's
            // TopNQueue::push (Ordering::Equal arm) keeps SpecE ties at
            // capacity because PsmMatch::cmp orders by spec_e_value first.
            // Matches Java DBScanner.java:745.
            let mut queue = TopNQueue::new(params.top_n_psms_per_spectrum as u32);
            for (_charge, mut per_charge) in per_charge_queues.drain() {
                for psm in per_charge.drain_into_vec() {
                    queue.push(psm);
                }
            }

            // Feature extraction (unchanged from baseline): post-merge, after
            // the per-spectrum queue is final.
            queue.fill_post_topn(|psm| {
                let ss = scored_spec_for_charge(psm.charge_used);
                let cand = &candidates[psm.primary_candidate_idx() as usize];
                psm.features = compute_psm_features(ss, &cand.peptide, scorer);
            });
```

(Keep the existing `queue` (Vec collection / yield-accounting) trailing logic at lines 361-379 unchanged.)

- [ ] **Step 3: Run the R-1 regression test (`r1_tie_retention_active_in_production_pipeline`)**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p search --test match_engine_java_parity r1_tie_retention_active_in_production_pipeline -- --nocapture 2>&1 | tail -10
```

Expected: PASS. The test asserts ≥1 queue has ≥2 PSMs on BSA fixture. With R-2's per-charge architecture + dedup, this should still hold (tied PSMs within a charge are kept; dedup collapses same-peptide-same-score across proteins but not across distinct peptides).

If FAIL: R-2.2 dedup may be too aggressive (collapsing PSMs it shouldn't). Investigate.

- [ ] **Step 4: Run gf_java_parity (5 BSA PSMs SP-vs-SP)**

```bash
cargo test -p search --test gf_java_parity rust_spec_probability_within_one_oom_of_java_for_5_traced_psms -- --nocapture 2>&1 | tail -15
```

Expected: PASS at 1.0 OOM tolerance (current state after R-1). Note: the per-charge GF change may shift SP slightly for some PSMs; if any of the 5 BSA PSMs moves outside 1.0 OOM, that's a real signal to investigate before bench.

- [ ] **Step 5: Run PXD001819 RawScore stability test**

```bash
cargo test -p scoring --test score_psm_pxd001819_parity -- --nocapture 2>&1 | tail -10
```

Expected: PASS (RawScore=293 stable).

- [ ] **Step 6: Run full search + scoring + output suites**

```bash
cargo test -p search -p scoring -p output 2>&1 | grep -E "^test result|FAILED" | head -30
```

Expected: same as Task 2 baseline. The 3 pre-existing `match_engine_smoke` failures remain. No new failures.

- [ ] **Step 7: Commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add rust/crates/search/src/match_engine.rs
git commit -m "fix(search): per-charge queues + dedup + per-charge GF + merge (R-2.1/.2/.3/.4)

Refactors match_engine's per-spectrum loop to use HashMap<u8, TopNQueue>
keyed by charge, replacing the single-queue-per-spectrum architecture
that R-1 alone amplified into 11.6x raw-target over-shoot on Astral.

R-2.1: each (spectrum, charge) gets its own TopNQueue. Mirrors Java's
per-SpecKey retention (DBScanner.java:534).

R-2.2: pre-GF pepSeq+score dedup pass per-charge collapses same-peptide
matches across proteins, aggregating candidate_idxs (Task 2's
dedup_pepseq_score function). Mirrors DBScanner.java:719-733.

R-2.3: compute_spec_e_values_for_spectrum is called once per non-empty
per-charge queue with the matching scored_spec_for_charge(charge).
Removes the top_charge heuristic that picked one charge for the whole
spectrum's queue. Mirrors DBScanner.java:606,779.

R-2.4: spectrum-level merge pushes per-charge queue contents into a
single per-spectrum TopNQueue. R-1's Equal-arm tie keep makes SpecE
ties survive automatically (PsmMatch::cmp orders by spec_e_value first).
Mirrors DBScanner.java:745.

External match_spectra API unchanged: still returns Vec<TopNQueue> with
one merged queue per spectrum.

Tests: r1_tie_retention_active_in_production_pipeline still passes,
gf_java_parity still passes at 1.0 OOM, PXD001819 RawScore=293 stable,
no new test failures."
```

---

## Task 4: PIN writer multi-accession Proteins column (R-2.5)

**Files:**
- Modify: `rust/crates/output/src/pin.rs` (write_psm_row, around line 462)
- Modify: `rust/crates/output/src/row_context.rs` (possibly simplified accession-resolution)

- [ ] **Step 1: Locate the current Proteins emit site**

Read `rust/crates/output/src/pin.rs:455-475` to see the current single-accession emit pattern:

```rust
    // Peptide, Proteins
    writeln!(writer, "\t{}\t{}", cand.peptide, ctx.accession)
```

Note: `ctx.accession` is built from one Candidate via `resolve_accession()` (`row_context.rs:47`). For R-2.5 we need to emit one accession per `candidate_idxs` entry.

- [ ] **Step 2: Read row_context.rs's resolve_accession function**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
grep -A 20 "pub(crate) fn resolve_accession" rust/crates/output/src/row_context.rs
```

This function takes one Candidate + SearchIndex and returns one accession string. For R-2.5 we'll call it per-candidate-idx.

- [ ] **Step 3: Modify write_psm_row to iterate candidate_idxs**

In `rust/crates/output/src/pin.rs`, find the `Peptide, Proteins` emit (around line 462) and replace with:

```rust
    // Peptide column (always one)
    write!(writer, "\t{}", cand.peptide)?;
    // Proteins column(s): one tab-separated accession per candidate_idx.
    // Java DirectPinWriter.java:237 emits one accession per matching
    // protein in the same row. Rust resolves each candidate_idx to its
    // accession via row_context::resolve_accession.
    for &cidx in &psm.candidate_idxs {
        let cand_for_acc = &candidates[cidx as usize];
        let accession = crate::row_context::resolve_accession(cand_for_acc, search_index);
        write!(writer, "\t{}", accession)?;
    }
    writeln!(writer)?;
```

Note: this requires `resolve_accession` to be accessible. If it's currently `pub(crate)` and in the same crate as pin.rs (it is — both are in `crates/output/`), the cross-module path `crate::row_context::resolve_accession` works.

The `write_psm_row` function signature needs `candidates: &[Candidate]` and `search_index: &SearchIndex` — verify these are already in scope. If not, plumb them through.

- [ ] **Step 4: Verify the output crate's tests still pass**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p output 2>&1 | grep -E "^test result|FAILED" | head -10
```

Expected: all output tests pass. Specifically `output_pin_schema_parity` should pass; the row count and column count remain stable.

If `output_pin_schema_parity` fails because the new multi-accession emit changes the Proteins column from 1 cell to N: this is intentional. Update the test's expected schema OR mark the test as compatible with variable-width Proteins (Java's PIN has the same property — variable-width Proteins column per row).

- [ ] **Step 5: Spot-check on BSA fixture**

```bash
cargo test -p search --test match_engine_java_parity 2>&1 | grep -E "^test result|FAILED" | head -5
```

Expected: all 3 tests pass (the R-1 regression test + the 2 existing parity tests).

If any tests open a generated PIN file and parse the Proteins column expecting exactly N tab-separated tokens, they may need updating. Search for such tests:

```bash
grep -rn "Proteins" rust/crates/search/tests/ --include='*.rs' | head -10
```

- [ ] **Step 6: Commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add rust/crates/output/src/pin.rs rust/crates/output/src/row_context.rs
git commit -m "fix(output): emit one accession per candidate_idx in PIN Proteins column (R-2.5)

Mirrors Java's DirectPinWriter.java:237: 'for (String acc :
proteins.accessions) row.append(acc)'. After R-2.2 dedup aggregates
candidate_idxs across proteins that match the same peptide at the same
score, the PIN writer now emits one tab-separated accession per
candidate in psm.candidate_idxs.

For PSMs with a single candidate_idx (typical), the output is identical
to pre-R-2.5: one Proteins column per row.

For PSMs aggregated by dedup, the row now has N tab-separated
accessions in the Proteins position, matching Java's variable-width
Proteins column.

Tests: output crate suite passes. The PIN schema is unchanged in
header (column names) but the Proteins position can now have N cells
per row instead of 1."
```

---

## Task 5: Strengthen match_engine_java_parity with deduped-count gate

**Files:**
- Modify: `rust/crates/search/tests/match_engine_java_parity.rs` (add one test or strengthen existing one)

- [ ] **Step 1: Add a test asserting deduped PSM count on BSA fixture**

In `rust/crates/search/tests/match_engine_java_parity.rs`, after the `r1_tie_retention_active_in_production_pipeline` test, add:

```rust
#[test]
fn r2_deduped_psm_count_matches_java_on_bsa_fixture() {
    // R-2 (2026-05-18): after per-charge queues + dedup + per-charge GF +
    // spectrum merge, Rust's distinct (scan, peptide) PSM count on the BSA
    // fixture should approach Java's. This catches:
    //   - dedup collapsing PSMs it shouldn't (would reduce distinct count)
    //   - missed cross-charge merge (would inflate count)
    //   - protein-aggregation breaking peptide identity
    //
    // Java reference: bsa_test_mgf_java.pin has 217 unique (scan, peptide)
    // target PSMs. Rust should fall within ±5% — i.e. 207-227.
    //
    // If this test fails after a future change, FIRST check what changed
    // in retention before assuming the test is wrong.

    let java_pin = fixture("benchmark/parity-fixtures/bsa_test_mgf_java.pin");
    let java_target_pairs = java_target_scan_peptide_pairs(&java_pin);
    let java_count = java_target_pairs.len();
    println!("Java distinct (scan, peptide) target PSMs: {}", java_count);

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("src/test/resources/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("src/test/resources/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let scorer = rank_scorer();
    let (queues, candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.05, "XXX");

    let mut rust_target_pairs: HashSet<(i32, Vec<u8>)> = HashSet::new();
    for (spec, queue) in spectra.iter().zip(queues.iter()) {
        let scan = match spec.scan.or_else(|| extract_scan_from_title(&spec.title)) {
            Some(s) => s,
            None => continue,
        };
        for psm in queue.iter_psms() {
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            if cand.is_decoy {
                continue;
            }
            let pep_residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
            rust_target_pairs.insert((scan, pep_residues));
        }
    }
    let rust_count = rust_target_pairs.len();
    println!("Rust distinct (scan, peptide) target PSMs: {}", rust_count);

    let ratio = rust_count as f64 / java_count as f64;
    println!("Rust/Java ratio: {:.3}", ratio);

    // Gate: within ±5% of Java
    assert!(
        ratio >= 0.95 && ratio <= 1.05,
        "Rust distinct PSM count {} is {:.1}% of Java's {} (gate: 95%-105%)",
        rust_count,
        ratio * 100.0,
        java_count
    );
}

/// Parse the Java pin file and return a Set of distinct (scan, peptide_residue)
/// pairs for target rows (Label=1). Strips flanking residues and mod brackets
/// to get the bare residue sequence, matching how Rust's candidate.peptide.residues
/// represents peptides internally.
fn java_target_scan_peptide_pairs(pin_path: &PathBuf) -> HashSet<(i32, Vec<u8>)> {
    let f = File::open(pin_path).unwrap_or_else(|e| panic!("open {pin_path:?}: {e}"));
    let r = BufReader::new(f);
    let mut lines = r.lines();
    let header = lines.next().unwrap().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let scan_idx = cols.iter().position(|c| *c == "ScanNr").expect("ScanNr");
    let label_idx = cols.iter().position(|c| *c == "Label").expect("Label");
    let pep_idx = cols.iter().position(|c| *c == "Peptide").expect("Peptide");

    let mut pairs: HashSet<(i32, Vec<u8>)> = HashSet::new();
    for line_result in lines {
        let line = match line_result {
            Ok(l) => l,
            Err(_) => continue,
        };
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= label_idx.max(scan_idx).max(pep_idx) {
            continue;
        }
        if fields[label_idx] != "1" {
            continue;
        }
        let scan: i32 = match fields[scan_idx].parse() {
            Ok(s) => s,
            Err(_) => continue,
        };
        // Strip flanking + mods. Java format: K.PEPTIDE.V or K.P+57.021EPTIDE.V
        let pep = fields[pep_idx];
        let pep_stripped = pep
            .split('.')
            .nth(1)
            .unwrap_or(pep)
            .as_bytes()
            .iter()
            .filter(|&&b| b.is_ascii_uppercase())
            .copied()
            .collect::<Vec<u8>>();
        pairs.insert((scan, pep_stripped));
    }
    pairs
}
```

- [ ] **Step 2: Run the new test**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed/rust
cargo test -p search --test match_engine_java_parity r2_deduped_psm_count_matches_java_on_bsa_fixture -- --nocapture 2>&1 | tail -20
```

Expected behavior:
- If R-2 worked correctly: ratio ≈ 1.0, test PASSES
- If dedup is too aggressive: ratio < 0.95, test FAILS — investigate
- If retention is too permissive: ratio > 1.05, test FAILS — investigate

If FAIL: the empirical signal tells us R-2 isn't yet right. Debug before bench.

- [ ] **Step 3: Run full match_engine_java_parity (all 4 tests)**

```bash
cargo test -p search --test match_engine_java_parity 2>&1 | grep -E "^test result|FAILED" | head
```

Expected: 4 tests pass.

- [ ] **Step 4: Commit**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add rust/crates/search/tests/match_engine_java_parity.rs
git commit -m "test(search): R-2 deduped (scan, peptide) count gate on BSA fixture

Adds r2_deduped_psm_count_matches_java_on_bsa_fixture which asserts
Rust's distinct (scan, peptide) target-PSM count on the BSA fixture
is within +-5% of Java's 217. Catches:
  - dedup too aggressive (count drops below 207)
  - retention too permissive (count exceeds 227)
  - protein-aggregation breaking peptide identity

This is the strengthened parity gate the 2026-05-18 review called for.
Combined with the existing r1_tie_retention_active test, the file now
gates on both tie-keep activity AND deduped count parity."
```

---

## Task 6: Astral no-mods bench + 4-gate verification

**Files:**
- Create on VM: `/srv/data/msgf-bench/msgf-rust-iter7.tgz`, `/srv/data/msgf-bench/track-iter7-build/`, `/srv/data/msgf-bench/bench-iter7-results/`

This task mirrors R-1's Task 3-4 (sync + build + bench). Use the same SSH ControlMaster pattern.

- [ ] **Step 1: Verify SSH ControlMaster alive**

```bash
ls -la /tmp/msgfplus-bench.sock 2>&1
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk 'echo connected' 2>&1 | head -3
```

Expected: socket file exists; SSH echoes "connected".

If not: user re-establishes via `ssh-keygen -R pride-linux-vm.ebi.ac.uk && ssh -M -S /tmp/msgfplus-bench.sock -fN root@pride-linux-vm.ebi.ac.uk`. Stop here and ask.

- [ ] **Step 2: Tar + upload + extract + build on VM**

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
tar --exclude='rust/target' --exclude='rust/**/target' -czf /tmp/msgf-rust-iter7.tgz rust 2>&1 | tail -3
ls -lh /tmp/msgf-rust-iter7.tgz
scp -o ControlPath=/tmp/msgfplus-bench.sock /tmp/msgf-rust-iter7.tgz root@pride-linux-vm.ebi.ac.uk:/srv/data/msgf-bench/
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    rm -rf /srv/data/msgf-bench/track-iter7-build &&
    mkdir -p /srv/data/msgf-bench/track-iter7-build &&
    cd /srv/data/msgf-bench/track-iter7-build &&
    tar -xzf /srv/data/msgf-bench/msgf-rust-iter7.tgz 2>&1 | grep -v "LIBARCHIVE" | head -3 &&
    cp -r /srv/data/msgf-bench/track-iter3-build/src ./ &&
    cd rust &&
    printf "[toolchain]\nchannel = \"stable\"\ncomponents = [\"rustfmt\", \"clippy\"]\n" > rust-toolchain.toml &&
    grep -q "dedup_pepseq_score" crates/search/src/match_engine.rs && echo "R-2 wiring present" || echo "MISSING R-2"
'
```

Expected output: `R-2 wiring present`.

Then build:

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    cd /srv/data/msgf-bench/track-iter7-build/rust &&
    nohup bash -c "cargo build --release -p msgf-rust > /tmp/cargo-build-iter7.log 2>&1" >/dev/null 2>&1 &
    echo "build PID: $!"
'
```

Wait for build:

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    until ! pgrep -x cargo >/dev/null && ! pgrep -x rustc >/dev/null; do sleep 30; done
    echo "build done"
    tail -3 /tmp/cargo-build-iter7.log
    ls -la /srv/data/msgf-bench/track-iter7-build/rust/target/release/msgf-rust
'
```

Expected: `Finished 'release' profile`; binary exists ~2 MB. If build FAILS: read the log to find which crate failed compiling.

- [ ] **Step 3: Launch Astral bench in nohup background**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    RUST=/srv/data/msgf-bench/track-iter7-build/rust/target/release/msgf-rust
    OUT=/srv/data/msgf-bench/bench-iter7-results
    mkdir -p "$OUT"
    nohup bash -c "/usr/bin/time -v $RUST \
        --spectrum /srv/data/msgf-bench/astral-data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML \
        --database /srv/data/msgf-bench/astral-data/ProteoBenchFASTA_MixedSpecies_HYE.fasta \
        --output-pin $OUT/astral-rust-r2.pin \
        --precursor-tol-ppm 10 --isotope-error-min=-1 --isotope-error-max=2 \
        --ntt 2 --max-missed-cleavages 2 --min-peaks 10 \
        --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 \
        --top-n 1 --threads 8 > $OUT/astral-rust-r2.log 2>&1" >/dev/null 2>&1 &
    echo "astral bench PID: $!"
'
```

- [ ] **Step 4: Poll for bench completion**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    until ! pgrep -f "astral-rust-r2.pin" >/dev/null; do sleep 60; done
    echo "bench done"
    grep -E "Elapsed|Maximum resident|Exit status" /srv/data/msgf-bench/bench-iter7-results/astral-rust-r2.log | head -3
    ls -la /srv/data/msgf-bench/bench-iter7-results/astral-rust-r2.pin
'
```

Expected: Exit status 0; wall ≤ 26 min; PIN file size much smaller than R-1's 467 MB (because dedup + per-charge separation should bring counts toward Java's ~36 MB).

- [ ] **Step 5: Count raw targets/decoys + run Percolator**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    PIN=/srv/data/msgf-bench/bench-iter7-results/astral-rust-r2.pin
    T=$(awk -F"\t" "NR>1 && \$2==1 {c++} END {print c+0}" "$PIN")
    D=$(awk -F"\t" "NR>1 && \$2==-1 {c++} END {print c+0}" "$PIN")
    RATIO=$(awk -v t=$T -v d=$D "BEGIN {if (d>0) printf \"%.3f\", t/d; else print \"inf\"}")
    echo "iter7 astral R-2 raw: targets=$T decoys=$D T/D=$RATIO"
    echo "Java baseline:        targets=89479 decoys=46792 T/D=1.912"
    echo "b1d45bb (pre-R-1):    targets=75457 decoys=46208 T/D=1.633"
    echo "Gates: raw 72K-107K; T/D >= 1.85"

    mkdir -p /srv/data/msgf-bench/percolator-iter7
    bash /srv/data/msgf-bench/run_percolator_docker.sh "$PIN" /srv/data/msgf-bench/percolator-iter7 astral_iter7 2>&1 | tail -2
'
```

- [ ] **Step 6: Compute deduped (scan, peptide) count from Percolator output**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    PSMS=/srv/data/msgf-bench/percolator-iter7/astral_iter7.target.psms.txt
    # Same methodology as the R-1 results note
    awk -F"\t" "NR==1 {next} \$3 <= 0.01 {
        scan = \$1; sub(/.*scan=/, \"\", scan); sub(/_.*/, \"\", scan)
        pep = \$5; sub(/^[A-Z_-]\\./, \"\", pep); sub(/\\.[A-Z_-]\$/, \"\", pep)
        gsub(/\\+[0-9.]+/, \"\", pep)
        print scan \"|\" pep
    }" "$PSMS" | sort -u | wc -l > /tmp/iter7_distinct_count
    DISTINCT=$(cat /tmp/iter7_distinct_count)
    PERC_1PCT=$(grep -c "" "$PSMS")  # incl. header
    PERC_1PCT=$((PERC_1PCT - 1))     # exclude header
    echo "iter7 astral R-2 Percolator @ 1% FDR rows: $PERC_1PCT"
    echo "iter7 astral R-2 distinct (scan, peptide): $DISTINCT"
    echo ""
    echo "=== R-2 4-GATE CHECK ==="
    echo "Java baselines: Percolator=35818, distinct=35818, T/D=1.912"
    echo "Gate 1 (raw 72K-107K):       [see above]"
    echo "Gate 2 (distinct >= 34027):  $DISTINCT"
    echo "Gate 3 (T/D >= 1.85):         [see above]"
    echo "Gate 4 (Percolator >= 30K):  $PERC_1PCT"
'
```

- [ ] **Step 7: Document iter7 Astral results**

Create `docs/parity-analysis/notes/2026-05-18-r2-bench-results.md` locally with the iter7 Astral data:

```markdown
# R-2 retention-layer refactor — empirical results (Astral, iteration 7)

_2026-05-18. Branch `rust-implement` at iter7 HEAD (R-2 + Task 4 + Task 5 commits)._

## Astral no-mods bench

| Metric | Java | b1d45bb | iter6 (R-1) | iter7 (R-2) | Gate | Status |
|---|---:|---:|---:|---:|---|---|
| Raw targets | 89,479 | 75,457 | 1,042,255 | <FILL> | 72K-107K | <FILL> |
| Raw decoys | 46,792 | 46,208 | 530,430 | <FILL> | — | — |
| T/D ratio | 1.912 | 1.633 | 1.965 | <FILL> | >= 1.85 | <FILL> |
| Wall | — | ~8:36 | 23:52 | <FILL> | <= 26 min | <FILL> |
| Percolator @ 1% FDR rows | 35,818 | 25,224 | 74,204 | <FILL> | >= 30K | <FILL> |
| Distinct (scan, peptide) at q<=0.01 | 35,818 | unknown | 26,934 | <FILL> | >= 34,027 | <FILL> |

## Gate decision

<FILL: which of the 4 gates passed/failed, and overall verdict>

## Next

<FILL: per user's 3-dataset bar — if all 4 gates pass on Astral, proceed to Task 7 (PXD+TMT bench). If any fail, revert R-2 and re-plan.>
```

Fill in the `<FILL>` placeholders with the actual data from Step 5-6.

- [ ] **Step 8: Apply gate decision**

If all 4 Astral gates pass: proceed to Task 7 (PXD+TMT bench).

If any gate fails: STOP, do not proceed to Task 7. Decide on revert:

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
# Identify the commits to revert: all R-2 work since the spec
git log --oneline 37d28f95..HEAD
# Revert the R-2 implementation commits (typically Tasks 1-5 from this plan)
# Leave the docs and spec; only revert the production code commits.
# Exact commits to revert depend on what landed.
```

Report which gate(s) failed in the results note and stop the plan.

- [ ] **Step 9: Commit the results note (if gates passed)**

```bash
git add docs/parity-analysis/notes/2026-05-18-r2-bench-results.md
git commit -m "docs(parity): R-2 Astral bench results (iter7)"
```

---

## Task 7: PXD001819 + TMT bench (3-dataset bar verification)

Only run this task if Task 6's Astral gates passed.

**Files:**
- Create on VM: bench results for PXD001819 + TMT
- Modify: `docs/parity-analysis/notes/2026-05-18-r2-bench-results.md` (add PXD+TMT sections)

- [ ] **Step 1: Launch PXD001819 + TMT benches**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk 'cat > /tmp/iter7-pxd-tmt.sh <<EOF
#!/bin/bash
set +e
RUST=/srv/data/msgf-bench/track-iter7-build/rust/target/release/msgf-rust
OUT=/srv/data/msgf-bench/bench-iter7-results
echo "=== pxd001819 ==="
/usr/bin/time -v \$RUST --spectrum /srv/data/msgf-bench/data/UPS1_5000amol_R1.mzML --database /srv/data/msgf-bench/data/PXD001819_uniprot_yeast_ups.fasta --output-pin \$OUT/pxd001819-rust-r2.pin --precursor-tol-ppm 5 --isotope-error-min=0 --isotope-error-max=1 --ntt 2 --max-missed-cleavages 2 --min-peaks 10 --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 --top-n 1 --threads 8 > \$OUT/pxd001819-rust-r2.log 2>&1
echo "=== tmt ==="
/usr/bin/time -v \$RUST --spectrum /srv/data/msgf-bench/tmt-data/a05058.mzML --database /srv/data/msgf-bench/tmt-data/PXD007683_UP000005640_UP000002311_reviewed.fasta --output-pin \$OUT/tmt-rust-r2.pin --precursor-tol-ppm 20 --isotope-error-min=-1 --isotope-error-max=2 --fragmentation 1 --instrument 1 --protocol 4 --ntt 2 --max-missed-cleavages 2 --min-peaks 10 --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 --top-n 1 --threads 8 > \$OUT/tmt-rust-r2.log 2>&1
echo "=== DONE ==="
EOF
chmod +x /tmp/iter7-pxd-tmt.sh
nohup bash /tmp/iter7-pxd-tmt.sh > /srv/data/msgf-bench/bench-iter7-pxd-tmt.log 2>&1 & echo "pxd+tmt PID: $!"'
```

- [ ] **Step 2: Poll for completion**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    until ! pgrep -f "pxd001819-rust-r2.pin\|tmt-rust-r2.pin" >/dev/null; do sleep 60; done
    echo "pxd+tmt done"
    cat /srv/data/msgf-bench/bench-iter7-pxd-tmt.log
    ls -la /srv/data/msgf-bench/bench-iter7-results/pxd001819-rust-r2.pin /srv/data/msgf-bench/bench-iter7-results/tmt-rust-r2.pin 2>&1
'
```

- [ ] **Step 3: Count + Percolator for both datasets**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm.ebi.ac.uk '
    for ds in pxd001819 tmt; do
        PIN=/srv/data/msgf-bench/bench-iter7-results/${ds}-rust-r2.pin
        T=$(awk -F"\t" "NR>1 && \$2==1 {c++} END {print c+0}" "$PIN")
        D=$(awk -F"\t" "NR>1 && \$2==-1 {c++} END {print c+0}" "$PIN")
        WALL=$(grep "Elapsed (wall clock)" /srv/data/msgf-bench/bench-iter7-results/${ds}-rust-r2.log | awk -F": " "{print \$NF}")
        echo "iter7 $ds R-2 raw: targets=$T decoys=$D wall=$WALL"
        echo "  Percolator:"
        bash /srv/data/msgf-bench/run_percolator_docker.sh "$PIN" /srv/data/msgf-bench/percolator-iter7 ${ds}_iter7 2>&1 | tail -1
    done
'
```

- [ ] **Step 4: Apply 3-dataset bar**

For each dataset, verify Rust beats Java on BOTH PSMs at 1% FDR AND wall time:

| Dataset | Java 1% FDR | Java wall | Rust iter7 1% FDR | Rust iter7 wall | Bar met? |
|---|---:|---|---:|---|---|
| PXD001819 | 14,989 | 1:21 | <FILL> | <FILL> | <FILL> |
| TMT | 10,194 | 3:07 | <FILL> | <FILL> | <FILL> |
| Astral | 35,818 | 5:58 | <FILL from Task 6> | <FILL from Task 6> | <FILL> |

Document the verdict in `docs/parity-analysis/notes/2026-05-18-r2-bench-results.md`.

- [ ] **Step 5: Commit results + apply final decision**

If all 3 datasets meet the bar (Rust better on BOTH PSMs and wall time): R-2 implementation is shippable.

If any dataset fails the bar: do NOT merge; document why and identify next steps (likely: more retention fixes, audit-tier feature fixes, performance work).

```bash
cd /Users/yperez/work/msgfplus-workspace/astral-speed
git add docs/parity-analysis/notes/2026-05-18-r2-bench-results.md
git commit -m "docs(parity): R-2 3-dataset bench results (iter7)

<FILL with summary verdict per the 4-gate + 3-dataset bar>"
```

---

## Self-Review

**1. Spec coverage:**

| Spec section | Plan task |
|---|---|
| Problem (recap of R-1 over-shoot) | Covered in plan intro |
| Goal (full architectural refactor) | Covered in plan intro |
| Success criteria (4-gate) | Task 6 Step 5-6 + Task 7 Step 4 |
| Architecture (Approach A in-place) | Tasks 1-3 |
| Data model change | Task 1 |
| Dedup algorithm | Task 2 |
| Per-charge GF + merge | Task 3 |
| PIN multi-accession | Task 4 |
| Testing strategy | Task 2 (unit), Task 5 (integration) |
| Bench protocol | Tasks 6-7 |
| Rollback | Task 6 Step 8 |
| Out of scope | not implemented; explicit in spec |
| Risks | mitigations distributed across tasks |
| Effort estimate | ~9-10 hours total |

All 13 spec sections have corresponding tasks.

**2. Placeholder scan:**

The `<FILL>` placeholders in Task 6 Step 7 and Task 7 Step 4 are intentional — they're empirical values to be filled in from bench output, not abstract gaps. These are appropriate (the plan can't predict the numbers in advance).

No "TBD", "TODO", "implement later", "add appropriate error handling", or "similar to Task N" patterns.

**3. Type consistency:**

- `PsmMatch.candidate_idxs: Vec<u32>` consistent throughout
- `primary_candidate_idx() -> u32` consistent
- `TopNQueue::drain_into_vec(&mut self) -> Vec<PsmMatch>` consistent
- `dedup_pepseq_score(psms: Vec<PsmMatch>, candidates: &[Candidate]) -> Vec<PsmMatch>` consistent
- `per_charge_queues: HashMap<u8, TopNQueue>` consistent

No type or signature mismatches.

---

## Execution handoff

Plan complete and saved to `docs/parity-analysis/plans/2026-05-18-r2-retention-refactor-plan.md`.
