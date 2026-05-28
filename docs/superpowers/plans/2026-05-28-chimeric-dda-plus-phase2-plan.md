# Chimeric DDA+ — Phase 2 Implementation Plan (MS1 targeted-XIC isotope refinement)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use `- [ ]` checkboxes.
>
> Phase 2 of `docs/superpowers/specs/2026-05-28-chimeric-dda-plus-integration-design.md`. Builds on Phase 1 (`feat/chimeric-dda-plus`, Tasks 1-5). Motivated by the Phase 1 bench finding (`docs/parity-analysis/notes/2026-05-28-chimeric-phase1-bench.md`): multi-PSM emission **inflates FDR** (Astral +94% on narrow windows; decoy fraction → 1:1) because spurious co-isolated PSMs aren't filtered. Phase 2 adds the MS1 isotope-envelope check that suppresses them.

**Goal:** For each chimeric PSM, verify its precursor's isotope envelope actually exists in the linked MS1 survey scan, and emit the match quality (isotope KL-divergence + precursor SNR) as **additive** Percolator columns so the rescorer can reject false co-IDs. Target: Astral chimeric drops back to a plausible gain; PXD/TMT retain genuine co-fragmentation gains.

**Architecture:** Under `--chimeric`, load MS1 scans (currently excluded), link each MS2 to its preceding MS1, and compute — per PSM — the observed isotope-envelope intensities at the peptide's theoretical precursor m/z, compared to an averagine theoretical envelope via KL-divergence. v1 uses the single linked MS1 scan (apex check); multi-scan XIC correlation is deferred to a Phase 2b. The feature is additive (new PIN columns); existing columns and the `--chimeric off` path are unchanged.

**Tech Stack:** Rust workspace; quick-xml mzML; Percolator PIN.

## Conventions

- Build: `cargo build --release -p msgf-rust`
- Bit-identical gate (must hold every task; chimeric off ⇒ unchanged): `cargo test --release -p msgf-rust precursor_cal_off_pin_tsv_match_golden_after_sort`
- Clippy: `cargo clippy --workspace --all-targets -- -D warnings`
- VM bench (decision): `--chimeric` on, 3 datasets, compare to the Phase-1 bench (PXD 17,015 / Astral 71,347 / TMT 9,608 @1% FDR with the inflated counts) — Phase 2 should pull Astral down toward plausibility while keeping a real PXD gain.

---

### Task 1: Averagine theoretical isotope distribution

**Files:** Create `crates/model/src/isotope.rs`; register `pub mod isotope;` in `crates/model/src/lib.rs`.

- [ ] **Step 1: Failing test** in `isotope.rs`:
```rust
#[test]
fn averagine_envelope_monotone_decreasing_after_apex_for_small_peptide() {
    // ~1000 Da peptide: monoisotope dominant, isotopes decay.
    let env = averagine_isotope_envelope(1000.0, 4);
    assert_eq!(env.len(), 4);
    let sum: f64 = env.iter().sum();
    assert!((sum - 1.0).abs() < 1e-9, "envelope must be normalized, got {sum}");
    assert!(env[0] > env[1] && env[1] > env[2], "small-peptide envelope should decay from m0");
}

#[test]
fn averagine_apex_shifts_up_for_large_peptide() {
    // ~3000 Da: the +1 isotope rivals/exceeds the monoisotope.
    let env = averagine_isotope_envelope(3000.0, 5);
    assert!(env[1] >= env[0] * 0.8, "large-peptide +1 isotope should be comparable to m0");
}
```
Run: `cargo test --release -p model averagine` → FAIL.

- [ ] **Step 2: Implement** `pub fn averagine_isotope_envelope(neutral_mass: f64, n_isotopes: usize) -> Vec<f64>`. Use the standard averagine approximation: average amino-acid residue C/H/N/O/S composition scaled to `neutral_mass`, then a Poisson/binomial isotope approximation. Acceptable v1: the Poisson-lambda averagine model — `lambda ≈ neutral_mass * 0.000594` (empirical C-count → 13C probability), envelope `p_k = exp(-lambda) * lambda^k / k!` for `k in 0..n_isotopes`, then normalize to sum 1.0. Document the constant's origin (averagine carbon count × natural 13C abundance 1.07%).

- [ ] **Step 3:** `cargo test --release -p model averagine` → PASS. Clippy clean.
- [ ] **Step 4: Commit** `feat(model): averagine theoretical isotope envelope`

---

### Task 2: Load MS1 + link each MS2 to its precursor MS1

**Files:** `crates/input/src/mzml.rs`; new `Ms1Link` side structure (in `crates/input/src/lib.rs` or a small new module).

Design: when chimeric, parse MS1 (ms_level 1) spectra too. Don't feed them to the search as scorable spectra; instead build a side index: `ms1_peaks: Vec<Vec<(f64,f32)>>` (unique MS1 scans, in order) and `ms2_to_ms1: Vec<Option<usize>>` indexed by the returned MS2 spectrum order, pointing at the most-recent preceding MS1. Memory: store MS1 peaks once; MS2→MS1 is an index, not a copy.

- [ ] **Step 1: Failing test** — a small mzML with one MS1 (ms level 1, a couple peaks) followed by two MS2s; assert the reader (in a new `--chimeric`-style MS1 mode) returns 2 MS2 spectra AND an `Ms1Link` where both MS2s link to MS1 index 0, and the MS1 peaks are retrievable. Mirror the existing mzml test harness.
- [ ] **Step 2: Implement.** Add a builder flag (e.g. `with_ms1_capture(true)`) that, instead of discarding ms-level-1 spectra, routes them into a separate `ms1_peaks` vec and records, for each subsequently emitted MS2, the current latest MS1 index into `ms2_to_ms1`. Return both the MS2 `Vec<Spectrum>` and the `Ms1Link { ms1_peaks, ms2_to_ms1 }` (extend the reader's return type or add a sibling method `read_with_ms1`). Keep the default (non-chimeric) path returning MS2-only, unchanged.
- [ ] **Step 3:** test PASS; existing input tests green; the non-MS1 path unchanged (bit-identical gate green).
- [ ] **Step 4: Commit** `feat(mzml): optional MS1 capture + MS2→MS1 linkage for chimeric`

---

### Task 3: Observed precursor isotope envelope + KL feature

**Files:** `crates/search/src/psm.rs` (`PsmFeatures` + new fields), `crates/search/src/match_engine.rs` (compute it), plumb `Ms1Link` into `PreparedSearch`.

- [ ] **Step 1: Add fields** to `PsmFeatures`:
```rust
/// KL-divergence between the observed precursor isotope envelope (from the
/// linked MS1) and the averagine theoretical envelope. High = poor isotope
/// match = likely spurious co-isolation. 0 when MS1/feature unavailable.
pub precursor_isotope_kl: f32,
/// Observed monoisotopic precursor intensity / local MS1 noise (SNR proxy).
pub precursor_snr: f32,
```
Default both to 0.0 in `PsmFeatures::default()`.

- [ ] **Step 2: Failing test** — a unit test of a new helper `precursor_isotope_kl(ms1_peaks, theo_mz, charge, neutral_mass, tol_da) -> (f32, f32)`: synthesize an MS1 peak list with a clean averagine-shaped envelope at the theoretical m/z → assert low KL; synthesize one with only noise (no envelope) → assert high KL. Run → FAIL.
- [ ] **Step 3: Implement the helper** (in `match_engine.rs` or a new `chimeric_features.rs`): for `k in 0..N`, find the nearest MS1 peak to `theo_mz + k*ISOTOPE/charge` within `tol_da`, take its intensity (0 if absent); normalize the observed vector; compute `KL(observed || averagine_isotope_envelope(neutral_mass, N))` with the standard epsilon-guard; SNR = m0 intensity / median MS1 intensity. Return `(kl, snr)`.
- [ ] **Step 4: Plumb + compute.** Thread `Ms1Link` into `PreparedSearch` (store `Option<&Ms1Link>`). In the post-top-N feature fill (or `compute_psm_features`), when chimeric and the MS2 has a linked MS1, call the helper at the PSM's theoretical precursor m/z and set the two feature fields. When not chimeric or no MS1 → leave 0.0.
- [ ] **Step 5:** helper test PASS; bit-identical gate green (off path sets nothing); clippy clean.
- [ ] **Step 6: Commit** `feat(search): precursor isotope-KL + SNR features for chimeric PSMs`

---

### Task 4: Emit the additive PIN columns

**Files:** `crates/output/src/pin.rs` (`write_header` + `write_psm_row`), `crates/output/tests/output_pin_schema_parity.rs`.

- [ ] **Step 1: Failing test** — assert the header contains `PrecursorIsotopeKL` and `PrecursorSNR` (positioned next to `EdgeScore`, before `Peptide`), and that row width grows by 2. Run → FAIL.
- [ ] **Step 2: Implement** — add both columns in `write_header` and write `features.precursor_isotope_kl` / `features.precursor_snr` in `write_psm_row` at the matching position.
- [ ] **Step 3:** `cargo test --release -p output` PASS; **bit-identical gate** — NOTE: adding columns changes the golden PIN width, so the golden must be regenerated (the columns are 0.0 in the non-chimeric golden run; regenerate per the test's documented procedure and confirm only the two new zero-columns are added). Commit the regenerated golden with the code.
- [ ] **Step 4: Commit** `feat(output): additive PrecursorIsotopeKL + PrecursorSNR PIN columns`

---

### Task 5: VM bench — does the MS1 filter deflate the FDR inflation?

**Files:** none (validation; I/the operator drive the VM).

- [ ] **Step 1:** workspace tests (CI skip list incl. `java_fixtures_load`) green; clippy clean; bit-identical golden green (regenerated).
- [ ] **Step 2:** ship + build chimeric-build on the VM; re-run `--chimeric` on all 3 datasets (cal=auto) producing the new PIN with the isotope columns.
- [ ] **Step 3:** Percolator @1% on each. Compare to the Phase-1 inflated baseline:
  - **Expected if working:** Astral chimeric drops from 71,347 toward a plausible small gain over 36,715 (the spurious co-IDs now have high KL → Percolator rejects them); PXD retains a real gain over 14,755; decoy fraction recovers from ~1.2:1 toward the off-mode ratio.
  - **If Astral stays inflated:** the KL feature isn't discriminative enough alone → Phase 3 (greedy shared-fragment rescoring) is required before chimeric is trustworthy; record and proceed to Phase 3 plan.
- [ ] **Step 4: Decide + document** in `docs/parity-analysis/notes/2026-05-28-chimeric-phase2-bench.md`; open/continue the chimeric PR with the off-vs-on + phase1-vs-phase2 tables.

---

## Self-review

- **Spec coverage:** Phase-2 spec items 2a→Task 2, 2b→Tasks 1+3, 2c→Task 4, validation→Task 5. Averagine (Task 1) is the theoretical model 2b needs.
- **Additive-safe:** the only PIN change is two new columns (Task 4); existing columns and the `--chimeric off` path are untouched; the golden is regenerated solely for the two zero-columns.
- **v1 scoping is explicit:** single-linked-MS1 envelope (apex), not multi-scan XIC — that's deferred (Phase 2b) and called out so the implementer doesn't over-build.
- **Decision gate (Task 5):** if the KL feature alone doesn't deflate Astral, that's the documented trigger for Phase 3 — the plan doesn't pretend Phase 2 is guaranteed sufficient.
- **Type consistency:** `precursor_isotope_kl`/`precursor_snr` (`f32`), `averagine_isotope_envelope(f64, usize) -> Vec<f64>`, `Ms1Link { ms1_peaks, ms2_to_ms1 }` used consistently across tasks.
