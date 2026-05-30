# Chimeric Two-Pass Cascade Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover most of chimeric's +116% Astral PSM gain at wall time < Java, by replacing the blind wide-window search with a narrow first pass plus an MS1-gated targeted second-peptide search.

**Architecture:** Pass 1 = the existing narrow search (primary peptide/scan, already fast). Pass 2 = for scans where MS1 shows a co-isolated precursor, score a handful of candidates at the MS1-detected co-isolated mass on the residual spectrum (primary peaks removed) and emit a secondary PSM. Both feed one PIN → Percolator. Targets the profiled bottleneck (per-candidate scoring ~65%) by scoring few candidates instead of thousands.

**Tech Stack:** Rust (crates: `search`, `scoring`, `model`, `input`, `output`, `msgf-rust`), cargo, Percolator + entrapment harness on VM.

**Spec:** `docs/superpowers/specs/2026-05-30-chimeric-two-pass-cascade-design.md`

**Reused machinery:** `input::Ms1Link { ms1_peaks: Vec<Vec<(f64,f32)>>, ms2_to_ms1: Vec<Option<usize>> }` (already loaded under `--chimeric`); `chimeric_features::precursor_isotope_match(ms1_peaks, theo_mono_mz, charge, neutral_mass, tol_da, n_isotopes) -> (kl, snr)` (envelope-quality scorer); `model::isotope::averagine_isotope_envelope`; `model::mass::{ISOTOPE, PROTON, H2O}`; `match_engine::{score_psm, psm_edge_score, matched_peak_keys, compute_spec_e_values_for_spectrum}`.

---

## File Structure

- **Create** `crates/search/src/coisolation.rs` — co-isolated-precursor detection from an MS1 isolation window (P1) + the targeted second-peptide search (P2). One responsibility: produce secondary PSMs for a scan.
- **Modify** `crates/search/src/lib.rs` — register `pub(crate) mod coisolation;`.
- **Modify** `crates/msgf-rust/src/bin/msgf-rust.rs` — the two-pass driver in `run()` (P3): after Pass-1 `match_spectra`, when `--chimeric`, run Pass 2 over MS1-gated scans and merge secondary PSMs into the queues before `write_pin`.

---

## Task 1: Co-isolated precursor detector

**Files:** Create `crates/search/src/coisolation.rs`; Modify `crates/search/src/lib.rs`.

- [ ] **Step 1: Register module.** In `crates/search/src/lib.rs` add near the other decls: `pub(crate) mod coisolation;`

- [ ] **Step 2: Write the failing test.** Create `crates/search/src/coisolation.rs`:

```rust
//! Chimeric two-pass cascade: detect co-isolated precursors in an MS2 scan's MS1
//! isolation window (excluding the selected precursor), then run a targeted
//! second-peptide search at each. This is the speed-correct chimeric path: it
//! scores few candidates at MS1-confirmed masses instead of thousands across the
//! blind window (see docs/parity-analysis/notes/2026-05-30-chimeric-cost-profile.md).

use model::mass::{ISOTOPE, PROTON};
use crate::chimeric_features::precursor_isotope_match;

/// A co-isolated precursor detected in the MS1 isolation window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CoIsolated {
    pub mono_mz: f64,
    pub charge: u8,
    pub neutral_mass: f64,
}

/// Detect co-isolated precursors in `ms1_peaks` (m/z-sorted) within the isolation
/// window `[win_lo, win_hi]`, EXCLUDING the envelope at `selected_mz` (the peptide
/// Pass 1 already searched). Tries charges in `charge_range`; accepts an envelope
/// whose averagine KL is below `max_kl`. Returns at most `max_n` (highest-intensity
/// monoisotopic peaks first).
pub(crate) fn detect_coisolated(
    ms1_peaks: &[(f64, f32)],
    win_lo: f64,
    win_hi: f64,
    selected_mz: f64,
    charge_range: std::ops::RangeInclusive<u8>,
    tol_da: f64,
    max_kl: f32,
    max_n: usize,
) -> Vec<CoIsolated> {
    // Candidate monoisotopic peaks = peaks inside the window, sorted by intensity desc.
    let lo_idx = ms1_peaks.partition_point(|&(mz, _)| mz < win_lo);
    let mut cands: Vec<(f64, f32)> = ms1_peaks[lo_idx..]
        .iter().take_while(|&&(mz, _)| mz <= win_hi).copied().collect();
    cands.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut out: Vec<CoIsolated> = Vec::new();
    for &(mz, _inten) in &cands {
        if (mz - selected_mz).abs() <= tol_da { continue; } // skip the selected precursor
        // Don't re-report a peak that's an isotope of an already-accepted envelope.
        if out.iter().any(|c| {
            let d = (mz - c.mono_mz).abs();
            (0..6).any(|k| (d - k as f64 * ISOTOPE / c.charge as f64).abs() <= tol_da)
        }) { continue; }
        // Try charges; accept the lowest-KL charge under max_kl.
        let mut best: Option<(f32, CoIsolated)> = None;
        for z in charge_range.clone() {
            if z == 0 { continue; }
            let neutral = (mz - PROTON) * z as f64;
            let (kl, _snr) = precursor_isotope_match(ms1_peaks, mz, z, neutral, tol_da, 4);
            if kl <= max_kl && best.as_ref().map_or(true, |(bk, _)| kl < *bk) {
                best = Some((kl, CoIsolated { mono_mz: mz, charge: z, neutral_mass: neutral }));
            }
        }
        if let Some((_, c)) = best { out.push(c); }
        if out.len() >= max_n { break; }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::isotope::averagine_isotope_envelope;

    /// Build a synthetic MS1 peak list (m/z-sorted) containing a 4-peak averagine
    /// envelope for `(mono_mz, charge, neutral_mass)` scaled by `scale`.
    fn envelope(mono_mz: f64, charge: u8, neutral: f64, scale: f32) -> Vec<(f64, f32)> {
        let env = averagine_isotope_envelope(neutral, 4);
        (0..4).map(|k| (mono_mz + k as f64 * ISOTOPE / charge as f64, (env[k] as f32) * scale)).collect()
    }

    #[test]
    fn detects_coisolated_excludes_selected() {
        let z = 2u8;
        let selected_mz = 600.0;
        let sel_neutral = (selected_mz - PROTON) * z as f64;
        let co_mz = 600.7; // a second precursor within a ~2 Da window
        let co_neutral = (co_mz - PROTON) * z as f64;
        let mut peaks = envelope(selected_mz, z, sel_neutral, 1000.0);
        peaks.extend(envelope(co_mz, z, co_neutral, 500.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
        assert_eq!(got.len(), 1, "exactly one co-isolated (selected excluded)");
        assert!((got[0].mono_mz - co_mz).abs() < 0.02);
        assert_eq!(got[0].charge, z);
    }

    #[test]
    fn no_coisolation_when_only_selected_present() {
        let z = 2u8;
        let selected_mz = 600.0;
        let peaks = envelope(selected_mz, z, (selected_mz - PROTON) * z as f64, 1000.0);
        let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
        assert!(got.is_empty(), "only the selected precursor -> no co-isolation");
    }
}
```

Verify `averagine_isotope_envelope` return type (slice/Vec of f64) and `precursor_isotope_match`'s exact arg order against the codebase before running; adjust the test harness to match.

- [ ] **Step 3: Run to verify fail.** `cargo test -p search --lib coisolation 2>&1 | tail` → FAIL (type/method missing). After pasting the impl it should pass; split if needed for a true red.

- [ ] **Step 4: Run to verify pass.** `cargo test -p search --lib coisolation 2>&1 | tail` → 2 passed. `cargo build -p search` clean. `cargo clippy -p search 2>&1 | tail -3` clean. (Add `#![allow(dead_code)]` if Task-2 items aren't present yet.)

- [ ] **Step 5: Commit.**
```bash
git add crates/search/src/coisolation.rs crates/search/src/lib.rs
git commit -m "feat(chimeric): MS1 co-isolated precursor detector (cascade P1)"
```

---

## Task 2: Targeted second-peptide search on the residual spectrum

**Files:** Modify `crates/search/src/coisolation.rs`.

This produces the best secondary `PsmMatch` for one `(scan, CoIsolated)` by scoring only candidates near the co-isolated mass on the residual spectrum. It reuses the existing scoring; the new logic is (a) residual peak removal, (b) restricting candidates to the co-isolated mass, (c) one targeted GF.

- [ ] **Step 1: Write the integration test.** Add to `coisolation.rs` tests a test `secondary_search_finds_planted_peptide` that builds a tiny in-memory `PreparedSearch` over 2–3 candidates (mirror the construction used in `crates/search/tests/match_engine_smoke.rs` — read that test for the exact builder), a synthetic `Spectrum` whose peaks are a planted secondary peptide's b/y ions, a `CoIsolated` at that peptide's mass, and asserts `search_secondary(...)` returns a `PsmMatch` whose `primary_candidate_idx()` resolves to the planted peptide. (Read `match_engine_smoke.rs` first; reuse its fixtures verbatim to avoid re-deriving the scorer/index setup.)

- [ ] **Step 2: Run to verify fail.** `cargo test -p search --lib coisolation::tests::secondary 2>&1 | tail` → FAIL (`search_secondary` missing).

- [ ] **Step 3: Implement `search_secondary`.** Add to `coisolation.rs`:

```rust
use model::spectrum::Spectrum;
use model::peptide::Peptide;
use scoring_crate::scoring::{RankScorer, ScoredSpectrum, score_psm, psm_edge_score};
use crate::candidate_gen::Candidate;
use crate::psm::{PsmMatch, PsmFeatures, TopNQueue};
use crate::match_engine::{matched_peak_keys, compute_spec_e_values_for_spectrum};
use crate::search_params::SearchParams;
use crate::search_index::SearchIndex;
use model::aa_set::AminoAcidSet;
use model::enzyme::Enzyme;
use std::collections::BTreeMap;

/// Best secondary PSM for `co` on `spec`, after removing the primary peptide's
/// matched charge-1 b/y peaks (residual). Scores ONLY candidates within
/// `params.precursor_tolerance` of `co.neutral_mass` (the candidate-count cut).
/// Returns `None` if no candidate clears scoring. `bucket_index` maps
/// `nominal(mass - H2O) -> candidate ids` (same as PreparedSearch).
#[allow(clippy::too_many_arguments)]
pub(crate) fn search_secondary(
    spec: &Spectrum,
    primary: &Peptide,
    co: CoIsolated,
    candidates: &[Candidate],
    bucket_index: &BTreeMap<i32, Vec<usize>>,
    scorer: &RankScorer,
    aa_set: &AminoAcidSet,
    enzyme: Option<Enzyme>,
    params: &SearchParams,
    search_index: &SearchIndex,
    fragment_tolerance_da: f64,
) -> Option<PsmMatch> {
    let z = co.charge;
    if z == 0 { return None; }
    // 1. residual spectrum: drop the primary's matched charge-1 b/y peaks.
    let full_ss = ScoredSpectrum::new(spec, scorer, z);
    let claimed = matched_peak_keys(&full_ss, primary, scorer);
    let mut residual = spec.clone();
    residual.peaks.retain(|&(mz, _)| !claimed.contains(&((mz * 1000.0).round() as i64)));
    let res_ss = ScoredSpectrum::new(&residual, scorer, z);

    // 2. candidates within precursor tol of the co-isolated neutral mass.
    let nominal = |m: f64| model::mass::nominal_from(m - model::mass::H2O);
    let tol = params.precursor_tolerance.left.as_da(co.neutral_mass).max(0.01);
    let lo = nominal(co.neutral_mass - tol) - 1;
    let hi = nominal(co.neutral_mass + tol) + 1;
    let mut queue = TopNQueue::new(1);
    for (_nm, idxs) in bucket_index.range(lo..=hi) {
        for &ci in idxs {
            let cand = &candidates[ci];
            // exact mass gate
            if (cand.peptide.mass() - co.neutral_mass).abs() > tol { continue; }
            let pin = score_psm(&res_ss, &cand.peptide, scorer, z, fragment_tolerance_da);
            let edge = psm_edge_score(&res_ss, &cand.peptide, scorer, z);
            let psm = PsmMatch {
                spectrum_idx: 0, candidate_idxs: vec![ci as u32], charge_used: z,
                mass_error_ppm: ((cand.peptide.mass() - co.neutral_mass) / co.neutral_mass * 1e6),
                score: pin, rank_score: pin + edge as f32, edge_score: edge,
                spec_e_value: 1.0, de_novo_score: i32::MIN,
                activation_method: Some(scorer.param().data_type.activation),
                e_value: 1.0, features: PsmFeatures::default(), isotope_offset: 0,
            };
            queue.push(psm);
        }
    }
    if queue.is_empty() { return None; }
    // 3. one targeted GF SpecEValue on the residual.
    compute_spec_e_values_for_spectrum(
        spec, params, &mut queue, aa_set, enzyme, scorer, &res_ss, z,
        fragment_tolerance_da, search_index, candidates,
    );
    queue.drain_into_vec().into_iter().next()
}
```

Confirm the exact `PsmMatch` field set + `precursor_tolerance.left.as_da` + `nominal_from`/`H2O` against the codebase (mirror how `match_engine.rs` builds a `PsmMatch` and computes the window). Adjust to match.

- [ ] **Step 4: Run to verify pass.** `cargo test -p search --lib coisolation 2>&1 | tail` → all pass. `cargo clippy -p search 2>&1 | tail -3` clean.

- [ ] **Step 5: Commit.**
```bash
git add crates/search/src/coisolation.rs
git commit -m "feat(chimeric): targeted second-peptide residual search (cascade P2)"
```

---

## Task 3: Two-pass driver (binary)

**Files:** Modify `crates/msgf-rust/src/bin/msgf-rust.rs`.

- [ ] **Step 1: After Pass-1 `match_spectra`, add the Pass-2 block.** In `run()`, the search currently produces `queues` (per-spectrum `TopNQueue`) then `write_pin`. When `params.chimeric` and an `Ms1Link` is present (`prepared.ms1_link`), after Pass 1, for each spectrum that has a primary PSM: get its linked MS1 (`link.ms2_to_ms1[spec_idx]` → `link.ms1_peaks[ms1_idx]`), compute the isolation window from `spec.isolation_lower_offset/upper_offset` (fallback `params.chimeric_isolation_halfwidth_da`), call `coisolation::detect_coisolated(...)`, and for each `CoIsolated` call `coisolation::search_secondary(...)` with the primary peptide (the queue's top PSM's candidate peptide). Push any returned secondary `PsmMatch` (fix its `spectrum_idx`) into a SECOND set of per-spectrum queues (or append to a secondary Vec keyed by spec_idx). Concretely, after the `let queues = match_spectra(...)` line:

```rust
    // Pass 2: MS1-gated targeted second-peptide search (chimeric cascade).
    if params.chimeric {
        if let Some(link) = prepared.ms1_link.as_ref() {
            let frag_tol = prepared.fragment_tolerance_da;
            let enz = if params.enzyme != model::enzyme::Enzyme::NoCleavage
                && params.enzyme != model::enzyme::Enzyme::NonSpecific
                { Some(params.enzyme) } else { None };
            for (spec_idx, q) in queues.iter_mut().enumerate() {
                if q.is_empty() { continue; }
                let spec = &spectra[spec_idx];
                let Some(Some(ms1_idx)) = link.ms2_to_ms1.get(spec_idx) else { continue };
                let Some(ms1) = link.ms1_peaks.get(*ms1_idx) else { continue };
                let lo = spec.precursor_mz - spec.isolation_lower_offset.unwrap_or(params.chimeric_isolation_halfwidth_da);
                let hi = spec.precursor_mz + spec.isolation_upper_offset.unwrap_or(params.chimeric_isolation_halfwidth_da);
                let tol = params.precursor_tolerance.left.as_da(spec.precursor_mz).max(0.01);
                let cos = search::coisolation::detect_coisolated(ms1, lo, hi, spec.precursor_mz,
                    *params.charge_range.start()..=*params.charge_range.end(), tol, 1.0, 2);
                if cos.is_empty() { continue; }
                let primary = prepared.candidates[q.iter_psms().next().unwrap().primary_candidate_idx() as usize].peptide.clone();
                for co in cos {
                    if let Some(mut psm) = search::coisolation::search_secondary(
                        spec, &primary, co, &prepared.candidates, &prepared.bucket_index,
                        prepared.scorer, &prepared.aa_set_for_gf, enz, &params, prepared.idx, frag_tol) {
                        psm.spectrum_idx = spec_idx;
                        q.push(psm); // secondary PSM joins the scan's queue (emitted as an extra row)
                    }
                }
            }
        }
    }
```
Note: `coisolation`/`search_secondary`/`detect_coisolated` and the `PreparedSearch` fields they need (`candidates`, `bucket_index`, `scorer`, `aa_set_for_gf`, `idx`, `fragment_tolerance_da`, `ms1_link`) must be `pub`/`pub(crate)` and re-exported as needed — make the minimal visibility changes in `search` to expose them to the binary (or add a thin `search::run_pass2(...)` wrapper in `match_engine.rs` that takes `&PreparedSearch` + `&[Spectrum]` + `&mut [TopNQueue]` and does the loop, keeping the binary thin and the internals `pub(crate)`). PREFER the wrapper: add `pub fn run_pass2_coisolation(prepared: &PreparedSearch, spectra: &[Spectrum], queues: &mut [TopNQueue], params: &SearchParams)` to `match_engine.rs` and call it from the binary. This keeps `coisolation` internals crate-private.

- [ ] **Step 2: Ensure the PIN emits the secondary rows distinctly.** The chimeric SpecId format already appends a per-row index (`pin.rs` `params.chimeric` branch), so a second PSM in the queue gets a unique SpecId. Confirm `write_pin` emits both. No pin.rs change expected; verify.

- [ ] **Step 3: Build + off-path bit-identity.** `cargo build -p msgf-rust 2>&1 | tail`. `cargo test -p search -p output 2>&1 | grep -E 'test result|FAILED'` → all pass except the pre-existing fixture test. `cargo clippy -p search -p output -p msgf-rust 2>&1 | tail`. Off path: `--chimeric off` does not enter the Pass-2 block → narrow PIN bit-identical (diff a baseline off-PIN if in doubt).

- [ ] **Step 4: BSA smoke.**
```bash
cargo build --release -p msgf-rust 2>&1 | tail -2
target/release/msgf-rust --spectrum test-fixtures/test.mgf --database test-fixtures/BSA.fasta --output-pin /tmp/casc.pin --chimeric 2>/tmp/casc.log; echo exit=$?
echo "rows=$(($(wc -l</tmp/casc.pin)-1))"
```
Expected: exit 0; rows ≥ the narrow row count (secondary PSMs add rows on co-isolated scans). (BSA is single-protein with little co-isolation, so few/no secondaries — mainly a no-crash check.)

- [ ] **Step 5: Commit.**
```bash
git add -A
git commit -m "feat(chimeric): two-pass cascade driver (narrow Pass 1 + MS1-gated Pass 2) (cascade P3)"
```

---

## Task 4: VM gates — recall + speed + entrapment

**Files:** none (bench on VM).

- [ ] **Step 1: Ship + rebuild.** scp `coisolation.rs`, `match_engine.rs`, `lib.rs`, `msgf-rust.rs` to `chimeric-build`; `cargo build --release -p msgf-rust`.
- [ ] **Step 2: Astral cascade vs Java + vs blind-chimeric.** Run Astral `--chimeric` (cascade) NO_RESCORE; Percolator @1%; capture wall + MaxRSS. Compare: wall **< Java 6:18**? @1% vs Java 35,818 and vs blind-chimeric 77,287 (what fraction of the +116% did MS1-gating recover?). Record `docs/parity-analysis/notes/2026-05-30-cascade-astral-gate.md`.
- [ ] **Step 3: Entrapment FDP.** Run cascade on `ASTRAL_entrapment.fasta`; Percolator; `compute_entrapment_fdp.py` → secondary PSMs must keep FDP ~nominal.
- [ ] **Step 4: PXD + decision.** Same on PXD. If recall (recovered fraction) is too low, the spec's fallback is all-scans-MS1-localized (a `coisolation.rs` gating change) — note and decide. Final table vs Java (PSMs + wall, 3 datasets). Commit the note; push the branch (parked per `merge-gate-beat-java` until it beats Java on PSMs AND speed).

---

## Self-Review notes (author)

- **Spec coverage:** §1 architecture → T3 driver; §2.1 detection → T1; §2.2 residual + §2.3 targeted scoring → T2; §3 speed model → measured T4; §4 gates → T3 off-bit-identity + T4 recall/speed/entrapment; §5 phases → Tasks 1–4.
- **Type consistency:** `CoIsolated{mono_mz,charge,neutral_mass}`; `detect_coisolated(...)→Vec<CoIsolated>`; `search_secondary(...)→Option<PsmMatch>`; `run_pass2_coisolation` wrapper. Reuses `precursor_isotope_match`, `matched_peak_keys`, `compute_spec_e_values_for_spectrum`, `score_psm`, `psm_edge_score` with their real signatures (verify at execution time — flagged in T1/T2 steps).
- **Off-path:** Pass-2 block gated on `params.chimeric` + `ms1_link.is_some()`; `--chimeric off` unreached → narrow bit-identical.
- **Known verification points (flagged inline):** the exact `PsmMatch` field set, `averagine_isotope_envelope` return type, the nominal↔neutral window (reuse `nominal_from`/`H2O` as in match_engine), and the `match_engine_smoke.rs` fixture for the T2 integration test. These are existing-code lookups, not placeholders.
