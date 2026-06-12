//! Bootstrap label adapter: run a seed search, compute target-decoy FDR, and
//! return the confident target PSMs as [`LabeledMatch`]es for the trainer.
//!
//! # Algorithm
//!
//! 1. Load the FASTA target database and build a [`SearchIndex`] (targets +
//!    reversed decoys), reusing the same construction as the production binary.
//! 2. Run the full search with the seed `RankScorer` via [`match_spectra`] —
//!    no search or scoring logic is re-implemented here.
//! 3. For each spectrum take the **best PSM** (highest `rank_score`) from the
//!    queue.
//! 4. Sort best-per-spectrum PSMs by `rank_score` descending (best first).
//! 5. Walk the sorted list computing a running target/decoy count and q-value
//!    at each position: `q = running_decoys / max(running_targets, 1)`.
//! 6. Convert to a proper monotone q-value by scanning from the bottom and
//!    propagating the minimum seen so far.
//! 7. Emit a [`LabeledMatch`] for every TARGET PSM whose q-value ≤ `train_fdr`.
//!
//! # Determinism
//!
//! The search engine (`match_spectra`) is deterministic for fixed inputs and a
//! fixed Rayon thread count. The Rayon thread pool is shared across the
//! process; the first caller establishes the pool size. For bit-identical
//! results across runs, the caller should ensure a fixed number of threads
//! (e.g. via `RAYON_NUM_THREADS`). The FDR computation and the returned Vec
//! are fully deterministic after the search.
//!
//! # Decoy prefix
//!
//! The default decoy prefix `"XXX_"` is used. Reversed decoys receive
//! accessions like `"XXX_P02769"`, and `Candidate::is_decoy` is set by
//! `starts_with("XXX_")` in `enumerate_candidates`. The prefix is consistent
//! with the production binary default.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use input::FastaReader;
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::scoring::rank_scorer::RankScorer;
use search::{match_spectra, SearchIndex, SearchParams};

use crate::TrainError;

/// The decoy prefix used by `bootstrap_labels` when building the reversed
/// decoy database and classifying candidates as target or decoy.
pub const BOOTSTRAP_DECOY_PREFIX: &str = "XXX_";

/// Fragment-mass tolerance (Da) for the GF node-scoring match window.
///
/// Intentionally 0.5 Da regardless of instrument, and MUST stay equal to the
/// production search binary's `fragment_tol_da` (also 0.5). The generating
/// function scores over an integer nominal-mass axis (~0.5 Da bins), so the
/// node-match window cannot be tighter than ~0.5 Da — accurate-mass (ppm)
/// discrimination is recaptured in the Percolator features, NOT here. Deriving
/// this from `seed.param().mme` (e.g. 20 ppm) would make training diverge from
/// the production search and is incorrect.
const FRAGMENT_TOL_DA: f64 = 0.5;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A confident target PSM accepted from the bootstrap seed search.
///
/// Carries everything the [`crate::accumulate::StatsAccumulator`] needs:
/// `spectrum`, `peptide`, and `charge`. The spectrum itself is not stored
/// here — the caller supplies the `spectra` slice to the trainer and looks
/// up `spectra[spectrum_index]`.
#[derive(Debug, Clone)]
pub struct LabeledMatch {
    /// Index into the `spectra` slice passed to [`bootstrap_labels`].
    pub spectrum_index: usize,
    /// The accepted (target) peptide for this spectrum.
    pub peptide: Peptide,
    /// The precursor charge state used when scoring.
    pub charge: u8,
    /// The q-value (monotone TDC FDR) at which this PSM was accepted.
    /// Always ≤ the `train_fdr` passed to [`bootstrap_labels`].
    pub confidence: f64,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Run a seed search over `spectra` against the FASTA `database`, generate
/// reversed decoys automatically, compute target-decoy competition FDR over
/// the best PSM per spectrum, and return the TARGET PSMs whose q-value ≤
/// `train_fdr`.
///
/// # Arguments
///
/// * `spectra`       — the MS2 spectra to search (already loaded by the caller).
/// * `database`      — path to the FASTA target database (decoys generated internally).
/// * `seed`          — the scoring model used for the seed search.
/// * `search_params` — search parameters (enzyme, mods, tolerances, …).
/// * `train_fdr`     — q-value threshold for accepting a PSM as a confident label.
///
/// # Returns
///
/// `Ok(Vec<LabeledMatch>)` — possibly empty when `train_fdr` is very strict or
/// the fixture is tiny with no confident identifications.
///
/// # Errors
///
/// Returns [`TrainError::Io`] if the FASTA cannot be opened or read, or
/// [`TrainError::Other`] for downstream parsing failures.
pub fn bootstrap_labels(
    spectra: &[Spectrum],
    database: &Path,
    seed: &RankScorer,
    search_params: &SearchParams,
    train_fdr: f64,
) -> Result<Vec<LabeledMatch>, TrainError> {
    // ── 1. Load FASTA + build SearchIndex (target + reversed decoys) ─────────
    let file = File::open(database)?;
    let target_db = FastaReader::load_all(BufReader::new(file))
        .map_err(|e| TrainError::Other(format!("FASTA parse error: {e}")))?;
    let idx = SearchIndex::from_target_db(&target_db, BOOTSTRAP_DECOY_PREFIX);

    // ── 2. Run the seed search ────────────────────────────────────────────────
    // Reuse the production `match_spectra` entry point — no search logic
    // is duplicated here.
    let (queues, candidates) = match_spectra(
        spectra,
        &idx,
        search_params,
        seed,
        FRAGMENT_TOL_DA,
        BOOTSTRAP_DECOY_PREFIX,
    );

    // ── 3. Collect best PSM per spectrum ──────────────────────────────────────
    // Each TopNQueue is already sorted best-first by `rank_score`. We take
    // the single best entry per spectrum.
    struct BestPsm {
        spectrum_index: usize,
        peptide: Peptide,
        charge: u8,
        rank_score: f32,
        is_decoy: bool,
    }

    let mut best_psms: Vec<BestPsm> = Vec::new();
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        // peek_top returns the best (highest rank_score).
        if let Some(psm) = queue.peek_top() {
            let cand_idx = psm.primary_candidate_idx() as usize;
            let cand = &candidates[cand_idx];
            best_psms.push(BestPsm {
                spectrum_index: spec_idx,
                peptide: cand.peptide.clone(),
                charge: psm.charge_used,
                rank_score: psm.rank_score,
                is_decoy: cand.is_decoy,
            });
        }
    }

    // ── 4. Sort by rank_score DESCENDING (best/most-confident first) ──────────
    // Tie-break by is_decoy ascending (target before decoy) for stability.
    best_psms.sort_by(|a, b| {
        let av = if a.rank_score.is_nan() { f32::NEG_INFINITY } else { a.rank_score };
        let bv = if b.rank_score.is_nan() { f32::NEG_INFINITY } else { b.rank_score };
        bv.partial_cmp(&av)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.is_decoy.cmp(&b.is_decoy))
    });

    // ── 5. Compute running q-values via TDC ──────────────────────────────────
    // Walk from best to worst, accumulating target/decoy counts.
    // q-value at position i = running_decoys / max(running_targets, 1).
    let n = best_psms.len();
    let mut raw_q = vec![1.0_f64; n];
    let mut n_targets = 0u64;
    let mut n_decoys = 0u64;
    for (i, psm) in best_psms.iter().enumerate() {
        if psm.is_decoy {
            n_decoys += 1;
        } else {
            n_targets += 1;
        }
        raw_q[i] = n_decoys as f64 / n_targets.max(1) as f64;
    }

    // ── 6. Monotone q-value: propagate minimum from the bottom ────────────────
    // Scan from worst (end) to best (start), ensuring q[i] <= q[i+1].
    let mut mono_q = raw_q;
    let mut min_q = 1.0_f64;
    for q in mono_q.iter_mut().rev() {
        if *q < min_q {
            min_q = *q;
        }
        *q = min_q;
    }

    // ── 6b. Conservative tie handling for equal `rank_score` buckets ──────────
    // Raw `rank_score` (RawScore) is coarse, so ties are common. The sort above
    // is not guaranteed to order target-vs-decoy within a tie favourably, and
    // target-before-decoy ordering would deflate the running decoy ratio inside
    // the bucket, accepting a target at the requested FDR even when a tied decoy
    // should fail it. Assign every PSM in an equal-`rank_score` bucket the WORST
    // (max) q-value seen anywhere in that bucket so the whole bucket passes or
    // fails together (pessimistic). `best_psms` is sorted by `rank_score`
    // descending, so equal scores are contiguous.
    assign_bucket_worst_q(&best_psms, &mut mono_q, |p| p.rank_score);

    // ── 7. Collect accepted TARGET PSMs ──────────────────────────────────────
    let mut labels: Vec<LabeledMatch> = Vec::new();
    for (psm, q) in best_psms.into_iter().zip(mono_q.into_iter()) {
        if !psm.is_decoy && q <= train_fdr {
            labels.push(LabeledMatch {
                spectrum_index: psm.spectrum_index,
                peptide: psm.peptide,
                charge: psm.charge,
                confidence: q,
            });
        }
    }

    Ok(labels)
}

/// Make q-values conservative across equal-`rank_score` ties.
///
/// `psms` MUST be sorted by `rank_score` descending (equal scores contiguous),
/// matching the order `q` was computed in. For each maximal run of PSMs with
/// the same `rank_score` (compared on exact `f32` bits, NaN-safe), every q in
/// that run is set to the WORST (maximum) q seen in the run. This prevents a
/// favourable within-tie ordering (e.g. target-before-decoy) from accepting a
/// target at a given FDR when a tied decoy should make the whole bucket fail.
///
/// Shared by [`bootstrap_labels`] and [`crate::gate::count_target_psms`].
pub(crate) fn assign_bucket_worst_q<T>(
    psms: &[T],
    q: &mut [f64],
    score_of: impl Fn(&T) -> f32,
) {
    debug_assert_eq!(psms.len(), q.len());
    let mut start = 0usize;
    while start < psms.len() {
        let s = score_of(&psms[start]);
        let mut end = start + 1;
        while end < psms.len() && scores_tie(score_of(&psms[end]), s) {
            end += 1;
        }
        // Worst (max) q within [start, end).
        let mut worst = q[start];
        for &qi in &q[start + 1..end] {
            if qi > worst {
                worst = qi;
            }
        }
        for qi in &mut q[start..end] {
            *qi = worst;
        }
        start = end;
    }
}

/// Two `rank_score`s tie iff they have identical bit patterns. Treats all NaNs
/// as a single tie class (they sort together as `NEG_INFINITY` above).
fn scores_tie(a: f32, b: f32) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    a.to_bits() == b.to_bits()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal stand-in for the per-spectrum best PSM used by the FDR walk.
    struct P {
        rank_score: f32,
        is_decoy: bool,
    }

    /// Reproduce Steps 4-6b of the FDR computation on a fixture and return the
    /// per-PSM q-values plus the number of accepted targets at `fdr`.
    fn fdr(mut psms: Vec<P>, fdr: f64) -> (Vec<f64>, usize) {
        // Step 4: target-before-decoy tie-break (the original, "optimistic"
        // ordering) — the conservative bucket pass must neutralise its bias.
        psms.sort_by(|a, b| {
            let av = if a.rank_score.is_nan() { f32::NEG_INFINITY } else { a.rank_score };
            let bv = if b.rank_score.is_nan() { f32::NEG_INFINITY } else { b.rank_score };
            bv.partial_cmp(&av)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.is_decoy.cmp(&b.is_decoy))
        });
        let n = psms.len();
        let mut q = vec![1.0_f64; n];
        let (mut t, mut d) = (0u64, 0u64);
        for (i, p) in psms.iter().enumerate() {
            if p.is_decoy { d += 1 } else { t += 1 }
            q[i] = d as f64 / t.max(1) as f64;
        }
        let mut min_q = 1.0;
        for qi in q.iter_mut().rev() {
            if *qi < min_q { min_q = *qi; }
            *qi = min_q;
        }
        assign_bucket_worst_q(&psms, &mut q, |p| p.rank_score);
        let accepted = psms
            .iter()
            .zip(q.iter())
            .filter(|(p, &qi)| !p.is_decoy && qi <= fdr)
            .count();
        (q, accepted)
    }

    #[test]
    fn tied_target_and_decoy_bucket_is_counted_conservatively() {
        // One confident target, then a target tied with a decoy at the same
        // raw rank_score. With target-before-decoy ordering the running ratio
        // would dip *before* the decoy is counted, accepting the tied target at
        // 1% FDR. The conservative bucket pass must give both tied PSMs the same
        // (worst) q, so the tied decoy makes the bucket fail and the tied target
        // is NOT accepted.
        let psms = vec![
            P { rank_score: 10.0, is_decoy: false }, // clearly good target
            P { rank_score: 5.0, is_decoy: false },  // tied with decoy below
            P { rank_score: 5.0, is_decoy: true },   // tied decoy
        ];
        let (q, accepted) = fdr(psms, 0.01);
        // The two tied PSMs share the worst q in their bucket (q = 0.5 here),
        // so neither passes 1% FDR. Only the first, clearly-good target does.
        assert_eq!(accepted, 1, "tied target must not be accepted when its bucket has a decoy");
        assert_eq!(q[1], q[2], "tied PSMs must share a q-value");
        assert!(q[1] > 0.01, "tied bucket q must reflect the decoy");
    }

    #[test]
    fn untied_scores_are_unaffected() {
        // With strictly distinct scores the bucket pass is a no-op: each PSM is
        // its own bucket, so q-values match the plain monotone TDC result.
        let psms = vec![
            P { rank_score: 10.0, is_decoy: false },
            P { rank_score: 9.0, is_decoy: false },
            P { rank_score: 8.0, is_decoy: true },
        ];
        let (q, accepted) = fdr(psms, 0.01);
        assert_eq!(accepted, 2, "both distinct-score targets pass at 1% FDR");
        assert!(q[0] <= 0.01 && q[1] <= 0.01);
    }
}
