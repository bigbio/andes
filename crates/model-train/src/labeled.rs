//! Bootstrap label adapter: run a seed search, compute target-decoy FDR, and
//! return the confident target PSMs as [`LabeledMatch`]es for the trainer.
//!
//! # Algorithm
//!
//! 1. Load the FASTA target database and build a [`SearchIndex`] (targets +
//!    reversed decoys), reusing the same construction as the production binary.
//! 2. Run the full search with the seed `RankScorer` via [`match_spectra`] —
//!    no search or scoring logic is re-implemented here.
//! 3. For each spectrum take the **best PSM** (lowest `spec_e_value`, then
//!    highest `score`) from the queue.
//! 4. Sort best-per-spectrum PSMs by `spec_e_value` ascending (best first).
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

/// Fragment-mass tolerance (Da) passed to the scoring engine.
/// Matches the production binary's default for HCD.
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
    // Each TopNQueue is already sorted best-first by `spec_e_value`. We take
    // the single best entry per spectrum.
    struct BestPsm {
        spectrum_index: usize,
        peptide: Peptide,
        charge: u8,
        spec_e_value: f64,
        is_decoy: bool,
    }

    let mut best_psms: Vec<BestPsm> = Vec::new();
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        // peek_top returns the best (lowest spec_e_value, then highest score).
        if let Some(psm) = queue.peek_top() {
            let cand_idx = psm.primary_candidate_idx() as usize;
            let cand = &candidates[cand_idx];
            best_psms.push(BestPsm {
                spectrum_index: spec_idx,
                peptide: cand.peptide.clone(),
                charge: psm.charge_used,
                spec_e_value: psm.spec_e_value,
                is_decoy: cand.is_decoy,
            });
        }
    }

    // ── 4. Sort by spec_e_value ascending (best/most-confident first) ─────────
    // Tie-break by is_decoy ascending (target before decoy) for stability.
    best_psms.sort_by(|a, b| {
        let av = if a.spec_e_value.is_nan() { f64::INFINITY } else { a.spec_e_value };
        let bv = if b.spec_e_value.is_nan() { f64::INFINITY } else { b.spec_e_value };
        av.partial_cmp(&bv)
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
