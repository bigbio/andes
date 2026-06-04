//! Acceptance gate for incremental model updates.
//!
//! [`evaluate_candidate`] runs the same TDC FDR search used by
//! [`crate::labeled::bootstrap_labels`] with both the *current* and the
//! *candidate* model and returns a [`YieldDelta`] that the caller can use to
//! decide whether to commit the update.
//!
//! # Reuse of search/FDR code
//!
//! Rather than re-implementing the search, we factor the inner search-and-count
//! logic into [`count_target_psms`], which is the same algorithm as
//! `bootstrap_labels` (Steps 1-7 from [`crate::labeled`]) but returns only
//! the count of accepted TARGET PSMs instead of the full `LabeledMatch` list.
//! `bootstrap_labels` itself is left unchanged.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use input::FastaReader;
use model::spectrum::Spectrum;
use scoring_crate::scoring::rank_scorer::RankScorer;
use search::{match_spectra, SearchIndex, SearchParams};

use crate::labeled::BOOTSTRAP_DECOY_PREFIX;
use crate::TrainError;

/// The fragment-mass tolerance used when running the acceptance gate search.
/// Matches the constant in [`crate::labeled`].
const FRAGMENT_TOL_DA: f64 = 0.5;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Number of TARGET PSMs at the requested FDR for the current model and the
/// candidate model.
///
/// A candidate is accepted when `candidate_count >= current_count`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct YieldDelta {
    /// Target PSMs at `fdr` for the current (stored) model.
    pub current_count: usize,
    /// Target PSMs at `fdr` for the candidate model.
    pub candidate_count: usize,
}

impl YieldDelta {
    /// Returns `true` iff the candidate model is at least as good as the
    /// current model (i.e. `candidate_count >= current_count`).
    pub fn is_accepted(&self) -> bool {
        self.candidate_count >= self.current_count
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Evaluate a candidate model against a validation dataset using TDC FDR.
///
/// Runs two searches (one with `current`, one with `candidate`) over
/// `spectra` against `database`, applies the same TDC FDR used in
/// [`crate::labeled::bootstrap_labels`], and returns the target-PSM count at
/// `fdr` for each scorer.
///
/// # Arguments
///
/// * `spectra`       — validation spectra (already loaded).
/// * `database`      — path to the FASTA target database.
/// * `current`       — the currently stored scoring model.
/// * `candidate`     — the proposed replacement model.
/// * `search_params` — search parameters (enzyme, mods, tolerances, …).
/// * `fdr`           — q-value threshold for counting accepted PSMs.
pub fn evaluate_candidate(
    spectra: &[Spectrum],
    database: &Path,
    current: &RankScorer,
    candidate: &RankScorer,
    search_params: &SearchParams,
    fdr: f64,
) -> Result<YieldDelta, TrainError> {
    // Build the SearchIndex once (shared between both searches).
    let file = File::open(database)?;
    let target_db = FastaReader::load_all(BufReader::new(file))
        .map_err(|e| TrainError::Other(format!("FASTA parse error: {e}")))?;
    let idx = SearchIndex::from_target_db(&target_db, BOOTSTRAP_DECOY_PREFIX);

    let current_count = count_target_psms(spectra, &idx, search_params, current, fdr)?;
    let candidate_count = count_target_psms(spectra, &idx, search_params, candidate, fdr)?;

    Ok(YieldDelta { current_count, candidate_count })
}

// ---------------------------------------------------------------------------
// Private helper: search + TDC FDR count
// ---------------------------------------------------------------------------

/// Run a single search and return the number of TARGET PSMs at `fdr`.
///
/// This is the same algorithm as [`crate::labeled::bootstrap_labels`] (Steps
/// 2-7) but returns `usize` instead of `Vec<LabeledMatch>`.
pub(crate) fn count_target_psms(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    search_params: &SearchParams,
    scorer: &RankScorer,
    fdr: f64,
) -> Result<usize, TrainError> {
    // Step 2: run the search.
    let (queues, candidates) = match_spectra(
        spectra,
        idx,
        search_params,
        scorer,
        FRAGMENT_TOL_DA,
        BOOTSTRAP_DECOY_PREFIX,
    );

    // Step 3: collect best PSM per spectrum.
    struct BestPsm {
        rank_score: f32,
        is_decoy: bool,
    }

    let mut best_psms: Vec<BestPsm> = Vec::new();
    for queue in queues.iter() {
        if queue.is_empty() {
            continue;
        }
        if let Some(psm) = queue.peek_top() {
            let cand_idx = psm.primary_candidate_idx() as usize;
            let cand = &candidates[cand_idx];
            best_psms.push(BestPsm {
                rank_score: psm.rank_score,
                is_decoy: cand.is_decoy,
            });
        }
    }

    // Step 4: sort by rank_score DESCENDING (highest RawScore = most confident
    // first), now that the generating function / SpecEValue is removed.
    best_psms.sort_by(|a, b| {
        let av = if a.rank_score.is_nan() { f32::NEG_INFINITY } else { a.rank_score };
        let bv = if b.rank_score.is_nan() { f32::NEG_INFINITY } else { b.rank_score };
        bv.partial_cmp(&av)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.is_decoy.cmp(&b.is_decoy))
    });

    // Step 5: running TDC q-values.
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

    // Step 6: monotone q-values.
    let mut mono_q = raw_q;
    let mut min_q = 1.0_f64;
    for q in mono_q.iter_mut().rev() {
        if *q < min_q {
            min_q = *q;
        }
        *q = min_q;
    }

    // Step 7: count accepted targets.
    let count = best_psms
        .iter()
        .zip(mono_q.iter())
        .filter(|(psm, &q)| !psm.is_decoy && q <= fdr)
        .count();

    Ok(count)
}
