//! `StatsAccumulator` — converts confident (spectrum, peptide) PSMs into
//! [`CountStats`] histograms by reusing the production peak-matching code.
//!
//! # Correctness guarantee
//!
//! The accumulator calls [`scoring_crate::ScoredSpectrum::ion_match_facts`],
//! which reuses **exactly** the same matcher as `directional_node_score_inner`:
//! highest-intensity-in-window peak selection, precursor-peak filtering, global
//! intensity ranking, optional deconvolution, and `param.mme`-based tolerance.
//! No matching logic is duplicated in this crate.

use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::param_model::IonType;
use scoring_crate::scoring::rank_scorer::RankScorer;
use scoring_crate::scoring::scored_spectrum::ScoredSpectrum;

use crate::counts::CountStats;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Accumulates per-ion match statistics from confident PSMs into a
/// [`CountStats`].
///
/// Build one per Rayon worker (it holds an immutable borrow of the `RankScorer`)
/// and merge results with [`merge`].
pub struct StatsAccumulator<'a> {
    scorer: &'a RankScorer,
}

impl<'a> StatsAccumulator<'a> {
    /// Construct a new accumulator backed by `scorer`.
    pub fn new(scorer: &'a RankScorer) -> Self {
        Self { scorer }
    }

    /// Build the `ScoredSpectrum` (applies precursor-peak filtering, global
    /// intensity ranking, and optional deconvolution from the scorer's `Param`)
    /// and accumulate this PSM's ion match facts into `stats`.
    ///
    /// For each theoretical ion of `peptide`:
    ///
    /// - **Matched** (`rank = Some(r)`): calls `stats.bump_rank(partition,
    ///   ion_type, rank_index)` where `rank_index` mirrors the rank-distribution
    ///   array indexing (`rank - 1`, clamped to `[0, max_rank - 1]`), and
    ///   `stats.bump_error(partition, error_bin)` when `error_bin` is `Some`.
    ///
    /// - **Unmatched** (`rank = None`): calls `stats.bump_rank(partition,
    ///   ion_type, max_rank_slot)` — the "missing ion" slot at index `max_rank`
    ///   (the last entry of the rank-distribution array, matching the semantics
    ///   of [`RankScorer::missing_ion_score`]).
    ///
    /// Also increments `stats.bump_charge(charge as i32)` once per PSM.
    pub fn accumulate(
        &self,
        stats: &mut CountStats,
        spec: &Spectrum,
        peptide: &Peptide,
        charge: u8,
    ) {
        // Build ScoredSpectrum via the production constructor: applies precursor
        // filtering, intensity ranking, deconvolution, prob_peak, main_ion, and
        // the segment-partition cache — identical to the search engine path.
        let scored_spec = ScoredSpectrum::new(spec, self.scorer, charge);

        // max_rank is the index of the "missing ion" slot in rank arrays.
        let max_rank = self.scorer.max_rank();

        // Call the production matcher: ion_match_facts reuses directional_node_score_inner's
        // loop body exactly (same partition lookup, same nearest_peak_rank_in, same mme).
        let facts = scored_spec.ion_match_facts(peptide, self.scorer);

        for fact in &facts {
            let rank_idx = match fact.rank {
                // Observed: rank index = (clamped_rank - 1), 0-based.
                // `fact.rank` is already clamped to [1, max_rank] by ion_match_facts.
                Some(r) => {
                    let clamped = r.min(max_rank).max(1);
                    clamped - 1
                }
                // Missing: use max_rank slot (last entry, as per RankScorer::missing_ion_score).
                None => max_rank,
            };
            stats.bump_rank(fact.partition, fact.ion_type, rank_idx);

            // Error bin: only for matched ions when error_scaling_factor > 0.
            if let (Some(_), Some(bin)) = (fact.rank, fact.error_bin) {
                stats.bump_error(fact.partition, bin);
            }
        }

        // Noise rank distribution: probe background m/z positions (same matcher,
        // same density as the ions) and record their ranks under IonType::Noise.
        // RankScorer needs a Noise entry per partition, and its SHAPE (dominated
        // by the "missing" slot) calibrates the ion-vs-noise likelihood ratio.
        // Without this the missing-ion penalty inverts and the model scores
        // target and decoy alike (0 PSMs at 1% FDR).
        for (partition, rank) in scored_spec.noise_match_facts(peptide, self.scorer) {
            let rank_idx = match rank {
                Some(r) => r.min(max_rank).max(1) - 1,
                None => max_rank,
            };
            stats.bump_rank(partition, IonType::Noise, rank_idx);
        }

        // Charge histogram: one bump per PSM.
        stats.bump_charge(charge as i32);
    }
}

/// Merge a collection of independently accumulated `CountStats` into one.
///
/// Designed for Rayon parallel accumulation: each worker fills its own
/// `CountStats`, then the results are merged with this function.
///
/// Returns an empty `CountStats` if `parts` is empty.
pub fn merge(parts: Vec<CountStats>) -> CountStats {
    parts.into_iter().fold(CountStats::new(), |mut acc, part| {
        acc.add(&part);
        acc
    })
}
