//! Per-ion intensity-rank log-likelihood-ratio (LLR) scoring.
//!
//! Implements the published intensity-rank node-scoring method (Frank,
//! *J. Proteome Res.* 2005; PMC2738854): peaks are ranked by intensity, and a
//! matched fragment ion observed at a given rank contributes a log-likelihood
//! ratio comparing how often a true ion of that type lands at that rank
//! (`ion_freq`) against how often a random noise peak does (`noise_freq`):
//!
//! ```text
//!   norm        = min(ion_charge, num_segments)
//!   llr[rank]    = ln( ion_freq[rank] / (noise_freq[rank] * norm) )
//! ```
//!
//! The `norm` term (denoted `charge_or_seg` below) caps the ion's charge at
//! the model's segment count; it normalizes the noise expectation so that
//! multiply-charged ions, which can land in several m/z segments, are not
//! over-credited. It is a model-derived normalization, not a tunable.
//!
//! The per-rank frequency tables come from the trained model and have length
//! `max_rank + 1`. The first `max_rank` entries (indices `0..max_rank`) hold
//! the LLR for observed ranks `1..=max_rank`; the final entry (index
//! `max_rank`) is a separate "absent" slot scoring an expected ion that is not
//! observed in the spectrum (see [`RankScorer::missing_ion_score`]).

use std::collections::HashMap;

use crate::param_model::{IonType, Param, Partition};

#[derive(Debug, Clone)]
pub struct RankScorer {
    /// The model this scorer was precomputed from. Held so callers (e.g.
    /// `match_engine`) can reach precursor-filter settings without threading a
    /// separate `Param` argument through `ScoredSpectrum::new`.
    param: Param,
    /// Precomputed LLR tables keyed by `(partition, ion_type)`, each of length
    /// `max_rank + 1` (entries `0..max_rank` are the observed-rank scores for
    /// ranks `1..=max_rank`; entry `max_rank` is the absent-ion score). Noise
    /// is never a key — it only supplies the denominator. This map backs the
    /// per-ion query API (`node_score`, `missing_ion_score`) used by tests and
    /// diagnostics; the DP hot path indexes `partition_ion_logs` instead.
    pub(crate) log_table: HashMap<(Partition, IonType), Vec<f32>>,
    /// Same LLR tables flattened per partition into a dense list of
    /// `(ion_type, table)` pairs, ordered to match
    /// `Param::ion_types_for_partition_slice`. The DP loop holds a partition
    /// fixed across many node evaluations, so iterating this slice replaces a
    /// `(partition, ion)` hash lookup per ion with a plain array walk — on a
    /// PXD001819 search that removes on the order of 200M map probes.
    pub(crate) partition_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>>,
    /// Per-partition LLR tables for **neutral-loss** ion types only
    /// (`IonType::is_loss()`), keyed exactly like `partition_ion_logs`. The
    /// mass-indexed DP hot path (`partition_ion_logs`) is intentionally
    /// intact-only because a loss ion's m/z shift (`−L/z`) is peptide-specific
    /// (it comes from the matched peptide's mod, not the model), so loss ions
    /// must be scored in a peptide-aware pass (`ScoredSpectrum::loss_node_score`)
    /// rather than the peptide-agnostic cache. Empty for every model without a
    /// trained loss table (all 39 bundled models) ⇒ no-loss scoring is
    /// byte-identical to the pre-feature engine.
    pub(crate) partition_loss_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>>,
    /// Highest observed rank the model distinguishes; ranks beyond this are
    /// clamped down to it before indexing.
    max_rank: u32,
    /// Optional CLI override for the fragment-matching tolerance (MGF input).
    /// `None` ⇒ derive from the model's instrument resolution class (the
    /// historical 20 ppm high-res / 0.5 Da low-res default).
    fragment_tol_override: Option<model::tolerance::Tolerance>,
}

impl RankScorer {
    pub fn new(param: &Param) -> Self {
        // Precompute, for every (partition, ion) pair, the per-rank LLR table
        // ln(ion_freq / (noise_freq * norm)). The noise frequencies are the
        // shared denominator for all ions within a partition.
        let mut log_table: HashMap<(Partition, IonType), Vec<f32>> = HashMap::new();

        for (partition, ion_table) in &param.rank_dist_table {
            // A partition without a noise distribution has no denominator, so
            // nothing in it can be scored — skip it entirely.
            let Some(noise_freqs) = ion_table.get(&IonType::Noise) else {
                continue;
            };

            for (ion_type, ion_freqs) in ion_table {
                let ion_charge = match ion_type {
                    IonType::Prefix { charge, .. } | IonType::Suffix { charge, .. } => *charge,
                    // Noise is the denominator, not a scored ion.
                    IonType::Noise => continue,
                };

                // Cap the ion charge at the model's segment count (see module
                // docs): a charge-z ion can fall into up to `num_segments`
                // segments, so its noise expectation is scaled accordingly.
                let charge_or_seg = (ion_charge as u32).min(param.num_segments as u32) as f32;

                // Tables may differ in length across releases; the shared
                // prefix is all that is jointly defined.
                let rank_count = ion_freqs.len().min(noise_freqs.len());
                let table: Vec<f32> = (0..rank_count)
                    .map(|r| (ion_freqs[r] / (noise_freqs[r] * charge_or_seg)).ln())
                    .collect();

                log_table.insert((*partition, *ion_type), table);
            }
        }

        // Flatten the per-pair tables into per-partition dense lists so the DP
        // hot path can iterate by index. Ion order mirrors
        // `Param::ion_types_for_partition_slice`; ions with no LLR table
        // (e.g. missing a noise denominator) are simply absent here.
        let mut partition_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>> = HashMap::new();
        for (&partition, ions) in &param.partition_ion_types_cache {
            let paired: Vec<(IonType, Vec<f32>)> = ions
                .iter()
                // INTACT ONLY: loss ions are excluded from the mass-indexed hot
                // path (their m/z shift is peptide-specific). This filter makes
                // the hot path byte-identical even for a model that does carry
                // loss ion types in its frag-offset cache.
                .filter(|ion| !ion.is_loss())
                .filter_map(|&ion| {
                    log_table
                        .get(&(partition, ion))
                        .map(|table| (ion, table.clone()))
                })
                .collect();
            partition_ion_logs.insert(partition, paired);
        }

        // Loss-ion tables, gathered directly from `log_table` (i.e. from
        // `rank_dist_table`) so a trained loss table is visible to the
        // peptide-aware loss pass even when no `FragmentOffsetFrequency` entry
        // exists for it. Keyed by partition; values are the loss `IonType`s with
        // their per-rank LLR tables.
        let mut partition_loss_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>> =
            HashMap::new();
        for (&(partition, ion), table) in &log_table {
            if ion.is_loss() {
                partition_loss_ion_logs
                    .entry(partition)
                    .or_default()
                    .push((ion, table.clone()));
            }
        }

        Self {
            param: param.clone(),
            log_table,
            partition_ion_logs,
            partition_loss_ion_logs,
            max_rank: param.max_rank as u32,
            fragment_tol_override: None,
        }
    }

    /// Borrow the dense `(IonType, log_table)` pairs for `partition`. Used by
    /// the GF DP hot path so per-ion scoring is array indexing, not HashMap
    /// lookup. Returns empty slice if the partition has no ions.
    pub fn partition_ion_logs(&self, partition: &Partition) -> &[(IonType, Vec<f32>)] {
        self.partition_ion_logs
            .get(partition)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Borrow the `(IonType, log_table)` pairs for **neutral-loss** ion types in
    /// `partition`. Used by the peptide-aware loss-scoring pass
    /// (`ScoredSpectrum::loss_node_score`). Returns an empty slice for any model
    /// without a trained loss table — i.e. for every standard search.
    pub fn partition_loss_ion_logs(&self, partition: &Partition) -> &[(IonType, Vec<f32>)] {
        self.partition_loss_ion_logs
            .get(partition)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// True if this model carries any trained neutral-loss ion table. When
    /// false, the peptide-aware loss pass is a no-op (byte-identical scoring).
    pub fn has_loss_tables(&self) -> bool {
        !self.partition_loss_ion_logs.is_empty()
    }

    /// Maximum rank used for clamping. Exposed so callers can apply
    /// rank-clamp / missing-ion semantics without going through `node_score`.
    pub fn max_rank(&self) -> u32 {
        self.max_rank
    }

    /// Return the `Param` this scorer was built from.
    pub fn param(&self) -> &Param {
        &self.param
    }

    /// Set (or clear) the CLI fragment-tolerance override. Used by the search
    /// binary for metadata-less (MGF) input; metadata-bearing input leaves it
    /// `None` so behavior is byte-identical to auto-detection.
    pub fn set_fragment_tol_override(&mut self, tol: Option<model::tolerance::Tolerance>) {
        self.fragment_tol_override = tol;
    }

    /// Effective fragment-matching tolerance for Percolator-feature and
    /// peak-match counting. Returns the CLI override when set, else the
    /// instrument-derived default: 20 ppm for high-resolution analyzers
    /// (Kim et al., Nat Commun 5:5277, 2014), 0.5 Da for ion-trap low-res.
    pub fn feature_match_tolerance(&self) -> model::tolerance::Tolerance {
        use model::tolerance::Tolerance;
        self.fragment_tol_override.unwrap_or_else(|| {
            if self.param.data_type.instrument.is_high_resolution() {
                Tolerance::Ppm(20.0)
            } else {
                Tolerance::Da(0.5)
            }
        })
    }

    /// LLR for a matched ion observed at intensity `rank` (1-based; rank 1 is
    /// the most intense peak). Ranks beyond `max_rank` are folded onto
    /// `max_rank`, so they index the last *observed-rank* entry
    /// (`max_rank - 1`) — never the absent-ion slot. Returns 0 (neutral) for an
    /// ion/partition the model doesn't cover.
    pub fn node_score(&self, partition: Partition, ion_type: IonType, rank: u32) -> f32 {
        let Some(table) = self.log_table.get(&(partition, ion_type)) else {
            return 0.0;
        };
        // Panic-safety: `clamp` requires min <= max. A degenerate model with
        // `max_rank == 0` would make `clamp(1, 0)` panic. Real models always
        // have `max_rank > 0`, so this is inert for any valid model.
        if self.max_rank == 0 {
            return 0.0;
        }
        let observed_rank = rank.clamp(1, self.max_rank);
        table
            .get((observed_rank - 1) as usize)
            .copied()
            .unwrap_or(0.0)
    }

    /// LLR penalty for an expected ion that does not appear in the spectrum.
    /// This is the dedicated absent slot at index `max_rank`, the final entry
    /// of the `max_rank + 1`-length table. Returns 0 if the ion/partition is
    /// uncovered.
    pub fn missing_ion_score(&self, partition: Partition, ion_type: IonType) -> f32 {
        let Some(table) = self.log_table.get(&(partition, ion_type)) else {
            return 0.0;
        };
        table.get(self.max_rank as usize).copied().unwrap_or(0.0)
    }

    /// Ion-pair existence LLR.
    ///
    /// Scores whether the prefix/suffix peak pair flanking a cleavage site is
    /// present, by comparing the model's learned existence probability for the
    /// pair against the probability that random noise would produce the same
    /// presence/absence pattern. With `prob_peak` the per-position chance of a
    /// noise peak, the four patterns indexed `0..=3` get the independent-noise
    /// baseline:
    ///
    /// - index 0 (neither present): `(1 - prob_peak)^2`
    /// - index 1, 2 (exactly one present): `prob_peak * (1 - prob_peak)`
    /// - index 3 (both present): `prob_peak^2`
    ///
    /// Returns `ln(ion_existence_prob / noise_baseline)`. Yields 0 if the
    /// partition has no existence table, or `index` is out of range.
    ///
    /// Degenerate-input behavior is load-bearing and intentionally left
    /// un-clamped on the denominator: for very peak-dense spectra at small
    /// parent mass, `prob_peak` can exceed 1, which makes the one-present
    /// baseline (`prob_peak * (1 - prob_peak)`) negative. `ln` of a
    /// positive/negative ratio is NaN, and the caller's `round() as i32` maps
    /// NaN to 0 — neutralizing that edge. Clamping the denominator to a tiny
    /// positive value instead would emit a large spurious positive score
    /// (e.g. `ln(0.028 / 1e-38) ≈ +84`) per affected edge and inflate DP
    /// maxima by roughly an order of magnitude on short charge-2 peptides, so
    /// we let NaN/±inf flow through to the rounding step unchanged.
    pub fn ion_existence_score(&self, partition: Partition, index: usize, prob_peak: f32) -> f32 {
        let Some(table) = self.param.ion_existence_table.get(&partition) else {
            return 0.0;
        };
        if index >= table.len() {
            return 0.0;
        }
        let noise_baseline = match index {
            0 => (1.0 - prob_peak) * (1.0 - prob_peak),
            3 => prob_peak * prob_peak,
            _ => prob_peak * (1.0 - prob_peak),
        };
        // Floor an exact-zero learned probability to 0.01 so ln stays finite;
        // a true 0 would otherwise force ln(0) = -inf for an observed pair.
        let ion_prob = if table[index] == 0.0 { 0.01 } else { table[index] };
        // Deliberately no denominator clamp (see the doc note): NaN/±inf are
        // expected on degenerate input and are resolved by the caller's round.
        (ion_prob / noise_baseline).ln()
    }

    /// Mass-error score.
    ///
    /// Converts `error` (in Da) to an index using `error_scaling_factor`,
    /// clamps to `[-esf, esf]`, then returns
    /// `log(ionErrHist[idx] / noiseErrHist[idx])`.
    ///
    /// Returns 0.0 if `error_scaling_factor == 0` or tables are missing.
    pub fn error_score(&self, partition: Partition, error: f32) -> f32 {
        let esf = self.param.error_scaling_factor;
        if esf == 0 {
            return 0.0;
        }
        let mut err_index = (error * esf as f32).round() as i32;
        if err_index > esf { err_index = esf; }
        else if err_index < -esf { err_index = -esf; }
        err_index += esf;
        let idx = err_index as usize;

        let ion_err = match self.param.ion_err_dist_table.get(&partition) {
            Some(v) => v,
            None => return 0.0,
        };
        let noise_err = match self.param.noise_err_dist_table.get(&partition) {
            Some(v) => v,
            None => return 0.0,
        };
        if idx >= ion_err.len() || idx >= noise_err.len() {
            return 0.0;
        }
        let ion_f = ion_err[idx];
        let noise_f = noise_err[idx];
        if ion_f <= 0.0 || noise_f <= 0.0 {
            return 0.0;
        }
        (ion_f / noise_f).ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::tiny_param;

    /// Build a minimal low-resolution `Param` (LowRes instrument) for
    /// `feature_match_tolerance` tests. Reuses `tiny_param()` and patches
    /// `data_type.instrument` to `InstrumentType::LowRes`.
    fn low_res_param() -> Param {
        use model::instrument::InstrumentType;
        let mut p = tiny_param();
        p.data_type.instrument = InstrumentType::LowRes;
        p
    }

    #[test]
    fn feature_match_tolerance_defaults_by_resolution() {
        use model::tolerance::Tolerance;
        // tiny_param() uses QExactive → high-resolution → 20 ppm
        let hi = RankScorer::new(&tiny_param());
        assert_eq!(hi.feature_match_tolerance(), Tolerance::Ppm(20.0));
        // low_res_param() uses LowRes → 0.5 Da
        let lo = RankScorer::new(&low_res_param());
        assert_eq!(lo.feature_match_tolerance(), Tolerance::Da(0.5));
    }

    #[test]
    fn feature_match_tolerance_honors_override() {
        use model::tolerance::Tolerance;
        let mut s = RankScorer::new(&low_res_param());
        s.set_fragment_tol_override(Some(Tolerance::Da(0.6)));
        assert_eq!(s.feature_match_tolerance(), Tolerance::Da(0.6));
        s.set_fragment_tol_override(None);
        assert_eq!(s.feature_match_tolerance(), Tolerance::Da(0.5));
    }

    #[test]
    fn node_score_log_formula() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };

        // Rank 1 → index 0. chargeOrSeg = min(1, 1) = 1. log(0.6 / (0.1 * 1)) = log(6.0).
        let s1 = scorer.node_score(part, ion, 1);
        assert!((s1 - 6.0_f32.ln()).abs() < 1e-5, "rank1: got {s1}, expected {}", 6.0_f32.ln());

        // Rank 2 → index 1. log(0.3 / 0.2) = log(1.5).
        let s2 = scorer.node_score(part, ion, 2);
        assert!((s2 - 1.5_f32.ln()).abs() < 1e-5);
    }

    #[test]
    fn rank_above_max_clamps() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };

        // rank > max_rank clamps to rank_index = max_rank - 1.
        // max_rank = 3 → rank_index = 2 → log(0.05 / 0.3).
        let s5 = scorer.node_score(part, ion, 5);
        let expected = (0.05_f32 / 0.3_f32).ln();
        assert!((s5 - expected).abs() < 1e-5);
    }

    #[test]
    fn missing_ion_score_uses_last_slot() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };

        // missing slot = index `maxRank` = 3 (the last entry in length-4 array).
        // log(0.001 / 0.4) = log(0.0025).
        let s_missing = scorer.missing_ion_score(part, ion);
        let expected = (0.001_f32 / 0.4_f32).ln();
        assert!((s_missing - expected).abs() < 1e-5);
    }

    #[test]
    fn chargeorseg_uses_min_of_ion_charge_and_num_segments() {
        // Build a param with num_segments=1 but an ion with charge 3.
        // charge_or_seg = min(3, 1) = 1.
        // Verify the log score uses 1 (not 3).
        let mut param = tiny_param();
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let ion3 = IonType::Prefix { charge: 3, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
        param.rank_dist_table.get_mut(&part).unwrap().insert(ion3, ion_freqs);

        let scorer = RankScorer::new(&param);
        let s1 = scorer.node_score(part, ion3, 1);
        // charge_or_seg = min(3, 1) = 1. log(0.6 / (0.1 * 1)) = log(6).
        assert!((s1 - 6.0_f32.ln()).abs() < 1e-5);
    }

    #[test]
    fn unknown_partition_returns_zero() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let unknown = Partition { charge: 99, parent_mass: 0.0, seg_num: 0 };
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        // Out-of-table partition → return 0 (neutral score).
        assert_eq!(scorer.node_score(unknown, ion, 1), 0.0);
        assert_eq!(scorer.missing_ion_score(unknown, ion), 0.0);
    }

    #[test]
    fn loss_ions_routed_to_loss_logs_not_intact_path() {
        // A model carrying a loss-class ion type must expose it ONLY through
        // partition_loss_ion_logs; the intact hot path must filter it out so
        // the mass-indexed DP stays peptide-agnostic and intact-only.
        let mut param = tiny_param();
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let intact = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let loss = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 1 };
        // Give the loss ion a rank distribution (so it gets an LLR table).
        param
            .rank_dist_table
            .get_mut(&part)
            .unwrap()
            .insert(loss, vec![0.6_f32, 0.3, 0.05, 0.001]);
        // Populate the ion-type cache with BOTH so we can prove the intact path
        // filters the loss ion out (an empty cache would pass trivially).
        param.partition_ion_types_cache.insert(part, vec![intact, loss]);

        let scorer = RankScorer::new(&param);

        // Hot path: intact only.
        assert!(scorer.partition_ion_logs(&part).iter().all(|(ion, _)| !ion.is_loss()));
        assert!(scorer.partition_ion_logs(&part).iter().any(|(ion, _)| *ion == intact));
        // Loss path: carries the loss ion.
        assert!(scorer.has_loss_tables());
        assert!(scorer.partition_loss_ion_logs(&part).iter().any(|(ion, _)| *ion == loss));
        // A standard model has no loss tables.
        assert!(!RankScorer::new(&tiny_param()).has_loss_tables());
    }

    #[test]
    fn unknown_ion_returns_zero() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let unknown_ion = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        // Suffix isn't in the table → return 0.
        assert_eq!(scorer.node_score(part, unknown_ion, 1), 0.0);
    }
}
