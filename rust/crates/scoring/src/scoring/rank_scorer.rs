//! Per-ion rank score lookup.
//!
//! Score formula:
//!   chargeOrSeg = min(ionType.charge, numSegments)
//!   log_score[i] = log(ion_freq[i] / (noise_freq[i] * chargeOrSeg))
//!
//! Rank-distribution arrays have length `maxRank + 1`. Indices `[0..maxRank-1]`
//! correspond to ranks 1..maxRank. Index `maxRank` (the last) is the
//! "missing ion" slot, used by `missing_ion_score`.

use std::collections::HashMap;

use crate::param_model::{IonType, Param, Partition};

#[derive(Debug, Clone)]
pub struct RankScorer {
    /// The `Param` this scorer was built from. Cloned at construction so
    /// that `match_engine` can forward precursor-filter information to
    /// `ScoredSpectrum::new` without a separate `Param` argument.
    param: Param,
    /// Cached log scores: `(partition, non-noise ion_type) → Vec<f32>` where
    /// the Vec has length `max_rank + 1` (indices 0..max_rank-1 for ranks
    /// 1..max_rank, index max_rank for the missing-ion slot).
    /// Retained for `node_score`/`missing_ion_score` API compatibility (tests,
    /// diagnostics). Hot-path callers should use `partition_ion_logs` instead.
    pub(crate) log_table: HashMap<(Partition, IonType), Vec<f32>>,
    /// Dense-indexed log tables for the hot path. Keyed by `Partition` →
    /// `Vec<(IonType, Vec<f32>)>` parallel to that partition's ion list
    /// (Noise excluded). The inner `(IonType, Vec<f32>)` pairs are in the
    /// same order as `Param::ion_types_for_partition_slice`. Replaces the
    /// per-ion HashMap lookup in `directional_node_score` with array
    /// indexing — eliminates ~200M HashMap operations per PXD001819 search.
    pub(crate) partition_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>>,
    /// Cached `min(rank - 1, max_rank - 1)` clamp constant.
    max_rank: u32,
}

impl RankScorer {
    pub fn new(param: &Param) -> Self {
        let mut log_table: HashMap<(Partition, IonType), Vec<f32>> = HashMap::new();

        for (partition, ion_table) in &param.rank_dist_table {
            // Noise frequencies come from the IonType::Noise entry in the
            // same partition's rank-dist table. Skip if absent.
            let noise_freqs = match ion_table.get(&IonType::Noise) {
                Some(v) => v,
                None => continue,
            };

            for (ion_type, ion_freqs) in ion_table {
                if matches!(ion_type, IonType::Noise) {
                    continue;
                }
                let charge = match ion_type {
                    IonType::Prefix { charge, .. } | IonType::Suffix { charge, .. } => *charge,
                    IonType::Noise => unreachable!(),
                };
                // chargeOrSeg = min(ion.charge, num_segments).
                let charge_or_seg = (charge as u32).min(param.num_segments as u32) as f32;
                let n = ion_freqs.len().min(noise_freqs.len());
                let mut logs = Vec::with_capacity(n);
                for i in 0..n {
                    let ion_f = ion_freqs[i];
                    let noise_f = noise_freqs[i] * charge_or_seg;
                    logs.push((ion_f / noise_f).ln());
                }
                log_table.insert((*partition, *ion_type), logs);
            }
        }

        // Build the dense partition_ion_logs cache. For each partition, walk
        // its ion list (in the same order as
        // `Param::ion_types_for_partition_slice`) and pair each ion with its
        // log table (cloned). Used by the hot path to avoid HashMap lookups
        // per-ion; the partition is constant per outer-segment iteration in
        // `directional_node_score`.
        let mut partition_ion_logs: HashMap<Partition, Vec<(IonType, Vec<f32>)>> = HashMap::new();
        for (&partition, ions) in &param.partition_ion_types_cache {
            let mut paired: Vec<(IonType, Vec<f32>)> = Vec::with_capacity(ions.len());
            for &ion in ions {
                if let Some(logs) = log_table.get(&(partition, ion)) {
                    paired.push((ion, logs.clone()));
                }
            }
            partition_ion_logs.insert(partition, paired);
        }

        Self {
            param: param.clone(),
            log_table,
            partition_ion_logs,
            max_rank: param.max_rank as u32,
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

    /// Maximum rank used for clamping. Exposed so callers can apply
    /// rank-clamp / missing-ion semantics without going through `node_score`.
    pub fn max_rank(&self) -> u32 {
        self.max_rank
    }

    /// Return the `Param` this scorer was built from.
    pub fn param(&self) -> &Param {
        &self.param
    }

    /// Score a peak-matched ion at rank `rank` (1-based, 1 = highest intensity).
    /// `rank > max_rank` clamps to `rank = max_rank` (so the rank index
    /// becomes `max_rank - 1`, the LAST observed-rank entry, NOT the
    /// missing-ion sentinel).
    pub fn node_score(&self, partition: Partition, ion_type: IonType, rank: u32) -> f32 {
        let logs = match self.log_table.get(&(partition, ion_type)) {
            Some(v) => v,
            None => return 0.0,
        };
        let rank_clamped = rank.min(self.max_rank).max(1);
        let idx = (rank_clamped - 1) as usize;
        if idx < logs.len() {
            logs[idx]
        } else {
            0.0
        }
    }

    /// Score for an ion that isn't observed in the spectrum. Uses the slot
    /// at index `max_rank` (the LAST entry in the `max_rank + 1`-length array).
    pub fn missing_ion_score(&self, partition: Partition, ion_type: IonType) -> f32 {
        let logs = match self.log_table.get(&(partition, ion_type)) {
            Some(v) => v,
            None => return 0.0,
        };
        let idx = self.max_rank as usize;
        if idx < logs.len() {
            logs[idx]
        } else {
            0.0
        }
    }

    /// Ion-existence score.
    ///
    /// Computes `log(ionExistenceProb[index] / noiseExistenceProb)` where:
    /// - `index == 0` (nn): `noiseProb = (1 - probPeak)^2`
    /// - `index == 3` (yy): `noiseProb = probPeak^2`
    /// - otherwise: `noiseProb = probPeak * (1 - probPeak)`
    ///
    /// Returns 0.0 if the `ion_existence_table` has no entry for `part`.
    pub fn ion_existence_score(&self, partition: Partition, index: usize, prob_peak: f32) -> f32 {
        let table = match self.param.ion_existence_table.get(&partition) {
            Some(t) => t,
            None => return 0.0,
        };
        if index >= table.len() {
            return 0.0;
        }
        let noise_existence_prob = match index {
            0 => (1.0 - prob_peak) * (1.0 - prob_peak),
            3 => prob_peak * prob_peak,
            _ => prob_peak * (1.0 - prob_peak),
        };
        let mut ion_prob = table[index];
        // Zero-probability slots are clamped to 0.01 to avoid log(0).
        if ion_prob == 0.0 {
            ion_prob = 0.01;
        }
        let denom = noise_existence_prob.max(f32::MIN_POSITIVE);
        (ion_prob / denom).ln()
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

    #[test]
    fn node_score_log_formula() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };

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
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };

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
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };

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
        let ion3 = IonType::Prefix { charge: 3, offset_bits: 0.0_f32.to_bits() };
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
        let ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        // Out-of-table partition → return 0 (neutral score).
        assert_eq!(scorer.node_score(unknown, ion, 1), 0.0);
        assert_eq!(scorer.missing_ion_score(unknown, ion), 0.0);
    }

    #[test]
    fn unknown_ion_returns_zero() {
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let unknown_ion = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        // Suffix isn't in the table → return 0.
        assert_eq!(scorer.node_score(part, unknown_ion, 1), 0.0);
    }
}
