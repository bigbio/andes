//! Per-spectrum precomputed state for Phase 5 scoring.
//!
//! Phase 5 Task 1 scope: peak ranking by intensity + nearest-peak-by-mz
//! lookup.
//!
//! Phase 5b Task 1: precursor-peak filtering before ranking. Mirrors Java's
//! `Spectrum.filterPrecursorPeaks(tolerance, reducedCharge, offset)`:
//!
//! ```text
//! // Java: edu.ucsd.msjava.msutil.Spectrum.filterPrecursorPeaks
//! public void filterPrecursorPeaks(Tolerance tolerance, int reducedCharge, float offset) {
//!     int c = this.getCharge() - reducedCharge;   // effective charge for the ion
//!     float mass = (this.getPrecursorMass() + c * ChargeCarrierMass()) / c + offset;
//!     for (Peak p : getPeakListByMass(mass, tolerance))
//!         p.setIntensity(0);
//! }
//! ```
//!
//! Where:
//! - `this.getPrecursorMass()` = `(precursor_mz - PROTON) * charge`  (neutral mass)
//! - `ChargeCarrierMass()` = `PROTON` = 1.00727649 Da
//! - `c = charge - reduced_charge`
//! - `filter_mz = (neutral_mass + c * PROTON) / c + offset`
//!   (offset is in m/z space, added after dividing by c)
//! - `getPeakListByMass` compares against each peak's m/z (not mass),
//!   so `filter_mz` is the m/z to match against
//!
//! The `precursor_off_map` (from `Param`) maps precursor charge → list of
//! `PrecursorOffsetFrequency { reduced_charge, offset, tolerance, frequency }`.
//! For each entry, any peak whose m/z is within `tolerance` Da of `filter_mz`
//! is excluded from ranking.
//!
//! Phase 6 Task 4: adds `prob_peak`, `main_ion`, `node_score`, `edge_score`,
//! and `observed_node_mass` for the GF DP graph traversal.

use crate::param_model::{IonType, Param, PrecursorOffsetFrequency};
use crate::scoring::rank_scorer::RankScorer;
use crate::spectrum::Spectrum;

const PROTON: f64 = 1.007_276_49;

#[derive(Debug, Clone)]
pub struct ScoredSpectrum<'a> {
    spec: &'a Spectrum,
    /// Per-peak rank (1 = highest intensity), aligned with `spec.peaks`
    /// indices. `ranks[i]` is the rank of the peak at index `i` in the
    /// original `spec.peaks` array. Ties broken by ascending m/z.
    /// Peaks filtered out by precursor-peak filtering receive rank `u32::MAX`.
    ranks: Vec<u32>,
    /// Number of peaks that survived precursor-peak filtering (used for
    /// `peak_count_after_filtering`).
    kept_count: usize,
    /// Probability that a random m/z bin contains a peak. Java:
    /// `probPeak = spec.size() / max(approxNumBins, 1)` where
    /// `approxNumBins = parentMass / (mme.getValue() * 2)`.
    ///
    /// For `new_without_filtering` (tests / unit use) this is set to a
    /// sentinel value of `1.0` — callers relying on `edge_score` accuracy
    /// should use the `new` constructor with a full `Param`.
    pub(crate) prob_peak: f32,
    /// The "main ion" for this spectrum's precursor partition. Used by
    /// `observed_node_mass` to look up the observed peak closest to a
    /// theoretical node mass. Set to a Prefix(charge=1, offset=0) fallback
    /// when `new_without_filtering` is used, or derived from the scorer's
    /// table when `new` is used.
    pub(crate) main_ion: IonType,
}

impl<'a> ScoredSpectrum<'a> {
    /// Construct, filtering precursor peaks at offsets from
    /// `param.precursor_off_map[charge]` before ranking. Also computes
    /// `prob_peak` and selects `main_ion` from the scorer.
    ///
    /// `charge` is the precursor charge of `spec`; if `spec.precursor_charge`
    /// is `Some(z)`, callers typically pass `z`; if `None`, pass the charge
    /// being tried by the search loop.
    ///
    /// Any peak whose m/z is within the tolerance of a precursor filter m/z
    /// gets rank `u32::MAX` and is effectively invisible to `nearest_peak_rank`.
    pub fn new(spec: &'a Spectrum, param: &Param, charge: u8) -> Self {
        let n = spec.peaks.len();

        // Collect filter m/z values from param.precursor_off_map for this charge.
        let filter_entries: &[PrecursorOffsetFrequency] = param
            .precursor_off_map
            .get(&(charge as i32))
            .map(Vec::as_slice)
            .unwrap_or(&[]);

        // Compute each filter m/z: mirror Java's filterPrecursorPeaks formula.
        // neutral_mass = (precursor_mz - PROTON) * charge
        // c = charge - reduced_charge
        // filter_mz = (neutral_mass + c * PROTON) / c + offset
        let neutral_mass = (spec.precursor_mz - PROTON) * (charge as f64);
        let filter_mzs: Vec<(f64, f64)> = filter_entries
            .iter()
            .filter_map(|pof| {
                let c = (charge as i32 - pof.reduced_charge) as f64;
                if c <= 0.0 {
                    // Would produce division by zero or negative charge; skip.
                    return None;
                }
                let filter_mz = (neutral_mass + c * PROTON) / c + (pof.offset as f64);
                let tol_da = pof.tolerance.as_da(filter_mz);
                Some((filter_mz, tol_da))
            })
            .collect();

        // Determine which peaks survive filtering.
        let ranks = vec![u32::MAX; n];
        let mut kept: Vec<(usize, f32, f64)> = Vec::with_capacity(n);
        for (i, &(mz, intensity)) in spec.peaks.iter().enumerate() {
            let filtered = filter_mzs
                .iter()
                .any(|&(fmz, tol)| (mz - fmz).abs() <= tol);
            if !filtered {
                kept.push((i, intensity, mz));
            }
        }

        let kept_count = kept.len();

        // Compute prob_peak. Java: spec.getPeptideMass() = (precursor_mz - PROTON) * charge.
        // approxNumBins = parentMass / (mme.getValue() * 2).
        // probPeak = (size if size>0 else 1) / max(approxNumBins, 1).
        //
        // Mirrors Java's mme.getValue() (raw, not Da-converted) for Java SpecEValue parity.
        // For Ppm(20), Java's getValue() returns 20.0 — NOT the Da-equivalent — so
        // approxNumBins is computed with the literal stored value.
        let parent_mass = neutral_mass; // = (precursor_mz - PROTON) * charge
        let mme_raw = param.mme.raw_value();
        let approx_num_bins = if mme_raw > 0.0 { parent_mass / (mme_raw * 2.0) } else { 1.0 };
        let peak_count = if kept_count == 0 { 1 } else { kept_count } as f64;
        let prob_peak = (peak_count / approx_num_bins.max(1.0)) as f32;

        // Select main_ion: per-partition main ion for (charge, parent_mass, last_seg).
        let last_seg = (param.num_segments - 1).max(0) as usize;
        let part = param.partition_for(charge, parent_mass, last_seg);
        let main_ion = main_ion_from_param(param, part);

        Self::rank_kept(spec, kept, kept_count, ranks, prob_peak, main_ion)
    }

    /// Constructor that skips precursor-peak filtering. Convenient for
    /// tests; preserves the simpler Phase 5 Task 1 API.
    ///
    /// Sets `prob_peak = 1.0` and `main_ion = Prefix(charge=1, offset=0)`
    /// as sentinels. For accurate `edge_score` computations, use `new`.
    pub fn new_without_filtering(spec: &'a Spectrum) -> Self {
        let n = spec.peaks.len();
        let kept: Vec<(usize, f32, f64)> = spec
            .peaks
            .iter()
            .enumerate()
            .map(|(i, &(mz, intensity))| (i, intensity, mz))
            .collect();
        let kept_count = kept.len();
        let ranks = vec![u32::MAX; n];
        let prob_peak = 1.0_f32;
        let main_ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        Self::rank_kept(spec, kept, kept_count, ranks, prob_peak, main_ion)
    }

    /// Shared ranking logic: sort `kept` by intensity DESC / mz ASC and
    /// write ranks back into the `ranks` vec. Returns the finished
    /// `ScoredSpectrum`.
    fn rank_kept(
        spec: &'a Spectrum,
        mut kept: Vec<(usize, f32, f64)>,
        kept_count: usize,
        mut ranks: Vec<u32>,
        prob_peak: f32,
        main_ion: IonType,
    ) -> Self {
        kept.sort_by(|a, b| {
            // Higher intensity first; if equal, lower m/z first.
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });
        for (rank_minus_one, &(orig_idx, _, _)) in kept.iter().enumerate() {
            ranks[orig_idx] = (rank_minus_one + 1) as u32;
        }
        Self { spec, ranks, kept_count, prob_peak, main_ion }
    }

    /// Returns `true` if the main ion is a prefix ion (b-ion direction),
    /// `false` if it is a suffix ion (y-ion direction).
    ///
    /// Mirrors Java `ScoredSpectrum.getMainIonDirection()` which returns
    /// `mainIon.isPrefixIon()`. Used by `PrimitiveAaGraph` to decide which
    /// end is the graph source.
    pub fn main_ion_direction(&self) -> bool {
        self.main_ion.is_prefix()
    }

    /// Total number of peaks in the original spectrum (before any filtering).
    pub fn peak_count(&self) -> usize {
        self.spec.peaks.len()
    }

    /// Number of peaks that survived precursor-peak filtering (and were ranked).
    pub fn peak_count_after_filtering(&self) -> usize {
        self.kept_count
    }

    /// Find the peak closest to `target_mz` within `tolerance_da`. Returns
    /// the peak's rank, or `None` if no peak falls within the window.
    ///
    /// Filtered-out peaks (rank == `u32::MAX`) are never returned.
    ///
    /// `spec.peaks` is sorted ascending by m/z (Phase 3a MGF reader
    /// guarantees this), so a binary search would be optimal; for
    /// Task 1 MVP we use a linear scan since spectrum sizes are small
    /// (typically < 2000 peaks).
    pub fn nearest_peak_rank(&self, target_mz: f64, tolerance_da: f64) -> Option<u32> {
        let mut best: Option<(usize, f64)> = None;
        for (i, &(mz, _intensity)) in self.spec.peaks.iter().enumerate() {
            // Skip filtered-out peaks.
            if self.ranks[i] == u32::MAX {
                continue;
            }
            let delta = (mz - target_mz).abs();
            if delta > tolerance_da {
                continue;
            }
            if best.as_ref().map_or(true, |(_, d)| delta < *d) {
                best = Some((i, delta));
            }
        }
        best.map(|(i, _)| self.ranks[i])
    }

    // -----------------------------------------------------------------------
    // Phase 6 / Task 4: GF DP scoring methods
    // -----------------------------------------------------------------------

    /// Mirror Java `NewScoredSpectrum.getNodeScore(prm, srm)`:
    /// `round(prefix_score(prefix_nominal) + suffix_score(suffix_nominal))`.
    ///
    /// `prefix_nominal` and `suffix_nominal` are the float node masses in Da
    /// (not integer nominal-mass indices). `parent_mass` is the precursor
    /// neutral mass. `fragment_tolerance_da` is the m/z window for peak lookup.
    pub fn node_score(
        &self,
        prefix_nominal: f64,
        suffix_nominal: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
    ) -> i32 {
        let pref = self.directional_node_score(
            prefix_nominal, true, scorer, charge, parent_mass, fragment_tolerance_da,
        );
        let suff = self.directional_node_score(
            suffix_nominal, false, scorer, charge, parent_mass, fragment_tolerance_da,
        );
        (pref + suff).round() as i32
    }

    /// Score for a single directional (prefix or suffix) node at `nominal_mass`.
    /// Mirrors the inner loop of Java's
    /// `NewScoredSpectrum.getNodeScore(nodeMass, isPrefix)`.
    fn directional_node_score(
        &self,
        nominal_mass: f64,
        is_prefix: bool,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
    ) -> f32 {
        use crate::scoring::fragment_ions::ions_for_node;
        let mut total = 0.0_f32;
        for (ion, theo_mz) in ions_for_node(nominal_mass, is_prefix, scorer.param(), parent_mass, charge) {
            let seg = scorer.param().segment_num(theo_mz, parent_mass);
            let part = scorer.param().partition_for(charge, parent_mass, seg);
            match self.nearest_peak_rank(theo_mz, fragment_tolerance_da) {
                Some(rank) => total += scorer.node_score(part, ion, rank),
                None => total += scorer.missing_ion_score(part, ion),
            }
        }
        total
    }

    /// Return the observed node mass for `node_nominal`, or `None` if no
    /// peak is near the theoretical m/z of the main ion.
    ///
    /// Java: `NewScoredSpectrum.getNodeMass(node)` — computes
    /// `theo_mz = main_ion.getMz(node.getMass())`, calls
    /// `spec.getPeakByMass(theoMass, scorer.getMME())` which returns the
    /// **highest-intensity** peak within the MME window (via
    /// `Collections.max(matchList, new IntensityComparator())`), then returns
    /// `main_ion.getMass(peak_mz)` if found, else `-1` (we use `None`).
    ///
    /// Single-pass: uses `scorer.param().mme.as_da(theo_mz)` as the window,
    /// selects the highest-intensity non-filtered peak in that window.
    pub fn observed_node_mass(
        &self,
        node_nominal: i32,
        scorer: &RankScorer,
        charge: u8,
        _parent_mass: f64,
    ) -> Option<f64> {
        let _ = charge; // not needed in formula; kept for API symmetry
        if node_nominal == 0 {
            // Source node mass is exactly 0 by convention (Java returns 0 when
            // `node.getNominalMass() == 0`).
            return Some(0.0);
        }
        let theo_mz = self.main_ion.mz(node_nominal as f64);
        let tol_da = scorer.param().mme.as_da(theo_mz);
        // Select the highest-intensity peak within [theo_mz - tol_da, theo_mz + tol_da].
        // Mirrors Java's Collections.max(matchList, new IntensityComparator()).
        // Skip filtered peaks (ranks[i] == u32::MAX).
        let mut best_peak_mz: Option<(f64, f32)> = None; // (mz, intensity)
        for (i, &(mz, intensity)) in self.spec.peaks.iter().enumerate() {
            if self.ranks[i] == u32::MAX {
                continue;
            }
            if (mz - theo_mz).abs() > tol_da {
                continue;
            }
            if best_peak_mz.as_ref().map_or(true, |&(_, best_int)| intensity > best_int) {
                best_peak_mz = Some((mz, intensity));
            }
        }
        best_peak_mz.map(|(peak_mz, _)| self.main_ion.mass_from_mz(peak_mz))
    }

    /// Mirror Java `NewScoredSpectrum.getEdgeScore(curNode, prevNode, theoMass)`.
    ///
    /// If `param.ion_existence_table` is empty (Java's `!scorer.supportEdgeScores()`),
    /// returns 0. Otherwise:
    ///   1. Look up observed node masses for `cur_nominal` and `prev_nominal`.
    ///   2. `ion_existence_index` = (cur observed?) + 2*(prev observed?).
    ///   3. `score = ion_existence_score(part, idx, prob_peak)`.
    ///   4. If `idx == 3` (both observed), also add `error_score(cur_mass - prev_mass - theo_aa_mass)`.
    ///   5. Return `round(score) as i32`.
    pub fn edge_score(
        &self,
        cur_nominal: i32,
        prev_nominal: i32,
        theo_aa_mass: f64,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
    ) -> i32 {
        // Java: if (!scorer.supportEdgeScores()) return 0;
        // supportEdgeScores() ↔ errorScalingFactor != 0.
        if scorer.param().error_scaling_factor == 0 {
            return 0;
        }
        if scorer.param().ion_existence_table.is_empty() {
            return 0;
        }

        // 1. Observed masses for cur and prev nodes.
        let cur_mass = self.observed_node_mass(cur_nominal, scorer, charge, parent_mass);
        let prev_mass = self.observed_node_mass(prev_nominal, scorer, charge, parent_mass);

        // 2. ion_existence_index: 1 if cur observed, +2 if prev observed.
        let mut idx = 0usize;
        if cur_mass.is_some() { idx += 1; }
        if prev_mass.is_some() { idx += 2; }

        // 3. Partition for this spectrum — Java uses the "last segment" partition
        //    stored at construction time.
        let last_seg = (scorer.param().num_segments - 1).max(0) as usize;
        let part = scorer.param().partition_for(charge, parent_mass, last_seg);

        // 4. Ion existence score.
        let mut s = scorer.ion_existence_score(part, idx, self.prob_peak);

        // 5. If both observed, add error score.
        if idx == 3 {
            let delta = cur_mass.unwrap() - prev_mass.unwrap() - theo_aa_mass;
            s += scorer.error_score(part, delta as f32);
        }

        s.round() as i32
    }
}

/// Select the main ion for `partition` from `param.rank_dist_table`.
/// Java: per-partition main ion = highest-frequency-at-rank-1 prefix ion.
/// We pick the Prefix ion with the highest freq at rank-1 index (index 0).
/// Falls back to `Prefix { charge: 1, offset_bits: 0 }` if table is empty.
///
/// TODO(Phase 6 followup): main_ion selection currently uses per-partition
/// rank-1 prefix-ion frequency from rank_dist_table. Java's
/// NewRankScorer.determineIonTypes (lines 611-640) aggregates fragOFFTable
/// across segments and considers all ion types. For HCD these agree; for
/// ETD/ECD they may diverge. Revisit if Task 9/10 parity tests fail.
fn main_ion_from_param(param: &Param, partition: crate::param_model::Partition) -> IonType {
    let fallback = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let table = match param.rank_dist_table.get(&partition) {
        Some(t) => t,
        None => return fallback,
    };
    let mut best_ion = None;
    let mut best_freq = f32::NEG_INFINITY;
    for (ion, freqs) in table {
        if !ion.is_prefix() {
            continue;
        }
        let freq_at_rank1 = freqs.first().copied().unwrap_or(0.0);
        if freq_at_rank1 > best_freq {
            best_freq = freq_at_rank1;
            best_ion = Some(*ion);
        }
    }
    best_ion.unwrap_or(fallback)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::param_model::{IonType, Partition};
    use crate::scoring::rank_scorer::RankScorer;

    fn spec(peaks: &[(f64, f32)]) -> Spectrum {
        Spectrum {
            title: "test".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: peaks.to_vec(),
        }
    }

    /// A richer `tiny_param` that has a Prefix(charge=1, offset=0) ion in the
    /// rank_dist_table under partition (charge=2, parent_mass=1000.0, seg_num=0),
    /// so that `ion_types_for_segment(0)` returns a non-empty list and
    /// `node_score` / `edge_score` can exercise the scoring paths.
    fn tiny_param_with_ions() -> Param {
        use crate::activation::ActivationMethod;
        use crate::instrument::InstrumentType;
        use crate::param_model::{FragmentOffsetFrequency, SpecDataType};
        use crate::protocol::Protocol;
        use crate::tolerance::Tolerance;
        use std::collections::HashMap;

        let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let noise = IonType::Noise;

        // max_rank=3 → 4 slots. Ion has higher freq at rank 1.
        let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
        let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

        let mut ion_table: HashMap<IonType, Vec<f32>> = HashMap::new();
        ion_table.insert(prefix1, ion_freqs);
        ion_table.insert(noise, noise_freqs);

        let mut rank_dist_table: HashMap<Partition, HashMap<IonType, Vec<f32>>> = HashMap::new();
        rank_dist_table.insert(part, ion_table);

        // frag_off_table: one prefix ion entry so ion_types_for_segment returns it.
        let mut frag_off_table = HashMap::new();
        frag_off_table.insert(part, vec![FragmentOffsetFrequency {
            ion_type: prefix1,
            frequency: 0.7,
        }]);

        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Da(0.5),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: HashMap::new(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
        }
    }

    // --- Phase 6 / Task 4 spec-review fix: prob_peak uses raw mme value ---

    /// Verify that `prob_peak` is computed using the raw stored mme value (Java
    /// parity), not the Da-converted form. For `Tolerance::Ppm(20.0)`:
    ///   Java formula: approxNumBins = parentMass / (mme.getValue() * 2)
    ///                               = parentMass / (20.0 * 2)
    ///   NOT:          parentMass / (as_da(parentMass) * 2)
    ///                               = parentMass / (parentMass * 20e-6 * 2)
    #[test]
    fn prob_peak_uses_raw_mme_value_not_da_converted() {
        use crate::activation::ActivationMethod;
        use crate::instrument::InstrumentType;
        use crate::param_model::SpecDataType;
        use crate::protocol::Protocol;
        use crate::tolerance::Tolerance;
        use std::collections::HashMap;

        // Spectrum: precursor_mz=501.00727649 → neutral_mass≈(501.007-PROTON)*2≈1000.0 Da,
        // charge=2.
        let precursor_mz = 501.007_276_49_f64; // ≈ (1000/2) + PROTON
        let s = Spectrum {
            title: "parity_test".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
        };

        let param = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![],
            num_precursor_off: 0,
            precursor_off_map: HashMap::new(),
            frag_off_table: HashMap::new(),
            max_rank: 3,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
        };

        let ss = ScoredSpectrum::new(&s, &param, 2);

        // Expected: raw_value = 20.0, parent_mass ≈ (501.007276 - PROTON) * 2.
        let parent_mass = (precursor_mz - PROTON) * 2.0;
        let raw_mme = 20.0_f64;
        let approx_num_bins = parent_mass / (raw_mme * 2.0);
        let expected_prob_peak = (3.0_f64 / approx_num_bins.max(1.0)) as f32;

        // The Da-converted form would be: parent_mass / (parent_mass * 20e-6 * 2) ≈ 25_000.0,
        // giving prob_peak ≈ 3/25000 = 0.00012, not the raw-value result ≈ 3/100 = 0.06.
        let wrong_approx_num_bins = parent_mass / (parent_mass * 20e-6 * 2.0);
        let wrong_prob_peak = (3.0_f64 / wrong_approx_num_bins.max(1.0)) as f32;

        // Sanity: raw and Da results must differ significantly for this to be a meaningful test.
        assert!(
            (expected_prob_peak - wrong_prob_peak).abs() > 0.001,
            "test precondition failed: Ppm raw vs Da-converted did not produce different prob_peak values"
        );

        assert!(
            (ss.prob_peak - expected_prob_peak).abs() < 1e-5,
            "prob_peak={} but expected={} (Java raw formula). Wrong Da-converted value would be {}",
            ss.prob_peak, expected_prob_peak, wrong_prob_peak
        );
    }

    // --- Phase 6 / Task 4 followup: observed_node_mass picks highest-intensity ---

    #[test]
    fn observed_node_mass_picks_highest_intensity_peak_in_window() {
        // Two peaks within the MME window of theo_mz; the higher-intensity one wins.
        // tiny_param_with_ions uses Tolerance::Da(0.5) → window ±0.5 Da.
        // main_ion = Prefix { charge: 1, offset_bits: 0 }
        // theo_mz = (node_nominal + PROTON) / 1 = node_nominal + PROTON ≈ node_nominal + 1.00728.
        // We use node_nominal = 100 → theo_mz ≈ 101.00728.
        // Place two peaks both within ±0.5:
        //   peak A at 101.1 (delta ≈ 0.093, low intensity 1.0) — CLOSER
        //   peak B at 101.4 (delta ≈ 0.393, high intensity 100.0) — FARTHER but HIGHER intensity
        // Java picks highest-intensity → peak B must win.
        use crate::mass::PROTON;
        let node_nominal = 100_i32;
        let theo_mz = node_nominal as f64 + PROTON;
        let closer_mz = theo_mz + 0.093; // delta 0.093 < 0.393
        let farther_mz = theo_mz + 0.393; // still within ±0.5
        let s = spec(&[(closer_mz, 1.0), (farther_mz, 100.0)]);
        let param = tiny_param_with_ions(); // mme = Da(0.5)
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let result = ss.observed_node_mass(node_nominal, &scorer, 2, 1000.0);
        let result_mass = result.expect("should find a peak in the window");
        // main_ion.mass_from_mz(peak_mz): farther_mz * 1 - PROTON - 0 = farther_mz - PROTON
        let expected_mass = farther_mz - PROTON;
        let wrong_mass = closer_mz - PROTON;
        assert!(
            (result_mass - expected_mass).abs() < 1e-6,
            "expected highest-intensity (farther) peak mass {expected_mass:.6}, \
             got {result_mass:.6} (closest/wrong would be {wrong_mass:.6})"
        );
    }

    // --- Phase 6 / Task 4 tests: node_score and edge_score ---

    #[test]
    fn node_score_does_not_panic_on_empty_spectrum() {
        // Spectrum with no peaks; every ion is missing → all contributions
        // come from missing_ion_score. With no matching peaks the missing
        // score for Prefix(charge=1) is log(0.001/0.4) < 0, but we also
        // include the suffix side which has no ions. Sum rounds to a small
        // negative.
        let s = spec(&[]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(100.0, 900.0, &scorer, 2, 1000.0, 0.5);
        // With empty ion_types_for_segment the suffix side contributes 0,
        // and no suffix ions are in the table → suffix score is 0.
        // The prefix missing-ion score is negative → total rounds negative or 0.
        assert!(n <= 0, "missing-ion score on empty spectrum should be non-positive, got {n}");
    }

    #[test]
    fn node_score_nonzero_when_peak_matches_prefix_ion() {
        // Place a high-intensity peak at the predicted b1 m/z for a node of
        // nominal mass = 100. Prefix ion: Prefix(charge=1, offset=0).
        // theo_mz = (100.0 + 0 + PROTON) / 1 = 100 + 1.00727649
        use crate::mass::PROTON;
        let nominal = 100.0_f64;
        let b1_mz = nominal + PROTON; // charge=1, offset=0
        let s = spec(&[(50.0, 1.0), (b1_mz, 100.0), (200.0, 2.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // prefix_nominal = 100, suffix_nominal = 900 (doesn't matter, no suffix ions in table).
        let n = ss.node_score(nominal, 900.0, &scorer, 2, 1000.0, 0.5);
        // Peak at b1_mz gets rank 1 (highest intensity = 100.0).
        // node_score(rank=1, Prefix) = log(0.6 / (0.1 * 1)) = log(6) > 0.
        // Total suffix = 0. Round(log(6)) = round(1.79) = 2.
        assert!(n > 0, "expected positive node_score when b-ion peak present, got {n}");
    }

    #[test]
    fn node_score_prefix_only_match() {
        // Only prefix ions in table; suffix side always contributes 0.
        use crate::mass::PROTON;
        let nominal = 57.0_f64; // roughly glycine residue mass
        let mz = (nominal + PROTON) / 1.0;
        let s = spec(&[(mz, 50.0), (300.0, 1.0)]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(nominal, 900.0, &scorer, 2, 1000.0, 0.5);
        // Peak at mz is rank 1. score = log(0.6 / 0.1) = log(6) ≈ 1.79 → rounds to 2.
        assert!(n > 0, "prefix-only match: expected positive score, got {n}");
    }

    #[test]
    fn node_score_no_matching_ions_returns_negative_or_zero() {
        // With a peak far from any ion, all ions are missing → negative score.
        let s = spec(&[(5000.0, 100.0)]); // peak far from any fragment ion
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(100.0, 900.0, &scorer, 2, 1000.0, 0.5);
        // missing_ion_score for Prefix(1) = log(0.001/0.4) < 0 → n <= 0.
        assert!(n <= 0, "missing ion should produce non-positive score, got {n}");
    }

    #[test]
    fn node_score_nominal_mass_zero_prefix_returns_zero() {
        // nominal_mass = 0 is the source node. Java's getNodeScore treats source
        // and sink nodes specially, but our Rust impl evaluates ions_for_node(0.0, …).
        // With prefix_nominal=0 and suffix_nominal=1000 (parent mass), and no peaks
        // in the spectrum, the missing-ion score for the Prefix ion governs.
        // The suffix nominal = 1000 > parent_mass → ions_for_node produces no suffix
        // ions for that degenerate case. Net result: non-positive score.
        let s = spec(&[]);
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let n = ss.node_score(0.0, 1000.0, &scorer, 2, 1000.0, 0.5);
        // Score is non-positive (missing-ion penalty applies).
        assert!(n <= 0, "source-node score with empty spectrum should be non-positive, got {n}");
    }

    #[test]
    fn edge_score_returns_zero_when_table_empty() {
        // No ion_existence_table → Java path returns 0.
        let s = spec(&[(100.0, 1.0)]);
        let mut param = tiny_param_with_ions();
        param.ion_existence_table.clear();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
        assert_eq!(e, 0);
    }

    #[test]
    fn edge_score_returns_zero_when_error_scaling_factor_zero() {
        // error_scaling_factor == 0 ↔ supportEdgeScores() == false → returns 0.
        let s = spec(&[(100.0, 1.0)]);
        let param = tiny_param_with_ions(); // error_scaling_factor defaults to 0
        assert_eq!(param.error_scaling_factor, 0);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
        assert_eq!(e, 0);
    }

    #[test]
    fn edge_score_nonzero_with_existence_table() {
        // Build a param with error_scaling_factor > 0 and a populated
        // ion_existence_table. Check that edge_score is computed (non-zero).
        use crate::activation::ActivationMethod;
        use crate::instrument::InstrumentType;
        use crate::param_model::{FragmentOffsetFrequency, SpecDataType};
        use crate::protocol::Protocol;
        use crate::tolerance::Tolerance;
        use std::collections::HashMap;

        let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let noise = IonType::Noise;

        let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
        let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

        let mut ion_table: HashMap<IonType, Vec<f32>> = HashMap::new();
        ion_table.insert(prefix1, ion_freqs);
        ion_table.insert(noise, noise_freqs);

        let mut rank_dist_table: HashMap<Partition, HashMap<IonType, Vec<f32>>> = HashMap::new();
        rank_dist_table.insert(part, ion_table);

        let mut frag_off_table = HashMap::new();
        frag_off_table.insert(part, vec![FragmentOffsetFrequency {
            ion_type: prefix1,
            frequency: 0.7,
        }]);

        // error_scaling_factor = 2 → dist_len = 5; ion_existence = 4 entries
        let error_scaling_factor = 2_i32;
        let dist_len = (error_scaling_factor as usize) * 2 + 1;

        let mut ion_err_dist_table: HashMap<Partition, Vec<f32>> = HashMap::new();
        ion_err_dist_table.insert(part, vec![0.1_f32, 0.2, 0.4, 0.2, 0.1]);

        let mut noise_err_dist_table: HashMap<Partition, Vec<f32>> = HashMap::new();
        noise_err_dist_table.insert(part, vec![0.05_f32, 0.1, 0.7, 0.1, 0.05]);

        let mut ion_existence_table: HashMap<Partition, Vec<f32>> = HashMap::new();
        // [nn, ?, ?, yy] = [0.1, 0.3, 0.3, 0.5]
        ion_existence_table.insert(part, vec![0.1_f32, 0.3, 0.3, 0.5]);

        let _ = dist_len; // used for documentation

        let param = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Da(0.5),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: HashMap::new(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor,
            ion_err_dist_table,
            noise_err_dist_table,
            ion_existence_table,
        };

        // No peaks in spectrum → cur_mass = None, prev_mass = None → idx = 0 (nn).
        let s = spec(&[]);
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
        // ion_existence_score(part, 0, prob_peak): ionExistenceProb[0]=0.1,
        // noiseExistenceProb = (1-p)^2. With many bins prob_peak ≈ 0.
        // log(0.1 / ~1.0) = ~log(0.1) ≈ -2.3 → rounds to -2.
        // Confirm the table is used (non-zero result).
        assert_ne!(e, 0, "edge_score should be nonzero with populated existence table");
    }

    #[test]
    fn empty_spectrum_yields_no_ranks() {
        let s = spec(&[]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        assert_eq!(ss.peak_count(), 0);
        assert!(ss.nearest_peak_rank(500.0, 0.1).is_none());
    }

    #[test]
    fn highest_intensity_gets_rank_1() {
        // Peaks sorted ascending by m/z (Phase 3a MGF reader does this).
        let s = spec(&[(100.0, 1.0), (200.0, 5.0), (300.0, 3.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        assert_eq!(ss.peak_count(), 3);
        // Peak at m/z 200 has the highest intensity (5.0) → rank 1.
        // The lookup window of 0.1 should find it.
        assert_eq!(ss.nearest_peak_rank(200.0, 0.1), Some(1));
        // Peak at m/z 300 has intensity 3.0 → rank 2.
        assert_eq!(ss.nearest_peak_rank(300.0, 0.1), Some(2));
        // Peak at m/z 100 has intensity 1.0 → rank 3 (lowest).
        assert_eq!(ss.nearest_peak_rank(100.0, 0.1), Some(3));
    }

    #[test]
    fn nearest_peak_within_tolerance() {
        let s = spec(&[(100.0, 1.0), (200.5, 5.0), (300.0, 3.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // Target 200.4 with tol 0.2 → finds peak at 200.5 (within 0.1).
        assert_eq!(ss.nearest_peak_rank(200.4, 0.2), Some(1));
        // Target 200.5 with tol 0.001 → exact match.
        assert_eq!(ss.nearest_peak_rank(200.5, 0.001), Some(1));
        // Target 200.4 with tol 0.05 → outside window, no match.
        assert_eq!(ss.nearest_peak_rank(200.4, 0.05), None);
    }

    #[test]
    fn ties_broken_deterministically() {
        // Two peaks with identical intensity — the lower m/z gets rank 1
        // (matching Java's behavior of sort stability + ties going to
        // earlier-indexed peaks).
        let s = spec(&[(100.0, 5.0), (200.0, 5.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // Both peaks should have a defined rank; the test asserts the
        // ranking is total (no two peaks share a rank).
        let r1 = ss.nearest_peak_rank(100.0, 0.1).unwrap();
        let r2 = ss.nearest_peak_rank(200.0, 0.1).unwrap();
        assert_ne!(r1, r2);
        assert!(r1 == 1 || r2 == 1);
        assert!(r1 == 2 || r2 == 2);
    }

    #[test]
    fn closest_among_multiple_in_tolerance() {
        // Multiple peaks within the tolerance window; the closest wins.
        let s = spec(&[(99.5, 1.0), (100.0, 5.0), (100.5, 2.0)]);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        // Target 100.1 with tol 0.6: all three are within. Closest is 100.0 → rank 1.
        assert_eq!(ss.nearest_peak_rank(100.1, 0.6), Some(1));
    }
}

#[cfg(test)]
mod precursor_filter_tests {
    use super::*;
    use crate::activation::ActivationMethod;
    use crate::instrument::InstrumentType;
    use crate::param_model::{Param, PrecursorOffsetFrequency, SpecDataType};
    use crate::protocol::Protocol;
    use crate::tolerance::Tolerance;
    use std::collections::HashMap;

    /// Build a Param with a single precursor offset entry: charge 2,
    /// reduced_charge 2, offset 0.0 Da (the precursor itself), tolerance 0.5 Da.
    fn param_with_precursor_filter() -> Param {
        let mut precursor_off_map: HashMap<i32, Vec<PrecursorOffsetFrequency>> = HashMap::new();
        precursor_off_map.insert(
            2,
            vec![PrecursorOffsetFrequency {
                reduced_charge: 2,
                offset: 0.0,
                tolerance: Tolerance::Da(0.5),
                frequency: 1.0,
            }],
        );

        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![],
            num_precursor_off: 1,
            precursor_off_map,
            frag_off_table: HashMap::new(),
            max_rank: 3,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
        }
    }

    fn make_spec(precursor_mz: f64, peaks: &[(f64, f32)]) -> Spectrum {
        Spectrum {
            title: "test".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: peaks.to_vec(),
        }
    }

    /// Verify the filter_mz formula for reduced_charge=2, offset=0:
    /// neutral_mass = (500.0 - PROTON) * 2 = 997.985450...
    /// c = 2 - 2 = 0 → filtered (c <= 0), so no filtering happens.
    ///
    /// Re-check: the task says "charge 2, reduced_charge 2" for the
    /// precursor itself. With c = charge - reduced_charge = 0, that
    /// would be division by zero. Real param files use reduced_charge < charge.
    ///
    /// Let's use reduced_charge=0 for the precursor filter test:
    /// c = 2 - 0 = 2; filter_mz = (neutral + 2*PROTON) / 2 + 0 = precursor_mz.
    fn param_with_precursor_filter_rc0() -> Param {
        let mut precursor_off_map: HashMap<i32, Vec<PrecursorOffsetFrequency>> = HashMap::new();
        precursor_off_map.insert(
            2,
            vec![PrecursorOffsetFrequency {
                reduced_charge: 0,
                offset: 0.0,
                tolerance: Tolerance::Da(0.5),
                frequency: 1.0,
            }],
        );

        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![],
            num_precursor_off: 1,
            precursor_off_map,
            frag_off_table: HashMap::new(),
            max_rank: 3,
            rank_dist_table: HashMap::new(),
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
        }
    }

    #[test]
    fn precursor_peak_is_filtered_out() {
        // precursor m/z = 500.0, charge 2, reduced_charge=0:
        // c = 2 - 0 = 2
        // neutral_mass = (500.0 - PROTON) * 2 ≈ 997.9855 Da
        // filter_mz = (997.9855 + 2 * PROTON) / 2 + 0.0 = 500.0 (the precursor m/z)
        //
        // A peak AT 500.0 (the precursor m/z itself, very high intensity) should be filtered.
        let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0), (300.0, 5.0)]);
        let param = param_with_precursor_filter_rc0();
        let ss = ScoredSpectrum::new(&s, &param, 2);

        // The precursor peak (500.0) should be filtered out (rank u32::MAX, not returned).
        assert!(
            ss.nearest_peak_rank(500.0, 0.1).is_none(),
            "precursor peak at 500.0 should be filtered, but a peak at that m/z was found"
        );

        // The other peaks should still be present and ranked.
        // (300.0, 5.0) is now rank 1 (highest among non-filtered);
        // (100.0, 1.0) is rank 2.
        assert_eq!(ss.nearest_peak_rank(300.0, 0.1), Some(1));
        assert_eq!(ss.nearest_peak_rank(100.0, 0.1), Some(2));
    }

    #[test]
    fn non_precursor_peaks_kept() {
        // Without filtering hitting any peak, all peaks should be present.
        // The filter is at precursor m/z = 500.0 ± 0.5, no peak in this set is there.
        let s = make_spec(500.0, &[(100.0, 1.0), (200.0, 50.0), (300.0, 5.0)]);
        let param = param_with_precursor_filter_rc0();
        let ss = ScoredSpectrum::new(&s, &param, 2);

        assert_eq!(ss.peak_count_after_filtering(), 3);
        assert_eq!(ss.nearest_peak_rank(200.0, 0.1), Some(1));
    }

    #[test]
    fn missing_precursor_off_map_falls_back_to_unfiltered() {
        // If param has no precursor offsets for this charge, all peaks
        // are kept and ranked normally.
        let mut param = param_with_precursor_filter_rc0();
        param.precursor_off_map.clear();
        let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0)]);
        let ss = ScoredSpectrum::new(&s, &param, 2);
        assert_eq!(ss.peak_count_after_filtering(), 2);
    }

    #[test]
    fn invalid_reduced_charge_skipped() {
        // reduced_charge >= charge → c = 0 → skip (no div-by-zero).
        // Using param_with_precursor_filter which has reduced_charge=2, charge=2.
        let param = param_with_precursor_filter();
        let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0)]);
        let ss = ScoredSpectrum::new(&s, &param, 2);
        // No filtering occurred (c <= 0 was skipped) → both peaks kept.
        assert_eq!(ss.peak_count_after_filtering(), 2);
    }
}
