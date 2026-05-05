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

use crate::param_model::{Param, PrecursorOffsetFrequency};
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
}

impl<'a> ScoredSpectrum<'a> {
    /// Construct, filtering precursor peaks at offsets from
    /// `param.precursor_off_map[charge]` before ranking.
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
        Self::rank_kept(spec, kept, kept_count, ranks)
    }

    /// Constructor that skips precursor-peak filtering. Convenient for
    /// tests; preserves the simpler Phase 5 Task 1 API.
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
        Self::rank_kept(spec, kept, kept_count, ranks)
    }

    /// Shared ranking logic: sort `kept` by intensity DESC / mz ASC and
    /// write ranks back into the `ranks` vec. Returns the finished
    /// `ScoredSpectrum`.
    fn rank_kept(
        spec: &'a Spectrum,
        mut kept: Vec<(usize, f32, f64)>,
        kept_count: usize,
        mut ranks: Vec<u32>,
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
        Self { spec, ranks, kept_count }
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
}

#[cfg(test)]
mod tests {
    use super::*;

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
