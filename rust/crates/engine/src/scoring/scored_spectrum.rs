//! Per-spectrum precomputed state for Phase 5 scoring.
//!
//! Phase 5 Task 1 scope: peak ranking by intensity + nearest-peak-by-mz
//! lookup. Tasks 2-7 will add fragment-ion prediction, rank scoring,
//! integration into match_engine, and (stretch) precursor-peak filtering
//! plus deconvolution.

use crate::spectrum::Spectrum;

#[derive(Debug, Clone)]
pub struct ScoredSpectrum<'a> {
    spec: &'a Spectrum,
    /// Per-peak rank (1 = highest intensity), aligned with `spec.peaks`
    /// indices. `ranks[i]` is the rank of the peak at index `i` in the
    /// original `spec.peaks` array. Ties broken by ascending m/z.
    ranks: Vec<u32>,
}

impl<'a> ScoredSpectrum<'a> {
    pub fn new(spec: &'a Spectrum) -> Self {
        let n = spec.peaks.len();
        // Build (orig_idx, intensity, mz) tuples, sort by intensity DESC,
        // tie-break by mz ASC, then assign rank.
        let mut indexed: Vec<(usize, f32, f64)> = spec.peaks.iter()
            .enumerate()
            .map(|(i, &(mz, intensity))| (i, intensity, mz))
            .collect();
        indexed.sort_by(|a, b| {
            // Higher intensity first; if equal, lower m/z first.
            b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });

        let mut ranks = vec![0u32; n];
        for (rank_minus_one, &(orig_idx, _, _)) in indexed.iter().enumerate() {
            ranks[orig_idx] = (rank_minus_one + 1) as u32;
        }

        Self { spec, ranks }
    }

    pub fn peak_count(&self) -> usize { self.spec.peaks.len() }

    /// Find the peak closest to `target_mz` within `tolerance_da`. Returns
    /// the peak's rank, or `None` if no peak falls within the window.
    /// `spec.peaks` is sorted ascending by m/z (Phase 3a MGF reader
    /// guarantees this), so a binary search would be optimal; for
    /// Task 1 MVP we use a linear scan since spectrum sizes are small
    /// (typically < 2000 peaks).
    pub fn nearest_peak_rank(&self, target_mz: f64, tolerance_da: f64) -> Option<u32> {
        let mut best: Option<(usize, f64)> = None;
        for (i, &(mz, _intensity)) in self.spec.peaks.iter().enumerate() {
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
        let ss = ScoredSpectrum::new(&s);
        assert_eq!(ss.peak_count(), 0);
        assert!(ss.nearest_peak_rank(500.0, 0.1).is_none());
    }

    #[test]
    fn highest_intensity_gets_rank_1() {
        // Peaks sorted ascending by m/z (Phase 3a MGF reader does this).
        let s = spec(&[(100.0, 1.0), (200.0, 5.0), (300.0, 3.0)]);
        let ss = ScoredSpectrum::new(&s);
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
        let ss = ScoredSpectrum::new(&s);
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
        let ss = ScoredSpectrum::new(&s);
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
        let ss = ScoredSpectrum::new(&s);
        // Target 100.1 with tol 0.6: all three are within. Closest is 100.0 → rank 1.
        assert_eq!(ss.nearest_peak_rank(100.1, 0.6), Some(1));
    }
}
