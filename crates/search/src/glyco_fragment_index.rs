//! Fragment-ion index for glyco candidate generation (peptide-first).
//!
//! The precursor-mass bucket index (`bucket_index`) that andes inherited from
//! MS-GF+ forces the glyco driver to either fully b/y-score every candidate
//! backbone (the ~40-min brute force) or truncate cheaply by core-Y intensity
//! (which loses weak/absent-core-Y spectra — the measured candidate-generation
//! ceiling). the reference engine avoids both with a FRAGMENT-ion index: peptides are
//! indexed by their fragment masses, so a spectrum's peaks look up the peptides
//! that actually have matching b/y ions in O(peaks), and all candidates can be
//! scored cheaply.
//!
//! This module provides that index. For glyco it enables a PEPTIDE-FIRST query:
//! spectrum peaks → sequon peptides with real b/y support → glycan by
//! subtraction (`glycan = precursor − peptide`), with no backbone enumeration
//! and no lossy core-Y truncation. Because candidates are selected by PEPTIDE
//! evidence (b/y), it works on weak-core-Y spectra where the core-Y selector
//! fails.

use rustc_hash::{FxHashMap, FxHashSet};

use model::peptide::Peptide;
use scoring_crate::scoring::fragment_ions::predict_by_ions;

/// Inverted index: fragment-m/z bin → `(candidate_index, theoretical_mz)`
/// postings. The theoretical m/z is stored so `query` can validate the real
/// fragment tolerance (bins are only a coarse bucket) and count DISTINCT matched
/// theoretical ions per candidate rather than raw postings.
///
/// Hashing uses `rustc-hash` (FxHash) throughout: the per-spectrum query hits
/// the bin map ~3×/peak and the match-dedup structures once per matched ion, so
/// the default SipHash was the dominant glyco-phase cost (profiled ~93%). FxHash
/// on the integer keys is behaviour-identical and much cheaper.
pub struct FragmentIndex {
    bin_width: f64,
    tol: f64,
    /// When set, the per-ion acceptance window is `theo_mz * tol_ppm * 1e-6`
    /// (floored at `PPM_FLOOR_DA`) instead of the fixed `tol` in Da. `bin_width`
    /// is then sized for the widest window in the indexed m/z range so the
    /// `b-1..=b+1` neighbourhood still covers every admissible match.
    tol_ppm: Option<f64>,
    bins: FxHashMap<i64, Vec<(u32, f64)>>,
}

/// Absolute floor (Da) on the ppm acceptance window, so a very light fragment is
/// not held to a sub-millidalton tolerance no instrument delivers.
pub const PPM_FLOOR_DA: f64 = 0.005;
/// Upper m/z the ppm bin width is sized for; a heavier fragment still matches
/// exactly because the distance check is per ion, only the bin neighbourhood
/// might miss it, and b/y ions above this m/z are rare in DDA glycopeptide scans.
const PPM_BIN_MAX_MZ: f64 = 5000.0;

impl FragmentIndex {
    #[inline]
    fn bin_of(mz: f64, bin_width: f64) -> i64 {
        (mz / bin_width).round() as i64
    }

    /// Identity of a theoretical ion for distinct-ion counting: high-resolution
    /// rounded m/z (1 mDa), fine enough to separate distinct b/y ions.
    #[inline]
    fn ion_id(mz: f64) -> i64 {
        (mz * 1000.0).round() as i64
    }

    /// Build the index over `(candidate_index, &Peptide)` entries, indexing each
    /// peptide's b/y ions at charges `1..=max_charge`. `tol` is the fragment match
    /// tolerance in Da; `bin_width` should be ≥ `tol` (bins bucket coarsely, the
    /// stored theoretical m/z is distance-checked at query time). A peptide is
    /// recorded at most once per bin so a bin collision within one peptide does not
    /// bloat the postings.
    ///
    /// CHARGE-AWARE (reviewer feedback): indexing ONLY +1 b/y misses the multiply-
    /// charged backbone fragments of LARGE / HIGH-CHARGE glycopeptides (z4/z5+),
    /// whose b/y ions land at +2/+3 — exactly the bucket peptide-first is meant to
    /// rescue. `max_charge` (clamped 1..=3) indexes those too, so a +2 backbone ion
    /// can select its peptide. `max_charge = 1` reproduces the legacy +1-only index.
    pub fn build<'a>(
        entries: impl IntoIterator<Item = (u32, &'a Peptide)>,
        tol: f64,
        max_charge: u8,
    ) -> Self {
        let bin_width = tol.max(0.001);
        let zmax = max_charge.clamp(1, 3);
        let mut bins: FxHashMap<i64, Vec<(u32, f64)>> = FxHashMap::default();
        for (idx, pep) in entries {
            for ion in predict_by_ions(pep, 1..=zmax) {
                let b = Self::bin_of(ion.mz, bin_width);
                bins.entry(b).or_default().push((idx, ion.mz));
            }
        }
        FragmentIndex {
            bin_width,
            tol,
            tol_ppm: None,
            bins,
        }
    }

    /// [`Self::build`] with a RELATIVE (ppm) acceptance window per ion. This is
    /// the retrieval-tolerance decoupling the roadmap's Stage S1A asks for: the
    /// fixed-Da `tol` the engine passes here is the rank-scoring model's window
    /// (0.5 Da for the low-res models), which on high-resolution data admits b/y
    /// matches ~50x wider than the 20 ppm every glycan-side matcher uses.
    /// `tol_ppm <= 0` falls back to [`Self::build`].
    pub fn build_ppm<'a>(
        entries: impl IntoIterator<Item = (u32, &'a Peptide)>,
        tol_ppm: f64,
        max_charge: u8,
    ) -> Self {
        Self::build_ppm_sized(entries, tol_ppm, max_charge, PPM_BIN_MAX_MZ)
    }

    /// [`Self::build_ppm`] with an explicit m/z the bin width is sized for. The
    /// query walks `ceil(window / bin_width)` neighbouring bins on each side, so a
    /// fragment ABOVE `bin_max_mz` still matches exactly (it just visits more bins);
    /// the parameter only trades bin count against per-bin posting length. Exposed
    /// so a test can force the multi-bin path at ordinary m/z.
    pub fn build_ppm_sized<'a>(
        entries: impl IntoIterator<Item = (u32, &'a Peptide)>,
        tol_ppm: f64,
        max_charge: u8,
        bin_max_mz: f64,
    ) -> Self {
        if tol_ppm.is_nan() || tol_ppm <= 0.0 {
            return Self::build(entries, 0.5, max_charge);
        }
        let widest = (bin_max_mz * tol_ppm * 1e-6).max(PPM_FLOOR_DA);
        let mut idx = Self::build(entries, widest, max_charge);
        idx.tol_ppm = Some(tol_ppm);
        idx
    }

    /// Bins to visit on each side of a peak's own bin so every theoretical ion
    /// within the acceptance window at that m/z is reachable. 1 for the fixed-Da
    /// index (bin_width >= tol) and for ppm queries below the sizing m/z.
    #[inline]
    fn neighbourhood(&self, mz: f64) -> i64 {
        ((self.window_da(mz) / self.bin_width).ceil() as i64).max(1)
    }

    /// Acceptance half-window (Da) for a theoretical ion at `theo` m/z.
    #[inline]
    fn window_da(&self, theo: f64) -> f64 {
        match self.tol_ppm {
            Some(p) => (theo * p * 1e-6).max(PPM_FLOOR_DA),
            None => self.tol,
        }
    }

    /// Return `(candidate_index, matched_ion_count)` for candidates whose
    /// theoretical b/y ions match at least `min_matches` spectrum peaks within
    /// the fragment tolerance, where each observed PEAK explains at most one
    /// theoretical ion and each theoretical ion is explained by at most one peak
    /// (a greedy 1:1 peak↔ion matching). So neither noise clustered near one ion
    /// NOR one peak sitting on a cross-charge m/z coincidence can inflate the
    /// count: the reported count never exceeds `min(#distinct ions, #peaks)`.
    pub fn query(&self, peaks: &[(f64, f32)], min_matches: u32) -> Vec<(u32, u32)> {
        // Greedy 1:1 peak↔ion matching per candidate. Two dedup directions must
        // both hold or the count inflates:
        //   (a) several peaks near ONE theoretical ion  → count it once
        //   (b) ONE peak near SEVERAL theoretical ions  → count it once
        // Direction (b) is why multi-charge indexing needs the peak-side guard:
        // a b_i⁺ and a y_j²⁺ of the SAME peptide can fall within tol of a single
        // observed peak, and counting both would let a peptide pass `min_matches`
        // on fewer real peaks than the threshold (an FDR-relevant inflation).
        //
        // `used_ion` enforces (a): a theoretical ion (keyed by candidate+ion_id)
        // is consumed once. `matched_this_peak` enforces (b): within one peak's
        // scan, a candidate accepts at most one of its ions. It is cleared (not
        // reallocated) per peak, so the profiled per-candidate-HashSet alloc the
        // old single-charge path removed does not return.
        // (A `(u32, i64)` key — not a packed `u64` — keeps the ion id lossless for
        // any accepted params: `--max-length` is uncapped and mod deltas reach
        // ±5000 Da, so `ion_id = round(mz·1000)` is not guaranteed to fit in 32 bits.)
        let mut used_ion: FxHashSet<(u32, i64)> = FxHashSet::default();
        let mut matched_this_peak: FxHashSet<u32> = FxHashSet::default();
        let mut count: FxHashMap<u32, u32> = FxHashMap::default();
        for &(mz, _) in peaks {
            matched_this_peak.clear();
            let b = Self::bin_of(mz, self.bin_width);
            let k = self.neighbourhood(mz);
            for nb in (b - k)..=(b + k) {
                if let Some(v) = self.bins.get(&nb) {
                    for &(idx, theo) in v {
                        // Accept this (peak, ion) match only if the peak has not
                        // yet been consumed by this candidate AND the ion has not
                        // yet been consumed by any peak. Check both before either
                        // insert, so a peak whose first-seen ion is already used
                        // can still match a still-free ion of the same candidate.
                        let ion_key = (idx, Self::ion_id(theo));
                        if (mz - theo).abs() <= self.window_da(theo)
                            && !matched_this_peak.contains(&idx)
                            && !used_ion.contains(&ion_key)
                        {
                            matched_this_peak.insert(idx);
                            used_ion.insert(ion_key);
                            *count.entry(idx).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
        count
            .into_iter()
            .filter(|&(_, c)| c >= min_matches)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    /// A fragment above the bin-sizing m/z must still match within its ppm window:
    /// the query widens its bin neighbourhood instead of silently missing it.
    /// Forced at ordinary m/z by sizing the bins for 300 m/z.
    #[test]
    fn ppm_query_reaches_matches_beyond_the_bin_sizing_mz() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDENK.R", &aa).expect("valid");
        let theo = predict_by_ions(&pep, 1..=1)
            .iter()
            .map(|i| i.mz)
            .fold(0.0f64, f64::max);
        assert!(
            theo > 300.0,
            "test premise: heaviest ion above the sizing m/z"
        );
        let idx = FragmentIndex::build_ppm_sized([(0u32, &pep)], 20.0, 1, 300.0);
        assert!(
            idx.neighbourhood(theo) > 1,
            "must need more than one bin per side here"
        );
        // +18 ppm: inside the window, but several bins away at this bin width.
        let peaks = vec![(theo * (1.0 + 18e-6), 100.0f32)];
        assert_eq!(idx.query(&peaks, 1), vec![(0, 1)]);
        // and still rejects +30 ppm
        let far = vec![(theo * (1.0 + 30e-6), 100.0f32)];
        assert!(idx.query(&far, 1).is_empty());
    }

    /// Stage S1A: a ppm index rejects a peak the 0.5 Da index accepts, and still
    /// accepts a peak inside its own window, for the same peptide.
    #[test]
    fn ppm_window_is_per_ion_and_much_narrower_than_the_da_window() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDENK.R", &aa).expect("valid");
        let ions = predict_by_ions(&pep, 1..=1);
        let theo = ions.iter().map(|i| i.mz).fold(0.0f64, f64::max); // heaviest y ion
        let da = FragmentIndex::build([(0u32, &pep)], 0.5, 1);
        let ppm = FragmentIndex::build_ppm([(0u32, &pep)], 20.0, 1);
        let off_by_0_3 = vec![(theo + 0.3, 100.0f32)];
        assert_eq!(
            da.query(&off_by_0_3, 1),
            vec![(0, 1)],
            "0.5 Da index must accept +0.3 Da"
        );
        assert!(
            ppm.query(&off_by_0_3, 1).is_empty(),
            "20 ppm index must reject +0.3 Da"
        );
        let inside = vec![(theo * (1.0 + 10e-6), 100.0f32)];
        assert_eq!(
            ppm.query(&inside, 1),
            vec![(0, 1)],
            "20 ppm index must accept +10 ppm"
        );
        // Floor: a light fragment at +4 mDa is inside the 5 mDa floor.
        let light = ions.iter().map(|i| i.mz).fold(f64::MAX, f64::min);
        let near_light = vec![(light + 0.004, 100.0f32)];
        assert_eq!(
            ppm.query(&near_light, 1),
            vec![(0, 1)],
            "floor must admit +4 mDa on a light ion"
        );
        // tol_ppm <= 0 degrades to the Da index.
        let zero = FragmentIndex::build_ppm([(0u32, &pep)], 0.0, 1);
        assert_eq!(zero.query(&off_by_0_3, 1), vec![(0, 1)]);
    }

    use super::*;
    use model::aa_set::AminoAcidSetBuilder;

    fn aa() -> model::aa_set::AminoAcidSet {
        AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .build()
            .unwrap()
    }

    /// The index must select the peptide whose b/y ions populate the spectrum,
    /// and rank it above an unrelated peptide.
    #[test]
    fn fragment_index_selects_peptide_with_matching_by_ions() {
        let aa = aa();
        let pep_a = Peptide::from_str("K.PEPTIDESK.R", &aa).expect("valid a");
        let pep_b = Peptide::from_str("K.ELVISLIVESR.K", &aa).expect("valid b");
        let idx = FragmentIndex::build([(0u32, &pep_a), (1u32, &pep_b)], 0.02, 1);

        // Spectrum = pep_a's own singly-charged b/y ions.
        let peaks: Vec<(f64, f32)> = predict_by_ions(&pep_a, 1..=1)
            .iter()
            .map(|i| (i.mz, 100.0))
            .collect();

        let hits = idx.query(&peaks, 3);
        let a = hits
            .iter()
            .find(|&&(i, _)| i == 0)
            .map(|&(_, c)| c)
            .unwrap_or(0);
        let b = hits
            .iter()
            .find(|&&(i, _)| i == 1)
            .map(|&(_, c)| c)
            .unwrap_or(0);
        assert!(a >= 3, "pep_a must match its own ladder (got {a})");
        assert!(a > b, "pep_a ({a}) must outscore the unrelated pep_b ({b})");
    }

    /// Codex re-review #2: match counts must be tolerance-checked DISTINCT ions.
    /// A peak in the bin neighbourhood but OUTSIDE the fragment tolerance must
    /// not count, and several peaks clustered around ONE theoretical ion must
    /// count that ion only once.
    #[test]
    fn fragment_index_counts_tolerance_checked_distinct_ions() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDESK.R", &aa).unwrap();
        let tol = 0.02;
        let idx = FragmentIndex::build([(0u32, &pep)], tol, 1);
        let ions = predict_by_ions(&pep, 1..=1);

        // Take the first real ion; add: (a) an in-tolerance duplicate, (b) an
        // out-of-tolerance neighbour (~0.5 Da off, same 0.02 bin neighbourhood
        // it is NOT, but well outside tol) — neither should raise the count above
        // the count of distinct real ions.
        let first = ions[0].mz;
        let mut peaks: Vec<(f64, f32)> = vec![
            (first, 100.0),
            (first + 0.005, 100.0), // within tol → same ion, must not double-count
            (first + 0.5, 100.0),   // far off → must not count
        ];
        // Add a couple more genuine ions so the peptide can pass a threshold of 3.
        peaks.push((ions[1].mz, 100.0));
        peaks.push((ions[2].mz, 100.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let hits = idx.query(&peaks, 1);
        let c = hits
            .iter()
            .find(|&&(i, _)| i == 0)
            .map(|&(_, c)| c)
            .unwrap_or(0);
        // Exactly 3 distinct ions matched (ions[0], [1], [2]); the duplicate and
        // the far peak add nothing.
        assert_eq!(c, 3, "must count 3 distinct tolerance-valid ions, got {c}");
        // And the far-off peak alone must not let it pass a threshold of 4.
        assert!(
            idx.query(&peaks, 4).is_empty(),
            "must not inflate past distinct-ion count"
        );
    }

    /// Charge-aware: a spectrum of a peptide's +2 b/y ions selects it only when the
    /// index was built with max_charge >= 2 — the high-charge glycopeptide bucket.
    #[test]
    fn fragment_index_is_charge_aware() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDESKELVISR.R", &aa).unwrap();
        // Spectrum = the peptide's DOUBLY-charged b/y ions.
        let peaks: Vec<(f64, f32)> = predict_by_ions(&pep, 2..=2)
            .iter()
            .map(|i| (i.mz, 100.0))
            .collect();
        let cnt = |z: u8| {
            FragmentIndex::build([(0u32, &pep)], 0.02, z)
                .query(&peaks, 1)
                .iter()
                .find(|&&(i, _)| i == 0)
                .map(|&(_, c)| c)
                .unwrap_or(0)
        };
        let c1 = cnt(1);
        let c2 = cnt(2);
        // The +2 index selects the peptide from its +2 ladder, and finds strictly
        // MORE +2 ions than the +1-only index (which can only match the rare +2/+1
        // m/z coincidence). This is the high-charge-glycopeptide capability.
        assert!(
            c2 >= 3,
            "+2 index must select the peptide from its +2 ladder (got {c2})"
        );
        assert!(
            c2 > c1,
            "+2 index must match more +2 ions than +1-only ({c2} > {c1})"
        );
    }

    /// Regression (adversarial review): with a multi-charge index, one observed
    /// peak can fall within tolerance of TWO theoretical ions of the SAME peptide
    /// (e.g. a bᵢ⁺ and a yⱼ²⁺ coinciding). The distinct-ion count must NOT credit
    /// both — the matched count can never exceed the number of observed peaks, or
    /// a peptide would pass `min_matches` on fewer real peaks than the threshold.
    #[test]
    fn fragment_index_count_never_exceeds_peak_count() {
        let aa = aa();
        // Short peptide → dense, low-m/z ladder where +1 and +2 ions crowd and
        // coincide within a 0.02 Da tolerance under a max_charge=2 index.
        let pep = Peptide::from_str("K.SNSSGGVG.R", &aa).unwrap();
        let idx = FragmentIndex::build([(0u32, &pep)], 0.02, 2);

        // Feed the union of the peptide's +1 and +2 ions, then collapse any
        // near-coincident theoretical m/z into a SINGLE observed peak (what a real
        // spectrum would show): distinct observed peaks < distinct theoretical ions.
        let mut theo: Vec<f64> = predict_by_ions(&pep, 1..=2).iter().map(|i| i.mz).collect();
        theo.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut peaks: Vec<(f64, f32)> = Vec::new();
        for m in theo {
            if peaks
                .last()
                .map(|&(p, _)| (m - p).abs() > 0.02)
                .unwrap_or(true)
            {
                peaks.push((m, 100.0));
            }
        }

        let n_peaks = peaks.len() as u32;
        let c = idx
            .query(&peaks, 1)
            .iter()
            .find(|&&(i, _)| i == 0)
            .map(|&(_, c)| c)
            .unwrap_or(0);
        assert!(
            c <= n_peaks,
            "matched count {c} must not exceed observed peak count {n_peaks} (cross-charge inflation)"
        );
        // And the peptide must NOT pass a threshold above the real peak count.
        assert!(
            idx.query(&peaks, n_peaks + 1).is_empty(),
            "must not pass a min_matches above the number of observed peaks"
        );
    }

    /// A spectrum of pure noise must select nothing at a sane threshold.
    #[test]
    fn fragment_index_rejects_noise() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDESK.R", &aa).unwrap();
        let idx = FragmentIndex::build([(0u32, &pep)], 0.02, 1);
        let noise: Vec<(f64, f32)> = vec![(123.456, 1.0), (777.111, 1.0), (1500.9, 1.0)];
        let hits = idx.query(&noise, 3);
        assert!(
            hits.is_empty(),
            "noise must not select the peptide, got {hits:?}"
        );
    }
}
