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
    bins: FxHashMap<i64, Vec<(u32, f64)>>,
}

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
    /// peptide's singly-charged b/y ions. `tol` is the fragment match tolerance
    /// in Da; `bin_width` should be ≥ `tol` (bins bucket coarsely, the stored
    /// theoretical m/z is distance-checked at query time). A peptide is recorded
    /// at most once per bin so a bin collision within one peptide does not bloat
    /// the postings.
    pub fn build<'a>(
        entries: impl IntoIterator<Item = (u32, &'a Peptide)>,
        tol: f64,
    ) -> Self {
        let bin_width = tol.max(0.001);
        let mut bins: FxHashMap<i64, Vec<(u32, f64)>> = FxHashMap::default();
        for (idx, pep) in entries {
            for ion in predict_by_ions(pep, 1..=1) {
                let b = Self::bin_of(ion.mz, bin_width);
                bins.entry(b).or_default().push((idx, ion.mz));
            }
        }
        FragmentIndex { bin_width, tol, bins }
    }

    /// Return `(candidate_index, matched_ion_count)` for candidates whose
    /// theoretical b/y ions match at least `min_matches` DISTINCT spectrum peaks
    /// within the fragment tolerance. Each theoretical ion is counted at most
    /// once per candidate (a peak within tol of the ion "explains" it), so
    /// noise clustered near one ion cannot inflate the count.
    pub fn query(&self, peaks: &[(f64, f32)], min_matches: u32) -> Vec<(u32, u32)> {
        // Count DISTINCT (candidate, theoretical-ion) matches. `seen` dedups a
        // theoretical ion matched by several peaks, keyed by the exact
        // `(candidate_idx, ion_id)` pair; `count` tallies the first sighting per
        // candidate. This replaces the previous `HashMap<u32, HashSet<i64>>`,
        // which allocated a HashSet per matched candidate on every query — the
        // profiled hotspot. Same distinct-ion semantics, no per-candidate alloc.
        // (A `(u32, i64)` key — not a packed `u64` — keeps the ion id lossless for
        // any accepted params: `--max-length` is uncapped and mod deltas reach
        // ±5000 Da, so `ion_id = round(mz·1000)` is not guaranteed to fit in 32 bits.)
        let mut seen: FxHashSet<(u32, i64)> = FxHashSet::default();
        let mut count: FxHashMap<u32, u32> = FxHashMap::default();
        for &(mz, _) in peaks {
            let b = Self::bin_of(mz, self.bin_width);
            for nb in [b - 1, b, b + 1] {
                if let Some(v) = self.bins.get(&nb) {
                    for &(idx, theo) in v {
                        if (mz - theo).abs() <= self.tol {
                            if seen.insert((idx, Self::ion_id(theo))) {
                                *count.entry(idx).or_insert(0) += 1;
                            }
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
        let idx = FragmentIndex::build([(0u32, &pep_a), (1u32, &pep_b)], 0.02);

        // Spectrum = pep_a's own singly-charged b/y ions.
        let peaks: Vec<(f64, f32)> = predict_by_ions(&pep_a, 1..=1)
            .iter()
            .map(|i| (i.mz, 100.0))
            .collect();

        let hits = idx.query(&peaks, 3);
        let a = hits.iter().find(|&&(i, _)| i == 0).map(|&(_, c)| c).unwrap_or(0);
        let b = hits.iter().find(|&&(i, _)| i == 1).map(|&(_, c)| c).unwrap_or(0);
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
        let idx = FragmentIndex::build([(0u32, &pep)], tol);
        let ions = predict_by_ions(&pep, 1..=1);

        // Take the first real ion; add: (a) an in-tolerance duplicate, (b) an
        // out-of-tolerance neighbour (~0.5 Da off, same 0.02 bin neighbourhood
        // it is NOT, but well outside tol) — neither should raise the count above
        // the count of distinct real ions.
        let first = ions[0].mz;
        let mut peaks: Vec<(f64, f32)> = vec![
            (first, 100.0),
            (first + 0.005, 100.0),  // within tol → same ion, must not double-count
            (first + 0.5, 100.0),    // far off → must not count
        ];
        // Add a couple more genuine ions so the peptide can pass a threshold of 3.
        peaks.push((ions[1].mz, 100.0));
        peaks.push((ions[2].mz, 100.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let hits = idx.query(&peaks, 1);
        let c = hits.iter().find(|&&(i, _)| i == 0).map(|&(_, c)| c).unwrap_or(0);
        // Exactly 3 distinct ions matched (ions[0], [1], [2]); the duplicate and
        // the far peak add nothing.
        assert_eq!(c, 3, "must count 3 distinct tolerance-valid ions, got {c}");
        // And the far-off peak alone must not let it pass a threshold of 4.
        assert!(idx.query(&peaks, 4).is_empty(), "must not inflate past distinct-ion count");
    }

    /// A spectrum of pure noise must select nothing at a sane threshold.
    #[test]
    fn fragment_index_rejects_noise() {
        let aa = aa();
        let pep = Peptide::from_str("K.PEPTIDESK.R", &aa).unwrap();
        let idx = FragmentIndex::build([(0u32, &pep)], 0.02);
        let noise: Vec<(f64, f32)> = vec![(123.456, 1.0), (777.111, 1.0), (1500.9, 1.0)];
        let hits = idx.query(&noise, 3);
        assert!(hits.is_empty(), "noise must not select the peptide, got {hits:?}");
    }
}
