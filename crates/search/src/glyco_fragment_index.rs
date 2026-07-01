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

use std::collections::{HashMap, HashSet};

use model::peptide::Peptide;
use scoring_crate::scoring::fragment_ions::predict_by_ions;

/// Inverted index: fragment-m/z bin → candidate indices with a b/y ion there.
pub struct FragmentIndex {
    bin_width: f64,
    bins: HashMap<i64, Vec<u32>>,
}

impl FragmentIndex {
    #[inline]
    fn bin_of(mz: f64, bin_width: f64) -> i64 {
        (mz / bin_width).round() as i64
    }

    /// Build the index over `(candidate_index, &Peptide)` entries, indexing each
    /// peptide's singly-charged b/y ions. `bin_width` is the m/z bin size (Da);
    /// pick it near the fragment tolerance (e.g. 0.02 for 20 ppm at ~1 kDa). A
    /// peptide is recorded at most once per bin so a bin collision within one
    /// peptide does not inflate its match count.
    pub fn build<'a>(
        entries: impl IntoIterator<Item = (u32, &'a Peptide)>,
        bin_width: f64,
    ) -> Self {
        let mut bins: HashMap<i64, Vec<u32>> = HashMap::new();
        for (idx, pep) in entries {
            let mut seen: HashSet<i64> = HashSet::new();
            for ion in predict_by_ions(pep, 1..=1) {
                let b = Self::bin_of(ion.mz, bin_width);
                if seen.insert(b) {
                    bins.entry(b).or_default().push(idx);
                }
            }
        }
        FragmentIndex { bin_width, bins }
    }

    /// Return `(candidate_index, match_count)` for candidates whose b/y ions
    /// match at least `min_matches` DISTINCT spectrum peaks. Each peak is matched
    /// against its bin ±1 (≈ ±1.5·bin_width tolerance) and contributes at most
    /// one count per candidate.
    pub fn query(&self, peaks: &[(f64, f32)], min_matches: u32) -> Vec<(u32, u32)> {
        let mut counts: HashMap<u32, u32> = HashMap::new();
        for &(mz, _) in peaks {
            let b = Self::bin_of(mz, self.bin_width);
            let mut hit: HashSet<u32> = HashSet::new();
            for nb in [b - 1, b, b + 1] {
                if let Some(v) = self.bins.get(&nb) {
                    for &idx in v {
                        hit.insert(idx);
                    }
                }
            }
            for idx in hit {
                *counts.entry(idx).or_default() += 1;
            }
        }
        counts
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
