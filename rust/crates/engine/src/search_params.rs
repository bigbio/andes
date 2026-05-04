//! Search parameters consumed by candidate enumeration + scoring.

use std::ops::RangeInclusive;

use crate::aa_set::AminoAcidSet;
use crate::enzyme::Enzyme;
use crate::tolerance::{PrecursorTolerance, Tolerance};

#[derive(Debug, Clone)]
pub struct SearchParams {
    pub aa_set: AminoAcidSet,
    pub enzyme: Enzyme,
    pub min_length: u32,
    pub max_length: u32,
    pub max_missed_cleavages: u32,
    pub max_variable_mods_per_peptide: u32,
    /// Precursor mass tolerance (default 20 ppm symmetric).
    pub precursor_tolerance: PrecursorTolerance,
    /// Charges to try for spectra without explicit charge (default 2..=3).
    pub charge_range: RangeInclusive<u8>,
    /// Top-N PSMs to keep per spectrum (default 10).
    pub top_n_psms_per_spectrum: u32,
}

impl SearchParams {
    /// Defaults matching MS-GF+ tryptic search:
    /// - enzyme: Trypsin
    /// - length: 6-40
    /// - missed cleavages: 1
    /// - variable mods per peptide: 3
    /// - precursor tolerance: 20 ppm symmetric
    /// - charge range: 2..=3
    /// - top-N PSMs: 10
    pub fn default_tryptic(aa_set: AminoAcidSet) -> Self {
        Self {
            aa_set,
            enzyme: Enzyme::Trypsin,
            min_length: 6,
            max_length: 40,
            max_missed_cleavages: 1,
            max_variable_mods_per_peptide: 3,
            precursor_tolerance: PrecursorTolerance::symmetric(Tolerance::Ppm(20.0)),
            charge_range: 2..=3,
            top_n_psms_per_spectrum: 10,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aa_set::AminoAcidSetBuilder;

    #[test]
    fn default_tryptic_has_expected_values() {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let params = SearchParams::default_tryptic(aa_set);
        assert_eq!(params.enzyme, Enzyme::Trypsin);
        assert_eq!(params.min_length, 6);
        assert_eq!(params.max_length, 40);
        assert_eq!(params.max_missed_cleavages, 1);
        assert_eq!(params.max_variable_mods_per_peptide, 3);
        assert_eq!(*params.charge_range.start(), 2);
        assert_eq!(*params.charge_range.end(), 3);
        assert_eq!(params.top_n_psms_per_spectrum, 10);
        match params.precursor_tolerance.left {
            Tolerance::Ppm(v) => assert_eq!(v, 20.0),
            _ => panic!("expected Ppm(20.0)"),
        }
    }
}
