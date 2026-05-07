//! Search parameters consumed by candidate enumeration + scoring.

use std::ops::RangeInclusive;

use model::aa_set::AminoAcidSet;
use model::enzyme::Enzyme;
use model::tolerance::{PrecursorTolerance, Tolerance};

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
    /// Isotope offsets to try when matching the precursor mass (default
    /// -1..=2, matching Java MS-GF+'s `-ti -1,2` flag). Each offset is
    /// a unit of `ISOTOPE` (~1.00335 Da) subtracted from the spectrum's
    /// observed neutral mass before comparison.
    pub isotope_error_range: RangeInclusive<i8>,
    /// Top-N PSMs to keep per spectrum (default 10).
    pub top_n_psms_per_spectrum: u32,
    /// Number of Tolerable Termini — mirrors Java's `-ntt N` flag.
    ///
    /// Controls how strictly enzymatic cleavage is enforced at the span boundaries:
    /// - `2` (default): both termini must be enzyme-cleavage sites (strict / fully tryptic).
    /// - `1`: at least one terminus must be a cleavage site (semi-specific). Generates
    ///   semi-tryptic peptides arising from non-canonical proteolysis (e.g., chymotrypsin
    ///   contamination, in-source fragmentation, signal-peptide cleavage).
    /// - `0`: neither terminus needs to be a cleavage site (non-specific). Equivalent to
    ///   using `Enzyme::NonSpecific` — all subsequences within length bounds are emitted.
    ///
    /// Values > 2 are treated identically to 2. Supported values: 0, 1, 2.
    pub num_tolerable_termini: u8,
}

impl SearchParams {
    /// Defaults matching MS-GF+ tryptic search:
    /// - enzyme: Trypsin
    /// - length: 6-40
    /// - missed cleavages: 1
    /// - variable mods per peptide: 3
    /// - precursor tolerance: 20 ppm symmetric
    /// - charge range: 2..=3
    /// - isotope error range: -1..=2 (matches Java's `-ti -1,2` default)
    /// - top-N PSMs: 10
    /// - num_tolerable_termini: 2 (strict tryptic)
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
            isotope_error_range: -1..=2,
            top_n_psms_per_spectrum: 10,
            num_tolerable_termini: 2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;

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
        assert_eq!(*params.isotope_error_range.start(), -1);
        assert_eq!(*params.isotope_error_range.end(), 2);
        assert_eq!(params.top_n_psms_per_spectrum, 10);
        match params.precursor_tolerance.left {
            Tolerance::Ppm(v) => assert_eq!(v, 20.0),
            _ => panic!("expected Ppm(20.0)"),
        }
        assert_eq!(params.num_tolerable_termini, 2);
    }
}
