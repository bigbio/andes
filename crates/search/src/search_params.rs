//! Search parameters consumed by candidate enumeration + scoring.

use std::ops::RangeInclusive;

use model::aa_set::AminoAcidSet;
use model::enzyme::Enzyme;
use model::tolerance::{PrecursorTolerance, Tolerance};

use crate::precursor_cal::PrecursorCalMode;

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
    /// -1..=2). Each offset is a unit of `ISOTOPE` (~1.00335 Da) subtracted
    /// from the spectrum's observed neutral mass before comparison.
    pub isotope_error_range: RangeInclusive<i8>,
    /// Top-N PSMs to keep per spectrum (default 10).
    pub top_n_psms_per_spectrum: u32,
    /// Number of Tolerable Termini.
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
    /// Minimum number of peaks required in an MS2 spectrum to attempt scoring.
    ///
    /// Spectra with fewer peaks than this threshold are skipped entirely.
    /// Default 10.
    pub min_peaks: u32,
    /// Precursor mass calibration mode (Java `-precursorCal`). Default `Off`
    /// (opt-in).
    pub precursor_cal_mode: PrecursorCalMode,
    /// Learned file-wide ppm shift applied to observed neutral masses in the
    /// main pass. Stays 0.0 until the pre-pass calibrator runs.
    pub precursor_mass_shift_ppm: f64,
    /// Full-isolation-window chimeric search (MSFragger-DDA+ style). Default false.
    pub chimeric: bool,
    /// Fallback isolation half-width (Da) used when the mzML lacks
    /// `<isolationWindow>` offsets. Only consulted when `chimeric` is true.
    pub chimeric_isolation_halfwidth_da: f64,
    /// GF-free scoring mode (opt-in, default false). When true, the patented
    /// generating function (SpecEValue DP) is SKIPPED entirely: candidate
    /// selection/ranking stays on `rank_score` (RawScore), PSMs are ordered
    /// for output by `rank_score` (descending), and the GF-derived PIN/TSV
    /// columns (`lnSpecEValue`, `lnEValue`, `lnDeltaSpecEValue`, `DeNovoScore`,
    /// `SpecEValue`, `EValue`) are omitted so Percolator calibrates FDR from
    /// RawScore + the remaining features. The `PsmMatch` GF fields keep their
    /// "not yet computed" sentinels (`spec_e_value = 1.0`, `de_novo_score =
    /// i32::MIN`, `e_value = 1.0`). MUST be bit-identical to the default path
    /// when false.
    pub gf_free: bool,
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
    /// - min_peaks: 10 (matches Java's `-minNumPeaks 10` default)
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
            min_peaks: 10,
            precursor_cal_mode: PrecursorCalMode::Off,
            precursor_mass_shift_ppm: 0.0,
            chimeric: false,
            chimeric_isolation_halfwidth_da: 1.5,
            gf_free: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precursor_cal::PrecursorCalMode;
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
        assert_eq!(params.precursor_cal_mode, PrecursorCalMode::Off);
        assert_eq!(params.precursor_mass_shift_ppm, 0.0);
        assert!(!params.chimeric);
        assert_eq!(params.chimeric_isolation_halfwidth_da, 1.5);
        assert!(!params.gf_free);
    }
}
