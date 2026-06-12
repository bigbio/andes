//! Search parameters consumed by candidate enumeration + scoring.

use std::ops::RangeInclusive;

use model::aa_set::AminoAcidSet;
use model::enzyme::Enzyme;
use model::tolerance::{PrecursorTolerance, Tolerance};

use crate::precursor_cal::PrecursorCalMode;

/// Rank-score pool retained in the hot loop when `--score strong` and user
/// `top_n` is smaller. Runners-up survive to `fill_post_topn` / strong re-rank.
///
/// Widened 10 -> 25 after the Phase-V gate: strong improved FDP but lost
/// PSMs@1%, and with K=10 a true peptide that rank-scores outside the top-10
/// can never be promoted by the strong re-rank (the pool is rank-gated). 25
/// gives the strong score a larger pool to recover those top-1 flips. Strong
/// mode only — default `rank` retention is unchanged. Re-run the Astral/TMT
/// A/B after changing this; raise further if it keeps recovering PSMs.
pub const STRONG_SCORE_RETENTION_K: u32 = 25;

/// Primary ranking mode for candidate selection and PIN `RawScore` emission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ScoreMode {
    /// Rank by inherited `node + cleavage + edge` (`rank_score`). Default; byte-identical path.
    #[default]
    Rank,
    /// Rank by fused `strong_score` and emit it as PIN `RawScore`.
    Strong,
}

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
    /// Precursor mass calibration mode (`auto`, `on`, `off`). Default `Off`
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
    /// Ranking / RawScore source: `Rank` (default) or `Strong` (S3 fused score).
    pub score_mode: ScoreMode,
}

impl SearchParams {
    /// Defaults for a standard tryptic search (Kim et al., Nat Commun 5:5277, 2014):
    /// - enzyme: Trypsin
    /// - length: 6-40
    /// - missed cleavages: 1
    /// - variable mods per peptide: 3
    /// - precursor tolerance: 20 ppm symmetric
    /// - charge range: 2..=3
    /// - isotope error range: -1..=2
    /// - top-N PSMs: 10
    /// - num_tolerable_termini: 2 (strict tryptic)
    /// - min_peaks: 10
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
            score_mode: ScoreMode::Rank,
        }
    }

    /// Top-N cap for `TopNQueue` retention in the candidate hot loop.
    /// Under `Strong`, widens to at least [`STRONG_SCORE_RETENTION_K`] so the
    /// post-loop strong re-rank can promote runners-up; PIN emission still
    /// trims to `top_n_psms_per_spectrum` afterward.
    pub fn hot_loop_retention_cap(&self) -> u32 {
        if self.score_mode == ScoreMode::Strong {
            self.top_n_psms_per_spectrum.max(STRONG_SCORE_RETENTION_K)
        } else {
            self.top_n_psms_per_spectrum
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precursor_cal::PrecursorCalMode;
    use model::aa_set::AminoAcidSetBuilder;

    #[test]
    fn hot_loop_retention_cap_widens_only_under_strong() {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.top_n_psms_per_spectrum = 1;
        assert_eq!(params.hot_loop_retention_cap(), 1);
        params.score_mode = ScoreMode::Strong;
        assert_eq!(
            params.hot_loop_retention_cap(),
            STRONG_SCORE_RETENTION_K
        );
        // User top_n larger than K → user value wins.
        params.top_n_psms_per_spectrum = STRONG_SCORE_RETENTION_K + 5;
        assert_eq!(params.hot_loop_retention_cap(), STRONG_SCORE_RETENTION_K + 5);
    }

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
        assert_eq!(params.score_mode, ScoreMode::Rank);
    }
}
