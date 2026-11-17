use std::ops::RangeInclusive;

use model::aa_set::AminoAcidSet;
use model::enzyme::Enzyme;
use model::tolerance::{PrecursorTolerance, Tolerance};

use crate::precursor_cal::PrecursorCalMode;

pub const STRONG_SCORE_RETENTION_K: u32 = 25;

pub const FRAGMENT_INDEX_TOP_K: u32 = 100;

pub const FRAGMENT_INDEX_MIN_MATCHED: u16 = 3;

pub const FRAGMENT_INDEX_LOWRES_VOTE_PEAKS: usize = 150;

pub const MMAP_WINDOW_CACHE_MAX_CANDIDATES: usize = 4_000_000;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpenSearch {
    pub delta_min: f64,
    pub delta_max: f64,
    pub top_k: u32,
}

impl OpenSearch {
    pub fn precursor_tolerance(&self) -> PrecursorTolerance {
        PrecursorTolerance::asymmetric(
            Tolerance::Da(self.delta_max.max(1e-6)),
            Tolerance::Da(self.delta_min.abs().max(1e-6)),
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ScoreMode {
    #[default]
    Rank,
    Strong,
}

#[derive(Debug, Clone)]
pub struct SearchParams {
    pub aa_set: AminoAcidSet,
    pub enzyme: Enzyme,
    pub extra_enzymes: Vec<Enzyme>,
    pub min_length: u32,
    pub max_length: u32,
    pub max_missed_cleavages: u32,
    pub max_variable_mods_per_peptide: u32,
    pub max_label_sites: u32,
    pub precursor_tolerance: PrecursorTolerance,
    pub charge_range: RangeInclusive<u8>,
    pub isotope_error_range: RangeInclusive<i8>,
    pub top_n_psms_per_spectrum: u32,
    pub num_tolerable_termini: u8,
    pub min_peaks: u32,
    pub precursor_cal_mode: PrecursorCalMode,
    pub cal_min_spec_keys: usize,
    pub fragment_index_top_k: u32,
    pub precursor_mass_shift_ppm: f64,
    pub chimeric: bool,
    pub dia: bool,
    pub chimeric_isolation_halfwidth_da: f64,
    pub chimeric_max_coisolated: usize,
    pub chimeric_max_kl: f32,
    pub score_mode: ScoreMode,
    pub refine_select_psm_fdr: f64,
    pub require_nxst_sequon: bool,
    pub open_search: Option<OpenSearch>,
}

impl SearchParams {
    pub fn default_tryptic(aa_set: AminoAcidSet) -> Self {
        Self {
            aa_set,
            enzyme: Enzyme::Trypsin,
            extra_enzymes: Vec::new(),
            min_length: 6,
            max_length: 50,
            max_missed_cleavages: 1,
            max_variable_mods_per_peptide: 3,
            max_label_sites: 0,
            precursor_tolerance: PrecursorTolerance::symmetric(Tolerance::Ppm(20.0)),
            charge_range: 2..=3,
            isotope_error_range: -1..=2,
            top_n_psms_per_spectrum: 10,
            num_tolerable_termini: 2,
            min_peaks: 10,
            precursor_cal_mode: PrecursorCalMode::Off,
            cal_min_spec_keys: crate::precursor_cal::constants::MIN_SPECKEYS_FOR_PREPASS,
            fragment_index_top_k: 0,
            precursor_mass_shift_ppm: 0.0,
            chimeric: false,
            dia: false,
            chimeric_isolation_halfwidth_da: 1.5,
            chimeric_max_coisolated: 4,
            chimeric_max_kl: 0.3,
            score_mode: ScoreMode::Rank,
            refine_select_psm_fdr: 0.01,
            require_nxst_sequon: false,
            open_search: None,
        }
    }

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
    use crate::precursor_matching::matches_precursor;
    use model::aa_set::AminoAcidSetBuilder;
    use model::mass::PROTON;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;

    #[test]
    fn hot_loop_retention_cap_widens_only_under_strong() {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.top_n_psms_per_spectrum = 1;
        assert_eq!(params.hot_loop_retention_cap(), 1);
        params.score_mode = ScoreMode::Strong;
        assert_eq!(params.hot_loop_retention_cap(), STRONG_SCORE_RETENTION_K);
        params.top_n_psms_per_spectrum = STRONG_SCORE_RETENTION_K + 5;
        assert_eq!(
            params.hot_loop_retention_cap(),
            STRONG_SCORE_RETENTION_K + 5
        );
    }

    #[test]
    fn default_tryptic_has_expected_values() {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let params = SearchParams::default_tryptic(aa_set);
        assert_eq!(params.enzyme, Enzyme::Trypsin);
        assert_eq!(params.min_length, 6);
        assert_eq!(params.max_length, 50);
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
        assert_eq!(params.chimeric_max_coisolated, 4);
        assert_eq!(params.chimeric_max_kl, 0.3);
        assert_eq!(params.score_mode, ScoreMode::Rank);
        assert_eq!(params.refine_select_psm_fdr, 0.01);
    }

    #[test]
    fn open_search_tolerance_covers_the_requested_offsets() {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let pep = Peptide::from_str("K.PEPTIDEK.R", &aa_set).unwrap();
        let mass = pep.mass();
        let os = OpenSearch {
            delta_min: -100.0,
            delta_max: 300.0,
            top_k: 100,
        };
        let tol = os.precursor_tolerance();
        let spec_at = |delta: f64| {
            let z = 2.0f64;
            Spectrum {
                precursor_mz: (mass + delta + z * PROTON) / z,
                precursor_charge: Some(2),
                ..Default::default()
            }
        };
        for inside in [250.0, 50.0, -50.0, 0.0] {
            assert!(
                matches_precursor(&spec_at(inside), &pep, 2, 0, &tol, 0.0).is_some(),
                "offset {inside} inside -100..300"
            );
        }
        for outside in [350.0, -150.0] {
            assert!(
                matches_precursor(&spec_at(outside), &pep, 2, 0, &tol, 0.0).is_none(),
                "offset {outside} outside -100..300"
            );
        }
    }
}
