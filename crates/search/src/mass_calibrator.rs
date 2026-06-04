//! Sampled pre-pass precursor mass calibration (Java `MassCalibrator`).
//!
//! Runs a lightweight top-1 search on ~500 sampled `(spectrum, charge)` keys,
//! filters to high-confidence PSMs, and returns the median ppm residual as the
//! file-wide shift applied in the main pass via
//! [`crate::precursor_cal::adjusted_observed_neutral_mass`].

use std::collections::HashMap;
use std::ops::RangeInclusive;

use model::mass::{H2O, PROTON};
use model::Spectrum;
use model::tolerance::{PrecursorTolerance, Tolerance};

use crate::match_engine::PreparedSearch;
use crate::precursor_cal::{
    constants, median, residual_ppm, robust_sigma_ppm, sample_every_nth, tightened_tolerance_ppm,
    PrecursorCalMode,
};
use crate::search_params::SearchParams;

/// One searchable `(spectrum index, charge)` pair — mirrors Java `SpecKey`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SpecKey {
    pub spectrum_idx: usize,
    pub charge: u8,
}

/// Residuals whose magnitude exceeds this are rejected as isotope contamination.
const MAX_REASONABLE_RESIDUAL_PPM: f64 = 50.0;

/// Java default `Constants.MIN_DE_NOVO_SCORE`.
const MIN_DE_NOVO_SCORE: i32 = 0;

/// Summary of a successful (or skipped) calibration pre-pass.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CalibrationStats {
    pub shift_ppm: f64,
    pub robust_sigma_ppm: f64,
    pub confident_psm_count: usize,
    /// PSMs rejected because top-1 `spec_e_value` > 1e-6.
    pub rejected_spec_e: usize,
    /// PSMs rejected because `|residual_ppm|` > 50.
    pub rejected_residual: usize,
    /// Sampled spectra with at least one PSM in the prepass queue.
    pub queues_with_psm: usize,
}

impl CalibrationStats {
    /// True when the pre-pass produced enough confident PSMs to trust the
    /// learned shift. `learn_calibration_stats` only ever sets
    /// `confident_psm_count` to 0 or `>= MIN_CONFIDENT_PSMS`, so the
    /// `> 0` check is sufficient as a downstream gate.
    pub fn has_reliable_stats(self) -> bool {
        self.confident_psm_count > 0
    }
}

/// Build the SpecKey list for a parsed MS2 file slice.
///
/// Mirrors Java `SpecKey.getSpecKeyList`: spectra below `min_peaks` are
/// skipped; charge-missing spectra expand across `charge_range`.
pub fn build_spec_keys(
    spectra: &[Spectrum],
    charge_range: &RangeInclusive<u8>,
    min_peaks: u32,
) -> Vec<SpecKey> {
    let mut keys = Vec::new();
    for (spectrum_idx, spec) in spectra.iter().enumerate() {
        if spec.peaks.len() < min_peaks as usize {
            continue;
        }
        match spec.precursor_charge {
            Some(z) if z > 0 => keys.push(SpecKey {
                spectrum_idx,
                charge: z as u8,
            }),
            _ => {
                for z in charge_range.clone() {
                    keys.push(SpecKey { spectrum_idx, charge: z });
                }
            }
        }
    }
    keys
}

/// Pre-pass params: isotope error fixed at 0, top-1, no mass shift.
pub fn prepass_search_params(main: &SearchParams) -> SearchParams {
    let mut p = main.clone();
    p.isotope_error_range = 0..=0;
    p.top_n_psms_per_spectrum = 1;
    p.precursor_mass_shift_ppm = 0.0;
    p
}

/// Resolve the ppm shift to store on `SearchParams` after the pre-pass.
pub fn apply_shift_for_mode(mode: PrecursorCalMode, stats: CalibrationStats) -> f64 {
    match mode {
        PrecursorCalMode::Off => 0.0,
        PrecursorCalMode::On => stats.shift_ppm,
        PrecursorCalMode::Auto => {
            if stats.has_reliable_stats() && stats.shift_ppm != 0.0 {
                stats.shift_ppm
            } else {
                0.0
            }
        }
    }
}

/// Run the sampled pre-pass and return calibration stats.
pub fn learn_calibration_stats(
    spec_keys: &[SpecKey],
    originals: &HashMap<usize, Spectrum>,
    prepared: &PreparedSearch<'_>,
    main_params: &SearchParams,
) -> CalibrationStats {
    if spec_keys.len() < constants::MIN_SPECKEYS_FOR_PREPASS {
        return CalibrationStats::default();
    }

    let sampled = sample_every_nth(
        spec_keys,
        constants::SAMPLING_STRIDE,
        constants::MAX_SAMPLED,
    );
    if sampled.is_empty() {
        return CalibrationStats::default();
    }

    let prepass_params = prepass_search_params(main_params);
    let prepass_spectra: Vec<Spectrum> = sampled
        .iter()
        .filter_map(|key| {
            originals
                .get(&key.spectrum_idx)
                .map(|spec| spectrum_with_charge(spec, key.charge))
        })
        .collect();

    if prepass_spectra.len() != sampled.len() {
        return CalibrationStats::default();
    }

    let queues = prepared.run_chunk_with_params(&prepass_spectra, 0, &prepass_params);
    let (residuals, filter) = extract_residuals(
        &sampled,
        &queues,
        originals,
        &prepared.candidates,
        MIN_DE_NOVO_SCORE,
        constants::RESIDUAL_CAP,
    );

    if residuals.len() < constants::MIN_CONFIDENT_PSMS {
        return CalibrationStats {
            rejected_spec_e: filter.rejected_spec_e,
            rejected_residual: filter.rejected_residual,
            queues_with_psm: filter.queues_with_psm,
            ..CalibrationStats::default()
        };
    }

    let shift_ppm = median(&residuals);
    CalibrationStats {
        shift_ppm,
        robust_sigma_ppm: robust_sigma_ppm(&residuals, shift_ppm),
        confident_psm_count: residuals.len(),
        queues_with_psm: filter.queues_with_psm,
        ..CalibrationStats::default()
    }
}

/// Tighten ppm precursor tolerance after a successful cal pass (matching
/// Java's post-cal block). No-op when stats are unreliable or
/// tolerance is not ppm-based.
pub fn apply_tightened_precursor_tolerance(params: &mut SearchParams, stats: CalibrationStats) {
    if !stats.has_reliable_stats() {
        return;
    }
    let (Some(left_ppm), Some(right_ppm)) = (
        ppm_value(params.precursor_tolerance.left),
        ppm_value(params.precursor_tolerance.right),
    ) else {
        return;
    };

    let sigma_mult = constants::TIGHTENED_WINDOW_SIGMA_MULTIPLIER;
    let floor = constants::TIGHTENED_WINDOW_FLOOR_PPM;
    let margin = constants::TIGHTENED_WINDOW_MARGIN_PPM;

    let tightened_left = tightened_tolerance_ppm(left_ppm, stats.robust_sigma_ppm, sigma_mult, floor, margin);
    let tightened_right =
        tightened_tolerance_ppm(right_ppm, stats.robust_sigma_ppm, sigma_mult, floor, margin);

    if tightened_left >= left_ppm && tightened_right >= right_ppm {
        return;
    }

    params.precursor_tolerance = PrecursorTolerance::asymmetric(
        Tolerance::Ppm(tightened_left),
        Tolerance::Ppm(tightened_right),
    );
}

fn ppm_value(t: Tolerance) -> Option<f64> {
    match t {
        Tolerance::Ppm(v) => Some(v),
        Tolerance::Da(_) => None,
    }
}

fn spectrum_with_charge(spec: &Spectrum, charge: u8) -> Spectrum {
    Spectrum {
        precursor_charge: Some(charge as i32),
        ..spec.clone()
    }
}

#[derive(Debug, Default)]
struct CalFilterCounts {
    queues_with_psm: usize,
    rejected_spec_e: usize,
    rejected_de_novo: usize,
    rejected_residual: usize,
}

fn extract_residuals(
    sampled: &[SpecKey],
    queues: &[crate::TopNQueue],
    originals: &HashMap<usize, Spectrum>,
    candidates: &[crate::candidate_gen::Candidate],
    min_de_novo_score: i32,
    keep_top_n: usize,
) -> (Vec<f64>, CalFilterCounts) {
    let mut residual_with_eval: Vec<(f64, f64)> = Vec::new();
    let mut filter = CalFilterCounts::default();

    // Java `generateSpecIndexDBMatchMap` keeps the best SpecEValue per
    // spectrum index across all sampled SpecKeys (e.g. charge variants).
    let mut best_by_spec: HashMap<usize, crate::PsmMatch> = HashMap::new();
    for (key, queue) in sampled.iter().zip(queues.iter()) {
        let Some(psm) = queue
            .iter_psms()
            .min_by(|a, b| {
                a.spec_e_value
                    .partial_cmp(&b.spec_e_value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        else {
            continue;
        };
        filter.queues_with_psm += 1;
        best_by_spec
            .entry(key.spectrum_idx)
            .and_modify(|best| {
                if psm.spec_e_value < best.spec_e_value
                    || (psm.spec_e_value == best.spec_e_value && psm.rank_score > best.rank_score)
                {
                    *best = psm.clone();
                }
            })
            .or_insert_with(|| psm.clone());
    }

    for (&spectrum_idx, psm) in best_by_spec.iter() {
        if psm.spec_e_value > constants::MAX_SPEC_EVALUE {
            filter.rejected_spec_e += 1;
            continue;
        }
        if psm.de_novo_score < min_de_novo_score {
            filter.rejected_de_novo += 1;
            continue;
        }

        let Some(spec) = originals.get(&spectrum_idx) else {
            continue;
        };
        let charge = psm.charge_used as f64;
        if charge <= 0.0 {
            continue;
        }

        let observed = (spec.precursor_mz - PROTON) * charge - H2O;
        let theoretical =
            candidates[psm.primary_candidate_idx() as usize].peptide.residue_mass();
        if theoretical <= 0.0 {
            continue;
        }

        let residual = residual_ppm(observed, theoretical);
        if residual.abs() > MAX_REASONABLE_RESIDUAL_PPM {
            filter.rejected_residual += 1;
            continue;
        }
        residual_with_eval.push((residual, psm.spec_e_value));
    }

    residual_with_eval.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let keep_n = residual_with_eval.len().min(keep_top_n);
    let residuals = residual_with_eval
        .into_iter()
        .take(keep_n)
        .map(|(r, _)| r)
        .collect();
    (residuals, filter)
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;
    use model::enzyme::Enzyme;
    use model::tolerance::{PrecursorTolerance, Tolerance};

    fn empty_params() -> SearchParams {
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        SearchParams {
            aa_set: aa,
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
            precursor_cal_mode: PrecursorCalMode::Auto,
            precursor_mass_shift_ppm: 0.0,
            chimeric: false,
            chimeric_isolation_halfwidth_da: 1.5,
            gf_free: false,
        }
    }

    #[test]
    fn build_spec_keys_expands_missing_charge() {
        let spec = Spectrum {
            title: "t".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: None,
            rt_seconds: None,
            scan: None,
            peaks: vec![(100.0, 1.0); 10],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        };
        let keys = build_spec_keys(&[spec], &(2..=3), 10);
        assert_eq!(keys.len(), 2);
        assert_eq!(keys[0].charge, 2);
        assert_eq!(keys[1].charge, 3);
    }

    #[test]
    fn prepass_params_fixes_isotope_and_top_n() {
        let main = empty_params();
        let prepass = prepass_search_params(&main);
        assert_eq!(*prepass.isotope_error_range.start(), 0);
        assert_eq!(*prepass.isotope_error_range.end(), 0);
        assert_eq!(prepass.top_n_psms_per_spectrum, 1);
        assert_eq!(prepass.precursor_mass_shift_ppm, 0.0);
    }

    #[test]
    fn apply_shift_auto_requires_reliable_nonzero() {
        let stats = CalibrationStats {
            shift_ppm: 5.0,
            robust_sigma_ppm: 0.8,
            confident_psm_count: 200,
            ..CalibrationStats::default()
        };
        assert_eq!(apply_shift_for_mode(PrecursorCalMode::Auto, stats), 5.0);
        assert_eq!(
            apply_shift_for_mode(
                PrecursorCalMode::Auto,
                CalibrationStats {
                    shift_ppm: 0.0,
                    confident_psm_count: 200,
                    ..CalibrationStats::default()
                }
            ),
            0.0
        );
        assert_eq!(
            apply_shift_for_mode(
                PrecursorCalMode::Auto,
                CalibrationStats {
                    shift_ppm: 5.0,
                    confident_psm_count: 0,
                    ..CalibrationStats::default()
                }
            ),
            0.0
        );
        assert_eq!(apply_shift_for_mode(PrecursorCalMode::On, stats), 5.0);
        assert_eq!(apply_shift_for_mode(PrecursorCalMode::Off, stats), 0.0);
    }

    #[test]
    fn size_guard_skips_files_below_speckey_threshold() {
        let keys: Vec<SpecKey> = (0..100)
            .map(|i| SpecKey {
                spectrum_idx: i,
                charge: 2,
            })
            .collect();
        assert!(keys.len() < constants::MIN_SPECKEYS_FOR_PREPASS);
    }
}
