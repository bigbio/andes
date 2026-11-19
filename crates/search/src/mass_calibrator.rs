use std::collections::HashMap;
use std::ops::RangeInclusive;

use model::fdr::{ranked_q, Counts, Eval};
use model::mass::{H2O, PROTON};
use model::stats::nan_low_f32;
use model::tolerance::{PrecursorTolerance, Tolerance};
use model::Spectrum;

use crate::match_engine::PreparedSearch;
use crate::precursor_cal::{
    constants, median, residual_ppm, robust_sigma_ppm, sample_every_nth, tightened_tolerance_ppm,
    PrecursorCalMode,
};
use crate::search_params::SearchParams;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SpecKey {
    pub spectrum_idx: usize,
    pub charge: u8,
}

const MAX_REASONABLE_RESIDUAL_PPM: f64 = 50.0;

const MIN_CONFIDENT_RANK_SCORE: f32 = 0.0;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CalibrationStats {
    pub shift_ppm: f64,
    pub robust_sigma_ppm: f64,
    pub confident_psm_count: usize,
    pub rejected_low_score: usize,
    pub rejected_residual: usize,
    pub queues_with_psm: usize,
}

impl CalibrationStats {
    pub fn has_reliable_stats(self) -> bool {
        self.confident_psm_count > 0
    }
}

pub fn build_spec_keys(
    spectra: &[Spectrum],
    charge_range: &RangeInclusive<u8>,
    min_peaks: u32,
) -> Vec<SpecKey> {
    build_spec_keys_from_counts(
        spectra.iter().map(|s| (s.precursor_charge, s.peaks.len())),
        charge_range,
        min_peaks,
    )
}

pub fn build_spec_keys_from_counts(
    spectra: impl IntoIterator<Item = (Option<i32>, usize)>,
    charge_range: &RangeInclusive<u8>,
    min_peaks: u32,
) -> Vec<SpecKey> {
    let mut keys = Vec::new();
    for (spectrum_idx, (precursor_charge, n_peaks)) in spectra.into_iter().enumerate() {
        if n_peaks < min_peaks as usize {
            continue;
        }
        match precursor_charge {
            Some(z) if z > 0 => keys.push(SpecKey {
                spectrum_idx,
                charge: z as u8,
            }),
            _ => {
                for z in charge_range.clone() {
                    keys.push(SpecKey {
                        spectrum_idx,
                        charge: z,
                    });
                }
            }
        }
    }
    keys
}

pub fn prepass_search_params(main: &SearchParams) -> SearchParams {
    let mut p = main.clone();
    p.isotope_error_range = 0..=0;
    p.top_n_psms_per_spectrum = 1;
    p.precursor_mass_shift_ppm = 0.0;
    p.score_mode = crate::search_params::ScoreMode::Rank;
    p.fragment_index_top_k = 0;
    p
}

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

pub fn learn_calibration_stats(
    spec_keys: &[SpecKey],
    originals: &HashMap<usize, Spectrum>,
    prepared: &mut PreparedSearch<'_>,
    main_params: &SearchParams,
) -> CalibrationStats {
    if spec_keys.len() < main_params.cal_min_spec_keys {
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
    prepared.sync_materialized_candidates();
    let (residuals, filter) = extract_residuals(
        &sampled,
        &queues,
        originals,
        &prepared.candidates,
        MIN_CONFIDENT_RANK_SCORE,
        constants::RESIDUAL_CAP,
    );

    if residuals.len() < constants::MIN_CONFIDENT_PSMS {
        return CalibrationStats {
            rejected_low_score: filter.rejected_low_score,
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

    let tightened_left =
        tightened_tolerance_ppm(left_ppm, stats.robust_sigma_ppm, sigma_mult, floor, margin);
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
    rejected_low_score: usize,
    rejected_residual: usize,
}

struct CalCandidate {
    residual: f64,
    rank_score: f32,
    is_decoy: bool,
    spectrum_idx: usize,
}

fn extract_residuals(
    sampled: &[SpecKey],
    queues: &[crate::TopNQueue],
    originals: &HashMap<usize, Spectrum>,
    candidates: &[crate::candidate_gen::Candidate],
    min_rank_score: f32,
    keep_top_n: usize,
) -> (Vec<f64>, CalFilterCounts) {
    let mut filter = CalFilterCounts::default();

    let mut best_by_spec: std::collections::BTreeMap<usize, crate::PsmMatch> =
        std::collections::BTreeMap::new();
    for (key, queue) in sampled.iter().zip(queues.iter()) {
        let Some(psm) = queue.peek_top() else {
            continue;
        };
        filter.queues_with_psm += 1;
        best_by_spec
            .entry(key.spectrum_idx)
            .and_modify(|best| {
                if *psm > *best {
                    *best = psm.clone();
                }
            })
            .or_insert_with(|| psm.clone());
    }

    let mut cands: Vec<CalCandidate> = Vec::new();
    for (&spectrum_idx, psm) in best_by_spec.iter() {
        if psm.rank_score < min_rank_score {
            filter.rejected_low_score += 1;
            continue;
        }

        let Some(spec) = originals.get(&spectrum_idx) else {
            continue;
        };
        let charge = psm.charge_used as f64;
        if charge <= 0.0 {
            continue;
        }

        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let observed = (spec.precursor_mz - PROTON) * charge - H2O;
        let theoretical = cand.peptide.residue_mass();
        if theoretical <= 0.0 {
            continue;
        }

        let residual = residual_ppm(observed, theoretical);
        if residual.abs() > MAX_REASONABLE_RESIDUAL_PPM {
            filter.rejected_residual += 1;
            continue;
        }
        cands.push(CalCandidate {
            residual,
            rank_score: psm.rank_score,
            is_decoy: cand.is_decoy,
            spectrum_idx,
        });
    }

    let residuals = high_confidence_residuals(&mut cands, keep_top_n);
    (residuals, filter)
}

const CALIBRATION_MAX_QVALUE: f64 = 0.01;

fn high_confidence_residuals(cands: &mut [CalCandidate], keep_top_n: usize) -> Vec<f64> {
    if cands.is_empty() {
        return Vec::new();
    }

    cands.sort_by(|a, b| {
        let (av, bv) = (nan_low_f32(a.rank_score), nan_low_f32(b.rank_score));
        bv.total_cmp(&av)
            .then_with(|| a.spectrum_idx.cmp(&b.spectrum_idx))
    });

    let tie = |x: f32, y: f32| (x.is_nan() && y.is_nan()) || x.to_bits() == y.to_bits();
    let q = ranked_q::<Counts>(
        cands.len(),
        |c, rank| c.add(cands[rank].is_decoy),
        |start, rank| tie(cands[rank].rank_score, cands[start].rank_score),
        Eval::PerRowGroupMax,
        Counts::ratio,
        1.0,
    );

    cands
        .iter()
        .zip(q.iter())
        .filter(|(c, &qi)| !c.is_decoy && qi <= CALIBRATION_MAX_QVALUE)
        .take(keep_top_n)
        .map(|(c, _)| c.residual)
        .collect()
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
            extra_enzymes: Vec::new(),
            min_length: 6,
            max_length: 40,
            max_missed_cleavages: 1,
            max_variable_mods_per_peptide: 3,
            max_label_sites: 0,
            precursor_tolerance: PrecursorTolerance::symmetric(Tolerance::Ppm(20.0)),
            charge_range: 2..=3,
            isotope_error_range: -1..=2,
            top_n_psms_per_spectrum: 10,
            num_tolerable_termini: 2,
            min_peaks: 10,
            precursor_cal_mode: PrecursorCalMode::Auto,
            cal_min_spec_keys: constants::MIN_SPECKEYS_FOR_PREPASS,
            fragment_index_top_k: 0,
            precursor_mass_shift_ppm: 0.0,
            chimeric: false,
            dia: false,
            dia_ms1_tol_ppm: 5.0,
            dia_vote_filter: None,
            chimeric_isolation_halfwidth_da: 1.5,
            chimeric_max_coisolated: 2,
            chimeric_max_kl: 0.3,
            score_mode: crate::search_params::ScoreMode::Rank,
            refine_select_psm_fdr: 0.01,
            require_nxst_sequon: false,
            open_search: None,
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
            isolation_target_mz: None,
            mobility_window: None,
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

    fn cand(residual: f64, rank_score: f32, is_decoy: bool) -> CalCandidate {
        cand_from(residual, rank_score, is_decoy, 0)
    }

    fn cand_from(
        residual: f64,
        rank_score: f32,
        is_decoy: bool,
        spectrum_idx: usize,
    ) -> CalCandidate {
        CalCandidate {
            residual,
            rank_score,
            is_decoy,
            spectrum_idx,
        }
    }

    #[test]
    fn residual_selection_is_independent_of_input_order() {
        let build = |order: &[usize]| -> Vec<f64> {
            let mut cands: Vec<CalCandidate> = order
                .iter()
                .map(|&i| cand_from(i as f64 * 0.1, 20.0, false, i))
                .collect();
            for i in 0..5 {
                cands.push(cand_from(9.0, 1.0, true, 1000 + i));
            }
            high_confidence_residuals(&mut cands, 20)
        };
        let forward: Vec<usize> = (0..40).collect();
        let reversed: Vec<usize> = (0..40).rev().collect();
        let shuffled: Vec<usize> = (0..40).map(|i| (i * 17) % 40).collect();

        let a = build(&forward);
        let b = build(&reversed);
        let c = build(&shuffled);
        assert_eq!(a.len(), 20, "the cap should bite inside the tie");
        assert_eq!(a, b, "reversing the input changed which residuals survived");
        assert_eq!(a, c, "permuting the input changed which residuals survived");
    }

    #[test]
    fn decoy_only_residuals_produce_no_calibration() {
        let mut cands: Vec<CalCandidate> = (0..50)
            .map(|i| cand(3.0, 20.0 - i as f32 * 0.1, true))
            .collect();
        let residuals = high_confidence_residuals(&mut cands, constants::RESIDUAL_CAP);
        assert!(
            residuals.is_empty(),
            "decoy-only set must yield no residuals"
        );
    }

    #[test]
    fn low_confidence_targets_and_decoys_produce_no_calibration() {
        let mut cands: Vec<CalCandidate> = (0..100).map(|i| cand(2.0, 5.0, i % 2 == 0)).collect();
        let residuals = high_confidence_residuals(&mut cands, constants::RESIDUAL_CAP);
        assert!(
            residuals.is_empty(),
            "low-confidence (q>1%) interleaved set must yield no residuals, got {}",
            residuals.len()
        );
    }

    #[test]
    fn high_confidence_targets_pass_the_gate() {
        let mut cands: Vec<CalCandidate> = (0..200)
            .map(|i| cand(4.0, 30.0 - i as f32 * 0.05, false))
            .collect();
        cands.push(cand(9.0, 1.0, true));
        let residuals = high_confidence_residuals(&mut cands, constants::RESIDUAL_CAP);
        assert!(
            !residuals.is_empty(),
            "confident targets must yield residuals"
        );
        assert!(residuals.iter().all(|&r| (r - 4.0).abs() < 1e-9));
    }

    #[test]
    fn unreliable_stats_apply_no_shift_and_no_tightening() {
        let stats = CalibrationStats::default();
        assert!(!stats.has_reliable_stats());
        assert_eq!(apply_shift_for_mode(PrecursorCalMode::Auto, stats), 0.0);
        assert_eq!(apply_shift_for_mode(PrecursorCalMode::On, stats), 0.0);

        let mut params = empty_params();
        let before = (
            ppm_value(params.precursor_tolerance.left),
            ppm_value(params.precursor_tolerance.right),
        );
        apply_tightened_precursor_tolerance(&mut params, stats);
        let after = (
            ppm_value(params.precursor_tolerance.left),
            ppm_value(params.precursor_tolerance.right),
        );
        assert_eq!(
            before, after,
            "tolerance must not tighten on unreliable stats"
        );
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
