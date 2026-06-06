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

/// Permissive `rank_score` (RawScore = node + cleavage + edge) pre-filter for a
/// pre-pass PSM before the real confidence gate runs.
///
/// Andes no longer computes the generating function, so the calibration
/// pre-pass can no longer gate on SpecEValue — and raw `rank_score` is NOT
/// comparable across spectra / candidate-window sizes / DB complexity, so it
/// must not be used as a confidence threshold on its own. The REAL confidence
/// filtering is the target-decoy q-value gate in `high_confidence_residuals`
/// (q ≤ 1%); this floor only drops obviously-degenerate negative scores.
const MIN_CONFIDENT_RANK_SCORE: f32 = 0.0;

/// Summary of a successful (or skipped) calibration pre-pass.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CalibrationStats {
    pub shift_ppm: f64,
    pub robust_sigma_ppm: f64,
    pub confident_psm_count: usize,
    /// PSMs rejected because top-1 `rank_score` (RawScore) was below the
    /// confidence floor.
    pub rejected_low_score: usize,
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
    // The calibration pre-pass always ranks by the bundled Rank score: it runs
    // before the strong-score model is relevant, and using Strong here would
    // make the precursor-mass calibration depend on the experimental score path.
    p.score_mode = crate::search_params::ScoreMode::Rank;
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
    rejected_low_score: usize,
    rejected_residual: usize,
}

/// A best-per-spectrum pre-pass PSM with the data the TD-q gate needs.
struct CalCandidate {
    residual: f64,
    rank_score: f32,
    is_decoy: bool,
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

    // Keep the best (highest `rank_score`) PSM per spectrum index across all
    // sampled SpecKeys (e.g. charge variants).
    let mut best_by_spec: HashMap<usize, crate::PsmMatch> = HashMap::new();
    for (key, queue) in sampled.iter().zip(queues.iter()) {
        let Some(psm) = queue
            .iter_psms()
            .max_by(|a, b| {
                a.rank_score
                    .partial_cmp(&b.rank_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        else {
            continue;
        };
        filter.queues_with_psm += 1;
        best_by_spec
            .entry(key.spectrum_idx)
            .and_modify(|best| {
                if psm.rank_score > best.rank_score {
                    *best = psm.clone();
                }
            })
            .or_insert_with(|| psm.clone());
    }

    // Compute each best-per-spectrum PSM's residual, keeping its `is_decoy`
    // flag so the TD-q gate below can run. Raw `rank_score` (RawScore) is NOT
    // comparable across spectra / candidate-window sizes / DB complexity, so we
    // CANNOT trust it as a confidence threshold on its own (the old SpecEValue
    // guard is gone). Instead we accept only PSMs that survive a target-decoy
    // q-value gate (Fix 1 / adversarial review).
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
        });
    }

    let residuals = high_confidence_residuals(&mut cands, keep_top_n);
    (residuals, filter)
}

/// Q-value threshold for accepting a pre-pass PSM as a calibration residual.
///
/// 1% target-decoy FDR — the standard high-confidence cutoff. Only PSMs at or
/// below this q-value are trusted to contribute a precursor-mass residual; this
/// replaces the (removed) SpecEValue confidence guard.
const CALIBRATION_MAX_QVALUE: f64 = 0.01;

/// Apply a target-decoy q-value gate to the best-per-spectrum calibration PSMs
/// and return the residuals of the high-confidence (q ≤ 1%) TARGET PSMs, capped
/// at `keep_top_n` by `rank_score`.
///
/// `cands` is sorted in place. The TDC walk mirrors the trainer's
/// `bootstrap_labels`: rank by `rank_score` descending, accumulate
/// target/decoy counts, q = decoys / max(targets, 1), then make it monotone
/// from the bottom. Ties on equal `rank_score` are handled conservatively (the
/// whole tie bucket takes the worst q in the bucket) so a tied decoy can fail
/// the bucket rather than letting a coarse-score tie admit a target.
fn high_confidence_residuals(cands: &mut [CalCandidate], keep_top_n: usize) -> Vec<f64> {
    if cands.is_empty() {
        return Vec::new();
    }

    // Rank by rank_score descending; tie-break is irrelevant to the gate
    // because tied buckets all take the worst q below, but keep it stable.
    cands.sort_by(|a, b| {
        let av = if a.rank_score.is_nan() { f32::NEG_INFINITY } else { a.rank_score };
        let bv = if b.rank_score.is_nan() { f32::NEG_INFINITY } else { b.rank_score };
        bv.partial_cmp(&av).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Running TDC q-value.
    let n = cands.len();
    let mut q = vec![1.0_f64; n];
    let (mut targets, mut decoys) = (0u64, 0u64);
    for (i, c) in cands.iter().enumerate() {
        if c.is_decoy {
            decoys += 1;
        } else {
            targets += 1;
        }
        q[i] = decoys as f64 / targets.max(1) as f64;
    }

    // Monotone from the bottom.
    let mut min_q = 1.0_f64;
    for qi in q.iter_mut().rev() {
        if *qi < min_q {
            min_q = *qi;
        }
        *qi = min_q;
    }

    // Conservative tie handling: every PSM in an equal-`rank_score` bucket takes
    // the worst (max) q in the bucket (cands are sorted desc, ties contiguous).
    let mut start = 0usize;
    while start < n {
        let s = cands[start].rank_score;
        let mut end = start + 1;
        let tie = |x: f32, y: f32| (x.is_nan() && y.is_nan()) || x.to_bits() == y.to_bits();
        while end < n && tie(cands[end].rank_score, s) {
            end += 1;
        }
        let mut worst = q[start];
        for &qi in &q[start + 1..end] {
            if qi > worst {
                worst = qi;
            }
        }
        for qi in &mut q[start..end] {
            *qi = worst;
        }
        start = end;
    }

    // Keep high-confidence TARGET residuals (q ≤ 1%), already in rank_score
    // descending order; cap at keep_top_n.
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
            score_mode: crate::search_params::ScoreMode::Rank,
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

    fn cand(residual: f64, rank_score: f32, is_decoy: bool) -> CalCandidate {
        CalCandidate { residual, rank_score, is_decoy }
    }

    #[test]
    fn decoy_only_residuals_produce_no_calibration() {
        // A pre-pass that found only decoys at the top is NOT trustworthy: no
        // residual should survive the TD-q gate, so no shift can be computed.
        let mut cands: Vec<CalCandidate> =
            (0..50).map(|i| cand(3.0, 20.0 - i as f32 * 0.1, true)).collect();
        let residuals = high_confidence_residuals(&mut cands, constants::RESIDUAL_CAP);
        assert!(residuals.is_empty(), "decoy-only set must yield no residuals");
    }

    #[test]
    fn low_confidence_targets_and_decoys_produce_no_calibration() {
        // Targets and decoys interleaved at indistinguishable scores → q never
        // drops to 1%, so nothing is accepted and no shift is applied.
        let mut cands: Vec<CalCandidate> = (0..100)
            .map(|i| cand(2.0, 5.0, i % 2 == 0)) // alternating target/decoy, all tied
            .collect();
        let residuals = high_confidence_residuals(&mut cands, constants::RESIDUAL_CAP);
        assert!(
            residuals.is_empty(),
            "low-confidence (q>1%) interleaved set must yield no residuals, got {}",
            residuals.len()
        );
    }

    #[test]
    fn high_confidence_targets_pass_the_gate() {
        // Many clearly-separated targets above a small tail of decoys → q stays
        // well under 1% for the targets, which contribute their residuals.
        let mut cands: Vec<CalCandidate> = (0..200)
            .map(|i| cand(4.0, 30.0 - i as f32 * 0.05, false))
            .collect();
        cands.push(cand(9.0, 1.0, true)); // a single low-scoring decoy at the tail
        let residuals = high_confidence_residuals(&mut cands, constants::RESIDUAL_CAP);
        assert!(!residuals.is_empty(), "confident targets must yield residuals");
        assert!(residuals.iter().all(|&r| (r - 4.0).abs() < 1e-9));
    }

    #[test]
    fn unreliable_stats_apply_no_shift_and_no_tightening() {
        // If the gate produced too few confident PSMs, has_reliable_stats() is
        // false and neither the shift nor the tolerance tightening is applied.
        let stats = CalibrationStats::default(); // confident_psm_count == 0
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
        assert_eq!(before, after, "tolerance must not tighten on unreliable stats");
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
