//! Chimeric two-pass cascade: detect co-isolated precursors in an MS2 scan's MS1
//! isolation window (excluding the selected precursor), then run a targeted
//! second-peptide search at each. Speed-correct chimeric path: scores few
//! candidates at MS1-confirmed masses instead of thousands across the blind
//! window.
#![allow(dead_code)]

use crate::candidate_gen::Candidate;
use crate::chimeric_features::precursor_isotope_match;
use crate::match_engine::compute_spec_e_values_for_spectrum;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::psm::{PsmFeatures, PsmMatch, TopNQueue};
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use model::aa_set::AminoAcidSet;
use model::enzyme::Enzyme;
use model::mass::{nominal_from, H2O, ISOTOPE, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::scoring::{
    predict_by_ions, psm_edge_score, score_psm, RankScorer, ScoredSpectrum,
};
use std::collections::BTreeMap;

/// A co-isolated precursor detected in the MS1 isolation window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CoIsolated {
    pub mono_mz: f64,
    pub charge: u8,
    pub neutral_mass: f64,
}

/// Detect co-isolated precursors in `ms1_peaks` (m/z-sorted) within the isolation
/// window `[win_lo, win_hi]`, EXCLUDING the envelope at `selected_mz` (the peptide
/// Pass 1 already searched). Tries charges in `charge_range`; accepts an envelope
/// whose averagine KL is below `max_kl`. Returns at most `max_n` (highest-intensity
/// monoisotopic peaks first).
#[allow(clippy::too_many_arguments)]
pub(crate) fn detect_coisolated(
    ms1_peaks: &[(f64, f32)],
    win_lo: f64,
    win_hi: f64,
    selected_mz: f64,
    charge_range: std::ops::RangeInclusive<u8>,
    tol_da: f64,
    max_kl: f32,
    max_n: usize,
) -> Vec<CoIsolated> {
    // Candidate monoisotopic peaks = peaks inside the window, sorted by intensity desc.
    let lo_idx = ms1_peaks.partition_point(|&(mz, _)| mz < win_lo);
    let mut cands: Vec<(f64, f32)> = ms1_peaks[lo_idx..]
        .iter()
        .take_while(|&&(mz, _)| mz <= win_hi)
        .copied()
        .collect();
    cands.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut out: Vec<CoIsolated> = Vec::new();
    for &(mz, _inten) in &cands {
        if (mz - selected_mz).abs() <= tol_da {
            continue; // skip the selected precursor (monoisotope)
        }
        // Skip the selected precursor's isotope peaks in BOTH directions: a peak at
        // selected_mz +/- k*ISOTOPE/z is part of the Pass-1 precursor's own envelope,
        // not a distinct species. The lower side matters when the instrument selected
        // M+1/M+2; without this guard Pass 2 re-discovers the primary as a fake
        // secondary (self-inflation). Selected charge is unknown, so reject any isotope
        // spacing in `charge_range` on either side (`d.abs()` makes it symmetric).
        if charge_range.clone().filter(|&z| z != 0).any(|z| {
            let d = (mz - selected_mz).abs();
            (1..6).any(|k| (d - k as f64 * ISOTOPE / z as f64).abs() <= tol_da)
        }) {
            continue;
        }
        // Don't re-report a peak that's an isotope of an already-accepted envelope.
        if out.iter().any(|c| {
            let d = (mz - c.mono_mz).abs();
            (0..6).any(|k| (d - k as f64 * ISOTOPE / c.charge as f64).abs() <= tol_da)
        }) {
            continue;
        }
        // Try charges; accept the lowest-KL charge under max_kl.
        let mut best: Option<(f32, CoIsolated)> = None;
        for z in charge_range.clone() {
            if z == 0 {
                continue;
            }
            let neutral = (mz - PROTON) * z as f64;
            let (kl, _snr) = precursor_isotope_match(ms1_peaks, mz, z, neutral, tol_da, 4);
            if kl <= max_kl && best.as_ref().is_none_or(|(bk, _)| kl < *bk) {
                best = Some((
                    kl,
                    CoIsolated {
                        mono_mz: mz,
                        charge: z,
                        neutral_mass: neutral,
                    },
                ));
            }
        }
        if let Some((_, c)) = best {
            out.push(c);
        }
        if out.len() >= max_n {
            break;
        }
    }
    out
}

/// Raw MS2 peaks claimed by `peptide`'s charge-1 b/y ions, as quantized m/z keys
/// (`round(mz·1000)`). Binary-searches the raw, m/z-sorted `spec.peaks` directly
/// (no `ScoredSpectrum` build). Used by `search_secondary` to strip the primary's
/// peaks before forming the residual; matches `nearest_peak_full`'s max-intensity
/// selection within the tolerance window.
fn primary_matched_peak_keys(
    spec: &Spectrum,
    peptide: &Peptide,
    scorer: &RankScorer,
) -> std::collections::HashSet<i64> {
    let mut keys = std::collections::HashSet::new();
    if peptide.length() < 2 {
        return keys;
    }
    let predicted = predict_by_ions(peptide, 1..=1);
    let tol_is_ppm = scorer.param().data_type.instrument.is_high_resolution();
    let tol = if tol_is_ppm { 20.0_f64 } else { 0.5_f64 };
    for p in &predicted {
        let tol_da = if tol_is_ppm { p.mz * tol / 1e6 } else { tol };
        let lo_mz = p.mz - tol_da;
        let hi_mz = p.mz + tol_da;
        // `spec.peaks` is m/z-sorted; binary-search the window start, scan to `hi_mz`.
        let start = spec.peaks.partition_point(|&(mz, _)| mz < lo_mz);
        let mut best: Option<(f64, f32)> = None; // (peak_mz, intensity)
        for &(mz, intensity) in &spec.peaks[start..] {
            if mz > hi_mz {
                break;
            }
            if best.is_none_or(|(_, best_int)| intensity > best_int) {
                best = Some((mz, intensity));
            }
        }
        if let Some((peak_mz, _)) = best {
            keys.insert((peak_mz * 1000.0).round() as i64);
        }
    }
    keys
}

/// Best secondary PSM for `co` on `spec`, after removing the primary peptide's
/// matched charge-1 b/y peaks (residual spectrum). Scores ONLY candidates within
/// `params.precursor_tolerance` of `co.neutral_mass` (the candidate-count cut that
/// makes the chimeric cascade cheap), then runs ONE targeted GF SpecEValue DP on
/// the residual. Returns `None` if no candidate clears scoring.
///
/// `bucket_index` maps `nominal(peptide.mass() - H2O) -> candidate ids`, identical
/// to `PreparedSearch.bucket_index` (built from `Peptide::nominal_residue_mass`).
/// Returns the winning secondary PSM AND the set of raw-peak keys it claimed (its
/// matched charge-1 b/y peaks). When a scan has multiple co-isolated precursors,
/// the caller threads these keys into `prior_claimed` for the NEXT call so the
/// secondaries COMPETE for residual evidence (a peak explained by one secondary is
/// removed before the next is scored) instead of double-counting shared leftovers.
#[allow(clippy::too_many_arguments)]
pub(crate) fn search_secondary(
    spec: &Spectrum,
    primary: &Peptide,
    prior_claimed: &std::collections::HashSet<i64>,
    co: CoIsolated,
    candidates: &[Candidate],
    bucket_index: &BTreeMap<i32, Vec<usize>>,
    scorer: &RankScorer,
    aa_set: &AminoAcidSet,
    enzyme: Option<Enzyme>,
    params: &SearchParams,
    search_index: &SearchIndex,
    fragment_tolerance_da: f64,
) -> Option<(PsmMatch, std::collections::HashSet<i64>)> {
    let z = co.charge;
    if z == 0 {
        return None;
    }

    // 1. Residual spectrum: drop the peaks already explained on this scan — the
    //    primary's matched charge-1 b/y peaks PLUS any peaks claimed by earlier
    //    secondaries (`prior_claimed`) — so this peptide is scored only against
    //    still-unexplained signal. Overwrite the precursor fields with `co`'s so
    //    the GF mass window (derived from the spectrum's precursor) centers on the
    //    co-isolated mass. `co_spec` feeds both `ScoredSpectrum::new` and the GF.
    let mut claimed = primary_matched_peak_keys(spec, primary, scorer);
    claimed.extend(prior_claimed.iter().copied());
    let mut co_spec = spec.clone();
    co_spec
        .peaks
        .retain(|&(mz, _)| !claimed.contains(&((mz * 1000.0).round() as i64)));
    co_spec.precursor_mz = co.mono_mz;
    co_spec.precursor_charge = Some(co.charge as i32);
    let res_ss = ScoredSpectrum::new(&co_spec, scorer, z);

    // 2. Candidates within precursor tol of the co-isolated neutral mass. Nominal
    //    bucket key matches PreparedSearch: `nominal_from(mass - H2O)`.
    //
    //    Apply the learned precursor calibration shift to the co-isolated neutral
    //    mass FIRST, so the candidate prefilter and reported mass error use the same
    //    calibrated scale as the rest of the search (`matches_precursor` /
    //    `candidate_nominal_bounds`). Without this, a non-zero `--precursor-cal`
    //    shift would exclude the true secondary before GF scoring. The GF below
    //    derives its own window from `co_spec`'s raw precursor m/z and applies the
    //    shift internally, so it stays consistent.
    let co_neutral = adjusted_observed_neutral_mass(co.neutral_mass, params.precursor_mass_shift_ppm);
    let nominal = |m: f64| nominal_from(m - H2O);
    let tol = params.precursor_tolerance.left.as_da(co_neutral).max(0.01);
    let lo = nominal(co_neutral - tol) - 1;
    let hi = nominal(co_neutral + tol) + 1;

    let mut queue = TopNQueue::new(1);
    for (_nm, idxs) in bucket_index.range(lo..=hi) {
        for &ci in idxs {
            let cand = &candidates[ci];
            // Exact-mass gate (the nominal range is integer-coarse).
            if (cand.peptide.mass() - co_neutral).abs() > tol {
                continue;
            }
            let pin = score_psm(&res_ss, &cand.peptide, scorer, z, fragment_tolerance_da);
            let edge = psm_edge_score(&res_ss, &cand.peptide, scorer, z);
            // Mirror the production RawScore scale: `score_psm(...) + cleavage_credit`.
            let cleavage = crate::match_engine::cleavage_credit_for(cand, params.enzyme, aa_set);
            let psm = PsmMatch {
                spectrum_idx: 0,
                candidate_idxs: vec![ci as u32],
                charge_used: z,
                mass_error_ppm: (cand.peptide.mass() - co_neutral) / co_neutral * 1e6,
                score: pin + cleavage as f32,
                rank_score: pin + cleavage as f32 + edge as f32,
                edge_score: edge,
                spec_e_value: 1.0,
                de_novo_score: i32::MIN,
                activation_method: Some(scorer.param().data_type.activation),
                e_value: 1.0,
                features: PsmFeatures::default(),
                isotope_offset: 0,
                // Set to co.mono_mz once the winner is chosen (end of fn).
                precursor_mz_override: None,
            };
            queue.push(psm);
        }
    }

    if queue.is_empty() {
        return None;
    }

    // 3. One targeted GF SpecEValue DP on the residual (fills spec_e_value /
    //    de_novo_score / e_value). The secondary's mass is KNOWN from the MS1
    //    monoisotopic peak, so clamp `isotope_error_range` to 0..=0 to build a
    //    single GF mass bin instead of 5-7 (cuts GF bins + SinkUnreachable retries).
    let mut p2 = params.clone();
    p2.isotope_error_range = 0..=0;
    compute_spec_e_values_for_spectrum(
        &co_spec,
        &p2,
        &mut queue,
        aa_set,
        enzyme,
        scorer,
        &res_ss,
        z,
        fragment_tolerance_da,
        search_index,
        candidates,
    );

    // Pick the winner by SCORE, not heap order: `drain_into_vec` is unordered and
    // `TopNQueue` keeps ties even at capacity 1, so selecting `.next()` would be
    // heap-order dependent (nondeterministic) in a user-visible ranking path.
    // Order by smallest spec_e_value, then largest rank_score, then smallest
    // candidate index as a deterministic final tiebreak.
    let mut best = queue
        .drain_into_vec()
        .into_iter()
        .min_by(|a, b| {
            a.spec_e_value
                .partial_cmp(&b.spec_e_value)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(
                    b.rank_score
                        .partial_cmp(&a.rank_score)
                        .unwrap_or(std::cmp::Ordering::Equal),
                )
                .then_with(|| a.primary_candidate_idx().cmp(&b.primary_candidate_idx()))
        })?;
    // Features on the RESIDUAL (the spectrum the secondary was scored against), so
    // they stay consistent with its RawScore / SpecEValue. The override makes the
    // PIN writer compute ExpMass/dm/absdm from the co-isolated mass.
    let cand_peptide = &candidates[best.primary_candidate_idx() as usize].peptide;
    let mut features = crate::match_engine::compute_psm_features(&res_ss, cand_peptide, scorer, z);
    features.edge_score = best.edge_score;
    best.features = features;
    best.precursor_mz_override = Some(co.mono_mz);
    // Peaks this secondary explained (its charge-1 b/y matches on the FULL spectrum),
    // returned so the caller can remove them before the next co-isolated mass on
    // this scan is searched (sequential competition; see fn doc).
    let winner_claimed = primary_matched_peak_keys(spec, cand_peptide, scorer);
    Some((best, winner_claimed))
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::isotope::averagine_isotope_envelope;

    /// Build a synthetic MS1 peak list (m/z-sorted) containing a 4-peak averagine
    /// envelope for `(mono_mz, charge, neutral_mass)` scaled by `scale`.
    fn envelope(mono_mz: f64, charge: u8, neutral: f64, scale: f32) -> Vec<(f64, f32)> {
        let env = averagine_isotope_envelope(neutral, 4);
        (0..4)
            .map(|k| {
                (
                    mono_mz + k as f64 * ISOTOPE / charge as f64,
                    (env[k] as f32) * scale,
                )
            })
            .collect()
    }

    #[test]
    fn detects_coisolated_excludes_selected() {
        let z = 2u8;
        let selected_mz = 600.0;
        let sel_neutral = (selected_mz - PROTON) * z as f64;
        let co_mz = 600.7; // a second precursor within a ~2 Da window
        let co_neutral = (co_mz - PROTON) * z as f64;
        let mut peaks = envelope(selected_mz, z, sel_neutral, 1000.0);
        peaks.extend(envelope(co_mz, z, co_neutral, 500.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
        assert_eq!(got.len(), 1, "exactly one co-isolated (selected excluded)");
        assert!((got[0].mono_mz - co_mz).abs() < 0.02);
        assert_eq!(got[0].charge, z);
    }

    #[test]
    fn no_coisolation_when_only_selected_present() {
        let z = 2u8;
        let selected_mz = 600.0;
        let peaks = envelope(selected_mz, z, (selected_mz - PROTON) * z as f64, 1000.0);
        let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
        assert!(got.is_empty(), "only the selected precursor -> no co-isolation");
    }

    // ── Task 2: targeted second-peptide residual search ────────────────────

    use model::{AminoAcid, AminoAcidSetBuilder, Protein, ProteinDb};
    use rustc_hash::FxHashMap;
    use scoring_crate::param_model::{IonType, Partition, SpecDataType};
    use scoring_crate::scoring::fragment_ions::predict_by_ions;
    use scoring_crate::Param;
    use crate::PreparedSearch;
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use model::protocol::Protocol;
    use model::Tolerance;

    /// Minimal RankScorer (mirrors `tests/match_engine_smoke.rs::tiny_scorer`):
    /// non-trivial prefix/suffix rank tables so b/y matches earn positive score.
    fn tiny_scorer() -> RankScorer {
        let part = Partition { charge: 2, parent_mass: 500.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let suffix1 = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let noise = IonType::Noise;

        let mut ion_table = FxHashMap::default();
        ion_table.insert(prefix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
        ion_table.insert(suffix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
        ion_table.insert(noise, vec![0.05_f32, 0.05, 0.05, 0.05]);

        let mut rank_dist_table = FxHashMap::default();
        rank_dist_table.insert(part, ion_table);

        let mut frag_off_table = FxHashMap::default();
        frag_off_table.insert(part, vec![]);

        let mut param = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: FxHashMap::default(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
        };
        param.rebuild_cache();
        RankScorer::new(&param)
    }

    fn peptide(residues: &[u8], pre: u8, post: u8) -> Peptide {
        let aas: Vec<AminoAcid> = residues
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(aas, pre, post)
    }

    #[test]
    fn secondary_search_finds_planted_peptide() {
        // Two distinct tryptic peptides from one protein:
        //   WVTFISLLR (positions 2..11) and DAFLGSFLYEYSR (positions 12..25).
        // Plant DAFLGSFLYEYSR's charge-1 b/y ions as the spectrum; pass
        // WVTFISLLR as the (different) primary so residual removal does NOT
        // delete the planted peaks. search_secondary must recover the planted
        // peptide at its own co-isolated mass.
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            }],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = tiny_scorer();
        let frag_tol = 0.02_f64;

        let prepared = PreparedSearch::prepare(&idx, &params, &scorer, frag_tol, "XXX");

        // Planted secondary peptide and its co-isolated precursor.
        let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
        let z = 2u8;
        let co_neutral = planted.mass();
        let co_mz = (co_neutral + z as f64 * PROTON) / z as f64;
        let co = CoIsolated { mono_mz: co_mz, charge: z, neutral_mass: co_neutral };

        // Spectrum peaks = planted peptide's predicted charge-1 b/y ions.
        let peaks: Vec<(f64, f32)> = predict_by_ions(&planted, 1..=1)
            .into_iter()
            .map(|p| (p.mz, 100.0_f32))
            .collect();
        let mut peaks = peaks;
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // Primary precursor is a DIFFERENT peptide (so its matched peaks, if any,
        // don't strip the planted ions). Precursor m/z is the primary's.
        let primary = peptide(b"WVTFISLLR", b'K', b'K');
        let prim_mz = (primary.mass() + z as f64 * PROTON) / z as f64;
        let spec = Spectrum {
            title: "chimeric-secondary".into(),
            precursor_mz: prim_mz,
            precursor_intensity: None,
            precursor_charge: Some(z as i32),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        };

        let got = search_secondary(
            &spec,
            &primary,
            &std::collections::HashSet::new(),
            co,
            &prepared.candidates,
            &prepared.bucket_index,
            &scorer,
            &prepared.aa_set_for_gf,
            Some(params.enzyme),
            &params,
            &idx,
            frag_tol,
        );

        let (psm, winner_claimed) = got.expect("secondary search should return a PSM at the co-isolated mass");
        let found = &prepared.candidates[psm.primary_candidate_idx() as usize].peptide;
        assert_eq!(
            found.residues, planted.residues,
            "secondary PSM should resolve to the planted peptide DAFLGSFLYEYSR"
        );

        // Competition check: searching the SAME co-isolated mass again with the
        // winner's peaks already claimed strips them from the residual, so the
        // peptide can no longer earn its fragment matches — it must NOT come back
        // with the same evidence. (Proves prior_claimed removes shared peaks so two
        // co-isolated precursors can't double-count the same leftover signal.)
        assert!(!winner_claimed.is_empty(), "winner should claim its matched peaks");
        let again = search_secondary(
            &spec,
            &primary,
            &winner_claimed,
            co,
            &prepared.candidates,
            &prepared.bucket_index,
            &scorer,
            &prepared.aa_set_for_gf,
            Some(params.enzyme),
            &params,
            &idx,
            frag_tol,
        );
        let matched_after = again
            .as_ref()
            .map(|(p, _)| p.features.num_matched_main_ions)
            .unwrap_or(0);
        assert!(
            matched_after < psm.features.num_matched_main_ions,
            "with its peaks already claimed, the secondary must match fewer ions \
             (was {}, now {})",
            psm.features.num_matched_main_ions,
            matched_after
        );
        // The planted peptide is in-window for the co-isolated mass, so the GF
        // must compute a real SpecEValue (< 1.0). Exactly 1.0 means the GF mass
        // window did not include the candidate (the `return 1.0` guard fired).
        assert!(
            psm.spec_e_value < 1.0,
            "secondary PSM SpecEValue must be a real probability < 1.0, got {}",
            psm.spec_e_value
        );
    }

    #[test]
    fn secondary_search_applies_precursor_calibration_shift() {
        // Regression for the Pass-2 calibration bug: the co-isolated neutral mass
        // derived from the raw MS1 m/z must be calibration-adjusted before the
        // candidate prefilter, exactly like the main search. Here the MS1
        // observation is biased high by a known ppm shift; with the shift applied,
        // the planted peptide is recovered and its reported error is ~0. Without
        // the fix the prefilter would search the biased mass and either drop the
        // peptide or report a ~`shift_ppm` error.
        let shift_ppm = 30.0_f64; // exaggerated vs a real ~1 ppm shift, to be decisive
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            }],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.precursor_mass_shift_ppm = shift_ppm;
        let scorer = tiny_scorer();
        let frag_tol = 0.02_f64;
        let prepared = PreparedSearch::prepare(&idx, &params, &scorer, frag_tol, "XXX");

        let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
        let z = 2u8;
        // Raw observed mass biased high by `shift_ppm` (so adjusted == true mass).
        let true_mass = planted.mass();
        let raw_neutral = true_mass / (1.0 - shift_ppm * 1e-6);
        let co_mz = (raw_neutral + z as f64 * PROTON) / z as f64;
        let co = CoIsolated { mono_mz: co_mz, charge: z, neutral_mass: raw_neutral };

        let mut peaks: Vec<(f64, f32)> = predict_by_ions(&planted, 1..=1)
            .into_iter()
            .map(|p| (p.mz, 100.0_f32))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let primary = peptide(b"WVTFISLLR", b'K', b'K');
        let prim_mz = (primary.mass() + z as f64 * PROTON) / z as f64;
        let spec = Spectrum {
            title: "chimeric-cal".into(),
            precursor_mz: prim_mz,
            precursor_intensity: None,
            precursor_charge: Some(z as i32),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        };

        let got = search_secondary(
            &spec, &primary, &std::collections::HashSet::new(), co, &prepared.candidates,
            &prepared.bucket_index, &scorer, &prepared.aa_set_for_gf, Some(params.enzyme),
            &params, &idx, frag_tol,
        );
        let (psm, _claimed) = got.expect("secondary must be found at the calibration-adjusted mass");
        let found = &prepared.candidates[psm.primary_candidate_idx() as usize].peptide;
        assert_eq!(found.residues, planted.residues, "should resolve to the planted peptide");
        // The candidate matches the ADJUSTED mass, so the reported error is ~0,
        // not ~-30 ppm (which is what the un-calibrated raw mass would report).
        assert!(
            psm.mass_error_ppm.abs() < 5.0,
            "mass error must be reported on the calibrated scale (~0), got {} ppm",
            psm.mass_error_ppm
        );
    }
}
