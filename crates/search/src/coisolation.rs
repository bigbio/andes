//! Chimeric two-pass cascade: detect co-isolated precursors in an MS2 scan's MS1
//! isolation window (excluding the selected precursor), then run a targeted
//! second-peptide search at each. This is the speed-correct chimeric path: it
//! scores few candidates at MS1-confirmed masses instead of thousands across the
//! blind window (see docs/parity-analysis/notes/2026-05-30-chimeric-cost-profile.md).
//!
//! Task 1 (this commit) ships only the detector. The targeted second-peptide
//! search (`search_secondary`) and the binary-level driver land in Tasks 2/3,
//! which consume `CoIsolated` / `detect_coisolated`. Until then they are
//! unreferenced outside tests, so allow dead_code at the module level.
#![allow(dead_code)]

use crate::candidate_gen::Candidate;
use crate::chimeric_features::precursor_isotope_match;
use crate::match_engine::{compute_spec_e_values_for_spectrum, matched_peak_keys};
use crate::psm::{PsmFeatures, PsmMatch, TopNQueue};
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use model::aa_set::AminoAcidSet;
use model::enzyme::Enzyme;
use model::mass::{nominal_from, H2O, ISOTOPE, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::scoring::{psm_edge_score, score_psm, RankScorer, ScoredSpectrum};
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
        // Skip the selected precursor's HIGHER isotope peaks too: a peak at
        // selected_mz + k*ISOTOPE/z (k >= 1) is part of the Pass-1 envelope, not
        // a distinct co-isolated species. The selected charge is unknown here, so
        // reject if the peak lines up with any isotope spacing in `charge_range`.
        if mz > selected_mz
            && charge_range.clone().filter(|&z| z != 0).any(|z| {
                let d = mz - selected_mz;
                (1..6).any(|k| (d - k as f64 * ISOTOPE / z as f64).abs() <= tol_da)
            })
        {
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

/// Best secondary PSM for `co` on `spec`, after removing the primary peptide's
/// matched charge-1 b/y peaks (residual spectrum). Scores ONLY candidates within
/// `params.precursor_tolerance` of `co.neutral_mass` (the candidate-count cut that
/// makes the chimeric cascade cheap), then runs ONE targeted GF SpecEValue DP on
/// the residual. Returns `None` if no candidate clears scoring.
///
/// `bucket_index` maps `nominal(peptide.mass() - H2O) -> candidate ids`, identical
/// to `PreparedSearch.bucket_index` (built from `Peptide::nominal_residue_mass`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn search_secondary(
    spec: &Spectrum,
    primary: &Peptide,
    co: CoIsolated,
    candidates: &[Candidate],
    bucket_index: &BTreeMap<i32, Vec<usize>>,
    scorer: &RankScorer,
    aa_set: &AminoAcidSet,
    enzyme: Option<Enzyme>,
    params: &SearchParams,
    search_index: &SearchIndex,
    fragment_tolerance_da: f64,
) -> Option<PsmMatch> {
    let z = co.charge;
    if z == 0 {
        return None;
    }

    // 1. Residual spectrum: drop the primary's matched charge-1 b/y peaks so the
    //    second peptide is scored against signal the primary did not explain.
    //    Build ONE `co_spec` = residual peaks + the CO-ISOLATED precursor fields,
    //    used for BOTH `ScoredSpectrum::new` and the GF below. The residual is a
    //    clone of `spec`, so its `precursor_mz` is the PRIMARY's; overwrite the
    //    precursor fields with `co`'s so the GF mass window (derived from the
    //    spectrum's precursor in `compute_spec_e_values_for_spectrum`) centers on
    //    the co-isolated mass rather than the primary mass (Bug-1 fix).
    let full_ss = ScoredSpectrum::new(spec, scorer, z);
    let claimed = matched_peak_keys(&full_ss, primary, scorer);
    let mut co_spec = spec.clone();
    co_spec
        .peaks
        .retain(|&(mz, _)| !claimed.contains(&((mz * 1000.0).round() as i64)));
    co_spec.precursor_mz = co.mono_mz;
    co_spec.precursor_charge = Some(co.charge as i32);
    let res_ss = ScoredSpectrum::new(&co_spec, scorer, z);

    // 2. Candidates within precursor tol of the co-isolated neutral mass. The
    //    nominal bucket key matches PreparedSearch's `nominal_residue_mass`
    //    convention: `nominal_from(mass - H2O)`.
    let nominal = |m: f64| nominal_from(m - H2O);
    let tol = params.precursor_tolerance.left.as_da(co.neutral_mass).max(0.01);
    let lo = nominal(co.neutral_mass - tol) - 1;
    let hi = nominal(co.neutral_mass + tol) + 1;

    let mut queue = TopNQueue::new(1);
    for (_nm, idxs) in bucket_index.range(lo..=hi) {
        for &ci in idxs {
            let cand = &candidates[ci];
            // Exact-mass gate (the nominal range is integer-coarse).
            if (cand.peptide.mass() - co.neutral_mass).abs() > tol {
                continue;
            }
            let pin = score_psm(&res_ss, &cand.peptide, scorer, z, fragment_tolerance_da);
            let edge = psm_edge_score(&res_ss, &cand.peptide, scorer, z);
            // Cleavage credit: production candidate loop emits RawScore as
            // `score_psm(...) + cleavage_credit`; mirror it so secondary PSMs
            // share the primaries' RawScore scale (Bug-2 fix).
            let cleavage = crate::match_engine::cleavage_credit_for(cand, params.enzyme, aa_set);
            let psm = PsmMatch {
                spectrum_idx: 0,
                candidate_idxs: vec![ci as u32],
                charge_used: z,
                mass_error_ppm: (cand.peptide.mass() - co.neutral_mass) / co.neutral_mass * 1e6,
                score: pin + cleavage as f32,
                rank_score: pin + cleavage as f32 + edge as f32,
                edge_score: edge,
                spec_e_value: 1.0,
                de_novo_score: i32::MIN,
                activation_method: Some(scorer.param().data_type.activation),
                e_value: 1.0,
                features: PsmFeatures::default(),
                isotope_offset: 0,
            };
            queue.push(psm);
        }
    }

    if queue.is_empty() {
        return None;
    }

    // 3. One targeted GF SpecEValue DP on the residual spectrum (fills
    //    spec_e_value / de_novo_score / e_value on the retained PSM). Pass
    //    `co_spec` (residual peaks + co precursor fields) so the GF mass window
    //    centers on the co-isolated mass; `res_ss` supplies the node scores.
    //
    //    Perf: the secondary's mass is KNOWN (co.neutral_mass, detected from the
    //    MS1 monoisotopic peak), so we only need the single GF mass bin at
    //    isotope offset 0. Clamp `isotope_error_range` to 0..=0 (vs the default
    //    -1..=2) so `compute_spec_e_values_for_spectrum` builds ~1 mass bin
    //    instead of 5-7, cutting GF bins and the associated SinkUnreachable
    //    retries. The candidate was enumerated within `precursor_tolerance` of
    //    `co.neutral_mass`, so with isotope 0 + the same precursor tol it stays
    //    in-window (see `secondary_search_finds_planted_peptide`).
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
    queue.drain_into_vec().into_iter().next()
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

        let psm = got.expect("secondary search should return a PSM at the co-isolated mass");
        let found = &prepared.candidates[psm.primary_candidate_idx() as usize].peptide;
        assert_eq!(
            found.residues, planted.residues,
            "secondary PSM should resolve to the planted peptide DAFLGSFLYEYSR"
        );
        // Bug-1 lock-in: the planted peptide is in-window for the co-isolated
        // mass, so the GF must compute a real SpecEValue (< 1.0). A value of
        // exactly 1.0 means the GF mass window did not include the candidate
        // mass (the `return 1.0` guard fired).
        assert!(
            psm.spec_e_value < 1.0,
            "secondary PSM SpecEValue must be a real probability < 1.0, got {}",
            psm.spec_e_value
        );
    }
}
