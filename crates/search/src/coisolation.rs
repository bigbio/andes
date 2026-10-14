//! Chimeric two-pass cascade: detect co-isolated precursors in an MS2 scan's MS1
//! isolation window (excluding the selected precursor), then run a targeted
//! second-peptide search at each. Speed-correct chimeric path: scores few
//! candidates at MS1-confirmed masses instead of thousands across the blind
//! window.
use crate::candidate_gen::Candidate;
use crate::chimeric_features::precursor_isotope_match;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::psm::{PsmFeatures, PsmMatch};
use crate::search_params::SearchParams;
use model::aa_set::AminoAcidSet;
use model::mass::{nominal_from, H2O, ISOTOPE, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::scoring::{
    fuse_strong_score, predict_by_ions, psm_edge_score, score_psm, score_psm_float,
    strong_score_calibrated, OnlineStats, RankScorer, ScoredSpectrum, StrongScoreInputs,
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

    let is_selected_isotope = |mz: f64| {
        charge_range.clone().filter(|&z| z != 0).any(|z| {
            let d = (mz - selected_mz).abs();
            (1..6).any(|k| (d - k as f64 * ISOTOPE / z as f64).abs() <= tol_da)
        })
    };
    let mut out: Vec<CoIsolated> = Vec::new();
    for &(mz, _inten) in &cands {
        // Cap check BEFORE accepting any candidate, so `max_n == 0` is a hard
        // disable (zero secondaries), not "one then stop".
        if out.len() >= max_n {
            break;
        }
        if (mz - selected_mz).abs() <= tol_da {
            continue; // skip the selected precursor (monoisotope)
        }
        // Skip the selected precursor's isotope peaks in BOTH directions: a peak at
        // selected_mz +/- k*ISOTOPE/z is part of the Pass-1 precursor's own envelope,
        // not a distinct species. The lower side matters when the instrument selected
        // M+1/M+2; without this guard Pass 2 re-discovers the primary as a fake
        // secondary (self-inflation). Selected charge is unknown, so reject any isotope
        // spacing in `charge_range` on either side (`d.abs()` makes it symmetric).
        if is_selected_isotope(mz) {
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
        if let Some((kl, c)) = best {
            out.push(walk_to_monoisotope(ms1_peaks, c, kl, tol_da, |m| {
                (m - selected_mz).abs() <= tol_da || is_selected_isotope(m)
            }));
        }
    }
    out
}

/// Most isotope steps [`walk_to_monoisotope`] takes below an accepted peak.
const MAX_MONO_WALK: usize = 3;

/// The monoisotopic form of an accepted envelope. Peaks are tried most intense
/// first, and above ~1.8 kDa the M+1 peak outweighs M, so the envelope is first
/// accepted at M+1 (its shifted shape still fits averagine) and the true M is then
/// skipped as its isotope: the secondary is searched one isotope spacing too
/// heavy. Step down by `ISOTOPE / z` while a peak sits there (inside the window
/// or not, since an M+1 in the window means M was co-isolated too), the envelope
/// from it fits averagine strictly better, and it is not the selected
/// precursor's own envelope (`is_selected`).
fn walk_to_monoisotope(
    ms1_peaks: &[(f64, f32)],
    mut c: CoIsolated,
    mut kl: f32,
    tol_da: f64,
    is_selected: impl Fn(f64) -> bool,
) -> CoIsolated {
    for _ in 0..MAX_MONO_WALK {
        let lower = c.mono_mz - ISOTOPE / c.charge as f64;
        let start = ms1_peaks.partition_point(|&(mz, _)| mz < lower - tol_da);
        let has_peak = ms1_peaks[start..]
            .first()
            .is_some_and(|&(mz, _)| mz <= lower + tol_da);
        if !has_peak || is_selected(lower) {
            break;
        }
        let neutral = (lower - PROTON) * c.charge as f64;
        let (lower_kl, _) = precursor_isotope_match(ms1_peaks, lower, c.charge, neutral, tol_da, 4);
        if lower_kl >= kl {
            break;
        }
        c = CoIsolated {
            mono_mz: lower,
            charge: c.charge,
            neutral_mass: neutral,
        };
        kl = lower_kl;
    }
    c
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
    for p in &predicted {
        // Use the SAME tolerance the scoring phase matched the primary's peaks
        // with (`param().mme`), NOT the PIN feature tolerance — otherwise pass-2
        // strips a different peak set than pass-1 actually scored, leaving stray
        // primary peaks in (or pulling extra peaks out of) the residual.
        let tol_da = scorer.param().mme.as_da(p.mz);
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

/// Largest `rank_score` first, then smallest candidate index: the deterministic
/// order secondaries are ranked in.
fn rank_desc_then_idx(a: &PsmMatch, b: &PsmMatch) -> std::cmp::Ordering {
    b.rank_score
        .partial_cmp(&a.rank_score)
        .unwrap_or(std::cmp::Ordering::Equal)
        .then_with(|| a.primary_candidate_idx().cmp(&b.primary_candidate_idx()))
}

/// Calibration-adjusted neutral mass of `co` and the precursor tolerance (Da)
/// around it that a secondary candidate must fall within.
pub(crate) fn secondary_mass_window(co: CoIsolated, params: &SearchParams) -> (f64, f64) {
    let co_neutral =
        adjusted_observed_neutral_mass(co.neutral_mass, params.precursor_mass_shift_ppm);
    let tol = params.precursor_tolerance.left.as_da(co_neutral).max(0.01);
    (co_neutral, tol)
}

/// Inclusive `nominal(peptide.mass() - H2O)` bucket range [`search_secondary`]
/// scans for `co`: the exact-mass window of [`secondary_mass_window`] widened by
/// one nominal unit on each side.
pub(crate) fn secondary_nominal_window(co: CoIsolated, params: &SearchParams) -> (i32, i32) {
    let (co_neutral, tol) = secondary_mass_window(co, params);
    let nominal = |m: f64| nominal_from(m - H2O);
    (nominal(co_neutral - tol) - 1, nominal(co_neutral + tol) + 1)
}

/// Best secondary PSM for `co` on `spec`, after removing the primary peptide's
/// matched charge-1 b/y peaks (residual spectrum). Scores ONLY candidates within
/// `params.precursor_tolerance` of `co.neutral_mass` (the candidate-count cut that
/// makes the chimeric cascade cheap), then runs ONE targeted GF SpecEValue DP on
/// the residual. Returns `None` if no candidate clears scoring.
///
/// `bucket_index` maps `nominal(peptide.mass() - H2O) -> candidate ids` (the
/// key of `Peptide::nominal_residue_mass`) over `candidates`, a Pass-2 pool.
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
    params: &SearchParams,
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
    //    the residual `ScoredSpectrum` is built at the co-isolated precursor mass.
    let mut claimed = primary_matched_peak_keys(spec, primary, scorer);
    claimed.extend(prior_claimed.iter().copied());
    // Build the residual peak list in a single pass (drop claimed peaks) instead
    // of cloning the whole Spectrum and `retain()`-ing it — avoids copying the
    // full peak vector per co-isolated secondary on dense Astral scans.
    let residual_peaks: Vec<(f64, f32)> = spec
        .peaks
        .iter()
        .copied()
        .filter(|&(mz, _)| !claimed.contains(&((mz * 1000.0).round() as i64)))
        .collect();
    let co_spec = Spectrum {
        title: spec.title.clone(),
        precursor_mz: co.mono_mz,
        precursor_charge: Some(co.charge as i32),
        precursor_intensity: spec.precursor_intensity,
        rt_seconds: spec.rt_seconds,
        scan: spec.scan,
        peaks: residual_peaks,
        activation_method: spec.activation_method,
        isolation_lower_offset: spec.isolation_lower_offset,
        isolation_upper_offset: spec.isolation_upper_offset,
        isolation_target_mz: spec.isolation_target_mz,
    };
    let res_ss = ScoredSpectrum::new(&co_spec, scorer, z);

    // 2. Candidates within precursor tol of the co-isolated neutral mass. Nominal
    //    bucket key: `nominal_from(mass - H2O)`.
    //
    //    Apply the learned precursor calibration shift to the co-isolated neutral
    //    mass FIRST, so the candidate prefilter and reported mass error use the same
    //    calibrated scale as the rest of the search (`matches_precursor` /
    //    `candidate_nominal_bounds`). Without this, a non-zero `--precursor-cal`
    //    shift would exclude the true secondary before GF scoring. The GF below
    //    derives its own window from `co_spec`'s raw precursor m/z and applies the
    //    shift internally, so it stays consistent.
    let (co_neutral, tol) = secondary_mass_window(co, params);
    let (lo, hi) = secondary_nominal_window(co, params);

    // Tailor + DeltaRawScore accumulation over this secondary's candidate score
    // distribution on the RESIDUAL spectrum, as `run_chunk_inner` does per spectrum.
    let mut tailor_hist: std::collections::HashMap<i32, u32> = std::collections::HashMap::new();
    let mut tailor_total: u32 = 0;
    let mut best_raw: i32 = i32::MIN;
    let mut best_pep_key: u64 = 0;
    let mut have_best = false;
    let mut second_raw: i32 = i32::MIN;
    let mut strong_null = OnlineStats::default();
    // Best PSM per distinct (peptidoform, decoy). The pool holds one candidate per
    // protein occurrence, so a shared target peptide appears several times with
    // the same score; ranked as-is, its listwise gap would be 0 against itself, a
    // value decoys (rarely shared) almost never take.
    let mut distinct: rustc_hash::FxHashMap<(u64, bool), PsmMatch> =
        rustc_hash::FxHashMap::default();
    for (_nm, idxs) in bucket_index.range(lo..=hi) {
        for &ci in idxs {
            let cand = &candidates[ci];
            // Exact-mass gate (the nominal range is integer-coarse).
            if (cand.peptide.mass() - co_neutral).abs() > tol {
                continue;
            }
            // Don't re-emit the primary peptide as a secondary on its own scan: a
            // co-isolated mass can fall within tol of the primary's mass, which
            // would yield a duplicate PSM on the same spectrum and inflate FDR.
            if cand.peptide == *primary {
                continue;
            }
            let pin = score_psm(&res_ss, &cand.peptide, scorer, z, fragment_tolerance_da);
            let edge = psm_edge_score(&res_ss, &cand.peptide, scorer, z);
            // Production RawScore scale: `score_psm(...) + cleavage_credit`.
            let cleavage = crate::match_engine::cleavage_credit_for(cand, params.enzyme, aa_set);
            let raw_score = pin + cleavage as f32;

            strong_null.push(raw_score);
            let s = raw_score.round() as i32;
            *tailor_hist.entry(s).or_insert(0) += 1;
            tailor_total += 1;
            // Best / second-best DISTINCT peptidoform, keyed on peptidoform identity:
            // the exact-mass gate admits only one nominal mass, so a nominal key
            // would never see a runner-up.
            let mkey = crate::match_engine::peptidoform_key(&cand.peptide);
            if have_best && mkey == best_pep_key {
                if s > best_raw {
                    best_raw = s;
                }
            } else if s > best_raw {
                second_raw = second_raw.max(best_raw);
                best_raw = s;
                best_pep_key = mkey;
                have_best = true;
            } else {
                second_raw = second_raw.max(s);
            }

            let psm = PsmMatch {
                spectrum_idx: 0,
                candidate_idxs: vec![ci as u32],
                charge_used: z,
                mass_error_ppm: (cand.peptide.mass() - co_neutral) / co_neutral * 1e6,
                score: raw_score,
                rank_score: raw_score + edge as f32,
                edge_score: edge,
                activation_method: Some(scorer.param().data_type.activation),
                features: PsmFeatures::default(),
                isotope_offset: 0,
                // Set to co.mono_mz once the winner is chosen (end of fn).
                precursor_mz_override: None,
            };
            distinct
                .entry((mkey, cand.is_decoy))
                .and_modify(|kept| {
                    if rank_desc_then_idx(&psm, kept) == std::cmp::Ordering::Less {
                        *kept = psm.clone();
                    }
                })
                .or_insert(psm);
        }
    }

    if distinct.is_empty() {
        return None;
    }

    // 3. The retained pool: the best distinct peptidoforms by `rank_score`, cut to
    // the same retention cap the primary queue uses, ordered deterministically
    // (largest rank_score, then smallest candidate index).
    let mut retained: Vec<PsmMatch> = distinct.into_values().collect();
    retained.sort_by(rank_desc_then_idx);
    retained.truncate(params.hot_loop_retention_cap().max(1) as usize);

    // Listwise terms over the retained pool's RawScores, as for primaries.
    let mut retained_scores: Vec<f32> = retained.iter().map(|p| p.score).collect();
    retained_scores.sort_by(|a, b| b.total_cmp(a)); // total order; release-safe
    let listwise_score_gap = scoring_crate::scoring::listwise_score_gap(&retained_scores);
    let candidate_rank_entropy = scoring_crate::scoring::candidate_rank_entropy(&retained_scores);
    let tailor_denom = crate::psm::tailor_denominator(&tailor_hist, tailor_total) as f32;
    let delta_trackers = crate::match_engine::DeltaTrackers {
        best_raw,
        best_pep_key,
        have_best,
        second_raw,
    };

    // Features on the RESIDUAL (the spectrum the secondary was scored against),
    // filled exactly as `run_chunk_inner` fills a primary's.
    let make_features = |psm: &PsmMatch, full: bool| -> PsmFeatures {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let mut f = crate::match_engine::compute_psm_features_with_edge(
            &res_ss,
            &cand.peptide,
            scorer,
            z,
            Some(psm.edge_score),
            full,
        );
        f.delta_raw_score = delta_trackers
            .delta_for(psm.score, crate::match_engine::peptidoform_key(&cand.peptide));
        f.listwise_score_gap = listwise_score_gap;
        f.candidate_rank_entropy = candidate_rank_entropy;
        f.tailor_score = psm.score / tailor_denom;
        let credit = crate::match_engine::cleavage_credit_for(cand, params.enzyme, aa_set) as f32;
        f.rank_score_float =
            score_psm_float(&res_ss, &cand.peptide, scorer, z, fragment_tolerance_da) + credit;
        f.strong_score = fuse_strong_score(&StrongScoreInputs {
            intensity_signal: f.intensity_signal,
            chance_match_surprise: f.chance_match_surprise,
            mass_competition_evidence: f.mass_competition_evidence,
            candidate_rank_entropy: f.candidate_rank_entropy,
            listwise_score_gap: f.listwise_score_gap,
        });
        f
    };
    let retained_strong: Vec<f32> = retained
        .iter()
        .map(|p| make_features(p, false).strong_score)
        .collect();

    // 4. The winner: the top of the retained pool, or under `--score strong` the
    // highest strong score (earlier = higher rank_score wins ties), mirroring the
    // primary's `reorder_by_strong_score`.
    let winner_pos = if params.score_mode == crate::search_params::ScoreMode::Strong {
        (0..retained.len())
            .reduce(|a, b| {
                if retained_strong[b] > retained_strong[a] {
                    b
                } else {
                    a
                }
            })
            .unwrap_or(0)
    } else {
        0
    };
    let mut best = retained.swap_remove(winner_pos);
    let mut features = make_features(&best, true);
    features.strong_score_cal =
        strong_score_calibrated(features.strong_score, &retained_strong, &strong_null);
    if params.score_mode == crate::search_params::ScoreMode::Strong {
        // `score` stays the RawScore so RankScore keeps the primaries' scale.
        best.rank_score = features.strong_score;
    }
    let cand_peptide = &candidates[best.primary_candidate_idx() as usize].peptide;
    best.features = features;
    // The override makes the PIN writer compute ExpMass/dm/absdm from the
    // co-isolated mass.
    best.precursor_mz_override = Some(co.mono_mz);
    // Peaks this secondary explained (its charge-1 b/y matches on the FULL spectrum),
    // returned so the caller can remove them before the next co-isolated mass on
    // this scan is searched (sequential competition; see fn doc).
    let winner_claimed = primary_matched_peak_keys(spec, cand_peptide, scorer);
    Some((best, winner_claimed))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::search_index::SearchIndex;
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
    fn heavy_coisolated_envelope_is_reported_at_its_monoisotope() {
        // At ~3 kDa the M+1 peak is the most intense, so it is tried first; the
        // envelope must still be reported at M, not one isotope spacing heavier.
        let z = 2u8;
        let selected_mz = 1200.0;
        let sel_neutral = (selected_mz - PROTON) * z as f64;
        let co_mz = 1501.2;
        let co_neutral = (co_mz - PROTON) * z as f64;
        let env = averagine_isotope_envelope(co_neutral, 4);
        assert!(env[1] > env[0], "fixture needs M+1 > M");
        let mut peaks = envelope(selected_mz, z, sel_neutral, 1000.0);
        peaks.extend(envelope(co_mz, z, co_neutral, 500.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let got = detect_coisolated(&peaks, 1499.0, 1504.0, selected_mz, 2..=3, 0.01, 0.5, 4);
        assert_eq!(got.len(), 1, "one co-isolated envelope: {got:?}");
        assert!(
            (got[0].mono_mz - co_mz).abs() < 0.01,
            "reported at {} instead of the monoisotope {co_mz}",
            got[0].mono_mz
        );
        assert_eq!(got[0].charge, z);
        assert!((got[0].neutral_mass - co_neutral).abs() < 0.02);
    }

    #[test]
    fn max_n_zero_is_a_hard_disable() {
        // Same setup as above, but max_n=0 must emit ZERO secondaries (the cap is
        // checked BEFORE pushing — not "one then stop").
        let z = 2u8;
        let selected_mz = 600.0;
        let sel_neutral = (selected_mz - PROTON) * z as f64;
        let co_mz = 600.7;
        let co_neutral = (co_mz - PROTON) * z as f64;
        let mut peaks = envelope(selected_mz, z, sel_neutral, 1000.0);
        peaks.extend(envelope(co_mz, z, co_neutral, 500.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 0);
        assert!(
            got.is_empty(),
            "max_n=0 must emit zero secondaries, got {}",
            got.len()
        );
        // And max_n=1 still caps at exactly one.
        let one = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 1);
        assert_eq!(one.len(), 1, "max_n=1 caps at one");
    }

    #[test]
    fn no_coisolation_when_only_selected_present() {
        let z = 2u8;
        let selected_mz = 600.0;
        let peaks = envelope(selected_mz, z, (selected_mz - PROTON) * z as f64, 1000.0);
        let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
        assert!(
            got.is_empty(),
            "only the selected precursor -> no co-isolation"
        );
    }

    // ── Targeted second-peptide residual search ─────────────────────────────

    use crate::PreparedSearch;
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use model::protocol::Protocol;
    use model::Tolerance;
    use model::{AminoAcid, AminoAcidSetBuilder, Protein, ProteinDb};
    use rustc_hash::FxHashMap;
    use scoring_crate::param_model::{IonType, Partition, SpecDataType};
    use scoring_crate::scoring::fragment_ions::predict_by_ions;
    use scoring_crate::Param;

    /// Minimal RankScorer (mirrors `tests/match_engine_smoke.rs::tiny_scorer`):
    /// non-trivial prefix/suffix rank tables so b/y matches earn positive score.
    pub(crate) fn tiny_scorer() -> RankScorer {
        let part = Partition {
            charge: 2,
            parent_mass: 500.0,
            seg_num: 0,
        };
        let prefix1 = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let suffix1 = IonType::Suffix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
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
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        };
        param.rebuild_cache();
        RankScorer::new(&param)
    }

    /// Every candidate of the prepared out-of-core index and its nominal-mass
    /// bucket index: the shape of one Pass-2 pool covering the whole database.
    fn whole_db_pool(prepared: &PreparedSearch<'_>) -> (Vec<Candidate>, BTreeMap<i32, Vec<usize>>) {
        let pool = prepared.out_of_core_candidates_where(|_, _| true);
        let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
        for (i, c) in pool.iter().enumerate() {
            bucket_index
                .entry(c.peptide.nominal_residue_mass())
                .or_default()
                .push(i);
        }
        (pool, bucket_index)
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

        let prepared =
            PreparedSearch::prepare_in_memory(&idx, &params, &scorer, frag_tol, "XXX").unwrap();
        let (candidates, bucket_index) = whole_db_pool(&prepared);

        // Planted secondary peptide and its co-isolated precursor.
        let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
        let z = 2u8;
        let co_neutral = planted.mass();
        let co_mz = (co_neutral + z as f64 * PROTON) / z as f64;
        let co = CoIsolated {
            mono_mz: co_mz,
            charge: z,
            neutral_mass: co_neutral,
        };

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
            isolation_target_mz: None,
        };

        let got = search_secondary(
            &spec,
            &primary,
            &std::collections::HashSet::new(),
            co,
            &candidates,
            &bucket_index,
            &scorer,
            &prepared.aa_set_for_scoring,
            &params,
            frag_tol,
        );

        let (psm, winner_claimed) =
            got.expect("secondary search should return a PSM at the co-isolated mass");
        let found = &candidates[psm.primary_candidate_idx() as usize].peptide;
        assert_eq!(
            found.residues, planted.residues,
            "secondary PSM should resolve to the planted peptide DAFLGSFLYEYSR"
        );

        // Competition check: searching the SAME co-isolated mass again with the
        // winner's peaks already claimed strips them from the residual, so the
        // peptide can no longer earn its fragment matches — it must NOT come back
        // with the same evidence. (Proves prior_claimed removes shared peaks so two
        // co-isolated precursors can't double-count the same leftover signal.)
        assert!(
            !winner_claimed.is_empty(),
            "winner should claim its matched peaks"
        );
        let again = search_secondary(
            &spec,
            &primary,
            &winner_claimed,
            co,
            &candidates,
            &bucket_index,
            &scorer,
            &prepared.aa_set_for_scoring,
            &params,
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
        // The planted peptide is in-window for the co-isolated mass and matched
        // its fragment ions, so it must earn a positive RawScore (rank_score).
        assert!(
            psm.rank_score > 0.0,
            "secondary PSM rank_score (RawScore) must be positive, got {}",
            psm.rank_score
        );
    }

    /// Spectrum whose peaks are `planted`'s charge-1 b/y ions, with `primary`'s
    /// precursor m/z, plus the co-isolated precursor at `planted`'s mass.
    fn planted_secondary_spectrum(planted: &Peptide, primary: &Peptide) -> (Spectrum, CoIsolated) {
        let z = 2u8;
        let mut peaks: Vec<(f64, f32)> = predict_by_ions(planted, 1..=1)
            .into_iter()
            .map(|p| (p.mz, 100.0_f32))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let spec = Spectrum {
            title: "chimeric-planted".into(),
            precursor_mz: (primary.mass() + z as f64 * PROTON) / z as f64,
            precursor_intensity: None,
            precursor_charge: Some(z as i32),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        };
        let co = CoIsolated {
            mono_mz: (planted.mass() + z as f64 * PROTON) / z as f64,
            charge: z,
            neutral_mass: planted.mass(),
        };
        (spec, co)
    }

    #[test]
    fn secondary_listwise_gap_ignores_duplicate_occurrences_of_the_winner() {
        // DAFLGSFLYEYSR occurs in two proteins with the same flanks, so the pool
        // holds it twice at the same RawScore. The listwise gap must be taken
        // against the best DIFFERENT peptidoform (here its isobaric decoy), not
        // against the winner's own second occurrence.
        let target = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: "".into(),
                    sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
                },
                Protein {
                    accession: "P2".into(),
                    description: "".into(),
                    sequence: b"GGGKDAFLGSFLYEYSRKGGG".to_vec(),
                },
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = tiny_scorer();
        let frag_tol = 0.02_f64;
        let prepared =
            PreparedSearch::prepare_in_memory(&idx, &params, &scorer, frag_tol, "XXX").unwrap();
        let (candidates, bucket_index) = whole_db_pool(&prepared);

        let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
        let occurrences = candidates
            .iter()
            .filter(|c| !c.is_decoy && c.peptide.residues == planted.residues)
            .count();
        assert!(occurrences >= 2, "fixture must hold the winner twice, got {occurrences}");

        let primary = peptide(b"WVTFISLLR", b'K', b'K');
        let (spec, co) = planted_secondary_spectrum(&planted, &primary);
        let (psm, _) = search_secondary(
            &spec,
            &primary,
            &std::collections::HashSet::new(),
            co,
            &candidates,
            &bucket_index,
            &scorer,
            &prepared.aa_set_for_scoring,
            &params,
            frag_tol,
        )
        .expect("planted secondary must be found");
        let winner = &candidates[psm.primary_candidate_idx() as usize];
        assert!(!winner.is_decoy);
        assert_eq!(winner.peptide.residues, planted.residues);
        assert!(
            psm.features.listwise_score_gap > 0.0,
            "gap must be measured against a different peptidoform, got {}",
            psm.features.listwise_score_gap
        );
    }

    #[test]
    fn secondary_features_share_the_primary_scales() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            }],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let scorer = tiny_scorer();
        let frag_tol = 0.02_f64;
        let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
        let primary = peptide(b"WVTFISLLR", b'K', b'K');
        let (spec, co) = planted_secondary_spectrum(&planted, &primary);

        let run = |mode: crate::search_params::ScoreMode| {
            let mut params = SearchParams::default_tryptic(aa_set.clone());
            params.score_mode = mode;
            let prepared =
                PreparedSearch::prepare_in_memory(&idx, &params, &scorer, frag_tol, "XXX").unwrap();
            let (candidates, bucket_index) = whole_db_pool(&prepared);
            search_secondary(
                &spec,
                &primary,
                &std::collections::HashSet::new(),
                co,
                &candidates,
                &bucket_index,
                &scorer,
                &prepared.aa_set_for_scoring,
                &params,
                frag_tol,
            )
            .expect("planted secondary must be found")
            .0
        };
        let rank = run(crate::search_params::ScoreMode::Rank);
        let f = &rank.features;
        // RankScoreFloat is the unrounded RawScore, never left at 0.
        assert!(f.rank_score_float > 0.0, "RankScoreFloat {}", f.rank_score_float);
        assert!((f.rank_score_float - rank.score).abs() < 5.0);
        // The planted peptide and its isobaric decoy are both retained, so
        // RawScoreCal is the leave-one-out z-score over the retained strong
        // scores, not the raw strong score.
        assert!(f.listwise_score_gap > 0.0);
        assert!(f.strong_score_cal.is_finite());
        assert_ne!(f.strong_score_cal, f.strong_score);

        // Under --score strong, RankScore (`score`) keeps the RawScore scale; only
        // the ranking field carries the strong score.
        let strong = run(crate::search_params::ScoreMode::Strong);
        assert_eq!(strong.score, rank.score);
        assert_eq!(strong.rank_score, strong.features.strong_score);
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
        let prepared =
            PreparedSearch::prepare_in_memory(&idx, &params, &scorer, frag_tol, "XXX").unwrap();
        let (candidates, bucket_index) = whole_db_pool(&prepared);

        let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
        let z = 2u8;
        // Raw observed mass biased high by `shift_ppm` (so adjusted == true mass).
        let true_mass = planted.mass();
        let raw_neutral = true_mass / (1.0 - shift_ppm * 1e-6);
        let co_mz = (raw_neutral + z as f64 * PROTON) / z as f64;
        let co = CoIsolated {
            mono_mz: co_mz,
            charge: z,
            neutral_mass: raw_neutral,
        };

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
            isolation_target_mz: None,
        };

        let got = search_secondary(
            &spec,
            &primary,
            &std::collections::HashSet::new(),
            co,
            &candidates,
            &bucket_index,
            &scorer,
            &prepared.aa_set_for_scoring,
            &params,
            frag_tol,
        );
        let (psm, _claimed) =
            got.expect("secondary must be found at the calibration-adjusted mass");
        let found = &candidates[psm.primary_candidate_idx() as usize].peptide;
        assert_eq!(
            found.residues, planted.residues,
            "should resolve to the planted peptide"
        );
        // The candidate matches the ADJUSTED mass, so the reported error is ~0,
        // not ~-30 ppm (which is what the un-calibrated raw mass would report).
        assert!(
            psm.mass_error_ppm.abs() < 5.0,
            "mass error must be reported on the calibrated scale (~0), got {} ppm",
            psm.mass_error_ppm
        );
    }
}
