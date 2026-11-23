//! Co-isolated precursors: detecting a second precursor in a spectrum's isolation window,
//! searching its mass window, and the second (co-isolation) pass that adds those matches.
//! Part of the match engine; `PreparedSearch` runs the pass.

use super::*;
use crate::candidates::Candidate;
use crate::chimeric::precursor_isotope_match;
use crate::precursor::calibration::adjusted_observed_neutral_mass;
use crate::psm::{PsmFeatures, PsmMatch};
use crate::params::SearchParams;
use model::amino_acid::set::AminoAcidSet;
use model::mass::{nominal_from, H2O, ISOTOPE, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::scoring::{
    fuse_strong_score, predict_by_ions, psm_edge_score, score_psm, score_psm_float,
    strong_score_calibrated, OnlineStats, RankScorer, ScoredSpectrum, StrongScoreInputs,
};
/// The isolation window around its target m/z (the precursor m/z when the scan records no
/// target), and the precursor tolerance in Da.
pub(super) fn coisolation_window(spec: &Spectrum, params: &SearchParams) -> (f64, f64, f64) {
    let (lo, hi) = crate::precursor::matching::isolation_window_mz(
        spec,
        params.chimeric_isolation_halfwidth_da,
    );
    let tol = params
        .precursor_tolerance
        .left
        .as_da(spec.precursor_mz)
        .max(0.01);
    (lo, hi, tol)
}

pub fn coisolation_ms1_bounds(spec: &Spectrum, params: &SearchParams) -> (f64, f64) {
    const MARGIN_DA: f64 = 0.01;
    let (lo, hi, tol) = coisolation_window(spec, params);
    (
        lo - tol - MARGIN_DA,
        hi + 3.0 * model::mass::ISOTOPE + tol + MARGIN_DA,
    )
}

pub(super) struct Pass2Scan {
    primary: Peptide,
    coisolated: Vec<CoIsolated>,
}

pub fn coisolation_pass(
    prepared: &PreparedSearch,
    spectra: &[Spectrum],
    queues: &mut [TopNQueue],
    params: &SearchParams,
    link: &Ms1Link,
    offset: usize,
) {
    if !params.chimeric {
        return;
    }
    let shifted;
    let spectra = match recentred_spectra(spectra, prepared.scorer.fragment_recentre_ppm()) {
        Some(v) => {
            shifted = v;
            &shifted[..]
        }
        None => spectra,
    };
    let scorer_set = prepared.scorer_set(spectra);

    let scans: Vec<Option<Pass2Scan>> = queues
        .par_iter()
        .enumerate()
        .map(|(spec_idx, q)| {
            if q.is_empty() {
                return None;
            }
            let spec = spectra.get(spec_idx)?;
            let ms1_idx = (*link.ms2_to_ms1.get(spec_idx)?)?;
            let ms1 = link.ms1_peaks.get(ms1_idx)?;
            let (lo, hi, tol) = coisolation_window(spec, params);
            let coisolated = detect_coisolated(
                ms1,
                lo,
                hi,
                spec.precursor_mz,
                *params.charge_range.start()..=*params.charge_range.end(),
                tol,
                params.chimeric_max_kl,
                params.chimeric_max_coisolated,
            );
            if coisolated.is_empty() {
                return None;
            }
            let primary = prepared.candidate_peptide(q.peek_top()?.primary_candidate_idx());
            Some(Pass2Scan {
                primary,
                coisolated,
            })
        })
        .collect();

    let found = mmap_secondaries(prepared, &scorer_set, spectra, &scans, params);

    let mut accum = prepared.mmap_accum.lock().unwrap();
    for (spec_idx, (q, secondaries)) in queues.iter_mut().zip(found).enumerate() {
        if secondaries.is_empty() {
            continue;
        }
        let mut emitted: std::collections::BTreeSet<(PepDedupKey, u8, bool)> = q
            .iter_psms()
            .filter_map(|p| {
                let i = p.primary_candidate_idx() as usize;
                let c = accum
                    .candidates
                    .get(i)
                    .or_else(|| prepared.candidates.get(i))?;
                Some((
                    PepDedupKey::from_peptide(&c.peptide),
                    p.charge_used,
                    c.is_decoy,
                ))
            })
            .collect();
        for (mut psm, cand) in secondaries {
            let key = (
                PepDedupKey::from_peptide(&cand.peptide),
                psm.charge_used,
                cand.is_decoy,
            );
            if !emitted.insert(key) {
                continue;
            }
            psm.candidate_idxs = smallvec::smallvec![accum.intern(&cand)];
            psm.spectrum_idx = offset + spec_idx;
            q.force_push(psm);
        }
    }
}

pub(super) const PASS2_POOL_SPAN_NOMINAL: i32 = 32;

pub(super) fn mmap_secondaries(
    prepared: &PreparedSearch,
    scorer_set: &ScorerSet<'_>,
    spectra: &[Spectrum],
    scans: &[Option<Pass2Scan>],
    params: &SearchParams,
) -> Vec<Vec<(PsmMatch, Candidate)>> {
    type Hit = (PsmMatch, Candidate, std::collections::HashSet<i64>);
    type Item = (usize, CoIsolated, (i32, i32));

    let mut claimed: Vec<std::collections::HashSet<i64>> =
        vec![std::collections::HashSet::new(); scans.len()];
    let mut found: Vec<Vec<(PsmMatch, Candidate)>> = vec![Vec::new(); scans.len()];
    let rounds = scans
        .iter()
        .flatten()
        .map(|s| s.coisolated.len())
        .max()
        .unwrap_or(0);
    for k in 0..rounds {
        let mut items: Vec<Item> = scans
            .iter()
            .enumerate()
            .filter_map(|(i, s)| {
                let co = *s.as_ref()?.coisolated.get(k)?;
                Some((i, co, secondary_nominal_window(co, params)))
            })
            .collect();
        items.sort_by_key(|&(i, _, w)| (w, i));
        let mut batches: Vec<&[Item]> = Vec::new();
        let mut rest = items.as_slice();
        while let Some(&(_, _, (lo, _))) = rest.first() {
            let n = rest
                .iter()
                .take_while(|&&(_, _, (_, hi))| hi - lo <= PASS2_POOL_SPAN_NOMINAL)
                .count()
                .max(1);
            let (batch, tail) = rest.split_at(n);
            batches.push(batch);
            rest = tail;
        }
        let claimed_ref = &claimed;
        let hits: Vec<(usize, Option<Hit>)> = batches
            .par_iter()
            .flat_map_iter(|batch| {
                let mut windows: Vec<(i32, i32)> = Vec::new();
                for &(_, _, (lo, hi)) in batch.iter() {
                    match windows.last_mut() {
                        Some(last) if lo <= last.1 + 1 => last.1 = last.1.max(hi),
                        _ => windows.push((lo, hi)),
                    }
                }
                let pool = prepared.expand_mmap_window_candidates(prepared.params, &windows, None);
                let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
                for (i, c) in pool.iter().enumerate() {
                    bucket_index
                        .entry(c.peptide.nominal_residue_mass())
                        .or_default()
                        .push(i);
                }
                batch
                    .iter()
                    .map(|&(spec_idx, co, _)| {
                        let scan = scans[spec_idx]
                            .as_ref()
                            .expect("only scans with co-isolated precursors are batched");
                        let hit = search_secondary(
                            &spectra[spec_idx],
                            &scan.primary,
                            &claimed_ref[spec_idx],
                            co,
                            &pool,
                            &bucket_index,
                            scorer_set.for_spec(&spectra[spec_idx]),
                            &prepared.aa_set_for_scoring,
                            params,
                            prepared.fragment_tolerance_da,
                        )
                        .map(|(psm, winner_claimed)| {
                            let cand = pool[psm.primary_candidate_idx() as usize].clone();
                            (psm, cand, winner_claimed)
                        });
                        (spec_idx, hit)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        for (spec_idx, hit) in hits {
            if let Some((psm, cand, winner_claimed)) = hit {
                claimed[spec_idx].extend(winner_claimed);
                found[spec_idx].push((psm, cand));
            }
        }
    }
    found
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CoIsolated {
    pub mono_mz: f64,
    pub charge: u8,
    pub neutral_mass: f64,
}

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
        if out.len() >= max_n {
            break;
        }
        if (mz - selected_mz).abs() <= tol_da {
            continue;
        }
        if is_selected_isotope(mz) {
            continue;
        }
        if out.iter().any(|c| {
            let d = (mz - c.mono_mz).abs();
            (0..6).any(|k| (d - k as f64 * ISOTOPE / c.charge as f64).abs() <= tol_da)
        }) {
            continue;
        }
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

const MAX_MONO_WALK: usize = 3;

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
        let tol_da = scorer.node_score_window().as_da(p.mz);
        let lo_mz = p.mz - tol_da;
        let hi_mz = p.mz + tol_da;
        let start = spec.peaks.partition_point(|&(mz, _)| mz < lo_mz);
        let mut best: Option<(f64, f32)> = None;
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

fn rank_desc_then_idx(a: &PsmMatch, b: &PsmMatch) -> std::cmp::Ordering {
    b.cmp(a)
}

pub(crate) fn secondary_mass_window(co: CoIsolated, params: &SearchParams) -> (f64, f64, f64) {
    let co_neutral =
        adjusted_observed_neutral_mass(co.neutral_mass, params.precursor_mass_shift_ppm);
    let left = params.precursor_tolerance.left.as_da(co_neutral).max(0.01);
    let right = params.precursor_tolerance.right.as_da(co_neutral).max(0.01);
    (co_neutral, left, right)
}

fn in_secondary_mass_window(peptide_mass: f64, co_neutral: f64, left: f64, right: f64) -> bool {
    let d = peptide_mass - co_neutral;
    d >= -left && d <= right
}

pub(crate) fn secondary_nominal_window(co: CoIsolated, params: &SearchParams) -> (i32, i32) {
    let (co_neutral, left, right) = secondary_mass_window(co, params);
    let nominal = |m: f64| nominal_from(m - H2O);
    (
        nominal(co_neutral - left) - 1,
        nominal(co_neutral + right) + 1,
    )
}

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

    let mut claimed = primary_matched_peak_keys(spec, primary, scorer);
    claimed.extend(prior_claimed.iter().copied());
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
        mobility_window: spec.mobility_window,
    };
    let res_ss = ScoredSpectrum::new(&co_spec, scorer, z);

    let (co_neutral, tol_left, tol_right) = secondary_mass_window(co, params);
    let (lo, hi) = secondary_nominal_window(co, params);

    let mut tailor_hist: std::collections::HashMap<i32, u32> = std::collections::HashMap::new();
    let mut tailor_total: u32 = 0;
    let mut best_raw: i32 = i32::MIN;
    let mut best_pep_key: u64 = 0;
    let mut have_best = false;
    let mut second_raw: i32 = i32::MIN;
    let mut strong_null = OnlineStats::default();
    let mut distinct: rustc_hash::FxHashMap<(u64, bool), PsmMatch> =
        rustc_hash::FxHashMap::default();
    for (_nm, idxs) in bucket_index.range(lo..=hi) {
        for &ci in idxs {
            let cand = &candidates[ci];
            if !in_secondary_mass_window(cand.peptide.mass(), co_neutral, tol_left, tol_right) {
                continue;
            }
            if cand.peptide == *primary {
                continue;
            }
            let pin = score_psm(&res_ss, &cand.peptide, scorer, z, fragment_tolerance_da);
            let edge = psm_edge_score(&res_ss, &cand.peptide, scorer, z);
            let cleavage = crate::engine::cleavage_credit_for(cand, params.enzyme, aa_set);
            let raw_score = pin + cleavage as f32;

            strong_null.push(raw_score);
            let s = raw_score.round() as i32;
            *tailor_hist.entry(s).or_insert(0) += 1;
            tailor_total += 1;
            let mkey = crate::engine::peptidoform_key(&cand.peptide);
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
                candidate_idxs: smallvec::smallvec![ci as u32],
                charge_used: z,
                mass_error_ppm: (cand.peptide.mass() - co_neutral) / co_neutral * 1e6,
                score: raw_score,
                rank_score: raw_score + edge as f32,
                edge_score: edge,
                activation_method: Some(scorer.param().data_type.activation),
                features: PsmFeatures::default(),
                isotope_offset: 0,
                precursor_mz_override: None,
                site_evidence: None,
                tie_key: crate::psm::tie_key_of(cand),
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

    let mut retained: Vec<PsmMatch> = distinct.into_values().collect();
    retained.sort_by(rank_desc_then_idx);
    retained.truncate(params.hot_loop_retention_cap().max(1) as usize);

    let mut retained_scores: Vec<f32> = retained.iter().map(|p| p.score).collect();
    retained_scores.sort_by(|a, b| b.total_cmp(a));
    let listwise_score_gap = scoring_crate::scoring::listwise_score_gap(&retained_scores);
    let candidate_rank_entropy = scoring_crate::scoring::candidate_rank_entropy(&retained_scores);
    let tailor_denom = crate::psm::tailor_denominator(&tailor_hist, tailor_total) as f32;
    let delta_raw_score = if second_raw > i32::MIN {
        (best_raw - second_raw) as f32
    } else {
        0.0
    };

    let make_features = |psm: &PsmMatch, full: bool| -> PsmFeatures {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let mut f = crate::engine::compute_psm_features_with_edge(
            &res_ss,
            &cand.peptide,
            scorer,
            z,
            Some(psm.edge_score),
            full,
            None,
        );
        f.delta_raw_score = delta_raw_score;
        f.listwise_score_gap = listwise_score_gap;
        f.candidate_rank_entropy = candidate_rank_entropy;
        f.tailor_score = psm.score / tailor_denom;
        let credit = crate::engine::cleavage_credit_for(cand, params.enzyme, aa_set) as f32;
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

    let winner_pos = if params.score_mode == crate::params::ScoreMode::Strong {
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
    if params.score_mode == crate::params::ScoreMode::Strong {
        best.rank_score = features.strong_score;
    }
    let cand_peptide = &candidates[best.primary_candidate_idx() as usize].peptide;
    best.features = features;
    best.precursor_mz_override = Some(co.mono_mz);
    let winner_claimed = primary_matched_peak_keys(spec, cand_peptide, scorer);
    Some((best, winner_claimed))
}

#[cfg(test)]
#[path = "coisolation/tests/coisolation.rs"]
pub(crate) mod tests;
