//! Searching a chunk of spectra: candidate windows, scoring, the kept features, and the inner loop
//! (`run_chunk_inner`) that fills each spectrum's top-N queue.

use std::sync::atomic::{AtomicUsize, Ordering};

use model::mass::PROTON;
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use scoring_crate::psm::{
    edge_scoring_enabled, fuse_strong_score, psm_edge_score, score_psm, score_psm_float,
    scored_for_charge, strong_score_calibrated, OnlineStats, ScoredSpectrum, StrongScoreInputs,
};
use smallvec::SmallVec;

use super::dedup::*;
use super::features::*;
use super::{
    charges_to_try, enumeration_pass, mmap_window_key, prefetch_window_candidates,
    recalibrate_strong_over_spectrum, recentred_spectra, windowed_top_k, MmapRecordCache,
    PreparedSearch, RankScorer, RefTracker, ScorerSet, KEPT_FEATURE_BATCH_SPECTRA,
};
use crate::candidates::{
    base_record_key, base_records_for_nominal_window, distinct_base_records_for_nominal_window,
    expand_base_record, retain_records_with_forms_in, BaseRecordKey, Candidate, FormWalkTables,
};
use crate::fragments::index::ChunkFragmentIndex;
use crate::params::{ScoreMode, SearchParams};
use crate::precursor::calibration::adjusted_observed_neutral_mass;
use crate::precursor::matching::{matches_precursor, MassError};
use crate::psm::{PsmFeatures, PsmMatch, TopNQueue};

/// What every spectrum of a chunk shares while it is scored: the scorers, the RT source, the
/// chunk's candidate records or fragment index, and the search parameters.
#[derive(Clone, Copy)]
struct ChunkScoring<'c> {
    scorer_set: &'c ScorerSet<'c>,
    rt_source: &'c Option<crate::fragments::index::FormRtSource<'c>>,
    record_cache: &'c Option<MmapRecordCache>,
    frag_index: Option<&'c ChunkFragmentIndex>,
    params: &'c SearchParams,
    spectrum_idx_offset: usize,
    complete_later: bool,
}

/// A spectrum's DIA isolation windows: (charge, lowest and highest neutral mass) per charge.
type DiaWindows = SmallVec<[(u8, f64, f64); 4]>;

/// The scored spectrum of charge `z`; every candidate charge has one.
fn scored_spec<'a, 's>(
    scored_per_charge: &'a [(u8, ScoredSpectrum<'s>)],
    z: u8,
) -> &'a ScoredSpectrum<'s> {
    scored_for_charge(scored_per_charge, z).expect("scored spectrum exists for candidate charge")
}

/// Moves the PSMs and pool scores of one charge search's units (`got`, unit `k` of spectrum
/// `owners[k]`) into their spectra's merged queues and pools, the PSMs numbered from
/// `offset`. When the units belong to distinct spectra (one spectrum's charges fall in
/// different searches) the spectra are filled in parallel, each in the same order as the
/// sequential fill; otherwise unit by unit.
fn merge_units(
    got: Vec<(TopNQueue, Vec<f32>)>,
    owners: &[usize],
    merged: &mut [TopNQueue],
    pools: &mut [Vec<f32>],
    offset: usize,
) {
    let fill = |(mut q, pool): (TopNQueue, Vec<f32>),
                i: usize,
                queue: &mut TopNQueue,
                kept: &mut Vec<f32>| {
        kept.extend(pool);
        for mut psm in q.drain_into_vec() {
            psm.spectrum_idx = offset + i;
            queue.push(psm);
        }
    };
    let mut seen = vec![false; merged.len()];
    if !owners
        .iter()
        .all(|&i| !std::mem::replace(&mut seen[i], true))
    {
        for (unit, &i) in got.into_iter().zip(owners) {
            fill(unit, i, &mut merged[i], &mut pools[i]);
        }
        return;
    }
    let mut slot: Vec<Option<(TopNQueue, Vec<f32>)>> = (0..merged.len()).map(|_| None).collect();
    for (unit, &i) in got.into_iter().zip(owners) {
        slot[i] = Some(unit);
    }
    merged
        .par_iter_mut()
        .zip(pools.par_iter_mut())
        .zip(slot.into_par_iter())
        .enumerate()
        .for_each(|(i, ((queue, kept), unit))| {
            if let Some(unit) = unit {
                fill(unit, i, queue, kept);
            }
        });
}

/// What finishing a spectrum's scoring carries over from its candidate loop: the scorer and scored
/// spectra, the candidates, the top-N queue and the per-charge queues, and the loop's statistics.
struct SpectrumTail<'a, 's> {
    scorer: &'a RankScorer,
    scored_per_charge: &'a [(u8, ScoredSpectrum<'s>)],
    cand_slice: &'a [Candidate],
    queue: TopNQueue,
    per_charge_queues: FxHashMap<u8, TopNQueue>,
    best_raw: i32,
    second_raw: i32,
    ref_tracker: RefTracker,
    track_ref: bool,
    tailor_hist: FxHashMap<i32, u32>,
    tailor_total: u32,
    strong_null_stats: OnlineStats,
    open_peaks: Option<Vec<(f64, f32)>>,
    open_attempts: usize,
    open_wins: usize,
}

impl<'a> PreparedSearch<'a> {
    pub fn run_chunk(&self, spectra: &[Spectrum], spectrum_idx_offset: usize) -> Vec<TopNQueue> {
        self.run_chunk_with_params(spectra, spectrum_idx_offset, self.params)
    }

    /// Scores a chunk of spectra; each returned queue holds exactly its PSMs, without the spare
    /// capacity reserved while scoring (the queues of a whole run are kept until the PIN is
    /// written).
    pub fn run_chunk_with_params(
        &self,
        spectra: &[Spectrum],
        spectrum_idx_offset: usize,
        params: &SearchParams,
    ) -> Vec<TopNQueue> {
        let mut queues = self.score_chunk(spectra, spectrum_idx_offset, params);
        for q in &mut queues {
            q.shrink_to_fit();
        }
        queues
    }

    fn score_chunk(
        &self,
        spectra: &[Spectrum],
        spectrum_idx_offset: usize,
        params: &SearchParams,
    ) -> Vec<TopNQueue> {
        let shifted;
        let spectra = match recentred_spectra(spectra, self.scorer.fragment_recentre_ppm()) {
            Some(v) => {
                shifted = v;
                &shifted[..]
            }
            None => spectra,
        };
        let has_charge = |s: &Spectrum| s.precursor_charge.is_some_and(|z| z > 0);
        if params.fragment_index_top_k == 0
            || params.charge_range.clone().count() < 2
            || spectra.iter().all(has_charge)
        {
            let queues = self
                .run_chunk_inner(spectra, spectrum_idx_offset, params, false)
                .into_iter()
                .map(|(queue, _pool)| queue)
                .collect();
            self.frag_pred_cache.clear();
            return queues;
        }
        let mass = crate::fragments::index::index_order_mass;
        let (lo, hi) = spectra
            .iter()
            .map(mass)
            .fold((f64::MAX, f64::MIN), |(lo, hi), m| (lo.min(m), hi.max(m)));
        let span_da = (hi - lo).max(1.0);
        let mut units: Vec<(usize, Spectrum)> = Vec::with_capacity(spectra.len());
        for (i, s) in spectra.iter().enumerate() {
            if has_charge(s) {
                units.push((i, s.clone()));
            } else {
                for z in params.charge_range.clone() {
                    let mut unit = s.clone();
                    unit.precursor_charge = Some(z as i32);
                    units.push((i, unit));
                }
            }
        }
        units.sort_by(|a, b| mass(&a.1).total_cmp(&mass(&b.1)).then(a.0.cmp(&b.0)));
        let mut merged: Vec<TopNQueue> = (0..spectra.len())
            .map(|_| TopNQueue::new(params.top_n_psms_per_spectrum))
            .collect();
        let mut pools: Vec<Vec<f32>> = vec![Vec::new(); spectra.len()];
        let complete_later = params.score_mode == ScoreMode::Strong && params.open_search.is_none();
        while !units.is_empty() {
            let take = crate::fragments::index::index_chunk_len(&units, |u| mass(&u.1), span_da);
            let (owners, chunk): (Vec<usize>, Vec<Spectrum>) = units.drain(..take).unzip();
            let got = self.run_chunk_inner(&chunk, 0, params, complete_later);
            merge_units(got, &owners, &mut merged, &mut pools, spectrum_idx_offset);
        }
        merged
            .par_iter_mut()
            .zip(pools.par_iter())
            .for_each(|(queue, pool)| {
                queue.trim_to_capacity(params.top_n_psms_per_spectrum);
                recalibrate_strong_over_spectrum(queue, pool);
            });
        if complete_later {
            self.complete_kept_features(spectra, &mut merged);
        }
        self.frag_pred_cache.clear();
        merged
    }

    /// Completes the features of the PSMs each spectrum keeps after its charges were merged
    /// (the charge searches of [`Self::score_chunk`] leave the full-pass-only features unset),
    /// so only the kept PSMs pay for them. Peptides are copied out of the candidate store in
    /// batches of spectra to keep the lock short.
    fn complete_kept_features(&self, spectra: &[Spectrum], queues: &mut [TopNQueue]) {
        let scorer_set = self.scorer_set(&[]);
        for (specs, batch) in spectra
            .chunks(KEPT_FEATURE_BATCH_SPECTRA)
            .zip(queues.chunks_mut(KEPT_FEATURE_BATCH_SPECTRA))
        {
            let peptides: Vec<Vec<Peptide>> = {
                let guard = self.mmap_accum.lock().unwrap();
                let accum = &*guard;
                batch
                    .par_iter()
                    .map(|q| {
                        q.iter_psms()
                            .map(|p| {
                                let idx = p.primary_candidate_idx() as usize;
                                accum
                                    .candidates
                                    .get(idx)
                                    .unwrap_or_else(|| &self.candidates[idx])
                                    .peptide
                                    .clone()
                            })
                            .collect()
                    })
                    .collect()
            };
            batch
                .par_iter_mut()
                .zip(specs.par_iter())
                .zip(peptides.par_iter())
                .for_each(|((queue, spec), peps)| {
                    if queue.is_empty() {
                        return;
                    }
                    let scorer = scorer_set.for_spec(spec);
                    let mut scored: SmallVec<[(u8, ScoredSpectrum<'_>); 4]> = SmallVec::new();
                    let mut k = 0;
                    queue.update_in_place(|psm| {
                        let z = psm.charge_used;
                        if scored.iter().all(|(c, _)| *c != z) {
                            scored.push((z, ScoredSpectrum::new(spec, scorer, z)));
                        }
                        let ss = &scored.iter().find(|(c, _)| *c == z).expect("scored").1;
                        complete_psm_features(
                            &mut psm.features,
                            ss,
                            &peps[k],
                            scorer,
                            z,
                            Some(&self.frag_pred_cache),
                        );
                        k += 1;
                    });
                });
        }
    }

    pub(super) fn expand_mmap_window_candidates(
        &self,
        params: &SearchParams,
        windows: &[(i32, i32)],
        cache: Option<&MmapRecordCache>,
    ) -> Vec<Candidate> {
        let mi = &self.mmap_index;

        let mut records: FxHashMap<BaseRecordKey, (crate::candidates::records::IndexRecord, u32)> =
            FxHashMap::default();
        for &(min_nominal, max_nominal) in windows {
            for (rec, mult) in base_records_for_nominal_window(mi, params, min_nominal, max_nominal)
            {
                records.entry(base_record_key(&rec)).or_insert((rec, mult));
            }
        }
        let mut records: Vec<(crate::candidates::records::IndexRecord, u32)> =
            records.into_values().collect();
        records.sort_unstable_by_key(|(r, _)| base_record_key(r));

        let in_any_window = |c: &Candidate| {
            let key = c.peptide.nominal_residue_mass();
            windows.iter().any(|&(lo, hi)| key >= lo && key <= hi)
        };

        let mut tagged: Vec<(u32, u32, u32, u16, Candidate)> = Vec::new();
        let mut owned: Vec<Candidate>;
        for (rec, mult) in &records {
            let forms: &[Candidate] = match cache.and_then(|c| c.get(&base_record_key(rec))) {
                Some(cached) => cached.as_slice(),
                None => {
                    owned = expand_base_record(self.idx, params, rec);
                    owned.as_slice()
                }
            };
            for cand in forms.iter().filter(|c| in_any_window(c)) {
                for k in 0..*mult {
                    let pass = enumeration_pass(cand, k);
                    tagged.push((
                        cand.protein_index as u32,
                        pass,
                        cand.start_offset_in_protein as u32,
                        cand.peptide.residues.len() as u16,
                        cand.clone(),
                    ));
                }
            }
        }
        tagged.sort_by(|a, b| (a.0, a.1, a.2, a.3).cmp(&(b.0, b.1, b.2, b.3)));
        let mut mmap_cands: Vec<Candidate> = tagged.into_iter().map(|t| t.4).collect();

        relabel_collision_decoys_with(&mut mmap_cands, &self.target_bare_seqs);
        if self.candidate_filter.is_some() {
            mmap_cands.retain(|c| self.keeps(c));
        }
        mmap_cands
    }

    /// Scores each spectrum of `spectra` at every charge it is searched with and returns, per
    /// spectrum, its retained PSMs and the strong scores of the pool its `strong_score_cal` was
    /// calibrated against (empty when the spectrum was skipped).
    /// With `complete_later` (strong mode) the full-pass-only features of the retained PSMs are
    /// left unset for the caller to complete with [`complete_psm_features`] once it knows which
    /// PSMs it keeps.
    fn run_chunk_inner(
        &self,
        spectra: &[Spectrum],
        spectrum_idx_offset: usize,
        params: &SearchParams,
        complete_later: bool,
    ) -> Vec<(TopNQueue, Vec<f32>)> {
        let scorer_set = self.scorer_set(spectra);

        let t_chunk_prep = std::time::Instant::now();
        let rt_prior = self.dia_rt_prior.read().expect("RT prior lock").clone();
        let rt_source = match (
            params.dia,
            self.form_rt_model.as_deref(),
            rt_prior.as_deref(),
        ) {
            (true, Some(model), Some(prior)) => Some(crate::fragments::index::FormRtSource {
                db: self.idx,
                params,
                model,
                cache: &self.form_rt_cache,
                prior,
            }),
            _ => None,
        };
        let owned_index: Option<ChunkFragmentIndex>;
        let (record_cache, frag_index): (Option<MmapRecordCache>, Option<&ChunkFragmentIndex>) = {
            let mi = &self.mmap_index;
            let mut windows: Vec<(i32, i32)> = spectra
                .iter()
                .filter(|s| s.peaks.len() >= params.min_peaks as usize)
                .flat_map(|s| mmap_window_key(s, params))
                .collect();
            windows.sort_unstable();
            windows.dedup();
            if let Some(open) = self.open_index.as_deref() {
                (Some(MmapRecordCache::default()), Some(open))
            } else if params.fragment_index_top_k > 0 {
                owned_index = Some({
                    let shift_ppm = params.precursor_mass_shift_ppm;
                    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
                    for s in spectra
                        .iter()
                        .filter(|s| s.peaks.len() >= params.min_peaks as usize)
                    {
                        for z in charges_to_try(s, params) {
                            if params.dia {
                                let (a, b) = crate::precursor::matching::dia_neutral_window(
                                    s,
                                    z,
                                    &params.precursor_tolerance,
                                    params.chimeric_isolation_halfwidth_da,
                                    shift_ppm,
                                );
                                lo = lo.min(a);
                                hi = hi.max(b);
                                continue;
                            }
                            let zf = z as f64;
                            let obs = adjusted_observed_neutral_mass(
                                s.precursor_mz * zf - zf * PROTON,
                                shift_ppm,
                            );
                            for o in params.isotope_error_range.clone() {
                                let c = obs - (o as f64) * model::mass::ISOTOPE;
                                lo = lo.min(c - params.precursor_tolerance.left.as_da(c));
                                hi = hi.max(c + params.precursor_tolerance.right.as_da(c));
                            }
                        }
                    }
                    let tables = FormWalkTables::new(params);
                    let recs = match (
                        windows.iter().map(|w| w.0).min(),
                        windows.iter().map(|w| w.1).max(),
                    ) {
                        (Some(wlo), Some(whi)) => {
                            distinct_base_records_for_nominal_window(mi, params, wlo, whi, |v| {
                                retain_records_with_forms_in(self.idx, &tables, v, lo, hi)
                            })
                        }
                        _ => Vec::new(),
                    };
                    let mut index = ChunkFragmentIndex::build(
                        recs,
                        self.idx,
                        params,
                        self.scorer.feature_match_tolerance(),
                        lo,
                        hi,
                    );
                    if let Some(src) = &rt_source {
                        index.enable_rt(src.prior);
                    }
                    index
                });
                (Some(MmapRecordCache::default()), owned_index.as_ref())
            } else {
                let records: FxHashMap<BaseRecordKey, crate::candidates::records::IndexRecord> =
                    windows
                        .par_iter()
                        .fold(FxHashMap::default, |mut m, &(lo, hi)| {
                            for (rec, _) in base_records_for_nominal_window(mi, params, lo, hi) {
                                m.entry(base_record_key(&rec)).or_insert(rec);
                            }
                            m
                        })
                        .reduce(FxHashMap::default, |mut a, b| {
                            a.extend(b);
                            a
                        });
                let budget = crate::params::MMAP_WINDOW_CACHE_MAX_CANDIDATES;
                let cached_total = AtomicUsize::new(0);
                let cache: MmapRecordCache = records
                    .into_par_iter()
                    .filter_map(|(key, rec)| {
                        let v = expand_base_record(self.idx, params, &rec);
                        let before = cached_total.fetch_add(v.len(), Ordering::Relaxed);
                        (before + v.len() <= budget).then_some((key, v))
                    })
                    .collect();
                (Some(cache), None)
            }
        };
        self.chunk_prep_nanos
            .fetch_add(t_chunk_prep.elapsed().as_nanos() as u64, Ordering::Relaxed);

        let chunk = ChunkScoring {
            scorer_set: &scorer_set,
            rt_source: &rt_source,
            record_cache: &record_cache,
            frag_index,
            params,
            spectrum_idx_offset,
            complete_later,
        };
        let queues: Vec<(TopNQueue, Vec<f32>)> = spectra
            .par_iter()
            .enumerate()
            .map(|(local_idx, spec)| self.score_spectrum(&chunk, local_idx, spec))
            .collect();

        queues
    }

    /// Scores one spectrum of a chunk against its candidate windows and fills its top-N queue;
    /// returns the queue and the strong scores kept for the chunk's calibration.
    fn score_spectrum(
        &self,
        chunk: &ChunkScoring<'_>,
        local_idx: usize,
        spec: &Spectrum,
    ) -> (TopNQueue, Vec<f32>) {
        let ChunkScoring {
            scorer_set,
            params,
            spectrum_idx_offset,
            ..
        } = *chunk;
        let fragment_tolerance_da = self.fragment_tolerance_da;
        let aa_set_for_scoring = &self.aa_set_for_scoring;
        let spec_idx = local_idx + spectrum_idx_offset;
        let scorer = scorer_set.for_spec(spec);
        let retention_cap = params.hot_loop_retention_cap();
        let queue = TopNQueue::new(retention_cap);

        if spec.peaks.len() < params.min_peaks as usize {
            return (queue, Vec::new());
        }

        let mut best_two = BestTwo::default();
        let mut ref_tracker = RefTracker::default();
        let track_ref = self.idx.has_alternative_entries;

        let mut tailor_hist: FxHashMap<i32, u32> = FxHashMap::default();
        let mut tailor_total: u32 = 0;
        let mut strong_null_stats = OnlineStats::default();

        let charges_to_try: SmallVec<[u8; 4]> = charges_to_try(spec, params);

        let mut scored_per_charge: SmallVec<[(u8, ScoredSpectrum<'_>); 4]> = SmallVec::new();
        for &z in &charges_to_try {
            if scored_per_charge.iter().all(|(charge, _)| *charge != z) {
                scored_per_charge.push((z, ScoredSpectrum::new(spec, scorer, z)));
            }
        }
        let scored_spec_for_charge = |z: u8| scored_spec(&scored_per_charge, z);

        let ms1_gate = if params.dia {
            self.dia_ms1_scans(spec.rt_seconds)
        } else {
            None
        };
        let (mmap_cands_owned, dia_windows) =
            self.spectrum_candidates(chunk, spec, scorer, &charges_to_try, &ms1_gate);
        let window_cand_indices: Vec<usize> = (0..mmap_cands_owned.len()).collect();
        let mmap_cands: &[Candidate] = &mmap_cands_owned;
        let cand_slice: &[Candidate] = mmap_cands;
        let resolve_cand = |slot: usize| -> (&Candidate, u32) { (&cand_slice[slot], slot as u32) };

        let enz_credit_neighboring = aa_set_for_scoring.neighboring_aa_cleavage_credit();
        let enz_penalty_neighboring = aa_set_for_scoring.neighboring_aa_cleavage_penalty();
        let enz_credit_peptide = aa_set_for_scoring.peptide_cleavage_credit();
        let enz_penalty_peptide = aa_set_for_scoring.peptide_cleavage_penalty();
        let enz_is_c_term = params.enzyme.is_c_term();
        let enz_is_n_term = params.enzyme.is_n_term();
        let enz = params.enzyme;
        let edge_enabled = edge_scoring_enabled(scorer);

        let mut per_charge_queues: FxHashMap<u8, TopNQueue> = FxHashMap::default();
        let open_peaks: Option<Vec<(f64, f32)>> = params
            .open_search
            .map(|_| crate::open::sorted_peaks(&spec.peaks));
        let mut open_attempts: usize = 0;
        let mut open_wins: usize = 0;

        for (pos, &cand_slot) in window_cand_indices.iter().enumerate() {
            prefetch_window_candidates(cand_slice, &window_cand_indices, pos);
            let (cand, cand_idx) = resolve_cand(cand_slot);
            let cleavage_credit = enzyme_cleavage_credit(
                cand,
                enz,
                enz_is_c_term,
                enz_is_n_term,
                enz_credit_neighboring,
                enz_penalty_neighboring,
                enz_credit_peptide,
                enz_penalty_peptide,
            ) as f32;
            let max_edge_bonus_per_edge: f32 = 10.0;
            let n_minus_1 = cand.peptide.length().saturating_sub(1) as f32;
            let max_edge_bonus = max_edge_bonus_per_edge * n_minus_1;
            for &z in &charges_to_try {
                let scored_spec = scored_spec_for_charge(z);
                let iso_errs =
                    self.precursor_errors(spec, cand, z, params, &dia_windows, &ms1_gate);
                if iso_errs.is_empty() {
                    continue;
                }

                let mut pin_score =
                    score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da)
                        + cleavage_credit;

                let mut mod_pep: Option<Peptide> = None;
                let mut open_iso: Option<i8> = None;
                if let (Some(open), Some(peaks)) = (params.open_search, open_peaks.as_deref()) {
                    let zf = f64::from(z);
                    let observed = adjusted_observed_neutral_mass(
                        spec.precursor_mz * zf - zf * PROTON,
                        params.precursor_mass_shift_ppm,
                    );
                    for off in params.isotope_error_range.clone() {
                        let corrected = observed - f64::from(off) * model::mass::ISOTOPE;
                        let delta = corrected - cand.peptide.mass();
                        if delta.abs() <= 0.5 || delta < open.delta_min || delta > open.delta_max {
                            continue;
                        }
                        let loc_tol = fragment_tolerance_da.max(2e-5 * cand.peptide.mass());
                        let Some(loc) = crate::open::localize_and_refine(
                            &cand.peptide,
                            peaks,
                            z,
                            delta,
                            loc_tol,
                            fragment_tolerance_da,
                        ) else {
                            continue;
                        };
                        if loc.best_site == 0 {
                            continue;
                        }
                        open_attempts += 1;
                        let site0 = usize::from(loc.best_site - 1);
                        let shifted = crate::open::shifted_peptide(&cand.peptide, site0, loc.delta);
                        let s2 = score_psm(scored_spec, &shifted, scorer, z, fragment_tolerance_da)
                            + cleavage_credit;
                        if s2 > pin_score {
                            pin_score = s2;
                            mod_pep = Some(shifted);
                            open_iso = Some(off);
                            open_wins += 1;
                        }
                    }
                }
                {
                    let s = pin_score.round() as i32;
                    if track_ref
                        && self.idx.db.proteins[cand.protein_index].entry_kind()
                            == model::EntryKind::Reference
                    {
                        ref_tracker.fold(s, peptidoform_key(&cand.peptide));
                    }
                    *tailor_hist.entry(s).or_insert(0) += 1;
                    tailor_total += 1;
                    strong_null_stats.push(pin_score);
                    best_two.fold(s, || peptidoform_key(&cand.peptide));
                }

                let could_win = match per_charge_queues.get(&z) {
                    Some(q) if q.len() >= q.capacity() as usize => {
                        q.worst_rank_score().is_none_or(|worst| {
                            if edge_enabled {
                                pin_score + max_edge_bonus > worst
                            } else {
                                pin_score >= worst
                            }
                        })
                    }
                    _ => true,
                };
                if !could_win {
                    continue;
                }

                let edge_i = psm_edge_score(
                    scored_spec,
                    mod_pep.as_ref().unwrap_or(&cand.peptide),
                    scorer,
                    z,
                );
                let rank_score = pin_score + edge_i as f32;

                let err = chosen_error(&iso_errs, open_iso, params.open_search.is_some());

                let features = PsmFeatures {
                    open_used_shift: mod_pep.is_some(),
                    ..PsmFeatures::default()
                };
                let psm = PsmMatch {
                    spectrum_idx: spec_idx,
                    candidate_idxs: smallvec::smallvec![cand_idx],
                    charge_used: z,
                    mass_error_ppm: err.mass_error_ppm,
                    score: pin_score,
                    rank_score,
                    edge_score: edge_i,
                    activation_method: Some(scorer.param().data_type.activation),
                    features,
                    isotope_offset: err.isotope_offset,
                    precursor_mz_override: params
                        .dia
                        .then(|| cand.peptide.mass() / z as f64 + PROTON),
                    site_evidence: None,
                    tie_key: crate::psm::tie_key_of(cand),
                };
                per_charge_queues
                    .entry(z)
                    .or_insert_with(|| TopNQueue::new(retention_cap))
                    .push(psm);
            }
        }

        self.finish_spectrum(
            chunk,
            spec,
            SpectrumTail {
                scorer,
                scored_per_charge: &scored_per_charge,
                cand_slice,
                queue,
                per_charge_queues,
                best_raw: best_two.best_raw,
                second_raw: best_two.second_raw,
                ref_tracker,
                track_ref,
                tailor_hist,
                tailor_total,
                strong_null_stats,
                open_peaks,
                open_attempts,
                open_wins,
            },
        )
    }

    /// The precursor mass errors at which candidate `cand` matches the spectrum at charge `z`:
    /// under `--dia` inside the charge's isolation window (and an MS1 peak when the MS1 gate is
    /// on), otherwise at each isotope offset of the range within the precursor tolerance.
    fn precursor_errors(
        &self,
        spec: &Spectrum,
        cand: &Candidate,
        z: u8,
        params: &SearchParams,
        dia_windows: &DiaWindows,
        ms1_gate: &Option<std::ops::Range<usize>>,
    ) -> SmallVec<[MassError; 4]> {
        let mut iso_errs: SmallVec<[MassError; 4]> = SmallVec::new();
        if params.dia {
            let (lo, hi) = dia_windows
                .iter()
                .find(|w| w.0 == z)
                .map_or((f64::MAX, f64::MIN), |w| (w.1, w.2));
            let mass = cand.peptide.mass();
            let mono_mz = mass / z as f64 + PROTON;
            if (lo..=hi).contains(&mass)
                && ms1_gate
                    .clone()
                    .is_none_or(|r| self.dia_ms1_has(r, mono_mz, params.dia_ms1_tol_ppm))
            {
                iso_errs.push(MassError {
                    mass_error_da: 0.0,
                    mass_error_ppm: 0.0,
                    isotope_offset: 0,
                });
            }
        } else {
            for offset in params.isotope_error_range.clone() {
                let matched = matches_precursor(
                    spec,
                    &cand.peptide,
                    z,
                    offset,
                    &params.precursor_tolerance,
                    params.precursor_mass_shift_ppm,
                );
                if let Some(err) = matched {
                    iso_errs.push(err);
                }
            }
        }
        iso_errs
    }

    /// The candidates of one spectrum and its DIA isolation windows per charge: from the
    /// fragment-ion index (gated by the windows and the DIA MS1 scans) or the precursor window of
    /// the out-of-core index, then the candidate filter.
    fn spectrum_candidates(
        &self,
        chunk: &ChunkScoring<'_>,
        spec: &Spectrum,
        scorer: &RankScorer,
        charges_to_try: &[u8],
        ms1_gate: &Option<std::ops::Range<usize>>,
    ) -> (Vec<Candidate>, DiaWindows) {
        let ChunkScoring {
            rt_source,
            record_cache,
            frag_index,
            params,
            ..
        } = *chunk;
        let shift_ppm = params.precursor_mass_shift_ppm;
        let dia_windows: SmallVec<[(u8, f64, f64); 4]> = if params.dia {
            charges_to_try
                .iter()
                .map(|&z| {
                    let (lo, hi) = crate::precursor::matching::dia_neutral_window(
                        spec,
                        z,
                        &params.precursor_tolerance,
                        params.chimeric_isolation_halfwidth_da,
                        shift_ppm,
                    );
                    (z, lo, hi)
                })
                .collect()
        } else {
            SmallVec::new()
        };
        let mut mmap_cands_owned: Vec<Candidate> = if let Some(fi) = &frag_index {
            let gate = |mass: f64| {
                ms1_gate.clone().is_none_or(|r| {
                    dia_windows.iter().any(|&(z, lo, hi)| {
                        (lo..=hi).contains(&mass)
                            && self.dia_ms1_has(
                                r.clone(),
                                mass / z as f64 + PROTON,
                                params.dia_ms1_tol_ppm,
                            )
                    })
                })
            };
            let voted = params
                .dia_vote_filter
                .filter(|_| params.dia)
                .map(|(w, k)| Spectrum {
                    peaks: windowed_top_k(&spec.peaks, w, k),
                    title: String::new(),
                    ..*spec
                });
            let sel = fi.query_where(
                voted.as_ref().unwrap_or(spec),
                charges_to_try,
                params,
                scorer.feature_match_tolerance(),
                params.fragment_index_top_k as usize,
                crate::params::FRAGMENT_INDEX_MIN_MATCHED,
                crate::fragments::index::vote_peaks_for(scorer.feature_match_tolerance()),
                &gate,
                rt_source.as_ref(),
            );
            let relabel =
                |v: &mut [Candidate]| relabel_collision_decoys_with(v, &self.target_bare_seqs);
            fi.materialise(&sel, self.idx, params, record_cache.as_ref(), &relabel)
        } else {
            let key = mmap_window_key(spec, params);
            self.expand_mmap_window_candidates(params, &key, record_cache.as_ref())
        };
        if frag_index.is_some() && self.candidate_filter.is_some() {
            mmap_cands_owned.retain(|c| self.keeps(c));
        }
        (mmap_cands_owned, dia_windows)
    }

    /// Ends scoring one spectrum once its candidates are scored: merges the per-charge queues,
    /// computes the spectrum-level features (delta, Tailor denominator, listwise gap, rank entropy),
    /// fills each retained PSM's features, calibrates the strong score, adds the site evidence and
    /// interns the candidates; returns the queue and the retained strong scores.
    fn finish_spectrum(
        &self,
        chunk: &ChunkScoring<'_>,
        spec: &Spectrum,
        tail: SpectrumTail<'_, '_>,
    ) -> (TopNQueue, Vec<f32>) {
        let SpectrumTail {
            scorer,
            scored_per_charge,
            cand_slice,
            mut queue,
            mut per_charge_queues,
            best_raw,
            second_raw,
            ref_tracker,
            track_ref,
            tailor_hist,
            tailor_total,
            strong_null_stats,
            open_peaks,
            open_attempts,
            open_wins,
        } = tail;
        let ChunkScoring {
            params,
            complete_later,
            ..
        } = *chunk;
        let fragment_tolerance_da = self.fragment_tolerance_da;
        let aa_set_for_scoring = &self.aa_set_for_scoring;
        let enz = params.enzyme;
        let mmap_cands: &[Candidate] = cand_slice;
        let scored_spec_for_charge = |z: u8| scored_spec(scored_per_charge, z);
        for queue in per_charge_queues.values_mut() {
            if queue.len() > 1 {
                let drained = queue.drain_into_vec();
                let mut deduped = dedup_pepseq_score(drained, cand_slice);
                deduped.sort_by(|a, b| b.cmp(a));
                for psm in deduped {
                    queue.push(psm);
                }
            }
        }

        let mut merged: Vec<PsmMatch> = per_charge_queues
            .drain()
            .flat_map(|(_charge, mut per_charge)| per_charge.drain_into_vec())
            .collect();
        merged.sort_by(|a, b| b.cmp(a));
        for psm in merged {
            queue.push(psm);
        }

        let delta_raw_score = if second_raw > i32::MIN {
            (best_raw - second_raw) as f32
        } else {
            0.0
        };

        let tailor_denom = crate::psm::tailor_denominator(&tailor_hist, tailor_total) as f32;

        let mut retained_scores: Vec<f32> = queue.iter_psms().map(|p| p.score).collect();
        retained_scores.sort_by(|a, b| b.total_cmp(a));
        let listwise_score_gap = scoring_crate::psm::listwise_score_gap(&retained_scores);
        let candidate_rank_entropy = scoring_crate::psm::candidate_rank_entropy(&retained_scores);

        let strong_mode = params.score_mode == ScoreMode::Strong;
        let make_features = |psm: &PsmMatch, full: bool| -> PsmFeatures {
            let ss = scored_spec_for_charge(psm.charge_used);
            let cand = &cand_slice[psm.primary_candidate_idx() as usize];
            let used_shift = psm.features.open_used_shift;
            let mut open_loc: Option<crate::open::Localization> = None;
            let mut shifted_pep: Option<Peptide> = None;
            let mut open_delta_ref: Option<f64> = None;
            {
                if let (Some(_), Some(peaks)) = (params.open_search, open_peaks.as_deref()) {
                    let zf = f64::from(psm.charge_used);
                    let mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
                    let observed = adjusted_observed_neutral_mass(
                        mz * zf - zf * PROTON,
                        params.precursor_mass_shift_ppm,
                    ) - f64::from(psm.isotope_offset) * model::mass::ISOTOPE;
                    let delta = observed - cand.peptide.mass();
                    if delta.abs() > 0.5 {
                        let loc_tol = fragment_tolerance_da.max(2e-5 * cand.peptide.mass());
                        if let Some(loc) = crate::open::localize_and_refine(
                            &cand.peptide,
                            peaks,
                            psm.charge_used,
                            delta,
                            loc_tol,
                            fragment_tolerance_da,
                        ) {
                            if loc.best_site != 0 {
                                if used_shift {
                                    shifted_pep = Some(crate::open::shifted_peptide(
                                        &cand.peptide,
                                        usize::from(loc.best_site - 1),
                                        loc.delta,
                                    ));
                                }
                                open_delta_ref = Some(loc.delta);
                            }
                            open_loc = Some(crate::open::Localization {
                                site: loc.site,
                                score: loc.score,
                                margin: loc.margin,
                            });
                        }
                    }
                }
            }
            let pep = shifted_pep.as_ref().unwrap_or(&cand.peptide);
            let mut features = compute_psm_features_with_edge(
                ss,
                pep,
                scorer,
                psm.charge_used,
                Some(psm.edge_score),
                full,
                Some(&self.frag_pred_cache),
            );

            features.delta_raw_score = delta_raw_score;
            if track_ref {
                if let Some(r) = ref_tracker.best_excluding(peptidoform_key(&cand.peptide)) {
                    features.ref_competitor = 1;
                    features.ref_best_raw = r as f32;
                }
            }
            features.listwise_score_gap = listwise_score_gap;
            features.candidate_rank_entropy = candidate_rank_entropy;
            features.tailor_score = psm.score / tailor_denom;
            let credit = cleavage_credit_for(cand, enz, aa_set_for_scoring) as f32;
            features.rank_score_float =
                score_psm_float(ss, pep, scorer, psm.charge_used, fragment_tolerance_da) + credit;
            features.strong_score = fuse_strong_score(&StrongScoreInputs {
                intensity_signal: features.intensity_signal,
                chance_match_surprise: features.chance_match_surprise,
                mass_competition_evidence: features.mass_competition_evidence,
                candidate_rank_entropy: features.candidate_rank_entropy,
                listwise_score_gap: features.listwise_score_gap,
            });
            if let Some(loc) = open_loc {
                features.open_site = loc.site;
                features.open_loc_score = loc.score;
                features.open_loc_margin = loc.margin;
                features.open_delta = open_delta_ref.map_or(0.0, |d| d as f32);
            }
            features.open_used_shift = used_shift;
            features
        };
        queue.fill_post_topn(|psm| {
            psm.features = make_features(psm, !strong_mode);
        });

        let retained_strong: Vec<f32> =
            queue.iter_psms().map(|p| p.features.strong_score).collect();
        queue.fill_post_topn(|psm| {
            psm.features.strong_score_cal = strong_score_calibrated(
                psm.features.strong_score,
                &retained_strong,
                &strong_null_stats,
            );
        });

        if strong_mode {
            queue.reorder_by_strong_score();
            queue.trim_to_capacity(params.top_n_psms_per_spectrum);
        }
        if strong_mode && params.open_search.is_some() {
            queue.update_in_place(|psm| {
                let strong = psm.features.strong_score;
                let strong_cal = psm.features.strong_score_cal;
                psm.features = make_features(psm, true);
                debug_assert_eq!(psm.features.strong_score.to_bits(), strong.to_bits());
                psm.features.strong_score_cal = strong_cal;
            });
        } else if strong_mode && !complete_later {
            queue.update_in_place(|psm| {
                let cand = &cand_slice[psm.primary_candidate_idx() as usize];
                complete_psm_features(
                    &mut psm.features,
                    scored_spec_for_charge(psm.charge_used),
                    &cand.peptide,
                    scorer,
                    psm.charge_used,
                    Some(&self.frag_pred_cache),
                );
            });
        }

        {
            let psms: Vec<&PsmMatch> = queue.iter_psms().collect();
            let rows: Vec<crate::sites::Row<'_>> = psms
                .iter()
                .map(|psm| {
                    let cand = &cand_slice[psm.primary_candidate_idx() as usize];
                    let z = psm.charge_used;
                    let mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
                    crate::sites::Row {
                        peptide: &cand.peptide,
                        is_decoy: cand.is_decoy,
                        charge: z,
                        precursor_mass: (mz - model::mass::PROTON) * z as f64,
                    }
                })
                .collect();
            let key = |p: &PsmMatch| {
                (
                    p.primary_candidate_idx(),
                    p.charge_used,
                    p.rank_score.to_bits(),
                )
            };
            let mut evidence: FxHashMap<_, _> = psms
                .iter()
                .map(|p| key(p))
                .zip(crate::sites::spectrum_evidence(
                    &spec.peaks,
                    &rows,
                    scorer.pin_match_tolerance(),
                ))
                .collect();
            drop(rows);
            drop(psms);
            queue.update_in_place(|psm| {
                psm.site_evidence = evidence.remove(&key(psm)).flatten();
            });
        }

        {
            let mut accum = self.mmap_accum.lock().unwrap();
            queue.remap_candidate_idxs(|local| accum.intern(&mmap_cands[local as usize]));
        }

        if params.open_search.is_some() {
            self.open_rescore_attempts
                .fetch_add(open_attempts, Ordering::Relaxed);
            self.open_rescore_wins
                .fetch_add(open_wins, Ordering::Relaxed);
        }
        (queue, retained_strong)
    }
}

/// The best and second-best raw scores among a spectrum's candidates; the second comes from a
/// peptidoform other than the best's.
struct BestTwo {
    best_raw: i32,
    best_pep_key: u64,
    have_best: bool,
    second_raw: i32,
}

impl Default for BestTwo {
    fn default() -> Self {
        BestTwo {
            best_raw: i32::MIN,
            best_pep_key: 0,
            have_best: false,
            second_raw: i32::MIN,
        }
    }
}

impl BestTwo {
    /// Folds in raw score `s`; `key` (the candidate's peptidoform key) is computed only when `s`
    /// beats the second-best score.
    #[inline]
    fn fold(&mut self, s: i32, key: impl FnOnce() -> u64) {
        if s > self.second_raw {
            let mkey = key();
            if self.have_best && mkey == self.best_pep_key {
                if s > self.best_raw {
                    self.best_raw = s;
                }
            } else if s > self.best_raw {
                self.second_raw = self.second_raw.max(self.best_raw);
                self.best_raw = s;
                self.best_pep_key = mkey;
                self.have_best = true;
            } else {
                self.second_raw = s;
            }
        }
    }
}

/// The precursor mass error a PSM reports: with `--open-search`, the error at the isotope
/// offset its shift was localized at (`open_iso`, else offset 0); otherwise, and as the
/// fallback, the smallest |ppm|.
fn chosen_error(iso_errs: &[MassError], open_iso: Option<i8>, open_search: bool) -> MassError {
    let smaller_ppm = |a: &MassError, b: &MassError| {
        a.mass_error_ppm
            .abs()
            .partial_cmp(&b.mass_error_ppm.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    };
    match open_iso {
        Some(want) => iso_errs
            .iter()
            .copied()
            .find(|e| e.isotope_offset == want)
            .or_else(|| iso_errs.iter().copied().min_by(smaller_ppm)),
        None if open_search => iso_errs
            .iter()
            .copied()
            .find(|e| e.isotope_offset == 0)
            .or_else(|| iso_errs.iter().copied().min_by(smaller_ppm)),
        None => iso_errs.iter().copied().min_by(smaller_ppm),
    }
    .unwrap()
}
