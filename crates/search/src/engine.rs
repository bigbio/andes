//! The search engine: a `PreparedSearch` over the fragment index and the database, run chunk by
//! chunk over the spectra; preparation, setup, DIA MS1, the search loop and candidate lookups are
//! child modules.

mod accumulator;
use accumulator::MmapAccumulator;
pub(crate) mod coisolation;
mod dedup;
mod features;
mod ms1;
mod prepare;
mod queries;
mod run;
mod setup;
mod slices;

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, AtomicUsize};
use std::sync::{Arc, Mutex, RwLock};

use input::Ms1Link;
use model::activation::ActivationMethod;
use model::amino_acid::set::AminoAcidSet;
use model::enzyme::Enzyme;
use model::mass::{nominal_from, H2O, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use scoring_crate::features::sites::{
    mod_site_features, FEAT_SITE_DET_COUNT, FEAT_SITE_INTENS_FRAC, FEAT_SITE_LOCALIZED,
    FEAT_SITE_SHIFTED_FRAC, FEAT_SITE_SHIFTED_MATCHED,
};
use scoring_crate::psm::ions::{predict_by_ions, predict_by_ions_with_losses, IonKind};
use scoring_crate::psm::{
    frag_llr_battery, fragment_ambiguity, intensity_signal, mass_competition_evidence,
    predict_frag_intensities, psm_edge_score, rich_ion_llr, score_psm, strong_score_calibrated_loo,
    theoretical_fragment_ions, FragmentIon, MatchedIon, RankScorer, ScoredSpectrum, DENSITY_HW,
};
use smallvec::{smallvec, SmallVec};

use crate::candidates::records::MmapCandidateIndex;
use crate::candidates::{expand_base_record, BaseRecordKey, Candidate};
use crate::fragments::index::ChunkFragmentIndex;
use crate::index::SearchIndex;
use crate::params::SearchParams;
use crate::precursor::calibration::adjusted_observed_neutral_mass;
use crate::psm::{PsmFeatures, PsmMatch, TopNQueue};
pub use coisolation::*;
pub(crate) use dedup::*;
pub(crate) use features::*;
pub use slices::*;

pub type CandidateFilter = Box<dyn Fn(&Candidate) -> bool + Send + Sync>;

const DIA_MS1_GATE_SCANS: usize = 1;
/// Spectra whose kept PSMs get their full-pass-only features per lock of the candidate store.
const KEPT_FEATURE_BATCH_SPECTRA: usize = 2048;

fn windowed_top_k(peaks: &[(f64, f32)], width: f64, k: usize) -> Vec<(f64, f32)> {
    let mut by_window: Vec<(i64, f32, usize)> = peaks
        .iter()
        .enumerate()
        .map(|(i, p)| ((p.0 / width).floor() as i64, p.1, i))
        .collect();
    by_window.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.total_cmp(&a.1)));
    let mut keep: Vec<usize> = Vec::with_capacity(peaks.len().min(k * 32));
    let (mut win, mut n) = (i64::MIN, 0usize);
    for &(w, _, i) in &by_window {
        if w != win {
            (win, n) = (w, 0);
        }
        if n < k {
            keep.push(i);
            n += 1;
        }
    }
    keep.sort_unstable();
    keep.into_iter().map(|i| peaks[i]).collect()
}

pub type ScorerLoader = Box<dyn Fn(ActivationMethod) -> Option<RankScorer> + Send + Sync>;

struct PerScanScorers {
    run_activation: ActivationMethod,
    loader: ScorerLoader,
    loaded: RwLock<Vec<ActivationRoute>>,
}

struct ActivationRoute {
    activation: ActivationMethod,
    scorer: Option<Arc<RankScorer>>,
    scans: AtomicU64,
}

pub(crate) struct ScorerSet<'a> {
    primary: &'a RankScorer,
    run_activation: Option<ActivationMethod>,
    alts: Vec<(ActivationMethod, Arc<RankScorer>)>,
}

impl ScorerSet<'_> {
    #[inline]
    pub(crate) fn for_spec(&self, spec: &Spectrum) -> &RankScorer {
        match spec.activation_method {
            Some(a) if Some(a) != self.run_activation => self
                .alts
                .iter()
                .find(|(m, _)| *m == a)
                .map_or(self.primary, |(_, s)| &**s),
            _ => self.primary,
        }
    }
}

pub struct IpfScoring<'a> {
    scorers: ScorerSet<'a>,
    fragment_tolerance_da: f64,
    declared: Vec<&'static model::modification::Modification>,
}

impl IpfScoring<'_> {
    pub fn set_fragment_error(&self, m: Option<scoring_crate::psm::FragmentErrorModel>) {
        self.scorers.primary.set_fragment_error(m);
        for (_, s) in &self.scorers.alts {
            s.set_fragment_error(m);
        }
    }

    pub fn fragment_error(&self) -> Option<scoring_crate::psm::FragmentErrorModel> {
        self.scorers.primary.fragment_error()
    }

    pub fn evidence(
        &self,
        spec: &Spectrum,
        rows: &[(&Peptide, u8)],
    ) -> Vec<Option<crate::ipf::IpfEvidence>> {
        let shifted = recentred_spectra(
            std::slice::from_ref(spec),
            self.scorers.primary.fragment_recentre_ppm(),
        );
        let spec = shifted.as_ref().map_or(spec, |v| &v[0]);
        let scorer = self.scorers.for_spec(spec);
        let tol = scorer.pin_match_tolerance();
        let tiers = crate::ipf::peak_tiers(&spec.peaks);
        let mut scored: SmallVec<[(u8, ScoredSpectrum<'_>); 4]> = SmallVec::new();
        rows.iter()
            .map(|&(peptide, z)| {
                if !scored.iter().any(|(c, _)| *c == z) {
                    scored.push((z, ScoredSpectrum::new(spec, scorer, z)));
                }
                let ss = &scored.iter().find(|(c, _)| *c == z).expect("built above").1;
                let score = |p: &Peptide| {
                    score_psm(ss, p, scorer, z, self.fragment_tolerance_da)
                        + psm_edge_score(ss, p, scorer, z) as f32
                };
                crate::ipf::psm_evidence(
                    &spec.peaks,
                    &tiers,
                    peptide,
                    z,
                    move |mz| tol.as_da(mz),
                    &self.declared,
                    Some(&score),
                )
            })
            .collect()
    }
}

pub(crate) fn recentred_spectra(spectra: &[Spectrum], offset_ppm: f64) -> Option<Vec<Spectrum>> {
    if offset_ppm == 0.0 || !offset_ppm.is_finite() {
        return None;
    }
    let scale = 1.0 / (1.0 + offset_ppm * 1e-6);
    Some(
        spectra
            .iter()
            .map(|s| {
                let mut c = s.clone();
                for p in &mut c.peaks {
                    p.0 *= scale;
                }
                c
            })
            .collect(),
    )
}

pub struct PreparedSearch<'a> {
    pub idx: &'a SearchIndex,
    pub params: &'a SearchParams,
    pub scorer: &'a RankScorer,
    pub fragment_tolerance_da: f64,
    pub candidates: Vec<Candidate>,
    pub mmap_index: MmapCandidateIndex,
    mmap_accum: Mutex<MmapAccumulator>,
    chunk_prep_nanos: AtomicU64,
    open_rescore_attempts: AtomicUsize,
    open_rescore_wins: AtomicUsize,
    pub aa_set_for_scoring: AminoAcidSet,
    pub ms1_link: Option<Ms1Link>,
    dia_ms1: Vec<model::scan::Ms1Scan>,
    dia_ms1_bins: Vec<ms1::Ms1MzBins>,
    dia_scan_mobility: input::DiaScanMobility,
    dia_mobility_margin: RwLock<Option<f32>>,
    target_bare_seqs: FoldedSet,
    candidate_filter: Option<CandidateFilter>,
    open_index: Option<Box<ChunkFragmentIndex>>,
    per_scan: Option<PerScanScorers>,
    form_rt_model: Option<Arc<crate::dia::FormRtModel>>,
    dia_rt_prior: RwLock<Option<Arc<crate::dia::DiaRtPrior>>>,
    form_rt_cache: crate::dia::FormRtCache,
    frag_pred_cache: crate::fragments::predictions::FragPredCache,
}

fn enumeration_pass(cand: &Candidate, copy: u32) -> u32 {
    if cand.is_protein_n_term && cand.start_offset_in_protein == 1 {
        1
    } else {
        copy
    }
}

pub(crate) fn build_mass_order(candidates: &[Candidate]) -> Vec<u32> {
    let mut order: Vec<u32> = (0..candidates.len() as u32).collect();
    order.par_sort_unstable_by(|&a, &b| {
        candidates[a as usize]
            .peptide
            .mass()
            .total_cmp(&candidates[b as usize].peptide.mass())
            .then(a.cmp(&b))
    });
    order
}

const PREFETCH_STRUCT_AHEAD: usize = 16;
const PREFETCH_RESIDUES_AHEAD: usize = 8;

const PREFETCH_RESIDUE_LINES: usize = 8;

#[inline(always)]
fn prefetch_window_candidates(cands: &[Candidate], window: &[usize], pos: usize) {
    #[cfg(target_arch = "x86_64")]
    {
        use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
        if let Some(c) = window
            .get(pos + PREFETCH_STRUCT_AHEAD)
            .and_then(|&slot| cands.get(slot))
        {
            let p = (c as *const Candidate).cast::<i8>();
            // SAFETY: prefetch never dereferences; `p` and the struct's last byte
            #[allow(unused_unsafe)]
            unsafe {
                _mm_prefetch(p, _MM_HINT_T0);
                _mm_prefetch(p.add(std::mem::size_of::<Candidate>() - 1), _MM_HINT_T0);
            }
        }
        if let Some(c) = window
            .get(pos + PREFETCH_RESIDUES_AHEAD)
            .and_then(|&slot| cands.get(slot))
        {
            let residues = c.peptide.residues.as_slice();
            let bytes = std::mem::size_of_val(residues);
            let base = residues.as_ptr().cast::<i8>();
            for line in 0..PREFETCH_RESIDUE_LINES.min(bytes.div_ceil(64)) {
                // SAFETY: prefetch never dereferences; every address is inside
                #[allow(unused_unsafe)]
                unsafe {
                    _mm_prefetch(base.add(line * 64), _MM_HINT_T0);
                }
            }
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = (cands, window, pos);
}

pub(crate) fn candidate_nominal_bounds(
    spec: &Spectrum,
    z: u8,
    params: &SearchParams,
    shift_ppm: f64,
) -> (i32, i32) {
    if params.dia {
        let (lo, hi) = crate::precursor::matching::dia_neutral_window(
            spec,
            z,
            &params.precursor_tolerance,
            params.chimeric_isolation_halfwidth_da,
            shift_ppm,
        );
        return (nominal_from(lo - H2O) - 1, nominal_from(hi - H2O) + 1);
    }
    let charge_f = z as f64;
    let iso_min = *params.isotope_error_range.start() as i32;
    let iso_max = *params.isotope_error_range.end() as i32;

    let neutral_mass =
        adjusted_observed_neutral_mass((spec.precursor_mz - PROTON) * charge_f - H2O, shift_ppm);
    let nominal_center = nominal_from(neutral_mass);
    let tol_da_left = params.precursor_tolerance.left.as_da(neutral_mass);
    let tol_da_right = params.precursor_tolerance.right.as_da(neutral_mass);
    let widen_left = (tol_da_left - 0.4999_f64).round() as i32;
    let widen_right = (tol_da_right - 0.4999_f64).round() as i32;
    let min_nominal = nominal_center - iso_max - widen_left;
    let max_nominal = nominal_center - iso_min + widen_right;
    let observed =
        adjusted_observed_neutral_mass(spec.precursor_mz * charge_f - charge_f * PROTON, shift_ppm);
    let (mut lightest, mut heaviest) = (f64::MAX, f64::MIN);
    for offset in params.isotope_error_range.clone() {
        let target = observed - offset as f64 * model::mass::ISOTOPE;
        lightest = lightest.min(target - params.precursor_tolerance.left.as_da(target));
        heaviest = heaviest.max(target + params.precursor_tolerance.right.as_da(target));
    }
    (
        min_nominal.min(nominal_from(lightest - H2O)),
        max_nominal.max(nominal_from(heaviest - H2O)),
    )
}

fn charges_to_try(spec: &Spectrum, params: &SearchParams) -> SmallVec<[u8; 4]> {
    match spec.precursor_charge.and_then(|z| u8::try_from(z).ok()) {
        Some(z) if z > 0 => smallvec![z],
        _ => params.charge_range.clone().collect(),
    }
}

fn mmap_window_key(spec: &Spectrum, params: &SearchParams) -> MmapWindowKey {
    let shift_ppm = params.precursor_mass_shift_ppm;
    charges_to_try(spec, params)
        .into_iter()
        .map(|z| candidate_nominal_bounds(spec, z, params, shift_ppm))
        .collect()
}

type MmapWindowKey = SmallVec<[(i32, i32); 4]>;

type MmapRecordCache = FxHashMap<BaseRecordKey, Vec<Candidate>>;

/// Re-calibrates `strong_score_cal` of every PSM of one spectrum against `pool`, the strong
/// scores retained over every charge the spectrum was searched with. With fewer than two pooled
/// scores the per-charge value is kept.
fn recalibrate_strong_over_spectrum(queue: &mut TopNQueue, pool: &[f32]) {
    if pool.len() < 2 {
        return;
    }
    queue.update_in_place(|psm| {
        psm.features.strong_score_cal =
            strong_score_calibrated_loo(pool, psm.features.strong_score);
    });
}

#[cfg(test)]
#[path = "engine/tests/accumulator.rs"]
mod accumulator_tests;

#[cfg(test)]
#[path = "engine/tests/per_scan_scorer.rs"]
mod per_scan_scorer_tests;

#[cfg(test)]
#[path = "engine/tests/chargeless.rs"]
mod chargeless_unit_tests;

#[cfg(test)]
#[path = "engine/tests/charges_to_try.rs"]
mod charges_to_try_tests;

#[cfg(test)]
#[path = "engine/tests/feature.rs"]
mod feature_tests;

#[cfg(test)]
#[path = "engine/tests/dedup.rs"]
mod dedup_tests;

#[derive(Debug, Default, Clone, Copy)]
struct RefTracker {
    best: Option<(i32, u64)>,
    second: Option<i32>,
}

impl RefTracker {
    fn fold(&mut self, score: i32, key: u64) {
        match self.best {
            None => self.best = Some((score, key)),
            Some((best, best_key)) if key == best_key => {
                self.best = Some((best.max(score), key));
            }
            Some((best, _)) if score > best => {
                self.second = Some(self.second.map_or(best, |s| s.max(best)));
                self.best = Some((score, key));
            }
            Some(_) => self.second = Some(self.second.map_or(score, |s| s.max(score))),
        }
    }

    fn best_excluding(&self, key: u64) -> Option<i32> {
        match self.best {
            Some((_, best_key)) if best_key == key => self.second,
            Some((best, _)) => Some(best),
            None => None,
        }
    }
}

#[cfg(test)]
#[path = "engine/tests/ref_tracker.rs"]
mod ref_tracker_tests;
