use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, RwLock};

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::{smallvec, SmallVec};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use crate::candidate_gen::{
    base_record_key, base_records_for_nominal_window, distinct_base_records_for_nominal_window,
    expand_base_record, retain_records_with_forms_in, BaseRecordKey, Candidate, FormWalkTables,
};
use crate::candidate_index::MmapCandidateIndex;
use crate::fragment_index::ChunkFragmentIndex;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::precursor_matching::{matches_precursor, MassError};
use crate::psm::{PsmFeatures, PsmMatch, TopNQueue};
use crate::search_index::SearchIndex;
use crate::search_params::{ScoreMode, SearchParams};
use input::Ms1Link;
use model::aa_set::AminoAcidSet;
use model::activation::ActivationMethod;
use model::enzyme::Enzyme;
use model::mass::{nominal_from, H2O, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::mod_site_features::{
    mod_site_features, FEAT_SITE_DET_COUNT, FEAT_SITE_INTENS_FRAC, FEAT_SITE_LOCALIZED,
    FEAT_SITE_SHIFTED_FRAC, FEAT_SITE_SHIFTED_MATCHED,
};
use scoring_crate::scoring::fragment_ions::{
    predict_by_ions, predict_by_ions_with_losses, IonKind,
};
use scoring_crate::scoring::{
    edge_scoring_enabled, frag_llr_battery, fragment_ambiguity, fuse_strong_score,
    intensity_signal, mass_competition_evidence, predict_frag_intensities, psm_edge_score,
    rich_ion_llr, score_psm, score_psm_float, scored_for_charge, strong_score_calibrated,
    strong_score_calibrated_loo, theoretical_fragment_ions, FragmentIon, MatchedIon, OnlineStats,
    RankScorer, ScoredSpectrum, StrongScoreInputs, DENSITY_HW,
};

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
    pub fn set_fragment_error(&self, m: Option<scoring_crate::scoring::FragmentErrorModel>) {
        self.scorers.primary.set_fragment_error(m);
        for (_, s) in &self.scorers.alts {
            s.set_fragment_error(m);
        }
    }

    pub fn fragment_error(&self) -> Option<scoring_crate::scoring::FragmentErrorModel> {
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
    pub dia_ms1: Vec<model::scan::Ms1Scan>,
    target_bare_seqs: FxHashSet<Box<[u8]>>,
    candidate_filter: Option<CandidateFilter>,
    open_index: Option<Box<ChunkFragmentIndex>>,
    per_scan: Option<PerScanScorers>,
    form_rt_model: Option<Arc<crate::dia_rt_rank::FormRtModel>>,
    dia_rt_prior: RwLock<Option<Arc<crate::dia_rt_rank::DiaRtPrior>>>,
    form_rt_cache: crate::dia_rt_rank::FormRtCache,
    frag_pred_cache: crate::frag_pred_cache::FragPredCache,
}

#[derive(Clone, PartialEq, Eq, Hash)]
struct GlobalCandKey {
    protein_index: u32,
    start_offset: u32,
    residues: Box<[u8]>,
    mod_units: Box<[i32]>,
    is_decoy: bool,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
}

impl GlobalCandKey {
    fn from_candidate(c: &Candidate) -> Self {
        let mut residues = Vec::with_capacity(c.peptide.residues.len());
        let mut mod_units = Vec::with_capacity(c.peptide.residues.len());
        for aa in &c.peptide.residues {
            residues.push(aa.residue);
            mod_units.push(
                aa.mod_
                    .as_ref()
                    .map(|m| (m.mass_delta * 100_000.0).round() as i32)
                    .unwrap_or(0),
            );
        }
        Self {
            protein_index: c.protein_index as u32,
            start_offset: c.start_offset_in_protein as u32,
            residues: residues.into_boxed_slice(),
            mod_units: mod_units.into_boxed_slice(),
            is_decoy: c.is_decoy,
            is_protein_n_term: c.is_protein_n_term,
            is_protein_c_term: c.is_protein_c_term,
        }
    }
}

#[derive(Default)]
struct MmapAccumulator {
    candidates: Vec<Candidate>,
    index_of: FxHashMap<GlobalCandKey, u32>,
}

impl MmapAccumulator {
    fn intern(&mut self, cand: &Candidate) -> u32 {
        let key = GlobalCandKey::from_candidate(cand);
        if let Some(&idx) = self.index_of.get(&key) {
            return idx;
        }
        let idx = self.candidates.len() as u32;
        self.candidates.push(cand.clone());
        self.index_of.insert(key, idx);
        idx
    }
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
        let (lo, hi) = crate::precursor_matching::dia_neutral_window(
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

impl<'a> PreparedSearch<'a> {
    pub fn prepare_mmap(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        decoy_prefix: &str,
        index_path: &std::path::Path,
    ) -> std::io::Result<Self> {
        let (mmap_index, was_built) =
            MmapCandidateIndex::open_or_build(index_path, idx, params, decoy_prefix)?;
        if mmap_index.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "the digest built no candidate peptides (enzyme {}, extra enzymes {}, \
                     length {}..{}, max missed cleavages {}, termini {}); check --enzyme, \
                     --min-length, --max-length, --max-missed-cleavages and the database",
                    params.enzyme.name(),
                    params.extra_enzymes.len(),
                    params.min_length,
                    params.max_length,
                    params.max_missed_cleavages,
                    params.num_tolerable_termini,
                ),
            ));
        }
        if was_built {
            log::info!(
                target: "candidate-index",
                "built and cached {} records → {}",
                mmap_index.len(),
                index_path.display()
            );
        } else {
            log::info!(
                target: "candidate-index",
                "reusing cached index ({} records) from {}",
                mmap_index.len(),
                index_path.display()
            );
        }
        let mut prepared = Self::over_index(idx, params, scorer, fragment_tolerance_da, mmap_index);
        prepared.open_index = Self::build_open_index(idx, params, scorer, &prepared.mmap_index);
        Ok(prepared)
    }

    pub fn prepare_in_memory(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        decoy_prefix: &str,
    ) -> std::io::Result<Self> {
        let mmap_index = MmapCandidateIndex::build_in_memory(idx, params, decoy_prefix)?;
        let mut prepared = Self::over_index(idx, params, scorer, fragment_tolerance_da, mmap_index);
        prepared.open_index = Self::build_open_index(idx, params, scorer, &prepared.mmap_index);
        Ok(prepared)
    }

    fn build_open_index(
        idx: &SearchIndex,
        params: &SearchParams,
        scorer: &RankScorer,
        mmap_index: &MmapCandidateIndex,
    ) -> Option<Box<ChunkFragmentIndex>> {
        params.open_search?;
        let t = std::time::Instant::now();
        let records = mmap_index.records();
        let (lo, hi) = match (records.first(), records.last()) {
            (Some(a), Some(b)) => (a.mass_milli as f64 / 1000.0, b.mass_milli as f64 / 1000.0),
            _ => return None,
        };
        let max_mod = params
            .aa_set
            .distinct_mods()
            .iter()
            .filter(|m| !m.fixed)
            .map(|m| m.mass_delta.max(0.0))
            .fold(0.0_f64, f64::max);
        let max_sites = params.max_variable_mods_per_peptide + params.max_label_sites;
        let hi = hi + max_mod * f64::from(max_sites) + 1.0;
        let index = ChunkFragmentIndex::build(
            records,
            idx,
            params,
            scorer.feature_match_tolerance(),
            lo,
            hi,
        );
        log::info!(
            target: "open-search-index",
            "{:.1}-{:.1} Da, {} forms, built in {:.1} s",
            lo,
            hi,
            index.n_forms(),
            t.elapsed().as_secs_f64()
        );
        Some(Box::new(index))
    }

    fn over_index(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        mmap_index: MmapCandidateIndex,
    ) -> Self {
        let mut aa_set_for_scoring: AminoAcidSet = params.aa_set.clone();
        if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
            aa_set_for_scoring.register_enzyme(params.enzyme, 0.95, 0.95);
        }
        let target_bare_seqs = build_target_bare_seqs(&mmap_index, idx);
        PreparedSearch {
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            candidates: Vec::new(),
            mmap_index,
            mmap_accum: Mutex::new(MmapAccumulator::default()),
            chunk_prep_nanos: AtomicU64::new(0),
            open_rescore_attempts: AtomicUsize::new(0),
            open_rescore_wins: AtomicUsize::new(0),
            aa_set_for_scoring,
            ms1_link: None,
            dia_ms1: Vec::new(),
            target_bare_seqs,
            candidate_filter: None,
            open_index: None,
            per_scan: None,
            form_rt_model: None,
            dia_rt_prior: RwLock::new(None),
            form_rt_cache: crate::dia_rt_rank::FormRtCache::default(),
            frag_pred_cache: crate::frag_pred_cache::FragPredCache::default(),
        }
    }

    /// Installs the per-form RT predictor of the DIA vote cut.
    pub fn set_form_rt_model(&mut self, model: Arc<crate::dia_rt_rank::FormRtModel>) {
        self.form_rt_model = Some(model);
    }

    pub fn form_rt_model(&self) -> Option<&crate::dia_rt_rank::FormRtModel> {
        self.form_rt_model.as_deref()
    }

    /// Sets (or with `None` clears) the run's RT prior; while set, every DIA chunk index
    /// carries calibrated per-form RTs and its vote cut uses the RT term.
    /// Clearing also drops the cached per-form RT predictions.
    pub fn set_dia_rt_prior(&self, prior: Option<crate::dia_rt_rank::DiaRtPrior>) {
        if prior.is_none() {
            self.form_rt_cache.clear();
        }
        *self.dia_rt_prior.write().expect("RT prior lock") = prior.map(Arc::new);
    }

    pub fn with_per_scan_scorers(
        mut self,
        run_activation: ActivationMethod,
        loader: ScorerLoader,
    ) -> Self {
        self.per_scan = Some(PerScanScorers {
            run_activation,
            loader,
            loaded: RwLock::new(Vec::new()),
        });
        self
    }

    pub fn set_fragment_error(&self, m: Option<scoring_crate::scoring::FragmentErrorModel>) {
        self.scorer.set_fragment_error(m);
        if let Some(ps) = &self.per_scan {
            for route in ps.loaded.read().unwrap().iter() {
                if let Some(s) = &route.scorer {
                    s.set_fragment_error(m);
                }
            }
        }
    }

    pub fn per_scan_routes(&self) -> Vec<(ActivationMethod, bool, u64)> {
        self.per_scan.as_ref().map_or_else(Vec::new, |ps| {
            ps.loaded
                .read()
                .unwrap()
                .iter()
                .map(|r| {
                    (
                        r.activation,
                        r.scorer.is_some(),
                        r.scans.load(Ordering::Relaxed),
                    )
                })
                .collect()
        })
    }

    pub(crate) fn scorer_set(&self, spectra: &[Spectrum]) -> ScorerSet<'a> {
        let Some(ps) = &self.per_scan else {
            return ScorerSet {
                primary: self.scorer,
                run_activation: None,
                alts: Vec::new(),
            };
        };
        let mut seen: SmallVec<[(ActivationMethod, u64); 4]> = SmallVec::new();
        for a in spectra.iter().filter_map(|s| s.activation_method) {
            if a == ps.run_activation {
                continue;
            }
            match seen.iter_mut().find(|(m, _)| *m == a) {
                Some((_, n)) => *n += 1,
                None => seen.push((a, 1)),
            }
        }
        for &(a, n) in &seen {
            let known = ps.loaded.read().unwrap().iter().any(|r| r.activation == a);
            if !known {
                let mut loaded = ps.loaded.write().unwrap();
                if !loaded.iter().any(|r| r.activation == a) {
                    let scorer = (ps.loader)(a).map(|s| {
                        s.set_fragment_error(self.scorer.fragment_error());
                        Arc::new(s)
                    });
                    loaded.push(ActivationRoute {
                        activation: a,
                        scorer,
                        scans: AtomicU64::new(0),
                    });
                }
            }
            if let Some(r) = ps.loaded.read().unwrap().iter().find(|r| r.activation == a) {
                r.scans.fetch_add(n, Ordering::Relaxed);
            }
        }
        let alts = ps
            .loaded
            .read()
            .unwrap()
            .iter()
            .filter_map(|r| r.scorer.as_ref().map(|s| (r.activation, Arc::clone(s))))
            .collect();
        ScorerSet {
            primary: self.scorer,
            run_activation: Some(ps.run_activation),
            alts,
        }
    }

    pub fn ipf_scoring<'s>(&self, scorer: &'s RankScorer) -> IpfScoring<'s> {
        assert!(
            std::ptr::eq(scorer, self.scorer),
            "ipf_scoring needs the scorer the search was prepared with"
        );
        let set = self.scorer_set(&[]);
        IpfScoring {
            scorers: ScorerSet {
                primary: scorer,
                run_activation: set.run_activation,
                alts: set.alts,
            },
            fragment_tolerance_da: self.fragment_tolerance_da,
            declared: self.params.aa_set.distinct_mods(),
        }
    }

    pub fn chunk_prep_seconds(&self) -> f64 {
        self.chunk_prep_nanos.load(Ordering::Relaxed) as f64 * 1e-9
    }

    pub fn open_search_stats(&self) -> (usize, usize) {
        (
            self.open_rescore_attempts.load(Ordering::Relaxed),
            self.open_rescore_wins.load(Ordering::Relaxed),
        )
    }
    pub fn sync_materialized_candidates(&mut self) {
        let mut accum = self.mmap_accum.lock().unwrap();
        let materialized = std::mem::take(&mut accum.candidates);
        accum.index_of = Default::default();
        drop(accum);
        self.candidates = materialized;
    }

    pub fn with_candidate_filter(mut self, keep: CandidateFilter) -> Self {
        self.target_bare_seqs =
            build_target_bare_seqs_where(&self.mmap_index, self.idx, self.params, &keep);
        self.candidate_filter = Some(keep);
        self
    }

    fn keeps(&self, cand: &Candidate) -> bool {
        self.candidate_filter.as_ref().is_none_or(|keep| keep(cand))
    }

    pub fn with_ms1_link(mut self, ms1_link: Option<Ms1Link>) -> Self {
        self.ms1_link = ms1_link;
        self
    }

    pub fn set_dia_ms1(&mut self, mut scans: Vec<model::scan::Ms1Scan>) {
        scans.retain(|s| s.rt.is_finite());
        scans.sort_by(|a, b| a.rt.total_cmp(&b.rt));
        self.dia_ms1 = scans;
    }

    pub fn take_dia_ms1(&mut self) -> Vec<model::scan::Ms1Scan> {
        std::mem::take(&mut self.dia_ms1)
    }

    pub fn shift_dia_ms1_ppm(&mut self, ppm: f64) {
        let f = 1.0 - ppm * 1e-6;
        for s in &mut self.dia_ms1 {
            for p in &mut s.peaks {
                p.0 *= f;
            }
        }
    }

    pub fn dia_ms1_error_ppm(&self, rt: Option<f64>, mz: f64, window_ppm: f64) -> Option<f64> {
        let tol = mz * window_ppm * 1e-6;
        self.dia_ms1_scans(rt)?
            .filter_map(|k| {
                let s = &self.dia_ms1[k];
                let i = s.peaks.partition_point(|p| p.0 < mz - tol);
                s.peaks[i..]
                    .iter()
                    .take_while(|p| p.0 <= mz + tol)
                    .max_by(|a, b| a.1.total_cmp(&b.1))
                    .copied()
            })
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|p| (p.0 - mz) / mz * 1e6)
    }

    fn dia_ms1_has(&self, scans: std::ops::Range<usize>, mz: f64, tol_ppm: f64) -> bool {
        let tol = mz * tol_ppm * 1e-6;
        self.dia_ms1[scans].iter().any(|s| {
            let i = s.peaks.partition_point(|p| p.0 < mz - tol);
            s.peaks.get(i).is_some_and(|p| p.0 <= mz + tol)
        })
    }

    fn dia_ms1_scans(&self, rt: Option<f64>) -> Option<std::ops::Range<usize>> {
        let rt = rt?;
        if self.dia_ms1.is_empty() {
            return None;
        }
        let i = self.dia_ms1.partition_point(|s| s.rt < rt);
        let near = match (i.checked_sub(1), (i < self.dia_ms1.len()).then_some(i)) {
            (Some(a), Some(b)) => {
                if rt - self.dia_ms1[a].rt <= self.dia_ms1[b].rt - rt {
                    a
                } else {
                    b
                }
            }
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => return None,
        };
        Some(
            near.saturating_sub(DIA_MS1_GATE_SCANS)
                ..(near + DIA_MS1_GATE_SCANS + 1).min(self.dia_ms1.len()),
        )
    }

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
        let mass = crate::fragment_index::index_order_mass;
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
            let take = crate::fragment_index::index_chunk_len(&units, |u| mass(&u.1), span_da);
            let (owners, chunk): (Vec<usize>, Vec<Spectrum>) = units.drain(..take).unzip();
            for ((mut q, pool), i) in self
                .run_chunk_inner(&chunk, 0, params, complete_later)
                .into_iter()
                .zip(owners)
            {
                pools[i].extend(pool);
                for mut psm in q.drain_into_vec() {
                    psm.spectrum_idx = spectrum_idx_offset + i;
                    merged[i].push(psm);
                }
            }
        }
        for (queue, pool) in merged.iter_mut().zip(&pools) {
            recalibrate_strong_over_spectrum(queue, pool);
        }
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
                let accum = self.mmap_accum.lock().unwrap();
                batch
                    .iter()
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

    fn expand_mmap_window_candidates(
        &self,
        params: &SearchParams,
        windows: &[(i32, i32)],
        cache: Option<&MmapRecordCache>,
    ) -> Vec<Candidate> {
        let mi = &self.mmap_index;

        let mut records: FxHashMap<BaseRecordKey, (crate::candidate_index::IndexRecord, u32)> =
            FxHashMap::default();
        for &(min_nominal, max_nominal) in windows {
            for (rec, mult) in base_records_for_nominal_window(mi, params, min_nominal, max_nominal)
            {
                records.entry(base_record_key(&rec)).or_insert((rec, mult));
            }
        }
        let mut records: Vec<(crate::candidate_index::IndexRecord, u32)> =
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

    pub fn out_of_core_candidates_where(
        &self,
        keep: impl Fn(&[u8], u8) -> bool + Sync,
    ) -> Vec<Candidate> {
        let db = self.idx;
        let params = self.params;
        let records = self.out_of_core_records_where(keep);
        let mut tagged: Vec<(u32, u32, u32, u16, Candidate)> = records
            .par_iter()
            .flat_map_iter(|(rec, mult)| {
                let mult = *mult;
                expand_base_record(db, params, rec)
                    .into_iter()
                    .flat_map(move |cand| {
                        (0..mult).map(move |k| {
                            (
                                cand.protein_index as u32,
                                enumeration_pass(&cand, k),
                                cand.start_offset_in_protein as u32,
                                cand.peptide.residues.len() as u16,
                                cand.clone(),
                            )
                        })
                    })
            })
            .collect();
        tagged.par_sort_by(|a, b| (a.0, a.1, a.2, a.3).cmp(&(b.0, b.1, b.2, b.3)));
        let mut cands: Vec<Candidate> = tagged.into_par_iter().map(|t| t.4).collect();
        relabel_collision_decoys_with(&mut cands, &self.target_bare_seqs);
        if self.candidate_filter.is_some() {
            cands.retain(|c| self.keeps(c));
        }
        cands
    }

    pub fn out_of_core_records_where(
        &self,
        keep: impl Fn(&[u8], u8) -> bool + Sync,
    ) -> Vec<(crate::candidate_index::IndexRecord, u32)> {
        let mi = &self.mmap_index;
        let db = self.idx;
        let mut kept: Vec<crate::candidate_index::IndexRecord> = (0..mi.len())
            .into_par_iter()
            .filter_map(|i| {
                let rec = mi.record_at(i);
                let seq = &db.db.proteins.get(rec.protein_index as usize)?.sequence;
                let start = rec.start_offset as usize;
                let end = start + rec.length as usize;
                let span = seq.get(start..end)?;
                let post = seq.get(end).copied().unwrap_or(b'-');
                keep(span, post).then_some(rec)
            })
            .collect();
        kept.par_sort_unstable_by_key(base_record_key);
        let mut records: Vec<(crate::candidate_index::IndexRecord, u32)> = Vec::new();
        for rec in kept {
            match records.last_mut() {
                Some((last, n)) if base_record_key(last) == base_record_key(&rec) => *n += 1,
                _ => records.push((rec, 1)),
            }
        }
        records
    }

    pub fn candidate_peptide(&self, global_idx: u32) -> Peptide {
        let accum = self.mmap_accum.lock().unwrap();
        match accum.candidates.get(global_idx as usize) {
            Some(c) => c.peptide.clone(),
            None => self.candidates[global_idx as usize].peptide.clone(),
        }
    }

    pub fn candidate_is_decoy(&self, global_idx: u32) -> bool {
        let accum = self.mmap_accum.lock().unwrap();
        match accum.candidates.get(global_idx as usize) {
            Some(c) => c.is_decoy,
            None => self.candidates[global_idx as usize].is_decoy,
        }
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
        let fragment_tolerance_da = self.fragment_tolerance_da;
        let aa_set_for_scoring = &self.aa_set_for_scoring;

        let t_chunk_prep = std::time::Instant::now();
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
                                let (a, b) = crate::precursor_matching::dia_neutral_window(
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
                    let prior = self.dia_rt_prior.read().expect("RT prior lock").clone();
                    if let (true, Some(model), Some(prior)) =
                        (params.dia, self.form_rt_model.as_deref(), prior)
                    {
                        index.attach_rt(self.idx, params, model, &self.form_rt_cache, &prior);
                    }
                    index
                });
                (Some(MmapRecordCache::default()), owned_index.as_ref())
            } else {
                let records: FxHashMap<BaseRecordKey, crate::candidate_index::IndexRecord> =
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
                let budget = crate::search_params::MMAP_WINDOW_CACHE_MAX_CANDIDATES;
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

        let queues: Vec<(TopNQueue, Vec<f32>)> = spectra
            .par_iter()
            .enumerate()
            .map(|(local_idx, spec)| {
                let spec_idx = local_idx + spectrum_idx_offset;
                let scorer = scorer_set.for_spec(spec);
                let retention_cap = params.hot_loop_retention_cap();
                let mut queue = TopNQueue::new(retention_cap);

                if spec.peaks.len() < params.min_peaks as usize {
                    return (queue, Vec::new());
                }

                let mut best_raw = i32::MIN;
                let mut best_pep_key: u64 = 0;
                let mut have_best = false;
                let mut second_raw = i32::MIN;
                let mut ref_tracker = RefTracker::default();
                let track_ref = self.idx.has_alternative_entries;

                let mut tailor_hist: FxHashMap<i32, u32> = FxHashMap::default();
                let mut tailor_total: u32 = 0;
                let mut strong_null_stats = OnlineStats::default();

                let charges_to_try: SmallVec<[u8; 4]> = charges_to_try(spec, params);

                let mut scored_per_charge: SmallVec<[(u8, ScoredSpectrum<'_>); 4]> =
                    SmallVec::new();
                for &z in &charges_to_try {
                    if scored_per_charge.iter().all(|(charge, _)| *charge != z) {
                        scored_per_charge.push((z, ScoredSpectrum::new(spec, scorer, z)));
                    }
                }
                let scored_spec_for_charge = |z: u8| {
                    scored_for_charge(&scored_per_charge, z)
                        .expect("scored spectrum exists for candidate charge")
                };

                let shift_ppm = params.precursor_mass_shift_ppm;
                let ms1_gate = if params.dia {
                    self.dia_ms1_scans(spec.rt_seconds)
                } else {
                    None
                };
                let dia_windows: SmallVec<[(u8, f64, f64); 4]> = if params.dia {
                    charges_to_try
                        .iter()
                        .map(|&z| {
                            let (lo, hi) = crate::precursor_matching::dia_neutral_window(
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
                    let voted =
                        params
                            .dia_vote_filter
                            .filter(|_| params.dia)
                            .map(|(w, k)| Spectrum {
                                peaks: windowed_top_k(&spec.peaks, w, k),
                                title: String::new(),
                                ..*spec
                            });
                    let sel = fi.query_where(
                        voted.as_ref().unwrap_or(spec),
                        &charges_to_try,
                        params,
                        scorer.feature_match_tolerance(),
                        params.fragment_index_top_k as usize,
                        crate::search_params::FRAGMENT_INDEX_MIN_MATCHED,
                        crate::fragment_index::vote_peaks_for(scorer.feature_match_tolerance()),
                        &gate,
                    );
                    let relabel = |v: &mut [Candidate]| {
                        relabel_collision_decoys_with(v, &self.target_bare_seqs)
                    };
                    fi.materialise(&sel, self.idx, params, record_cache.as_ref(), &relabel)
                } else {
                    let key = mmap_window_key(spec, params);
                    self.expand_mmap_window_candidates(params, &key, record_cache.as_ref())
                };
                if frag_index.is_some() && self.candidate_filter.is_some() {
                    mmap_cands_owned.retain(|c| self.keeps(c));
                }
                let window_cand_indices: Vec<usize> = (0..mmap_cands_owned.len()).collect();
                let mmap_cands: &[Candidate] = &mmap_cands_owned;
                let cand_slice: &[Candidate] = mmap_cands;
                let resolve_cand =
                    |slot: usize| -> (&Candidate, u32) { (&cand_slice[slot], slot as u32) };

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
                    .map(|_| crate::open_search::sorted_peaks(&spec.peaks));
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
                        let mut iso_errs: SmallVec<[MassError; 4]> = SmallVec::new();
                        if params.dia {
                            let (lo, hi) = dia_windows
                                .iter()
                                .find(|w| w.0 == z)
                                .map_or((f64::MAX, f64::MIN), |w| (w.1, w.2));
                            let mass = cand.peptide.mass();
                            let mono_mz = mass / z as f64 + PROTON;
                            if (lo..=hi).contains(&mass)
                                && ms1_gate.clone().is_none_or(|r| {
                                    self.dia_ms1_has(r, mono_mz, params.dia_ms1_tol_ppm)
                                })
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
                                    shift_ppm,
                                );
                                if let Some(err) = matched {
                                    iso_errs.push(err);
                                }
                            }
                        }
                        if iso_errs.is_empty() {
                            continue;
                        }

                        let mut pin_score =
                            score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da)
                                + cleavage_credit;

                        let mut mod_pep: Option<Peptide> = None;
                        let mut open_iso: Option<i8> = None;
                        if let (Some(open), Some(peaks)) =
                            (params.open_search, open_peaks.as_deref())
                        {
                            let zf = f64::from(z);
                            let observed = adjusted_observed_neutral_mass(
                                spec.precursor_mz * zf - zf * PROTON,
                                params.precursor_mass_shift_ppm,
                            );
                            for off in params.isotope_error_range.clone() {
                                let corrected = observed
                                    - f64::from(off) * model::mass::ISOTOPE;
                                let delta = corrected - cand.peptide.mass();
                                if delta.abs() <= 0.5
                                    || delta < open.delta_min
                                    || delta > open.delta_max
                                {
                                    continue;
                                }
                                let loc_tol =
                                    fragment_tolerance_da.max(2e-5 * cand.peptide.mass());
                                let Some(loc) = crate::open_search::localize_and_refine(
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
                                let shifted = crate::open_search::shifted_peptide(
                                    &cand.peptide,
                                    site0,
                                    loc.delta,
                                );
                                let s2 = score_psm(
                                    scored_spec,
                                    &shifted,
                                    scorer,
                                    z,
                                    fragment_tolerance_da,
                                ) + cleavage_credit;
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
                            if s > second_raw {
                                let mkey = peptidoform_key(&cand.peptide);
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
                                    second_raw = s;
                                }
                            }
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

                        let smaller_ppm = |a: &MassError, b: &MassError| {
                            a.mass_error_ppm
                                .abs()
                                .partial_cmp(&b.mass_error_ppm.abs())
                                .unwrap_or(std::cmp::Ordering::Equal)
                        };
                        let err = match open_iso {
                            Some(want) => iso_errs
                                .iter()
                                .copied()
                                .find(|e| e.isotope_offset == want)
                                .or_else(|| iso_errs.iter().copied().min_by(smaller_ppm)),
                            None if params.open_search.is_some() => iso_errs
                                .iter()
                                .copied()
                                .find(|e| e.isotope_offset == 0)
                                .or_else(|| iso_errs.iter().copied().min_by(smaller_ppm)),
                            None => iso_errs.iter().copied().min_by(smaller_ppm),
                        }
                        .unwrap();

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

                let tailor_denom =
                    crate::psm::tailor_denominator(&tailor_hist, tailor_total) as f32;

                let mut retained_scores: Vec<f32> = queue.iter_psms().map(|p| p.score).collect();
                retained_scores.sort_by(|a, b| b.total_cmp(a));
                let listwise_score_gap =
                    scoring_crate::scoring::listwise_score_gap(&retained_scores);
                let candidate_rank_entropy =
                    scoring_crate::scoring::candidate_rank_entropy(&retained_scores);

                let strong_mode = params.score_mode == ScoreMode::Strong;
                let make_features = |psm: &PsmMatch, full: bool| -> PsmFeatures {
                    let ss = scored_spec_for_charge(psm.charge_used);
                    let cand = &cand_slice[psm.primary_candidate_idx() as usize];
                    let used_shift = psm.features.open_used_shift;
                    let mut open_loc: Option<crate::open_search::Localization> = None;
                    let mut shifted_pep: Option<Peptide> = None;
                    let mut open_delta_ref: Option<f64> = None;
                    {
                        if let (Some(_), Some(peaks)) =
                            (params.open_search, open_peaks.as_deref())
                        {
                            let zf = f64::from(psm.charge_used);
                            let mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
                            let observed = adjusted_observed_neutral_mass(
                                mz * zf - zf * PROTON,
                                params.precursor_mass_shift_ppm,
                            ) - f64::from(psm.isotope_offset) * model::mass::ISOTOPE;
                            let delta = observed - cand.peptide.mass();
                            if delta.abs() > 0.5 {
                                let loc_tol =
                                    fragment_tolerance_da.max(2e-5 * cand.peptide.mass());
                                if let Some(loc) = crate::open_search::localize_and_refine(
                                    &cand.peptide,
                                    peaks,
                                    psm.charge_used,
                                    delta,
                                    loc_tol,
                                    fragment_tolerance_da,
                                ) {
                                    if loc.best_site != 0 {
                                        if used_shift {
                                            shifted_pep =
                                                Some(crate::open_search::shifted_peptide(
                                                    &cand.peptide,
                                                    usize::from(loc.best_site - 1),
                                                    loc.delta,
                                                ));
                                        }
                                        open_delta_ref = Some(loc.delta);
                                    }
                                    open_loc = Some(crate::open_search::Localization {
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
                        if let Some(r) = ref_tracker.best_excluding(peptidoform_key(&cand.peptide))
                        {
                            features.ref_competitor = 1;
                            features.ref_best_raw = r as f32;
                        }
                    }
                    features.listwise_score_gap = listwise_score_gap;
                    features.candidate_rank_entropy = candidate_rank_entropy;
                    features.tailor_score = psm.score / tailor_denom;
                    let credit = cleavage_credit_for(cand, enz, aa_set_for_scoring) as f32;
                    features.rank_score_float =
                        score_psm_float(ss, pep, scorer, psm.charge_used, fragment_tolerance_da)
                            + credit;
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
                    let rows: Vec<crate::site_evidence::Row<'_>> = psms
                        .iter()
                        .map(|psm| {
                            let cand = &cand_slice[psm.primary_candidate_idx() as usize];
                            let z = psm.charge_used;
                            let mz = psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
                            crate::site_evidence::Row {
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
                        .zip(crate::site_evidence::spectrum_evidence(
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
                    queue.remap_candidate_idxs(|local| {
                        accum.intern(&mmap_cands[local as usize])
                    });
                }

                if params.open_search.is_some() {
                    self.open_rescore_attempts
                        .fetch_add(open_attempts, Ordering::Relaxed);
                    self.open_rescore_wins
                        .fetch_add(open_wins, Ordering::Relaxed);
                }
                (queue, retained_strong)
            })
            .collect();

        queues
    }
}

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

pub const OUT_OF_CORE_SLICE_DA: f64 = 80.0;

pub fn match_spectra_out_of_core(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    scorer: &RankScorer,
    fragment_tolerance_da: f64,
    decoy_prefix: &str,
    cache_dir: &std::path::Path,
) -> std::io::Result<(Vec<TopNQueue>, Vec<Candidate>)> {
    let mut params = params.clone();
    params.fragment_index_top_k = crate::search_params::FRAGMENT_INDEX_TOP_K;
    let cache = crate::candidate_index::index_cache_path(idx, &params, cache_dir);
    let mut prepared = PreparedSearch::prepare_mmap(
        idx,
        &params,
        scorer,
        fragment_tolerance_da,
        decoy_prefix,
        &cache,
    )?;
    let mut order: Vec<usize> = (0..spectra.len()).collect();
    order.sort_by(|&a, &b| {
        crate::fragment_index::index_order_mass(&spectra[a])
            .total_cmp(&crate::fragment_index::index_order_mass(&spectra[b]))
            .then(a.cmp(&b))
    });
    let mut queues: Vec<Option<TopNQueue>> = (0..spectra.len()).map(|_| None).collect();
    let mut rest: &[usize] = &order;
    while !rest.is_empty() {
        let take = crate::fragment_index::index_chunk_len(
            rest,
            |&i| crate::fragment_index::index_order_mass(&spectra[i]),
            OUT_OF_CORE_SLICE_DA,
        );
        let (head, tail) = rest.split_at(take);
        let chunk: Vec<Spectrum> = head.iter().map(|&i| spectra[i].clone()).collect();
        for (mut q, &orig) in prepared.run_chunk(&chunk, 0).into_iter().zip(head) {
            q.update_in_place(|psm| psm.spectrum_idx = orig);
            queues[orig] = Some(q);
        }
        rest = tail;
    }
    prepared.sync_materialized_candidates();
    let queues = queues
        .into_iter()
        .map(|q| q.expect("every spectrum is in exactly one chunk"))
        .collect();
    Ok((queues, prepared.candidates))
}

fn coisolation_window(spec: &Spectrum, params: &SearchParams) -> (f64, f64, f64) {
    let lo = spec.precursor_mz
        - spec
            .isolation_lower_offset
            .unwrap_or(params.chimeric_isolation_halfwidth_da);
    let hi = spec.precursor_mz
        + spec
            .isolation_upper_offset
            .unwrap_or(params.chimeric_isolation_halfwidth_da);
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

struct Pass2Scan {
    primary: Peptide,
    coisolated: Vec<crate::coisolation::CoIsolated>,
}

pub fn run_pass2_coisolation(
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
            let coisolated = crate::coisolation::detect_coisolated(
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

const PASS2_POOL_SPAN_NOMINAL: i32 = 32;

fn mmap_secondaries(
    prepared: &PreparedSearch,
    scorer_set: &ScorerSet<'_>,
    spectra: &[Spectrum],
    scans: &[Option<Pass2Scan>],
    params: &SearchParams,
) -> Vec<Vec<(PsmMatch, Candidate)>> {
    use crate::coisolation::{search_secondary, secondary_nominal_window, CoIsolated};
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

#[inline(always)]
#[allow(
    clippy::too_many_arguments,
    reason = "hot-loop driver; all args are orthogonal cleavage parameters"
)]
pub(crate) fn enzyme_cleavage_credit(
    cand: &Candidate,
    enz: Enzyme,
    enz_is_c_term: bool,
    enz_is_n_term: bool,
    credit_neighboring: i32,
    penalty_neighboring: i32,
    credit_peptide: i32,
    penalty_peptide: i32,
) -> i32 {
    let mut score: i32 = 0;
    let pre = cand.peptide.pre;
    let post = cand.peptide.post;
    if enz_is_c_term {
        score += if cand.is_protein_n_term || enz.is_cleavable(pre) {
            credit_neighboring
        } else {
            penalty_neighboring
        };
        let last = match cand.peptide.residues.last() {
            Some(aa) => aa.residue,
            None => 0,
        };
        score += if enz.is_cleavable(last) {
            credit_peptide
        } else {
            penalty_peptide
        };
    } else if enz_is_n_term {
        let first = match cand.peptide.residues.first() {
            Some(aa) => aa.residue,
            None => 0,
        };
        score += if enz.is_cleavable(first) {
            credit_peptide
        } else {
            penalty_peptide
        };
        score += if cand.is_protein_c_term || enz.is_cleavable(post) {
            credit_neighboring
        } else {
            penalty_neighboring
        };
    }
    score
}

pub(crate) fn cleavage_credit_for(cand: &Candidate, enz: Enzyme, aa_set: &AminoAcidSet) -> i32 {
    enzyme_cleavage_credit(
        cand,
        enz,
        enz.is_c_term(),
        enz.is_n_term(),
        aa_set.neighboring_aa_cleavage_credit(),
        aa_set.neighboring_aa_cleavage_penalty(),
        aa_set.peptide_cleavage_credit(),
        aa_set.peptide_cleavage_penalty(),
    )
}

pub(crate) fn compute_psm_features(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
) -> PsmFeatures {
    compute_psm_features_with_edge(scored_spec, peptide, scorer, charge, None, true, None)
}

/// [`predict_frag_intensities`], through `cache` when one is given.
fn predict_frag_with(
    model: &scoring_crate::gbdt_eval::GbdtPeakModel,
    peptide: &Peptide,
    charge: u8,
    cache: Option<&crate::frag_pred_cache::FragPredCache>,
) -> Vec<f64> {
    match cache {
        Some(c) => c.intensities(model, peptide, charge),
        None => predict_frag_intensities(model, peptide, charge),
    }
}

/// The PSM features that only a full feature pass computes: the fragment-intensity LLR
/// battery, the rich-ion LLR and the modification-site features.
#[derive(Default)]
struct KeptPsmFeatureValues {
    frag_pred_explained: f32,
    frag_pred_chance_llr: f32,
    frag_topk_observed: f32,
    rich_ion_llr: f32,
    mod_site: [f32; scoring_crate::mod_site_features::N_MOD_SITE_FEATURES],
}

/// Computes the full-pass-only features of `peptide` at `charge`. `frag_pred` is the
/// fragment-intensity prediction of the scorer's model for this peptide and charge (`None`
/// without a model).
fn kept_psm_feature_values(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    frag_pred: Option<&[f64]>,
) -> KeptPsmFeatureValues {
    let model_tol = scorer.feature_match_tolerance();
    let model_tol_is_ppm = matches!(model_tol, model::tolerance::Tolerance::Ppm(_));
    let model_tol = model_tol.raw_value();
    let (frag_pred_explained, frag_pred_chance_llr, frag_topk_observed) = frag_llr_battery(
        scorer.param().frag_intensity_model.as_deref(),
        scored_spec,
        peptide,
        charge,
        model_tol,
        model_tol_is_ppm,
        frag_pred,
    );
    let rich_ion_llr = rich_ion_llr(
        scorer.param().rich_ion_model.as_deref(),
        scored_spec,
        peptide,
        charge,
        model_tol,
        model_tol_is_ppm,
    );
    let mod_site = mod_site_features(peptide, scored_spec, model_tol, model_tol_is_ppm);
    KeptPsmFeatureValues {
        frag_pred_explained,
        frag_pred_chance_llr,
        frag_topk_observed,
        rich_ion_llr,
        mod_site,
    }
}

/// Completes `features`, computed by a pass with `full = false` for `peptide` at `charge`,
/// with the full-pass-only features, giving the same values as a full pass. `frag_cache`
/// shares fragment-intensity predictions within a chunk.
pub(crate) fn complete_psm_features(
    features: &mut PsmFeatures,
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    frag_cache: Option<&crate::frag_pred_cache::FragPredCache>,
) {
    if peptide.length() < 2 {
        return;
    }
    let frag_pred: Option<Vec<f64>> = scorer
        .param()
        .frag_intensity_model
        .as_deref()
        .map(|g| predict_frag_with(g, peptide, charge, frag_cache));
    let kept = kept_psm_feature_values(scored_spec, peptide, scorer, charge, frag_pred.as_deref());
    features.frag_pred_explained = kept.frag_pred_explained;
    features.frag_pred_chance_llr = kept.frag_pred_chance_llr;
    features.frag_topk_observed = kept.frag_topk_observed;
    features.rich_ion_llr = kept.rich_ion_llr;
    features.mod_site_shifted_matched = kept.mod_site[FEAT_SITE_SHIFTED_MATCHED];
    features.mod_site_shifted_frac = kept.mod_site[FEAT_SITE_SHIFTED_FRAC];
    features.mod_site_intens_frac = kept.mod_site[FEAT_SITE_INTENS_FRAC];
    features.mod_site_localized = kept.mod_site[FEAT_SITE_LOCALIZED];
    features.mod_site_det_count = kept.mod_site[FEAT_SITE_DET_COUNT];
}

pub(crate) fn compute_psm_features_with_edge(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    known_edge_score: Option<i32>,
    full: bool,
    frag_cache: Option<&crate::frag_pred_cache::FragPredCache>,
) -> PsmFeatures {
    let n = peptide.length();
    if n < 2 {
        return PsmFeatures::default();
    }

    let edge_score =
        known_edge_score.unwrap_or_else(|| psm_edge_score(scored_spec, peptide, scorer, charge));

    let predict_losses = scorer
        .param()
        .data_type
        .activation
        .predicts_neutral_losses();
    let predicted = predict_by_ions_with_losses(peptide, 1..=1, predict_losses);
    let mut b_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut y_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut b_rank: SmallVec<[u32; 64]> = smallvec![u32::MAX; n - 1];
    let mut y_rank: SmallVec<[u32; 64]> = smallvec![u32::MAX; n - 1];

    let mut matched_ions: SmallVec<[MatchedIon; 96]> = SmallVec::new();
    let mut matched_ion_ranks: SmallVec<[u32; 96]> = SmallVec::new();

    let feat_tol = scorer.pin_match_tolerance();
    let feature_tol_is_ppm = matches!(feat_tol, model::tolerance::Tolerance::Ppm(_));
    let feature_tol = feat_tol.raw_value();

    let model_tol = scorer.feature_match_tolerance();
    let model_tol_is_ppm = matches!(model_tol, model::tolerance::Tolerance::Ppm(_));
    let model_tol = model_tol.raw_value();

    const H2O_LOSS_DA: f64 = model::mass::H2O;
    const NH3_LOSS_DA: f64 = model::mass::NH3;
    let mut neutral_loss_ion_count: u32 = 0;

    for p in &predicted {
        let tol_da = if feature_tol_is_ppm {
            p.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        if let Some((rank, intensity, peak_mz)) = scored_spec.nearest_peak_full(p.mz, tol_da) {
            let is_b = matches!(p.kind, IonKind::B | IonKind::C);
            matched_ions.push((intensity, peak_mz, p.mz, is_b, FragmentIon::of_predicted(p)));
            matched_ion_ranks.push(rank);

            let pos = (p.position - 1) as usize;
            match p.kind {
                IonKind::B | IonKind::C => {
                    if pos < b_matched.len() {
                        b_matched[pos] = true;
                        if !p.is_loss() && rank < b_rank[pos] {
                            b_rank[pos] = rank;
                        }
                    }
                }
                IonKind::Y | IonKind::Z => {
                    if pos < y_matched.len() {
                        y_matched[pos] = true;
                        if !p.is_loss() && rank < y_rank[pos] {
                            y_rank[pos] = rank;
                        }
                    }
                }
            }

            let ion_charge = 1.0_f64;
            let h2o_target = p.mz - H2O_LOSS_DA / ion_charge;
            let nh3_target = p.mz - NH3_LOSS_DA / ion_charge;
            if scored_spec.nearest_peak_full(h2o_target, tol_da).is_some()
                || scored_spec.nearest_peak_full(nh3_target, tol_da).is_some()
            {
                neutral_loss_ion_count += 1;
            }
        }
    }

    let num_matched: u32 = (b_matched.iter().filter(|&&m| m).count()
        + y_matched.iter().filter(|&&m| m).count()) as u32;

    let predicted_z2 = predict_by_ions_with_losses(peptide, 2..=2, predict_losses);
    let mut doubly_charged_matched_ion_count: u32 = 0;
    for p in &predicted_z2 {
        let tol_da = if feature_tol_is_ppm {
            p.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        if scored_spec.nearest_peak_full(p.mz, tol_da).is_some() {
            doubly_charged_matched_ion_count += 1;
        }
    }

    fn longest_run(matched: &[bool]) -> u32 {
        let mut best = 0u32;
        let mut cur = 0u32;
        for &m in matched {
            if m {
                cur += 1;
                if cur > best {
                    best = cur;
                }
            } else {
                cur = 0;
            }
        }
        best
    }

    let parent_mass = scored_spec.parent_mass();
    let num_segments = scorer.param().num_segments.max(1) as usize;

    let mut b_any_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut y_any_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut sum_prefix_intensity: f64 = 0.0;
    let mut sum_suffix_intensity: f64 = 0.0;

    let mut prm_accurate: f64 = 0.0;
    let mut srm_accurate: f64 = 0.0;

    let segment_ions: SmallVec<[&[scoring_crate::param_model::IonType]; 8]> = (0..num_segments)
        .map(|seg| {
            scorer
                .param()
                .ion_types_for_partition_slice(charge, parent_mass, seg)
        })
        .collect();

    for i in 0..(n - 1) {
        let aa_n = &peptide.residues[i];
        let aa_c = &peptide.residues[n - 1 - i];
        prm_accurate += aa_n.mass + aa_n.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        srm_accurate += aa_c.mass + aa_c.mod_.as_ref().map_or(0.0, |m| m.mass_delta);

        let mut b_any_this = false;
        let mut y_any_this = false;

        for seg in 0..num_segments {
            let ions = segment_ions[seg];
            for &ion in ions {
                let (is_prefix, residue_mass) = match ion {
                    scoring_crate::param_model::IonType::Prefix {
                        charge: ic,
                        offset_bits,
                        ..
                    } => {
                        let offset = f32::from_bits(offset_bits) as f64;
                        let z = ic as f64;
                        (true, (prm_accurate / z + offset, ion))
                    }
                    scoring_crate::param_model::IonType::Suffix {
                        charge: ic,
                        offset_bits,
                        ..
                    } => {
                        let offset = f32::from_bits(offset_bits) as f64;
                        let z = ic as f64;
                        (false, (srm_accurate / z + offset, ion))
                    }
                    scoring_crate::param_model::IonType::Noise => continue,
                };
                let theo_mz = residue_mass.0;
                if scorer.param().segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = if feature_tol_is_ppm {
                    theo_mz * feature_tol / 1e6
                } else {
                    feature_tol
                };
                if let Some((_rank, intensity, _peak_mz)) =
                    scored_spec.nearest_peak_full(theo_mz, tol_da)
                {
                    if is_prefix {
                        sum_prefix_intensity += intensity as f64;
                        b_any_this = true;
                    } else {
                        sum_suffix_intensity += intensity as f64;
                        y_any_this = true;
                    }
                }
            }
        }

        b_any_matched[i] = b_any_this;
        y_any_matched[i] = y_any_this;
    }

    let longest_b = longest_run(&b_any_matched);
    let longest_y = longest_run(&y_any_matched);

    let total_intensity = scored_spec.total_intensity();
    let matched_b_intensity: f64 = sum_prefix_intensity;
    let matched_y_intensity: f64 = sum_suffix_intensity;
    let matched_total = matched_b_intensity + matched_y_intensity;

    let safe_div = |num: f64, denom: f64| -> f32 {
        if denom > 0.0 {
            (num / denom) as f32
        } else {
            0.0
        }
    };

    let explained_ion_current_ratio = safe_div(matched_total, total_intensity);
    let n_term_ion_current_ratio = safe_div(matched_b_intensity, total_intensity);
    let c_term_ion_current_ratio = safe_div(matched_y_intensity, total_intensity);
    let ms2_ion_current = if total_intensity > 0.0 {
        total_intensity as f32
    } else {
        0.0
    };
    let isolation_window_efficiency = 0.0_f32;

    matched_ions.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let top7 = &matched_ions[..matched_ions.len().min(7)];

    let abs_ppm_errors: Vec<f64> = top7
        .iter()
        .filter(|&&(_, _, pred, _, _)| pred > 0.0)
        .map(|&(_, obs, pred, _, _)| ((obs - pred) / pred * 1e6).abs())
        .collect();
    let rel_ppm_errors: Vec<f64> = top7
        .iter()
        .filter(|&&(_, _, pred, _, _)| pred > 0.0)
        .map(|&(_, obs, pred, _, _)| (obs - pred) / pred * 1e6)
        .collect();

    fn mean_and_pop_stdev(values: &[f64]) -> (f32, f32) {
        if values.is_empty() {
            return (0.0, 0.0);
        }
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let sum_sq: f64 = values.iter().map(|v| v * v).sum();
        let var = (sum_sq / n - mean * mean).max(0.0);
        (mean as f32, var.sqrt() as f32)
    }

    let (mean_error_top7, stdev_error_top7) = mean_and_pop_stdev(&abs_ppm_errors);
    let (mean_rel_error_top7, stdev_rel_error_top7) = mean_and_pop_stdev(&rel_ppm_errors);

    let (ppm_centre, ppm_sigma) = scorer.ppm_gaussian();
    let ppm_gaussian_score: f32 = matched_ions
        .iter()
        .filter(|&&(_, _, pred, _, _)| pred > 0.0)
        .map(|&(_, obs, pred, _, _)| {
            let z = ((obs - pred) / pred * 1e6 - ppm_centre) / ppm_sigma;
            (-0.5 * z * z).exp()
        })
        .sum::<f64>() as f32;

    let chance_match_surprise: f32 = matched_ions
        .iter()
        .filter(|&&(_, _, pred, _, _)| pred > 0.0)
        .map(|&(_, obs, pred, _, _)| {
            let tol_da = if feature_tol_is_ppm {
                pred * feature_tol / 1e6
            } else {
                feature_tol
            };
            let rho = scored_spec.local_peak_density(obs, DENSITY_HW);
            let p_chance = rho * (2.0 * tol_da);
            if p_chance > 0.0 {
                (-p_chance.ln()).max(0.0)
            } else {
                0.0
            }
        })
        .sum::<f64>() as f32;

    let mut complementary_sites: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    for i in 1..n {
        complementary_sites[i - 1] = b_any_matched[i - 1] && y_any_matched[n - i - 1];
    }
    let longest_complementary_ladder = longest_run(&complementary_sites);

    let mut complementary_ion_balance = 0.0_f32;
    for i in 1..n {
        let rb = b_rank[i - 1];
        let ry = y_rank[n - i - 1];
        if rb != u32::MAX && ry != u32::MAX {
            let diff = (rb as i64 - ry as i64).unsigned_abs() as f32;
            complementary_ion_balance += 1.0 / (1.0 + diff);
        }
    }

    let mean_matched_intensity_rank = if matched_ion_ranks.is_empty() {
        0.0
    } else {
        matched_ion_ranks.iter().map(|&r| r as f32).sum::<f32>() / matched_ion_ranks.len() as f32
    };

    let charge1_predicted = predict_by_ions(peptide, 1..=1);
    let theoretical_ions = theoretical_fragment_ions(&charge1_predicted);
    let unique_match_fraction = if num_matched == 0 {
        0.0
    } else {
        let evidence_sum: f64 = matched_ions
            .iter()
            .map(|&(_, obs, _, _, ion)| {
                let ambiguity = fragment_ambiguity(
                    &theoretical_ions,
                    obs,
                    feature_tol,
                    feature_tol_is_ppm,
                    ion,
                );
                1.0 / (1.0 + ambiguity as f64)
            })
            .sum();
        (evidence_sum / num_matched as f64) as f32
    };
    let mass_competition_evidence_val = mass_competition_evidence(
        scored_spec,
        &matched_ions,
        &theoretical_ions,
        feature_tol,
        feature_tol_is_ppm,
    );

    let frag_intensity_model = scorer.param().frag_intensity_model.as_deref();
    let frag_pred: Option<Vec<f64>> = frag_intensity_model
        .filter(|_| peptide.length() >= 2)
        .map(|g| predict_frag_with(g, peptide, charge, frag_cache));
    let frag_pred_slice = frag_pred.as_deref();
    let intensity_signal_val = intensity_signal(
        frag_intensity_model,
        scored_spec,
        peptide,
        charge,
        model_tol,
        model_tol_is_ppm,
        frag_pred_slice,
    );

    let kept = if full {
        kept_psm_feature_values(scored_spec, peptide, scorer, charge, frag_pred_slice)
    } else {
        KeptPsmFeatureValues::default()
    };
    let KeptPsmFeatureValues {
        frag_pred_explained,
        frag_pred_chance_llr,
        frag_topk_observed,
        rich_ion_llr: rich_ion_llr_val,
        mod_site,
    } = kept;

    PsmFeatures {
        num_matched_main_ions: num_matched,
        longest_b,
        longest_y,
        longest_y_pct: longest_y as f32 / n as f32,
        matched_ion_ratio: num_matched as f32 / n as f32,
        explained_ion_current_ratio,
        n_term_ion_current_ratio,
        c_term_ion_current_ratio,
        ms2_ion_current,
        isolation_window_efficiency,
        mean_error_top7,
        stdev_error_top7,
        mean_rel_error_top7,
        stdev_rel_error_top7,
        edge_score,
        precursor_isotope_kl: 0.0,
        precursor_snr: 0.0,
        delta_raw_score: 0.0,
        tailor_score: 0.0,
        rank_score_float: 0.0,
        ppm_gaussian_score,
        longest_complementary_ladder,
        complementary_ion_balance,
        neutral_loss_ion_count,
        mean_matched_intensity_rank,
        doubly_charged_matched_ion_count,
        chance_match_surprise,
        unique_match_fraction,
        intensity_signal: intensity_signal_val,
        frag_pred_explained,
        frag_pred_chance_llr,
        frag_topk_observed,
        rich_ion_llr: rich_ion_llr_val,
        is_refinement: 0,
        num_mods: 0,
        refine_mod_class: 0,
        mod_site_shifted_matched: mod_site[FEAT_SITE_SHIFTED_MATCHED],
        mod_site_shifted_frac: mod_site[FEAT_SITE_SHIFTED_FRAC],
        mod_site_intens_frac: mod_site[FEAT_SITE_INTENS_FRAC],
        mod_site_localized: mod_site[FEAT_SITE_LOCALIZED],
        mod_site_det_count: mod_site[FEAT_SITE_DET_COUNT],
        mass_competition_evidence: mass_competition_evidence_val,
        candidate_rank_entropy: 0.0,
        listwise_score_gap: 0.0,
        strong_score: 0.0,
        strong_score_cal: 0.0,
        delta_rt: 0.0,
        abs_delta_rt: 0.0,
        delta_rt_norm: 0.0,
        ref_competitor: 0,
        ref_best_raw: 0.0,
        is_variant: 0,
        is_isoform: 0,
        is_chain: 0,
        delta_ref_score: 0.0,
        isotope_mimic: 0,
        predicted_rt_min: 0.0,
        isobaric_rt_margin: 0.0,
        open_site: 0,
        open_loc_score: 0.0,
        open_loc_margin: 0.0,
        open_delta: 0.0,
        open_used_shift: false,
    }
}

#[cfg(test)]
mod accumulator_tests {
    use super::*;
    use model::{AminoAcidSetBuilder, Protein, ProteinDb};

    #[test]
    fn sync_resets_the_identity_map() {
        let db = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
                annotation: None,
            }],
        };
        let idx = SearchIndex::from_target_db(&db, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = crate::coisolation::tests::tiny_scorer();
        let mut prepared =
            PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
        let cands = prepared.out_of_core_candidates_where(|_, _| true);
        assert!(cands.len() >= 2);
        {
            let mut accum = prepared.mmap_accum.lock().unwrap();
            assert_eq!(accum.intern(&cands[0]), 0);
            assert_eq!(accum.intern(&cands[1]), 1);
        }
        prepared.sync_materialized_candidates();
        assert_eq!(prepared.candidates.len(), 2);
        let mut accum = prepared.mmap_accum.lock().unwrap();
        assert!(accum.index_of.is_empty());
        assert_eq!(accum.intern(&cands[1]), 0);
        assert_eq!(accum.candidates.len(), 1);
    }

    #[test]
    fn dia_ms1_gate_reads_neighbouring_survey_scans() {
        let db = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRK".to_vec(),
                annotation: None,
            }],
        };
        let idx = SearchIndex::from_target_db(&db, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = crate::coisolation::tests::tiny_scorer();
        let mut prepared =
            PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
        let scan = |rt: f64, mz: f64| model::scan::Ms1Scan {
            rt,
            peaks: vec![(mz, 1.0)],
            ..Default::default()
        };
        prepared.set_dia_ms1(vec![
            scan(30.0, 700.0),
            scan(10.0, 500.0),
            scan(20.0, 600.0),
        ]);
        assert_eq!(prepared.dia_ms1_scans(None), None);
        let r = prepared.dia_ms1_scans(Some(9.0)).unwrap();
        assert_eq!(r, 0..2);
        assert!(prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 4e-6), 5.0));
        assert!(!prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 6e-6), 5.0));
        assert!(prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 6e-6), 25.0));
        assert!(!prepared.dia_ms1_has(r, 700.0, 25.0));
        assert_eq!(prepared.dia_ms1_scans(Some(26.0)).unwrap(), 1..3);

        let e = prepared
            .dia_ms1_error_ppm(Some(20.0), 600.0 * (1.0 - 7e-6), 25.0)
            .unwrap();
        assert!((e - 7.0).abs() < 1e-3, "{e}");
        assert_eq!(prepared.dia_ms1_error_ppm(Some(20.0), 650.0, 25.0), None);
        prepared.shift_dia_ms1_ppm(7.0);
        let e = prepared
            .dia_ms1_error_ppm(Some(20.0), 600.0 * (1.0 - 7e-6), 25.0)
            .unwrap();
        assert!(e.abs() < 1e-3, "{e}");
    }

    #[test]
    fn windowed_top_k_keeps_the_k_most_intense_per_window_in_mz_order() {
        let peaks = [
            (101.0, 1.0),
            (150.0, 5.0),
            (180.0, 3.0),
            (199.0, 4.0),
            (250.0, 2.0),
        ];
        assert_eq!(
            windowed_top_k(&peaks, 100.0, 2),
            vec![(150.0, 5.0), (199.0, 4.0), (250.0, 2.0)]
        );
        assert_eq!(windowed_top_k(&peaks, 100.0, 9), peaks.to_vec());
    }
}

#[cfg(test)]
mod per_scan_scorer_tests {
    use super::*;
    use model::{AminoAcidSetBuilder, Protein, ProteinDb};
    use std::sync::atomic::AtomicUsize;

    fn spectrum(activation: Option<ActivationMethod>) -> Spectrum {
        Spectrum {
            precursor_mz: 500.0,
            peaks: vec![(200.0, 10.0)],
            activation_method: activation,
            ..Spectrum::default()
        }
    }

    fn prepared_with<'a>(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
    ) -> PreparedSearch<'a> {
        PreparedSearch::prepare_in_memory(idx, params, scorer, 0.02, "XXX").unwrap()
    }

    #[test]
    fn scans_route_by_their_own_activation() {
        let db = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
                annotation: None,
            }],
        };
        let idx = SearchIndex::from_target_db(&db, "XXX");
        let params =
            SearchParams::default_tryptic(AminoAcidSetBuilder::new_standard().build().unwrap());
        let scorer = crate::coisolation::tests::tiny_scorer();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let prepared = prepared_with(&idx, &params, &scorer).with_per_scan_scorers(
            ActivationMethod::HCD,
            Box::new(move |a| {
                seen.fetch_add(1, Ordering::Relaxed);
                (a == ActivationMethod::EThcD).then(crate::coisolation::tests::tiny_scorer)
            }),
        );
        let chunk = [
            spectrum(Some(ActivationMethod::HCD)),
            spectrum(Some(ActivationMethod::EThcD)),
            spectrum(None),
            spectrum(Some(ActivationMethod::ETD)),
            spectrum(Some(ActivationMethod::EThcD)),
        ];
        let set = prepared.scorer_set(&chunk);
        let own = |s: &Spectrum| std::ptr::eq(set.for_spec(s), &scorer);
        assert!(own(&chunk[0]));
        assert!(!own(&chunk[1]), "EThcD scan has its own scorer");
        assert!(own(&chunk[2]));
        assert!(
            own(&chunk[3]),
            "a declined activation keeps the run's scorer"
        );
        assert!(!own(&chunk[4]));
        let _ = prepared.scorer_set(&chunk);
        assert_eq!(calls.load(Ordering::Relaxed), 2);
        let mut routes = prepared.per_scan_routes();
        routes.sort_by_key(|r| r.0 as u8);
        assert_eq!(
            routes.iter().map(|r| (r.0, r.1, r.2)).collect::<Vec<_>>(),
            vec![
                (ActivationMethod::ETD, false, 2),
                (ActivationMethod::EThcD, true, 4)
            ]
        );
    }

    #[test]
    fn no_routing_keeps_the_run_scorer() {
        let db = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
                annotation: None,
            }],
        };
        let idx = SearchIndex::from_target_db(&db, "XXX");
        let params =
            SearchParams::default_tryptic(AminoAcidSetBuilder::new_standard().build().unwrap());
        let scorer = crate::coisolation::tests::tiny_scorer();
        let prepared = prepared_with(&idx, &params, &scorer);
        let chunk = [spectrum(Some(ActivationMethod::EThcD))];
        let set = prepared.scorer_set(&chunk);
        assert!(std::ptr::eq(set.for_spec(&chunk[0]), &scorer));
        assert!(prepared.per_scan_routes().is_empty());
    }
}

#[cfg(test)]
mod chargeless_unit_tests {
    use super::*;
    use model::{AminoAcid, AminoAcidSetBuilder, Protein, ProteinDb};
    use scoring_crate::scoring::fragment_ions::predict_by_ions;

    fn peptide(residues: &[u8]) -> Peptide {
        let aas: Vec<AminoAcid> = residues
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(aas, b'K', b'K')
    }

    fn spectrum(precursor_mz: f64, charge: Option<i32>, peptides: &[&Peptide]) -> Spectrum {
        let mut peaks: Vec<(f64, f32)> = peptides
            .iter()
            .flat_map(|p| predict_by_ions(p, 1..=1))
            .map(|ion| (ion.mz, 100.0_f32))
            .collect();
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        Spectrum {
            title: "chargeless".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: charge,
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        }
    }

    #[test]
    fn chargeless_spectrum_is_searched_per_charge() {
        let db = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"GKDAFLGSFLYEYSRKSAFLWEYSPLLTNNGVRK".to_vec(),
                annotation: None,
            }],
        };
        let idx = SearchIndex::from_target_db(&db, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.charge_range = 2..=5;
        params.fragment_index_top_k = 1;
        let scorer = crate::coisolation::tests::tiny_scorer();
        let prepared =
            PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();

        let truth = peptide(b"DAFLGSFLYEYSR");
        let bait = peptide(b"SAFLWEYSPLLTNNGVR");
        let mz3 = (truth.mass() + 3.0 * PROTON) / 3.0;
        let spectra = [
            spectrum(mz3, Some(3), &[&truth, &bait]),
            spectrum(mz3, None, &[&truth, &bait]),
        ];
        let mut queues = prepared.run_chunk(&spectra, 7);
        assert_eq!(queues.len(), 2);
        for (i, q) in queues.iter_mut().enumerate() {
            let psms = q.drain_into_vec();
            let best = psms
                .iter()
                .max_by(|a, b| a.rank_score.total_cmp(&b.rank_score))
                .unwrap_or_else(|| panic!("spectrum {i} came back empty"));
            assert_eq!(best.spectrum_idx, 7 + i);
            assert_eq!(best.charge_used, 3, "spectrum {i}");
            let accum = prepared.mmap_accum.lock().unwrap();
            let pep = &accum.candidates[best.primary_candidate_idx() as usize].peptide;
            assert_eq!(pep.residues, truth.residues, "spectrum {i}");
        }
    }
}

#[cfg(test)]
mod charges_to_try_tests {
    use super::charges_to_try;
    use crate::search_params::SearchParams;
    use model::aa_set::AminoAcidSetBuilder;
    use model::Spectrum;

    fn params() -> SearchParams {
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        SearchParams::default_tryptic(aa)
    }

    fn spec(charge: Option<i32>) -> Spectrum {
        Spectrum {
            precursor_charge: charge,
            ..Default::default()
        }
    }

    #[test]
    fn reported_charge_is_searched_alone() {
        let p = params();
        assert_eq!(charges_to_try(&spec(Some(5)), &p).as_slice(), &[5]);
    }

    #[test]
    fn out_of_range_charge_falls_back_to_the_range() {
        let p = params();
        assert_eq!(charges_to_try(&spec(Some(300)), &p).as_slice(), &[2, 3]);
        assert_eq!(charges_to_try(&spec(Some(-2)), &p).as_slice(), &[2, 3]);
        assert_eq!(charges_to_try(&spec(Some(255)), &p).as_slice(), &[255]);
    }

    #[test]
    fn missing_charge_still_falls_back_to_the_range() {
        let p = params();
        assert_eq!(charges_to_try(&spec(None), &p).as_slice(), &[2, 3]);
    }
}

#[cfg(test)]
mod feature_tests {
    use super::*;
    use model::activation::ActivationMethod;
    use model::amino_acid::AminoAcid;
    use model::instrument::InstrumentType;
    use model::mass::PROTON;
    use model::peptide::Peptide;
    use model::protocol::Protocol;
    use model::spectrum::Spectrum;
    use model::tolerance::Tolerance;
    use rustc_hash::FxHashMap;
    use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Partition, SpecDataType};
    use scoring_crate::scoring::fragment_ions::predict_by_ions;
    use scoring_crate::scoring::ScoredSpectrum;

    fn make_scorer(tol_da: f64) -> RankScorer {
        use model::mass::{H2O, PROTON};
        let part = Partition {
            charge: 2,
            parent_mass: 0.0,
            seg_num: 0,
        };
        let prefix1 = IonType::Prefix {
            charge: 1,
            offset_bits: (PROTON as f32).to_bits(),
            loss_class: 0,
        };
        let suffix1 = IonType::Suffix {
            charge: 1,
            offset_bits: ((H2O + PROTON) as f32).to_bits(),
            loss_class: 0,
        };
        let noise = IonType::Noise;
        let mut ion_table = FxHashMap::default();
        ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
        ion_table.insert(suffix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
        ion_table.insert(noise, vec![0.1_f32, 0.2, 0.3, 0.4]);
        let mut rank_dist_table = FxHashMap::default();
        rank_dist_table.insert(part, ion_table);
        let mut frag_off_table = FxHashMap::default();
        frag_off_table.insert(
            part,
            vec![
                FragmentOffsetFrequency {
                    ion_type: prefix1,
                    frequency: 0.7,
                },
                FragmentOffsetFrequency {
                    ion_type: suffix1,
                    frequency: 0.7,
                },
            ],
        );
        let mut param = scoring_crate::Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Da(tol_da),
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

    fn ala_peptide(len: usize) -> Peptide {
        let aa = AminoAcid::standard(b'A').unwrap();
        Peptide::new(vec![aa; len], b'_', b'-')
    }

    fn make_spectrum(peaks: Vec<(f64, f32)>) -> Spectrum {
        Spectrum {
            title: "test".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        }
    }

    #[test]
    fn recentred_spectra_removes_the_measured_offset() {
        let spec = make_spectrum(vec![(300.0, 1.0), (1000.0, 2.0)]);
        assert!(recentred_spectra(std::slice::from_ref(&spec), 0.0).is_none());
        let theo = [300.0, 1000.0];
        let low = make_spectrum(theo.iter().map(|m| (m * (1.0 - 8e-6), 1.0)).collect());
        let c = recentred_spectra(std::slice::from_ref(&low), -8.0).unwrap();
        for (p, t) in c[0].peaks.iter().zip(theo) {
            assert!((p.0 - t).abs() / t < 1e-9, "{} vs {t}", p.0);
        }
        assert_eq!(low.peaks[0].0, 300.0 * (1.0 - 8e-6));
        assert!(c[0].peaks.windows(2).all(|w| w[0].0 < w[1].0));
    }

    #[test]
    fn compute_psm_features_top7_error_stats_zero_when_no_matches() {
        let pep = ala_peptide(4);
        let spec = make_spectrum(vec![]);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.5), 2);
        assert_eq!(
            f.mean_error_top7, 0.0,
            "mean_error_top7 should be 0 with no matches"
        );
        assert_eq!(
            f.stdev_error_top7, 0.0,
            "stdev_error_top7 should be 0 with no matches"
        );
        assert_eq!(
            f.mean_rel_error_top7, 0.0,
            "mean_rel_error_top7 should be 0 with no matches"
        );
        assert_eq!(
            f.stdev_rel_error_top7, 0.0,
            "stdev_rel_error_top7 should be 0 with no matches"
        );
        assert_eq!(
            f.explained_ion_current_ratio, 0.0,
            "ratio should be 0 with no peaks"
        );
        assert_eq!(
            f.ms2_ion_current, 0.0,
            "ms2_ion_current should be 0 with no peaks"
        );
    }

    #[test]
    fn compute_psm_features_populates_ion_current_ratios() {
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);

        let mut peaks: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        peaks.push((1500.0, 5.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = make_spectrum(peaks);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.01), 2);

        assert!(
            f.explained_ion_current_ratio > 0.0,
            "explained_ion_current_ratio should be > 0 when ions match, got {}",
            f.explained_ion_current_ratio
        );
        assert!(
            f.n_term_ion_current_ratio > 0.0,
            "n_term_ion_current_ratio should be > 0 when b-ions match"
        );
        assert!(
            f.c_term_ion_current_ratio > 0.0,
            "c_term_ion_current_ratio should be > 0 when y-ions match"
        );

        let sum = f.n_term_ion_current_ratio + f.c_term_ion_current_ratio;
        assert!(
            (sum - f.explained_ion_current_ratio).abs() < 1e-5,
            "n_term + c_term should == explained ({} + {} != {})",
            f.n_term_ion_current_ratio,
            f.c_term_ion_current_ratio,
            f.explained_ion_current_ratio
        );

        let total: f32 = ss.total_intensity() as f32;
        assert!(
            (f.ms2_ion_current - total).abs() < 1.0,
            "ms2_ion_current {} should match total spectrum intensity {}",
            f.ms2_ion_current,
            total
        );

        assert_eq!(f.isolation_window_efficiency, 0.0);
    }

    #[test]
    fn compute_psm_features_unique_match_fraction() {
        let scorer = make_scorer(0.01);

        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);
        let mut peaks: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_unique = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(peaks)),
            &pep,
            &scorer,
            2,
        );
        assert!(
            f_unique.unique_match_fraction > 0.7,
            "separated ions should yield a high fraction, got {}",
            f_unique.unique_match_fraction
        );

        const H2O: f64 = model::mass::H2O;
        const NH3: f64 = model::mass::NH3;
        let mut contested_case: Option<(Peptide, f64, f64)> = None;
        'search: for seq in ["PFAH", "PFSR", "PFLM", "PFML"] {
            let residues: Vec<AminoAcid> = seq
                .bytes()
                .map(|r| AminoAcid::standard(r).unwrap())
                .collect();
            let pep2 = Peptide::new(residues, b'_', b'-');
            let pred2 = predict_by_ions(&pep2, 1..=1);
            let mut theo: Vec<f64> = Vec::new();
            for p in &pred2 {
                theo.push(p.mz);
                theo.push(p.mz - H2O);
                theo.push(p.mz - NH3);
            }
            for p in &pred2 {
                let tol = p.mz * 20e-6;
                let near_real = pred2.iter().filter(|q| (q.mz - p.mz).abs() < 0.05).count();
                let others: Vec<f64> = theo
                    .iter()
                    .copied()
                    .filter(|t| (t - p.mz).abs() > 1e-6 && (t - p.mz).abs() <= tol)
                    .collect();
                if near_real == 1 && others.len() == 1 {
                    contested_case = Some((pep2, p.mz, others[0]));
                    break 'search;
                }
            }
        }
        let (pep2, peak_mz, other_mz) =
            contested_case.expect("should find a contested theo pair among the candidate peptides");
        let tol = peak_mz * 20e-6;
        assert!(
            (other_mz - peak_mz).abs() <= tol,
            "contested pair should be within feature tolerance"
        );
        let contested = vec![(peak_mz, 100.0_f32)];
        let f_contested = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(contested)),
            &pep2,
            &scorer,
            2,
        );
        assert!(
            f_contested.unique_match_fraction < f_unique.unique_match_fraction,
            "contested peak should lower fraction (unique={}, contested={})",
            f_unique.unique_match_fraction,
            f_contested.unique_match_fraction
        );
        assert!(
            (f_contested.unique_match_fraction - 0.5).abs() < 0.15,
            "two ions on one peak → fraction ~0.5, got {}",
            f_contested.unique_match_fraction
        );
    }

    #[test]
    fn compute_psm_features_doubly_charged_matched_ion_count() {
        let pep = ala_peptide(4);
        let predicted_z2 = predict_by_ions(&pep, 2..=2);
        assert!(
            !predicted_z2.is_empty(),
            "charge-2 ions should exist for 4-mer"
        );

        let f_empty = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(vec![])),
            &pep,
            &make_scorer(0.01),
            3,
        );
        assert_eq!(f_empty.doubly_charged_matched_ion_count, 0);

        let mut peaks: Vec<(f64, f32)> = predicted_z2
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(peaks)),
            &pep,
            &make_scorer(0.01),
            3,
        );
        assert!(
            f.doubly_charged_matched_ion_count > 0,
            "charge-2 peaks should match, got {}",
            f.doubly_charged_matched_ion_count
        );
        assert_eq!(
            f.doubly_charged_matched_ion_count,
            predicted_z2.len() as u32,
            "all charge-2 ions placed → count should equal {}",
            predicted_z2.len()
        );
    }

    #[test]
    fn compute_psm_features_mean_matched_intensity_rank() {
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);

        let mut top: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, 1000.0 - i as f32))
            .collect();
        top.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_top = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(top)),
            &pep,
            &make_scorer(0.01),
            2,
        );
        assert!(
            (f_top.mean_matched_intensity_rank - 2.5).abs() < 0.1,
            "ions as the only peaks → mean rank should be ~2.5, got {}",
            f_top.mean_matched_intensity_rank
        );

        let mut buried: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, 10.0 + i as f32))
            .collect();
        for j in 0..10 {
            buried.push((1500.0 + j as f64 * 10.0, 500.0 + j as f32 * 50.0));
        }
        buried.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_buried = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(buried)),
            &pep,
            &make_scorer(0.01),
            2,
        );
        assert!(
            f_buried.mean_matched_intensity_rank > f_top.mean_matched_intensity_rank,
            "buried ions should have higher mean rank (top={}, buried={})",
            f_top.mean_matched_intensity_rank,
            f_buried.mean_matched_intensity_rank
        );
    }

    #[test]
    fn compute_psm_features_longest_complementary_ladder() {
        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);
        let n = pep.length();

        let mut all: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_all = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(all)),
            &pep,
            &make_scorer(0.01),
            2,
        );
        let expected = (n - 1) as u32;
        assert_eq!(
            f_all.longest_complementary_ladder, expected,
            "all ions matched → ladder should equal n-1={expected}, got {}",
            f_all.longest_complementary_ladder
        );

        let mut split: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .filter(|(_, p)| !(matches!(p.kind, IonKind::Y) && p.position == 3))
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        split.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_split = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(split)),
            &pep,
            &make_scorer(0.01),
            2,
        );
        assert_eq!(
            f_split.longest_complementary_ladder, 2,
            "missing middle bond should leave a length-2 ladder, got {}",
            f_split.longest_complementary_ladder
        );
        assert!(
            f_split.longest_complementary_ladder < f_all.longest_complementary_ladder,
            "split ladder ({}) should be shorter than full ({})",
            f_split.longest_complementary_ladder,
            f_all.longest_complementary_ladder
        );
    }

    #[test]
    fn compute_psm_features_complementary_ion_balance() {
        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);

        let mut all: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_all = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(all)),
            &pep,
            &make_scorer(0.01),
            2,
        );
        assert!(
            f_all.complementary_ion_balance > 0.0,
            "all ions matched → balance should be > 0, got {}",
            f_all.complementary_ion_balance
        );

        let mut b_only: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .filter(|(_, p)| matches!(p.kind, IonKind::B))
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        b_only.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_b = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(b_only)),
            &pep,
            &make_scorer(0.01),
            2,
        );
        assert_eq!(
            f_b.complementary_ion_balance, 0.0,
            "b-only spectrum → no complementary pair → balance 0, got {}",
            f_b.complementary_ion_balance
        );
    }

    #[test]
    fn compute_psm_features_neutral_loss_ion_count() {
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);
        let b2 = predicted
            .iter()
            .find(|p| matches!(p.kind, IonKind::B) && p.position == 2)
            .expect("b2 should exist");

        let mut peaks = vec![(b2.mz, 100.0_f32), (b2.mz - 18.0106, 50.0_f32)];
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(peaks)),
            &pep,
            &make_scorer(0.01),
            2,
        );
        assert_eq!(
            f.neutral_loss_ion_count, 1,
            "b2 + b2−H2O partner should yield neutral_loss_ion_count=1, got {}",
            f.neutral_loss_ion_count
        );
    }

    #[test]
    fn compute_psm_features_strong_score_stage1_bolt_ons() {
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);

        let mut exact: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        exact.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_exact = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(exact)),
            &pep,
            &make_scorer(0.01),
            2,
        );

        assert!(
            (f_exact.ppm_gaussian_score - 4.0).abs() < 0.05,
            "ppm_gaussian_score should be ~4 for 4 exact matches, got {}",
            f_exact.ppm_gaussian_score
        );
        assert_eq!(
            f_exact.longest_complementary_ladder, 2,
            "both cleavage sites should form a length-2 ladder, got {}",
            f_exact.longest_complementary_ladder
        );

        let mut shifted: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz * (1.0 + 15e-6), (i + 1) as f32 * 10.0))
            .collect();
        shifted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f_shift = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(shifted)),
            &pep,
            &make_scorer(0.01),
            2,
        );
        assert!(
            f_shift.ppm_gaussian_score < f_exact.ppm_gaussian_score * 0.5,
            "ppm_gaussian_score should fall sharply when matches scatter (exact={}, shifted={})",
            f_exact.ppm_gaussian_score,
            f_shift.ppm_gaussian_score
        );
    }

    #[test]
    fn compute_psm_features_error_stats_nonzero_when_ions_match_with_offset() {
        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);

        let offset_da = 0.0005_f64;
        let mut peaks: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz + offset_da, (i + 1) as f32 * 10.0))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = make_spectrum(peaks);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.05), 2);

        assert!(
            f.mean_error_top7 > 0.0,
            "mean_error_top7 should be > 0 when peaks are systematically offset, got {}",
            f.mean_error_top7
        );
        assert!(
            f.stdev_error_top7 < 20.0,
            "stdev_error_top7 should be small (single-digit ppm) for identical-Da offset, got {}",
            f.stdev_error_top7
        );
        assert!(
            f.mean_rel_error_top7 != 0.0,
            "mean_rel_error_top7 should be nonzero when peaks are offset"
        );
    }

    #[test]
    fn ms2_ion_current_equals_total_intensity() {
        let pep = ala_peptide(3);
        let peaks = vec![(100.0, 50.0_f32), (200.0, 30.0), (300.0, 20.0)];
        let spec = make_spectrum(peaks.clone());
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.5), 2);

        let expected: f32 = peaks.iter().map(|&(_, i)| i).sum();
        assert_eq!(
            f.ms2_ion_current, expected,
            "ms2_ion_current {} should equal sum of peak intensities {}",
            f.ms2_ion_current, expected
        );
    }

    #[test]
    fn b1_mz_for_alanine_is_proton_plus_residue_mass() {
        use model::amino_acid::AminoAcid;
        let aa = AminoAcid::standard(b'A').unwrap();
        let residue_mass = aa.mass;
        let expected_b1_mz = residue_mass + PROTON;
        let pep = ala_peptide(2);
        let predicted = predict_by_ions(&pep, 1..=1);
        let b1 = predicted
            .iter()
            .find(|p| matches!(p.kind, IonKind::B) && p.position == 1)
            .expect("b1 ion should exist");
        assert!(
            (b1.mz - expected_b1_mz).abs() < 1e-6,
            "b1 mz {} expected {}",
            b1.mz,
            expected_b1_mz
        );
    }
}

pub(crate) fn dedup_pepseq_score(psms: Vec<PsmMatch>, candidates: &[Candidate]) -> Vec<PsmMatch> {
    use std::collections::btree_map::Entry;
    use std::collections::BTreeMap;

    let mut pep_key_cache: FxHashMap<u32, Arc<PepDedupKey>> = FxHashMap::default();
    let mut groups: BTreeMap<DedupMapKey, PsmMatch> = BTreeMap::new();

    let mut psms = psms;
    psms.sort_by(|a, b| b.cmp(a));
    for psm in psms {
        let primary = psm.primary_candidate_idx();
        let pep_key = pep_key_cache
            .entry(primary)
            .or_insert_with(|| {
                Arc::new(PepDedupKey::from_peptide(
                    &candidates[primary as usize].peptide,
                ))
            })
            .clone();
        let key = DedupMapKey {
            pep: pep_key,
            score: psm.rank_score.round() as i32,
            is_decoy: candidates[primary as usize].is_decoy,
        };

        match groups.entry(key) {
            Entry::Vacant(slot) => {
                slot.insert(psm);
            }
            Entry::Occupied(mut slot) => {
                let existing = slot.get_mut();
                merge_unique_candidate_idxs(&mut existing.candidate_idxs, &psm.candidate_idxs);
                if psm.rank_score > existing.rank_score {
                    let merged_idxs = std::mem::take(&mut existing.candidate_idxs);
                    let mut survivor = psm;
                    merge_unique_candidate_idxs(&mut survivor.candidate_idxs, &merged_idxs);
                    *existing = survivor;
                }
            }
        }
    }

    let mut out = Vec::with_capacity(groups.len());
    out.extend(groups.into_values());
    out
}

#[inline]
pub(crate) fn peptidoform_key(peptide: &model::Peptide) -> u64 {
    use std::hash::Hasher;
    let mut h = rustc_hash::FxHasher::default();
    for aa in &peptide.residues {
        h.write_u8(aa.residue);
        let units = aa
            .mod_
            .as_ref()
            .map(|m| (m.mass_delta * 100_000.0).round() as i64)
            .unwrap_or(0);
        h.write_i64(units);
    }
    h.finish()
}

/// I and L are isobaric, so a decoy peptide that equals a target peptide with I and L folded is
/// the same precursor and fragment ladder as that target.
fn fold_isoleucine(r: u8) -> u8 {
    if r == b'I' {
        b'L'
    } else {
        r
    }
}

/// The residues of a candidate with I folded to L, the key of the target peptide set.
fn bare_residues(cand: &Candidate) -> Box<[u8]> {
    cand.peptide
        .residues
        .iter()
        .map(|aa| fold_isoleucine(aa.residue))
        .collect()
}

fn folded_set(spans: Vec<&[u8]>) -> FxHashSet<Box<[u8]>> {
    let mut set = FxHashSet::with_capacity_and_hasher(spans.len(), Default::default());
    set.extend(
        spans
            .into_iter()
            .map(|s| s.iter().map(|&r| fold_isoleucine(r)).collect::<Box<[u8]>>()),
    );
    set
}

fn relabel_collision_decoys_with(candidates: &mut [Candidate], target_seqs: &FxHashSet<Box<[u8]>>) {
    if target_seqs.is_empty() {
        return;
    }
    for cand in candidates.iter_mut().filter(|c| c.is_decoy) {
        if target_seqs.contains(&bare_residues(cand)) {
            cand.is_decoy = false;
        }
    }
}

fn build_target_bare_seqs(mi: &MmapCandidateIndex, db: &SearchIndex) -> FxHashSet<Box<[u8]>> {
    use crate::candidate_index::flags as rec_flags;
    let mut spans: Vec<&[u8]> = (0..mi.len())
        .into_par_iter()
        .filter_map(|i| {
            let rec = mi.record_at(i);
            if rec.flags & rec_flags::IS_DECOY != 0 {
                return None;
            }
            let seq = &db.db.proteins.get(rec.protein_index as usize)?.sequence;
            let start = rec.start_offset as usize;
            seq.get(start..start + rec.length as usize)
        })
        .collect();
    spans.par_sort_unstable();
    spans.dedup();
    folded_set(spans)
}

fn build_target_bare_seqs_where(
    mi: &MmapCandidateIndex,
    db: &SearchIndex,
    params: &SearchParams,
    keep: &CandidateFilter,
) -> FxHashSet<Box<[u8]>> {
    use crate::candidate_index::flags as rec_flags;
    let mut spans: Vec<&[u8]> = (0..mi.len())
        .into_par_iter()
        .filter_map(|i| {
            let rec = mi.record_at(i);
            if rec.flags & rec_flags::IS_DECOY != 0 {
                return None;
            }
            let seq = &db.db.proteins.get(rec.protein_index as usize)?.sequence;
            let start = rec.start_offset as usize;
            let span = seq.get(start..start + rec.length as usize)?;
            expand_base_record(db, params, &rec)
                .iter()
                .any(keep)
                .then_some(span)
        })
        .collect();
    spans.par_sort_unstable();
    spans.dedup();
    folded_set(spans)
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct PepDedupKey {
    residues: Vec<u8>,
    mod_units: Vec<i32>,
}

impl PepDedupKey {
    fn from_peptide(peptide: &model::Peptide) -> Self {
        let len = peptide.residues.len();
        let mut residues = Vec::with_capacity(len);
        let mut mod_units = Vec::with_capacity(len);
        for aa in &peptide.residues {
            residues.push(aa.residue);
            mod_units.push(
                aa.mod_
                    .as_ref()
                    .map(|m| (m.mass_delta * 100_000.0).round() as i32)
                    .unwrap_or(0),
            );
        }
        Self {
            residues,
            mod_units,
        }
    }
}

#[derive(Clone, PartialEq, Eq)]
struct DedupMapKey {
    pep: Arc<PepDedupKey>,
    score: i32,
    is_decoy: bool,
}

impl PartialOrd for DedupMapKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DedupMapKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.pep
            .residues
            .cmp(&other.pep.residues)
            .then_with(|| self.pep.mod_units.cmp(&other.pep.mod_units))
            .then(self.score.cmp(&other.score))
            .then(self.is_decoy.cmp(&other.is_decoy))
    }
}

fn merge_unique_candidate_idxs(into: &mut crate::psm::CandidateIdxs, from: &[u32]) {
    if from.is_empty() {
        return;
    }
    let mut seen: FxHashSet<u32> = into.iter().copied().collect();
    for &idx in from {
        if seen.insert(idx) {
            into.push(idx);
        }
    }
}

#[cfg(test)]
mod dedup_tests {
    use super::*;
    use crate::psm::PsmMatch;
    use model::amino_acid::AminoAcid;
    use model::modification::{ModLocation, Modification};
    use model::peptide::Peptide;
    use model::ResidueSpec;

    fn seq_peptide(bytes: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = bytes
            .iter()
            .filter_map(|&b| AminoAcid::standard(b))
            .collect();
        Peptide::new(residues, b'R', b'K')
    }

    fn cand_with_peptide(peptide: Peptide) -> Candidate {
        Candidate {
            peptide,
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn credit(cand: &Candidate, enz: Enzyme) -> i32 {
        enzyme_cleavage_credit(
            cand,
            enz,
            enz.is_c_term(),
            enz.is_n_term(),
            10,
            -100,
            1,
            -1000,
        )
    }

    fn flanked(seq: &[u8], pre: u8, post: u8) -> Candidate {
        let residues: Vec<AminoAcid> = seq.iter().filter_map(|&b| AminoAcid::standard(b)).collect();
        cand_with_peptide(Peptide::new(residues, pre, post))
    }

    #[test]
    fn n_term_enzyme_peptide_credit_reads_the_first_residue() {
        assert_eq!(credit(&flanked(b"DAAAG", b'G', b'D'), Enzyme::AspN), 11);
        assert_eq!(credit(&flanked(b"AAAAG", b'D', b'D'), Enzyme::AspN), -990);
        assert_eq!(credit(&flanked(b"DAAAG", b'G', b'A'), Enzyme::AspN), -99);
        assert_eq!(credit(&flanked(b"KAAAG", b'G', b'K'), Enzyme::LysN), 11);
    }

    #[test]
    fn c_term_enzyme_credit_reads_the_last_residue() {
        assert_eq!(credit(&flanked(b"AAAAR", b'K', b'G'), Enzyme::Trypsin), 11);
        assert_eq!(
            credit(&flanked(b"AAAAG", b'G', b'G'), Enzyme::Trypsin),
            -1100
        );
    }

    fn psm(primary: u32, rank: f32, pin: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: smallvec::smallvec![primary],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: pin,
            rank_score: rank,
            edge_score: (rank - pin) as i32,
            activation_method: None,
            features: Default::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        }
    }

    #[test]
    fn collision_decoy_relabeled_as_target_genuine_decoy_kept() {
        let cand = |seq: &[u8], is_decoy: bool| Candidate {
            peptide: seq_peptide(seq),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        };
        let mut cands = vec![
            cand(b"PEPTIDEK", false),
            cand(b"PEPTIDEK", true),
            cand(b"EDITPEPK", true),
            cand(b"PEPTLDEK", true),
        ];
        let targets = folded_set(vec![b"PEPTIDEK".as_slice()]);
        relabel_collision_decoys_with(&mut cands, &targets);
        assert!(!cands[0].is_decoy, "target unchanged");
        assert!(!cands[1].is_decoy, "collision decoy relabeled to target");
        assert!(cands[2].is_decoy, "genuine decoy left as decoy");
        assert!(
            !cands[3].is_decoy,
            "a decoy equal to a target with I and L folded is relabeled to target"
        );
    }

    #[test]
    fn dedup_uses_rank_score_not_pin_score() {
        let pep = seq_peptide(b"PEPTK");
        let cands = vec![cand_with_peptide(pep.clone())];
        let psms = vec![psm(0, 100.4, 99.6), psm(0, 120.0, 99.6)];
        let out = dedup_pepseq_score(psms, &cands);
        assert_eq!(out.len(), 2, "different rank_score keys must not merge");
    }

    #[test]
    fn dedup_distinguishes_mod_state() {
        let mut ox = seq_peptide(b"PEPMK");
        ox.residues[3].mod_ = Some(model::modification::leak_mod(Modification {
            name: "Ox".into(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        }));
        let cands = vec![
            cand_with_peptide(seq_peptide(b"PEPMK")),
            cand_with_peptide(ox),
        ];
        let psms = vec![psm(0, 50.0, 50.0), psm(1, 50.0, 50.0)];
        let out = dedup_pepseq_score(psms, &cands);
        assert_eq!(out.len(), 2, "mod-aware pepSeq keys must differ");
    }

    #[test]
    fn dedup_keeps_highest_rank_score_survivor() {
        let pep = seq_peptide(b"PEPTK");
        let cands = vec![cand_with_peptide(pep)];
        let psms = vec![psm(0, 59.6, 50.0), psm(0, 60.4, 50.0)];
        let out = dedup_pepseq_score(psms, &cands);
        assert_eq!(out.len(), 1);
        assert!((out[0].rank_score - 60.4).abs() < f32::EPSILON);
    }
}

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
mod ref_tracker_tests {
    use super::RefTracker;

    #[test]
    fn competitor_excludes_the_rows_own_peptidoform() {
        let mut t = RefTracker::default();
        t.fold(30, 1);
        t.fold(25, 2);
        t.fold(28, 1);
        assert_eq!(t.best_excluding(1), Some(25));
        assert_eq!(t.best_excluding(2), Some(30));
        assert_eq!(t.best_excluding(3), Some(30));
        t.fold(40, 3);
        assert_eq!(t.best_excluding(3), Some(30));
        assert_eq!(RefTracker::default().best_excluding(1), None);
        let mut one = RefTracker::default();
        one.fold(30, 1);
        assert_eq!(one.best_excluding(1), None);
        assert_eq!(one.best_excluding(2), Some(30));
    }
}
