//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::{smallvec, SmallVec};
use std::sync::atomic::{AtomicUsize, Ordering};

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
    edge_scoring_enabled, frag_llr_battery, fuse_strong_score, intensity_signal,
    mass_competition_evidence, predict_frag_intensities, psm_edge_score, rich_ion_llr, score_psm,
    score_psm_float, scored_for_charge, strong_score_calibrated, OnlineStats, RankScorer,
    ScoredSpectrum, StrongScoreInputs, DENSITY_HW,
};

// One-time-built state shared across every chunk of a streamed search.
//
// `PreparedSearch::prepare_mmap` opens (or builds) the out-of-core base-peptide
// index once; `PreparedSearch::run_chunk` scores any slice of `Spectrum`s against
// it, materialising each spectrum's candidates lazily from the index.

/// Keeps or drops a candidate after retrieval; see
/// [`PreparedSearch::with_candidate_filter`].
pub type CandidateFilter = Box<dyn Fn(&Candidate) -> bool + Send + Sync>;

pub struct PreparedSearch<'a> {
    pub idx: &'a SearchIndex,
    pub params: &'a SearchParams,
    pub scorer: &'a RankScorer,
    pub fragment_tolerance_da: f64,
    /// The candidates the emitted PSMs reference, indexed by
    /// `PsmMatch::candidate_idxs`. Starts empty and is filled from the scan's
    /// accumulator by [`Self::sync_materialized_candidates`].
    pub candidates: Vec<Candidate>,
    /// Out-of-core base-peptide index the per-spectrum candidates are
    /// materialised from.
    pub mmap_index: MmapCandidateIndex,
    /// The per-spectrum materialized candidates the surviving PSMs reference,
    /// appended (under lock) during the parallel scan and keyed for
    /// de-duplication so a candidate scored across several chunks / spectra gets
    /// one stable global index. Drained into `self.candidates` by
    /// [`Self::sync_materialized_candidates`].
    mmap_accum: Mutex<MmapAccumulator>,
    /// `params.aa_set` with the search enzyme registered for cleavage
    /// scoring. Cheap to clone, but we keep one shared copy here.
    pub aa_set_for_scoring: AminoAcidSet,
    /// Optional MS1 linkage for the chimeric precursor isotope features.
    /// `None` unless the binary supplied one via
    /// [`Self::with_ms1_link`] under `--chimeric`. When `None` (the default
    /// and the entire `--chimeric off` path), the feature fill skips the
    /// precursor-envelope computation and the two new `PsmFeatures` fields
    /// stay 0.0 — keeping the off path bit-identical.
    pub ms1_link: Option<Ms1Link>,
    /// The GLOBAL set of target bare-peptide sequences (every non-decoy base
    /// record in the index), built once at prepare time. The per-spectrum
    /// collision-decoy relabel uses THIS set so a decoy is relabeled iff its bare
    /// sequence matches ANY DB target. Relabeling from only the window candidates
    /// would be WRONG: a decoy whose target twin lies outside the current
    /// precursor window would keep its `is_decoy` flag. Bare sequences only, so
    /// this is bounded (far smaller than the full candidate list).
    target_bare_seqs: FxHashSet<Box<[u8]>>,
    /// Optional post-retrieval candidate filter (see
    /// [`Self::with_candidate_filter`]).
    candidate_filter: Option<CandidateFilter>,
}

/// Global candidate identity for the scan accumulator: enough to dedup a
/// peptidoform reached from multiple spectra/charges/isotope-offsets to one
/// stable global index. The peptidoform mod state is folded into `mod_units`
/// (per-residue mod-mass in 1e-5 Da units) so two placements of the same residues
/// at the same span but different mods get distinct indices.
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

/// Global candidate store of a scan: maps a candidate's identity to its stable
/// global index in `candidates`, growing as the scan materializes new candidates.
#[derive(Default)]
struct MmapAccumulator {
    candidates: Vec<Candidate>,
    index_of: std::collections::HashMap<GlobalCandKey, u32>,
}

impl MmapAccumulator {
    /// Intern a candidate, returning its stable global index. If an identical
    /// candidate was already seen, returns the existing index without pushing.
    fn intern(&mut self, cand: Candidate) -> u32 {
        let key = GlobalCandKey::from_candidate(&cand);
        if let Some(&idx) = self.index_of.get(&key) {
            return idx;
        }
        let idx = self.candidates.len() as u32;
        self.candidates.push(cand);
        self.index_of.insert(key, idx);
        idx
    }
}

/// Which of `enumerate_protein`'s walks emits copy `copy` of `cand`'s span: the
/// whole protein (pass 0) or, for a protein starting with Met, the walk over
/// `sequence[1..]` (pass 1). Copy `k` of a span reachable from both walks is
/// pass `k`; a span at offset 1 carrying the protein N-terminus flag is reached
/// only by the Met-cleaved walk (the flag is set iff the span starts at that
/// walk's offset), so it is pass 1 even though it has a single copy.
fn enumeration_pass(cand: &Candidate, copy: u32) -> u32 {
    if cand.is_protein_n_term && cand.start_offset_in_protein == 1 {
        1
    } else {
        copy
    }
}

/// Candidate indices sorted by exact neutral mass, ties by index.
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

/// How many window slots ahead [`prefetch_window_candidates`] requests a
/// candidate's struct, and how many ahead its residue array (whose address is read
/// from the struct, so it is requested once the struct has had time to arrive).
const PREFETCH_STRUCT_AHEAD: usize = 16;
const PREFETCH_RESIDUES_AHEAD: usize = 8;
/// Cache lines of a residue array requested per candidate (8 lines cover a
/// 21-residue peptide).
const PREFETCH_RESIDUE_LINES: usize = 8;

/// Ask the CPU to start loading the candidates the window scan reaches a few
/// iterations from now. The window visits candidates in index order, which is
/// scattered across the whole candidate pool, so without this every candidate
/// costs a cache miss on its struct and another on its heap residue array. A
/// prefetch is only a hint: it never changes what the loop reads.
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
            // both lie inside the borrowed `Candidate`.
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
                // the residue array (`line * 64 < bytes`).
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

/// Derive the inclusive `[min_nominal, max_nominal]` nominal-mass bucket bounds
/// for one spectrum at one charge state `z`: the candidates whose nominal residue
/// mass lies inside are the ones retrieved for it.
///
/// The window is centered on the precursor m/z and widened by the precursor
/// tolerance plus the isotope-error range, then extended to cover every peptide
/// mass `matches_precursor` accepts at some isotope offset. The nominal centre
/// alone can miss one: the scaled, rounded nominal mass drifts from the residue
/// mass by more than a bucket over an isotope shift on heavy peptides (a 2887 Da
/// residue mass observed at isotope +2 centres the window on 2890 and starts it
/// at 2888).
pub(crate) fn candidate_nominal_bounds(
    spec: &Spectrum,
    z: u8,
    params: &SearchParams,
    shift_ppm: f64,
) -> (i32, i32) {
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
    // `matches_precursor` uses `tolerance.left` when the peptide is LIGHTER than
    // the observed mass and `tolerance.right` when HEAVIER. So the lower nominal
    // bound (lightest acceptable peptide) widens by `left`, and the upper bound
    // (heaviest) by `right`. (Inert when left == right; matters for the
    // calibrator's asymmetric tolerance and any Da-asymmetric config — using the
    // wrong side there prunes real candidates before the exact filter.)
    let min_nominal = nominal_center - iso_max - widen_left;
    let max_nominal = nominal_center - iso_min + widen_right;
    // The lightest and heaviest peptide masses `matches_precursor` accepts over
    // the isotope range. `nominal_from` is monotone, so their nominal residue
    // masses bound every accepted peptide's.
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

/// The charge states a spectrum is searched at: the reported precursor charge is
/// trusted as the ONLY charge, else the configured range. A charge that does not
/// fit a `u8` (or is not positive) counts as unknown.
fn charges_to_try(spec: &Spectrum, params: &SearchParams) -> SmallVec<[u8; 4]> {
    match spec.precursor_charge.and_then(|z| u8::try_from(z).ok()) {
        Some(z) if z > 0 => smallvec![z],
        _ => params.charge_range.clone().collect(),
    }
}

/// Cache key for the out-of-core per-chunk candidate cache: the ordered list of
/// nominal-mass windows this spectrum enumerates, one per charge state.
///
/// `shift_ppm` is `params.precursor_mass_shift_ppm`, a run-level constant, so it
/// is fully captured by `params` and needs no separate key component.
fn mmap_window_key(spec: &Spectrum, params: &SearchParams) -> MmapWindowKey {
    let shift_ppm = params.precursor_mass_shift_ppm;
    charges_to_try(spec, params)
        .into_iter()
        .map(|z| candidate_nominal_bounds(spec, z, params, shift_ppm))
        .collect()
}

/// One spectrum's nominal-mass windows, one per charge state tried.
type MmapWindowKey = SmallVec<[(i32, i32); 4]>;

/// Per-chunk map from a base record to ALL of its peptidoforms (no window
/// filter, no multiplicity copies, decoy label as stored). Keyed on the base
/// record rather than on the window because windows overlap heavily (±tolerance
/// × isotope offsets × charge states) while a base record's expansion is the
/// same in every window that contains it — so one expansion serves every
/// spectrum in the chunk that reaches the record. Bounded by total candidates
/// ([`crate::search_params::MMAP_WINDOW_CACHE_MAX_CANDIDATES`]).
type MmapRecordCache = FxHashMap<BaseRecordKey, Vec<Candidate>>;

impl<'a> PreparedSearch<'a> {
    /// Build the per-search state over the out-of-core base-peptide index.
    ///
    /// Opens the mass-sorted base-peptide index cached at `index_path`, building
    /// it (via
    /// [`build_base_peptide_index`](crate::candidate_index::build_base_peptide_index))
    /// on a cache miss, or in memory when the cache cannot be written.
    /// `self.candidates` starts empty and is grown lazily during the scan so the
    /// PIN/TSV writers resolve `candidate_idxs` against a real candidate slice.
    ///
    /// `decoy_prefix` MUST match the prefix the `SearchIndex` was built with so
    /// the index's per-record decoy flags align with the database.
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
        if mmap_index.is_in_memory() {
            eprintln!(
                "candidate-index: built {} records in memory (not cached)",
                mmap_index.len()
            );
        } else if was_built {
            eprintln!(
                "candidate-index: built and cached {} records → {}",
                mmap_index.len(),
                index_path.display()
            );
        } else {
            eprintln!(
                "candidate-index: reusing cached index ({} records) from {}",
                mmap_index.len(),
                index_path.display()
            );
        }
        Ok(Self::over_index(
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            mmap_index,
        ))
    }

    /// As [`Self::prepare_mmap`], with the index built in memory and never
    /// cached: for small, run-specific databases such as the refinement
    /// cascade's peptide-anchored index.
    pub fn prepare_in_memory(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        decoy_prefix: &str,
    ) -> std::io::Result<Self> {
        let mmap_index = MmapCandidateIndex::build_in_memory(idx, params, decoy_prefix)?;
        Ok(Self::over_index(
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            mmap_index,
        ))
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
        // The GLOBAL target bare-sequence set, built once from the full index.
        // Bare residues are reconstructed from each non-decoy record's protein
        // span — bounded, far smaller than the candidate list.
        let target_bare_seqs = build_target_bare_seqs(&mmap_index, idx);
        PreparedSearch {
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            candidates: Vec::new(),
            mmap_index,
            mmap_accum: Mutex::new(MmapAccumulator::default()),
            aa_set_for_scoring,
            ms1_link: None,
            target_bare_seqs,
            candidate_filter: None,
        }
    }

    /// Drain the scan accumulator into `self.candidates` so downstream
    /// consumers (PIN/TSV writers) resolve `candidate_idxs` against the final
    /// materialized candidate slice. The accumulator is left empty, identity map
    /// included, so a later scan interns from index 0 instead of returning
    /// indices into the drained store.
    pub fn sync_materialized_candidates(&mut self) {
        let mut accum = self.mmap_accum.lock().unwrap();
        let materialized = std::mem::take(&mut accum.candidates);
        accum.index_of = Default::default();
        drop(accum);
        self.candidates = materialized;
    }

    /// Score only the retrieved candidates `keep` accepts. The filter runs on
    /// each spectrum's retrieved candidates (after the fragment-index shortlist
    /// when that retrieval is on), and on [`Self::out_of_core_candidates_where`]
    /// and the chimeric Pass-2 pools, so every scored candidate passed it.
    ///
    /// The collision-decoy target set is rebuilt from the targets the filter
    /// keeps: a span none of whose peptidoforms passes is not a target of this
    /// search, so it must not relabel a decoy with the same sequence.
    pub fn with_candidate_filter(mut self, keep: CandidateFilter) -> Self {
        self.target_bare_seqs =
            build_target_bare_seqs_where(&self.mmap_index, self.idx, self.params, &keep);
        self.candidate_filter = Some(keep);
        self
    }

    /// Whether `cand` passes the optional candidate filter.
    fn keeps(&self, cand: &Candidate) -> bool {
        self.candidate_filter.as_ref().is_none_or(|keep| keep(cand))
    }

    /// Attach an [`Ms1Link`] for the chimeric precursor isotope features.
    /// Builder-style; only the binary's `--chimeric` path calls this. Has NO effect unless
    /// `params.chimeric` is also set (the feature fill double-guards on it).
    pub fn with_ms1_link(mut self, ms1_link: Option<Ms1Link>) -> Self {
        self.ms1_link = ms1_link;
        self
    }

    /// Score one chunk of spectra in parallel using the prepared candidate
    /// state. Returns one `TopNQueue` per input spectrum, in input order.
    ///
    /// The `spectrum_idx_offset` is the index of `spectra[0]` in the overall
    /// stream of spectra being searched. It is written into every emitted
    /// `PsmMatch::spectrum_idx` so the downstream PIN/TSV writers can still
    /// look up the right spectrum metadata in the concatenated metadata
    /// vector.
    pub fn run_chunk(&self, spectra: &[Spectrum], spectrum_idx_offset: usize) -> Vec<TopNQueue> {
        self.run_chunk_with_params(spectra, spectrum_idx_offset, self.params)
    }

    /// Like [`Self::run_chunk`] but uses `params` for precursor matching and
    /// queue sizing instead of the params reference stored at prepare time.
    /// Candidate enumeration and mass buckets are unchanged.
    ///
    /// In fragment-index mode a spectrum without a precursor charge is searched
    /// as one unit per charge in `params.charge_range`, each at its own
    /// precursor mass, and the units' PSMs are merged into that spectrum's
    /// queue. Searched whole, its windows span every charge (m/z x 2 to m/z x 5),
    /// so the chunk's index covered that whole mass range and the index's top-K
    /// was spent across all of it rather than on any one precursor window.
    pub fn run_chunk_with_params(
        &self,
        spectra: &[Spectrum],
        spectrum_idx_offset: usize,
        params: &SearchParams,
    ) -> Vec<TopNQueue> {
        let has_charge = |s: &Spectrum| s.precursor_charge.is_some_and(|z| z > 0);
        if params.fragment_index_top_k == 0
            || params.charge_range.clone().count() < 2
            || spectra.iter().all(has_charge)
        {
            return self.run_chunk_inner(spectra, spectrum_idx_offset, params);
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
        while !units.is_empty() {
            let take = crate::fragment_index::index_chunk_len(&units, |u| mass(&u.1), span_da);
            let (owners, chunk): (Vec<usize>, Vec<Spectrum>) = units.drain(..take).unzip();
            for (mut q, i) in self
                .run_chunk_inner(&chunk, 0, params)
                .into_iter()
                .zip(owners)
            {
                for mut psm in q.drain_into_vec() {
                    psm.spectrum_idx = spectrum_idx_offset + i;
                    merged[i].push(psm);
                }
            }
        }
        merged
    }

    /// The fully-prepared candidate list for one set of nominal-mass `windows`
    /// (one entry per charge state to try): every peptidoform whose nominal
    /// residue mass lies in a window, each replicated by its enumeration
    /// multiplicity, in database enumeration order, with collision decoys
    /// relabelled and the candidate filter applied.
    ///
    /// The result is a pure function of `(windows, params, index)` — it depends
    /// on the spectrum ONLY through those nominal bounds — which is what makes
    /// the per-chunk window cache in `run_chunk_inner` result-identical.
    /// `cache` holds base-record expansions shared across the chunk; records
    /// missing from it are expanded here.
    fn expand_mmap_window_candidates(
        &self,
        params: &SearchParams,
        windows: &[(i32, i32)],
        cache: Option<&MmapRecordCache>,
    ) -> Vec<Candidate> {
        let mi = &self.mmap_index;

        // Union of base records over the spectrum's (charge, isotope-offset)
        // windows, each with its enumeration multiplicity. Deduplicating at
        // the record level is exact: a record's copies all sit at one mass, so
        // every window either holds all of them or none, and two records never
        // yield the same peptidoform (the key includes protein and offset).
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
        // `enumerate_candidates` order: protein ASC, start ASC, length ASC, flags ASC.
        records.sort_unstable_by_key(|(r, _)| base_record_key(r));

        // A spectrum scores exactly the candidates whose nominal residue mass
        // lands in the coarse integer window `candidate_nominal_bounds(spec, z,
        // ..)`. The fetch above over-fetches base records over a padded mass
        // window; acceptance is gated per peptidoform on the exact predicate
        // `nominal_from(peptide.mass() − H2O) ∈ [min, max]` for ANY of the
        // spectrum's windows.
        let in_any_window = |c: &Candidate| {
            let key = c.peptide.nominal_residue_mass();
            windows.iter().any(|&(lo, hi)| key >= lo && key <= hi)
        };

        // Replicate each peptidoform by its base record's multiplicity, so the
        // per-spectrum candidate set matches the database enumeration
        // copy-for-copy (collapsed by the later `dedup_pepseq_score`; the PIN
        // `Proteins` column repeats the accession once per copy).
        //
        // The copies are NOT adjacent in the enumeration: `enumerate_protein`
        // walks the WHOLE protein at `seq_offset = 0`, then walks it AGAIN at
        // `seq_offset = 1` (N-terminal-Met cleavage). So the order inside one
        // protein is `[every pass-0 span, offset ASC] ++ [every pass-1 span,
        // offset ASC]`. Tag each copy with its PASS index (`enumeration_pass`)
        // and sort on it, so the two walks land in that block order.
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
        // Canonical enumeration order: protein_index ASC, then PASS (see above),
        // then start offset ASC, then span length ASC — exactly the order
        // `enumerate_protein` emits spans in. Within each such group candidates
        // keep the `expand_mod_combinations` order they were pushed in
        // (`sort_by` is stable), which is the order `enumerate_candidates`
        // uses. DO NOT add residues/mod_units tie-breaking here: lexicographic
        // residue order diverges from the recursive expansion order, and the
        // queue's tie-breaks follow candidate order.
        tagged.sort_by(|a, b| (a.0, a.1, a.2, a.3).cmp(&(b.0, b.1, b.2, b.3)));
        let mut mmap_cands: Vec<Candidate> = tagged.into_iter().map(|t| t.4).collect();

        // Collision-decoy relabeling against the GLOBAL target bare-sequence set
        // (built once at prepare time from the full index): a decoy is relabeled
        // iff its bare sequence matches ANY DB target. Order-independent, so it
        // runs after the sort.
        relabel_collision_decoys_with(&mut mmap_cands, &self.target_bare_seqs);
        if self.candidate_filter.is_some() {
            mmap_cands.retain(|c| self.keeps(c));
        }
        mmap_cands
    }

    /// Every candidate of the out-of-core index whose bare span passes `keep`,
    /// in database enumeration order (protein, enumeration pass, start, length,
    /// then mod expansion order), with collision decoys relabelled against the
    /// targets of the WHOLE database and the candidate filter applied. `keep`
    /// sees a span's residues and the residue after it (`b'-'` at the protein
    /// C-terminus), so it cannot depend on variable mods; every peptidoform of
    /// a kept span is returned.
    pub fn out_of_core_candidates_where(
        &self,
        keep: impl Fn(&[u8], u8) -> bool + Sync,
    ) -> Vec<Candidate> {
        let mi = &self.mmap_index;
        let db = self.idx;
        let params = self.params;
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
        // Identical copies of a span (the N-terminal-Met re-enumeration) share a
        // key; their count is the enumeration multiplicity.
        kept.par_sort_unstable_by_key(base_record_key);
        let mut records: Vec<(crate::candidate_index::IndexRecord, u32)> = Vec::new();
        for rec in kept {
            match records.last_mut() {
                Some((last, n)) if base_record_key(last) == base_record_key(&rec) => *n += 1,
                _ => records.push((rec, 1)),
            }
        }
        // Same ordering as `expand_mmap_window_candidates`: protein, enumeration
        // pass, start, length; a span's peptidoforms keep their expansion order
        // because the sort is stable and records arrive in key order.
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

    /// The peptide of the candidate at GLOBAL index `global_idx`: an entry of
    /// the scan accumulator (read under its lock) before
    /// [`Self::sync_materialized_candidates`], of `self.candidates` after it.
    pub fn candidate_peptide(&self, global_idx: u32) -> Peptide {
        let accum = self.mmap_accum.lock().unwrap();
        match accum.candidates.get(global_idx as usize) {
            Some(c) => c.peptide.clone(),
            None => self.candidates[global_idx as usize].peptide.clone(),
        }
    }

    /// Whether the candidate at GLOBAL index `global_idx` is a decoy (same
    /// lookup as [`Self::candidate_peptide`]).
    pub fn candidate_is_decoy(&self, global_idx: u32) -> bool {
        let accum = self.mmap_accum.lock().unwrap();
        match accum.candidates.get(global_idx as usize) {
            Some(c) => c.is_decoy,
            None => self.candidates[global_idx as usize].is_decoy,
        }
    }

    fn run_chunk_inner(
        &self,
        spectra: &[Spectrum],
        spectrum_idx_offset: usize,
        params: &SearchParams,
    ) -> Vec<TopNQueue> {
        let scorer = self.scorer;
        let fragment_tolerance_da = self.fragment_tolerance_da;
        let aa_set_for_scoring = &self.aa_set_for_scoring;

        // Out-of-core per-chunk expansion cache, keyed by BASE RECORD. Every
        // (spectrum, charge) search unit re-derives its candidate list from the
        // on-disk index, and consecutive windows overlap heavily, so the same
        // base peptide would otherwise be re-expanded for nearly every spectrum
        // that reaches it — on a PTM-rich search (phospho S/T/Y, four variable
        // mods) that expansion was ~95% of all CPU (issue #72). Fetch the
        // distinct base records of every window in the chunk (index range
        // scans, cheap), expand each ONCE in parallel, and let the per-spectrum
        // assembly below look them up. The cache is a pure memo, so scored
        // output never changes; it is bounded by total candidates — records
        // admitted after the budget is spent are dropped and expanded on
        // demand, and which ones land in the cache under parallel contention
        // is not deterministic.
        let (record_cache, frag_index): (Option<MmapRecordCache>, Option<ChunkFragmentIndex>) = {
            let mi = &self.mmap_index;
            let mut windows: Vec<(i32, i32)> = spectra
                .iter()
                .filter(|s| s.peaks.len() >= params.min_peaks as usize)
                .flat_map(|s| mmap_window_key(s, params))
                .collect();
            windows.sort_unstable();
            windows.dedup();
            // Fragment-ion index over the chunk's records (issue #76): the
            // forms are enumerated once here, and each spectrum then scores
            // only the forms its peaks vote for, instead of every form in its
            // precursor windows. Index mode indexes every form in the chunk's
            // whole mass interval and each spectrum queries only its own
            // window, so one range scan over the union of the windows yields
            // every record it can select.
            if params.fragment_index_top_k > 0 {
                let frag_index = {
                    // The chunk's precursor mass interval over every spectrum,
                    // charge and isotope offset; only forms inside it are indexed.
                    let shift_ppm = params.precursor_mass_shift_ppm;
                    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
                    for s in spectra
                        .iter()
                        .filter(|s| s.peaks.len() >= params.min_peaks as usize)
                    {
                        for z in charges_to_try(s, params) {
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
                    // Records none of whose forms can reach the interval are
                    // dropped during the fetch, before the records are sorted.
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
                    ChunkFragmentIndex::build(
                        recs,
                        self.idx,
                        params,
                        self.scorer.feature_match_tolerance(),
                        lo,
                        hi,
                    )
                };
                // The record cache serves the enumeration path only; index
                // mode materialises through the bounded walk and never reads it.
                (Some(MmapRecordCache::default()), Some(frag_index))
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

        // Parallel per-spectrum search. All inputs above are `&` immutable; the
        // closure owns its TopNQueue, scored_per_charge cache, and per-bin GF state.
        let queues: Vec<TopNQueue> = spectra
            .par_iter()
            .enumerate()
            .map(|(local_idx, spec)| {
                let spec_idx = local_idx + spectrum_idx_offset;
                let retention_cap = params.hot_loop_retention_cap();
                let mut queue = TopNQueue::new(retention_cap);

                // Cascade Pass 1 is a NARROW precursor-window search, identical to the
                // non-chimeric path; Pass 2 (`run_pass2_coisolation`) is the only
                // chimeric-specific enumeration. MS1 load + the
                // `precursor_isotope_match` feature fill still run for Pass 2.

                // Skip spectra with too few peaks.
                if spec.peaks.len() < params.min_peaks as usize {
                    return queue;
                }

                // DeltaRawScore capture (additive feature): track the best and the
                // second-best *distinct-peptide* rounded RawScore across ALL
                // mass-matching candidates for this spectrum — independent of the
                // TopNQueue, so it survives `top_n = 1` (where the runner-up is
                // evicted) AND never feeds the GF `min_score` (so no emitted PSM's
                // SpecEValue is perturbed). Distinctness is keyed on the
                // PEPTIDOFORM (residues + per-residue mod deltas), NOT on nominal mass.
                //
                // It was previously keyed on `nominal_residue_mass()`, on the stated
                // assumption that isobaric different-sequence peptides at the same
                // nominal mass are "rare". That is inverted: the candidate window IS a
                // nominal-mass window (`candidate_nominal_bounds` —
                // four buckets at the default 20 ppm / -1..=2), so EVERY genuine sequence
                // competitor shares the winner's nominal mass. They all took the
                // `mkey == best_mass_key` branch, which only ever raises `best_raw` and
                // never writes `second_raw` — so the emitted delta was the lead over a
                // candidate one or two Da away (a different isotope hypothesis), not over
                // the runner-up peptide. DeltaRankScore is one of Percolator's strongest
                // classical features and it reported maximum confidence on genuine ties.
                //
                // The peptidoform key still collapses the same peptide matched in multiple
                // proteins or at multiple charges, which is what the mass key was reaching
                // for. A 64-bit hash collision would merge two unrelated peptides, leaving
                // the delta unchanged or slightly overstated; the probability against the
                // single stored `best_pep_key` is negligible.
                let mut best_raw = i32::MIN;
                let mut best_pep_key: u64 = 0;
                let mut have_best = false;
                let mut second_raw = i32::MIN;
                // Best RawScores among candidates from Reference entries, for
                // `DeltaRefScore`; only tracked when the database holds
                // isoform, variant or chain entries.
                let mut ref_tracker = RefTracker::default();
                let track_ref = self.idx.has_alternative_entries;

                // Tailor per-spectrum score calibration (additive PIN feature):
                // accumulate a histogram of EVERY scored candidate's rounded RawScore
                // for this spectrum. RawScore is a bounded integer so the map is
                // small and cheap; it never influences which candidates are scored
                // or kept (purely a side-channel for the post-loop denominator).
                // Each spectrum's closure owns its own histogram → parallel-safe.
                let mut tailor_hist: FxHashMap<i32, u32> = FxHashMap::default();
                let mut tailor_total: u32 = 0;
                // S4 null pool: every scored candidate's pin RawScore on this spectrum.
                let mut strong_null_stats = OnlineStats::default();

                // Determine which charge states to try for this spectrum.
                // For charge-explicit spectra this is a single entry; for charge-missing,
                // typically 2-3 entries (small overhead, correct behavior).
                let charges_to_try: SmallVec<[u8; 4]> = charges_to_try(spec, params);

                // Build (and cache) a ScoredSpectrum per charge to evaluate.
                //
                // A single ScoredSpectrum keyed off `spec.precursor_charge.unwrap_or(2)`
                // would force charge-missing spectra to use z=2 even when evaluating
                // z=3 candidates — wrong precursor filtering, wrong partition, wrong
                // main_ion.
                //
                // For charge-explicit spectra the cache has exactly 1 entry (no overhead).
                // For charge-missing spectra, typically 2-3 entries per spectrum.
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

                // This spectrum's candidates, retrieved from the out-of-core index:
                // the fragment-index shortlist when that retrieval is on, else every
                // peptidoform in its nominal precursor windows, in database
                // enumeration order. The precursor gate is still
                // `matches_precursor` inside the loop; retrieval only decides which
                // candidates are looked at. PSM `candidate_idxs` start as LOCAL slot
                // indices into this list and are remapped to global indices once, at
                // the end of the closure.
                let shift_ppm = params.precursor_mass_shift_ppm;
                let mut mmap_cands_owned: Vec<Candidate> = if let Some(fi) = &frag_index {
                    let sel = fi.query(
                        spec,
                        &charges_to_try,
                        params,
                        scorer.feature_match_tolerance(),
                        params.fragment_index_top_k as usize,
                        crate::search_params::FRAGMENT_INDEX_MIN_MATCHED,
                        crate::fragment_index::vote_peaks_for(scorer.feature_match_tolerance()),
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

                // Hoist the loop-invariant cleavage-credit constants out of the
                // per-candidate hot path: resolving them once here avoids
                // re-invoking four small `aa_set` accessor methods (each a
                // HashMap field deref) for every candidate.
                //
                // The four credit/penalty values are SearchParams-constant; we
                // resolve them ONCE here. The per-candidate logic becomes four
                // branches over precomputed i32 constants.
                let enz_credit_neighboring = aa_set_for_scoring.neighboring_aa_cleavage_credit();
                let enz_penalty_neighboring = aa_set_for_scoring.neighboring_aa_cleavage_penalty();
                let enz_credit_peptide = aa_set_for_scoring.peptide_cleavage_credit();
                let enz_penalty_peptide = aa_set_for_scoring.peptide_cleavage_penalty();
                let enz_is_c_term = params.enzyme.is_c_term();
                let enz_is_n_term = params.enzyme.is_n_term();
                let enz = params.enzyme;
                let edge_enabled = edge_scoring_enabled(scorer);

                // Per-charge queue keyed by charge state. Retains top-N PSMs
                // independently per precursor charge (Kim et al., Nat Commun 5:5277, 2014).
                let mut per_charge_queues: FxHashMap<u8, TopNQueue> = FxHashMap::default();

                // Cascade Pass 1 is a NARROW brute-force window scan: iterate every
                // candidate whose nominal mass falls in the precursor window.
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
                    // Conservative per-peptide bound on the cumulative
                    // edge_score for two-stage gating. `psm_edge_score` returns
                    // `sum of n-1 per-edge scores`, each clamped to roughly [-4, +4]
                    // (log probability ratios). 10 per edge is a very loose upper
                    // bound; we only need it to never UNDER-estimate the max so
                    // we don't skip a candidate that could win.
                    let max_edge_bonus_per_edge: f32 = 10.0;
                    let n_minus_1 = cand.peptide.length().saturating_sub(1) as f32;
                    let max_edge_bonus = max_edge_bonus_per_edge * n_minus_1;
                    for &z in &charges_to_try {
                        let scored_spec = scored_spec_for_charge(z);
                        // Track (pin_score, edge, rank_score) for the
                        // best isotope offset. `pin_score` (= node + cleavage)
                        // remains the PIN RawScore distribution Percolator
                        // was trained on. `rank_score` (= node + cleavage + edge)
                        // is the queue-ordering key (node + cleavage + edge;
                        // Kim et al., Nat Commun 5:5277, 2014).
                        //
                        // `score_psm` and `psm_edge_score` are BOTH
                        // iso-offset independent (they take `(scored_spec,
                        // peptide, scorer, charge)` — no iso parameter), so the
                        // iso loop only finds which offsets match (cheap
                        // precursor-mass check), then we compute
                        // pin_score + edge_score ONCE.
                        //
                        // Two-stage gate: if `pin_score + max_edge_bonus` can't
                        // exceed the queue's worst retained rank_score, skip the
                        // edge_score call entirely. For top-N=1 (Astral) this
                        // gates ~99% of candidates after the queue fills.
                        let mut iso_errs: SmallVec<[MassError; 4]> = SmallVec::new();
                        // Cascade Pass 1 is NARROW: the precursor match uses the tight
                        // `matches_precursor` path against the selected m/z, not the
                        // wide isolation window.
                        for offset in
                            *params.isotope_error_range.start()..=*params.isotope_error_range.end()
                        {
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
                        if iso_errs.is_empty() {
                            continue;
                        }

                        // Compute pin_score ONCE (iso-independent).
                        let pin_score =
                            score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da)
                                + cleavage_credit;

                        // DeltaRawScore capture: fold this candidate's RawScore into
                        // the per-spectrum best / second-best-distinct-peptide
                        // trackers BEFORE the `could_win` gate, so a strong runner-up
                        // that can't beat the top-1 still contributes to the delta.
                        {
                            let s = pin_score.round() as i32;
                            if track_ref
                                && self.idx.db.proteins[cand.protein_index].entry_kind()
                                    == model::EntryKind::Reference
                            {
                                ref_tracker.fold(s, peptidoform_key(&cand.peptide));
                            }
                            // Tailor histogram: fold this candidate's rounded RawScore
                            // in (same `s` the DeltaRawScore tracker uses). One u32
                            // increment, no per-candidate allocation.
                            *tailor_hist.entry(s).or_insert(0) += 1;
                            tailor_total += 1;
                            strong_null_stats.push(pin_score);
                            // A score at or below the runner-up changes neither
                            // tracker (the runner-up never exceeds the best), so
                            // the peptidoform hash is only taken above it.
                            if s > second_raw {
                                let mkey = peptidoform_key(&cand.peptide);
                                if have_best && mkey == best_pep_key {
                                    // Same peptidoform (another protein / another
                                    // charge): not a competitor, so it may only raise
                                    // the best.
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

                        // Gate against the queue's current worst rank_score
                        // before invoking edge_score.
                        let could_win = match per_charge_queues.get(&z) {
                            // With edges disabled the bound is 0, so the gate must admit
                            // a candidate tied with the worst retained one: the queue keeps
                            // ties.
                            Some(q) if q.len() >= q.capacity() as usize => {
                                q.worst_rank_score().is_none_or(|worst| {
                                    if edge_enabled {
                                        pin_score + max_edge_bonus > worst
                                    } else {
                                        pin_score >= worst
                                    }
                                })
                            }
                            // Queue below capacity (or doesn't exist yet): accept
                            // everything until it fills up.
                            _ => true,
                        };
                        if !could_win {
                            continue;
                        }

                        // Stage 2: compute edge_score ONCE (also iso-independent).
                        let edge_i = psm_edge_score(scored_spec, &cand.peptide, scorer, z);
                        let rank_score = pin_score + edge_i as f32;

                        // Pick the iso-offset with the smallest |mass_error_ppm|
                        // for the PIN row (preserves the tie-break where
                        // the first-matched iso wins when scores are equal). Since
                        // score is iso-independent, the iso choice only affects
                        // the pin `isotope_error` / `dm` columns.
                        let err = iso_errs
                            .into_iter()
                            .min_by(|a, b| {
                                a.mass_error_ppm
                                    .abs()
                                    .partial_cmp(&b.mass_error_ppm.abs())
                                    .unwrap_or(std::cmp::Ordering::Equal)
                            })
                            .unwrap();

                        let features = PsmFeatures::default();
                        let psm = PsmMatch {
                            spectrum_idx: spec_idx,
                            candidate_idxs: vec![cand_idx],
                            charge_used: z,
                            mass_error_ppm: err.mass_error_ppm,
                            score: pin_score,
                            rank_score,
                            edge_score: edge_i,
                            activation_method: Some(scorer.param().data_type.activation),
                            features,
                            isotope_offset: err.isotope_offset,
                            precursor_mz_override: None,
                            site_evidence: None,
                            tie_key: crate::psm::tie_key_of(cand),
                        };
                        per_charge_queues
                            .entry(z)
                            .or_insert_with(|| TopNQueue::new(retention_cap))
                            .push(psm);
                    }
                }

                // pepSeq + score dedup per-charge BEFORE GF compute.
                // Same peptide matched against multiple proteins collapses to one
                // PsmMatch with aggregated candidate_idxs (peptide-sequence dedup).
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

                // Candidate selection/ranking already happened on `rank_score`
                // (RawScore = node + cleavage + edge). The generating function has
                // been removed, so there is no per-charge SpecEValue DP — the
                // per-charge queues are merged directly into the spectrum queue.
                // Spectrum-level merge. The TopNQueue::push (Ordering::Equal arm)
                // keeps rank_score ties at capacity, matching Kim et al. (Nat Commun
                // 5:5277, 2014) tie-retention on the spectrum-level merge. Pushed
                // best first, so the queue keeps the top N plus the ties at the N-th
                // whatever order the charge states were searched in.
                let mut merged: Vec<PsmMatch> = per_charge_queues
                    .drain()
                    .flat_map(|(_charge, mut per_charge)| per_charge.drain_into_vec())
                    .collect();
                merged.sort_by(|a, b| b.cmp(a));
                for psm in merged {
                    queue.push(psm);
                }

                // Feature extraction (unchanged from baseline): post-merge, after
                // the per-spectrum queue is final.
                //
                // Pre-computed `psm.edge_score` from the candidate loop
                // is moved into `features.edge_score` to avoid the per-PSM
                // recomputation that `compute_psm_features` would otherwise do.
                // Spectrum-level DeltaRawScore: positive lead of the best RawScore
                // over the second-best distinct peptide. 0.0 when no distinct
                // runner-up was scored (uncontested top-1, e.g. top_n=1 with a
                // single candidate). The PIN writer emits it on the rank-1 row only.
                let delta_raw_score = if second_raw > i32::MIN {
                    (best_raw - second_raw) as f32
                } else {
                    0.0
                };

                // Tailor denominator: the RawScore at the top-1% quantile of this
                // spectrum's candidate-score distribution. Falls back to 1.0 (no
                // calibration) for small-N spectra or a non-positive quantile.
                // Computed ONCE per spectrum from the histogram (purely additive).
                let tailor_denom =
                    crate::psm::tailor_denominator(&tailor_hist, tailor_total) as f32;

                // S2 listwise terms over the retained top-K queue (not the full scored pool).
                let mut retained_scores: Vec<f32> = queue.iter_psms().map(|p| p.score).collect();
                // total_cmp (total order), not partial_cmp+unwrap_or — the latter can
                // sort incorrectly under release optimisation on the pinned toolchain.
                retained_scores.sort_by(|a, b| b.total_cmp(a));
                let listwise_score_gap =
                    scoring_crate::scoring::listwise_score_gap(&retained_scores);
                let candidate_rank_entropy =
                    scoring_crate::scoring::candidate_rank_entropy(&retained_scores);

                // Strong mode retains a wider pool, ranks it by the strong score and
                // trims it to the user's top-N. The strong score needs only the cheap
                // features, so the pool gets a light pass and only the survivors get
                // the full one (rich-ion and fragment-LLR models, site features).
                let strong_mode = params.score_mode == ScoreMode::Strong;
                let make_features = |psm: &PsmMatch, full: bool| -> PsmFeatures {
                    let ss = scored_spec_for_charge(psm.charge_used);
                    let cand = &cand_slice[psm.primary_candidate_idx() as usize];
                    let mut features = compute_psm_features_with_edge(
                        ss,
                        &cand.peptide,
                        scorer,
                        psm.charge_used,
                        Some(psm.edge_score),
                        full,
                    );

                    // The PIN `precursor_isotope_kl` / `precursor_snr` columns are
                    // left at their 0.0 defaults here. They were a per-PSM MS1
                    // isotope-envelope feature of the old single-pass chimeric mode;
                    // the two-pass cascade does not consume them (computing them ran
                    // ~121k times on Astral for ~2:40 of wall). The cascade uses MS1
                    // only in Pass 2's `run_pass2_coisolation` / `detect_coisolated`
                    // co-isolation detection — `ms1_link` and its loading stay.

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
                    // Tailor calibrated score: this PSM's RawScore divided by the
                    // spectrum's top-1% quantile RawScore. Additive — does not
                    // affect ranking. Uses `psm.score` (the PIN RawScore =
                    // node + cleavage), the same value emitted in the RawScore column.
                    features.tailor_score = psm.score / tailor_denom;
                    // Float-precision rank score (additive; ranking unchanged): the
                    // emitted RawScore (`psm.score` = node + cleavage) with the
                    // per-split node sum left UNROUNDED, recovering the discrimination
                    // integer rounding compresses on short low-evidence peptides.
                    // `psm.score - score_psm(...)` is the cleavage credit (everything
                    // RawScore adds beyond the rounded node sum); swapping the rounded
                    // node sum for the unrounded `score_psm_float` yields the float
                    // RawScore. Computed only for emitted PSMs. (b/y-only: for the rare
                    // loss-bearing peptide it omits RawScore's neutral-loss term.)
                    // The credit is an i32 and the rounded node sum is integer-valued, so
                    // `psm.score = node + credit` is exact in f32 and the credit is
                    // recoverable from the candidate alone; re-walking `score_psm` here
                    // only to subtract it out again was one of the two full cache walks
                    // per emitted PSM (competitive plan item 4.3). `cleavage_credit_for`
                    // resolves the enzyme constants and calls the same `enzyme_cleavage_credit`.
                    let credit = cleavage_credit_for(cand, enz, aa_set_for_scoring) as f32;
                    features.rank_score_float = score_psm_float(
                        ss,
                        &cand.peptide,
                        scorer,
                        psm.charge_used,
                        fragment_tolerance_da,
                    ) + credit;
                    features.strong_score = fuse_strong_score(&StrongScoreInputs {
                        intensity_signal: features.intensity_signal,
                        chance_match_surprise: features.chance_match_surprise,
                        mass_competition_evidence: features.mass_competition_evidence,
                        candidate_rank_entropy: features.candidate_rank_entropy,
                        listwise_score_gap: features.listwise_score_gap,
                    });
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
                    // Emit only user top-N rows (gate uses top_n=1); retention pool was wider.
                    queue.trim_to_capacity(params.top_n_psms_per_spectrum);
                }
                if strong_mode {
                    // Full features for the emitted rows. The strong score depends only
                    // on features both passes compute identically; the calibrated score
                    // was taken over the whole pool, so keep it.
                    queue.update_in_place(|psm| {
                        let strong = psm.features.strong_score;
                        let strong_cal = psm.features.strong_score_cal;
                        psm.features = make_features(psm, true);
                        debug_assert_eq!(psm.features.strong_score.to_bits(), strong.to_bits());
                        psm.features.strong_score_cal = strong_cal;
                    });
                }

                // Site evidence while the peaks are at hand: which modification
                // each row's residues carry, judged against the spectrum's other rows.
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

                // The surviving PSMs carry LOCAL `candidate_idxs` into this
                // spectrum's `mmap_cands`. Intern each referenced candidate into the
                // shared accumulator (de-duped by peptidoform identity, so the same
                // candidate seen from another spectrum reuses one global index) and
                // rewrite `candidate_idxs` to the global indices the PIN/TSV writers
                // resolve against `self.candidates` after the scan.
                {
                    let mut accum = self.mmap_accum.lock().unwrap();
                    queue.remap_candidate_idxs(|local| {
                        accum.intern(mmap_cands[local as usize].clone())
                    });
                }

                queue
            })
            .collect();

        queues
    }
}

/// Precursor-mass span (Da) of one chunk in [`match_spectra_out_of_core`].
pub const OUT_OF_CORE_SLICE_DA: f64 = 80.0;

/// Single-shot search with fragment-ion retrieval, the path `andes search`
/// takes. Spectra are scored in precursor-mass-ordered
/// chunks of at most [`OUT_OF_CORE_SLICE_DA`], so the per-chunk index stays
/// bounded however many spectra are passed.
///
/// Returns one queue per input spectrum, in input order (`spectrum_idx` is the
/// input position), and the candidates the queues reference. The candidate
/// index is cached under the system temp directory and reused by later calls
/// with the same database and parameters, or built in memory when the cache
/// cannot be written.
pub fn match_spectra_out_of_core(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    scorer: &RankScorer,
    fragment_tolerance_da: f64,
    decoy_prefix: &str,
) -> std::io::Result<(Vec<TopNQueue>, Vec<Candidate>)> {
    let mut params = params.clone();
    params.fragment_index_top_k = crate::search_params::FRAGMENT_INDEX_TOP_K;
    let cache = crate::candidate_index::index_cache_path(idx, &params);
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

/// Isolation window `[lo, hi]` (m/z) Pass 2 searches for co-isolated precursors
/// of `spec`, and the m/z tolerance used to match their MS1 isotope peaks. Uses
/// the per-scan isolation offsets when the parser recorded them, else the
/// configured chimeric half-width.
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

/// MS1 m/z range Pass 2 reads for `spec`: every peak co-isolation detection can
/// consult (a monoisotopic candidate inside the isolation window and its first
/// three isotope peaks at charge >= 1, each within the match tolerance), plus a
/// 0.01 Da margin. Restricting a linked MS1 scan to this range leaves
/// [`run_pass2_coisolation`]'s detections unchanged.
pub fn coisolation_ms1_bounds(spec: &Spectrum, params: &SearchParams) -> (f64, f64) {
    const MARGIN_DA: f64 = 0.01;
    let (lo, hi, tol) = coisolation_window(spec, params);
    (
        lo - tol - MARGIN_DA,
        hi + 3.0 * model::mass::ISOTOPE + tol + MARGIN_DA,
    )
}

/// One scan's Pass-2 input: its primary peptide and co-isolated precursors.
struct Pass2Scan {
    primary: Peptide,
    coisolated: Vec<crate::coisolation::CoIsolated>,
}

/// Pass 2 of the chimeric cascade. For each non-empty top-N queue (filled by
/// Pass 1 with its PRIMARY peptide), detects MS1 co-isolated precursors in the
/// isolation window, strips the primary's matched peaks, and runs a targeted
/// second-peptide GF search at each co-isolated mass; any secondary PSM is pushed
/// into the same queue as an extra PIN row.
///
/// Off-path bit-identity: no-op when `params.chimeric` is false.
///
/// `spectra`, `queues`, and `link` are a single CHUNK (the streaming chimeric
/// path scores one chunk at a time): `spectra[i]` ↔ `queues[i]` ↔
/// `link.ms2_to_ms1[i]`, all chunk-local. `offset` is the chunk's global start
/// index, added to each emitted secondary's `spectrum_idx` so it aligns with the
/// accumulated `all_spectra`. Peaks must still be present (the residual needs
/// them) — call BEFORE peaks are dropped.
///
/// Candidates come from the out-of-core index; see [`mmap_secondaries`].
pub fn run_pass2_coisolation(
    prepared: &PreparedSearch,
    spectra: &[Spectrum],
    queues: &mut [TopNQueue],
    params: &SearchParams,
    link: &Ms1Link,
    offset: usize,
) {
    // Bit-identical guard: no chimeric mode → no-op.
    if !params.chimeric {
        return;
    }

    // (A) Co-isolated precursors and the primary peptide of every scan Pass 1
    // identified.
    let scans: Vec<Option<Pass2Scan>> = queues
        .par_iter()
        .enumerate()
        .map(|(spec_idx, q)| {
            if q.is_empty() {
                return None;
            }
            let spec = spectra.get(spec_idx)?;
            // Linked MS1 scan for this MS2 (most-recent preceding MS1).
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
                // max_kl + max_n: the chimeric-N lever. Raising max_coisolated
                // searches deeper co-fragments on the residual; the KL gate keeps
                // secondaries clean.
                params.chimeric_max_kl,
                params.chimeric_max_coisolated,
            );
            if coisolated.is_empty() {
                return None;
            }
            // Primary peptide = the queue's best PSM (smallest SpecEValue).
            let primary = prepared.candidate_peptide(q.peek_top()?.primary_candidate_idx());
            Some(Pass2Scan {
                primary,
                coisolated,
            })
        })
        .collect();

    // (B) Secondaries per scan, each with its winning candidate (its
    // `candidate_idxs` index a pass-local pool).
    let found = mmap_secondaries(prepared, spectra, &scans, params);

    // (C) Emit. Winners are interned into the accumulator (scan order) and
    // their indices rewritten to the global ones. A winner whose peptidoform,
    // charge and target/decoy label the scan already emits (a retained Pass-1
    // candidate, typically at another isotope offset, or an earlier secondary)
    // is dropped: its distinct ExpMass would put it in its own competition
    // group and count the identification twice.
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
            psm.candidate_idxs = vec![accum.intern(cand)];
            // `spec_idx` is chunk-local (indexes `spectra` + `link`); the emitted
            // spectrum_idx must be global to align with `all_spectra`.
            psm.spectrum_idx = offset + spec_idx;
            // Distinct co-isolated peptide — a legitimate EXTRA emission, not a
            // competitor for the primary's slot. force_push skips capacity-based
            // eviction (plain `push` at capacity 1 would evict the primary or
            // drop the secondary).
            q.force_push(psm);
        }
    }
}

/// Widest nominal-mass span (Da) of one Pass-2 candidate pool. Bounds the
/// pool's memory: a pool holds every peptidoform of its secondaries' windows.
const PASS2_POOL_SPAN_NOMINAL: i32 = 32;

/// Pass 2: the winning secondary of every co-isolated precursor, per scan in
/// detection order, each with its candidate.
///
/// The k-th co-isolated precursors of all scans form round k (a scan's
/// secondaries compete for residual peaks, so round k needs the peaks claimed
/// in rounds `< k`). Within a round, secondaries are sorted by nominal window
/// and grouped into batches spanning at most [`PASS2_POOL_SPAN_NOMINAL`]; each
/// batch expands the union of its windows from the out-of-core index ONCE
/// (`expand_mmap_window_candidates`: enumeration multiset, order and
/// collision-decoy labels) and scans it through a pool-local bucket index. A
/// pool's order is the database enumeration order restricted to the pool, so
/// `search_secondary`'s smallest-index tie-break picks the same winner whichever
/// batch a secondary lands in.
fn mmap_secondaries(
    prepared: &PreparedSearch,
    spectra: &[Spectrum],
    scans: &[Option<Pass2Scan>],
    params: &SearchParams,
) -> Vec<Vec<(PsmMatch, Candidate)>> {
    use crate::coisolation::{search_secondary, secondary_nominal_window, CoIsolated};
    type Hit = (PsmMatch, Candidate, std::collections::HashSet<i64>);
    // (scan index, co-isolated precursor, nominal window)
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
                // Union of the batch's windows as disjoint intervals.
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
                            prepared.scorer,
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

/// Per-candidate enzymatic cleavage credit: the RawScore term for whether the
/// peptide's two termini sit at cleavable positions. `cand.pre`/`post` are the
/// flanking residues; protein termini always earn the neighboring credit.
///
/// Credit/penalty constants come from the ENZYME-REGISTERED `AminoAcidSet`
/// (the one populated by `register_enzyme`), never from an unregistered set.
///
/// `fn` + `#[inline(always)]` so LLVM inlines it into the candidate loop of
/// `run_chunk_inner`; a closure went through `FnMut::call_mut` dispatch.
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
        // N-term cleavage (neighboring)
        score += if cand.is_protein_n_term || enz.is_cleavable(pre) {
            credit_neighboring
        } else {
            penalty_neighboring
        };
        // C-term cleavage (peptide). Inline residues.last() to avoid
        // the Option::map call_mut dispatch that perf flagged.
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
        // N-term cleavage (peptide)
        score += if enz.is_cleavable(pre) {
            credit_peptide
        } else {
            penalty_peptide
        };
        // C-term cleavage (neighboring)
        score += if cand.is_protein_c_term || enz.is_cleavable(post) {
            credit_neighboring
        } else {
            penalty_neighboring
        };
    }
    score
}

/// [`enzyme_cleavage_credit`] with the constants resolved from the enzyme-registered
/// `aa_set`. `search_secondary` calls this so secondaries share the production
/// RawScore scale (`score_psm(...) + credit`).
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

/// Compute fragment-ion feature columns for a single PSM.
///
/// Uses charge-1 b/y ions only (the `NumMatchedMainIons` convention).
/// A peptide position counts at most once per ion series;
/// a position can contribute 1 from b AND 1 from y (so the maximum
/// `num_matched_main_ions` is `2 * (n - 1)` for a peptide of length n).
///
/// Returns `PsmFeatures::default()` for peptides shorter than 2 residues
/// (no cleavable fragment ions exist).
///
/// # Ion-current + error-stat features
///
/// The ion-current and error-stat PIN columns:
/// - Ion-current ratios use raw peak intensities vs total MS2 ion current.
/// - `MS2IonCurrent` is the raw sum (NOT log10); the PIN emitter emits it as-is.
/// - `IsolationWindowEfficiency` is always 0.0 (no isolation-window data
///   in the Spectrum object).
/// - Top-7 error stats: errors are collected for all matched b+y ions,
///   sorted descending by intensity, top-7 taken; absolute Da error for
///   mean/stdev, signed ppm for rel-mean/rel-stdev. Population stdev
///   formula: `sqrt(E[x²] - mean²)`.
pub(crate) fn compute_psm_features(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
) -> PsmFeatures {
    compute_psm_features_with_edge(scored_spec, peptide, scorer, charge, None, true)
}

/// As [`compute_psm_features`], reusing an `edge_score` the caller already has
/// (`PsmMatch::edge_score`, computed in the candidate loop for the same spectrum,
/// peptide and charge) instead of recomputing it.
pub(crate) fn compute_psm_features_with_edge(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    known_edge_score: Option<i32>,
    full: bool,
) -> PsmFeatures {
    let n = peptide.length();
    if n < 2 {
        return PsmFeatures::default();
    }

    // ADDITIVE edge-score feature (new PIN column). Computed
    // here so it shares the per-PSM ScoredSpectrum + scorer references that
    // the existing feature-extraction code already has on hand.
    let edge_score =
        known_edge_score.unwrap_or_else(|| psm_edge_score(scored_spec, peptide, scorer, charge));

    // Predict charge-1 b/y ions; one bool per fragment position.
    //
    // Stack-allocate b/y_matched on a 64-slot SmallVec (max
    // peptide length is 40 → n-1 ≤ 39). The prior `vec![false; n-1]` heap
    // allocations fired ~150k × 4 / PSM batch and were a measurable hot-path
    // cost. SmallVec inlines for n ≤ 64.
    // Activation-gated neutral-loss prediction. ETD preserves labile mods
    // (no loss ions); collisional/photon activation predicts them. Inert
    // unless a residue's mod declares `neutral_losses`, so output is
    // byte-identical to the canonical b/y set for every non-loss peptide.
    let predict_losses = scorer
        .param()
        .data_type
        .activation
        .predicts_neutral_losses();
    let predicted = predict_by_ions_with_losses(peptide, 1..=1, predict_losses);
    let mut b_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut y_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    // Per-position charge-1 INTACT ion intensity-ranks (u32::MAX = unmatched),
    // backing the additive complementary-ion-balance feature below.
    let mut b_rank: SmallVec<[u32; 64]> = smallvec![u32::MAX; n - 1];
    let mut y_rank: SmallVec<[u32; 64]> = smallvec![u32::MAX; n - 1];

    // Collect matched-ion details for ion-current ratio and error-stat features.
    // Each entry: (intensity, observed_mz, predicted_mz, is_b_ion).
    // SmallVec inlines for up to ~96 matched ions (b+y at n positions, with
    // some headroom for partition multi-ion-type matches at long peptides).
    let mut matched_ions: SmallVec<[(f32, f64, f64, bool); 96]> = SmallVec::new();
    let mut matched_ion_ranks: SmallVec<[u32; 96]> = SmallVec::new();

    // Percolator feature counting uses a fragment tolerance of its own, NOT
    // param.mme (Kim et al., Nat Commun 5:5277, 2014): on high-res data the
    // tolerance derived from the run's measured fragment error (20 ppm before
    // a measurement), 0.5 Da on low-res data. The param.mme value (0.5 Da for
    // the hcd_qexactive_tryp model) is the coarser binning tolerance used by
    // the rank-distribution tables — appropriate for node-score lookup but
    // ~50× too wide for feature counting at m/z 500.
    let feat_tol = scorer.pin_match_tolerance();
    let feature_tol_is_ppm = matches!(feat_tol, model::tolerance::Tolerance::Ppm(_));
    let feature_tol = feat_tol.raw_value(); // numeric value in the unit (ppm or Da)

    // The trained fragment-intensity and rich-ion models (and the mod-site
    // features beside them) match at the window they were trained at: the
    // class default, not the measured tolerance. Serving them at the measured
    // window lost PSMs on Astral and phospho in an A/B against this split.
    let model_tol = scorer.feature_match_tolerance();
    let model_tol_is_ppm = matches!(model_tol, model::tolerance::Tolerance::Ppm(_));
    let model_tol = model_tol.raw_value();

    // Neutral-loss masses (charge-1 fragment partners at −H2O / −NH3).
    const H2O_LOSS_DA: f64 = 18.0106;
    const NH3_LOSS_DA: f64 = 17.0265;
    let mut neutral_loss_ion_count: u32 = 0;

    for p in &predicted {
        let tol_da = if feature_tol_is_ppm {
            p.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        if let Some((rank, intensity, peak_mz)) = scored_spec.nearest_peak_full(p.mz, tol_da) {
            // c-ions are N-terminal like b; z•-ions are C-terminal like y. NOTE: this
            // b/y feature path only ever sees `predict_by_ions` output (B/Y), so the
            // C/Z arms here (and in the sibling b/y feature fns) are currently
            // unreachable — kept for exhaustiveness and to stay correct if a future
            // caller routes c/z through these features. ETD c/z scoring today uses the
            // dedicated `cz_hyperscore_psm` path, not this one.
            let is_b = matches!(p.kind, IonKind::B | IonKind::C);
            matched_ions.push((intensity, peak_mz, p.mz, is_b));
            matched_ion_ranks.push(rank);

            // position is 1-based (b1/y1 = index 0 in the matched arrays)
            let pos = (p.position - 1) as usize;
            match p.kind {
                IonKind::B | IonKind::C => {
                    if pos < b_matched.len() {
                        b_matched[pos] = true;
                        // Record the best (lowest) INTACT-ion rank per position
                        // for the complementary-balance feature; skip loss ions.
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

            // Matched ion with a neutral-loss partner at −H2O or −NH3 (charge 1).
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

    // NumMatchedMainIons: each (bond, direction) tuple contributes 1 if at
    // least one charge-1 prefix/suffix ion matched (Kim et al., Nat Commun
    // 5:5277, 2014 Percolator features).
    let num_matched: u32 = (b_matched.iter().filter(|&&m| m).count()
        + y_matched.iter().filter(|&&m| m).count()) as u32;

    // ── Strong-score Stage-1: charge-2 matched-ion count ─────────────────────
    // High-charge precursors fragment into +2 ions ignored by the charge-1
    // loop above. Separate probe — does not alter existing feature values.
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

    // Ion-current ratio features (partition-ion-list).
    //
    // Kim et al. (Nat Commun 5:5277, 2014): explained ion current sums
    // matched peak intensities across the full partition ion list (b, y, and
    // partition-specific variants). The b/y-only buffer above under-counts
    // intensity for Percolator features; we replace it with a partition-wide
    // sum and use partition-wide matches for longest_b/y. NumMatchedMainIons
    // continues to count charge-1 b/y matches.
    let parent_mass = scored_spec.parent_mass();
    let num_segments = scorer.param().num_segments.max(1) as usize;

    // Stack-allocate (same rationale as b/y_matched above).
    let mut b_any_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut y_any_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut sum_prefix_intensity: f64 = 0.0;
    let mut sum_suffix_intensity: f64 = 0.0;

    // Use accurate residue mass for theoretical m/z (Kim et al., Nat Commun
    // 5:5277, 2014). IonType::mz divides nominal mass by INTEGER_MASS_SCALER
    // (0.999497) to recover an approximate accurate mass — that
    // approximation can drift ~0.014 Da from the true accurate mass per
    // residue (NEEQSR's N: nominal 114 → 114.057 vs accurate 114.043),
    // which is way outside the 20 ppm feature-matching window for high-res
    // instruments. We bypass that conversion by computing theo_mz directly
    // from accurate residue mass + ion offset.
    let mut prm_accurate: f64 = 0.0;
    let mut srm_accurate: f64 = 0.0;

    // Cache the per-segment ion list ONCE per spectrum (constant
    // for fixed `(charge, parent_mass)`), avoiding the `partition_for` binary
    // search + HashMap lookup that fired for every (split × segment) pair.
    // On Astral with ~150k PSMs × ~12 splits × 2 segments = ~3.6M lookups
    // saved per run. SmallVec<[&[IonType]; 8]> inlines (num_segments is
    // typically 1-2; clamp at 8 to be safe).
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

        // Each segment's ion list is checked separately so per-bond ion sums
        // align with the partition model (Kim et al., Nat Commun 5:5277, 2014).
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

    let total_intensity = scored_spec.total_intensity(); // raw sum, all peaks
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
    // MS2 ion current is the raw sum (no log10 transform).
    let ms2_ion_current = if total_intensity > 0.0 {
        total_intensity as f32
    } else {
        0.0
    };
    // Isolation-window efficiency (the precursor's share of the ion current in
    // its isolation window) needs the MS1 scan around each precursor. The
    // standard search streams MS2 only (MS1 is read under `--chimeric` and
    // `--precursor-mono`), so the column is a constant 0.0 for every PSM and
    // carries no information to the rescorer.
    let isolation_window_efficiency = 0.0_f32;

    // ── Top-7 mass-error statistics ───────────────────────────────────────────

    // Sort matched ions descending by intensity.
    matched_ions.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let top7 = &matched_ions[..matched_ions.len().min(7)];

    // All four *ErrorTop7 columns are in ppm (Kim et al., Nat Commun 5:5277,
    // 2014 Percolator features). MeanErrorTop7 = mean of |ppm error|;
    // MeanRelErrorTop7 = mean of signed ppm error. The "Rel" suffix
    // distinguishes signed vs absolute, not Da-vs-ppm.
    //
    // Population stdev formula: sqrt(sum_sq/n - mean²).
    let abs_ppm_errors: Vec<f64> = top7
        .iter()
        .filter(|&&(_, _, pred, _)| pred > 0.0)
        .map(|&(_, obs, pred, _)| ((obs - pred) / pred * 1e6).abs())
        .collect();
    let rel_ppm_errors: Vec<f64> = top7
        .iter()
        .filter(|&&(_, _, pred, _)| pred > 0.0)
        .map(|&(_, obs, pred, _)| (obs - pred) / pred * 1e6)
        .collect();

    fn mean_and_pop_stdev(values: &[f64]) -> (f32, f32) {
        if values.is_empty() {
            return (0.0, 0.0);
        }
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let sum_sq: f64 = values.iter().map(|v| v * v).sum();
        let var = (sum_sq / n - mean * mean).max(0.0); // clamp negative rounding noise
        (mean as f32, var.sqrt() as f32)
    }

    let (mean_error_top7, stdev_error_top7) = mean_and_pop_stdev(&abs_ppm_errors);
    let (mean_rel_error_top7, stdev_rel_error_top7) = mean_and_pop_stdev(&rel_ppm_errors);

    // ── Strong-score Stage-1: ppm-Gaussian "tight-match evidence" ─────────────
    // Σ exp(-½ ((ppmᵢ - c)/σ)²) over ALL matched ions (not just top-7): c and σ
    // are the run's measured fragment offset and 1.5 × its spread on high-res
    // data, 0 and 7 ppm otherwise. Each ion contributes ~1 when matched near
    // the run's error centre and ~0 near the tolerance edge, so the sum rewards
    // PSMs whose fragments cluster tightly in mass accuracy — the high-res
    // signal the intensity-rank model throws away. For low-res spectra every
    // ppm error is large, so the term is ~0 for all candidates and Percolator
    // simply learns a near-zero weight (harmless).
    let (ppm_centre, ppm_sigma) = scorer.ppm_gaussian();
    let ppm_gaussian_score: f32 = matched_ions
        .iter()
        .filter(|&&(_, _, pred, _)| pred > 0.0)
        .map(|&(_, obs, pred, _)| {
            let z = ((obs - pred) / pred * 1e6 - ppm_centre) / ppm_sigma;
            (-0.5 * z * z).exp()
        })
        .sum::<f64>() as f32;

    // ── Strong-score Stage-2: per-peak chance-match "surprise" (null moat) ───
    // For each matched ion: p_chance = ρ·Δ, ρ = local peak density (peaks/Da)
    // in a ±50 m/z window, Δ = the match window (2·tol_da). surprise =
    // max(0, -ln(p_chance)): a tight match in a sparse region is improbable by
    // chance (high surprise); a wide match in a crowded region is nearly free
    // (~0). Summed over matched ions = total deterministic evidence that the
    // matches are NOT coincidental — the denominator no ML rescorer computes.
    let chance_match_surprise: f32 = matched_ions
        .iter()
        .filter(|&&(_, _, pred, _)| pred > 0.0)
        .map(|&(_, obs, pred, _)| {
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

    // ── Strong-score Stage-1: complementary ladder ─────────────────────────
    // Bond i (1-based) is reported by both b_i and y_{n-i}; track cleavage
    // sites where BOTH halves are observed. `b_any_matched`/`y_any_matched` are
    // 0-indexed by ion number, so b_i = [i-1] and y_{n-i} = [n-i-1].
    let mut complementary_sites: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    for i in 1..n {
        complementary_sites[i - 1] = b_any_matched[i - 1] && y_any_matched[n - i - 1];
    }
    let longest_complementary_ladder = longest_run(&complementary_sites);

    // ADDITIVE complementary-ion BALANCE: Σ over bonds where both the charge-1
    // intact b_i and complementary y_{n-i} matched, weighted by intensity-rank
    // agreement 1/(1+|rank_b − rank_y|). Orthogonal to the longest-RUN ladder
    // (no intensity) — captures whether the two halves of a bond appear at
    // correlated ranks (true peptide) vs lopsided one-series matches (decoy).
    let mut complementary_ion_balance = 0.0_f32;
    for i in 1..n {
        let rb = b_rank[i - 1];
        let ry = y_rank[n - i - 1];
        if rb != u32::MAX && ry != u32::MAX {
            let diff = (rb as i64 - ry as i64).unsigned_abs() as f32;
            complementary_ion_balance += 1.0 / (1.0 + diff);
        }
    }

    // ── Strong-score Stage-1: mean matched intensity rank ──────────────────
    // Average intensity-rank of matched b/y ions (1 = most intense). Lower =
    // better — real PSMs explain dominant peaks.
    let mean_matched_intensity_rank = if matched_ion_ranks.is_empty() {
        0.0
    } else {
        matched_ion_ranks.iter().map(|&r| r as f32).sum::<f32>() / matched_ion_ranks.len() as f32
    };

    // ── Strong-score Stage-2: within-peptide ambiguity + mass competition ───
    let charge1_predicted = predict_by_ions(peptide, 1..=1);
    let mut theo_mz_list: SmallVec<[f64; 128]> = SmallVec::new();
    for p in &charge1_predicted {
        theo_mz_list.push(p.mz);
        theo_mz_list.push(p.mz - H2O_LOSS_DA);
        theo_mz_list.push(p.mz - NH3_LOSS_DA);
    }
    let unique_match_fraction = if num_matched == 0 {
        0.0
    } else {
        let evidence_sum: f64 = matched_ions
            .iter()
            .map(|&(_, obs, pred, _)| {
                let tol_da = if feature_tol_is_ppm {
                    obs * feature_tol / 1e6
                } else {
                    feature_tol
                };
                let ambiguity = theo_mz_list
                    .iter()
                    .filter(|&&theo| (theo - obs).abs() <= tol_da && (theo - pred).abs() > 1e-9)
                    .count();
                1.0 / (1.0 + ambiguity as f64)
            })
            .sum();
        (evidence_sum / num_matched as f64) as f32
    };
    let mass_competition_evidence_val = mass_competition_evidence(
        scored_spec,
        &matched_ions,
        &theo_mz_list,
        feature_tol,
        feature_tol_is_ppm,
    );

    // ── Strong-score S1: intensity-model signal (additive PIN column) ───────
    // Driven by the v3 frag-intensity regressor on the param; 0.0 without one.
    let frag_intensity_model = scorer.param().frag_intensity_model.as_deref();
    // Walk the GBDT ONCE for this candidate. `intensity_signal` and
    // `frag_llr_battery` below both need the identical per-ion predictions
    // (same peptide, same charge, same model, same ion order), and each used to
    // recompute the whole 300-tree ensemble for every fragment — the single
    // hottest thing in the search, done twice for the same numbers.
    let frag_pred: Option<Vec<f64>> = frag_intensity_model
        .filter(|_| peptide.length() >= 2)
        .map(|g| predict_frag_intensities(g, peptide, charge));
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

    // The remaining features feed only the PIN, never the strong score, so a
    // light pass (`full == false`, used to rank strong-mode candidates before the
    // trim) leaves them at 0.0 and the survivors get a full pass.
    let (frag_pred_explained, frag_pred_chance_llr, frag_topk_observed, rich_ion_llr_val, mod_site) =
        if full {
            // Tier-2 frag-intensity LLR battery (additive PIN features; 0.0 when no
            // frag model). Discriminative alternative to the cosine intensity_signal.
            let (e, c, t) = frag_llr_battery(
                frag_intensity_model,
                scored_spec,
                peptide,
                charge,
                model_tol,
                model_tol_is_ppm,
                frag_pred_slice,
            );
            // Decoy-aware rich-ion LLR (additive PIN feature; 0.0 when no rich-ion model).
            let r = rich_ion_llr(
                scorer.param().rich_ion_model.as_deref(),
                scored_spec,
                peptide,
                charge,
                model_tol,
                model_tol_is_ppm,
            );
            // Mod-localization site-determining-ion features (additive PIN columns; all
            // 0.0 for an unmodified peptide). Uses the SAME tolerance as the
            // RichIonLLR / ion-feature path above (train/serve parity).
            let m = mod_site_features(peptide, scored_spec, model_tol, model_tol_is_ppm);
            (e, c, t, r, m)
        } else {
            (0.0, 0.0, 0.0, 0.0, Default::default())
        };

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
        // Chimeric precursor-envelope features default 0.0; populated only
        // under --chimeric by the feature fill.
        precursor_isotope_kl: 0.0,
        precursor_snr: 0.0,
        // Set by the caller (run_chunk_inner) from the per-spectrum capture;
        // `compute_psm_features` has no cross-candidate view, so default here.
        delta_raw_score: 0.0,
        // Set by the caller from the per-spectrum Tailor histogram + denominator;
        // `compute_psm_features` has no cross-candidate view, so default to the
        // uncalibrated RawScore-equivalent here (overwritten in fill_post_topn).
        tailor_score: 0.0,
        rank_score_float: 0.0, // set per-PSM at emission (fill_post_topn)
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
        // Engine-wide RT features: populated post-search (after per-run
        // calibration) by `output::rt_wiring::populate_rt_features`; neutral 0.0
        // here so a run without observed RT stays baseline-identical.
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
        isobaric_rt_margin: 0.0, // glyco-only; set by populate_glyco_rt_features
    }
}

// ── Unit tests for feature columns ───────────────────────────────────────────

#[cfg(test)]
mod accumulator_tests {
    use super::*;
    use model::{AminoAcidSetBuilder, Protein, ProteinDb};

    /// After the accumulator is drained, a candidate interned again gets a fresh
    /// index into the new store, not its index in the drained one.
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
            assert_eq!(accum.intern(cands[0].clone()), 0);
            assert_eq!(accum.intern(cands[1].clone()), 1);
        }
        prepared.sync_materialized_candidates();
        assert_eq!(prepared.candidates.len(), 2);
        let mut accum = prepared.mmap_accum.lock().unwrap();
        assert!(accum.index_of.is_empty());
        assert_eq!(accum.intern(cands[1].clone()), 0);
        assert_eq!(accum.candidates.len(), 1);
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
        }
    }

    /// A spectrum without a charge is searched at each charge's own window. The
    /// peaks also hold every ion of a longer peptide whose mass lies between the
    /// 3+ and 4+ windows: with one index slot over the union of all charges'
    /// windows, that out-of-window peptide out-voted the true 3+ peptide and the
    /// spectrum came back empty.
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

    /// A charge that does not fit a `u8` is unknown, not truncated (300 as u8
    /// would be 44).
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

    /// Minimal RankScorer for feature tests, with mme = Da(tol_da).
    ///
    /// Uses realistic prefix/suffix offsets so the partition-ion-list
    /// intensity-ratio path matches peaks placed at `predict_by_ions`'s
    /// standard b/y m/z values (b_neutral + PROTON; y_neutral = suffix +
    /// H2O + PROTON). The partition ion list is read directly, so the
    /// offsets matter (offset=0.0 with no suffix ion would not work here).
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

    /// Build a minimal peptide of `len` alanine residues with flanks `_-`.
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
        }
    }

    // ── Test: empty spectrum → all new features are 0 ───────────────────────

    #[test]
    fn compute_psm_features_top7_error_stats_zero_when_no_matches() {
        let pep = ala_peptide(4);
        let spec = make_spectrum(vec![]); // no peaks
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

    // ── Test: ion-current ratios populate and satisfy arithmetic invariant ───

    #[test]
    fn compute_psm_features_populates_ion_current_ratios() {
        // Use a 3-residue peptide (ALA-ALA-ALA). predict_by_ions(charge=1) gives:
        //   b1, y1, b2, y2 at definite m/z values.
        // We place spectrum peaks at exactly those m/z values so all ions match,
        // then verify explained_ratio > 0 and n + c == explained.
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);

        // Place peaks exactly at every predicted m/z with increasing intensities.
        let mut peaks: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        // Add some unmatched background intensity so total_intensity > matched.
        peaks.push((1500.0, 5.0)); // far from any ion
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = make_spectrum(peaks);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.01), 2); // tight tolerance

        // All ratios should be positive since all predicted ions match.
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

        // Invariant: n_term + c_term == explained (within float precision)
        let sum = f.n_term_ion_current_ratio + f.c_term_ion_current_ratio;
        assert!(
            (sum - f.explained_ion_current_ratio).abs() < 1e-5,
            "n_term + c_term should == explained ({} + {} != {})",
            f.n_term_ion_current_ratio,
            f.c_term_ion_current_ratio,
            f.explained_ion_current_ratio
        );

        // ms2_ion_current should equal total peak intensity sum.
        let total: f32 = ss.total_intensity() as f32;
        assert!(
            (f.ms2_ion_current - total).abs() < 1.0,
            "ms2_ion_current {} should match total spectrum intensity {}",
            f.ms2_ion_current,
            total
        );

        // isolation_window_efficiency always 0.0.
        assert_eq!(f.isolation_window_efficiency, 0.0);
    }

    // ── Test: unique-match fraction ──────────────────────────────────────────

    #[test]
    fn compute_psm_features_unique_match_fraction() {
        let scorer = make_scorer(0.01);

        // (a) Well-separated ions → fraction ~1.0.
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
        // Neutral-loss offsets in the theo list can add minor ambiguity even
        // when peaks are well separated; still expect a high fraction.
        assert!(
            f_unique.unique_match_fraction > 0.7,
            "separated ions should yield a high fraction, got {}",
            f_unique.unique_match_fraction
        );

        // (b) Two theoretical ions within tol of one peak → evidence 0.5 each.
        // Scan alanine peptides for a charge-1 ion pair whose m/z values fall
        // within the 20 ppm feature window (including neutral-loss offsets).
        const H2O: f64 = 18.0106;
        const NH3: f64 = 17.0265;
        let mut contested_case: Option<(Peptide, f64, f64)> = None;
        'search: for len in 3..=12 {
            let pep2 = ala_peptide(len);
            let pred2 = predict_by_ions(&pep2, 1..=1);
            let mut theo: Vec<f64> = Vec::new();
            for p in &pred2 {
                theo.push(p.mz);
                theo.push(p.mz - H2O);
                theo.push(p.mz - NH3);
            }
            for i in 0..theo.len() {
                let tol = theo[i] * 20e-6;
                for j in (i + 1)..theo.len() {
                    if (theo[i] - theo[j]).abs() <= 2.0 * tol {
                        contested_case = Some((pep2, theo[i], theo[j]));
                        break 'search;
                    }
                }
            }
        }
        let (pep2, peak_mz, other_mz) =
            contested_case.expect("should find a contested theo pair in alanine scan");
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
        // Both ions match one peak → each evidence 0.5 → fraction 0.5.
        assert!(
            (f_contested.unique_match_fraction - 0.5).abs() < 0.15,
            "two ions on one peak → fraction ~0.5, got {}",
            f_contested.unique_match_fraction
        );
    }

    // ── Test: doubly-charged matched-ion count ───────────────────────────────

    #[test]
    fn compute_psm_features_doubly_charged_matched_ion_count() {
        let pep = ala_peptide(4);
        let predicted_z2 = predict_by_ions(&pep, 2..=2);
        assert!(
            !predicted_z2.is_empty(),
            "charge-2 ions should exist for 4-mer"
        );

        // Empty spectrum → count 0.
        let f_empty = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(vec![])),
            &pep,
            &make_scorer(0.01),
            3,
        );
        assert_eq!(f_empty.doubly_charged_matched_ion_count, 0);

        // Peaks at charge-2 b/y m/z → count > 0.
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

    // ── Test: mean matched intensity rank ────────────────────────────────────

    #[test]
    fn compute_psm_features_mean_matched_intensity_rank() {
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);

        // (a) Ion peaks are the most intense → mean rank near 1.
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
        // Only ion peaks present → ranks 1..4 → mean 2.5.
        assert!(
            (f_top.mean_matched_intensity_rank - 2.5).abs() < 0.1,
            "ions as the only peaks → mean rank should be ~2.5, got {}",
            f_top.mean_matched_intensity_rank
        );

        // (b) Bury ion peaks under many bigger noise peaks → mean rank rises.
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

    // ── Test: longest complementary ladder ─────────────────────────────────

    #[test]
    fn compute_psm_features_longest_complementary_ladder() {
        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);
        let n = pep.length();

        // (a) All predicted ions matched → ladder == n-1.
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

        // (b) Drop the middle bond's complementary y (y3 for bond 2) → ladder
        //     splits; longest run is the longer of the two halves.
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
        // Bonds 1 and 2-4 split at bond 2: runs of 1 and 2 → longest = 2.
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

    // ── Test: complementary-ion balance (additive feature #campaign-10) ──────

    #[test]
    fn compute_psm_features_complementary_ion_balance() {
        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);

        // (a) All ions matched → bonds have both b_i and y_{n-i} → balance > 0.
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

        // (b) Only b ions matched (no y) → no complementary pair → balance == 0.
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

    // ── Test: neutral-loss ion count ─────────────────────────────────────────

    #[test]
    fn compute_psm_features_neutral_loss_ion_count() {
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);
        let b2 = predicted
            .iter()
            .find(|p| matches!(p.kind, IonKind::B) && p.position == 2)
            .expect("b2 should exist");

        // Peak at b2 and its −H2O partner → count == 1.
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

    // ── Test: strong-score Stage-1 bolt-ons (ppm-Gaussian + complementary) ──

    #[test]
    fn compute_psm_features_strong_score_stage1_bolt_ons() {
        // ALA-ALA-ALA → predict_by_ions(charge 1) yields b1,b2,y1,y2.
        // Bonds 1 and 2 are each reported by a (b, complementary-y) pair:
        //   bond 1 = b1 + y2, bond 2 = b2 + y1.
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);

        // (a) Peaks exactly at every predicted m/z (0 ppm error).
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

        // All 4 ions matched at 0 ppm → each Gaussian kernel == 1 → sum ≈ 4.
        assert!(
            (f_exact.ppm_gaussian_score - 4.0).abs() < 0.05,
            "ppm_gaussian_score should be ~4 for 4 exact matches, got {}",
            f_exact.ppm_gaussian_score
        );
        // Both bonds have b and complementary y observed → full ladder.
        assert_eq!(
            f_exact.longest_complementary_ladder, 2,
            "both cleavage sites should form a length-2 ladder, got {}",
            f_exact.longest_complementary_ladder
        );

        // (b) Same peaks shifted by ~15 ppm: still matched (20 ppm feature
        //     window) but each kernel exp(-½(15/7)²) ≈ 0.10, so the sum drops
        //     far below the exact case — tight matches are worth more evidence.
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

    // ── Test: top-7 error stats are nonzero when ions match ─────────────────

    #[test]
    fn compute_psm_features_error_stats_nonzero_when_ions_match_with_offset() {
        // Build a peptide and shift every peak by a fixed offset so errors are known.
        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);

        // 0.0005 Da offset = ~6 ppm at m/z 89 (Ala b1) — within the
        // hardcoded 20 ppm window that compute_psm_features uses for
        // high-resolution instruments (Kim et al., Nat Commun 5:5277, 2014).
        // The previous 0.01 Da offset assumed Rust used param.mme (~0.05 Da
        // in this fixture's make_scorer), but feature counting now uses
        // 20 ppm regardless of param.mme.
        let offset_da = 0.0005_f64;
        let mut peaks: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz + offset_da, (i + 1) as f32 * 10.0))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = make_spectrum(peaks);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        // make_scorer still accepts a tol arg for legacy compatibility, but
        // compute_psm_features uses the instrument-based hardcoded tolerance.
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.05), 2);

        // Mean error should be nonzero when peaks are systematically offset.
        // MeanErrorTop7 is in PPM, not Da. PPM error =
        // (Δm / mz) × 1e6 varies per-ion because mz differs across b1, y1,
        // b2, y2, … of the test peptide, so stdev is no longer ~0 (it's a
        // small but non-zero spread). Just verify mean is positive.
        assert!(
            f.mean_error_top7 > 0.0,
            "mean_error_top7 should be > 0 when peaks are systematically offset, got {}",
            f.mean_error_top7
        );
        // Stdev varies with m/z when offset is constant in Da and reported in
        // ppm. Just bound to "small" (PPM at typical fragment m/z 100-500 is
        // ~1-5 ppm for 0.0005 Da offset).
        assert!(
            f.stdev_error_top7 < 20.0,
            "stdev_error_top7 should be small (single-digit ppm) for identical-Da offset, got {}",
            f.stdev_error_top7
        );
        // Relative error should also be nonzero.
        assert!(
            f.mean_rel_error_top7 != 0.0,
            "mean_rel_error_top7 should be nonzero when peaks are offset"
        );
    }

    // ── Test: ms2_ion_current mirrors total_intensity exactly ───────────────

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

    // ── Test: PROTON mass sanity — b1 ion for alanine at charge 1 ───────────
    // This verifies the predict_by_ions formula aligns with our test setup.
    #[test]
    fn b1_mz_for_alanine_is_proton_plus_residue_mass() {
        use model::amino_acid::AminoAcid;
        let aa = AminoAcid::standard(b'A').unwrap();
        let residue_mass = aa.mass; // monoisotopic residue mass
        let expected_b1_mz = residue_mass + PROTON; // charge 1
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

/// Pre-merge dedup pass: collapse PSMs sharing the same peptide sequence
/// and rounded rank_score before per-charge merge.
pub(crate) fn dedup_pepseq_score(psms: Vec<PsmMatch>, candidates: &[Candidate]) -> Vec<PsmMatch> {
    use std::collections::btree_map::Entry;
    use std::collections::BTreeMap;

    let mut pep_key_cache: FxHashMap<u32, Arc<PepDedupKey>> = FxHashMap::default();
    let mut groups: BTreeMap<DedupMapKey, PsmMatch> = BTreeMap::new();

    // Best first in the `PsmMatch` total order, so each group's survivor and the
    // order its candidates (proteins) are listed in never depend on the order the
    // queue handed the PSMs over.
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
            // Never merge a target with a reversed-decoy that coincidentally shares
            // residues+rounded score (palindromic / reversal-coincident peptides):
            // a merged row would take a +1/-1 Label by heap order, silently
            // corrupting Percolator's target-decoy competition.
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

/// Allocation-free peptidoform identity for the per-spectrum DeltaRankScore tracker:
/// residues plus per-residue modification deltas, quantised to 1e-5 Da exactly as
/// [`PepDedupKey`] quantises them, so the two agree on what "the same peptidoform"
/// means. `PepDedupKey` allocates two `Vec`s and runs once per RETAINED PSM; this
/// runs once per SCORED CANDIDATE, so it must not allocate.
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

/// Bare residue sequence (no mods, no flanks) of a candidate's peptide.
fn bare_residues(cand: &Candidate) -> Box<[u8]> {
    cand.peptide.residues.iter().map(|aa| aa.residue).collect()
}

/// Relabel as targets the decoy candidates whose bare peptide sequence is ALSO a
/// target peptide. Such "decoys" are reversal/shuffle coincidences (or
/// palindromes / short peptides) that reproduce a real target sequence; counting
/// them as decoys over-counts the decoy space (TDC conservative) and, for
/// decoy-aware training, mislabels true-sequence ions as negatives. Relabeling to
/// target lets the per-spectrum dedup merge the collision with its target twin.
/// Bare residues only (mods don't change sequence identity).
///
/// `target_seqs` must be the GLOBAL set of all DB targets (built once from the
/// full index by [`build_target_bare_seqs`]), NOT the targets that happen to land
/// in the current precursor window.
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

/// Build the GLOBAL set of target bare-peptide sequences from the full
/// base-peptide index. For every NON-decoy record, reconstruct the residue span
/// from the protein sequence (`protein.sequence[start..start+length]`), so the
/// relabel decision is global without materializing all candidates.
///
/// Bare residues only (no mods, no flanks) → bounded memory: one `Box<[u8]>` per
/// distinct target backbone, far smaller than the full candidate list.
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
    let mut set = FxHashSet::with_capacity_and_hasher(spans.len(), Default::default());
    set.extend(spans.into_iter().map(|s| s.to_vec().into_boxed_slice()));
    set
}

/// [`build_target_bare_seqs`] restricted to the target spans with at least one
/// peptidoform that `keep` accepts.
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
    let mut set = FxHashSet::with_capacity_and_hasher(spans.len(), Default::default());
    set.extend(spans.into_iter().map(|s| s.to_vec().into_boxed_slice()));
    set
}

/// Mod-aware dedup key: bare residues plus per-position mod mass (1e-5 Da units).
/// Avoids string formatting on the hot path.
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

fn merge_unique_candidate_idxs(into: &mut Vec<u32>, from: &[u32]) {
    // O(k) dedup: track membership in an FxHashSet seeded from `into`'s
    // current contents while preserving the SAME output order as the prior
    // `Vec::contains` (O(k^2)) implementation — existing entries first, then
    // new unique `from` entries in `from` order.
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

    fn psm(primary: u32, rank: f32, pin: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![primary],
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
            cand(b"PEPTIDEK", false), // a target peptide
            cand(b"PEPTIDEK", true),  // a decoy that coincidentally equals the target
            cand(b"EDITPEPK", true),  // a genuine decoy (no target twin)
        ];
        let targets: FxHashSet<Box<[u8]>> = [b"PEPTIDEK".to_vec().into_boxed_slice()]
            .into_iter()
            .collect();
        relabel_collision_decoys_with(&mut cands, &targets);
        assert!(!cands[0].is_decoy, "target unchanged");
        assert!(!cands[1].is_decoy, "collision decoy relabeled to target");
        assert!(cands[2].is_decoy, "genuine decoy left as decoy");
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
        // Same rounded score bucket (60) but different float rank_score — merge to best.
        let psms = vec![psm(0, 59.6, 50.0), psm(0, 60.4, 50.0)];
        let out = dedup_pepseq_score(psms, &cands);
        assert_eq!(out.len(), 1);
        assert!((out[0].rank_score - 60.4).abs() < f32::EPSILON);
    }
}

/// Best RawScore among a spectrum's Reference-entry candidates, kept so a row
/// can be compared with the best Reference candidate of a *different*
/// peptidoform: a peptide shared by a canonical entry and an isoform or chain
/// must not compete with its own canonical copy.
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

    /// The best score among peptidoforms other than `key`.
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
