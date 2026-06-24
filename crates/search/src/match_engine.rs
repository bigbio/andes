//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::{smallvec, SmallVec};

use model::aa_set::AminoAcidSet;
use input::Ms1Link;
use crate::candidate_gen::{
    enumerate_candidates, lazy_candidates_for_nominal_window, Candidate,
};
use crate::candidate_index::MmapCandidateIndex;
use model::enzyme::Enzyme;
use model::mass::{nominal_from, H2O, PROTON};
use model::peptide::Peptide;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::precursor_matching::{matches_precursor, MassError};
use crate::psm::{PsmFeatures, PsmMatch, TopNQueue};
use scoring_crate::scoring::fragment_ions::{IonKind, predict_by_ions, predict_by_ions_with_losses};
use crate::search_index::SearchIndex;
use crate::search_params::{ScoreMode, SearchParams};
use scoring_crate::intensity_model::IntensityModel;
use scoring_crate::mod_site_features::{
    mod_site_features, FEAT_SITE_DET_COUNT, FEAT_SITE_INTENS_FRAC, FEAT_SITE_LOCALIZED,
    FEAT_SITE_SHIFTED_FRAC, FEAT_SITE_SHIFTED_MATCHED,
};
use scoring_crate::scoring::{
    frag_llr_battery, fuse_strong_score, intensity_signal, mass_competition_evidence,
    psm_edge_score, rich_ion_llr, score_psm,
    strong_score_calibrated, RankScorer, OnlineStats, ScoredSpectrum, StrongScoreInputs,
    DENSITY_HW,
};
use model::spectrum::Spectrum;

// One-time-built state shared across every chunk of a streamed search.
//
// `match_spectra` materializes its full set of candidates, bucket index,
// distinct-peptide counts, and enzyme-registered aa_set in a single pass at
// startup. For chunked / streaming spectrum loading we want to reuse that
// state instead of rebuilding it per chunk. `PreparedSearch::prepare` does
// the setup once; `PreparedSearch::run_chunk` runs the per-spectrum scoring
// loop on any slice of `Spectrum`s using that prepared state.
//
// The two-pass split mirrors the original `match_spectra` body — there is
// no algorithmic change. Pre-existing single-call callers can still use
// `match_spectra(...)` which is now a thin wrapper around
// `prepare` + a single `run_chunk` call.

/// Selects how `PreparedSearch` resolves the per-spectrum candidate set.
///
/// `Ram` (default) materializes the full `Vec<Candidate>` + mass-bucket index at
/// `prepare` time and resolves each spectrum's window through `bucket_index`
/// (the original, byte-for-byte unchanged path). `Mmap` keeps only the
/// out-of-core base-peptide index resident and resolves each spectrum's
/// candidates lazily via [`crate::candidate_gen::lazy_candidates_for_nominal_window`]
/// (the RAM nominal-bucket predicate, applied per lazily-fetched candidate); the materialized
/// candidates are accumulated into `self.candidates` during the scan so the
/// PIN/TSV writers resolve `candidate_idxs` exactly as in `Ram` mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CandidateBacking {
    /// In-RAM `Vec<Candidate>` + `bucket_index`. Default; byte-identical path.
    #[default]
    Ram,
    /// Out-of-core mmap'd base-peptide index + lazy per-spectrum mod enumeration.
    ///
    /// # Identity guarantee (result-equivalent, not byte-identical at scale)
    ///
    /// `Mmap` is verified **identification-equivalent** to `Ram`: the scored
    /// candidate **set** per spectrum is identical (it reproduces `Ram`'s exact
    /// `candidate_nominal_bounds` window via `lazy_candidates_for_nominal_window`),
    /// and on real data (PXD001468 b1931 + human_entrap, cam_only) the **top-1
    /// peptide matches on 99.97% of scans** with identical `RankScore` and FDR
    /// decisions.
    ///
    /// It is **NOT byte-identical at scale.** `Ram`'s `strong_score` (the
    /// `RawScore` PIN column) depends on **order-dependent listwise features**
    /// (`CandidateRankEntropy`, `listwise_score_gap`) computed over a
    /// capacity-limited top-K queue whose `could_win` pruning depends on the
    /// candidate *visitation order* — an `enumerate_candidates` emission-order
    /// artifact, not a spectral property. `Mmap` visits candidates in a different
    /// (lazy, mass-window) order, so on a minority of scans it retains a slightly
    /// different top-K → `strong_score` can differ by ≤~0.2, occasionally flipping
    /// a shared-peptide primary or a near-tie. This is **FDR-neutral in practice**
    /// (primary identification is preserved). Cosmetic PIN-order diffs
    /// (SpecId row index, `Proteins` column order) also occur.
    ///
    /// Byte-identity holds on small DBs / fixtures where the orders coincide. Full
    /// bit-identity at scale would require either reproducing `Ram`'s exact
    /// visitation order in the lazy path OR making the listwise features /
    /// `could_win` pruning order-independent in BOTH backings (the latter removes
    /// the latent `Ram` order-dependence) — both deferred follow-ups.
    Mmap,
}

pub struct PreparedSearch<'a> {
    pub idx: &'a SearchIndex,
    pub params: &'a SearchParams,
    pub scorer: &'a RankScorer,
    pub fragment_tolerance_da: f64,
    /// Final, deduplicated candidate list (target + decoy).
    ///
    /// In `Ram` mode this is the full enumeration built at `prepare` time. In
    /// `Mmap` mode it starts empty and is grown lazily during `run_chunk` from
    /// the per-spectrum materialized candidates (so downstream resolution of
    /// `PsmMatch::candidate_idxs` is identical to the `Ram` path).
    pub candidates: Vec<Candidate>,
    /// `nominal(peptide.mass() - H2O)` → indices into `candidates`. Empty in
    /// `Mmap` mode (candidate resolution goes through `mmap_index` instead).
    pub bucket_index: BTreeMap<i32, Vec<usize>>,
    /// Candidate-resolution mode. `Ram` (default) leaves every code path
    /// byte-identical to before this field existed.
    pub backing_mode: CandidateBacking,
    /// Out-of-core base-peptide index, present only in `Mmap` mode.
    pub mmap_index: Option<MmapCandidateIndex>,
    /// `Mmap`-mode accumulator: the per-spectrum materialized candidates appended
    /// (under lock) during the parallel scan, keyed for de-duplication so a
    /// candidate scored across several chunks / spectra gets one stable global
    /// index. Drained into `self.candidates` after each chunk's scan. `None` in
    /// `Ram` mode (the field is never touched, so `Ram` stays unchanged).
    mmap_accum: Option<Mutex<MmapAccumulator>>,
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
    /// Optional context-intensity model for the strong-score signal numerator
    /// (`IntensitySignal` PIN column). `None` unless the binary supplied
    /// `--intensity-model`; when absent the column is 0.0 and ranking is unchanged.
    pub intensity_model: Option<Arc<IntensityModel>>,
    /// `Mmap`-mode ONLY: the GLOBAL set of target bare-peptide sequences (every
    /// non-decoy base record in the index), built once at `prepare_mmap`. The
    /// per-spectrum collision-decoy relabel uses THIS set so a decoy is relabeled
    /// iff its bare sequence matches ANY DB target — identical to the in-RAM path,
    /// which runs `relabel_collision_decoys` once globally over all candidates.
    /// Relabeling per-spectrum from only the window candidates is WRONG: a decoy
    /// whose target twin is absent from the current precursor window would keep its
    /// `is_decoy` flag and diverge the TDC label from RAM. Bare sequences only, so
    /// this is bounded (far smaller than the full candidate `Vec` we avoid in
    /// out-of-core mode). `None` in `Ram` mode.
    mmap_target_bare_seqs: Option<std::collections::HashSet<Box<[u8]>>>,
}

/// Global candidate identity for the `Mmap`-mode accumulator: enough to dedup a
/// peptidoform reached from multiple spectra/charges/isotope-offsets to one
/// stable global index (mirrors the per-candidate identity the in-RAM `candidates`
/// Vec carries). The peptidoform mod state is folded into `mod_units` (per-residue
/// mod-mass in 1e-5 Da units) so two placements of the same residues at the same
/// span but different mods get distinct indices, exactly as the in-RAM enumeration
/// produces distinct candidates.
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

/// `Mmap`-mode global candidate store: maps a candidate's identity to its stable
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

/// Owned, precursor-tolerance-independent products of [`PreparedSearch::prepare`]
/// (candidates, bucket index, GF aa_set). Built once during the
/// calibration pre-pass and reused for the tightened-tolerance main pass, avoiding
/// a second full candidate enumeration (~15s on the 16.8M-candidate Astral search).
pub struct PreparedParts {
    candidates: Vec<Candidate>,
    bucket_index: BTreeMap<i32, Vec<usize>>,
    aa_set_for_scoring: AminoAcidSet,
}

/// Derive the inclusive `[min_nominal, max_nominal]` nominal-mass bucket bounds
/// for one spectrum at one charge state `z`, used to enumerate candidate
/// peptides from the mass-bucket index.
///
/// The window is centered on the precursor m/z and widened by the precursor
/// tolerance plus the isotope-error range. This is the original, byte-for-byte
/// unchanged derivation used by both the standard search and the cascade's
/// NARROW Pass 1.
pub(crate) fn candidate_nominal_bounds(
    spec: &Spectrum,
    z: u8,
    params: &SearchParams,
    shift_ppm: f64,
) -> (i32, i32) {
    let charge_f = z as f64;
    let iso_min = *params.isotope_error_range.start() as i32;
    let iso_max = *params.isotope_error_range.end() as i32;

    let neutral_mass = adjusted_observed_neutral_mass(
        (spec.precursor_mz - PROTON) * charge_f - H2O,
        shift_ppm,
    );
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
    (min_nominal, max_nominal)
}

impl<'a> PreparedSearch<'a> {
    /// Build the per-search state once. Enumerates candidates, builds the
    /// mass-bucket index, and clones+registers the aa_set for cleavage
    /// scoring.
    pub fn prepare(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        decoy_prefix: &str,
    ) -> Self {
        // Collect the production candidate list.
        let mut candidates: Vec<Candidate> =
            enumerate_candidates(idx, params, decoy_prefix).collect();

        // Decoy↔target peptide-collision removal: a decoy peptide whose bare
        // sequence coincides with a real target peptide (palindromes, reversal-
        // coincident, short peptides) is not a valid decoy — it IS a true target
        // sequence. Relabel it as target so per-spectrum dedup merges it with its
        // target twin instead of inflating the decoy count (which biases TDC
        // conservative and injects label noise into decoy-aware training).
        relabel_collision_decoys(&mut candidates);

        // Build mass-bucket index: nominal(peptide.mass() - H2O) → Vec<candidate_idx>.
        //
        // Uses the same nominal_from convention as the node-scoring DP's
        // nominal-mass binning so that bucket keys align with that lookup.
        // Stores only indices into `candidates` — no cloning, tiny memory overhead.
        let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
        for (cand_idx, cand) in candidates.iter().enumerate() {
            let nominal = cand.peptide.nominal_residue_mass();
            bucket_index.entry(nominal).or_default().push(cand_idx);
        }

        // Build an aa_set clone with enzyme registered (for cleavage scoring).
        // Defaults: peptide_eff = 0.95, neighboring_eff = 0.95.
        // Cloning is cheap (AminoAcidSet is a HashMap of ~20 entries).
        // This avoids mutating the shared SearchParams.aa_set borrow.
        let mut aa_set_for_scoring: AminoAcidSet = params.aa_set.clone();
        if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
            aa_set_for_scoring.register_enzyme(params.enzyme, 0.95, 0.95);
        }

        PreparedSearch {
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            candidates,
            bucket_index,
            backing_mode: CandidateBacking::Ram,
            mmap_index: None,
            mmap_accum: None,
            aa_set_for_scoring,
            ms1_link: None,
            intensity_model: None,
            mmap_target_bare_seqs: None,
        }
    }

    /// Build the per-search state for the OUT-OF-CORE (`Mmap`) candidate path.
    ///
    /// Builds the mass-sorted base-peptide index file (via
    /// [`build_base_peptide_index`](crate::candidate_index::build_base_peptide_index))
    /// at `index_path`, opens it as an [`MmapCandidateIndex`], and arms the lazy
    /// per-spectrum resolution path. Unlike [`Self::prepare`] this does NOT
    /// materialize the full `Vec<Candidate>`; `self.candidates` starts empty and
    /// is grown lazily during the scan so the PIN/TSV writers still resolve
    /// `candidate_idxs` against a real candidate slice.
    ///
    /// `decoy_prefix` MUST match the prefix the `SearchIndex` was built with so
    /// the index's per-record decoy flags align with the in-RAM enumeration.
    ///
    /// # Byte-identity scope
    ///
    /// The resulting PIN is byte-identical to [`Self::prepare`] (Ram mode) **only**
    /// for fully-tryptic searches (`num_tolerable_termini = 2`) with at most one
    /// variable mod per residue. For semi-tryptic or multi-mod-per-residue searches
    /// the candidate set and scores are identical but the per-scan candidate order
    /// (and therefore PIN row order / primary-peptide / Proteins-column order) may
    /// differ — see [`CandidateBacking::Mmap`] for details.
    pub fn prepare_mmap(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        decoy_prefix: &str,
        index_path: &std::path::Path,
    ) -> std::io::Result<Self> {
        // Out-of-core mmap mode cannot serve the chimeric Pass-2 (it needs the
        // in-RAM candidate index); reject the combination before building anything.
        reject_unsupported_mmap_combo(params)?;
        let (mmap_index, was_built) =
            MmapCandidateIndex::open_or_build(index_path, idx, params, decoy_prefix)?;
        if was_built {
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

        let mut aa_set_for_scoring: AminoAcidSet = params.aa_set.clone();
        if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
            aa_set_for_scoring.register_enzyme(params.enzyme, 0.95, 0.95);
        }

        // RC1 fix: build the GLOBAL target bare-sequence set ONCE from the full
        // base-peptide index (mirrors the in-RAM `relabel_collision_decoys`, which
        // runs once globally over all candidates). Per-spectrum relabel from only
        // the window candidates is WRONG (a collision decoy whose target twin is
        // outside the current window would keep its decoy flag → TDC label
        // diverges from RAM). Bare residues are reconstructed from each non-decoy
        // record's protein span — bounded, far smaller than the candidate Vec.
        let target_bare_seqs = build_target_bare_seqs(&mmap_index, idx);

        Ok(PreparedSearch {
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            candidates: Vec::new(),
            bucket_index: BTreeMap::new(),
            backing_mode: CandidateBacking::Mmap,
            mmap_index: Some(mmap_index),
            mmap_accum: Some(Mutex::new(MmapAccumulator::default())),
            aa_set_for_scoring,
            ms1_link: None,
            intensity_model: None,
            mmap_target_bare_seqs: Some(target_bare_seqs),
        })
    }

    /// Drain the `Mmap`-mode accumulator into `self.candidates` so downstream
    /// consumers (PIN/TSV writers, chimeric Pass 2) resolve `candidate_idxs`
    /// against the final materialized candidate slice. No-op in `Ram` mode (the
    /// candidates already live in `self.candidates`). Idempotent.
    pub fn sync_materialized_candidates(&mut self) {
        if let Some(accum) = self.mmap_accum.as_ref() {
            let materialized = std::mem::take(&mut accum.lock().unwrap().candidates);
            self.candidates = materialized;
        }
    }

    /// Consume this prepared search into owned [`PreparedParts`], dropping the
    /// `idx`/`params`/`scorer` borrows so the caller can mutate `params` (e.g.
    /// tighten the precursor tolerance) before rebuilding via [`Self::from_parts`].
    pub fn into_parts(self) -> PreparedParts {
        PreparedParts {
            candidates: self.candidates,
            bucket_index: self.bucket_index,
            aa_set_for_scoring: self.aa_set_for_scoring,
        }
    }

    /// Rebuild a [`PreparedSearch`] from previously-enumerated [`PreparedParts`]
    /// and a (possibly mutated) `params`, reusing the candidate enumeration. The
    /// parts are precursor-tolerance independent, so this is correct after
    /// calibration tightens the tolerance. `ms1_link` starts `None`.
    pub fn from_parts(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        parts: PreparedParts,
    ) -> Self {
        PreparedSearch {
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            candidates: parts.candidates,
            bucket_index: parts.bucket_index,
            backing_mode: CandidateBacking::Ram,
            mmap_index: None,
            mmap_accum: None,
            aa_set_for_scoring: parts.aa_set_for_scoring,
            ms1_link: None,
            intensity_model: None,
            mmap_target_bare_seqs: None,
        }
    }

    /// Attach an optional [`IntensityModel`] for the S1 signal PIN column.
    pub fn with_intensity_model(mut self, model: Option<Arc<IntensityModel>>) -> Self {
        self.intensity_model = model;
        self
    }

    /// Attach an [`Ms1Link`] for the chimeric precursor isotope features.
    /// Builder-style so existing `prepare` callers are unchanged;
    /// only the binary's `--chimeric` path calls this. Has NO effect unless
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
    pub fn run_chunk(
        &self,
        spectra: &[Spectrum],
        spectrum_idx_offset: usize,
    ) -> Vec<TopNQueue> {
        self.run_chunk_with_params(spectra, spectrum_idx_offset, self.params)
    }

    /// Like [`Self::run_chunk`] but uses `params` for precursor matching and
    /// queue sizing instead of the params reference stored at prepare time.
    /// Candidate enumeration and mass buckets are unchanged.
    pub fn run_chunk_with_params(
        &self,
        spectra: &[Spectrum],
        spectrum_idx_offset: usize,
        params: &SearchParams,
    ) -> Vec<TopNQueue> {
        self.run_chunk_inner(spectra, spectrum_idx_offset, params)
    }

    fn run_chunk_inner(
        &self,
        spectra: &[Spectrum],
        spectrum_idx_offset: usize,
        params: &SearchParams,
    ) -> Vec<TopNQueue> {
        let scorer = self.scorer;
        let fragment_tolerance_da = self.fragment_tolerance_da;
        let candidates = &self.candidates;
        let bucket_index = &self.bucket_index;
        let aa_set_for_scoring = &self.aa_set_for_scoring;
        // Out-of-core (`Mmap`) backing handles, `None` on the default `Ram` path.
        let mmap_index = self.mmap_index.as_ref();
        let mmap_accum = self.mmap_accum.as_ref();

        // Yield-accounting counters.
        // Aggregated across all worker threads via Relaxed atomics — exact counts
        // don't require ordering with other memory ops.
        let skipped_min_peaks = AtomicU64::new(0);
        let candidates_visited = AtomicU64::new(0);
        let psms_pushed = AtomicU64::new(0);
        let spectra_with_psms = AtomicU64::new(0);

        // Research diagnostic (env-gated, chimeric only): measure the
        // shared-fragment overlap between the top-2 co-emitted distinct peptides
        // per scan. Tests the "fragment theft" hypothesis behind chimeric FDR
        // inflation. Zero cost unless ANDES_CHIMERIC_OVERLAP is set AND --chimeric.
        let chim_overlap = params.chimeric
            && std::env::var("ANDES_CHIMERIC_OVERLAP").is_ok();

        // Parallel per-spectrum search. All inputs above are `&` immutable; the
        // closure owns its TopNQueue, scored_per_charge cache, and per-bin GF state.
        let queues: Vec<TopNQueue> = spectra
            .par_iter()
            .enumerate()
            .map(
                |(local_idx, spec)| {
                let spec_idx = local_idx + spectrum_idx_offset;
                let retention_cap = params.hot_loop_retention_cap();
                let mut queue = TopNQueue::new(retention_cap);

            // Cascade Pass 1 is a NARROW precursor-window search, identical to the
            // non-chimeric path; Pass 2 (`run_pass2_coisolation`) is the only
            // chimeric-specific enumeration. MS1 load + the
            // `precursor_isotope_match` feature fill still run for Pass 2.

            // Skip spectra with too few peaks.
            if spec.peaks.len() < params.min_peaks as usize {
                skipped_min_peaks.fetch_add(1, Ordering::Relaxed);
                return queue;
            }

            // DeltaRawScore capture (additive feature): track the best and the
            // second-best *distinct-peptide* rounded RawScore across ALL
            // mass-matching candidates for this spectrum — independent of the
            // TopNQueue, so it survives `top_n = 1` (where the runner-up is
            // evicted) AND never feeds the GF `min_score` (so no emitted PSM's
            // SpecEValue is perturbed). Distinct-peptide keyed by nominal
            // residue mass: collapses the same peptide scored at multiple
            // charges / matched in multiple proteins (which share a neutral
            // mass) so a shared peptide doesn't zero the delta. Isobaric
            // different-sequence peptides at the same nominal mass are rare and
            // only make the delta conservative (slightly understated).
            let mut best_raw = i32::MIN;
            let mut best_mass_key = i32::MIN;
            let mut second_raw = i32::MIN;

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
            let charges_to_try: SmallVec<[u8; 4]> = match spec.precursor_charge {
                Some(z) if z > 0 => smallvec![z as u8],
                _ => params.charge_range.clone().collect(),
            };

            // Build (and cache) a ScoredSpectrum per charge to evaluate.
            //
            // A single ScoredSpectrum keyed off `spec.precursor_charge.unwrap_or(2)`
            // would force charge-missing spectra to use z=2 even when evaluating
            // z=3 candidates — wrong precursor filtering, wrong partition, wrong
            // main_ion.
            //
            // For charge-explicit spectra the cache has exactly 1 entry (no overhead).
            // For charge-missing spectra, typically 2-3 entries per spectrum.
            let mut scored_per_charge: SmallVec<[(u8, ScoredSpectrum<'_>); 4]> = SmallVec::new();
            for &z in &charges_to_try {
                if scored_per_charge.iter().all(|(charge, _)| *charge != z) {
                    scored_per_charge.push((z, ScoredSpectrum::new(spec, scorer, z)));
                }
            }
            let scored_spec_for_charge = |z: u8| {
                scored_per_charge
                    .iter()
                    .find(|(charge, _)| *charge == z)
                    .map(|(_, spec)| spec)
                    .expect("scored spectrum exists for candidate charge")
            };

            // Compute per-charge candidate windows and union them into a deduplicated
            // set of candidate indices. Window derivation mirrors the precursor
            // mass-bin logic so any candidate admitted by matches_precursor is
            // guaranteed to be in at least one charge's window.
            //
            // Vec + sort_unstable + dedup is faster than BTreeSet for the typical
            // 1k-3k indices per spectrum: better cache locality, no tree pointer
            // chasing, single sort pass at end. Iteration order matches BTreeSet
            // (ascending), preserving downstream parity / determinism.
            let shift_ppm = params.precursor_mass_shift_ppm;

            // `Ram` path (default, byte-identical): coarse mass-bucket lookup.
            // Produces a sorted+deduped list of indices into `self.candidates`.
            //
            // `Mmap` path: lazily materialize the per-spectrum candidate set from
            // the out-of-core base-peptide index (Task 3 `lazy_candidates_for_precursor`)
            // for each (charge, isotope-offset) the in-RAM path would consider,
            // then sort into the SAME enumeration order the in-RAM `candidates`
            // Vec uses — `(protein_index, start_offset, mod placement)` — so the
            // candidate iteration order (and therefore queue tie-breaking and
            // `candidate_idxs` aggregation order) is identical. The actual
            // precursor gate is still `matches_precursor` inside the loop, exactly
            // as on the `Ram` path; the lazy retrieval only decides which
            // candidates are *looked at*. Collision-decoy relabeling is applied to
            // the per-spectrum set — a collision decoy and its target twin share a
            // bare sequence (hence mass), so both land in the same window.
            let mut window_cand_indices: Vec<usize> = Vec::new();
            let mut mmap_cands: Vec<Candidate> = Vec::new();
            match self.backing_mode {
                CandidateBacking::Ram => {
                    window_cand_indices.reserve(2048);
                    for &z in &charges_to_try {
                        let (min_nominal, max_nominal) =
                            candidate_nominal_bounds(spec, z, params, shift_ppm);
                        for (_nm, idxs) in bucket_index.range(min_nominal..=max_nominal) {
                            window_cand_indices.extend_from_slice(idxs);
                        }
                    }
                    window_cand_indices.sort_unstable();
                    window_cand_indices.dedup();
                }
                CandidateBacking::Mmap => {
                    let mi = mmap_index.expect("Mmap backing requires an mmap index");
                    // Union of candidates over the same (charge, isotope-offset)
                    // grid the in-RAM precursor match scans, deduped by peptidoform
                    // identity (residues+mods+protein+offset+termini+decoy).
                    //
                    // `base_mult` reproduces the in-RAM enumeration MULTIPLICITY: a
                    // base peptide can appear in `enumerate_candidates` more than
                    // once at the SAME coordinates (notably the N-terminal-Met
                    // re-enumeration emits an identical-coordinate span), and the
                    // in-RAM `candidates` Vec keeps each copy as a distinct index.
                    // `dedup_pepseq_score` later aggregates those indices into one
                    // PSM's `candidate_idxs` — so the PIN `Proteins` column repeats
                    // the accession once per copy. `lazy_candidates_for_precursor`
                    // de-duplicates base records, so to stay byte-identical we count
                    // each base record's raw multiplicity here and replicate the
                    // expanded peptidoforms that many times. The on-disk index
                    // preserves the duplicates (it is built from the same
                    // `enumerate_candidates`), so the count is exact.
                    let mut seen: FxHashSet<GlobalCandKey> = FxHashSet::default();
                    let mut deduped: Vec<Candidate> = Vec::new();
                    for &z in &charges_to_try {
                        // RESULT-IDENTITY: the in-RAM (`Ram`) backing scores
                        // exactly the candidates whose nominal residue mass lands
                        // in the coarse integer bucket window
                        // `candidate_nominal_bounds(spec, z, ..)` (via
                        // `bucket_index.range(min..=max)`). The earlier f64
                        // directional milli-Da window used here was a DIFFERENT set
                        // at bucket boundaries (proven divergence on scan 33203:
                        // RAM kept an iso=+2 candidate, mmap an iso=−1). To stay
                        // result-identical to RAM we derive the SAME nominal bounds
                        // (reusing `candidate_nominal_bounds`, the single source of
                        // truth) and gate per-candidate ACCEPTANCE on the EXACT RAM
                        // nominal-bucket predicate
                        // `nominal_from(peptide.mass() − H2O) ∈ [min, max]`.
                        // `lazy_candidates_for_nominal_window` over-fetches base
                        // records over a padded mass window (bounded-RSS preserved)
                        // and applies that predicate — making the scored multiset
                        // identical to RAM's `bucket_index.range(..)` set.
                        let (min_nominal, max_nominal) =
                            candidate_nominal_bounds(spec, z, params, shift_ppm);
                        for cand in lazy_candidates_for_nominal_window(
                            mi, self.idx, params, min_nominal, max_nominal,
                        ) {
                            if seen.insert(GlobalCandKey::from_candidate(&cand)) {
                                deduped.push(cand);
                            }
                        }
                    }
                    // Replicate each deduped peptidoform by its base record's raw
                    // multiplicity, so the per-spectrum candidate set matches the
                    // in-RAM enumeration copy-for-copy (collapsed identically by the
                    // later `dedup_pepseq_score`). The raw count is a fixed property
                    // of the index: query it at the candidate's EXACT base mass
                    // (peptide mass minus its mod deltas) and count records sharing
                    // the same (protein, start, length, flags) key.
                    for cand in deduped {
                        let mult = base_record_multiplicity(mi, &cand);
                        for _ in 0..mult {
                            mmap_cands.push(cand.clone());
                        }
                    }
                    // Collision-decoy relabeling against the GLOBAL target
                    // bare-sequence set (built once at `prepare_mmap` from the full
                    // index). This is IDENTICAL to the in-RAM path, which relabels
                    // once globally over all candidates: a decoy is relabeled iff
                    // its bare sequence matches ANY DB target. Relabeling against
                    // only this window's targets would be WRONG — a collision decoy
                    // whose target twin is outside the current precursor window
                    // would keep its decoy flag and diverge the TDC label from RAM.
                    let target_bare_seqs = self
                        .mmap_target_bare_seqs
                        .as_ref()
                        .expect("Mmap backing requires the global target bare-seq set");
                    relabel_collision_decoys_with(&mut mmap_cands, target_bare_seqs);
                    // Canonical enumeration order: protein_index ASC, then start
                    // offset ASC.  Within each (protein, offset) group, candidates
                    // retain the order produced by `lazy_candidates_for_nominal_window`
                    // → `expand_mod_combinations`, which is the SAME order
                    // `enumerate_candidates` uses.  `sort_by` is stable, so the
                    // within-group ordering from `lazy_candidates_for_nominal_window`
                    // is preserved — reproducing the RAM path's ascending global-index
                    // iteration order exactly.
                    //
                    // DO NOT add residues/mod_units tie-breaking here: lexicographic
                    // residue order diverges from `expand_mod_combinations` recursive
                    // output order and causes tie-breaking differences vs RAM.
                    mmap_cands.sort_by(|a, b| {
                        a.protein_index
                            .cmp(&b.protein_index)
                            .then(a.start_offset_in_protein.cmp(&b.start_offset_in_protein))
                    });
                    window_cand_indices = (0..mmap_cands.len()).collect();
                }
            }
            // Per-spectrum candidate resolution slice, uniform across both
            // backings. In `Ram` mode this borrows `self.candidates` and the PSM
            // `candidate_idxs` are already GLOBAL indices. In `Mmap` mode it
            // borrows the per-spectrum `mmap_cands` and the PSM `candidate_idxs`
            // start as LOCAL slot indices into `mmap_cands`; dedup + feature
            // extraction resolve against this same slice exactly as `Ram` resolves
            // against `self.candidates`. The local→global remap happens once at the
            // END of the closure (see `mmap_global_idx` below), so everything in
            // between is identical to the in-RAM path.
            let cand_slice: &[Candidate] = match self.backing_mode {
                CandidateBacking::Ram => candidates,
                CandidateBacking::Mmap => &mmap_cands,
            };
            // Resolve a per-spectrum candidate slot to its `&Candidate` and the
            // value to store in `candidate_idxs`. For both backings this is the
            // slot's index into `cand_slice` (global for `Ram`, local for `Mmap`).
            let resolve_cand = |slot: usize| -> (&Candidate, u32) {
                (&cand_slice[slot], slot as u32)
            };

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

            // Per-candidate cleavage credit:
            //   `cleavage_score = n_term_cleavage_score + c_term_cleavage_score`
            // added to the raw PSM score before queue insertion.
            //
            // Use the ENZYME-REGISTERED aa_set (cleavage credit/penalty are
            // populated by register_enzyme — params.aa_set is unregistered).
            //
            // `fn` (not closure) + `#[inline(always)]` ensures LLVM
            // monomorphizes + inlines into the candidate loop. Closure form
            // was not being inlined and went through FnMut::call_mut dispatch.
            #[inline(always)]
            #[allow(clippy::too_many_arguments, reason = "private inner driver for the per-chunk search loop; all args are orthogonal cleavage parameters")]
            fn compute_cleavage_credit(
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

            // Per-charge queue keyed by charge state. Retains top-N PSMs
            // independently per precursor charge (Kim et al., Nat Commun 5:5277, 2014).
            let mut per_charge_queues: FxHashMap<u8, TopNQueue> = FxHashMap::default();

            // Cascade Pass 1 is a NARROW brute-force window scan: iterate every
            // candidate whose nominal mass falls in the precursor window.
            for &cand_slot in &window_cand_indices {
                let (cand, cand_idx) = resolve_cand(cand_slot);
                let cleavage_credit = compute_cleavage_credit(
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
                            spec, &cand.peptide, z, offset,
                            &params.precursor_tolerance, shift_ppm,
                        );
                        if let Some(err) = matched {
                            iso_errs.push(err);
                        }
                    }
                    if iso_errs.is_empty() {
                        continue;
                    }

                    // Compute pin_score ONCE (iso-independent).
                    let pin_score = score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da)
                        + cleavage_credit;

                    // DeltaRawScore capture: fold this candidate's RawScore into
                    // the per-spectrum best / second-best-distinct-peptide
                    // trackers BEFORE the `could_win` gate, so a strong runner-up
                    // that can't beat the top-1 still contributes to the delta.
                    {
                        let s = pin_score.round() as i32;
                        // Tailor histogram: fold this candidate's rounded RawScore
                        // in (same `s` the DeltaRawScore tracker uses). One u32
                        // increment, no per-candidate allocation.
                        *tailor_hist.entry(s).or_insert(0) += 1;
                        tailor_total += 1;
                        strong_null_stats.push(pin_score);
                        let mkey = cand.peptide.nominal_residue_mass();
                        if mkey == best_mass_key {
                            if s > best_raw {
                                best_raw = s;
                            }
                        } else if s > best_raw {
                            second_raw = second_raw.max(best_raw);
                            best_raw = s;
                            best_mass_key = mkey;
                        } else {
                            second_raw = second_raw.max(s);
                        }
                    }

                    // Gate against the queue's current worst rank_score
                    // before invoking edge_score.
                    let could_win = match per_charge_queues.get(&z) {
                        Some(q) if q.len() >= q.capacity() as usize => {
                            q.worst_rank_score()
                                .is_none_or(|worst| pin_score + max_edge_bonus > worst)
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
                    let err = iso_errs.into_iter()
                        .min_by(|a, b| a.mass_error_ppm.abs().partial_cmp(&b.mass_error_ppm.abs()).unwrap_or(std::cmp::Ordering::Equal))
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
                    };
                    per_charge_queues
                        .entry(z)
                        .or_insert_with(|| TopNQueue::new(retention_cap))
                        .push(psm);
                    psms_pushed.fetch_add(1, Ordering::Relaxed);
                }
            }
            candidates_visited.fetch_add(window_cand_indices.len() as u64, Ordering::Relaxed);

            // pepSeq + score dedup per-charge BEFORE GF compute.
            // Same peptide matched against multiple proteins collapses to one
            // PsmMatch with aggregated candidate_idxs (peptide-sequence dedup).
            for queue in per_charge_queues.values_mut() {
                if queue.len() > 1 {
                    let drained = queue.drain_into_vec();
                    let deduped = dedup_pepseq_score(drained, cand_slice);
                    for psm in deduped {
                        queue.push(psm);
                    }
                }
            }

            // Candidate selection/ranking already happened on `rank_score`
            // (RawScore = node + cleavage + edge). The generating function has
            // been removed, so there is no per-charge SpecEValue DP — the
            // per-charge queues are merged directly into the spectrum queue.
            let any_queue_nonempty = per_charge_queues.values().any(|q| !q.is_empty());
            if any_queue_nonempty {
                spectra_with_psms.fetch_add(1, Ordering::Relaxed);
            }

            // Spectrum-level merge. The TopNQueue::push (Ordering::Equal arm)
            // keeps rank_score ties at capacity, matching Kim et al. (Nat Commun
            // 5:5277, 2014) tie-retention on the spectrum-level merge.
            for (_charge, mut per_charge) in per_charge_queues.drain() {
                for psm in per_charge.drain_into_vec() {
                    queue.push(psm);
                }
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
            let tailor_denom = crate::psm::tailor_denominator(&tailor_hist, tailor_total) as f32;

            // S2 listwise terms over the retained top-K queue (not the full scored pool).
            let mut retained_scores: Vec<f32> =
                queue.iter_psms().map(|p| p.score).collect();
            retained_scores.sort_by(|a, b| {
                b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
            });
            let listwise_score_gap =
                scoring_crate::scoring::listwise_score_gap(&retained_scores);
            let candidate_rank_entropy =
                scoring_crate::scoring::candidate_rank_entropy(&retained_scores);

            queue.fill_post_topn(|psm| {
                let ss = scored_spec_for_charge(psm.charge_used);
                let cand = &cand_slice[psm.primary_candidate_idx() as usize];
                let mut features = compute_psm_features(
                    ss,
                    &cand.peptide,
                    scorer,
                    psm.charge_used,
                    self.intensity_model.as_deref(),
                );
                features.edge_score = psm.edge_score; // reuse per-candidate value

                // The PIN `precursor_isotope_kl` / `precursor_snr` columns are
                // left at their 0.0 defaults here. They were a per-PSM MS1
                // isotope-envelope feature of the old single-pass chimeric mode;
                // the two-pass cascade does not consume them (computing them ran
                // ~121k times on Astral for ~2:40 of wall). The cascade uses MS1
                // only in Pass 2's `run_pass2_coisolation` / `detect_coisolated`
                // co-isolation detection — `ms1_link` and its loading stay.

                features.delta_raw_score = delta_raw_score;
                features.listwise_score_gap = listwise_score_gap;
                features.candidate_rank_entropy = candidate_rank_entropy;
                // Tailor calibrated score: this PSM's RawScore divided by the
                // spectrum's top-1% quantile RawScore. Additive — does not
                // affect ranking. Uses `psm.score` (the PIN RawScore =
                // node + cleavage), the same value emitted in the RawScore column.
                features.tailor_score = psm.score / tailor_denom;
                features.strong_score = fuse_strong_score(&StrongScoreInputs {
                    intensity_signal: features.intensity_signal,
                    chance_match_surprise: features.chance_match_surprise,
                    mass_competition_evidence: features.mass_competition_evidence,
                    candidate_rank_entropy: features.candidate_rank_entropy,
                    listwise_score_gap: features.listwise_score_gap,
                });
                psm.features = features;
            });

            let retained_strong: Vec<f32> = queue
                .iter_psms()
                .map(|p| p.features.strong_score)
                .collect();
            queue.fill_post_topn(|psm| {
                psm.features.strong_score_cal = strong_score_calibrated(
                    psm.features.strong_score,
                    &retained_strong,
                    &strong_null_stats,
                );
            });

            if params.score_mode == ScoreMode::Strong {
                queue.reorder_by_strong_score();
                // Emit only user top-N rows (gate uses top_n=1); retention pool was wider.
                queue.trim_to_capacity(params.top_n_psms_per_spectrum);
            }

            // Chimeric fragment-overlap diagnostic (env-gated). For scans that
            // emit ≥2 distinct peptides, measure how many MS2 peaks the runner-up
            // claims that the top peptide also claims (the "fragment theft" the
            // chimeric FDR inflation is hypothesized to come from).
            if chim_overlap {
                let sorted = queue.clone().into_sorted_vec(); // best-first
                let mut picks: Vec<&PsmMatch> = Vec::new();
                'outer: for psm in &sorted {
                    let seq: Vec<u8> = cand_slice[psm.primary_candidate_idx() as usize]
                        .peptide.residues.iter().map(|a| a.residue).collect();
                    for p in &picks {
                        let pseq: Vec<u8> = cand_slice[p.primary_candidate_idx() as usize]
                            .peptide.residues.iter().map(|a| a.residue).collect();
                        if pseq == seq { continue 'outer; }
                    }
                    picks.push(psm);
                    if picks.len() == 2 { break; }
                }
                if picks.len() == 2 {
                    let pa = &cand_slice[picks[0].primary_candidate_idx() as usize].peptide;
                    let pb = &cand_slice[picks[1].primary_candidate_idx() as usize].peptide;
                    let ka = matched_peak_keys(scored_spec_for_charge(picks[0].charge_used), pa, scorer);
                    let kb = matched_peak_keys(scored_spec_for_charge(picks[1].charge_used), pb, scorer);
                    let shared = ka.intersection(&kb).count();
                    let uni = ka.union(&kb).count();
                    let minlen = ka.len().min(kb.len());
                    eprintln!(
                        "CHIM_OVERLAP spec_idx={} nA={} nB={} shared={} jacc={:.3} fracmin={:.3}",
                        spec_idx, ka.len(), kb.len(), shared,
                        if uni > 0 { shared as f64 / uni as f64 } else { 0.0 },
                        if minlen > 0 { shared as f64 / minlen as f64 } else { 0.0 },
                    );
                }
            }

            // `Mmap` mode: the surviving PSMs carry LOCAL `candidate_idxs` into
            // this spectrum's `mmap_cands`. Intern each referenced candidate into
            // the shared accumulator (de-duped by peptidoform identity, so the same
            // candidate seen from another spectrum reuses one global index) and
            // rewrite `candidate_idxs` to the global indices the PIN/TSV writers
            // resolve against `self.candidates` after the scan. `Ram` mode is
            // untouched (indices are already global).
            if self.backing_mode == CandidateBacking::Mmap {
                let accum = mmap_accum.expect("Mmap backing requires an accumulator");
                let mut accum = accum.lock().unwrap();
                queue.remap_candidate_idxs(|local| accum.intern(mmap_cands[local as usize].clone()));
            }

                queue
            },
            )
            .collect();

        // Yield-accounting summary.
        // Helps disambiguate whether a PSM-yield gap comes from:
        //   - filtering (skipped_min_peaks)
        //   - enumeration (candidates_visited)
        //   - scoring (psms_pushed)
        //   - top-N retention (spectra_with_psms)
        eprintln!(
            "Yield (chunk): {} spectra in, {} skipped by min_peaks, {} candidates visited, \
             {} PSMs pushed, {} spectra with non-empty queue",
            spectra.len(),
            skipped_min_peaks.load(Ordering::Relaxed),
            candidates_visited.load(Ordering::Relaxed),
            psms_pushed.load(Ordering::Relaxed),
            spectra_with_psms.load(Ordering::Relaxed),
        );

        queues
    }
}

/// Match every spectrum against every candidate from the SearchIndex.
/// Returns one top-N PSM queue per spectrum (in input order) PLUS the
/// enumerated `Vec<Candidate>` that backs the `PsmMatch::candidate_idxs`
/// handles inside each queue.
///
/// Callers that need to resolve a PSM's peptide / protein info must hold
/// on to the returned candidates vector and look up by
/// `psm.primary_candidate_idx() as usize`. The previous API embedded a cloned
/// `Candidate` directly in every PsmMatch; that allocation cost is now
/// gone but the resolution responsibility shifts to the caller.
///
/// A `ScoredSpectrum` is built once per spectrum and reused across all
/// candidates; candidates are bucketed by mass for sub-linear precursor
/// lookup. After per-candidate scoring, SpecEValue is computed via the
/// generating-function DP across the precursor tolerance window in nominal
/// mass space and assigned to every PSM in the queue.
///
/// This is a thin wrapper around [`PreparedSearch::prepare`] +
/// [`PreparedSearch::run_chunk`] preserved for single-shot callers (tests
/// and the historic single-pass binary path).
pub fn match_spectra(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    scorer: &RankScorer,
    fragment_tolerance_da: f64,
    decoy_prefix: &str,
) -> (Vec<TopNQueue>, Vec<Candidate>) {
    let prepared = PreparedSearch::prepare(
        idx,
        params,
        scorer,
        fragment_tolerance_da,
        decoy_prefix,
    );
    let queues = prepared.run_chunk(spectra, 0);
    (queues, prepared.candidates)
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

    queues.par_iter_mut().enumerate().for_each(|(spec_idx, q)| {
        if q.is_empty() {
            return;
        }
        let Some(spec) = spectra.get(spec_idx) else {
            return;
        };

        // Linked MS1 scan for this MS2 (most-recent preceding MS1).
        let Some(Some(ms1_idx)) = link.ms2_to_ms1.get(spec_idx) else {
            return;
        };
        let Some(ms1) = link.ms1_peaks.get(*ms1_idx) else {
            return;
        };

        // Isolation window: prefer the per-scan offsets if the parser recorded
        // them, else fall back to the configured chimeric half-width.
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

        let cos = crate::coisolation::detect_coisolated(
            ms1,
            lo,
            hi,
            spec.precursor_mz,
            *params.charge_range.start()..=*params.charge_range.end(),
            tol,
            // max_kl + max_n: the chimeric-N lever (params; default 0.3 / 2 = the
            // proven +101%-Astral setting). Raising max_coisolated searches deeper
            // co-fragments on the residual; the KL gate keeps secondaries clean.
            params.chimeric_max_kl,
            params.chimeric_max_coisolated,
        );
        if cos.is_empty() {
            return;
        }

        // Primary peptide = the queue's best PSM (smallest SpecEValue).
        let primary = match q.peek_top() {
            Some(best) => {
                prepared.candidates[best.primary_candidate_idx() as usize]
                    .peptide
                    .clone()
            }
            None => return,
        };

        // Secondaries on the SAME scan compete for residual evidence: each accepted
        // secondary's matched peaks are added to `claimed` so the next co-isolated
        // mass is scored against still-unexplained signal (no double-counting of
        // shared leftover peaks across multiple co-isolated precursors).
        let mut claimed: std::collections::HashSet<i64> = std::collections::HashSet::new();
        for co in cos {
            if let Some((mut psm, winner_claimed)) = crate::coisolation::search_secondary(
                spec,
                &primary,
                &claimed,
                co,
                &prepared.candidates,
                &prepared.bucket_index,
                prepared.scorer,
                &prepared.aa_set_for_scoring,
                params,
                prepared.fragment_tolerance_da,
                prepared.intensity_model.as_deref(),
            ) {
                // `spec_idx` is chunk-local (indexes `spectra` + `link`); the
                // emitted spectrum_idx must be global to align with `all_spectra`.
                psm.spectrum_idx = offset + spec_idx;
                // Distinct co-isolated peptide — a legitimate EXTRA emission, not a
                // competitor for the primary's slot. force_push skips capacity-based
                // eviction (plain `push` at capacity 1 would evict the primary or
                // drop the secondary).
                q.force_push(psm);
                claimed.extend(winner_claimed);
            }
        }
    });
}

/// Per-candidate cleavage credit, module-level mirror of the nested
/// `compute_cleavage_credit` in `run_chunk_inner`. `search_secondary` calls this
/// so secondaries share the production RawScore scale (`score_psm(...) + credit`).
/// Credit/penalty constants come from the ENZYME-REGISTERED `aa_set`; keep the
/// branch structure bit-identical to `compute_cleavage_credit`.
pub(crate) fn cleavage_credit_for(cand: &Candidate, enz: Enzyme, aa_set: &AminoAcidSet) -> i32 {
    let credit_neighboring = aa_set.neighboring_aa_cleavage_credit();
    let penalty_neighboring = aa_set.neighboring_aa_cleavage_penalty();
    let credit_peptide = aa_set.peptide_cleavage_credit();
    let penalty_peptide = aa_set.peptide_cleavage_penalty();
    let enz_is_c_term = enz.is_c_term();
    let enz_is_n_term = enz.is_n_term();

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
        // C-term cleavage (peptide)
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


/// Research diagnostic: the set of observed MS2 peaks claimed by `peptide`'s
/// charge-1 b/y ions, as quantized m/z keys (round(mz·1000)). Mirrors the
/// matching in `compute_psm_features`. Used only by the env-gated chimeric
/// fragment-overlap diagnostic; not on any production path.
pub(crate) fn matched_peak_keys(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
) -> std::collections::HashSet<i64> {
    let mut keys = std::collections::HashSet::new();
    let n = peptide.length();
    if n < 2 {
        return keys;
    }
    let predicted = predict_by_ions(peptide, 1..=1);
    let feat_tol = scorer.feature_match_tolerance();
    for p in &predicted {
        let tol_da = feat_tol.as_da(p.mz);
        if let Some((_rank, _intensity, peak_mz)) = scored_spec.nearest_peak_full(p.mz, tol_da) {
            keys.insert((peak_mz * 1000.0).round() as i64);
        }
    }
    keys
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
    intensity_model: Option<&IntensityModel>,
) -> PsmFeatures {
    let n = peptide.length();
    if n < 2 {
        return PsmFeatures::default();
    }

    // ADDITIVE edge-score feature (new PIN column). Computed
    // here so it shares the per-PSM ScoredSpectrum + scorer references that
    // the existing feature-extraction code already has on hand.
    let edge_score = psm_edge_score(scored_spec, peptide, scorer, charge);

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
    let predict_losses = scorer.param().data_type.activation.predicts_neutral_losses();
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

    // Percolator feature counting uses a hardcoded fragment tolerance, NOT
    // param.mme (Kim et al., Nat Commun 5:5277, 2014). High-res instruments
    // (HighRes / TOF / QExactive) get 20 ppm; low-res LTQ gets 0.5 Da.
    // The param.mme value (0.5 Da for the hcd_qexactive_tryp model) is the
    // coarser binning tolerance used by the rank-distribution tables —
    // appropriate for node-score lookup but ~50× too wide for feature
    // counting at m/z 500.
    let feat_tol = scorer.feature_match_tolerance();
    let feature_tol_is_ppm = matches!(feat_tol, model::tolerance::Tolerance::Ppm(_));
    let feature_tol = feat_tol.raw_value(); // numeric value in the unit (ppm or Da)

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
        if let Some((rank, intensity, peak_mz)) =
            scored_spec.nearest_peak_full(p.mz, tol_da)
        {
            let is_b = matches!(p.kind, IonKind::B);
            matched_ions.push((intensity, peak_mz, p.mz, is_b));
            matched_ion_ranks.push(rank);

            // position is 1-based (b1/y1 = index 0 in the matched arrays)
            let pos = (p.position - 1) as usize;
            match p.kind {
                IonKind::B => {
                    if pos < b_matched.len() {
                        b_matched[pos] = true;
                        // Record the best (lowest) INTACT-ion rank per position
                        // for the complementary-balance feature; skip loss ions.
                        if !p.is_loss() && rank < b_rank[pos] {
                            b_rank[pos] = rank;
                        }
                    }
                }
                IonKind::Y => {
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
    let segment_ions: SmallVec<[&[scoring_crate::param_model::IonType]; 8]> =
        (0..num_segments)
            .map(|seg| scorer.param().ion_types_for_partition_slice(charge, parent_mass, seg))
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
                    scoring_crate::param_model::IonType::Prefix { charge: ic, offset_bits, .. } => {
                        let offset = f32::from_bits(offset_bits) as f64;
                        let z = ic as f64;
                        (true, (prm_accurate / z + offset, ion))
                    }
                    scoring_crate::param_model::IonType::Suffix { charge: ic, offset_bits, .. } => {
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
        if denom > 0.0 { (num / denom) as f32 } else { 0.0 }
    };

    let explained_ion_current_ratio = safe_div(matched_total, total_intensity);
    let n_term_ion_current_ratio    = safe_div(matched_b_intensity, total_intensity);
    let c_term_ion_current_ratio    = safe_div(matched_y_intensity, total_intensity);
    // MS2 ion current is the raw sum (no log10 transform).
    let ms2_ion_current = if total_intensity > 0.0 { total_intensity as f32 } else { 0.0 };
    // Isolation-window efficiency is not available → emit 0.0.
    let isolation_window_efficiency = 0.0_f32;

    // ── Top-7 mass-error statistics ───────────────────────────────────────────

    // Sort matched ions descending by intensity.
    matched_ions.sort_by(|a, b| {
        b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
    });
    let top7 = &matched_ions[..matched_ions.len().min(7)];

    // All four *ErrorTop7 columns are in ppm (Kim et al., Nat Commun 5:5277,
    // 2014 Percolator features). MeanErrorTop7 = mean of |ppm error|;
    // MeanRelErrorTop7 = mean of signed ppm error. The "Rel" suffix
    // distinguishes signed vs absolute, not Da-vs-ppm.
    //
    // Population stdev formula: sqrt(sum_sq/n - mean²).
    let abs_ppm_errors: Vec<f64> = top7.iter()
        .filter(|&&(_, _, pred, _)| pred > 0.0)
        .map(|&(_, obs, pred, _)| ((obs - pred) / pred * 1e6).abs())
        .collect();
    let rel_ppm_errors: Vec<f64> = top7.iter()
        .filter(|&&(_, _, pred, _)| pred > 0.0)
        .map(|&(_, obs, pred, _)| (obs - pred) / pred * 1e6)
        .collect();

    fn mean_and_pop_stdev(values: &[f64]) -> (f32, f32) {
        if values.is_empty() { return (0.0, 0.0); }
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let sum_sq: f64 = values.iter().map(|v| v * v).sum();
        let var = (sum_sq / n - mean * mean).max(0.0); // clamp negative rounding noise
        (mean as f32, var.sqrt() as f32)
    }

    let (mean_error_top7, stdev_error_top7)         = mean_and_pop_stdev(&abs_ppm_errors);
    let (mean_rel_error_top7, stdev_rel_error_top7) = mean_and_pop_stdev(&rel_ppm_errors);

    // ── Strong-score Stage-1: ppm-Gaussian "tight-match evidence" ─────────────
    // Σ exp(-½ (ppmᵢ/σ)²) over ALL matched ions (not just top-7), σ = 7 ppm.
    // Each ion contributes ~1 when matched within a few ppm and ~0 near the
    // tolerance edge, so the sum rewards PSMs whose fragments cluster tightly in
    // mass accuracy — the high-res signal the intensity-rank model throws away.
    // For low-res spectra every ppm error is large, so the term is ~0 for all
    // candidates and Percolator simply learns a near-zero weight (harmless).
    const PPM_SIGMA: f64 = 7.0;
    let ppm_gaussian_score: f32 = matched_ions
        .iter()
        .filter(|&&(_, _, pred, _)| pred > 0.0)
        .map(|&(_, obs, pred, _)| {
            let ppm = (obs - pred) / pred * 1e6;
            (-0.5 * (ppm / PPM_SIGMA) * (ppm / PPM_SIGMA)).exp()
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
            let tol_da = if feature_tol_is_ppm { pred * feature_tol / 1e6 } else { feature_tol };
            let rho = scored_spec.local_peak_density(obs, DENSITY_HW);
            let p_chance = rho * (2.0 * tol_da);
            if p_chance > 0.0 { (-p_chance.ln()).max(0.0) } else { 0.0 }
        })
        .sum::<f64>() as f32;

    // ── Strong-score Stage-1: complementary ladder ─────────────────────────
    // Bond i (1-based) is reported by both b_i and y_{n-i}; track cleavage
    // sites where BOTH halves are observed. `b_any_matched`/`y_any_matched` are
    // 0-indexed by ion number, so b_i = [i-1] and y_{n-i} = [n-i-1].
    let mut complementary_sites: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    for i in 1..n {
        complementary_sites[i - 1] =
            b_any_matched[i - 1] && y_any_matched[n - i - 1];
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
        matched_ion_ranks.iter().map(|&r| r as f32).sum::<f32>()
            / matched_ion_ranks.len() as f32
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
                    .filter(|&&theo| {
                        (theo - obs).abs() <= tol_da && (theo - pred).abs() > 1e-9
                    })
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
    // NCE is not carried on Spectrum today; use "unknown" (model backs off).
    // When a v3 frag-intensity regressor is present on the param it takes
    // precedence; the coarse IntensityModel table is the fallback.
    let frag_intensity_model = scorer.param().frag_intensity_model.as_deref();
    let intensity_signal_val = intensity_signal(
        intensity_model,
        frag_intensity_model,
        scored_spec,
        peptide,
        charge,
        "unknown",
        feature_tol,
        feature_tol_is_ppm,
    );

    // Tier-2 frag-intensity LLR battery (additive PIN features; 0.0 when no
    // frag model). Discriminative alternative to the cosine intensity_signal.
    let (frag_pred_explained, frag_pred_chance_llr, frag_topk_observed) = frag_llr_battery(
        frag_intensity_model,
        scored_spec,
        peptide,
        charge,
        feature_tol,
        feature_tol_is_ppm,
    );

    // Decoy-aware rich-ion LLR (additive PIN feature; 0.0 when no rich-ion model).
    let rich_ion_llr_val = rich_ion_llr(
        scorer.param().rich_ion_model.as_deref(),
        scored_spec,
        peptide,
        charge,
        feature_tol,
        feature_tol_is_ppm,
    );

    // Mod-localization site-determining-ion features (additive PIN columns; all
    // 0.0 for an unmodified peptide). Uses the SAME feature tolerance as the
    // RichIonLLR / ion-feature path above (train/serve parity).
    let mod_site = mod_site_features(peptide, scored_spec, feature_tol, feature_tol_is_ppm);

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
    }
}

/// Reject combinations the mmap (out-of-core) candidate backing cannot serve.
///
/// The chimeric two-pass cascade (`run_pass2_coisolation` → `search_secondary`)
/// resolves co-isolated peptides against the IN-RAM `candidates`/`bucket_index`,
/// which are empty under mmap mode (candidate resolution goes through
/// `mmap_index` instead). Running both together silently yields no secondaries,
/// so the combination is rejected up front in [`PreparedSearch::prepare_mmap`].
fn reject_unsupported_mmap_combo(params: &SearchParams) -> std::io::Result<()> {
    if params.chimeric {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "--chimeric is not supported together with --candidate-index mmap: the \
             chimeric Pass-2 search resolves co-isolated peptides against the \
             in-memory candidate index, which is not populated in out-of-core mmap \
             mode. Use either --chimeric or --candidate-index mmap, not both.",
        ));
    }
    Ok(())
}

// ── Unit tests for feature columns ───────────────────────────────────────────

#[cfg(test)]
mod mmap_combo_tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;

    #[test]
    fn mmap_rejects_chimeric() {
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut p = SearchParams::default_tryptic(aa);
        p.chimeric = true;
        let err = reject_unsupported_mmap_combo(&p).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
    }

    #[test]
    fn mmap_allows_non_chimeric() {
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        let p = SearchParams::default_tryptic(aa);
        assert!(reject_unsupported_mmap_combo(&p).is_ok());
    }
}

#[cfg(test)]
mod feature_tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::mass::PROTON;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;
    use scoring_crate::scoring::fragment_ions::predict_by_ions;
    use scoring_crate::scoring::ScoredSpectrum;
    use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Partition, SpecDataType};
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use model::protocol::Protocol;
    use model::tolerance::Tolerance;
    use rustc_hash::FxHashMap;

    /// Minimal RankScorer for feature tests, with mme = Da(tol_da).
    ///
    /// Uses realistic prefix/suffix offsets so the partition-ion-list
    /// intensity-ratio path matches peaks placed at `predict_by_ions`'s
    /// standard b/y m/z values (b_neutral + PROTON; y_neutral = suffix +
    /// H2O + PROTON). The partition ion list is read directly, so the
    /// offsets matter (offset=0.0 with no suffix ion would not work here).
    fn make_scorer(tol_da: f64) -> RankScorer {
        use model::mass::{H2O, PROTON};
        let part = Partition { charge: 2, parent_mass: 0.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: (PROTON as f32).to_bits(), loss_class: 0 };
        let suffix1 = IonType::Suffix { charge: 1, offset_bits: ((H2O + PROTON) as f32).to_bits(), loss_class: 0 };
        let noise = IonType::Noise;
        let mut ion_table = FxHashMap::default();
        ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
        ion_table.insert(suffix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
        ion_table.insert(noise, vec![0.1_f32, 0.2, 0.3, 0.4]);
        let mut rank_dist_table = FxHashMap::default();
        rank_dist_table.insert(part, ion_table);
        let mut frag_off_table = FxHashMap::default();
        frag_off_table.insert(part, vec![
            FragmentOffsetFrequency { ion_type: prefix1, frequency: 0.7 },
            FragmentOffsetFrequency { ion_type: suffix1, frequency: 0.7 },
        ]);
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
        }
    }

    // ── Test: empty spectrum → all new features are 0 ───────────────────────

    #[test]
    fn compute_psm_features_top7_error_stats_zero_when_no_matches() {
        let pep = ala_peptide(4);
        let spec = make_spectrum(vec![]); // no peaks
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.5), 2, None);
        assert_eq!(f.mean_error_top7,     0.0, "mean_error_top7 should be 0 with no matches");
        assert_eq!(f.stdev_error_top7,    0.0, "stdev_error_top7 should be 0 with no matches");
        assert_eq!(f.mean_rel_error_top7,  0.0, "mean_rel_error_top7 should be 0 with no matches");
        assert_eq!(f.stdev_rel_error_top7, 0.0, "stdev_rel_error_top7 should be 0 with no matches");
        assert_eq!(f.explained_ion_current_ratio, 0.0, "ratio should be 0 with no peaks");
        assert_eq!(f.ms2_ion_current, 0.0, "ms2_ion_current should be 0 with no peaks");
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
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.01), 2, None); // tight tolerance

        // All ratios should be positive since all predicted ions match.
        assert!(f.explained_ion_current_ratio > 0.0,
            "explained_ion_current_ratio should be > 0 when ions match, got {}",
            f.explained_ion_current_ratio);
        assert!(f.n_term_ion_current_ratio > 0.0,
            "n_term_ion_current_ratio should be > 0 when b-ions match");
        assert!(f.c_term_ion_current_ratio > 0.0,
            "c_term_ion_current_ratio should be > 0 when y-ions match");

        // Invariant: n_term + c_term == explained (within float precision)
        let sum = f.n_term_ion_current_ratio + f.c_term_ion_current_ratio;
        assert!(
            (sum - f.explained_ion_current_ratio).abs() < 1e-5,
            "n_term + c_term should == explained ({} + {} != {})",
            f.n_term_ion_current_ratio, f.c_term_ion_current_ratio, f.explained_ion_current_ratio
        );

        // ms2_ion_current should equal total peak intensity sum.
        let total: f32 = ss.total_intensity() as f32;
        assert!((f.ms2_ion_current - total).abs() < 1.0,
            "ms2_ion_current {} should match total spectrum intensity {}",
            f.ms2_ion_current, total);

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
            &pep, &scorer, 2, None,
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
            &pep2, &scorer, 2, None,
        );
        assert!(
            f_contested.unique_match_fraction < f_unique.unique_match_fraction,
            "contested peak should lower fraction (unique={}, contested={})",
            f_unique.unique_match_fraction, f_contested.unique_match_fraction
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
        assert!(!predicted_z2.is_empty(), "charge-2 ions should exist for 4-mer");

        // Empty spectrum → count 0.
        let f_empty = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(vec![])),
            &pep, &make_scorer(0.01), 3, None,
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
            &pep, &make_scorer(0.01), 3, None,
        );
        assert!(
            f.doubly_charged_matched_ion_count > 0,
            "charge-2 peaks should match, got {}",
            f.doubly_charged_matched_ion_count
        );
        assert_eq!(
            f.doubly_charged_matched_ion_count, predicted_z2.len() as u32,
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
            &pep, &make_scorer(0.01), 2, None,
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
            &pep, &make_scorer(0.01), 2, None,
        );
        assert!(
            f_buried.mean_matched_intensity_rank > f_top.mean_matched_intensity_rank,
            "buried ions should have higher mean rank (top={}, buried={})",
            f_top.mean_matched_intensity_rank, f_buried.mean_matched_intensity_rank
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
            &pep, &make_scorer(0.01), 2, None,
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
            &pep, &make_scorer(0.01), 2, None,
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
            f_split.longest_complementary_ladder, f_all.longest_complementary_ladder
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
            &pep, &make_scorer(0.01), 2, None,
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
            &pep, &make_scorer(0.01), 2, None,
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
        let mut peaks = vec![
            (b2.mz, 100.0_f32),
            (b2.mz - 18.0106, 50.0_f32),
        ];
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let f = compute_psm_features(
            &ScoredSpectrum::new_without_filtering(&make_spectrum(peaks)),
            &pep, &make_scorer(0.01), 2, None,
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
            &pep, &make_scorer(0.01), 2, None,
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
            &pep, &make_scorer(0.01), 2, None,
        );
        assert!(
            f_shift.ppm_gaussian_score < f_exact.ppm_gaussian_score * 0.5,
            "ppm_gaussian_score should fall sharply when matches scatter (exact={}, shifted={})",
            f_exact.ppm_gaussian_score, f_shift.ppm_gaussian_score
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
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.05), 2, None);

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
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.5), 2, None);

        let expected: f32 = peaks.iter().map(|&(_, i)| i).sum();
        assert_eq!(f.ms2_ion_current, expected,
            "ms2_ion_current {} should equal sum of peak intensities {}",
            f.ms2_ion_current, expected);
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
        let b1 = predicted.iter().find(|p| matches!(p.kind, IonKind::B) && p.position == 1)
            .expect("b1 ion should exist");
        assert!(
            (b1.mz - expected_b1_mz).abs() < 1e-6,
            "b1 mz {} expected {}", b1.mz, expected_b1_mz
        );
    }
}

/// Pre-merge dedup pass: collapse PSMs sharing the same peptide sequence
/// and rounded rank_score before per-charge merge.
pub(crate) fn dedup_pepseq_score(
    psms: Vec<PsmMatch>,
    candidates: &[Candidate],
) -> Vec<PsmMatch> {
    use std::collections::btree_map::Entry;
    use std::collections::BTreeMap;

    let mut pep_key_cache: FxHashMap<u32, Arc<PepDedupKey>> = FxHashMap::default();
    let mut groups: BTreeMap<DedupMapKey, PsmMatch> = BTreeMap::new();

    for psm in psms {
        let primary = psm.primary_candidate_idx();
        let pep_key = pep_key_cache
            .entry(primary)
            .or_insert_with(|| Arc::new(PepDedupKey::from_peptide(&candidates[primary as usize].peptide)))
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

/// Number of raw base-peptide records in the mmap index that share `cand`'s
/// `(protein_index, start_offset, length, decoy/n-term/c-term flags)` identity —
/// i.e. how many times the in-RAM `enumerate_candidates` emits this base peptide
/// at these exact coordinates (≥1; >1 only for the N-terminal-Met
/// re-enumeration that produces identical-coordinate spans). Used by the `Mmap`
/// path to replicate each lazily-expanded peptidoform so the per-spectrum
/// candidate set (and thus the `dedup_pepseq_score`-aggregated `candidate_idxs`,
/// hence the PIN `Proteins` column) matches the in-RAM path copy-for-copy.
///
/// The base mass (unmodified residue masses + water) is recovered from the
/// candidate by subtracting each residue's mod delta; the index is queried at
/// that exact `mass_milli`.
fn base_record_multiplicity(mi: &MmapCandidateIndex, cand: &Candidate) -> u32 {
    use crate::candidate_index::flags as rec_flags;
    // Base mass = peptide neutral mass with every mod delta removed.
    let mut base_mass = cand.peptide.mass();
    for aa in &cand.peptide.residues {
        if let Some(m) = aa.mod_.as_ref() {
            base_mass -= m.mass_delta;
        }
    }
    let mass_milli = (base_mass * 1000.0).round() as u64;
    let mut want_flags: u16 = 0;
    if cand.is_decoy {
        want_flags |= rec_flags::IS_DECOY;
    }
    if cand.is_protein_n_term {
        want_flags |= rec_flags::IS_PROTEIN_N_TERM;
    }
    if cand.is_protein_c_term {
        want_flags |= rec_flags::IS_PROTEIN_C_TERM;
    }
    let pidx = cand.protein_index as u32;
    let start = cand.start_offset_in_protein as u32;
    let len = cand.peptide.residues.len() as u16;
    // Tolerate ±1 milli rounding between `peptide.mass()` here and the index's
    // stored `mass_milli` (built from the same formula, but de-mod arithmetic can
    // drift one unit at the rounding boundary).
    let count = mi
        .mass_window(mass_milli.saturating_sub(1), mass_milli + 1)
        .into_iter()
        .filter(|r| {
            r.protein_index == pidx
                && r.start_offset == start
                && r.length == len
                && r.flags == want_flags
        })
        .count() as u32;
    count.max(1)
}

/// Bare residue sequence (no mods, no flanks) of a candidate's peptide.
fn bare_residues(cand: &Candidate) -> Box<[u8]> {
    cand.peptide.residues.iter().map(|aa| aa.residue).collect()
}

/// Relabel decoy candidates whose bare peptide sequence is ALSO a target peptide
/// as targets. Such "decoys" are reversal/shuffle coincidences (or palindromes /
/// short peptides) that reproduce a real target sequence; counting them as decoys
/// over-counts the decoy space (TDC conservative) and, for decoy-aware training,
/// mislabels true-sequence ions as negatives. Relabeling to target lets the
/// per-spectrum dedup merge the collision with its target twin — no EXTRA target
/// is created, because the sequence is a real target peptide so the twin already
/// exists. Bare residues only (mods don't change sequence identity).
///
/// Cost: one HashSet of target peptide sequences, built once at index setup.
fn relabel_collision_decoys(candidates: &mut [Candidate]) {
    use std::collections::HashSet;
    let target_seqs: HashSet<Box<[u8]>> = candidates
        .iter()
        .filter(|c| !c.is_decoy)
        .map(bare_residues)
        .collect();
    relabel_collision_decoys_with(candidates, &target_seqs);
}

/// Relabel collision decoys against a PRE-COMPUTED target bare-sequence set.
///
/// Used by the out-of-core (`Mmap`) path, where the target set must be the GLOBAL
/// set of all DB targets (built once from the full index), NOT the targets that
/// happen to land in the current precursor window. This makes the relabel decision
/// identical to the in-RAM path, which relabels once globally over all candidates.
/// (`relabel_collision_decoys` is the in-RAM convenience wrapper that derives the
/// set from the candidate slice itself — correct only when that slice IS the full
/// candidate list.)
fn relabel_collision_decoys_with(
    candidates: &mut [Candidate],
    target_seqs: &std::collections::HashSet<Box<[u8]>>,
) {
    if target_seqs.is_empty() {
        return;
    }
    for cand in candidates.iter_mut().filter(|c| c.is_decoy) {
        if target_seqs.contains(&bare_residues(cand)) {
            cand.is_decoy = false;
        }
    }
}

/// Build the GLOBAL set of target bare-peptide sequences from the full mmap
/// base-peptide index. For every NON-decoy record, reconstruct the residue span
/// from the protein sequence (`protein.sequence[start..start+length]`). This is
/// the set the in-RAM `relabel_collision_decoys` derives from the complete
/// candidate enumeration — built here directly from the index so the `Mmap` path
/// makes the SAME global relabel decision without materializing all candidates.
///
/// Bare residues only (no mods, no flanks) → bounded memory: one `Box<[u8]>` per
/// distinct target backbone, far smaller than the candidate `Vec` the out-of-core
/// mode exists to avoid.
fn build_target_bare_seqs(
    mi: &MmapCandidateIndex,
    db: &SearchIndex,
) -> std::collections::HashSet<Box<[u8]>> {
    use crate::candidate_index::flags as rec_flags;
    let mut set: std::collections::HashSet<Box<[u8]>> = std::collections::HashSet::new();
    for rec in mi.records() {
        if rec.flags & rec_flags::IS_DECOY != 0 {
            continue;
        }
        let pidx = rec.protein_index as usize;
        if pidx >= db.db.proteins.len() {
            continue;
        }
        let seq = &db.db.proteins[pidx].sequence;
        let start = rec.start_offset as usize;
        let end = start + rec.length as usize;
        if end > seq.len() {
            continue;
        }
        set.insert(seq[start..end].to_vec().into_boxed_slice());
    }
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
        Self { residues, mod_units }
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
    use std::sync::Arc;
    use model::amino_acid::AminoAcid;
    use model::modification::{ModLocation, Modification};
    use model::peptide::Peptide;
    use model::ResidueSpec;
    use crate::psm::PsmMatch;

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
        relabel_collision_decoys(&mut cands);
        assert!(!cands[0].is_decoy, "target unchanged");
        assert!(!cands[1].is_decoy, "collision decoy relabeled to target");
        assert!(cands[2].is_decoy, "genuine decoy left as decoy");
    }

    #[test]
    fn dedup_uses_rank_score_not_pin_score() {
        let pep = seq_peptide(b"PEPTK");
        let cands = vec![cand_with_peptide(pep.clone())];
        let psms = vec![
            psm(0, 100.4, 99.6),
            psm(0, 120.0, 99.6),
        ];
        let out = dedup_pepseq_score(psms, &cands);
        assert_eq!(out.len(), 2, "different rank_score keys must not merge");
    }

    #[test]
    fn dedup_distinguishes_mod_state() {
        let mut ox = seq_peptide(b"PEPMK");
        ox.residues[3].mod_ = Some(Arc::new(Modification {
            name: "Ox".into(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        }));
        let cands = vec![
            cand_with_peptide(seq_peptide(b"PEPMK")),
            cand_with_peptide(ox),
        ];
        let psms = vec![
            psm(0, 50.0, 50.0),
            psm(1, 50.0, 50.0),
        ];
        let out = dedup_pepseq_score(psms, &cands);
        assert_eq!(out.len(), 2, "mod-aware pepSeq keys must differ");
    }

    #[test]
    fn dedup_keeps_highest_rank_score_survivor() {
        let pep = seq_peptide(b"PEPTK");
        let cands = vec![cand_with_peptide(pep)];
        // Same rounded score bucket (60) but different float rank_score — merge to best.
        let psms = vec![
            psm(0, 59.6, 50.0),
            psm(0, 60.4, 50.0),
        ];
        let out = dedup_pepseq_score(psms, &cands);
        assert_eq!(out.len(), 1);
        assert!((out[0].rank_score - 60.4).abs() < f32::EPSILON);
    }
}
