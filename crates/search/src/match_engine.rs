//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use std::collections::{BTreeMap, HashMap};
use std::hash::Hasher;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

// GF failure-mode diagnostics. Module-level atomics
// incremented per-bin from compute_spec_e_values_for_spectrum and
// reported in the yield-accounting summary. Used to characterise the
// ~4.7% of Astral PSMs where GF compute fails. Module-level rather than
// per-PreparedSearch because we want cumulative counts across all
// chunks and the per-call wiring would be invasive.
//
// These are diagnostics-only; behavior is unchanged. They are reset at
// the start of each run_chunk invocation so per-bench numbers don't
// accumulate across calls.
static GF_EMPTY_SCORE_RANGE: AtomicU64 = AtomicU64::new(0);
static GF_SINK_UNREACHABLE: AtomicU64 = AtomicU64::new(0);
static GF_SINK_RETRY_OK: AtomicU64 = AtomicU64::new(0);
static GF_BIN_ATTEMPTS: AtomicU64 = AtomicU64::new(0);
static GF_SPECTRA_NO_GROUP: AtomicU64 = AtomicU64::new(0);

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet, FxHasher};
use smallvec::{smallvec, SmallVec};

use model::aa_set::AminoAcidSet;
use input::Ms1Link;
use crate::candidate_gen::{enumerate_candidates, Candidate};
use model::enzyme::Enzyme;
use scoring_crate::gf::generating_function::GeneratingFunction;
use scoring_crate::gf::group::GeneratingFunctionGroup;
use scoring_crate::gf::primitive_graph::PrimitiveAaGraph;
use model::mass::{nominal_from, H2O, PROTON};
use model::peptide::Peptide;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::precursor_matching::{matches_precursor, MassError};
use crate::psm::{PsmFeatures, PsmMatch, TopNQueue};
use scoring_crate::scoring::fragment_ions::{IonKind, predict_by_ions};
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use scoring_crate::scoring::{psm_edge_score, score_psm, RankScorer, ScoredSpectrum};
use model::spectrum::Spectrum;

/// One-time-built state shared across every chunk of a streamed search.
///
/// `match_spectra` materializes its full set of candidates, bucket index,
/// distinct-peptide counts, and enzyme-registered aa_set in a single pass at
/// startup. For chunked / streaming spectrum loading we want to reuse that
/// state instead of rebuilding it per chunk. `PreparedSearch::prepare` does
/// the setup once; `PreparedSearch::run_chunk` runs the per-spectrum scoring
/// loop on any slice of `Spectrum`s using that prepared state.
///
/// The two-pass split mirrors the original `match_spectra` body — there is
/// no algorithmic change. Pre-existing single-call callers can still use
/// `match_spectra(...)` which is now a thin wrapper around
/// `prepare` + a single `run_chunk` call.
pub struct PreparedSearch<'a> {
    pub idx: &'a SearchIndex,
    pub params: &'a SearchParams,
    pub scorer: &'a RankScorer,
    pub fragment_tolerance_da: f64,
    /// Final, deduplicated candidate list (target + decoy).
    pub candidates: Vec<Candidate>,
    /// `nominal(peptide.mass() - H2O)` → indices into `candidates`.
    pub bucket_index: BTreeMap<i32, Vec<usize>>,
    /// `params.aa_set` with the search enzyme registered for GF cleavage
    /// scoring. Cheap to clone, but we keep one shared copy here.
    pub aa_set_for_gf: AminoAcidSet,
    /// Optional MS1 linkage for the chimeric precursor isotope features.
    /// `None` unless the binary supplied one via
    /// [`Self::with_ms1_link`] under `--chimeric`. When `None` (the default
    /// and the entire `--chimeric off` path), the feature fill skips the
    /// precursor-envelope computation and the two new `PsmFeatures` fields
    /// stay 0.0 — keeping the off path bit-identical.
    pub ms1_link: Option<Ms1Link>,
}

/// Owned, precursor-tolerance-independent products of [`PreparedSearch::prepare`]
/// (candidates, bucket index, GF aa_set). Built once during the
/// calibration pre-pass and reused for the tightened-tolerance main pass, avoiding
/// a second full candidate enumeration (~15s on the 16.8M-candidate Astral search).
pub struct PreparedParts {
    candidates: Vec<Candidate>,
    bucket_index: BTreeMap<i32, Vec<usize>>,
    aa_set_for_gf: AminoAcidSet,
}

/// Derive the inclusive `[min_nominal, max_nominal]` nominal-mass bucket bounds
/// for one spectrum at one charge state `z`, used to enumerate candidate
/// peptides from the mass-bucket index.
///
/// The window is centered on the precursor m/z and widened by the precursor
/// tolerance plus the isotope-error range. This is the original, byte-for-byte
/// unchanged derivation used by both the standard search and the cascade's
/// NARROW Pass 1.
fn candidate_nominal_bounds(
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
    // Convention: max widens by tol_da_left, min widens by tol_da_right.
    let min_nominal = nominal_center - iso_max - widen_right;
    let max_nominal = nominal_center - iso_min + widen_left;
    (min_nominal, max_nominal)
}

impl<'a> PreparedSearch<'a> {
    /// Build the per-search state once. Enumerates candidates, builds the
    /// mass-bucket index, seeds the `SearchIndex` distinct-peptide counts,
    /// and clones+registers the aa_set for GF cleavage scoring.
    pub fn prepare(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        decoy_prefix: &str,
    ) -> Self {
        // Collect the production candidate list AND seed the per-length
        // distinct-peptide counts in a single pass. This avoids a second full
        // `enumerate_candidates(...)` walk just to populate the E-value
        // denominator map.
        let mut candidates: Vec<Candidate> = Vec::new();
        let mut seen_per_length: HashMap<usize, FxHashSet<u64>> = HashMap::new();
        for cand in enumerate_candidates(idx, params, decoy_prefix) {
            let residues = &cand.peptide.residues;
            let mut h = FxHasher::default();
            for aa in residues {
                h.write_u8(aa.residue);
            }
            seen_per_length
                .entry(residues.len())
                .or_default()
                .insert(h.finish());
            candidates.push(cand);
        }
        let distinct_counts: HashMap<usize, usize> = seen_per_length
            .into_iter()
            .map(|(len, set)| (len, set.len()))
            .collect();
        idx.set_distinct_peptide_counts_if_absent(distinct_counts);

        // Build mass-bucket index: nominal(peptide.mass() - H2O) → Vec<candidate_idx>.
        //
        // Uses the same nominal_from convention as the GF mass-bin loop so that
        // bucket keys align with the GF's mass-bin lookup.
        // Stores only indices into `candidates` — no cloning, tiny memory overhead.
        let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
        for (cand_idx, cand) in candidates.iter().enumerate() {
            let nominal = cand.peptide.nominal_residue_mass();
            bucket_index.entry(nominal).or_default().push(cand_idx);
        }

        // Build an aa_set clone with enzyme registered (for GF cleavage scoring).
        // Defaults: peptide_eff = 0.95, neighboring_eff = 0.95.
        // Cloning is cheap (AminoAcidSet is a HashMap of ~20 entries).
        // This avoids mutating the shared SearchParams.aa_set borrow.
        let mut aa_set_for_gf: AminoAcidSet = params.aa_set.clone();
        if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
            aa_set_for_gf.register_enzyme(params.enzyme, 0.95, 0.95);
        }

        PreparedSearch {
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            candidates,
            bucket_index,
            aa_set_for_gf,
            ms1_link: None,
        }
    }

    /// Consume this prepared search into owned [`PreparedParts`], dropping the
    /// `idx`/`params`/`scorer` borrows so the caller can mutate `params` (e.g.
    /// tighten the precursor tolerance) before rebuilding via [`Self::from_parts`].
    pub fn into_parts(self) -> PreparedParts {
        PreparedParts {
            candidates: self.candidates,
            bucket_index: self.bucket_index,
            aa_set_for_gf: self.aa_set_for_gf,
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
            aa_set_for_gf: parts.aa_set_for_gf,
            ms1_link: None,
        }
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
        let idx = self.idx;
        let fragment_tolerance_da = self.fragment_tolerance_da;
        let candidates = &self.candidates;
        let bucket_index = &self.bucket_index;
        let aa_set_for_gf = &self.aa_set_for_gf;

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
        // inflation. Zero cost unless MSGF_CHIMERIC_OVERLAP is set AND --chimeric.
        let chim_overlap = params.chimeric && std::env::var("MSGF_CHIMERIC_OVERLAP").is_ok();

        // Parallel per-spectrum search. All inputs above are `&` immutable; the
        // closure owns its TopNQueue, scored_per_charge cache, and per-bin GF state.
        let queues: Vec<TopNQueue> = spectra
            .par_iter()
            .enumerate()
            .map(
                |(local_idx, spec)| {
                let spec_idx = local_idx + spectrum_idx_offset;
                let mut queue = TopNQueue::new(params.top_n_psms_per_spectrum);

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
            // set of candidate indices. Window derivation mirrors
            // compute_spec_e_values_for_spectrum's logic so any candidate admitted by
            // matches_precursor is guaranteed to be in at least one charge's window.
            //
            // Vec + sort_unstable + dedup is faster than BTreeSet for the typical
            // 1k-3k indices per spectrum: better cache locality, no tree pointer
            // chasing, single sort pass at end. Iteration order matches BTreeSet
            // (ascending), preserving downstream parity / determinism.
            let mut window_cand_indices: Vec<usize> = Vec::with_capacity(2048);
            let shift_ppm = params.precursor_mass_shift_ppm;
            for &z in &charges_to_try {
                // Cascade Pass 1 stays NARROW: the bounds use the
                // precursor-tolerance mode, not the wide isolation window. Off path
                // and `--chimeric` Pass 1 agree → byte-identical.
                let (min_nominal, max_nominal) =
                    candidate_nominal_bounds(spec, z, params, shift_ppm);
                for (_nm, idxs) in bucket_index.range(min_nominal..=max_nominal) {
                    window_cand_indices.extend_from_slice(idxs);
                }
            }
            window_cand_indices.sort_unstable();
            window_cand_indices.dedup();

            // Hoist the loop-invariant cleavage-credit constants out of the
            // per-candidate hot path: resolving them once here avoids
            // re-invoking four small `aa_set` accessor methods (each a
            // HashMap field deref) for every candidate.
            //
            // The four credit/penalty values are SearchParams-constant; we
            // resolve them ONCE here. The per-candidate logic becomes four
            // branches over precomputed i32 constants.
            let enz_credit_neighboring = aa_set_for_gf.neighboring_aa_cleavage_credit();
            let enz_penalty_neighboring = aa_set_for_gf.neighboring_aa_cleavage_penalty();
            let enz_credit_peptide = aa_set_for_gf.peptide_cleavage_credit();
            let enz_penalty_peptide = aa_set_for_gf.peptide_cleavage_penalty();
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

            // Per-charge queue keyed by charge state. Mirrors Java's
            // per-SpecKey raw-score retention (Java parity).
            let mut per_charge_queues: FxHashMap<u8, TopNQueue> = FxHashMap::default();

            // Cascade Pass 1 is a NARROW brute-force window scan: iterate every
            // candidate whose nominal mass falls in the precursor window.
            for &cand_idx in &window_cand_indices {
                let cand = &candidates[cand_idx];
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
                    // is the Java-aligned queue-ordering key.
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
                        candidate_idxs: vec![cand_idx as u32],
                        charge_used: z,
                        mass_error_ppm: err.mass_error_ppm,
                        score: pin_score,
                        rank_score,
                        edge_score: edge_i,
                        spec_e_value: 1.0,
                        de_novo_score: i32::MIN,
                        activation_method: Some(scorer.param().data_type.activation),
                        e_value: 1.0,
                        features,
                        isotope_offset: err.isotope_offset,
                        precursor_mz_override: None,
                    };
                    per_charge_queues
                        .entry(z)
                        .or_insert_with(|| TopNQueue::new(params.top_n_psms_per_spectrum))
                        .push(psm);
                    psms_pushed.fetch_add(1, Ordering::Relaxed);
                }
            }
            candidates_visited.fetch_add(window_cand_indices.len() as u64, Ordering::Relaxed);

            // pepSeq + score dedup per-charge BEFORE GF compute.
            // Same peptide matched against multiple proteins collapses to one
            // PsmMatch with aggregated candidate_idxs (Java parity for pepSeq dedup).
            for queue in per_charge_queues.values_mut() {
                if queue.len() > 1 {
                    let drained = queue.drain_into_vec();
                    let deduped = dedup_pepseq_score(drained, candidates);
                    for psm in deduped {
                        queue.push(psm);
                    }
                }
            }

            // Per-charge GF / SpecEValue compute. Each per-charge queue
            // gets SpecE calibrated against its OWN charge's GF distribution
            // (Java parity: getRankScorer per SpecKey).
            let enzyme_opt = if params.enzyme != Enzyme::NoCleavage
                && params.enzyme != Enzyme::NonSpecific
            {
                Some(params.enzyme)
            } else {
                None
            };
            let mut any_queue_nonempty = false;
            for (&charge, queue) in per_charge_queues.iter_mut() {
                if queue.is_empty() {
                    continue;
                }
                any_queue_nonempty = true;
                let scored_spec_charge = scored_spec_for_charge(charge);
                compute_spec_e_values_for_spectrum(
                    spec,
                    params,
                    queue,
                    aa_set_for_gf,
                    enzyme_opt,
                    scorer,
                    scored_spec_charge,
                    charge,
                    fragment_tolerance_da,
                    idx,
                    candidates,
                );
            }
            if any_queue_nonempty {
                spectra_with_psms.fetch_add(1, Ordering::Relaxed);
            }

            // Spectrum-level merge with SpecE tie keep. The
            // TopNQueue::push (Ordering::Equal arm) keeps SpecE ties at
            // capacity because PsmMatch::cmp orders by spec_e_value first.
            // Matches Java parity: SpecE tie-keep on spectrum-level merge.
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
            queue.fill_post_topn(|psm| {
                let ss = scored_spec_for_charge(psm.charge_used);
                let cand = &candidates[psm.primary_candidate_idx() as usize];
                let mut features = compute_psm_features(ss, &cand.peptide, scorer, psm.charge_used);
                features.edge_score = psm.edge_score; // reuse per-candidate value

                // The PIN `precursor_isotope_kl` / `precursor_snr` columns are
                // left at their 0.0 defaults here. They were a per-PSM MS1
                // isotope-envelope feature of the old single-pass chimeric mode;
                // the two-pass cascade does not consume them (computing them ran
                // ~121k times on Astral for ~2:40 of wall). The cascade uses MS1
                // only in Pass 2's `run_pass2_coisolation` / `detect_coisolated`
                // co-isolation detection — `ms1_link` and its loading stay.

                features.delta_raw_score = delta_raw_score;
                psm.features = features;
            });

            // Chimeric fragment-overlap diagnostic (env-gated). For scans that
            // emit ≥2 distinct peptides, measure how many MS2 peaks the runner-up
            // claims that the top peptide also claims (the "fragment theft" the
            // chimeric FDR inflation is hypothesized to come from).
            if chim_overlap {
                let sorted = queue.clone().into_sorted_vec(); // best-first
                let mut picks: Vec<&PsmMatch> = Vec::new();
                'outer: for psm in &sorted {
                    let seq: Vec<u8> = candidates[psm.primary_candidate_idx() as usize]
                        .peptide.residues.iter().map(|a| a.residue).collect();
                    for p in &picks {
                        let pseq: Vec<u8> = candidates[p.primary_candidate_idx() as usize]
                            .peptide.residues.iter().map(|a| a.residue).collect();
                        if pseq == seq { continue 'outer; }
                    }
                    picks.push(psm);
                    if picks.len() == 2 { break; }
                }
                if picks.len() == 2 {
                    let pa = &candidates[picks[0].primary_candidate_idx() as usize].peptide;
                    let pb = &candidates[picks[1].primary_candidate_idx() as usize].peptide;
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
        // GF DP failure-mode diagnostics.
        // Cumulative across all chunks in this run; not reset between
        // chunks. Helps localize the ~4.7% Astral PSMs with sentinel
        // DeNovoScore / lnSpecEValue=0 (GF failed for that spectrum's
        // entire precursor-mass window).
        eprintln!(
            "GF diagnostics (cumulative): {} bin attempts, {} EmptyScoreRange, \
             {} SinkUnreachable, {} of those recovered by unthresholded retry, \
             {} spectra with no successful bin",
            GF_BIN_ATTEMPTS.load(Ordering::Relaxed),
            GF_EMPTY_SCORE_RANGE.load(Ordering::Relaxed),
            GF_SINK_UNREACHABLE.load(Ordering::Relaxed),
            GF_SINK_RETRY_OK.load(Ordering::Relaxed),
            GF_SPECTRA_NO_GROUP.load(Ordering::Relaxed),
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

    // The targeted secondary search needs the enzyme only when it actually
    // constrains cleavage; NoCleavage / NonSpecific carry no cleavage credit.
    let enzyme = if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
        Some(params.enzyme)
    } else {
        None
    };

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
            // max_kl: averagine-envelope KL gate. 0.3 requires a clean MS1 envelope
            // → fewer spurious secondaries → FDP toward nominal (small PSM cost).
            0.3,
            2,
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
                &prepared.aa_set_for_gf,
                enzyme,
                params,
                prepared.idx,
                prepared.fragment_tolerance_da,
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

/// For a single spectrum, compute the GF across the precursor tolerance
/// window in nominal mass space, then assign `spec_e_value` to every PSM
/// in `queue` whose nominal_peptide_mass falls within the window.
///
/// # Arguments
/// * `spec` — the spectrum (used for precursor m/z).
/// * `params` — search params (precursor_tolerance, isotope_error_range).
/// * `queue` — the PSM queue for this spectrum (mutated in place).
/// * `aa_set` — amino acid set with enzyme already registered via `register_enzyme`.
/// * `enzyme` — the search enzyme (passed to PrimitiveAaGraph; may be None).
/// * `scorer` — RankScorer.
/// * `scored_spec` — ScoredSpectrum built with `top_charge` (per-charge cache).
/// * `top_charge` — charge of the top PSM in the queue; used for GF mass window.
///   For charge-explicit spectra this equals `spec.precursor_charge.unwrap()`.
///   For charge-missing spectra, using the top PSM's charge ensures the GF
///   reflects the dominant scoring context.
/// * `fragment_tolerance_da` — fragment mass tolerance in Da.
/// * `search_index` — database (target+decoy); used to look up protein sequences
///   for protein-terminal flag derivation.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_spec_e_values_for_spectrum(
    spec: &Spectrum,
    params: &SearchParams,
    queue: &mut TopNQueue,
    aa_set: &AminoAcidSet,
    enzyme: Option<Enzyme>,
    scorer: &RankScorer,
    scored_spec: &ScoredSpectrum<'_>,
    top_charge: u8,
    fragment_tolerance_da: f64,
    search_index: &SearchIndex,
    candidates: &[Candidate],
) {
    // 1. Determine the peptide neutral mass and its tolerance window.
    // For charge-explicit spectra, `top_charge` == spec.precursor_charge.unwrap().
    // For charge-missing spectra, `top_charge` is the top PSM's charge (B3 fix).
    let charge = top_charge;
    if charge == 0 {
        return;
    }

    // peptide_neutral_mass = (precursor_mz - H) * charge - H2O
    // This matches Java: scoredSpec.getPrecursorPeak().getMass() - H2O
    // where getPrecursorPeak().getMass() = (mz - H) * charge.
    let shift_ppm = params.precursor_mass_shift_ppm;
    let peptide_neutral_mass = adjusted_observed_neutral_mass(
        (spec.precursor_mz - PROTON) * (charge as f64) - H2O,
        shift_ppm,
    );
    let nominal_peptide_mass = nominal_from(peptide_neutral_mass);

    // Isotope error convention: range [min_iso, max_iso] is applied as
    //   minNominalPeptideMass = nominalPeptideMass - maxIsotopeError
    //   maxNominalPeptideMass = nominalPeptideMass - minIsotopeError
    let iso_min = *params.isotope_error_range.start() as i32;
    let iso_max = *params.isotope_error_range.end() as i32;
    let min_iso_nominal = nominal_peptide_mass - iso_max;
    let max_iso_nominal = nominal_peptide_mass - iso_min;

    // Tolerance widening: round(tol_da - 0.4999).
    // tol_da_left governs the upper bound; tol_da_right governs the lower bound.
    let tol_da_left = params.precursor_tolerance.left.as_da(peptide_neutral_mass);
    let tol_da_right = params.precursor_tolerance.right.as_da(peptide_neutral_mass);
    let widen_left = (tol_da_left - 0.4999_f64).round() as i32;
    let widen_right = (tol_da_right - 0.4999_f64).round() as i32;

    let max_peptide_mass_idx = max_iso_nominal + widen_left;
    let min_peptide_mass_idx = min_iso_nominal - widen_right;

    if max_peptide_mass_idx < min_peptide_mass_idx {
        return;
    }

    // 2. Compute the minimum score across all PSMs (used as GF score threshold).
    //
    // Use `rank_score` (= node + cleavage + edge), not `score`
    // (= node + cleavage only). Java parity: `match.score` is
    // `cleavageScore + rawScore` where `rawScore` is `DBScanScorer.getScore`'s
    // `node + edge` return — i.e. Rust's `rank_score`. Using `score` here
    // seeds the GF threshold below Java's level by the per-PSM edge_score
    // value (~+20 typical), widening the score distribution and biasing
    // SpecEValue.
    let min_score = queue
        .iter_psms()
        .map(|p| p.rank_score.round() as i32)
        .min()
        .unwrap_or(i32::MIN);

    // parent_mass = (mz - PROTON) * charge  (precursor peak mass + proton, as in NewScoredSpectrum).
    let parent_mass = (spec.precursor_mz - PROTON) * (charge as f64);

    // 3. Derive protein-terminal flags by OR-ing across ALL PSMs in the queue.
    //
    // Aggregates `use_protein_n_term` / `use_protein_c_term` across all
    // candidates before GF construction. Iterates the full queue and sets
    // either flag the moment any PSM is at a protein N- or C-terminus,
    // short-circuiting once both are set.
    let (use_protein_n_term, use_protein_c_term) = {
        let mut any_n = false;
        let mut any_c = false;
        for psm in queue.iter_psms() {
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            if cand.is_protein_n_term { any_n = true; }
            if cand.is_protein_c_term { any_c = true; }
            if any_n && any_c { break; }
        }
        (any_n, any_c)
    };

    // 3b. Build the GF group across the nominal mass range.
    let mut group = GeneratingFunctionGroup::new();

    for nominal_mass_idx in min_peptide_mass_idx..=max_peptide_mass_idx {
        if nominal_mass_idx <= 0 {
            continue;
        }
        // Use the thread-local arena-pooled constructor: eliminates 11
        // Vec allocations per call (~4.4M allocs per PXD001819 run) by
        // recycling the buffers between graph builds. Output is bit-
        // identical to `new` (gated by primitive_graph_arena_parity tests).
        let graph = PrimitiveAaGraph::new_pooled(
            aa_set,
            nominal_mass_idx,
            enzyme,
            scored_spec,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
            use_protein_n_term,
            use_protein_c_term,
        );
        GF_BIN_ATTEMPTS.fetch_add(1, Ordering::Relaxed);
        match GeneratingFunction::with_score_threshold(&graph, min_score, aa_set) {
            Ok(gf) => group.accept(gf),
            Err(scoring_crate::gf::generating_function::GfError::EmptyScoreRange { .. }) => {
                GF_EMPTY_SCORE_RANGE.fetch_add(1, Ordering::Relaxed);
                continue;
            }
            Err(scoring_crate::gf::generating_function::GfError::SinkUnreachable) => {
                // SinkUnreachable from the thresholded DP means the
                // score-threshold pre-pass (`setup_score_threshold`) pruned
                // every path from source to sink because no AA-path could
                // theoretically reach the queue's `min_score`. This is a
                // pruning artifact, not a real reachability problem: the
                // unthresholded DP (`GeneratingFunction::compute`) still has
                // valid paths to compute a complete distribution from. Retry
                // without the threshold to recover ~10% of bin attempts that
                // would otherwise emit sentinel DeNovoScore / lnSpecEValue=0
                // and leave Percolator with broken features on ~5K Astral PSMs.
                GF_SINK_UNREACHABLE.fetch_add(1, Ordering::Relaxed);
                if let Ok(gf) = GeneratingFunction::compute(&graph, aa_set) {
                    GF_SINK_RETRY_OK.fetch_add(1, Ordering::Relaxed);
                    group.accept(gf);
                }
                continue;
            }
            Err(_) => continue,
        }
    }

    if !group.is_computed() {
        GF_SPECTRA_NO_GROUP.fetch_add(1, Ordering::Relaxed);
        return;
    }

    // 4. For each PSM in the queue, compute spec_e_value from its score.
    //
    // Use `rank_score` (Java-aligned `node + cleavage + edge`),
    // not `score` (Rust pin-only `node + cleavage`). Java parity:
    // `gf.getSpectralProbability(match.getScore())` where `match.getScore()`
    // is `node + cleavage + edge`. Using
    // `score` here would look up the wrong tail of the GF score distribution
    // (lower by the per-PSM edge contribution ~+20), giving inflated
    // SpecEValue values for PSMs whose top-1 was chosen via edge contribution.
    let max_score = group.max_score();

    queue.update_spec_e_values(|psm| {
        // Nominal peptide mass: residue masses sum + no water (mass-index convention).
        // Use nominal_from() (INTEGER_MASS_SCALER-aware) to match how graph nodes are indexed.
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let psm_nominal_mass = cand.peptide.nominal_residue_mass();
        if psm_nominal_mass < min_peptide_mass_idx || psm_nominal_mass > max_peptide_mass_idx {
            return 1.0;
        }
        let score_int = psm.rank_score.round() as i32;
        if score_int >= max_score {
            // Score exceeds GF range — return the probability at max_score - 1
            // (which already has the underflow guard applied by the GF DP).
            // Avoids returning a grossly inflated value (1/max_score ≈ 0.01)
            // that would invert ranking of the best PSMs.
            return group.spectral_probability(max_score - 1)
                .unwrap_or(f32::from_bits(1) as f64);
        }
        group.spectral_probability(score_int).unwrap_or(1.0)
    });

    // 5. Enrichment: set de_novo_score and e_value for output writers.
    //
    // de_novo_score = group.max_score() - 1.
    //
    // e_value = spec_e_value * num_distinct_peptides_at_length.
    //
    // Align lookup index with Java parity.
    //     `sa.getNumDistinctPeptides(enzyme == null ? length - 2 : length - 1)`
    // where `match.getLength() = pepLength + 2` (flanking residues included in
    // the stored length). So Java effectively queries
    //   - with enzyme: `numDistinctPeptides[pepLength + 1]`
    //   - without enzyme: `numDistinctPeptides[pepLength]`
    let de_novo_score = max_score - 1;
    let lookup_offset = match params.enzyme {
        Enzyme::NoCleavage | Enzyme::NonSpecific => 0,
        _ => 1,
    };
    queue.update_psm_enrichment(|psm| {
        psm.de_novo_score = de_novo_score;
        let len = candidates[psm.primary_candidate_idx() as usize].peptide.length();
        let num_distinct = search_index
            .num_distinct_peptides_at_length(len + lookup_offset)
            .max(1);
        psm.e_value = psm.spec_e_value * num_distinct as f64;
    });
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
    let tol_is_ppm = scorer.param().data_type.instrument.is_high_resolution();
    let tol = if tol_is_ppm { 20.0_f64 } else { 0.5_f64 };
    for p in &predicted {
        let tol_da = if tol_is_ppm { p.mz * tol / 1e6 } else { tol };
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
) -> PsmFeatures {
    let n = peptide.length();
    if n < 2 {
        return PsmFeatures::default();
    }

    // ADDITIVE Java-parity edge-score feature (new PIN column). Computed
    // here so it shares the per-PSM ScoredSpectrum + scorer references that
    // the existing feature-extraction code already has on hand.
    let edge_score = psm_edge_score(scored_spec, peptide, scorer, charge);

    // Predict charge-1 b/y ions; one bool per fragment position.
    //
    // Stack-allocate b/y_matched on a 64-slot SmallVec (max
    // peptide length is 40 → n-1 ≤ 39). The prior `vec![false; n-1]` heap
    // allocations fired ~150k × 4 / PSM batch and were a measurable hot-path
    // cost. SmallVec inlines for n ≤ 64.
    let predicted = predict_by_ions(peptide, 1..=1);
    let mut b_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut y_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];

    // Collect matched-ion details for ion-current ratio and error-stat features.
    // Each entry: (intensity, observed_mz, predicted_mz, is_b_ion).
    // SmallVec inlines for up to ~96 matched ions (b+y at n positions, with
    // some headroom for partition multi-ion-type matches at long peptides).
    let mut matched_ions: SmallVec<[(f32, f64, f64, bool); 96]> = SmallVec::new();

    // Java parity: feature-counting uses a
    // HARDCODED fragment tolerance, NOT param.mme. High-res instruments
    // (HighRes / TOF / QExactive) get 20 ppm; low-res LTQ gets 0.5 Da.
    // The param.mme value (0.5 Da for HCD_QExactive_Tryp.param) is the
    // coarser binning tolerance used by the rank-distribution tables —
    // appropriate for node-score lookup but ~50× too wide for feature
    // counting at m/z 500. Pre-fix Rust used param.mme for both, which
    // inflated NumMatchedMainIons by ~+3, longest_b by ~+2 vs Java, and
    // compressed all intensity ratios (more low-intensity noise matched
    // into the matched-ion sum). Confirmed by the Rust-vs-Java pin-diff
    // harness.
    let feature_tol = if scorer.param().data_type.instrument.is_high_resolution() {
        20.0_f64 // ppm
    } else {
        0.5_f64 // Da
    };
    let feature_tol_is_ppm = scorer.param().data_type.instrument.is_high_resolution();

    for p in &predicted {
        let tol_da = if feature_tol_is_ppm {
            p.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        if let Some((_rank, intensity, peak_mz)) =
            scored_spec.nearest_peak_full(p.mz, tol_da)
        {
            let is_b = matches!(p.kind, IonKind::B);
            matched_ions.push((intensity, peak_mz, p.mz, is_b));

            // position is 1-based (b1/y1 = index 0 in the matched arrays)
            let pos = (p.position - 1) as usize;
            match p.kind {
                IonKind::B => {
                    if pos < b_matched.len() {
                        b_matched[pos] = true;
                    }
                }
                IonKind::Y => {
                    if pos < y_matched.len() {
                        y_matched[pos] = true;
                    }
                }
            }
        }
    }

    // NumMatchedMainIons mirrors Java's PSMFeatureFinder count: each (bond, direction)
    // tuple contributes 1 if at least one charge-1 prefix/suffix ion matched.
    // Rust's b/y-charge-1 path above is a faithful subset of Java's
    // `getMassErrorWithIntensity`-driven count (which iterates the partition
    // ion list filtered to charge 1; for HCD_QExactive_Tryp the dominant
    // charge-1 prefix/suffix ions ARE b/y plus a few low-impact variants).
    let num_matched: u32 = (b_matched.iter().filter(|&&m| m).count()
        + y_matched.iter().filter(|&&m| m).count()) as u32;

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

    // ── Ion-current ratio features (partition-ion-list) ────────────────────────
    //
    // Java parity: `NewScoredSpectrum.getExplainedIonCurrent`
    // iterates the FULL partition ion list across all segments (b, y, plus
    // partition-specific variants like a-ion, b-H2O, etc.) and sums matched
    // peak intensities. The current Rust matched-ion buffer above only
    // contains b/y at charge 1, so it systematically UNDER-counts the
    // intensity sum. Rust-vs-Java pin-diff confirms: ExplainedIonCurrentRatio
    // median -0.026, NTerm -0.005, CTerm -0.018 — all compressed.
    //
    // We replace the b/y-only sum with a partition-wide sum AND use
    // partition-wide matches to drive longest_b/y (matches Java's "bIC > 0"
    // test). NumMatchedMainIons continues to count charge-1 b/y matches.
    let parent_mass = scored_spec.parent_mass();
    let num_segments = scorer.param().num_segments.max(1) as usize;

    // Stack-allocate (same rationale as b/y_matched above).
    let mut b_any_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut y_any_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut sum_prefix_intensity: f64 = 0.0;
    let mut sum_suffix_intensity: f64 = 0.0;

    // Use ACCURATE residue mass for theo m/z computation (matches Java's
    // PSMFeatureFinder which passes `peptide.get(i).getAccurateMass()`).
    // IonType::mz internally divides nominal mass by INTEGER_MASS_SCALER
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

        // Java iterates each segment's ion list separately and checks that
        // the computed theoMass falls into that segment (line 271-273). We
        // mirror that exactly so per-bond ion sums match Java's bIC / yIC.
        for seg in 0..num_segments {
            let ions = segment_ions[seg];
            for &ion in ions {
                let (is_prefix, residue_mass) = match ion {
                    scoring_crate::param_model::IonType::Prefix { charge: ic, offset_bits } => {
                        let offset = f32::from_bits(offset_bits) as f64;
                        let z = ic as f64;
                        (true, (prm_accurate / z + offset, ion))
                    }
                    scoring_crate::param_model::IonType::Suffix { charge: ic, offset_bits } => {
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

    // All four *ErrorTop7 columns are in PPM (matching Java
    // `NewScoredSpectrum.getMassErrorWithIntensity`, which always returns
    // `(p.getMz() - theoMass) / theoMass * 1e6f`). The Java column naming
    // is misleading: `MeanErrorTop7` = mean of |ppm error| (absolute),
    // `MeanRelErrorTop7` = mean of signed ppm error. Both are ppm; the
    // "Rel" suffix in Java distinguishes signed vs absolute, NOT
    // Da-vs-ppm. MeanErrorTop7/StdevErrorTop7 are therefore emitted as
    // absolute ppm (not Da) to match Java's units.
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
    }
}

// ── Unit tests for feature columns ───────────────────────────────────────────

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
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: (PROTON as f32).to_bits() };
        let suffix1 = IonType::Suffix { charge: 1, offset_bits: ((H2O + PROTON) as f32).to_bits() };
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
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.5), 2);
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
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.01), 2); // tight tolerance

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

    // ── Test: top-7 error stats are nonzero when ions match ─────────────────

    #[test]
    fn compute_psm_features_error_stats_nonzero_when_ions_match_with_offset() {
        // Build a peptide and shift every peak by a fixed offset so errors are known.
        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);

        // 0.0005 Da offset = ~6 ppm at m/z 89 (Ala b1) — within the
        // hardcoded 20 ppm window that compute_psm_features now uses for
        // high-resolution instruments (Java parity).
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

/// Pre-merge dedup pass: collapse PSMs sharing the same Java
/// `(pepSeq, score)` key before per-charge GF compute.
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

/// Mod-aware dedup key: bare residues plus per-position mod mass (1e-5 Da units).
/// Matches Java pepSeq semantics without string formatting on the hot path.
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
            spec_e_value: 1.0,
            de_novo_score: 0,
            activation_method: None,
            e_value: 1.0,
            features: Default::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
        }
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
