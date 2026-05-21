//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use std::collections::{BTreeMap, HashMap};
use std::hash::Hasher;
use std::sync::atomic::{AtomicU64, Ordering};

// GF failure-mode diagnostics (2026-05-19). Module-level atomics
// incremented per-bin from compute_spec_e_values_for_spectrum and
// reported in the yield-accounting summary. Used to characterise the
// ~4.7% of Astral PSMs where GF compute fails (docs/parity-analysis/
// notes/2026-05-19-gf-compute-failures.md). Module-level rather than
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
use rustc_hash::{FxHashSet, FxHasher};
use smallvec::{smallvec, SmallVec};

use model::aa_set::AminoAcidSet;
use crate::candidate_gen::{enumerate_candidates, Candidate};
use model::enzyme::Enzyme;
use scoring_crate::gf::generating_function::GeneratingFunction;
use scoring_crate::gf::group::GeneratingFunctionGroup;
use scoring_crate::gf::primitive_graph::PrimitiveAaGraph;
use model::mass::{nominal_from, H2O, PROTON};
use model::peptide::Peptide;
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
        // bucket keys align with the GF's mass-bin lookup (commit b89779a fix).
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
        }
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
        let params = self.params;
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

        // Parallel per-spectrum search. All inputs above are `&` immutable; the
        // closure owns its TopNQueue, scored_per_charge cache, and per-bin GF state.
        let queues: Vec<TopNQueue> = spectra
            .par_iter()
            .enumerate()
            .map(|(local_idx, spec)| {
                let spec_idx = local_idx + spectrum_idx_offset;
                let mut queue = TopNQueue::new(params.top_n_psms_per_spectrum);

            // Skip spectra with too few peaks.
            if spec.peaks.len() < params.min_peaks as usize {
                skipped_min_peaks.fetch_add(1, Ordering::Relaxed);
                return queue;
            }

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
            for &z in &charges_to_try {
                let charge_f = z as f64;
                let neutral_mass = (spec.precursor_mz - PROTON) * charge_f - H2O;
                let nominal_center = nominal_from(neutral_mass);
                let iso_min = *params.isotope_error_range.start() as i32;
                let iso_max = *params.isotope_error_range.end() as i32;
                let tol_da_left  = params.precursor_tolerance.left.as_da(neutral_mass);
                let tol_da_right = params.precursor_tolerance.right.as_da(neutral_mass);
                let widen_left  = (tol_da_left  - 0.4999_f64).round() as i32;
                let widen_right = (tol_da_right - 0.4999_f64).round() as i32;
                // Convention: max widens by tol_da_left, min widens by tol_da_right.
                let min_nominal = nominal_center - iso_max - widen_right;
                let max_nominal = nominal_center - iso_min + widen_left;
                for (_nm, idxs) in bucket_index.range(min_nominal..=max_nominal) {
                    window_cand_indices.extend_from_slice(idxs);
                }
            }
            window_cand_indices.sort_unstable();
            window_cand_indices.dedup();

            // Per-candidate cleavage credit:
            //   `cleavage_score = n_term_cleavage_score + c_term_cleavage_score`
            // added to the raw PSM score before queue insertion.
            // For C-term enzymes (Trypsin, default):
            //   - N-term: credit if isProteinNTerm OR enzyme.is_cleavable(prevAA),
            //     else penalty
            //   - C-term: credit if enzyme.is_cleavable(lastResidue), else penalty
            // Omitting this offsets every PSM score by ≈ -(credit + credit) = -4
            // for fully tryptic ntt=2 candidates.
            //
            // Use the ENZYME-REGISTERED aa_set (cleavage credit/penalty are
            // populated by register_enzyme — params.aa_set is unregistered).
            let compute_cleavage_credit = |cand: &Candidate| -> i32 {
                let aa_set = aa_set_for_gf;
                let enz = params.enzyme;
                let mut score: i32 = 0;
                let pre = cand.peptide.pre;
                let last = cand.peptide.residues.last().map(|aa| aa.residue).unwrap_or(0);
                let post = cand.peptide.post;
                if enz.is_c_term() {
                    // N-term cleavage
                    score += if cand.is_protein_n_term || enz.is_cleavable(pre) {
                        aa_set.neighboring_aa_cleavage_credit()
                    } else {
                        aa_set.neighboring_aa_cleavage_penalty()
                    };
                    // C-term cleavage
                    score += if enz.is_cleavable(last) {
                        aa_set.peptide_cleavage_credit()
                    } else {
                        aa_set.peptide_cleavage_penalty()
                    };
                } else if enz.is_n_term() {
                    // N-term cleavage (peptide N-term)
                    score += if enz.is_cleavable(pre) {
                        aa_set.peptide_cleavage_credit()
                    } else {
                        aa_set.peptide_cleavage_penalty()
                    };
                    // C-term cleavage (neighbor)
                    score += if cand.is_protein_c_term || enz.is_cleavable(post) {
                        aa_set.neighboring_aa_cleavage_credit()
                    } else {
                        aa_set.neighboring_aa_cleavage_penalty()
                    };
                }
                score
            };

            // R-2.1: per-charge queue keyed by charge state. Mirrors Java's
            // per-SpecKey raw-score retention (DBScanner.java:534).
            let mut per_charge_queues: HashMap<u8, TopNQueue> = HashMap::new();

            for &cand_idx in &window_cand_indices {
                let cand = &candidates[cand_idx];
                let cleavage_credit = compute_cleavage_credit(cand) as f32;
                for &z in &charges_to_try {
                    let scored_spec = scored_spec_for_charge(z);
                    let mut best_for_charge: Option<(MassError, f32)> = None;
                    for offset in params.isotope_error_range.clone() {
                        if let Some(err) = matches_precursor(spec, &cand.peptide, z, offset, &params.precursor_tolerance) {
                            let score = score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da)
                                + cleavage_credit;
                            if best_for_charge.as_ref().map_or(true, |(_, s)| score > *s) {
                                best_for_charge = Some((err, score));
                            }
                        }
                    }
                    if let Some((err, score)) = best_for_charge {
                        let features = PsmFeatures::default();
                        let psm = PsmMatch {
                            spectrum_idx: spec_idx,
                            candidate_idxs: vec![cand_idx as u32],
                            charge_used: z,
                            mass_error_ppm: err.mass_error_ppm,
                            score,
                            spec_e_value: 1.0,
                            de_novo_score: i32::MIN,
                            activation_method: Some(scorer.param().data_type.activation),
                            e_value: 1.0,
                            features,
                            isotope_offset: err.isotope_offset,
                        };
                        per_charge_queues
                            .entry(z)
                            .or_insert_with(|| TopNQueue::new(params.top_n_psms_per_spectrum))
                            .push(psm);
                        psms_pushed.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
            candidates_visited.fetch_add(window_cand_indices.len() as u64, Ordering::Relaxed);

            // R-2.2: pepSeq + score dedup per-charge BEFORE GF compute.
            // Same peptide matched against multiple proteins collapses to one
            // PsmMatch with aggregated candidate_idxs (Java DBScanner.java:719-733).
            for queue in per_charge_queues.values_mut() {
                if queue.len() > 1 {
                    let drained = queue.drain_into_vec();
                    let deduped = dedup_pepseq_score(drained, candidates);
                    for psm in deduped {
                        queue.push(psm);
                    }
                }
            }

            // R-2.3: per-charge GF / SpecEValue compute. Each per-charge queue
            // gets SpecE calibrated against its OWN charge's GF distribution
            // (Java DBScanner.java:606,779 — getRankScorer per SpecKey).
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

            // R-2.4: spectrum-level merge with SpecE tie keep. R-1's
            // TopNQueue::push (Ordering::Equal arm) keeps SpecE ties at
            // capacity because PsmMatch::cmp orders by spec_e_value first.
            // Matches Java DBScanner.java:745.
            for (_charge, mut per_charge) in per_charge_queues.drain() {
                for psm in per_charge.drain_into_vec() {
                    queue.push(psm);
                }
            }

            // Feature extraction (unchanged from baseline): post-merge, after
            // the per-spectrum queue is final.
            queue.fill_post_topn(|psm| {
                let ss = scored_spec_for_charge(psm.charge_used);
                let cand = &candidates[psm.primary_candidate_idx() as usize];
                psm.features = compute_psm_features(ss, &cand.peptide, scorer, psm.charge_used);
            });

                queue
            })
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
        // GF DP failure-mode diagnostics (2026-05-19; see
        // docs/parity-analysis/notes/2026-05-19-gf-compute-failures.md).
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
fn compute_spec_e_values_for_spectrum(
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
    let peptide_neutral_mass = (spec.precursor_mz - PROTON) * (charge as f64) - H2O;
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

    // 2. Compute the minimum score across all PSMs (used as score threshold).
    let min_score = queue
        .iter_psms()
        .map(|p| p.score.round() as i32)
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
            if let Some(prot) = search_index.protein_at(cand.protein_index) {
                let start = cand.start_offset_in_protein;
                let pep_len = cand.peptide.length();
                if start == 0 { any_n = true; }
                if start + pep_len >= prot.sequence.len() { any_c = true; }
                if any_n && any_c { break; }
            }
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
                // 2026-05-20: SinkUnreachable from the thresholded DP means the
                // score-threshold pre-pass (`setup_score_threshold`) pruned
                // every path from source to sink because no AA-path could
                // theoretically reach the queue's `min_score`. This is a
                // pruning artifact, not a real reachability problem: the
                // unthresholded DP (`GeneratingFunction::compute`) still has
                // valid paths to compute a complete distribution from. Retry
                // without the threshold to recover ~10% of bin attempts that
                // would otherwise emit sentinel DeNovoScore / lnSpecEValue=0
                // and leave Percolator with broken features on ~5K Astral PSMs.
                // See docs/parity-analysis/notes/2026-05-19-gf-compute-failures.md.
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
    let max_score = group.max_score();

    queue.update_spec_e_values(|psm| {
        // Nominal peptide mass: residue masses sum + no water (mass-index convention).
        // Use nominal_from() (INTEGER_MASS_SCALER-aware) to match how graph nodes are indexed.
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let psm_nominal_mass = cand.peptide.nominal_residue_mass();
        if psm_nominal_mass < min_peptide_mass_idx || psm_nominal_mass > max_peptide_mass_idx {
            return 1.0;
        }
        let score_int = psm.score.round() as i32;
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
    // HIGH-2 (2026-05-18): align lookup index with Java. Java's
    // `DirectPinWriter.java:165` does
    //     `sa.getNumDistinctPeptides(enzyme == null ? length - 2 : length - 1)`
    // where `match.getLength() = pepLength + 2` (DBScanner.java:521 includes the
    // two flanking residues in the stored length). So Java effectively queries
    //   - with enzyme: `numDistinctPeptides[pepLength + 1]`
    //   - without enzyme: `numDistinctPeptides[pepLength]`
    //
    // Rust previously queried `num_distinct(pepLength)` for both cases, which
    // was the right semantics for the "without enzyme" branch and an
    // off-by-one for the typical tryptic case.
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
/// All 9 previously zero-stubbed PIN columns are now filled:
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
    let predicted = predict_by_ions(peptide, 1..=1);
    let mut b_matched = vec![false; n - 1];
    let mut y_matched = vec![false; n - 1];

    // Collect matched-ion details for ion-current ratio and error-stat features.
    // Each entry: (intensity, observed_mz, predicted_mz, is_b_ion)
    let mut matched_ions: Vec<(f32, f64, f64, bool)> = Vec::new();

    // Java parity (PSMFeatureFinder.java:51-54): feature-counting uses a
    // HARDCODED fragment tolerance, NOT param.mme. High-res instruments
    // (HighRes / TOF / QExactive) get 20 ppm; low-res LTQ gets 0.5 Da.
    // The param.mme value (0.5 Da for HCD_QExactive_Tryp.param) is the
    // coarser binning tolerance used by the rank-distribution tables —
    // appropriate for node-score lookup but ~50× too wide for feature
    // counting at m/z 500. Pre-fix Rust used param.mme for both, which
    // inflated NumMatchedMainIons by ~+3, longest_b by ~+2 vs Java, and
    // compressed all intensity ratios (more low-intensity noise matched
    // into the matched-ion sum). Confirmed by iter16-vs-Java pin-diff
    // harness (docs/parity-analysis/notes/2026-05-19-pin-diff-findings.md).
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

    let longest_b = longest_run(&b_matched);
    let longest_y = longest_run(&y_matched);

    // ── Ion-current ratio features ────────────────────────────────────────────

    let total_intensity = scored_spec.total_intensity(); // raw sum, all peaks

    let matched_b_intensity: f64 = matched_ions.iter()
        .filter(|&&(_, _, _, is_b)| is_b)
        .map(|&(int, _, _, _)| int as f64)
        .sum();
    let matched_y_intensity: f64 = matched_ions.iter()
        .filter(|&&(_, _, _, is_b)| !is_b)
        .map(|&(int, _, _, _)| int as f64)
        .sum();
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

    // Absolute Da errors for mean7/sd7;
    // signed errors (no abs) for r_mean7/r_sd7 (ppm).
    // Population stdev formula: sqrt(sum_sq/n - mean²).
    let abs_da_errors: Vec<f64> = top7.iter()
        .map(|&(_, obs, pred, _)| (obs - pred).abs())
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

    let (mean_error_top7, stdev_error_top7)         = mean_and_pop_stdev(&abs_da_errors);
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
    use std::collections::HashMap;

    /// Minimal RankScorer for feature tests, with mme = Da(tol_da).
    fn make_scorer(tol_da: f64) -> RankScorer {
        let part = Partition { charge: 2, parent_mass: 0.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let noise = IonType::Noise;
        let mut ion_table = HashMap::new();
        ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
        ion_table.insert(noise, vec![0.1_f32, 0.2, 0.3, 0.4]);
        let mut rank_dist_table = HashMap::new();
        rank_dist_table.insert(part, ion_table);
        let mut frag_off_table = HashMap::new();
        frag_off_table.insert(part, vec![FragmentOffsetFrequency { ion_type: prefix1, frequency: 0.7 }]);
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
            precursor_off_map: HashMap::new(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
            partition_ion_types_cache: HashMap::new(),
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
        // high-resolution instruments (Java parity, PSMFeatureFinder.java:51-54).
        // The previous 0.01 Da offset assumed Rust used param.mme (~0.05 Da
        // in this fixture's make_scorer), but the iter20 fix makes feature
        // counting use 20 ppm regardless of param.mme.
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

        // All absolute Da errors should be ~offset_da.
        assert!(
            f.mean_error_top7 > 0.0,
            "mean_error_top7 should be > 0 when peaks are systematically offset, got {}",
            f.mean_error_top7
        );
        // With identical errors, stdev should be near 0.
        assert!(
            f.stdev_error_top7 < 1e-4,
            "stdev_error_top7 should be ~0 for identical errors, got {}",
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

/// Pre-merge dedup pass (R-2.2): collapse PSMs that share the same
/// (peptide_residue, rounded_score) key into a single entry, aggregating
/// their `candidate_idxs` into a unified Vec. Mirrors Java's
/// `DBScanner.java:719-733` `pepSeqMap` dedup.
///
/// Called by the per-spectrum loop after the per-candidate scoring loop,
/// before per-charge GF compute (so SpecE is computed on the deduped set).
///
/// Inputs:
/// - `psms`: drained from a per-charge `TopNQueue` via `drain_into_vec`
/// - `candidates`: the search's enumerated candidate slice; used to resolve
///   each PSM's peptide residue sequence for the dedup key
///
/// Returns: deduped `Vec<PsmMatch>`. The caller re-pushes these into the
/// per-charge queue via `queue.push()` for each entry.
///
/// Not yet called (Task 3); suppressing dead_code warning until integration.
#[allow(dead_code)]
pub(crate) fn dedup_pepseq_score(
    psms: Vec<PsmMatch>,
    candidates: &[Candidate],
) -> Vec<PsmMatch> {
    use std::collections::HashMap;

    // Key: (peptide_residue_bytes, rounded_score_i32)
    // The residue sequence is the unmodified bare AA string, matching Java's
    // `m.getPepSeq()` used as the dedup key (DBScanner.java:721).
    let mut groups: HashMap<(Vec<u8>, i32), PsmMatch> = HashMap::new();

    for psm in psms {
        let cand = &candidates[psm.primary_candidate_idx() as usize];
        let pep_residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
        let score_rounded = psm.score.round() as i32;
        let key = (pep_residues, score_rounded);

        groups
            .entry(key)
            .and_modify(|existing| {
                // Aggregate this PSM's indices into the surviving entry.
                // Avoid duplicates if the same idx somehow appears twice.
                for &idx in &psm.candidate_idxs {
                    if !existing.candidate_idxs.contains(&idx) {
                        existing.candidate_idxs.push(idx);
                    }
                }
            })
            .or_insert(psm);
    }

    groups.into_values().collect()
}
