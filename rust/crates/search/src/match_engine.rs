//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use std::collections::{BTreeMap, HashMap};
use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

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
use scoring_crate::scoring::{score_psm, RankScorer, ScoredSpectrum};
use model::spectrum::Spectrum;

/// Match every spectrum against every candidate from the SearchIndex.
/// Returns one top-N PSM queue per spectrum, in input order.
///
/// A `ScoredSpectrum` is built once per spectrum and reused across all
/// candidates; candidates are bucketed by mass for sub-linear precursor
/// lookup. After per-candidate scoring, SpecEValue is computed via the
/// generating-function DP across the precursor tolerance window in nominal
/// mass space and assigned to every PSM in the queue.
pub fn match_spectra(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    scorer: &RankScorer,
    fragment_tolerance_da: f64,
    decoy_prefix: &str,
) -> Vec<TopNQueue> {
    // Populate the per-length distinct-peptide counts on the SearchIndex.
    // Idempotent + lock-free (OnceLock); tests that pre-populate via
    // `with_distinct_peptide_counts` keep their populated map. Required so
    // `idx.num_distinct_peptides_at_length(...)` returns real values during
    // production search (the future Phase 7 e_value swap consumes it).
    idx.ensure_distinct_peptide_counts(params, decoy_prefix);

    let candidates: Vec<Candidate> = enumerate_candidates(idx, params, decoy_prefix).collect();

    // Build mass-bucket index: nominal(peptide.mass() - H2O) → Vec<candidate_idx>.
    //
    // Uses the same nominal_from convention as the GF mass-bin loop so that
    // bucket keys align with the GF's mass-bin lookup (commit b89779a fix).
    // Stores only indices into `candidates` — no cloning, tiny memory overhead.
    let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (cand_idx, cand) in candidates.iter().enumerate() {
        let nominal = nominal_from(cand.peptide.mass() - H2O);
        bucket_index.entry(nominal).or_default().push(cand_idx);
    }

    // Build an aa_set clone with enzyme registered (for GF cleavage scoring).
    // We use Java MS-GF+ defaults: peptide_eff = 0.95, neighboring_eff = 0.95.
    // Cloning is cheap (AminoAcidSet is a HashMap of ~20 entries).
    // This avoids mutating the shared SearchParams.aa_set borrow.
    let mut aa_set_for_gf: AminoAcidSet = params.aa_set.clone();
    if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
        aa_set_for_gf.register_enzyme(params.enzyme, 0.95, 0.95);
    }

    // Yield-accounting counters.
    // Aggregated across all worker threads via Relaxed atomics — exact counts
    // don't require ordering with other memory ops.
    let skipped_min_peaks = AtomicU64::new(0);
    let candidates_visited = AtomicU64::new(0);
    let psms_pushed = AtomicU64::new(0);
    let spectra_with_psms = AtomicU64::new(0);

    // Parallel per-spectrum search.
    //
    // Mirrors Java DBScanner.computeSpecEValues which runs per-spectrum work on
    // `-thread N` workers. All inputs above are `&` immutable; the closure owns
    // its TopNQueue, scored_per_charge cache, and per-bin GF state.
    let queues: Vec<TopNQueue> = spectra
        .par_iter()
        .enumerate()
        .map(|(spec_idx, spec)| {
            let mut queue = TopNQueue::new(params.top_n_psms_per_spectrum);

            // Skip spectra with too few peaks (mirrors Java's `-minNumPeaks` filter).
            if spec.peaks.len() < params.min_peaks as usize {
                skipped_min_peaks.fetch_add(1, Ordering::Relaxed);
                return queue;
            }

            // Determine which charge states to try for this spectrum.
            // For charge-explicit spectra this is a single entry; for charge-missing,
            // typically 2-3 entries (small overhead, correct behavior).
            let charges_to_try: Vec<u8> = match spec.precursor_charge {
                Some(z) if z > 0 => vec![z as u8],
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
            let mut scored_per_charge: HashMap<u8, ScoredSpectrum<'_>> = HashMap::new();
            for &z in &charges_to_try {
                scored_per_charge.entry(z)
                    .or_insert_with(|| ScoredSpectrum::new(spec, scorer, z));
            }

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
                // Java convention: max widens by tol_da_left, min widens by tol_da_right.
                let min_nominal = nominal_center - iso_max - widen_right;
                let max_nominal = nominal_center - iso_min + widen_left;
                for (_nm, idxs) in bucket_index.range(min_nominal..=max_nominal) {
                    window_cand_indices.extend_from_slice(idxs);
                }
            }
            window_cand_indices.sort_unstable();
            window_cand_indices.dedup();

            // Per-candidate cleavage credit. Mirrors Java DBScanner.java:441
            //   `cleavageScore = nTermCleavageScore + cTermCleavageScore`
            // added to `scorer.getScore(...)` before queue insertion.
            // For C-term enzymes (Trypsin, default):
            //   - N-term: credit if isProteinNTerm OR enzyme.isCleavable(prevAA),
            //     else penalty
            //   - C-term: credit if enzyme.isCleavable(lastResidue), else penalty
            // Omitting this offsets every PSM score by ≈ -(credit + credit) = -4
            // vs Java for fully tryptic ntt=2 candidates.
            //
            // Use the ENZYME-REGISTERED aa_set (cleavage credit/penalty are
            // populated by register_enzyme — params.aa_set is unregistered).
            let compute_cleavage_credit = |cand: &Candidate| -> i32 {
                let aa_set = &aa_set_for_gf;
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

            for &cand_idx in &window_cand_indices {
                let cand = &candidates[cand_idx];
                let cleavage_credit = compute_cleavage_credit(cand) as f32;
                for &z in &charges_to_try {
                    let scored_spec = &scored_per_charge[&z];
                    let mut best_for_charge: Option<(MassError, f32)> = None;
                    for offset in params.isotope_error_range.clone() {
                        if let Some(err) = matches_precursor(spec, &cand.peptide, z, offset, &params.precursor_tolerance) {
                            // Add cleavage credit (Java DBScanner.java:513:
                            // `score = cleavageScore + scorer.getScore(...)`).
                            let score = score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da)
                                + cleavage_credit;
                            if best_for_charge.as_ref().map_or(true, |(_, s)| score > *s) {
                                best_for_charge = Some((err, score));
                            }
                        }
                    }
                    if let Some((err, score)) = best_for_charge {
                        let features = compute_psm_features(scored_spec, &cand.peptide, scorer);
                        queue.push(PsmMatch {
                            spectrum_idx: spec_idx,
                            candidate: cand.clone(),
                            charge_used: z,
                            mass_error_ppm: err.mass_error_ppm,
                            score,
                            spec_e_value: 1.0,  // set by compute_spec_e_values_for_spectrum
                            de_novo_score: i32::MIN,  // set by compute_spec_e_values_for_spectrum
                            activation_method: Some(scorer.param().data_type.activation),
                            e_value: 1.0,  // set by compute_spec_e_values_for_spectrum
                            features,
                            isotope_offset: err.isotope_offset,
                        });
                        psms_pushed.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
            candidates_visited.fetch_add(window_cand_indices.len() as u64, Ordering::Relaxed);

            // Compute SpecEValue for the PSMs in this queue.
            if !queue.is_empty() {
                spectra_with_psms.fetch_add(1, Ordering::Relaxed);
                let enzyme_opt = if params.enzyme != Enzyme::NoCleavage
                    && params.enzyme != Enzyme::NonSpecific
                {
                    Some(params.enzyme)
                } else {
                    None
                };
                // Pick the ScoredSpectrum for the top PSM's charge.
                let top_charge = queue
                    .iter_psms()
                    .max_by(|a, b| a.cmp(b))
                    .map(|p| p.charge_used)
                    .unwrap_or(charges_to_try[0]);
                let scored_spec_for_gf = &scored_per_charge[&top_charge];
                compute_spec_e_values_for_spectrum(
                    spec,
                    params,
                    &mut queue,
                    &aa_set_for_gf,
                    enzyme_opt,
                    scorer,
                    scored_spec_for_gf,
                    top_charge,
                    fragment_tolerance_da,
                    idx,
                );
            }

            queue
        })
        .collect();

    // Yield-accounting summary.
    // Helps disambiguate whether a PSM-yield gap vs Java is from:
    //   - filtering (skipped_min_peaks)
    //   - enumeration (candidates_visited)
    //   - scoring (psms_pushed)
    //   - top-N retention (spectra_with_psms)
    eprintln!(
        "Yield: {} spectra in, {} skipped by min_peaks, {} candidates visited, \
         {} PSMs pushed, {} spectra with non-empty queue",
        spectra.len(),
        skipped_min_peaks.load(Ordering::Relaxed),
        candidates_visited.load(Ordering::Relaxed),
        psms_pushed.load(Ordering::Relaxed),
        spectra_with_psms.load(Ordering::Relaxed),
    );

    queues
}

/// For a single spectrum, compute the GF across the precursor tolerance
/// window in nominal mass space, then assign `spec_e_value` to every PSM
/// in `queue` whose nominal_peptide_mass falls within the window.
///
/// Mirrors Java DBScanner.java:597-650.
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

    // Java isotope error convention: range [min_iso, max_iso] is applied as
    //   minNominalPeptideMass = nominalPeptideMass - maxIsotopeError
    //   maxNominalPeptideMass = nominalPeptideMass - minIsotopeError
    let iso_min = *params.isotope_error_range.start() as i32;
    let iso_max = *params.isotope_error_range.end() as i32;
    let min_iso_nominal = nominal_peptide_mass - iso_max;
    let max_iso_nominal = nominal_peptide_mass - iso_min;

    // Tolerance widening: Java uses Math.round(tol_da - 0.4999f).
    // tolDaLeft governs the upper bound; tolDaRight governs the lower bound.
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

    // parent_mass = (mz - H) * charge  (precursor peak mass, with H added back in Java).
    let parent_mass = (spec.precursor_mz - PROTON) * (charge as f64);

    // 3. Derive protein-terminal flags by OR-ing across ALL PSMs in the queue.
    //
    // Java reference: DBScanner.java:592-602 aggregates useProteinNTerm /
    // useProteinCTerm across all candidates before GF construction. We mirror
    // this by iterating the full queue and setting either flag the moment any
    // PSM is at a protein N- or C-terminus, short-circuiting once both are set.
    let (use_protein_n_term, use_protein_c_term) = {
        let mut any_n = false;
        let mut any_c = false;
        for psm in queue.iter_psms() {
            if let Some(prot) = search_index.protein_at(psm.candidate.protein_index) {
                let start = psm.candidate.start_offset_in_protein;
                let pep_len = psm.candidate.peptide.length();
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
        match GeneratingFunction::with_score_threshold(&graph, min_score, aa_set) {
            Ok(gf) => group.accept(gf),
            Err(_) => continue, // skip degenerate / unreachable bins
        }
    }

    if !group.is_computed() {
        return;
    }

    // 4. For each PSM in the queue, compute spec_e_value from its score.
    let max_score = group.max_score();

    queue.update_spec_e_values(|psm| {
        // Nominal peptide mass: residue masses sum + no water (Java convention for mass index).
        // Use nominal_from() (INTEGER_MASS_SCALER-aware) to match how graph nodes are indexed.
        let psm_nominal_mass = nominal_from(psm.candidate.peptide.mass() - H2O);
        if psm_nominal_mass < min_peptide_mass_idx || psm_nominal_mass > max_peptide_mass_idx {
            return 1.0;
        }
        let score_int = psm.score.round() as i32;
        if score_int >= max_score {
            // Score exceeds GF range — return the probability at max_score - 1
            // (which already has the underflow guard applied by the GF DP).
            // Mirrors Java behavior; avoids returning a grossly inflated value
            // (1/max_score ≈ 0.01) that would invert ranking of the best PSMs.
            return group.spectral_probability(max_score - 1)
                .unwrap_or(f32::from_bits(1) as f64);
        }
        group.spectral_probability(score_int).unwrap_or(1.0)
    });

    // 5. Enrichment: set de_novo_score and e_value for output writers.
    //
    // de_novo_score = group.max_score() - 1  (mirrors Java's getDeNovoScore()).
    //
    // e_value = spec_e_value * num_distinct_peptides_at_length.
    // num_distinct sourced from SearchIndex (mirrors Java
    // PeptideEnumerator.getNumDistinctPeptides), replacing the prior
    // top-N-queue-derived proxy.
    let de_novo_score = max_score - 1;
    queue.update_psm_enrichment(|psm| {
        psm.de_novo_score = de_novo_score;
        let len = psm.candidate.peptide.length();
        let num_distinct = search_index.num_distinct_peptides_at_length(len).max(1);
        psm.e_value = psm.spec_e_value * num_distinct as f64;
    });
}

/// Compute fragment-ion feature columns for a single PSM.
///
/// Uses charge-1 b/y ions only (matching Java's `NumMatchedMainIons`
/// convention).  A peptide position counts at most once per ion series;
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
///   Mirrors Java `PSMFeatureFinder.computeExplainedIonCurrent()`.
/// - `MS2IonCurrent` is the raw sum (NOT log10).  Java `getMS2IonCurrent()`
///   returns the raw sum; the PIN emitter emits it as-is.
/// - `IsolationWindowEfficiency` is always 0.0; Java returns `null` here
///   (no isolation-window data in the Spectrum object).
/// - Top-7 error stats mirror Java `MassErrorStat`: errors are collected for
///   all matched b+y ions, sorted descending by intensity, top-7 taken;
///   absolute Da error for mean/stdev, signed ppm for rel-mean/rel-stdev.
///   Population stdev formula: `sqrt(E[x²] - mean²)` — matches Java.
pub(crate) fn compute_psm_features(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
) -> PsmFeatures {
    let n = peptide.length();
    if n < 2 {
        return PsmFeatures::default();
    }

    // Predict charge-1 b/y ions; one bool per fragment position.
    let predicted = predict_by_ions(peptide, 1..=1);
    let mut b_matched = vec![false; n - 1];
    let mut y_matched = vec![false; n - 1];

    // Collect matched-ion details for ion-current ratio and error-stat features.
    // Each entry: (intensity, observed_mz, predicted_mz, is_b_ion)
    let mut matched_ions: Vec<(f32, f64, f64, bool)> = Vec::new();

    for p in &predicted {
        let tol_da = scorer.param().mme.as_da(p.mz);
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
    // Java `getMS2IonCurrent()` returns the raw sum (no log10 transform).
    let ms2_ion_current = if total_intensity > 0.0 { total_intensity as f32 } else { 0.0 };
    // Java `getIsolationWindowEfficiency()` always returns null → emit 0.0.
    let isolation_window_efficiency = 0.0_f32;

    // ── Top-7 mass-error statistics ───────────────────────────────────────────

    // Sort matched ions descending by intensity (mirrors Java MassErrorStat
    // which sorts errorList by intensity via PairReverseComparator).
    let mut by_intensity = matched_ions.clone();
    by_intensity.sort_by(|a, b| {
        b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
    });
    let top7: Vec<(f32, f64, f64, bool)> = by_intensity.into_iter().take(7).collect();

    // Java MassErrorStat: absolute Da errors for mean7/sd7;
    //                     signed errors (no abs) for rMean7/rSd7 (ppm).
    // Population stdev formula (Java): sqrt(sumSq/n - mean²).
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
        }
    }

    // ── Test: empty spectrum → all new features are 0 ───────────────────────

    #[test]
    fn compute_psm_features_top7_error_stats_zero_when_no_matches() {
        let pep = ala_peptide(4);
        let spec = make_spectrum(vec![]); // no peaks
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.5));
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
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.01)); // tight tolerance

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

        let offset_da = 0.01_f64;  // 10 mDa error on every peak
        let mut peaks: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz + offset_da, (i + 1) as f32 * 10.0))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = make_spectrum(peaks);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        // tolerance = 0.05 Da so all offset peaks are still within window.
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.05));

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
        let f = compute_psm_features(&ss, &pep, &make_scorer(0.5));

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
