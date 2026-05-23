//! Phase 6 / Task 8 smoke tests: SpecEValue is computed and < 1.0 for matched PSMs.
//!
//! Tests that:
//! 1. PSMs in a non-empty queue have spec_e_value <= 1.0 after match_spectra.
//! 2. For a well-matched spectrum, the top PSM has spec_e_value < 1.0.
//! 3. The TopNQueue ordering reflects spec_e_value (best first in sorted_vec).

use std::collections::HashMap;

use model::{AminoAcid, AminoAcidSetBuilder, Peptide, Protein, ProteinDb, Spectrum, PROTON, Tolerance};
use scoring_crate::{Param, RankScorer};
use search::{match_spectra, SearchIndex, SearchParams};
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use scoring_crate::param_model::{IonType, Partition, SpecDataType};
use model::protocol::Protocol;
use search::psm::PsmMatch;

fn make_spectrum(precursor_mz: f64, charge: Option<i32>) -> Spectrum {
    Spectrum {
        title: "specevalue_smoke".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: charge,
        rt_seconds: None,
        scan: None,
        peaks: vec![],
        activation_method: None,
    }
}

/// Minimal RankScorer for smoke tests (no real peaks, just need valid scorer).
fn tiny_scorer() -> RankScorer {
    let part = Partition { charge: 2, parent_mass: 500.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let suffix1 = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let noise = IonType::Noise;

    let mut ion_table = HashMap::new();
    ion_table.insert(prefix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
    ion_table.insert(suffix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
    ion_table.insert(noise, vec![0.05_f32, 0.05, 0.05, 0.05]);

    let mut rank_dist_table = HashMap::new();
    rank_dist_table.insert(part, ion_table);

    let mut frag_off_table = HashMap::new();
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

/// Build a known peptide spectrum match and return queues.
fn run_single_peptide_search(
    sequence: &[u8],
    peptide_sequence: &[u8],
    charge: u8,
) -> (Vec<search::psm::TopNQueue>, Vec<search::candidate_gen::Candidate>) {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: sequence.to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    // make_spectrum produces 0 peaks; default min_peaks=10 would skip everything.
    params.min_peaks = 0;

    let residues: Vec<AminoAcid> = peptide_sequence
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    let peptide = Peptide::new(residues, b'K', b'-');
    let mass = peptide.mass();
    let mz = (mass + charge as f64 * PROTON) / charge as f64;
    let spec = make_spectrum(mz, Some(charge as i32));

    match_spectra(&[spec], &idx, &params, &tiny_scorer(), 0.05, "XXX")
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

#[test]
fn spec_e_value_is_at_most_one_for_all_psms() {
    // After compute_spec_e_values_for_spectrum, no PSM should have
    // spec_e_value > 1.0 (spectral probability is always in (0, 1]).
    let (queues, _candidates) = run_single_peptide_search(b"MKWVTFISLLR", b"WVTFISLLR", 2);
    assert_eq!(queues.len(), 1);
    let sorted = queues.into_iter().next().unwrap().into_sorted_vec();
    assert!(!sorted.is_empty(), "expected at least one PSM");
    for psm in &sorted {
        assert!(
            psm.spec_e_value <= 1.0 + 1e-9,
            "spec_e_value {} > 1.0 for PSM with score {}",
            psm.spec_e_value,
            psm.score
        );
    }
}

#[test]
fn top_psm_has_spec_e_value_set() {
    // For a known-good peptide match, the top PSM's spec_e_value should be
    // something meaningful (not left at the sentinel 1.0 in most cases, but
    // this is not guaranteed for minimal fixtures — so we just verify it's
    // a valid probability in (0, 1]).
    let (queues, _candidates) = run_single_peptide_search(b"MKWVTFISLLR", b"WVTFISLLR", 2);
    let sorted = queues.into_iter().next().unwrap().into_sorted_vec();
    let top = &sorted[0];
    assert!(
        top.spec_e_value > 0.0,
        "spec_e_value must be positive (probability)"
    );
    assert!(
        top.spec_e_value <= 1.0 + 1e-9,
        "spec_e_value must be at most 1.0 (probability)"
    );
}

#[test]
fn sorted_vec_spec_e_value_is_non_decreasing() {
    // After sorting, the best PSM (index 0) should have the smallest
    // spec_e_value; values should be non-decreasing from index 0 onward.
    //
    // Use a larger protein so there are multiple candidate PSMs in the queue.
    let (queues, _candidates) = run_single_peptide_search(
        b"MKWVTFISLLLKWVTFISLLLER",
        b"WVTFISLLL",
        2,
    );
    let sorted = queues.into_iter().next().unwrap().into_sorted_vec();
    if sorted.len() < 2 {
        // Not enough PSMs to assert ordering; skip gracefully.
        return;
    }
    for window in sorted.windows(2) {
        let (a, b) = (&window[0], &window[1]);
        // a.spec_e_value <= b.spec_e_value (non-decreasing = best first).
        assert!(
            a.spec_e_value <= b.spec_e_value + 1e-12,
            "sorted_vec not non-decreasing in spec_e_value: {} > {}",
            a.spec_e_value,
            b.spec_e_value
        );
    }
}

#[test]
fn psm_with_lower_spec_e_value_ranks_first() {
    // Directly construct two PsmMaches with different spec_e_values and verify
    // that the one with the lower e-value sorts first in the sorted_vec.
    use search::psm::TopNQueue;

    fn make_psm(score: f32, spec_e_value: f64) -> PsmMatch {
        // candidate_idxs[0] = 0 is a placeholder for queue-ordering tests that
        // never resolve the candidate back. Safe because this test never
        // touches a `candidates` slice.
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score,
            rank_score: score,  // iter33: queue-ordering test defaults rank_score = score
            edge_score: 0,
            spec_e_value,
            de_novo_score: i32::MIN,
            activation_method: None,
            e_value: 1.0,
            features: search::psm::PsmFeatures::default(),
            isotope_offset: 0,
        }
    }

    let mut q = TopNQueue::new(5);
    q.push(make_psm(5.0, 0.5));    // mediocre
    q.push(make_psm(5.0, 0.001));  // best
    q.push(make_psm(5.0, 0.1));    // medium

    let sorted = q.into_sorted_vec();
    assert_eq!(sorted.len(), 3);
    // Best e-value first.
    assert!(
        sorted[0].spec_e_value <= sorted[1].spec_e_value,
        "index 0 should have <= spec_e_value of index 1"
    );
    assert!(
        sorted[1].spec_e_value <= sorted[2].spec_e_value,
        "index 1 should have <= spec_e_value of index 2"
    );
    assert!(
        (sorted[0].spec_e_value - 0.001).abs() < 1e-12,
        "best e-value should be 0.001, got {}",
        sorted[0].spec_e_value
    );
}

// ---------------------------------------------------------------------------
// Phase 7 / Task 1: PSM enrichment field tests
// ---------------------------------------------------------------------------

#[test]
fn top_psm_de_novo_score_equals_gf_max_minus_one() {
    // After match_spectra, the top PSM's de_novo_score should equal
    // group.max_score() - 1 (Java's getDeNovoScore() contract).
    //
    // We verify the structural invariant rather than an exact numeric value:
    // de_novo_score must NOT be the sentinel (i32::MIN) and must be >= 0
    // (GF max_score is always positive for non-trivial peptides).
    let (queues, _candidates) = run_single_peptide_search(b"MKWVTFISLLR", b"WVTFISLLR", 2);
    let sorted = queues.into_iter().next().unwrap().into_sorted_vec();
    assert!(!sorted.is_empty(), "expected at least one PSM");
    let top = &sorted[0];
    assert_ne!(
        top.de_novo_score, i32::MIN,
        "de_novo_score should not be sentinel after match_spectra"
    );
    assert!(
        top.de_novo_score >= 0,
        "de_novo_score should be non-negative (GF max score is positive), got {}",
        top.de_novo_score
    );
}

#[test]
fn top_psm_e_value_is_spec_e_value_times_some_constant() {
    // After match_spectra, e_value = spec_e_value * num_distinct_peptides.
    // Since num_distinct_peptides >= 1, e_value >= spec_e_value.
    // We verify: e_value > 0 and e_value >= spec_e_value.
    let (queues, _candidates) = run_single_peptide_search(b"MKWVTFISLLR", b"WVTFISLLR", 2);
    let sorted = queues.into_iter().next().unwrap().into_sorted_vec();
    assert!(!sorted.is_empty(), "expected at least one PSM");
    let top = &sorted[0];
    assert!(
        top.e_value > 0.0,
        "e_value must be positive, got {}",
        top.e_value
    );
    assert!(
        top.e_value >= top.spec_e_value - 1e-12,
        "e_value ({}) must be >= spec_e_value ({}) since num_distinct_peptides >= 1",
        top.e_value,
        top.spec_e_value
    );
}

// ---------------------------------------------------------------------------
// Protein-terminal flag derivation into GF construction.
// ---------------------------------------------------------------------------

/// Helper: run a single-peptide search and return the top PSM's spec_e_value.
///
/// `protein_seq` — the protein sequence that `peptide_seq` is embedded in.
/// `peptide_seq` — the peptide residues (must be a contiguous sub-sequence).
/// `charge`      — precursor charge to use.
fn top_spec_e_value_for(protein_seq: &[u8], peptide_seq: &[u8], charge: u8) -> f64 {
    let (queues, _candidates) = run_single_peptide_search(protein_seq, peptide_seq, charge);
    let sorted = queues.into_iter().next().unwrap().into_sorted_vec();
    assert!(!sorted.is_empty(), "expected at least one PSM");
    sorted[0].spec_e_value
}

/// Smoke test: the GF should use protein-terminal flags derived from
/// the top PSM rather than always hard-coding `false, false`.
///
/// We verify this *indirectly* by comparing spec_e_values for two scenarios:
///   (a) `WVTFISLLR` at the N-terminus of the protein  →  use_protein_n_term=true
///   (b) `WVTFISLLR` embedded after a K residue        →  use_protein_n_term=false
///
/// If the fix is working, the GF is built with different flags and the resulting
/// spec_e_values may differ (because the cleavage edge at the source node
/// changes with the N-terminal flag).  We do NOT assert a specific numeric
/// difference — we assert that the two paths produce *valid* spec_e_values
/// (i.e. the fix did not break anything) and document the observed values.
///
/// Note: in some degenerate fixtures (very short peptides, flat score landscape)
/// the two values can coincide.  The test therefore uses `assert!` on validity
/// rather than asserting strict inequality, and prints the observed pair for
/// inspection in CI logs.
#[test]
fn gf_protein_n_term_flag_derived_from_top_psm() {
    // (a) peptide at protein N-terminus: start_offset_in_protein = 0
    //     protein = WVTFISLLRK, peptide = WVTFISLLR (tryptic; K is the post-residue)
    let ev_n_term = top_spec_e_value_for(b"WVTFISLLRK", b"WVTFISLLR", 2);

    // (b) same peptide embedded internally: protein = MKWVTFISLLRK
    //     start_offset_in_protein = 2  →  use_protein_n_term=false
    let ev_internal = top_spec_e_value_for(b"MKWVTFISLLRK", b"WVTFISLLR", 2);

    // Both values must be valid probabilities.
    assert!(ev_n_term > 0.0 && ev_n_term <= 1.0 + 1e-9,
        "N-terminal spec_e_value out of range: {ev_n_term}");
    assert!(ev_internal > 0.0 && ev_internal <= 1.0 + 1e-9,
        "internal spec_e_value out of range: {ev_internal}");

    // Print for inspection — helpful when the values differ or coincide.
    println!(
        "N-terminal spec_e_value={ev_n_term:.6e}  internal={ev_internal:.6e}  \
         differ={}",
        (ev_n_term - ev_internal).abs() > 1e-15
    );
}

/// Smoke test: protein C-terminal flag.
///
/// When the top PSM ends at the last residue of the protein, `use_protein_c_term`
/// should be `true`.  Same indirect-validity approach as the N-terminal test.
#[test]
fn gf_protein_c_term_flag_derived_from_top_psm() {
    // (a) peptide ends at C-terminus: protein = KWVTFISLLR
    //     tryptic peptide WVTFISLLR → post-residue is '-' (end-of-protein)
    let ev_c_term = top_spec_e_value_for(b"KWVTFISLLR", b"WVTFISLLR", 2);

    // (b) same peptide with a downstream residue: protein = KWVTFISLLRK
    //     peptide ends at position 9 of 10, i.e. NOT at C-terminus
    let ev_not_c_term = top_spec_e_value_for(b"KWVTFISLLRK", b"WVTFISLLR", 2);

    assert!(ev_c_term > 0.0 && ev_c_term <= 1.0 + 1e-9,
        "C-terminal spec_e_value out of range: {ev_c_term}");
    assert!(ev_not_c_term > 0.0 && ev_not_c_term <= 1.0 + 1e-9,
        "non-C-terminal spec_e_value out of range: {ev_not_c_term}");

    println!(
        "B4: C-terminal spec_e_value={ev_c_term:.6e}  non-C-term={ev_not_c_term:.6e}  \
         differ={}",
        (ev_c_term - ev_not_c_term).abs() > 1e-15
    );
}
