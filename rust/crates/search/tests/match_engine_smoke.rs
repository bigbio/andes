//! match_engine smoke tests.

use std::collections::HashMap;

use model::{AminoAcid, AminoAcidSetBuilder, Peptide, Protein, ProteinDb, Spectrum, PROTON, Tolerance};
use scoring_crate::{Param, RankScorer};
use search::{match_spectra, SearchIndex, SearchParams};
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use scoring_crate::param_model::{IonType, Partition, SpecDataType};
use model::protocol::Protocol;

fn make_spectrum(precursor_mz: f64, charge: Option<i32>) -> Spectrum {
    Spectrum {
        title: "smoke".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: charge,
        rt_seconds: None,
        scan: None,
        peaks: vec![],
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

#[test]
fn known_peptide_appears_in_top_n() {
    // Protein "MKWVTFISLLR" — Trypsin cleaves after K (pos 1) and R (pos 10).
    // Peptide "WVTFISLLR" (positions 2..11, length 9) is a perfect cleavage.
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MKWVTFISLLR".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);

    let target_residues: Vec<AminoAcid> = b"WVTFISLLR".iter()
        .map(|&r| AminoAcid::standard(r).unwrap()).collect();
    let target_peptide = Peptide::new(target_residues, b'K', b'-');
    let target_mass = target_peptide.mass();
    let charge = 2u8;
    let mz = (target_mass + charge as f64 * PROTON) / charge as f64;

    let spec = make_spectrum(mz, Some(charge as i32));
    let (queues, candidates) = match_spectra(&[spec], &idx, &params, &tiny_scorer(), 0.05, "XXX");

    assert_eq!(queues.len(), 1);
    let top = queues.into_iter().next().unwrap().into_sorted_vec();
    assert!(!top.is_empty(), "expected at least one match");
    let best = &top[0];
    assert_eq!(candidates[best.candidate_idx as usize].peptide.length(), 9);
    assert!(!candidates[best.candidate_idx as usize].is_decoy);
    assert!(best.mass_error_ppm.abs() < 1.0);
}

#[test]
fn top_n_capacity_respected() {
    // NoCleavage gives exactly 1 candidate per protein. Top-N cap at 1.
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"AAAAAAAAAA".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.enzyme = model::Enzyme::NoCleavage;
    params.top_n_psms_per_spectrum = 1;
    params.max_variable_mods_per_peptide = 0;

    let target_residues: Vec<AminoAcid> = b"AAAAAAAAAA".iter()
        .map(|&r| AminoAcid::standard(r).unwrap()).collect();
    let target_peptide = Peptide::new(target_residues, b'_', b'-');
    let mass = target_peptide.mass();
    let charge = 2u8;
    let mz = (mass + charge as f64 * PROTON) / charge as f64;

    let spec = make_spectrum(mz, Some(charge as i32));
    let (queues, candidates) = match_spectra(&[spec], &idx, &params, &tiny_scorer(), 0.05, "XXX");
    assert!(queues[0].len() <= 1);
}

#[test]
fn spectrum_without_charge_tries_charge_range() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MKWVTFISLLR".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);

    let target_residues: Vec<AminoAcid> = b"WVTFISLLR".iter()
        .map(|&r| AminoAcid::standard(r).unwrap()).collect();
    let target_peptide = Peptide::new(target_residues, b'K', b'-');
    let mass = target_peptide.mass();
    let charge = 2u8;
    let mz = (mass + charge as f64 * PROTON) / charge as f64;

    let spec = make_spectrum(mz, None);  // no charge!
    let (queues, candidates) = match_spectra(&[spec], &idx, &params, &tiny_scorer(), 0.05, "XXX");
    let top = queues.into_iter().next().unwrap().into_sorted_vec();
    assert!(!top.is_empty(), "expected charge_range to find a match");
    assert_eq!(top[0].charge_used, 2);
}

/// B3 correctness: for charge-missing spectra, each candidate is scored
/// against a ScoredSpectrum built with its own charge (not a fixed z=2).
///
/// We set up a peptide whose precursor m/z at z=3 matches the spectrum
/// but at z=2 does not.  With the pre-B3 code (single scored_spec at z=2)
/// the candidate would still be found but with a mismatched charge.
/// With the B3 fix (per-charge cache), each charge sees its own ScoredSpectrum
/// and the PSM's charge_used matches the charge that actually satisfied the
/// precursor-mass check.
#[test]
fn charge_missing_spectrum_uses_per_charge_scored_spec() {
    // Peptide "WVTFISLLR", a tryptic fragment from BSA-related sequences.
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MKWVTFISLLR".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    // charge_range 2..=3; spectrum has no charge.
    params.charge_range = 2..=3;

    let target_residues: Vec<AminoAcid> = b"WVTFISLLR".iter()
        .map(|&r| AminoAcid::standard(r).unwrap()).collect();
    let target_peptide = Peptide::new(target_residues, b'K', b'-');
    let mass = target_peptide.mass();

    // Set the precursor m/z at z=3 so only z=3 satisfies precursor matching.
    let charge = 3u8;
    let mz = (mass + charge as f64 * PROTON) / charge as f64;

    let spec = make_spectrum(mz, None);  // charge-missing
    let (queues, candidates) = match_spectra(&[spec], &idx, &params, &tiny_scorer(), 0.05, "XXX");
    let top = queues.into_iter().next().unwrap().into_sorted_vec();

    // The only match must be at charge 3 (the precursor m/z is z=3-exact).
    assert!(!top.is_empty(), "expected a charge-3 match for charge-missing spectrum");
    assert!(
        top.iter().all(|p| p.charge_used == 3),
        "all PSMs should be at z=3; found charges: {:?}",
        top.iter().map(|p| p.charge_used).collect::<Vec<_>>()
    );
}
