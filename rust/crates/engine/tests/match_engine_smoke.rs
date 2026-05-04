//! match_engine smoke tests.

use engine::{
    match_spectra, AminoAcid, AminoAcidSetBuilder, Peptide, Protein, ProteinDb,
    SearchIndex, SearchParams, Spectrum, PROTON,
};

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
    let queues = match_spectra(&[spec], &idx, &params, "XXX");

    assert_eq!(queues.len(), 1);
    let top = queues.into_iter().next().unwrap().into_sorted_vec();
    assert!(!top.is_empty(), "expected at least one match");
    let best = &top[0];
    assert_eq!(best.candidate.peptide.length(), 9);
    assert!(!best.candidate.is_decoy);
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
    params.enzyme = engine::Enzyme::NoCleavage;
    params.top_n_psms_per_spectrum = 1;
    params.max_variable_mods_per_peptide = 0;

    let target_residues: Vec<AminoAcid> = b"AAAAAAAAAA".iter()
        .map(|&r| AminoAcid::standard(r).unwrap()).collect();
    let target_peptide = Peptide::new(target_residues, b'_', b'-');
    let mass = target_peptide.mass();
    let charge = 2u8;
    let mz = (mass + charge as f64 * PROTON) / charge as f64;

    let spec = make_spectrum(mz, Some(charge as i32));
    let queues = match_spectra(&[spec], &idx, &params, "XXX");
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
    let queues = match_spectra(&[spec], &idx, &params, "XXX");
    let top = queues.into_iter().next().unwrap().into_sorted_vec();
    assert!(!top.is_empty(), "expected charge_range to find a match");
    assert_eq!(top[0].charge_used, 2);
}
