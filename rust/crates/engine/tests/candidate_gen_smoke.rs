//! Handcrafted candidate-enumeration tests.

use engine::{
    enumerate_candidates, AminoAcidSet, AminoAcidSetBuilder, Enzyme,
    ModLocation, Modification, Protein, ProteinDb, ResidueSpec,
    SearchIndex, SearchParams,
};

fn aa_set() -> AminoAcidSet {
    AminoAcidSetBuilder::new_standard().build().unwrap()
}

fn make_index(seq: &[u8]) -> SearchIndex {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: seq.to_vec(),
        }],
    };
    SearchIndex::from_target_db(&target, "XXX")
}

fn params(min: u32, max: u32, missed: u32) -> SearchParams {
    let mut p = SearchParams::default_tryptic(aa_set());
    p.min_length = min;
    p.max_length = max;
    p.max_missed_cleavages = missed;
    p.max_variable_mods_per_peptide = 0;
    p
}

#[test]
fn single_tryptic_peptide_no_missed() {
    let idx = make_index(b"MKWVTFISLLR");
    let p = params(6, 40, 0);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX").collect();
    let target_candidates: Vec<_> = candidates.iter().filter(|c| !c.is_decoy).collect();
    assert_eq!(target_candidates.len(), 1, "expected 1 target candidate, got {}", target_candidates.len());
    let cand = target_candidates[0];
    assert_eq!(cand.peptide.length(), 9);
    assert_eq!(cand.start_offset_in_protein, 2);
}

#[test]
fn protein_shorter_than_min_yields_nothing() {
    let idx = make_index(b"AB");
    let p = params(6, 40, 0);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX").collect();
    assert!(candidates.is_empty());
}

#[test]
fn each_candidate_is_decoy_or_target() {
    let idx = make_index(b"MKWVTFISLLR");
    let p = params(6, 40, 0);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX").collect();
    assert!(candidates.iter().any(|c| !c.is_decoy));
    assert!(candidates.iter().any(|c| c.is_decoy));
}

#[test]
fn no_cleavage_enzyme_emits_full_protein_only() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLR".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.enzyme = Enzyme::NoCleavage;
    p.min_length = 6;
    p.max_length = 40;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 0;
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX").collect();
    // 2 candidates: target whole + decoy whole
    assert_eq!(candidates.len(), 2);
    assert_eq!(candidates[0].peptide.length(), 11);
    assert_eq!(candidates[0].start_offset_in_protein, 0);
}

#[test]
fn nonspecific_enzyme_emits_every_length_valid_span() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"AAAAAA".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.enzyme = Enzyme::NonSpecific;
    p.min_length = 3;
    p.max_length = 6;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 0;
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX").collect();
    let _target_candidates: Vec<_> = candidates.iter().filter(|c| !c.is_decoy).collect();
    // For NonSpecific, every cleavage position can pair. With seq length 6
    // and missed=0, only ADJACENT cleavage positions form candidates.
    // Cleavage positions = [0, 1, 2, 3, 4, 5, 6]; adjacent spans have length 1.
    // None match length range 3-6, so 0 candidates with missed=0.
    // Wait — that's wrong. Re-read the spec: missed cleavages means count
    // of cleavage positions strictly between start and end. For NonSpecific
    // every position is cleavable, so a length-3 span (start, start+3) has
    // 2 internal cleavage positions, requiring missed_cleavages >= 2.
    //
    // So with missed=0 and NonSpecific, no length>1 spans are valid.
    // Re-do: change params to missed=5 (high enough to allow any).
    p.max_missed_cleavages = 5;
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX").collect();
    let target_candidates: Vec<_> = candidates.iter().filter(|c| !c.is_decoy).collect();
    // length 3: 4 starts; length 4: 3; length 5: 2; length 6: 1; total 10.
    assert_eq!(target_candidates.len(), 10);
}

#[test]
fn missed_cleavages_increase_candidate_count() {
    // Sequence "AKMKCKDK" — Trypsin cleaves after K at positions 2, 4, 6, 8.
    // Cleavage positions: [0, 2, 4, 6, 8].
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"AKMKCKDK".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.min_length = 2;
    p.max_length = 8;
    p.max_variable_mods_per_peptide = 0;

    p.max_missed_cleavages = 0;
    let c0_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();

    p.max_missed_cleavages = 1;
    let c1_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();

    p.max_missed_cleavages = 2;
    let c2_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();

    assert!(c0_count < c1_count, "missed=0 ({c0_count}) should be less than missed=1 ({c1_count})");
    assert!(c1_count < c2_count, "missed=1 ({c1_count}) should be less than missed=2 ({c2_count})");
}

#[test]
fn missed_cleavages_zero_emits_only_perfectly_cleaved() {
    // "AKMKLR" — Trypsin cleaves after positions 1 (K), 3 (K), 5 (R).
    // Cleavage positions: [0, 2, 4, 6].
    // missed=0, length 2-6: spans (0,2)="AK", (2,4)="MK", (4,6)="LR" — 3 spans.
    // (Note: 'B' is not standard so we use 'L' which IS standard.)
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"AKMKLR".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.min_length = 2;
    p.max_length = 6;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 0;
    let target_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert_eq!(target_count, 3, "expected 3 perfectly-cleaved peptides, got {target_count}");
}

fn aa_set_with_oxidation() -> engine::AminoAcidSet {
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
    };
    engine::AminoAcidSetBuilder::new_standard()
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

#[test]
fn one_variable_mod_site_doubles_candidates() {
    // "MKAR" — Trypsin spans (0,2)="MK" + (2,4)="AR".
    // With Oxidation-M variable: "MK" → 2 versions (unmod + Mox); "AR" → 1.
    // Total target = 3.
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MKAR".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set_with_oxidation());
    p.min_length = 2;
    p.max_length = 4;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 3;
    let target_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert_eq!(target_count, 3, "expected 3 target candidates (MK + MKox + AR)");
}

#[test]
fn two_variable_mod_sites_quadruple_candidates() {
    // "MMK" — single span (0,3) with 2 M positions.
    // Combos: {none, M0_ox, M1_ox, both_ox} = 4.
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MMK".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set_with_oxidation());
    p.min_length = 2;
    p.max_length = 5;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 3;
    let target_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert_eq!(target_count, 4);
}

#[test]
fn max_variable_mods_caps_combinations() {
    // "MMMK" — 3 M sites. With max_mods=1: {none, M0_ox, M1_ox, M2_ox} = 4.
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MMMK".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set_with_oxidation());
    p.min_length = 2;
    p.max_length = 5;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 1;
    let target_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert_eq!(target_count, 4);
}
