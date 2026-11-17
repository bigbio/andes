use model::{
    AminoAcidSet, AminoAcidSetBuilder, Enzyme, ModLocation, Modification, Protein, ProteinDb,
    ResidueSpec,
};
use search::{enumerate_candidates, SearchIndex, SearchParams};

fn aa_set() -> AminoAcidSet {
    AminoAcidSetBuilder::new_standard().build().unwrap()
}

fn make_index(seq: &[u8]) -> SearchIndex {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: seq.to_vec(),
            annotation: None,
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
    assert_eq!(
        target_candidates.len(),
        1,
        "expected 1 target candidate, got {}",
        target_candidates.len()
    );
    for cand in &target_candidates {
        assert_eq!(cand.peptide.length(), 9);
        assert_eq!(cand.start_offset_in_protein, 2);
    }
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
            annotation: None,
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
    assert_eq!(candidates.len(), 4);
    let target_candidates: Vec<_> = candidates.iter().filter(|c| !c.is_decoy).collect();
    assert_eq!(target_candidates.len(), 2);
    let full = target_candidates
        .iter()
        .find(|c| c.start_offset_in_protein == 0)
        .unwrap();
    assert_eq!(full.peptide.length(), 11);
    let met_cleaved = target_candidates
        .iter()
        .find(|c| c.start_offset_in_protein == 1)
        .unwrap();
    assert_eq!(met_cleaved.peptide.length(), 10);
}

#[test]
fn nonspecific_enzyme_emits_every_length_valid_span() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"AAAAAA".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.enzyme = Enzyme::NonSpecific;
    p.min_length = 3;
    p.max_length = 6;
    p.max_variable_mods_per_peptide = 0;
    for missed in [0, 1, 5] {
        p.max_missed_cleavages = missed;
        let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX").collect();
        let target_candidates: Vec<_> = candidates.iter().filter(|c| !c.is_decoy).collect();
        assert_eq!(target_candidates.len(), 10, "max_missed_cleavages {missed}");
    }
}

#[test]
fn nonspecific_extra_enzyme_ignores_missed_cleavage_cap() {
    let idx = make_index(b"AAAAAAAAAAAA");
    let mut p = params(8, 11, 1);
    p.enzyme = Enzyme::Trypsin;
    p.extra_enzymes = vec![Enzyme::NonSpecific];
    let n = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert_eq!(n, 14);
}

#[test]
fn nonspecific_default_missed_cap_builds_candidates() {
    let idx = make_index(b"AAAAAAAAAAAAAAAAAAAA");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.enzyme = Enzyme::NonSpecific;
    p.min_length = 8;
    p.max_length = 11;
    p.max_variable_mods_per_peptide = 0;
    assert_eq!(p.max_missed_cleavages, 1);
    let n = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert!(n > 0, "no candidates for a non-specific digest");
}

#[test]
fn missed_cleavages_increase_candidate_count() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"AKMKCKDK".to_vec(),
            annotation: None,
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

    assert!(
        c0_count < c1_count,
        "missed=0 ({c0_count}) should be less than missed=1 ({c1_count})"
    );
    assert!(
        c1_count < c2_count,
        "missed=1 ({c1_count}) should be less than missed=2 ({c2_count})"
    );
}

#[test]
fn missed_cleavages_zero_emits_only_perfectly_cleaved() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"AKMKLR".to_vec(),
            annotation: None,
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
    assert_eq!(
        target_count, 3,
        "expected 3 perfectly-cleaved peptides, got {target_count}"
    );
}

fn aa_set_with_oxidation() -> model::AminoAcidSet {
    let ox = Modification {
        name: "Oxidation".into(),
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
    };
    model::AminoAcidSetBuilder::new_standard()
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

#[test]
fn one_variable_mod_site_doubles_candidates() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKAR".to_vec(),
            annotation: None,
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
    assert_eq!(
        target_count, 3,
        "expected 3 target candidates (MK + MKox + AR)"
    );
}

#[test]
fn two_variable_mod_sites_quadruple_candidates() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MMK".to_vec(),
            annotation: None,
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
    assert_eq!(
        target_count, 6,
        "expected 6 (MMK×4 standard + MK×2 met-cleaved)"
    );
}

#[test]
fn max_variable_mods_caps_combinations() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MMMK".to_vec(),
            annotation: None,
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
    assert_eq!(
        target_count, 7,
        "expected 7 (MMMK×4 standard + MMK×3 met-cleaved)"
    );
}

fn aa_set_with_protein_nterm_acetyl() -> AminoAcidSet {
    let acetyl = Modification {
        name: "ProtNTermAcetyl".into(),
        mass_delta: 42.010565,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::ProtNTerm,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    AminoAcidSetBuilder::new_standard()
        .add_variable_mod(acetyl)
        .build()
        .unwrap()
}

fn aa_set_with_nterm_acetyl() -> AminoAcidSet {
    let acetyl = Modification {
        name: "NTermAcetyl".into(),
        mass_delta: 42.010565,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::NTerm,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    AminoAcidSetBuilder::new_standard()
        .add_variable_mod(acetyl)
        .build()
        .unwrap()
}

fn aa_set_with_both_cterm_mods() -> AminoAcidSet {
    let cterm = Modification {
        name: "Amide_CT".into(),
        mass_delta: -0.984016,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::CTerm,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    let prot_cterm = Modification {
        name: "GlyGly_PCT".into(),
        mass_delta: 114.042927,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::ProtCTerm,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    AminoAcidSetBuilder::new_standard()
        .add_variable_mod(cterm)
        .add_variable_mod(prot_cterm)
        .build()
        .unwrap()
}

#[test]
fn protein_n_term_mod_only_at_protein_start() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MAAAAKMAAAAAK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set_with_protein_nterm_acetyl());
    p.min_length = 6;
    p.max_length = 40;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 1;

    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();

    assert_eq!(
        candidates.len(),
        3,
        "expected 3 candidates (2 for protein-start peptide, 1 for offset-6 peptide), got {}",
        candidates.len()
    );

    for cand in &candidates {
        let has_mod = cand.peptide.residues[0].is_modified();
        if has_mod {
            assert_eq!(
                cand.start_offset_in_protein, 0,
                "ProtNTerm mod appeared on peptide starting at offset {} (should only be at 0)",
                cand.start_offset_in_protein
            );
        }
    }

    let mod_count = candidates
        .iter()
        .filter(|c| c.peptide.residues[0].is_modified())
        .count();
    assert_eq!(
        mod_count, 1,
        "exactly 1 candidate should have the ProtNTerm mod"
    );
}

#[test]
fn nterm_mod_applies_to_every_peptide_n_terminus() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"AAAAAAKMAAAAAK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set_with_nterm_acetyl());
    p.min_length = 7;
    p.max_length = 40;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 1;

    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();

    assert_eq!(
        candidates.len(), 4,
        "expected 4 candidates (unmod + NTerm-acetyl for both the offset-0 and offset-7 peptides), got {}",
        candidates.len()
    );

    let modified_offsets: std::collections::BTreeSet<usize> = candidates
        .iter()
        .filter(|c| c.peptide.residues[0].is_modified())
        .map(|c| c.start_offset_in_protein)
        .collect();
    assert_eq!(
        modified_offsets, std::collections::BTreeSet::from([0usize, 7usize]),
        "NTerm acetyl should appear on both peptide N-termini (offsets 0 and 7), got {modified_offsets:?}"
    );

    for cand in &candidates {
        let residues = &cand.peptide.residues;
        for (i, aa) in residues.iter().enumerate().skip(1) {
            assert!(
                !aa.is_modified(),
                "NTerm acetyl leaked to internal position {i} in peptide at offset {}",
                cand.start_offset_in_protein
            );
        }
    }
}

#[test]
fn c_term_and_protein_c_term_distinguished() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MAAAAKR".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set_with_both_cterm_mods());
    p.min_length = 1;
    p.max_length = 40;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 1;

    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();

    assert_eq!(
        candidates.len(),
        7,
        "expected 7 candidates, got {}",
        candidates.len()
    );

    let protein_len = 7usize;
    for cand in &candidates {
        let span_end = cand.start_offset_in_protein + cand.peptide.length();
        let is_prot_c_term = span_end == protein_len;
        let residues = &cand.peptide.residues;
        if let Some(last) = residues.last() {
            if let Some(m) = &last.mod_ {
                if is_prot_c_term {
                    assert!(
                        (m.mass_delta - (-0.984016)).abs() < 1e-4
                            || (m.mass_delta - 114.042927).abs() < 1e-4,
                        "protein C-term peptide got an unexpected mod delta ({}); expected CTerm Amide (-0.984016) or ProtCTerm GlyGly (114.042927)",
                        m.mass_delta
                    );
                } else {
                    assert!(
                        (m.mass_delta - (-0.984016)).abs() < 1e-4,
                        "non-protein-C-term peptide got an unexpected delta mod ({}); expected CTerm Amide (-0.984016)",
                        m.mass_delta
                    );
                }
            }
        }
    }
}

#[test]
fn met_cleavage_generates_alternative_candidates() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MAGER".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.enzyme = Enzyme::NoCleavage;
    p.min_length = 1;
    p.max_length = 40;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 0;

    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();

    assert_eq!(
        candidates.len(),
        2,
        "expected 2 target candidates (standard + Met-cleaved), got {}",
        candidates.len()
    );

    let has_full = candidates
        .iter()
        .any(|c| c.start_offset_in_protein == 0 && c.peptide.length() == 5);
    let has_met_cleaved = candidates
        .iter()
        .any(|c| c.start_offset_in_protein == 1 && c.peptide.length() == 4);

    assert!(has_full, "missing standard candidate at offset 0 (MAGER)");
    assert!(
        has_met_cleaved,
        "missing Met-cleaved candidate at offset 1 (AGER)"
    );
}

#[test]
fn non_met_first_residue_does_not_trigger_cleavage() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"KAGER".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.enzyme = Enzyme::NoCleavage;
    p.min_length = 1;
    p.max_length = 40;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 0;

    let target_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();

    assert_eq!(
        target_count, 1,
        "expected 1 candidate for non-M protein, got {}",
        target_count
    );
}

const NTT_PROTEIN: &[u8] = b"AAAAAKAAAAR";

fn ntt_protein_index() -> SearchIndex {
    make_index(NTT_PROTEIN)
}

fn ntt_params(ntt: u8) -> SearchParams {
    let mut p = params(5, 11, 0);
    p.num_tolerable_termini = ntt;
    p
}

#[test]
fn ntt_2_emits_only_strict_tryptic_spans() {
    let idx = ntt_protein_index();
    let p = ntt_params(2);
    let count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert_eq!(
        count, 2,
        "ntt=2 should emit exactly 2 strict tryptic spans, got {count}"
    );
}

#[test]
fn ntt_1_emits_strict_plus_semi_spans() {
    let idx = ntt_protein_index();
    let ntt2_count = enumerate_candidates(&idx, &ntt_params(2), "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    let ntt1_count = enumerate_candidates(&idx, &ntt_params(1), "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert!(
        ntt1_count > ntt2_count,
        "ntt=1 ({ntt1_count}) should generate more candidates than ntt=2 ({ntt2_count})"
    );
    assert_eq!(
        ntt1_count, 4,
        "expected 4 ntt=1 candidates, got {ntt1_count}"
    );
}

#[test]
fn ntt_1_includes_free_c_term_span() {
    let idx = ntt_protein_index();
    let p = ntt_params(1);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();
    let has_free_c = candidates.iter().any(|c| {
        c.start_offset_in_protein == 0 && c.peptide.length() == 5
    });
    assert!(
        has_free_c,
        "ntt=1 should include (start=0, end=5): tryptic N-term, free C-term"
    );
}

#[test]
fn ntt_1_includes_free_n_term_span() {
    let idx = ntt_protein_index();
    let p = ntt_params(1);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();
    let has_free_n = candidates
        .iter()
        .any(|c| c.start_offset_in_protein == 1 && c.peptide.length() == 5);
    assert!(
        has_free_n,
        "ntt=1 should include (start=1, end=6): free N-term, tryptic C-term"
    );
}

#[test]
fn ntt_1_no_dedup_for_strict_spans() {
    let idx = ntt_protein_index();
    let p = ntt_params(1);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();
    let count_strict = candidates
        .iter()
        .filter(|c| c.start_offset_in_protein == 0 && c.peptide.length() == 6)
        .count();
    assert_eq!(
        count_strict, 1,
        "strict span (0,6) should appear exactly once under ntt=1, got {count_strict}"
    );
}

#[test]
fn ntt_0_emits_all_spans() {
    let idx = ntt_protein_index();
    let ntt1_count = enumerate_candidates(&idx, &ntt_params(1), "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    let ntt0_count = enumerate_candidates(&idx, &ntt_params(0), "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    assert!(
        ntt0_count > ntt1_count,
        "ntt=0 ({ntt0_count}) should generate more candidates than ntt=1 ({ntt1_count})"
    );
    assert_eq!(
        ntt0_count, 28,
        "ntt=0 should emit all 28 valid-length spans, got {ntt0_count}"
    );
}

#[test]
fn ntt_0_trypsin_matches_nonspecific_high_missed() {
    let seq = b"AAAAAAAAAAAA";
    let idx = make_index(seq);

    let mut p_ntt0 = params(3, 8, 10);
    p_ntt0.enzyme = Enzyme::Trypsin;
    p_ntt0.num_tolerable_termini = 0;

    let mut p_ns = params(3, 8, 10);
    p_ns.enzyme = Enzyme::NonSpecific;
    p_ns.num_tolerable_termini = 2;

    let ntt0_count = enumerate_candidates(&idx, &p_ntt0, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    let ns_count = enumerate_candidates(&idx, &p_ns, "XXX")
        .filter(|c| !c.is_decoy)
        .count();

    assert_eq!(
        ntt0_count, 45,
        "ntt=0 + Trypsin should emit 45 spans for AAAAAAAAAAAA min=3 max=8, got {ntt0_count}"
    );
    assert_eq!(
        ns_count, 45,
        "NonSpecific + ntt=2 high missed should also emit 45 spans, got {ns_count}"
    );
}

#[test]
fn default_ntt_is_2() {
    let p = SearchParams::default_tryptic(aa_set());
    assert_eq!(p.num_tolerable_termini, 2, "default ntt should be 2");
}

#[test]
fn met_alone_does_not_trigger_cleavage() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"M".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = SearchParams::default_tryptic(aa_set());
    p.enzyme = Enzyme::NoCleavage;
    p.min_length = 1;
    p.max_length = 40;
    p.max_missed_cleavages = 0;
    p.max_variable_mods_per_peptide = 0;

    let target_count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();

    assert_eq!(
        target_count, 1,
        "expected 1 candidate for M-only protein, got {}",
        target_count
    );
}

#[test]
fn empty_decoy_prefix_labels_only_normalized_decoys() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLR".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "");
    let p = params(6, 40, 0);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "").collect();

    assert!(
        candidates.iter().any(|c| !c.is_decoy),
        "target proteins must not be labeled decoy when --decoy-prefix is empty"
    );
    assert!(
        candidates.iter().any(|c| c.is_decoy),
        "decoy proteins must still be labeled decoy when --decoy-prefix is empty"
    );
    assert!(
        candidates
            .iter()
            .all(|c| !c.is_decoy || c.protein_index >= target.proteins.len()),
        "only decoy half of the index may carry is_decoy=true"
    );
}
