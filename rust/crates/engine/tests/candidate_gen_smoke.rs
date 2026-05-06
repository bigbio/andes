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

// ─── Track B2: terminal-mod expansion tests ───────────────────────────────────
//
// Terminal-location semantics in expand_mod_combinations:
//   - Peptide at protein start (start_offset == 0): position 0 gets ProtNTerm variants.
//   - Peptide NOT at protein start: position 0 gets NTerm variants.
//   - Peptide at protein end (end == protein_len): last position gets ProtCTerm variants.
//   - Peptide NOT at protein end: last position gets CTerm variants.
//
// This mirrors Java CandidatePeptideGrid.java:43 which routes the first residue to
// addProtNTermResidue (protein start) or addNTermResidue (non-protein start), and the
// last residue to addProtCTermResidue (protein end) or addCTermResidue (non-protein end).

/// Build an AminoAcidSet with a Protein_N_Term-only variable mod (+42.0106 Acetyl on *).
fn aa_set_with_protein_nterm_acetyl() -> AminoAcidSet {
    let acetyl = Modification {
        name: "ProtNTermAcetyl".into(),
        mass_delta: 42.010565,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::ProtNTerm,
        fixed: false,
        accession: None,
    };
    AminoAcidSetBuilder::new_standard()
        .add_variable_mod(acetyl)
        .build()
        .unwrap()
}

/// Build an AminoAcidSet with an N-Term-only variable mod (+42.0106 Acetyl on *).
fn aa_set_with_nterm_acetyl() -> AminoAcidSet {
    let acetyl = Modification {
        name: "NTermAcetyl".into(),
        mass_delta: 42.010565,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::NTerm,
        fixed: false,
        accession: None,
    };
    AminoAcidSetBuilder::new_standard()
        .add_variable_mod(acetyl)
        .build()
        .unwrap()
}

/// Build an AminoAcidSet with both a C-Term and a Protein_C_Term variable mod.
fn aa_set_with_both_cterm_mods() -> AminoAcidSet {
    let cterm = Modification {
        name: "Amide_CT".into(),
        mass_delta: -0.984016,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::CTerm,
        fixed: false,
        accession: None,
    };
    let prot_cterm = Modification {
        name: "GlyGly_PCT".into(),
        mass_delta: 114.042927,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::ProtCTerm,
        fixed: false,
        accession: None,
    };
    AminoAcidSetBuilder::new_standard()
        .add_variable_mod(cterm)
        .add_variable_mod(prot_cterm)
        .build()
        .unwrap()
}

/// Protein_N_Term mod appears on the peptide starting at protein index 0.
///
/// Protein: "MAAAAKMAAAAAK" (length 13).
/// Trypsin + missed=0 → (0..6)="MAAAAK" (protein N-term start) + (6..13)="MAAAAAK" (not at start).
/// With ProtNTerm Acetyl variable mod and max_mods=1:
/// - "MAAAAK" (protein start): gets Anywhere (unmod M) + ProtNTerm (Acetyl-M) → 2 candidates.
/// - "MAAAAAK" (offset 6, not protein start): gets only Anywhere (unmod M) → 1 candidate.
///
/// Total: 3. This proves B2's fix: without the terminal-loc expansion, we'd get 2 (both unmod).
#[test]
fn protein_n_term_mod_only_at_protein_start() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MAAAAKMAAAAAK".to_vec(),
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

    // Peptide at offset 0: 2 (unmod + Protein_N_Term Acetyl).
    // Peptide at offset 6: 1 (unmod only — ProtNTerm does NOT apply here; this position gets NTerm).
    assert_eq!(
        candidates.len(), 3,
        "expected 3 candidates (2 for protein-start peptide, 1 for offset-6 peptide), got {}",
        candidates.len()
    );

    // Only candidates starting at protein offset 0 may have the ProtNTerm mod.
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

    // Exactly 1 candidate has the Protein_N_Term mod.
    let mod_count = candidates.iter()
        .filter(|c| c.peptide.residues[0].is_modified())
        .count();
    assert_eq!(mod_count, 1, "exactly 1 candidate should have the ProtNTerm mod");
}

/// N-Term mod applies to peptides NOT at the protein N-terminus.
///
/// Protein: "AAAAAAKMAAAAAK" (length 14).
/// Trypsin + missed=0 → (0..7)="AAAAAAK" (protein N-term) + (7..14)="MAAAAAK" (not at start).
/// With NTerm Acetyl variable mod and max_mods=1:
/// - "AAAAAAK" (protein start, offset=0): ProtNTerm lookup → NTerm mod does NOT apply → 1 unmod.
/// - "MAAAAAK" (offset=7): NTerm lookup → NTerm Acetyl applies to position 0 → 2 variants.
///
/// Total: 3.
#[test]
fn nterm_mod_applies_to_non_protein_start_peptides() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"AAAAAAKMAAAAAK".to_vec(),
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

    // "AAAAAAK" (protein start): no NTerm mod (gets ProtNTerm which is empty) → 1.
    // "MAAAAAK" (offset 7): NTerm Acetyl applies → 2.
    // Total: 3.
    assert_eq!(
        candidates.len(), 3,
        "expected 3 candidates (1 for protein-start, 2 for offset-7 with NTerm mod), got {}",
        candidates.len()
    );

    // The modified candidate must be at offset 7 (non-protein-start).
    let modified: Vec<_> = candidates.iter()
        .filter(|c| c.peptide.residues[0].is_modified())
        .collect();
    assert_eq!(modified.len(), 1, "exactly 1 candidate should have the NTerm mod");
    assert_eq!(
        modified[0].start_offset_in_protein, 7,
        "NTerm mod should appear on the offset-7 peptide, not at offset 0"
    );

    // The NTerm mod must NOT appear at any internal position.
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

/// C-Term and Protein_C_Term mods are routed to the correct peptide.
///
/// Protein: "MAAAAKR" (length 7).
/// Trypsin cleaves after K(5): spans (0..6)="MAAAAK" (not protein C-term) and (6..7)="R" (protein C-term).
/// With both CTerm Amide (-0.984) and ProtCTerm GlyGly (+114.04) variable mods:
/// - "MAAAAK" (end < protein_len): CTerm Amide applies to last position → 2 variants.
/// - "R" (end == protein_len): ProtCTerm GlyGly applies to last position → 2 variants.
///
/// Total: 4.
///
/// This also verifies the C-Term mod does NOT bleed into the protein-C-term peptide, and vice versa.
#[test]
fn c_term_and_protein_c_term_distinguished() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MAAAAKR".to_vec(),
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

    // "MAAAAK" → 2 (unmod + CTerm Amide on last residue K).
    // "R"      → 2 (unmod + ProtCTerm GlyGly on last residue R).
    assert_eq!(
        candidates.len(), 4,
        "expected 4 candidates, got {}",
        candidates.len()
    );

    // Verify the right mod appears on the right peptide.
    let protein_len = 7usize;
    for cand in &candidates {
        let span_end = cand.start_offset_in_protein + cand.peptide.length();
        let is_prot_c_term = span_end == protein_len;
        let residues = &cand.peptide.residues;
        if let Some(last) = residues.last() {
            if let Some(m) = &last.mod_ {
                if is_prot_c_term {
                    // Protein-C-term peptide "R": should get ProtCTerm GlyGly (+114.04).
                    assert!(
                        m.mass_delta > 0.0,
                        "protein C-term peptide 'R' got a negative delta mod ({}); expected ProtCTerm GlyGly",
                        m.mass_delta
                    );
                } else {
                    // Non-protein-C-term peptide "MAAAAK": should get CTerm Amide (-0.984).
                    assert!(
                        m.mass_delta < 0.0,
                        "non-protein-C-term peptide 'MAAAAK' got a positive delta mod ({}); expected CTerm Amide",
                        m.mass_delta
                    );
                }
            }
        }
    }
}
