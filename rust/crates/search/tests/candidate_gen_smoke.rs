//! Handcrafted candidate-enumeration tests.

use model::{AminoAcidSet, AminoAcidSetBuilder, Enzyme, ModLocation, Modification, Protein, ProteinDb, ResidueSpec};
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
    // Protein "MKWVTFISLLR": trypsin cleaves after K (pos 1) → spans "MK" (too short) + "WVTFISLLR".
    // Standard pass: 1 candidate "WVTFISLLR" at offset 2.
    // Met-cleavage pass (sub_seq="KWVTFISLLR"): trypsin cleaves after K (sub_pos 0) →
    //   sub-spans "K" (too short) + "WVTFISLLR" at abs_offset=2. Adds 1 more candidate.
    // Total target candidates: 2.
    let idx = make_index(b"MKWVTFISLLR");
    let p = params(6, 40, 0);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX").collect();
    let target_candidates: Vec<_> = candidates.iter().filter(|c| !c.is_decoy).collect();
    assert_eq!(target_candidates.len(), 2, "expected 2 target candidates (standard + Met-cleaved), got {}", target_candidates.len());
    // Both candidates are "WVTFISLLR" at offset 2 — one from each enumeration pass.
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
    // Protein starts with M, so Met-cleaved pass also runs.
    // Standard pass: target "MKWVTFISLLR" (len=11, offset=0) + decoy "RLLSIFTFVKM" (len=11, offset=0).
    // Met-cleaved pass (target only, since decoy "RLLSIFTFVKM" starts with R):
    //   sub_seq "KWVTFISLLR" (len=10) → 1 candidate at offset=1.
    // Total: 3 (2 standard + 1 met-cleaved target).
    assert_eq!(candidates.len(), 3);
    let target_candidates: Vec<_> = candidates.iter().filter(|c| !c.is_decoy).collect();
    assert_eq!(target_candidates.len(), 2);
    // Standard target: full protein at offset 0, length 11.
    let full = target_candidates.iter().find(|c| c.start_offset_in_protein == 0).unwrap();
    assert_eq!(full.peptide.length(), 11);
    // Met-cleaved target: sequence[1..] at offset 1, length 10.
    let met_cleaved = target_candidates.iter().find(|c| c.start_offset_in_protein == 1).unwrap();
    assert_eq!(met_cleaved.peptide.length(), 10);
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

fn aa_set_with_oxidation() -> model::AminoAcidSet {
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
    };
    model::AminoAcidSetBuilder::new_standard()
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

#[test]
fn one_variable_mod_site_doubles_candidates() {
    // "MKAR" — Trypsin spans (0,2)="MK" + (2,4)="AR".
    // Standard pass: "MK" → 2 (unmod + Mox); "AR" → 1. Total = 3.
    // Met-cleavage pass (sub_seq="KAR"): spans "K" (too short) + "AR" at abs_offset=2.
    //   "AR" has no M residue → 1 extra candidate.
    // Total target = 4.
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
    assert_eq!(target_count, 4, "expected 4 target candidates (MK + MKox + AR + AR[met-cleaved])");
}

#[test]
fn two_variable_mod_sites_quadruple_candidates() {
    // "MMK" — standard pass: single span (0,3) "MMK" with 2 M positions.
    // Standard combos: {none, M0_ox, M1_ox, both_ox} = 4.
    // Met-cleavage pass (sub_seq="MK"): single span "MK" (abs_offset=1) with 1 M position.
    // Met-cleaved combos: {none, Mox} = 2.
    // Total target = 4 + 2 = 6.
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
    assert_eq!(target_count, 6, "expected 6 (MMK×4 standard + MK×2 met-cleaved)");
}

#[test]
fn max_variable_mods_caps_combinations() {
    // "MMMK" — 3 M sites. Standard pass with max_mods=1: {none, M0_ox, M1_ox, M2_ox} = 4.
    // Met-cleavage pass (sub_seq="MMK"): 2 M sites, max_mods=1: {none, M0_ox, M1_ox} = 3.
    // Total target = 4 + 3 = 7.
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
    assert_eq!(target_count, 7, "expected 7 (MMMK×4 standard + MMK×3 met-cleaved)");
}

// ─── Terminal-mod expansion tests ────────────────────────────────────────────
//
// Terminal-location semantics in expand_mod_combinations:
//   - Peptide at protein start (start_offset == 0): position 0 gets ProtNTerm variants.
//   - Peptide NOT at protein start: position 0 gets NTerm variants.
//   - Peptide at protein end (end == protein_len): last position gets ProtCTerm variants.
//   - Peptide NOT at protein end: last position gets CTerm variants.

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
/// Met-cleavage pass (sub_seq="AAAAKMAAAAAK"):
/// - "AAAAK" (sub_seq 0..5): length=5 < min=6, skipped.
/// - "MAAAAAK" (sub_seq 5..12, abs_offset=6): is_protein_n_term=false, NTerm lookup empty → 1 candidate.
///
/// Total target: 3 + 1 = 4. The ProtNTerm mod still appears exactly once (on offset-0 peptide).
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

    // Standard pass: 2 (offset-0 "MAAAAK": unmod + ProtNTerm Acetyl) + 1 (offset-6 "MAAAAAK": unmod).
    // B5 Met-cleavage pass: 1 extra "MAAAAAK" at offset-6 (no ProtNTerm mod, NTerm lookup empty).
    // Total: 4.
    assert_eq!(
        candidates.len(), 4,
        "expected 4 candidates (2 for protein-start peptide, 1+1 for offset-6 peptide), got {}",
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
/// Standard pass:
/// - "MAAAAK" (end < protein_len): CTerm Amide applies → 2 variants.
/// - "R" (end == protein_len): ProtCTerm GlyGly applies → 2 variants.
///
/// Met-cleavage pass (sub_seq="AAAAKR"):
/// - "AAAA" (abs_end=5, not protein C-term): CTerm Amide → 2 variants.
/// - "KR" (abs_end=7, protein C-term): ProtCTerm GlyGly → 2 variants.
///
/// Total: 4 + 4 = 8.
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

    // Standard pass: "MAAAAK"×2 + "R"×2 = 4.
    // B5 Met-cleavage pass (sub_seq="AAAAKR"): "AAAA"×2 + "KR"×2 = 4.
    // Total: 8.
    assert_eq!(
        candidates.len(), 8,
        "expected 8 candidates, got {}",
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
                    // Protein-C-term peptide "R" or Met-cleaved "KR": should get ProtCTerm GlyGly (+114.04).
                    assert!(
                        m.mass_delta > 0.0,
                        "protein C-term peptide got a negative delta mod ({}); expected ProtCTerm GlyGly",
                        m.mass_delta
                    );
                } else {
                    // Non-protein-C-term peptide "MAAAAK" or Met-cleaved "AAAA": should get CTerm Amide (-0.984).
                    assert!(
                        m.mass_delta < 0.0,
                        "non-protein-C-term peptide got a positive delta mod ({}); expected CTerm Amide",
                        m.mass_delta
                    );
                }
            }
        }
    }
}

// ─── N-terminal Met cleavage tests ───────────────────────────────────────────

/// Met-cleavage generates alternative protein-N-term candidates for M-leading proteins.
///
/// Protein: "MAGER" (5 residues). With NoCleavage + min=1, the standard pass
/// emits the full protein as a single peptide at offset 0 (is_protein_n_term=true).
/// The Met-cleavage pass emits sub_seq="AGER" at offset 1 (is_protein_n_term=true,
/// since it starts at sub_seq index 0).
/// Both must be present in the candidate set.
#[test]
fn met_cleavage_generates_alternative_candidates() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"MAGER".to_vec(),
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

    // Standard: "MAGER" at offset 0, length 5.
    // Met-cleaved: "AGER" at offset 1, length 4.
    assert_eq!(candidates.len(), 2, "expected 2 target candidates (standard + Met-cleaved), got {}", candidates.len());

    let has_full = candidates.iter().any(|c| c.start_offset_in_protein == 0 && c.peptide.length() == 5);
    let has_met_cleaved = candidates.iter().any(|c| c.start_offset_in_protein == 1 && c.peptide.length() == 4);

    assert!(has_full, "missing standard candidate at offset 0 (MAGER)");
    assert!(has_met_cleaved, "missing Met-cleaved candidate at offset 1 (AGER)");
}

/// Non-M first residue does not trigger Met-cleavage enumeration.
///
/// Protein: "KAGER". Standard pass emits tryptic peptides. No second pass.
#[test]
fn non_met_first_residue_does_not_trigger_cleavage() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"KAGER".to_vec(),
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

    // Only 1 candidate: full sequence "KAGER". No Met-cleaved pass since first residue != M.
    assert_eq!(target_count, 1, "expected 1 candidate for non-M protein, got {}", target_count);
}

// ─── Phase 5: num_tolerable_termini (NTT) tests ──────────────────────────────
//
// Test protein: "AAAKBBBBBBR" (length 11)
//   - Trypsin cleaves after K(pos 3) and R(pos 10).
//   - Cleavage positions: [0, 4, 11].
//   - Strict spans (ntt=2, missed=0): (0,4)="AAAK" (too short at min=6), (4,11)="BBBBBBR" → 1 span.
//     With min=4: (0,4) and (4,11) → 2 spans.
//   - Semi-specific additional spans (ntt=1) with free-C from start=0:
//       end in [4,11] not at cleavage position → ends 5,6,7,8,9,10 → "AAAAK.." lengths 5-10.
//       With min=4: ends 4..=11, non-cleavage → 4,5,6,7,8,9,10 → 7 spans. But end=4 IS cleavage → skip. end=11 IS cleavage → skip. → ends 5,6,7,8,9,10 → 6 spans.
//       Actually let's use a simpler protein for clarity.
//
// Simpler test protein: "AAAAAKAAAAR" (length 11)
//   - Trypsin cleaves after K(4) and R(10).
//   - Cleavage positions: [0, 5, 11].
//   - Strict spans (ntt=2): (0,5)="AAAAK"(5), (5,11)="AAAAR"(6) → lengths 5 and 6.
//     With min=5, max=11: both qualify → 2 spans.
//   - Semi (ntt=1): free C from start=0: ends 5..=11 not cleavage → 6,7,8,9,10 → 5 spans.
//                   free C from start=5: ends 10..=11 not cleavage → 10 → 1 span.
//                   free N for end=5: starts 0..=0 not cleavage → (none, since 0 is cleavage pos) → 0.
//                   free N for end=11: starts 0..=6 not cleavage → 1,2,3,4,6 → 5 spans.
//   Total new semi spans = 5 + 1 + 0 + 5 = 11. Total ntt=1 = 2 (strict) + 11 = 13.
//
// Use "AAAAAKAAAAR" with min=5, max=11, missed=0, no mods.

const NTT_PROTEIN: &[u8] = b"AAAAAKAAAAR";
//   Trypsin cleavage positions: [0, 6, 11] (cleavage AFTER K at idx 5 → next pos = 6;
//   cleavage AFTER R at idx 10 → next pos = 11).
//   Let me recompute: for C-term enzyme, position i is in cleavage_positions if
//   enzyme.is_cleavable_after(seq[i-1]). K is at index 5 → position 6 (since i=6, seq[5]=K).
//   R is at index 10 → position 11. Plus 0 and 11.
//   Cleavage positions: [0, 6, 11].
//   Strict (ntt=2, min=5, max=11, missed=0): spans (0,6)=len6, (6,11)=len5 → 2.
//   Free-C from tryptic starts:
//     start=0: ends in [5,11] not in {0,6,11} → 5,7,8,9,10 → 5 spans.
//     start=6: ends in [11,11] not in {0,6,11} → none (11 is cleavage) → 0 spans.
//   Free-N for tryptic ends:
//     end=6: starts in [0,1] not in {0,6,11} → 1 → 1 span.
//     end=11: starts in [0,6] not in {0,6,11} → 1,2,3,4,5 → 5 spans. But start=6 is cleavage → {0} at start: 1,2,3,4,5 → 5 spans.
//   New semi spans = 5 + 0 + 1 + 5 = 11. Total ntt=1 = 2 + 11 = 13.

fn ntt_protein_index() -> SearchIndex {
    make_index(NTT_PROTEIN)
}

fn ntt_params(ntt: u8) -> SearchParams {
    let mut p = params(5, 11, 0);
    p.num_tolerable_termini = ntt;
    p
}

/// ntt=2 emits only strict tryptic spans (baseline).
#[test]
fn ntt_2_emits_only_strict_tryptic_spans() {
    let idx = ntt_protein_index();
    let p = ntt_params(2);
    let count = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    // Cleavage positions [0,6,11], min=5, max=11, missed=0:
    // Spans: (0,6)=len6 ✓, (6,11)=len5 ✓ → 2 strict spans.
    // NTT_PROTEIN does not start with M, so no Met-cleavage pass.
    assert_eq!(count, 2, "ntt=2 should emit exactly 2 strict tryptic spans, got {count}");
}

/// ntt=1 emits strictly more candidates than ntt=2.
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
    // Expected: 2 strict + 11 semi = 13.
    assert_eq!(ntt1_count, 13, "expected 13 ntt=1 candidates, got {ntt1_count}");
}

/// ntt=1 includes spans with a tryptic N-term but non-tryptic C-term.
#[test]
fn ntt_1_includes_free_c_term_span() {
    let idx = ntt_protein_index();
    let p = ntt_params(1);
    // A span starting at a tryptic position (0 or 6) with a non-tryptic end.
    // Example: start=0, end=5 (length 5) — start IS cleavage, end 5 is NOT cleavage.
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();
    let has_free_c = candidates.iter().any(|c| {
        // start at protein offset 0 (tryptic N-term), end at non-cleavage position.
        // end = start_offset + peptide.length() = 0 + 5 = 5 (not in {0,6,11}).
        c.start_offset_in_protein == 0 && c.peptide.length() == 5
    });
    assert!(has_free_c, "ntt=1 should include (start=0, end=5): tryptic N-term, free C-term");
}

/// ntt=1 includes spans with a non-tryptic N-term but tryptic C-term.
#[test]
fn ntt_1_includes_free_n_term_span() {
    let idx = ntt_protein_index();
    let p = ntt_params(1);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();
    // span with start=1 (non-cleavage), end=6 (tryptic C-term): length=5.
    let has_free_n = candidates.iter().any(|c| {
        c.start_offset_in_protein == 1 && c.peptide.length() == 5
    });
    assert!(has_free_n, "ntt=1 should include (start=1, end=6): free N-term, tryptic C-term");
}

/// A span where BOTH ends are tryptic should appear exactly once under ntt=1
/// (not twice from the strict + semi union).
#[test]
fn ntt_1_no_dedup_for_strict_spans() {
    let idx = ntt_protein_index();
    let p = ntt_params(1);
    let candidates: Vec<_> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .collect();
    // Count candidates with start=0, length=6 (span (0,6), both ends tryptic).
    let count_strict = candidates.iter()
        .filter(|c| c.start_offset_in_protein == 0 && c.peptide.length() == 6)
        .count();
    assert_eq!(
        count_strict, 1,
        "strict span (0,6) should appear exactly once under ntt=1, got {count_strict}"
    );
}

/// ntt=0 emits all valid-length spans regardless of cleavage sites,
/// and produces strictly more candidates than ntt=1.
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
    // For "AAAAAKAAAAR" (length 11), min=5, max=11:
    // All (start, end) pairs: start in 0..=6, end in (start+5)..=(start+11).min(11).
    // start=0: ends 5,6,7,8,9,10,11 → 7
    // start=1: ends 6,7,8,9,10,11 → 6
    // start=2: ends 7,8,9,10,11 → 5
    // start=3: ends 8,9,10,11 → 4
    // start=4: ends 9,10,11 → 3
    // start=5: ends 10,11 → 2
    // start=6: ends 11 → 1
    // Total = 7+6+5+4+3+2+1 = 28
    assert_eq!(ntt0_count, 28, "ntt=0 should emit all 28 valid-length spans, got {ntt0_count}");
}

/// ntt=0 with Trypsin should produce the same candidates as Enzyme::NonSpecific
/// with ntt=2 — WHEN missed_cleavages is set high enough to allow all spans.
///
/// Note: NonSpecific with ntt=2 routes through the cleavage-position loop where
/// every position is a cleavage site, so missed_cleavages acts as a filter.
/// For the spans to match, set missed_cleavages >= max_length so all spans pass.
#[test]
fn ntt_0_trypsin_matches_nonspecific_high_missed() {
    // Use a protein with no K/R (so trypsin has only [0, n] as cleavage positions).
    // With ntt=0 + Trypsin, we emit all (start, end) pairs — no missed-cleavage filter.
    // With NonSpecific + ntt=2 + high missed_cleavages, we also emit all pairs.
    let seq = b"AAAAAAAAAAAA"; // 12 residues, no K/R
    let idx = make_index(seq);

    let mut p_ntt0 = params(3, 8, 10); // high missed
    p_ntt0.enzyme = Enzyme::Trypsin;
    p_ntt0.num_tolerable_termini = 0;

    let mut p_ns = params(3, 8, 10); // same missed budget
    p_ns.enzyme = Enzyme::NonSpecific;
    p_ns.num_tolerable_termini = 2;

    let ntt0_count = enumerate_candidates(&idx, &p_ntt0, "XXX")
        .filter(|c| !c.is_decoy)
        .count();
    let ns_count = enumerate_candidates(&idx, &p_ns, "XXX")
        .filter(|c| !c.is_decoy)
        .count();

    // Both should emit all valid-length spans (start in 0..=9, lengths 3..=8).
    // The NonSpecific path counts internal cleavage positions as missed, but with
    // high missed budget all pass. The ntt=0 path has no cleavage constraint at all.
    // For a protein with no K/R, Trypsin has cleavage positions [0, 12].
    // ntt=0 + Trypsin: all (start, end) pairs, no filter.
    // NonSpecific: every position is cleavage, missed = end - start - 1.
    //   With missed_cleavages=10 and max_length=8: max missed = 7 → all length-8 spans pass.
    // Both should yield: sum of (n - len + 1) for len in 3..=8 = 10+9+8+7+6+5 = 45.
    assert_eq!(ntt0_count, 45, "ntt=0 + Trypsin should emit 45 spans for AAAAAAAAAAAA min=3 max=8, got {ntt0_count}");
    assert_eq!(ns_count, 45, "NonSpecific + ntt=2 high missed should also emit 45 spans, got {ns_count}");
}

/// ntt field in SearchParams defaults to 2 for default_tryptic.
#[test]
fn default_ntt_is_2() {
    let p = SearchParams::default_tryptic(aa_set());
    assert_eq!(p.num_tolerable_termini, 2, "default ntt should be 2");
}

/// A single-residue M-only protein does not trigger Met-cleavage (sequence.len() == 1).
#[test]
fn met_alone_does_not_trigger_cleavage() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(), description: "".into(),
            sequence: b"M".to_vec(),
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

    // Only 1 candidate: "M" at offset 0. Met-cleavage guard `len > 1` prevents empty sub_seq.
    assert_eq!(target_count, 1, "expected 1 candidate for M-only protein, got {}", target_count);
}
