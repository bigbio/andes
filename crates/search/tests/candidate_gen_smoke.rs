//! Handcrafted candidate-enumeration tests.

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
    assert_eq!(
        target_candidates.len(),
        2,
        "expected 2 target candidates (standard + Met-cleaved), got {}",
        target_candidates.len()
    );
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
    let full = target_candidates
        .iter()
        .find(|c| c.start_offset_in_protein == 0)
        .unwrap();
    assert_eq!(full.peptide.length(), 11);
    // Met-cleaved target: sequence[1..] at offset 1, length 10.
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
            accession: "P1".into(),
            description: "".into(),
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
    // "AKMKLR" — Trypsin cleaves after positions 1 (K), 3 (K), 5 (R).
    // Cleavage positions: [0, 2, 4, 6].
    // missed=0, length 2-6: spans (0,2)="AK", (2,4)="MK", (4,6)="LR" — 3 spans.
    // (Note: 'B' is not standard so we use 'L' which IS standard.)
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
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
            accession: "P1".into(),
            description: "".into(),
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
    assert_eq!(
        target_count, 4,
        "expected 4 target candidates (MK + MKox + AR + AR[met-cleaved])"
    );
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
            accession: "P1".into(),
            description: "".into(),
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
    assert_eq!(
        target_count, 6,
        "expected 6 (MMK×4 standard + MK×2 met-cleaved)"
    );
}

#[test]
fn max_variable_mods_caps_combinations() {
    // "MMMK" — 3 M sites. Standard pass with max_mods=1: {none, M0_ox, M1_ox, M2_ox} = 4.
    // Met-cleavage pass (sub_seq="MMK"): 2 M sites, max_mods=1: {none, M0_ox, M1_ox} = 3.
    // Total target = 4 + 3 = 7.
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
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
    assert_eq!(
        target_count, 7,
        "expected 7 (MMMK×4 standard + MMK×3 met-cleaved)"
    );
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
        neutral_losses: Vec::new(),
        loss_class: 0,
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
        neutral_losses: Vec::new(),
        loss_class: 0,
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
        neutral_losses: Vec::new(),
        loss_class: 0,
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
            accession: "P1".into(),
            description: "".into(),
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
        candidates.len(),
        4,
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
    let mod_count = candidates
        .iter()
        .filter(|c| c.peptide.residues[0].is_modified())
        .count();
    assert_eq!(
        mod_count, 1,
        "exactly 1 candidate should have the ProtNTerm mod"
    );
}

/// A peptide-N-term (NTerm) mod applies to EVERY peptide N-terminus, INCLUDING
/// the protein-start peptide — the protein's N-terminus IS a peptide N-terminus
/// (finding 2.3: the old ProtNTerm-XOR-NTerm lookup wrongly excluded it).
///
/// Protein: "AAAAAAKMAAAAAK" (length 14).
/// Trypsin + missed=0 → (0..7)="AAAAAAK" (protein N-term) + (7..14)="MAAAAAK".
/// With an NTerm Acetyl variable mod and max_mods=1, BOTH peptides get
/// unmod + NTerm-acetyl = 2 each → total 4.
#[test]
fn nterm_mod_applies_to_every_peptide_n_terminus() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
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

    // Both "AAAAAAK" (offset 0) and "MAAAAAK" (offset 7): unmod + NTerm-acetyl = 2 each.
    // Total: 4.
    assert_eq!(
        candidates.len(), 4,
        "expected 4 candidates (unmod + NTerm-acetyl for both the offset-0 and offset-7 peptides), got {}",
        candidates.len()
    );

    // The NTerm mod now appears on BOTH peptide N-termini (offsets 0 and 7).
    let modified_offsets: std::collections::BTreeSet<usize> = candidates
        .iter()
        .filter(|c| c.peptide.residues[0].is_modified())
        .map(|c| c.start_offset_in_protein)
        .collect();
    assert_eq!(
        modified_offsets, std::collections::BTreeSet::from([0usize, 7usize]),
        "NTerm acetyl should appear on both peptide N-termini (offsets 0 and 7), got {modified_offsets:?}"
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
/// A protein-C-term peptide is ALSO a peptide C-terminus, so it carries BOTH the
/// peptide-C-term (CTerm) and protein-C-term (ProtCTerm) mods (finding 2.3).
/// Standard pass:
/// - "MAAAAK" (end < protein_len): CTerm Amide → unmod + Amide = 2.
/// - "R" (end == protein_len): CTerm Amide ∪ ProtCTerm GlyGly → unmod + Amide + GlyGly = 3.
///
/// Met-cleavage pass (sub_seq="AAAAKR"):
/// - "AAAA" (abs_end=5, not protein C-term): CTerm Amide → 2.
/// - "KR" (abs_end=7, protein C-term): CTerm Amide ∪ ProtCTerm GlyGly → 3.
///
/// Total: (2 + 3) + (2 + 3) = 10.
///
/// The CTerm mod does NOT bleed into NON-protein-C-term peptides (those get only Amide).
#[test]
fn c_term_and_protein_c_term_distinguished() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
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

    // Standard pass: "MAAAAK"×2 + "R"×3 = 5.
    // B5 Met-cleavage pass (sub_seq="AAAAKR"): "AAAA"×2 + "KR"×3 = 5.
    // Total: 10 (protein-C-term peptides "R"/"KR" each carry unmod + CTerm Amide + ProtCTerm GlyGly).
    assert_eq!(
        candidates.len(),
        10,
        "expected 10 candidates, got {}",
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
                    // Protein-C-term peptide "R"/"KR": it is BOTH a peptide and a
                    // protein C-terminus, so either the CTerm Amide (-0.984016) OR
                    // the ProtCTerm GlyGly (+114.042927) is valid here (finding 2.3).
                    assert!(
                        (m.mass_delta - (-0.984016)).abs() < 1e-4
                            || (m.mass_delta - 114.042927).abs() < 1e-4,
                        "protein C-term peptide got an unexpected mod delta ({}); expected CTerm Amide (-0.984016) or ProtCTerm GlyGly (114.042927)",
                        m.mass_delta
                    );
                } else {
                    // Non-protein-C-term peptide "MAAAAK" or Met-cleaved "AAAA": only CTerm Amide (-0.984016).
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
            accession: "P1".into(),
            description: "".into(),
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

/// Non-M first residue does not trigger Met-cleavage enumeration.
///
/// Protein: "KAGER". Standard pass emits tryptic peptides. No second pass.
#[test]
fn non_met_first_residue_does_not_trigger_cleavage() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
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
    assert_eq!(
        target_count, 1,
        "expected 1 candidate for non-M protein, got {}",
        target_count
    );
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
//   (Semi spans are additionally bounded by max_missed_cleavages — see the
//   per-span derivation below the NTT_PROTEIN const for the missed=0 count.)
//
// Use "AAAAAKAAAAR" with min=5, max=11, missed=0, no mods.

const NTT_PROTEIN: &[u8] = b"AAAAAKAAAAR";
//   Trypsin cleavage positions: [0, 6, 11] (cleavage AFTER K at idx 5 → pos 6;
//   AFTER R at idx 10 → pos 11; plus the protein ends 0 and 11).
//   Strict (ntt=2, min=5, max=11, missed=0): (0,6)=len6, (6,11)=len5 → 2 spans.
//   Semi-tryptic spans (ntt=1) ALSO respect max_missed_cleavages (=0 here), so any
//   span whose [start,end] straddles the internal cleavage site (pos 6) is pruned —
//   the same internal-site bound the strict path applies (finding 2.4):
//     Free-C from start=0: non-cleavage ends 5,7,8,9,10; only end=5 has 0 internal
//       sites (5<6) → 1 span (7,8,9,10 each cross site 6 → pruned).
//     Free-C from start=6: no non-cleavage end in [11,11] → 0.
//     Free-N for end=6:  start=1 (non-cleavage), no site in (1,6) → 1 span.
//     Free-N for end=11: starts 1..5 all cross site 6 → pruned → 0.
//   New semi spans = 1 + 0 + 1 + 0 = 2. Total ntt=1 = 2 strict + 2 semi = 4.

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
    assert_eq!(
        count, 2,
        "ntt=2 should emit exactly 2 strict tryptic spans, got {count}"
    );
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
    // Expected: 2 strict + 2 semi = 4 (semi spans crossing the internal K/R site
    // are pruned by the max_missed_cleavages=0 bound — finding 2.4).
    assert_eq!(
        ntt1_count, 4,
        "expected 4 ntt=1 candidates, got {ntt1_count}"
    );
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
    assert!(
        has_free_c,
        "ntt=1 should include (start=0, end=5): tryptic N-term, free C-term"
    );
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
    let has_free_n = candidates
        .iter()
        .any(|c| c.start_offset_in_protein == 1 && c.peptide.length() == 5);
    assert!(
        has_free_n,
        "ntt=1 should include (start=1, end=6): free N-term, tryptic C-term"
    );
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
    let count_strict = candidates
        .iter()
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
    assert_eq!(
        ntt0_count, 28,
        "ntt=0 should emit all 28 valid-length spans, got {ntt0_count}"
    );
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
    assert_eq!(
        ntt0_count, 45,
        "ntt=0 + Trypsin should emit 45 spans for AAAAAAAAAAAA min=3 max=8, got {ntt0_count}"
    );
    assert_eq!(
        ns_count, 45,
        "NonSpecific + ntt=2 high missed should also emit 45 spans, got {ns_count}"
    );
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
            accession: "P1".into(),
            description: "".into(),
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
