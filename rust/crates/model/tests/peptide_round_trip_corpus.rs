//! Phase 1 Display ↔ from_str round-trip stress test. Validates that for
//! every constructible `Peptide` in our representative corpus,
//! `Peptide::from_str(&p.to_string(), &aa_set) == Ok(p)` (structural).
//!
//! This is NOT a Java byte-parity gate — Java's `Peptide.toString()`
//! uses a different format (3-decimal, fixed mods invisible, no
//! flanking). Phase 7's output crate owns PIN/TSV byte parity.

use model::{
    AminoAcid, AminoAcidSet, AminoAcidSetBuilder, ModLocation, Modification,
    Peptide, ResidueSpec,
};

fn corpus_aa_set() -> AminoAcidSet {
    let cam = Modification {
        name: "Carbamidomethyl".to_string(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
    };
    let ox = Modification {
        name: "Oxidation".to_string(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
    };
    let pyro_glu = Modification {
        name: "Pyro-glu".to_string(),
        mass_delta: -17.02655,
        residue: ResidueSpec::Specific(b'Q'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
    };
    AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .add_variable_mod(pyro_glu)
        .build()
        .unwrap()
}

/// Build a peptide from a sequence with optional `(index, mod_name)` annotations.
fn build_peptide(
    seq: &[u8],
    pre: u8,
    post: u8,
    mods: &[(usize, &str)],
    aa_set: &AminoAcidSet,
) -> Peptide {
    let mut residues: Vec<AminoAcid> = seq.iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    for &(idx, mod_name) in mods {
        let r = seq[idx];
        let variant = aa_set
            .variants_for(r, ModLocation::Anywhere)
            .iter()
            .find(|aa| aa.mod_.as_ref().map(|m| m.name == mod_name).unwrap_or(false))
            .cloned()
            .unwrap_or_else(|| panic!("mod {mod_name:?} not found for residue {}", r as char));
        residues[idx] = variant;
    }
    Peptide::new(residues, pre, post)
}

#[test]
fn round_trip_unmodified_corpus() {
    let aa_set = corpus_aa_set();
    let cases: &[(&[u8], u8, u8)] = &[
        (b"PEPTIDE",      b'_', b'-'),
        (b"PEPTIDE",      b'K', b'R'),
        (b"GAVL",         b'_', b'A'),
        (b"AAAAA",        b'A', b'A'),
        (b"WYRFLMHK",     b'R', b'P'),
        (b"GG",           b'_', b'-'),  // shortest realistic
        (b"M",            b'_', b'-'),  // single residue
    ];
    for &(seq, pre, post) in cases {
        let p = build_peptide(seq, pre, post, &[], &aa_set);
        let serialized = p.to_string();
        let parsed = Peptide::from_str(&serialized, &aa_set)
            .unwrap_or_else(|e| panic!("from_str failed on {serialized:?}: {e}"));
        assert_eq!(parsed.to_string(), serialized,
            "Display→from_str→Display drift on {serialized:?}");
        assert_eq!(parsed, p,
            "Structural mismatch on {serialized:?}");
    }
}

#[test]
fn round_trip_with_carbamidomethyl_c() {
    let aa_set = corpus_aa_set();
    // Carbamidomethyl is a FIXED mod — every C in the AA set is already
    // modified. The build_peptide helper picks up that variant.
    let p = build_peptide(b"PEC", b'K', b'R', &[(2, "Carbamidomethyl")], &aa_set);
    let serialized = p.to_string();
    assert_eq!(serialized, "K.PEC+57.02146.R");
    let parsed = Peptide::from_str(&serialized, &aa_set).unwrap();
    assert_eq!(parsed, p);
}

#[test]
fn round_trip_with_oxidation_m() {
    let aa_set = corpus_aa_set();
    let p = build_peptide(b"MEMD", b'_', b'-', &[(2, "Oxidation")], &aa_set);
    let serialized = p.to_string();
    assert_eq!(serialized, "_.MEM+15.99491D.-");
    let parsed = Peptide::from_str(&serialized, &aa_set).unwrap();
    assert_eq!(parsed, p);
}

#[test]
fn round_trip_with_negative_mass_mod() {
    let aa_set = corpus_aa_set();
    let p = build_peptide(b"QPEPT", b'_', b'-', &[(0, "Pyro-glu")], &aa_set);
    let serialized = p.to_string();
    assert_eq!(serialized, "_.Q-17.02655PEPT.-");
    let parsed = Peptide::from_str(&serialized, &aa_set).unwrap();
    assert_eq!(parsed, p);
}

#[test]
fn round_trip_with_multi_mod() {
    let aa_set = corpus_aa_set();
    let p = build_peptide(b"MCEM", b'K', b'R',
        &[(1, "Carbamidomethyl"), (3, "Oxidation")], &aa_set);
    let serialized = p.to_string();
    let parsed = Peptide::from_str(&serialized, &aa_set).unwrap();
    assert_eq!(parsed, p);
    assert_eq!(parsed.to_string(), serialized);
}

#[test]
fn from_str_then_display_is_identity() {
    let aa_set = corpus_aa_set();
    let inputs = [
        "_.PEPTIDE.-",
        "K.PEPTIDE.R",
        "_.M.-",
        "_.MEM+15.99491DE.-",
        "K.PEC+57.02146PM+15.99491DE.R",
        "_.Q-17.02655PEPT.-",
    ];
    for s in &inputs {
        let p = Peptide::from_str(s, &aa_set)
            .unwrap_or_else(|e| panic!("from_str failed on {s:?}: {e}"));
        assert_eq!(p.to_string(), *s, "from_str→Display drift on {s:?}");
    }
}
