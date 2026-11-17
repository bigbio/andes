#![allow(dead_code)]

use std::path::PathBuf;

use model::{AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec};

pub fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

pub fn aa_set() -> model::AminoAcidSet {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
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
    AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

pub fn strip_flanking_and_mods(pin_pep: &str) -> String {
    let bytes = pin_pep.as_bytes();
    if bytes.len() < 5 {
        return String::new();
    }
    if bytes[1] != b'.' || bytes[bytes.len() - 2] != b'.' {
        return String::new();
    }
    let middle = &pin_pep[2..pin_pep.len() - 2];
    let mut out = String::with_capacity(middle.len());
    let mut chars = middle.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '+' || c == '-' {
            while let Some(&nc) = chars.peek() {
                if nc.is_ascii_digit() || nc == '.' {
                    chars.next();
                } else {
                    break;
                }
            }
        } else if c.is_ascii_uppercase() {
            out.push(c);
        }
    }
    out
}

#[cfg(test)]
mod parser_tests {
    use super::strip_flanking_and_mods;

    #[test]
    fn strips_flanking_only() {
        assert_eq!(strip_flanking_and_mods("R.PEPTIDE.K"), "PEPTIDE");
    }

    #[test]
    fn strips_one_mod_mass() {
        assert_eq!(strip_flanking_and_mods("K.PEPTM+15.995DE.R"), "PEPTMDE");
    }

    #[test]
    fn strips_multiple_mod_masses() {
        assert_eq!(
            strip_flanking_and_mods("K.GAC+57.021LLPKIETM+15.995R.E"),
            "GACLLPKIETMR"
        );
    }

    #[test]
    fn strips_negative_mod_mass() {
        assert_eq!(strip_flanking_and_mods("K.PEPM-18.0R.E"), "PEPMR");
    }

    #[test]
    fn handles_protein_terminal_dash_flanking() {
        assert_eq!(strip_flanking_and_mods("-.PEPTIDE.R"), "PEPTIDE");
        assert_eq!(strip_flanking_and_mods("R.PEPTIDE.-"), "PEPTIDE");
    }
}
