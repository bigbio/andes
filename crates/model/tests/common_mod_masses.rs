//! Pin commonly-used modification monoisotopic mass deltas, as the bundled
//! UNIMOD table gives them to a modification list, to UniMod
//! (https://www.unimod.org) reference values.

use model::mod_spec::parse_mod_list;
use model::modification::{ModLocation, ResidueSpec};

fn bit_eq(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

/// (list item, expected_name, expected_mass_delta).
fn unimod_reference_mods() -> Vec<(&'static str, &'static str, f64)> {
    vec![
        ("Carbamidomethyl (C)", "Carbamidomethyl", 57.021464),
        ("Oxidation (M)", "Oxidation", 15.994915),
        ("Phospho (S)", "Phospho", 79.966331),
        ("Acetyl (Protein N-term)", "Acetyl", 42.010565),
        ("TMT6plex (K)", "TMT6plex", 229.162932),
        ("TMT6plex (N-term)", "TMT6plex", 229.162932),
        ("iTRAQ4plex (K)", "iTRAQ4plex", 144.102063),
        ("iTRAQ8plex (K)", "iTRAQ8plex", 304.205360),
        ("Methyl (K)", "Methyl", 14.015650),
        ("Dimethyl (K)", "Dimethyl", 28.031300),
        ("Trimethyl (K)", "Trimethyl", 42.046950),
    ]
}

#[test]
fn resolves_to_expected_name_and_mass() {
    for (item, expected_name, expected_mass) in unimod_reference_mods() {
        let m = parse_mod_list(item, false)
            .unwrap_or_else(|e| panic!("parse failed for {item:?}: {e}"))
            .remove(0);
        assert_eq!(m.name, expected_name, "name drift on {item:?}");
        assert!(
            bit_eq(m.mass_delta, expected_mass),
            "mass drift on {:?}: rust={}, expected={}",
            item,
            m.mass_delta,
            expected_mass
        );
    }
}

#[test]
fn nterm_tmt_uses_wildcard_residue() {
    let m = parse_mod_list("TMT6plex (N-term)", true).unwrap().remove(0);
    assert_eq!(m.residue, ResidueSpec::Wildcard);
    assert_eq!(m.location, ModLocation::NTerm);
    assert!(m.fixed);
}

#[test]
fn prot_nterm_acetyl_is_variable_wildcard() {
    let m = parse_mod_list("Acetyl (Protein N-term)", false)
        .unwrap()
        .remove(0);
    assert_eq!(m.residue, ResidueSpec::Wildcard);
    assert_eq!(m.location, ModLocation::ProtNTerm);
    assert!(!m.fixed);
}
