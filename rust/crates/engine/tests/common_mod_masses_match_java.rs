//! Pin ~10 commonly-used modification monoisotopic mass deltas to the
//! values used by Java MS-GF+'s default `MSGFPlus_Mods.txt` and
//! `Modification.java` factory methods. Source-of-truth values copied
//! from those files. Each value is verifiable against UniMod
//! (https://www.unimod.org).

use engine::modification::{Modification, ModLocation, ResidueSpec};

fn bit_eq(a: f64, b: f64) -> bool { a.to_bits() == b.to_bits() }

/// (mods_txt_line, expected_name, expected_mass_delta).
/// Lines are mass-based (not composition-based) since the Phase 1
/// parser only accepts numeric mass deltas. Multi-residue mods like
/// Phospho are tested with a single-residue substitute.
fn java_common_mods() -> Vec<(&'static str, &'static str, f64)> {
    vec![
        ("57.021464,C,fix,any,Carbamidomethyl",      "Carbamidomethyl",  57.021464),
        ("15.994915,M,opt,any,Oxidation",            "Oxidation",        15.994915),
        ("79.966331,S,opt,any,Phospho",              "Phospho",          79.966331),
        ("42.010565,*,opt,Prot-N-term,Acetyl",       "Acetyl",           42.010565),
        ("229.162932,K,fix,any,TMT6plex",            "TMT6plex",         229.162932),
        ("229.162932,*,fix,N-term,TMT6plex",         "TMT6plex",         229.162932),
        ("144.102063,K,fix,any,iTRAQ4plex",          "iTRAQ4plex",       144.102063),
        ("304.205360,K,fix,any,iTRAQ8plex",          "iTRAQ8plex",       304.205360),
        ("14.015650,K,opt,any,Methyl",               "Methyl",           14.015650),
        ("28.031300,K,opt,any,Dimethyl",             "Dimethyl",         28.031300),
        ("42.046950,K,opt,any,Trimethyl",            "Trimethyl",        42.046950),
    ]
}

#[test]
fn parses_to_expected_name_and_mass() {
    for (line, expected_name, expected_mass) in java_common_mods() {
        let m = Modification::from_mods_txt_line(line)
            .unwrap_or_else(|e| panic!("parse failed for {line:?}: {e:?}"));
        assert_eq!(m.name, expected_name, "name drift on {line:?}");
        assert!(
            bit_eq(m.mass_delta, expected_mass),
            "mass drift on {:?}: rust={}, expected={}",
            line, m.mass_delta, expected_mass
        );
    }
}

#[test]
fn nterm_tmt_uses_wildcard_residue() {
    let m = Modification::from_mods_txt_line("229.162932,*,fix,N-term,TMT6plex").unwrap();
    assert_eq!(m.residue, ResidueSpec::Wildcard);
    assert_eq!(m.location, ModLocation::NTerm);
    assert!(m.fixed);
}

#[test]
fn prot_nterm_acetyl_is_variable_wildcard() {
    let m = Modification::from_mods_txt_line("42.010565,*,opt,Prot-N-term,Acetyl").unwrap();
    assert_eq!(m.residue, ResidueSpec::Wildcard);
    assert_eq!(m.location, ModLocation::ProtNTerm);
    assert!(!m.fixed);
}
