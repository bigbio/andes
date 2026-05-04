//! Pin `ActivationMethod` variants to Java
//! `edu.ucsd.msjava.msutil.ActivationMethod` (lines 125-131).
//! Source-of-truth strings copied by hand.

use engine::ActivationMethod;

#[test]
fn java_canonical_names_resolve() {
    let java: &[(ActivationMethod, &str)] = &[
        (ActivationMethod::CID,  "CID"),
        (ActivationMethod::ETD,  "ETD"),
        (ActivationMethod::HCD,  "HCD"),
        (ActivationMethod::PQD,  "PQD"),
        (ActivationMethod::UVPD, "UVPD"),
    ];
    for &(variant, name) in java {
        assert_eq!(variant.name(), name);
        assert_eq!(ActivationMethod::from_name(name), Some(variant));
    }
}

#[test]
fn no_extra_variants() {
    let names: Vec<_> = [
        ActivationMethod::CID,  ActivationMethod::ETD,
        ActivationMethod::HCD,  ActivationMethod::PQD,
        ActivationMethod::UVPD,
    ].iter().map(|m| m.name()).collect();
    let mut sorted = names.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(names.len(), sorted.len(), "duplicate name(s) in ActivationMethod");
}
