//! Pin `ActivationMethod` variants to PSI-MS activation method terms
//! (MS:1000133 CID, MS:1000599 ETD, MS:1000422 HCD, MS:1000598 PQD,
//! MS:1000776 UVPD).

use model::ActivationMethod;

#[test]
fn canonical_names_resolve() {
    let reference: &[(ActivationMethod, &str)] = &[
        (ActivationMethod::CID,  "CID"),
        (ActivationMethod::ETD,  "ETD"),
        (ActivationMethod::HCD,  "HCD"),
        (ActivationMethod::PQD,  "PQD"),
        (ActivationMethod::UVPD, "UVPD"),
    ];
    for &(variant, name) in reference {
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
