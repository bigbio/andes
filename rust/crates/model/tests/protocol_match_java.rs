//! Pin `Protocol` variants to Java `edu.ucsd.msjava.msutil.Protocol`
//! (lines 56-61).

use model::Protocol;

#[test]
fn java_canonical_names_resolve() {
    let java: &[(Protocol, &str)] = &[
        (Protocol::Automatic,       "Automatic"),
        (Protocol::Phosphorylation, "Phosphorylation"),
        (Protocol::ITRAQ,           "iTRAQ"),
        (Protocol::ITRAQPhospho,    "iTRAQPhospho"),
        (Protocol::TMT,             "TMT"),
        (Protocol::Standard,        "Standard"),
    ];
    for &(variant, name) in java {
        assert_eq!(variant.name(), name);
        assert_eq!(Protocol::from_name(name), Some(variant));
    }
}
