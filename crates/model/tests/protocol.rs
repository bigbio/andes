//! Pin `Protocol` variant names to the andes search-protocol enum.

use model::Protocol;

#[test]
fn canonical_names_resolve() {
    let reference: &[(Protocol, &str)] = &[
        (Protocol::Automatic,       "Automatic"),
        (Protocol::Phosphorylation, "Phosphorylation"),
        (Protocol::ITRAQ,           "iTRAQ"),
        (Protocol::ITRAQPhospho,    "iTRAQPhospho"),
        (Protocol::TMT,             "TMT"),
        (Protocol::Standard,        "Standard"),
    ];
    for &(variant, name) in reference {
        assert_eq!(variant.name(), name);
        assert_eq!(Protocol::from_name(name), Some(variant));
    }
}
