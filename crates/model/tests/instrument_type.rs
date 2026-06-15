//! Pin `InstrumentType` variant names to the andes instrument-class enum.

use model::InstrumentType;

#[test]
fn canonical_names_resolve() {
    let reference: &[(InstrumentType, &str)] = &[
        (InstrumentType::LowRes,    "LowRes"),
        (InstrumentType::HighRes,   "HighRes"),
        (InstrumentType::TOF,       "TOF"),
        (InstrumentType::QExactive, "QExactive"),
    ];
    for &(variant, name) in reference {
        assert_eq!(variant.name(), name);
        assert_eq!(InstrumentType::from_name(name), Some(variant));
    }
}
