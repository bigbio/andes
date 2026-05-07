//! Pin `InstrumentType` variants to Java
//! `edu.ucsd.msjava.msutil.InstrumentType` (lines 73-76).

use model::InstrumentType;

#[test]
fn java_canonical_names_resolve() {
    let java: &[(InstrumentType, &str)] = &[
        (InstrumentType::LowRes,    "LowRes"),
        (InstrumentType::HighRes,   "HighRes"),
        (InstrumentType::TOF,       "TOF"),
        (InstrumentType::QExactive, "QExactive"),
    ];
    for &(variant, name) in java {
        assert_eq!(variant.name(), name);
        assert_eq!(InstrumentType::from_name(name), Some(variant));
    }
}
