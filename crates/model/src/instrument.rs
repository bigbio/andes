#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InstrumentType {
    LowRes,
    HighRes,
    TOF,
    QExactive,
    OrbitrapAstral,
    TimsTOF,
}

impl InstrumentType {
    pub fn name(self) -> &'static str {
        match self {
            InstrumentType::LowRes => "LowRes",
            InstrumentType::HighRes => "HighRes",
            InstrumentType::TOF => "TOF",
            InstrumentType::QExactive => "QExactive",
            InstrumentType::OrbitrapAstral => "OrbitrapAstral",
            InstrumentType::TimsTOF => "TimsTOF",
        }
    }

    pub fn is_high_resolution(self) -> bool {
        matches!(
            self,
            InstrumentType::HighRes
                | InstrumentType::TOF
                | InstrumentType::QExactive
                | InstrumentType::OrbitrapAstral
                | InstrumentType::TimsTOF
        )
    }

    pub fn from_name(s: &str) -> Option<Self> {
        match s {
            "LowRes" => Some(InstrumentType::LowRes),
            "HighRes" => Some(InstrumentType::HighRes),
            "TOF" => Some(InstrumentType::TOF),
            "QExactive" => Some(InstrumentType::QExactive),
            "OrbitrapAstral" => Some(InstrumentType::OrbitrapAstral),
            "TimsTOF" => Some(InstrumentType::TimsTOF),
            _ => None,
        }
    }

    pub fn family_fallback(self) -> Self {
        match self {
            InstrumentType::OrbitrapAstral => InstrumentType::QExactive,
            InstrumentType::TimsTOF => InstrumentType::TOF,
            other => other,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_round_trips() {
        for i in [
            InstrumentType::LowRes,
            InstrumentType::HighRes,
            InstrumentType::TOF,
            InstrumentType::QExactive,
            InstrumentType::OrbitrapAstral,
            InstrumentType::TimsTOF,
        ] {
            assert_eq!(InstrumentType::from_name(i.name()), Some(i));
        }
    }

    #[test]
    fn from_name_known_variants() {
        assert_eq!(
            InstrumentType::from_name("LowRes"),
            Some(InstrumentType::LowRes)
        );
        assert_eq!(
            InstrumentType::from_name("HighRes"),
            Some(InstrumentType::HighRes)
        );
        assert_eq!(InstrumentType::from_name("TOF"), Some(InstrumentType::TOF));
        assert_eq!(
            InstrumentType::from_name("QExactive"),
            Some(InstrumentType::QExactive)
        );
        assert_eq!(
            InstrumentType::from_name("OrbitrapAstral"),
            Some(InstrumentType::OrbitrapAstral)
        );
        assert_eq!(
            InstrumentType::from_name("TimsTOF"),
            Some(InstrumentType::TimsTOF)
        );
    }

    #[test]
    fn from_name_case_sensitive() {
        assert_eq!(InstrumentType::from_name("lowres"), None);
        assert_eq!(InstrumentType::from_name("tof"), None);
        assert_eq!(InstrumentType::from_name("orbitrapastral"), None);
        assert_eq!(InstrumentType::from_name("timstof"), None);
    }

    #[test]
    fn from_name_unknown() {
        assert_eq!(InstrumentType::from_name("Astral"), None);
        assert_eq!(InstrumentType::from_name(""), None);
    }

    #[test]
    fn family_fallback_mapping() {
        assert_eq!(
            InstrumentType::OrbitrapAstral.family_fallback(),
            InstrumentType::QExactive
        );
        assert_eq!(
            InstrumentType::TimsTOF.family_fallback(),
            InstrumentType::TOF
        );
        assert_eq!(
            InstrumentType::QExactive.family_fallback(),
            InstrumentType::QExactive
        );
        assert_eq!(InstrumentType::TOF.family_fallback(), InstrumentType::TOF);
        assert_eq!(
            InstrumentType::LowRes.family_fallback(),
            InstrumentType::LowRes
        );
        assert_eq!(
            InstrumentType::HighRes.family_fallback(),
            InstrumentType::HighRes
        );
    }

    #[test]
    fn is_high_resolution_includes_new_variants() {
        assert!(InstrumentType::OrbitrapAstral.is_high_resolution());
        assert!(InstrumentType::TimsTOF.is_high_resolution());
    }
}
