//! Mass spectrometer instrument categories.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InstrumentType {
    LowRes,
    HighRes,
    TOF,
    QExactive,
}

impl InstrumentType {
    pub fn name(self) -> &'static str {
        match self {
            InstrumentType::LowRes    => "LowRes",
            InstrumentType::HighRes   => "HighRes",
            InstrumentType::TOF       => "TOF",
            InstrumentType::QExactive => "QExactive",
        }
    }

    /// Whether the instrument produces high-resolution MS/MS spectra.
    ///
    /// Mirrors Java's `InstrumentType.isHighResolution()`: HighRes,
    /// TOF, and QExactive return `true`; LowRes returns `false`. Used by
    /// `compute_psm_features` to mirror Java's `PSMFeatureFinder` hardcoded
    /// 20 ppm (high-res) / 0.5 Da (low-res) fragment tolerance for
    /// feature counting, independent of `param.mme` (which the rank-based
    /// scoring tables use at a coarser resolution for binning).
    pub fn is_high_resolution(self) -> bool {
        matches!(
            self,
            InstrumentType::HighRes | InstrumentType::TOF | InstrumentType::QExactive
        )
    }

    /// Case-sensitive lookup.
    pub fn from_name(s: &str) -> Option<Self> {
        match s {
            "LowRes"    => Some(InstrumentType::LowRes),
            "HighRes"   => Some(InstrumentType::HighRes),
            "TOF"       => Some(InstrumentType::TOF),
            "QExactive" => Some(InstrumentType::QExactive),
            _           => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_round_trips() {
        for i in [
            InstrumentType::LowRes, InstrumentType::HighRes,
            InstrumentType::TOF,    InstrumentType::QExactive,
        ] {
            assert_eq!(InstrumentType::from_name(i.name()), Some(i));
        }
    }

    #[test]
    fn from_name_known_variants() {
        assert_eq!(InstrumentType::from_name("LowRes"),    Some(InstrumentType::LowRes));
        assert_eq!(InstrumentType::from_name("HighRes"),   Some(InstrumentType::HighRes));
        assert_eq!(InstrumentType::from_name("TOF"),       Some(InstrumentType::TOF));
        assert_eq!(InstrumentType::from_name("QExactive"), Some(InstrumentType::QExactive));
    }

    #[test]
    fn from_name_case_sensitive() {
        assert_eq!(InstrumentType::from_name("lowres"), None);
        assert_eq!(InstrumentType::from_name("tof"), None);
    }

    #[test]
    fn from_name_unknown() {
        assert_eq!(InstrumentType::from_name("Astral"), None);
        assert_eq!(InstrumentType::from_name(""), None);
    }
}
