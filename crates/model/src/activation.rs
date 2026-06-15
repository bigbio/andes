//! Activation methods used by tandem MS spectrum acquisition. The five
//! canonical variants (CID/ETD/HCD/PQD/UVPD) are pinned by
//! `tests/activation_method_match_java.rs`.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ActivationMethod {
    CID,
    ETD,
    HCD,
    PQD,
    UVPD,
}

impl ActivationMethod {
    pub fn name(self) -> &'static str {
        match self {
            ActivationMethod::CID  => "CID",
            ActivationMethod::ETD  => "ETD",
            ActivationMethod::HCD  => "HCD",
            ActivationMethod::PQD  => "PQD",
            ActivationMethod::UVPD => "UVPD",
        }
    }

    /// Whether fragment prediction should emit neutral-loss-shifted ions for
    /// this activation method.
    ///
    /// Collisional / vibrational and photon-based activation (CID, HCD, PQD,
    /// UVPD) cleave labile groups, so neutral-loss fragment ions are expected
    /// and should be predicted. Electron-based dissociation (ETD, and any
    /// electron-transfer/-capture variant) preserves labile modifications on
    /// the backbone fragments, so loss ions are NOT predicted.
    pub fn predicts_neutral_losses(self) -> bool {
        match self {
            // Electron-based: labile mods survive → no loss ions.
            ActivationMethod::ETD => false,
            // Collisional / photon-based: labile mods fragment → loss ions.
            ActivationMethod::CID
            | ActivationMethod::HCD
            | ActivationMethod::PQD
            | ActivationMethod::UVPD => true,
        }
    }

    /// Case-sensitive lookup. Returns `None` for unknown names, including the
    /// runtime sentinels `ASWRITTEN` and `FUSION` which never appear in
    /// stored `.param` files.
    pub fn from_name(s: &str) -> Option<Self> {
        match s {
            "CID"  => Some(ActivationMethod::CID),
            "ETD"  => Some(ActivationMethod::ETD),
            "HCD"  => Some(ActivationMethod::HCD),
            "PQD"  => Some(ActivationMethod::PQD),
            "UVPD" => Some(ActivationMethod::UVPD),
            _      => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_round_trips() {
        for m in [
            ActivationMethod::CID, ActivationMethod::ETD,
            ActivationMethod::HCD, ActivationMethod::PQD,
            ActivationMethod::UVPD,
        ] {
            assert_eq!(ActivationMethod::from_name(m.name()), Some(m));
        }
    }

    #[test]
    fn from_name_known_variants() {
        assert_eq!(ActivationMethod::from_name("CID"),  Some(ActivationMethod::CID));
        assert_eq!(ActivationMethod::from_name("ETD"),  Some(ActivationMethod::ETD));
        assert_eq!(ActivationMethod::from_name("HCD"),  Some(ActivationMethod::HCD));
        assert_eq!(ActivationMethod::from_name("PQD"),  Some(ActivationMethod::PQD));
        assert_eq!(ActivationMethod::from_name("UVPD"), Some(ActivationMethod::UVPD));
    }

    #[test]
    fn from_name_case_sensitive() {
        assert_eq!(ActivationMethod::from_name("cid"), None);
        assert_eq!(ActivationMethod::from_name("hcd"), None);
    }

    #[test]
    fn from_name_runtime_sentinels_unknown() {
        // ASWRITTEN and FUSION are runtime metadata strings that should
        // never appear in stored .param files; we omit them and return
        // None so the param loader can surface BadEnum.
        assert_eq!(ActivationMethod::from_name("As written in the spectrum or CID if no info"), None);
        assert_eq!(ActivationMethod::from_name("Merge spectra from the same precursor"), None);
    }

    #[test]
    fn from_name_unknown() {
        assert_eq!(ActivationMethod::from_name("garbage"), None);
        assert_eq!(ActivationMethod::from_name(""), None);
    }

    #[test]
    fn etd_does_not_predict_neutral_losses() {
        // Electron-based dissociation preserves labile mods on fragments.
        assert!(!ActivationMethod::ETD.predicts_neutral_losses());
    }

    #[test]
    fn collisional_and_photon_methods_predict_neutral_losses() {
        assert!(ActivationMethod::CID.predicts_neutral_losses());
        assert!(ActivationMethod::HCD.predicts_neutral_losses());
        assert!(ActivationMethod::PQD.predicts_neutral_losses());
        assert!(ActivationMethod::UVPD.predicts_neutral_losses());
    }
}
