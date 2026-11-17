use model::{ActivationMethod, InstrumentType};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ResolutionSource {
    Metadata,
    FragmentTolFlag,
    Peaks,
    Default,
}

impl ResolutionSource {
    pub(crate) fn tag(self) -> &'static str {
        match self {
            ResolutionSource::Metadata => "detected",
            ResolutionSource::FragmentTolFlag => "from --fragment-tol-*",
            ResolutionSource::Peaks => "inferred from peaks",
            ResolutionSource::Default => "default",
        }
    }
}

const PREDICTOR_NCE_RANGE: (f64, f64) = (10.0, 60.0);

#[derive(Debug, Clone)]
pub(crate) struct Acquisition {
    pub(crate) activation: Option<ActivationMethod>,
    pub(crate) analyzer: Option<InstrumentType>,
    pub(crate) high_res: bool,
    pub(crate) resolution_source: ResolutionSource,
    pub(crate) peak_evidence: Option<input::ResolutionEvidence>,
    pub(crate) meta: input::Ms2MetaSample,
}

impl Acquisition {
    pub(crate) fn resolve(
        activation: Option<ActivationMethod>,
        analyzer: Option<InstrumentType>,
        fragment_tol_ppm: Option<f64>,
        fragment_tol_da: Option<f64>,
        peak_evidence: Option<input::ResolutionEvidence>,
        meta: input::Ms2MetaSample,
    ) -> Self {
        let inferred = peak_evidence.as_ref().and_then(|e| e.high_res);
        let (high_res, resolution_source) = if let Some(i) = analyzer {
            (i.is_high_resolution(), ResolutionSource::Metadata)
        } else if fragment_tol_ppm.is_some() {
            (true, ResolutionSource::FragmentTolFlag)
        } else if fragment_tol_da.is_some() {
            (false, ResolutionSource::FragmentTolFlag)
        } else if let Some(h) = inferred {
            (h, ResolutionSource::Peaks)
        } else {
            (false, ResolutionSource::Default)
        };
        Acquisition {
            activation,
            analyzer,
            high_res,
            resolution_source,
            peak_evidence,
            meta,
        }
    }

    pub(crate) fn inferred_high_res(&self) -> Option<bool> {
        self.peak_evidence.as_ref().and_then(|e| e.high_res)
    }

    fn is_tof(&self) -> bool {
        matches!(
            self.analyzer,
            Some(InstrumentType::TOF | InstrumentType::TimsTOF)
        )
    }

    pub(crate) fn predictor_nce(&self) -> Option<f64> {
        if self.activation != Some(ActivationMethod::HCD) || self.is_tof() {
            return None;
        }
        let nce = self.meta.median_collision_energy()?;
        (PREDICTOR_NCE_RANGE.0..=PREDICTOR_NCE_RANGE.1)
            .contains(&nce)
            .then_some(nce)
    }

    pub(crate) fn ms2_activation(&self) -> &'static str {
        if self.is_tof() {
            return "unknown";
        }
        let dissociation = match self.activation {
            Some(ActivationMethod::HCD) => "HCD",
            Some(ActivationMethod::CID) => "CID",
            Some(_) => return "unknown",
            None if self.high_res => "HCD",
            None => "CID",
        };
        match (dissociation, self.analyzer, self.high_res) {
            ("HCD", _, true) => "HCD-FTMS",
            ("HCD", _, false) => "HCD-ITMS",
            ("CID", _, true) => "CID-FTMS",
            _ => "CID-ITMS",
        }
    }

    pub(crate) fn ms2_instrument(&self) -> &'static str {
        match self.analyzer {
            Some(InstrumentType::QExactive) => "Q Exactive",
            Some(InstrumentType::OrbitrapAstral) => "Orbitrap Exploris 480",
            _ => "",
        }
    }

    pub(crate) fn predictor_instrument(&self, model_instrument: InstrumentType) -> InstrumentType {
        self.analyzer.unwrap_or(model_instrument)
    }

    pub(crate) fn describe(&self) -> String {
        let mut parts = vec![format!(
            "{} MS2 [{}]",
            if self.high_res {
                "high-resolution"
            } else {
                "low-resolution"
            },
            self.resolution_source.tag()
        )];
        if let Some(i) = self.analyzer {
            parts.push(format!("analyzer {}", i.name()));
        }
        if let Some((lo, hi)) = self.meta.collision_energy_range() {
            let med = self.meta.median_collision_energy().unwrap_or(lo);
            if (hi - lo).abs() < 1e-6 {
                parts.push(format!("collision energy {med:.1}"));
            } else {
                parts.push(format!(
                    "collision energy median {med:.1} (range {lo:.1}-{hi:.1})"
                ));
            }
        }
        if let Some(w) = self.meta.median_isolation_width() {
            parts.push(format!("isolation width {w:.2} Th"));
        }
        let cvs = self.meta.distinct_faims_cvs();
        if !cvs.is_empty() {
            let list: Vec<String> = cvs.iter().map(|v| format!("{v}")).collect();
            parts.push(format!("FAIMS CV {} V", list.join("/")));
        }
        if self.meta.n_ms2 > 0 {
            parts.push(format!("from {} MS2 scans", self.meta.n_ms2));
        }
        parts.join(", ")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn evidence(high_res: Option<bool>) -> input::ResolutionEvidence {
        input::ResolutionEvidence {
            spectra: 1000,
            peaks: 100_000,
            isotope_pairs: 20_000,
            chance_pairs: 400.0,
            enrichment: 50.0,
            high_res,
        }
    }

    fn meta(ce: &[f64]) -> input::Ms2MetaSample {
        input::Ms2MetaSample {
            n_ms2: ce.len(),
            collision_energies: ce.to_vec(),
            ..Default::default()
        }
    }

    #[test]
    fn metadata_analyzer_decides_the_class_whatever_the_model() {
        let a = Acquisition::resolve(
            Some(ActivationMethod::ETD),
            Some(InstrumentType::QExactive),
            None,
            Some(0.5),
            Some(evidence(Some(false))),
            meta(&[]),
        );
        assert!(a.high_res);
        assert_eq!(a.resolution_source, ResolutionSource::Metadata);
        let a = Acquisition::resolve(
            Some(ActivationMethod::CID),
            Some(InstrumentType::LowRes),
            Some(20.0),
            None,
            None,
            meta(&[]),
        );
        assert!(!a.high_res);
    }

    #[test]
    fn analyzer_without_activation_still_decides_the_class() {
        let a = Acquisition::resolve(
            None,
            Some(InstrumentType::OrbitrapAstral),
            None,
            None,
            Some(evidence(Some(false))),
            meta(&[]),
        );
        assert!(a.high_res);
        assert_eq!(a.resolution_source, ResolutionSource::Metadata);
        assert_eq!(a.activation, None);
        assert_eq!(a.analyzer, Some(InstrumentType::OrbitrapAstral));
    }

    #[test]
    fn flags_then_peaks_then_default() {
        let a = Acquisition::resolve(
            None,
            None,
            Some(10.0),
            None,
            Some(evidence(Some(false))),
            meta(&[]),
        );
        assert!(a.high_res);
        assert_eq!(a.resolution_source, ResolutionSource::FragmentTolFlag);
        let a = Acquisition::resolve(
            None,
            None,
            None,
            None,
            Some(evidence(Some(true))),
            meta(&[]),
        );
        assert!(a.high_res);
        assert_eq!(a.resolution_source, ResolutionSource::Peaks);
        assert_eq!(a.inferred_high_res(), Some(true));
        let a = Acquisition::resolve(None, None, None, None, Some(evidence(None)), meta(&[]));
        assert!(!a.high_res);
        assert_eq!(a.resolution_source, ResolutionSource::Default);
    }

    #[test]
    fn predictor_nce_comes_from_an_hcd_run_outside_a_tof() {
        let orbi = Acquisition::resolve(
            Some(ActivationMethod::HCD),
            Some(InstrumentType::OrbitrapAstral),
            None,
            None,
            None,
            meta(&[25.0, 25.0, 27.0]),
        );
        assert_eq!(orbi.predictor_nce(), Some(25.0));
        assert_eq!(
            orbi.predictor_instrument(InstrumentType::QExactive),
            InstrumentType::OrbitrapAstral
        );
        let tims = Acquisition::resolve(
            Some(ActivationMethod::CID),
            Some(InstrumentType::TimsTOF),
            None,
            None,
            None,
            meta(&[42.0]),
        );
        assert_eq!(tims.predictor_nce(), None);
        assert_eq!(
            tims.predictor_instrument(InstrumentType::QExactive),
            InstrumentType::TimsTOF
        );
        let odd = Acquisition::resolve(
            Some(ActivationMethod::HCD),
            Some(InstrumentType::QExactive),
            None,
            None,
            None,
            meta(&[150.0]),
        );
        assert_eq!(odd.predictor_nce(), None);
        let cid = Acquisition::resolve(
            Some(ActivationMethod::CID),
            Some(InstrumentType::LowRes),
            None,
            None,
            None,
            meta(&[35.0]),
        );
        assert_eq!(cid.predictor_nce(), None);
    }

    #[test]
    fn describe_names_the_class_and_its_source() {
        let a = Acquisition::resolve(
            Some(ActivationMethod::HCD),
            Some(InstrumentType::QExactive),
            None,
            None,
            None,
            input::Ms2MetaSample {
                n_ms2: 3,
                collision_energies: vec![28.0, 28.0, 28.0],
                isolation_widths: vec![1.4, 1.4, 1.4],
                faims_cvs: vec![-45.0, -65.0],
            },
        );
        let d = a.describe();
        assert!(d.starts_with("high-resolution MS2 [detected]"), "{d}");
        assert!(d.contains("collision energy 28.0"), "{d}");
        assert!(d.contains("isolation width 1.40 Th"), "{d}");
        assert!(d.contains("FAIMS CV -65/-45 V"), "{d}");
    }

    #[test]
    fn ms2_conditions_follow_the_acquisition() {
        let resolve = |act, analyzer, tol_da| {
            Acquisition::resolve(Some(act), analyzer, None, tol_da, None, meta(&[]))
        };
        let astral = resolve(
            ActivationMethod::HCD,
            Some(InstrumentType::OrbitrapAstral),
            None,
        );
        assert_eq!(astral.ms2_activation(), "HCD-FTMS");
        assert_eq!(astral.ms2_instrument(), "Orbitrap Exploris 480");
        let qe = resolve(ActivationMethod::HCD, Some(InstrumentType::QExactive), None);
        assert_eq!(qe.ms2_instrument(), "Q Exactive");
        let trap = resolve(ActivationMethod::CID, Some(InstrumentType::LowRes), None);
        assert_eq!(trap.ms2_activation(), "CID-ITMS");
        assert_eq!(trap.ms2_instrument(), "");
        let tims = resolve(ActivationMethod::CID, Some(InstrumentType::TimsTOF), None);
        assert_eq!(tims.ms2_activation(), "unknown");
        assert_eq!(
            resolve(ActivationMethod::ETD, Some(InstrumentType::QExactive), None).ms2_activation(),
            "unknown"
        );
    }
}
