//! Tests of the retrieval-mode choice.

use super::retrieval_choice;
use crate::cli::FragmentIndexFlag;
use model::tolerance::Tolerance;

const HI: Tolerance = Tolerance::Ppm(20.0);
const LO: Tolerance = Tolerance::Da(0.5);

#[test]
fn high_res_uses_the_index() {
    let c = retrieval_choice(HI, None, FragmentIndexFlag::Auto);
    assert!(c.use_index);
    assert!(c.report.unwrap().contains("high-resolution"));
}

#[test]
fn low_res_uses_the_capped_index() {
    let c = retrieval_choice(LO, None, FragmentIndexFlag::Auto);
    assert!(c.use_index);
    assert!(c.report.unwrap().contains("most intense peaks vote"));
}

#[test]
fn off_is_honoured_on_low_res() {
    let c = retrieval_choice(LO, None, FragmentIndexFlag::Off);
    assert!(!c.use_index);
}

#[test]
fn glyco_never_uses_the_fragment_index() {
    const GLYCO: Option<&str> = Some("--glyco looks candidates up by backbone mass");
    for flag in [FragmentIndexFlag::Auto, FragmentIndexFlag::On] {
        let c = retrieval_choice(HI, GLYCO, flag);
        assert!(!c.use_index, "glyco looks candidates up by backbone mass");
    }
    let c = retrieval_choice(HI, GLYCO, FragmentIndexFlag::On);
    assert_eq!(
        c.refused_because,
        Some("--glyco looks candidates up by backbone mass")
    );
}

#[test]
fn off_is_honoured_where_the_index_would_otherwise_run() {
    let c = retrieval_choice(HI, None, FragmentIndexFlag::Off);
    assert!(!c.use_index);
    assert!(c.report.unwrap().contains("switched off"));
}

#[test]
fn activations_outside_the_indexed_ion_series_keep_enumeration() {
    use model::activation::ActivationMethod;
    for a in [
        ActivationMethod::CID,
        ActivationMethod::HCD,
        ActivationMethod::PQD,
        ActivationMethod::EThcD,
    ] {
        assert!(search::fragments::index::covers_activation(a), "{a:?}");
    }
    for a in [ActivationMethod::ETD, ActivationMethod::UVPD] {
        assert!(!search::fragments::index::covers_activation(a), "{a:?}");
    }
    const ETD: Option<&str> = Some("b/y only");
    for tol in [HI, LO] {
        assert!(!retrieval_choice(tol, ETD, FragmentIndexFlag::Auto).use_index);
        let forced = retrieval_choice(tol, ETD, FragmentIndexFlag::On);
        assert!(!forced.use_index);
        assert_eq!(forced.refused_because, ETD);
    }
}

#[test]
fn a_tight_dalton_tolerance_counts_as_high_resolution() {
    let tight = retrieval_choice(Tolerance::Da(0.02), None, FragmentIndexFlag::Auto);
    assert!(tight.report.unwrap().contains("high-resolution"));
    let wide = retrieval_choice(Tolerance::Da(0.2), None, FragmentIndexFlag::Auto);
    assert!(wide.report.unwrap().contains("low-resolution"));
}
