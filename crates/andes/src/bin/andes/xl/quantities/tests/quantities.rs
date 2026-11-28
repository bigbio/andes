//! Tests of which cross-link quantities are accepted and of their q-values.

use super::*;

fn settings() -> Settings {
    Settings {
        params: LfqParams::default(),
        rt_window_s: None,
        xl_fdr: 0.01,
        transfer_fdr: 0.05,
        mbr: true,
    }
}

fn feature(cosine: f32, apex_intensity: f32) -> FeatureQuant {
    FeatureQuant {
        apex_rt: 100.0,
        rt_start: 90.0,
        rt_stop: 110.0,
        area: 1e6,
        apex_intensity,
        isotope_areas: vec![1.0; 4],
        cosine,
        n_isotopes_found: 4,
        rt_delta_s: 0.0,
        n_scans: 20,
    }
}

fn row(evidence: Evidence, target: Option<FeatureQuant>, decoy: Option<FeatureQuant>) -> Row {
    Row {
        precursor: 0,
        run: 0,
        evidence,
        anchor_rt: 100.0,
        anchor_sd_s: 5.0,
        window_s: 40.0,
        scan: 0,
        csm_score: 1.0,
        n_csms: 1,
        feature: target,
        decoy,
        score: None,
        q: None,
    }
}

#[test]
fn an_identification_is_accepted_when_its_isotope_pattern_fits() {
    let s = settings();
    let mut rows = vec![
        row(Evidence::Identified, Some(feature(0.95, 1e5)), None),
        row(Evidence::Identified, Some(feature(0.40, 1e5)), None),
        row(Evidence::Identified, None, Some(feature(0.95, 1e5))),
    ];
    score_rows(&mut rows, &s);
    assert!(accepted(&rows[0], &s));
    assert!(!accepted(&rows[1], &s), "below the minimum cosine");
    assert!(!accepted(&rows[2], &s), "no feature");
    assert_eq!(rows[1].score, None);
}

#[test]
fn a_transfer_is_accepted_at_the_transfer_fdr_of_its_competition() {
    let s = settings();
    let mut rows: Vec<Row> = (0..40)
        .map(|_| {
            row(
                Evidence::Transferred,
                Some(feature(0.99, 1e6)),
                Some(feature(0.75, 1e3)),
            )
        })
        .collect();
    rows.push(row(
        Evidence::Transferred,
        Some(feature(0.80, 1e3)),
        Some(feature(0.99, 1e6)),
    ));
    rows.push(row(Evidence::Transferred, Some(feature(0.99, 1e6)), None));
    score_rows(&mut rows, &s);
    assert!(
        rows[..40].iter().all(|r| accepted(r, &s)),
        "{:?}",
        rows[0].q
    );
    assert!(!accepted(&rows[40], &s), "its decoy wins: {:?}", rows[40].q);
    assert!(accepted(&rows[41], &s));
    let strict = Settings {
        transfer_fdr: 0.0,
        ..settings()
    };
    assert!(rows.iter().all(|r| !accepted(r, &strict)));
}

#[test]
fn the_evidence_is_written_as_the_benchmark_tables_read_it() {
    assert_eq!(Evidence::Identified.as_str(), "identified");
    assert_eq!(Evidence::Transferred.as_str(), "transferred");
}
