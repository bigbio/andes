//! Tests of the cross-link MS1 measurement: the elution width of a run, what follows from it,
//! and the decoy each row competes with.

use model::mass::ISOTOPE;
use model::Tolerance;
use quant::ms1::Ms1Scan;

use super::super::library::{Precursor, Side};
use super::*;

const ENVELOPE: [f64; 5] = [0.15, 0.30, 0.28, 0.17, 0.10];
const MONO: f64 = 845.4321;
const CHARGE: u8 = 4;
const STEP_S: f64 = 0.5;
const PEAK_SD_S: f64 = 4.0;

fn settings() -> Settings {
    Settings {
        params: LfqParams {
            tol: Tolerance::Ppm(10.0),
            ..LfqParams::default()
        },
        rt_window_s: None,
        xl_fdr: 0.01,
        transfer_fdr: 0.05,
        mbr: true,
    }
}

fn side(peptide: &str) -> Side {
    Side {
        peptide: peptide.into(),
        bare: peptide.into(),
        site: 2,
        site_kind: "K".into(),
    }
}

fn library() -> Library {
    Library {
        precursors: vec![Precursor {
            key: "DSS|AKA@2|GKG@2/4".into(),
            linker: "DSS".into(),
            sides: [side("AKA"), side("GKG")],
            charge: CHARGE,
            mono_mz: MONO,
            envelope: ENVELOPE.to_vec(),
            isotope_share: 0.9,
            residue_pair: "P:2 -- P:9".into(),
            protein_pair: "P -- P".into(),
            link: "intra",
            q: 0.0,
        }],
        identifications: Vec::new(),
    }
}

/// A run of 600 s with the precursor eluting at each of `apexes`.
fn run(apexes: &[f64]) -> Ms1RunIndex {
    let scans = (0..(600.0 / STEP_S) as usize)
        .map(|i| {
            let rt = i as f64 * STEP_S;
            let g: f64 = apexes
                .iter()
                .map(|apex| (-0.5 * ((rt - apex) / PEAK_SD_S).powi(2)).exp())
                .sum();
            let peaks = ENVELOPE
                .iter()
                .enumerate()
                .map(|(k, share)| {
                    let mz = MONO + k as f64 * ISOTOPE / f64::from(CHARGE);
                    (mz, (1e6 * share * g) as f32)
                })
                .filter(|p| p.1 > 1.0)
                .collect();
            Ms1Scan {
                rt,
                peaks,
                ..Ms1Scan::default()
            }
        })
        .collect();
    Ms1RunIndex::new(scans)
}

fn row(evidence: Evidence, precursor: usize, run: usize, anchor_rt: f64) -> Row {
    Row {
        precursor,
        run,
        evidence,
        anchor_rt,
        anchor_sd_s: 0.0,
        window_s: 0.0,
        scan: 0,
        csm_score: f64::NAN,
        n_csms: 0,
        feature: None,
        decoy: None,
        score: None,
        q: None,
    }
}

#[test]
fn the_peak_scale_of_a_run_is_the_width_of_its_elution_peaks() {
    let ms1 = run(&[100.0, 300.0]);
    let targets = [
        target_of(MONO, CHARGE, &ENVELOPE, 1.0, 97.0),
        target_of(MONO, CHARGE, &ENVELOPE, 1.0, 304.0),
        target_of(MONO, CHARGE, &ENVELOPE, 1.0, 500.0),
    ];
    let scale = peak_scale(&ms1, &targets, &settings()).expect("scale");
    assert!(
        (scale.sd_s - PEAK_SD_S).abs() / PEAK_SD_S < 0.05,
        "{scale:?}"
    );
    assert!((scale.step_s - STEP_S).abs() < 1e-9);
    assert_eq!(peak_scale(&ms1, &targets[2..], &settings()), None);
    assert_eq!(peak_scale(&ms1, &[], &settings()), None);
}

#[test]
fn the_smoothing_and_the_climb_follow_the_peak_width() {
    let scale = PeakScale {
        sd_s: 6.0,
        step_s: 0.4,
    };
    let p = params_at(&settings(), scale, 50.0);
    assert_eq!(p.rt_window_s, 50.0);
    assert_eq!(p.gaussian_width_s, Some(8.0 * 0.5 * 6.0));
    assert_eq!(p.max_apex_climb, 15);
    let given = Settings {
        params: LfqParams {
            gaussian_width_s: Some(12.0),
            ..settings().params
        },
        ..settings()
    };
    assert_eq!(params_at(&given, scale, 50.0).gaussian_width_s, Some(12.0));
}

#[test]
fn an_identification_is_looked_for_within_a_peak_width_of_its_spectrum() {
    let scale = PeakScale {
        sd_s: 5.0,
        step_s: 0.5,
    };
    let (sd, window) = identified_prior(&settings(), scale);
    assert!((sd - FWHM_PER_SD * 5.0).abs() < 1e-12);
    assert!((window - (4.0 * sd + 4.0 * 5.0)).abs() < 1e-9);
    let given = Settings {
        rt_window_s: Some(33.0),
        ..settings()
    };
    assert_eq!(identified_prior(&given, scale), (sd, 33.0));
}

#[test]
fn a_transfer_is_looked_for_within_the_prediction_error_and_a_peak_width() {
    let scale = PeakScale {
        sd_s: 3.0,
        step_s: 0.5,
    };
    let t = Transfer {
        precursor: 7,
        run: 2,
        rt: 1234.0,
        sd_s: 4.0,
    };
    let r = transfer_row(&t, scale);
    assert_eq!(
        (r.precursor, r.run, r.evidence),
        (7, 2, Evidence::Transferred)
    );
    assert_eq!(r.anchor_rt, 1234.0);
    assert!((r.anchor_sd_s - 5.0).abs() < 1e-12);
    assert!((r.window_s - (4.0 * 5.0 + 4.0 * 3.0)).abs() < 1e-9);
    let exact = transfer_row(&Transfer { sd_s: 0.0, ..t }, scale);
    assert!((exact.anchor_sd_s - 3.0).abs() < 1e-12);
}

#[test]
fn an_identification_competes_with_a_shifted_mass_and_transfers_with_either_decoy() {
    let all = |r: &Row| {
        probes(
            r.precursor as u32,
            r.run as u16,
            MONO,
            CHARGE,
            r.anchor_rt,
            5.0,
            4.0,
            (0.0, 6000.0),
        )
    };
    let (mut by_mass, mut by_rt) = (0, 0);
    for precursor in 0..40 {
        for run in 0..5 {
            let id = row(Evidence::Identified, precursor, run, 3000.0);
            let probes = all(&id);
            let decoy = decoy_of(&id, &probes).expect("decoy");
            assert!(matches!(decoy.kind, Hypothesis::MassShiftDecoy(_)));
            let transfer = row(Evidence::Transferred, precursor, run, 3000.0);
            let decoy = decoy_of(&transfer, &probes).expect("decoy");
            assert_eq!(decoy_of(&transfer, &probes), Some(decoy));
            match decoy.kind {
                Hypothesis::MassShiftDecoy(_) => {
                    by_mass += 1;
                    assert_eq!(decoy.t_hat, 3000.0);
                    assert!((decoy.mz - MONO).abs() > 1.0);
                }
                Hypothesis::RtShiftDecoy => {
                    by_rt += 1;
                    assert_eq!(decoy.mz, MONO);
                    assert!((decoy.t_hat - 3000.0).abs() >= 8.0 * 5.0 + 4.0 * 4.0);
                }
                other => panic!("{other:?} is not a decoy"),
            }
        }
    }
    assert!(
        by_mass >= 60 && by_rt >= 60,
        "{by_mass} by mass, {by_rt} by retention time"
    );
}

#[test]
fn a_row_is_measured_at_its_peak_and_its_decoy_finds_nothing() {
    let ms1 = run(&[100.0]);
    let lib = library();
    let scale = PeakScale {
        sd_s: PEAK_SD_S,
        step_s: STEP_S,
    };
    let s = settings();
    let (sd, window) = identified_prior(&s, scale);
    let mut id = row(Evidence::Identified, 0, 0, 94.0);
    (id.anchor_sd_s, id.window_s) = (sd, window);
    let (feature, decoy) = measure_row(&ms1, &lib, &id, &s, scale);
    let feature = feature.expect("feature");
    assert!((feature.apex_rt - 100.0).abs() <= 1.0, "{feature:?}");
    assert!(feature.cosine > 0.999);
    let expected = 1e6 * PEAK_SD_S * (2.0 * std::f64::consts::PI).sqrt() / 0.9;
    assert!(
        (feature.area - expected).abs() / expected < 0.03,
        "{} vs {expected}",
        feature.area
    );
    assert_eq!(decoy, None);
    let far = transfer_row(
        &Transfer {
            precursor: 0,
            run: 0,
            rt: 400.0,
            sd_s: 3.0,
        },
        scale,
    );
    assert_eq!(measure_row(&ms1, &lib, &far, &s, scale), (None, None));
}

#[test]
fn a_linear_anchor_enters_with_the_apex_of_its_feature() {
    let ms1 = run(&[100.0]);
    let scale = PeakScale {
        sd_s: PEAK_SD_S,
        step_s: STEP_S,
    };
    let anchor = |rt: f64, mono_mz: f64| LinearAnchor {
        key: (1, CHARGE),
        rt,
        mono_mz,
        charge: CHARGE,
        envelope: ENVELOPE.to_vec(),
    };
    let apex = linear_apex(&ms1, &anchor(105.0, MONO), &settings(), scale).expect("apex");
    assert!((apex - 100.0).abs() <= 1.0, "{apex}");
    assert_eq!(
        linear_apex(&ms1, &anchor(105.0, MONO + 3.0), &settings(), scale),
        None
    );
}
