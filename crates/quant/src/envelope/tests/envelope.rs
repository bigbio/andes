use super::*;
use crate::ms1::Ms1Scan;
use model::Tolerance;

const ENVELOPE: [f64; 7] = [0.08, 0.20, 0.27, 0.22, 0.13, 0.07, 0.03];
const MONO: f64 = 912.4567;
const CHARGE: u8 = 4;
const HEIGHT: f64 = 1e6;
const SIGMA_S: f64 = 4.0;
/// The prior on the apex around an anchor that is within the elution peak.
const PRIOR_SD_S: f64 = 2.0 * SIGMA_S;

struct Elution {
    apex_s: f64,
    height: f64,
    isotopes: Vec<(usize, f64)>,
}

fn pattern(apex_s: f64, height: f64) -> Elution {
    Elution {
        apex_s,
        height,
        isotopes: ENVELOPE.iter().copied().enumerate().collect(),
    }
}

fn run(elutions: &[Elution]) -> Ms1RunIndex {
    let scans = (0..240)
        .map(|i| {
            let rt = i as f64;
            let mut peaks: Vec<(f64, f32)> = vec![(400.0, 50.0), (1500.0, 50.0)];
            for k in 0..ENVELOPE.len() {
                let mz = MONO + k as f64 * ISOTOPE / f64::from(CHARGE);
                let mut inten = 0.0;
                for e in elutions {
                    let g = (-0.5 * ((rt - e.apex_s) / SIGMA_S).powi(2)).exp();
                    for &(iso, share) in &e.isotopes {
                        if iso == k {
                            inten += e.height * share * g;
                        }
                    }
                }
                if inten > 1.0 {
                    peaks.push((mz, inten as f32));
                }
            }
            peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
            Ms1Scan {
                rt,
                peaks,
                ..Ms1Scan::default()
            }
        })
        .collect();
    Ms1RunIndex::new(scans)
}

fn target(anchor_rt: f64) -> LfqTarget {
    LfqTarget {
        mono_mz: MONO,
        charge: CHARGE,
        anchor_rt,
        envelope: ENVELOPE.to_vec(),
        isotope_share: 1.0,
    }
}

fn params() -> LfqParams {
    LfqParams {
        tol: Tolerance::Ppm(10.0),
        rt_window_s: 60.0,
        ..LfqParams::default()
    }
}

/// The area under a Gaussian elution of unit envelope: height * sigma * sqrt(2 pi).
fn gaussian_area(height: f64) -> f64 {
    height * SIGMA_S * (2.0 * std::f64::consts::PI).sqrt()
}

#[test]
fn a_wide_envelope_with_a_weak_monoisotopic_peak_is_quantified_on_all_isotopes() {
    let index = run(&[pattern(100.0, HEIGHT)]);
    let t = target(97.0);
    let f = quantify_envelope(&index, &t, t.mono_mz, &params(), PRIOR_SD_S).expect("feature");
    assert!((f.apex_rt - 100.0).abs() <= 1.0, "{f:?}");
    assert!(f.cosine > 0.999, "{f:?}");
    assert_eq!(f.n_isotopes_found, 7);
    let expected = gaussian_area(HEIGHT);
    assert!(
        (f.area - expected).abs() / expected < 0.03,
        "{} vs {expected}",
        f.area
    );
    assert!((f.rt_delta_s - 3.0).abs() <= 1.0, "{f:?}");
}

#[test]
fn the_area_is_the_amplitude_over_the_isotope_share() {
    let index = run(&[pattern(100.0, HEIGHT)]);
    let whole = target(100.0);
    let part = LfqTarget {
        isotope_share: 0.8,
        ..whole.clone()
    };
    let a = quantify_envelope(&index, &whole, MONO, &params(), PRIOR_SD_S).expect("feature");
    let b = quantify_envelope(&index, &part, MONO, &params(), PRIOR_SD_S).expect("feature");
    assert!((b.area - a.area / 0.8).abs() / a.area < 1e-12);
}

#[test]
fn an_isotope_another_ion_adds_to_does_not_inflate_the_quantity() {
    let clean = run(&[pattern(100.0, HEIGHT)]);
    let interfered = run(&[
        pattern(100.0, HEIGHT),
        Elution {
            apex_s: 100.0,
            height: 3.0 * HEIGHT,
            isotopes: vec![(3, ENVELOPE[3])],
        },
    ]);
    let t = target(100.0);
    let a = quantify_envelope(&clean, &t, MONO, &params(), PRIOR_SD_S).expect("feature");
    let b = quantify_envelope(&interfered, &t, MONO, &params(), PRIOR_SD_S).expect("feature");
    let summed: f64 = b.isotope_areas.iter().sum();
    assert!(
        summed > 1.5 * a.area,
        "the plain sum is inflated: {summed} vs {}",
        a.area
    );
    assert!(
        (b.area - a.area).abs() / a.area < 0.05,
        "{} vs {}",
        b.area,
        a.area
    );
    assert!(b.cosine < a.cosine);
}

#[test]
fn a_predicted_anchor_finds_the_pattern_near_it_and_not_a_larger_one_far_away() {
    let index = run(&[pattern(100.0, HEIGHT), pattern(150.0, 20.0 * HEIGHT)]);
    let t = target(108.0);
    let wide = LfqParams {
        rt_window_s: 80.0,
        ..params()
    };
    let near = quantify_envelope(&index, &t, MONO, &wide, 5.0).expect("feature");
    assert!((near.apex_rt - 100.0).abs() <= 1.0, "{near:?}");
    let far = target(145.0);
    let other = quantify_envelope(&index, &far, MONO, &wide, 5.0).expect("feature");
    assert!((other.apex_rt - 150.0).abs() <= 1.0, "{other:?}");
}

#[test]
fn a_scan_where_one_isotope_alone_has_signal_is_not_taken_for_the_pattern() {
    let index = run(&[
        pattern(110.0, HEIGHT),
        Elution {
            apex_s: 95.0,
            height: HEIGHT,
            isotopes: vec![(2, 1.0)],
        },
    ]);
    let t = target(102.0);
    let f = quantify_envelope(&index, &t, MONO, &params(), 10.0).expect("feature");
    assert!((f.apex_rt - 110.0).abs() <= 1.0, "{f:?}");
}

#[test]
fn a_shifted_decoy_mass_finds_nothing() {
    let index = run(&[pattern(100.0, HEIGHT)]);
    let t = target(100.0);
    let p = params();
    assert!(quantify_envelope(&index, &t, MONO + p.decoy_mz_shift, &p, PRIOR_SD_S).is_none());
}

#[test]
fn the_apex_of_a_noisy_peak_is_found_from_an_anchor_on_its_flank() {
    let scans = (0..240)
        .map(|i| {
            let rt = i as f64;
            let g = (-0.5 * ((rt - 100.0) / SIGMA_S).powi(2)).exp();
            let jitter = if i % 2 == 0 { 1.5 } else { 0.5 };
            let peaks = ENVELOPE
                .iter()
                .enumerate()
                .map(|(k, share)| {
                    let mz = MONO + k as f64 * ISOTOPE / f64::from(CHARGE);
                    (mz, (HEIGHT * share * g * jitter) as f32)
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
    let index = Ms1RunIndex::new(scans);
    let smoothed = LfqParams {
        gaussian_width_s: Some(8.0 * 0.5 * SIGMA_S),
        ..params()
    };
    let f =
        quantify_envelope(&index, &target(94.0), MONO, &smoothed, 2.0 * SIGMA_S).expect("feature");
    assert!((f.apex_rt - 100.0).abs() <= 2.0, "{f:?}");
    assert!(f.rt_start <= 93.0 && f.rt_stop >= 107.0, "{f:?}");
    let expected = gaussian_area(HEIGHT);
    assert!(
        (f.area - expected).abs() / expected < 0.1,
        "{} vs {expected}",
        f.area
    );
}

fn gaussian_trace(sd_s: f64, step_s: f64, noise: bool) -> (Vec<f64>, Vec<f32>) {
    let n = (120.0 / step_s) as usize;
    (0..n)
        .map(|i| {
            let rt = i as f64 * step_s;
            let g = (-0.5 * ((rt - 60.0) / sd_s).powi(2)).exp();
            let jitter = if !noise {
                1.0
            } else if i % 2 == 0 {
                1.6
            } else {
                0.4
            };
            (rt, (1e5 * g * jitter) as f32)
        })
        .unzip()
}

#[test]
fn the_elution_width_is_read_from_the_autocorrelation() {
    for sd in [2.0, 5.0, 9.0] {
        let (rt, signal) = gaussian_trace(sd, 0.5, false);
        let got = elution_sd_s(&rt, &signal, 0.5).expect("width");
        assert!((got - sd).abs() / sd < 0.05, "{got} vs {sd}");
    }
}

#[test]
fn the_elution_width_does_not_depend_on_scan_to_scan_noise_or_on_the_scan_spacing() {
    let (rt, signal) = gaussian_trace(5.0, 0.5, true);
    let noisy = elution_sd_s(&rt, &signal, 0.5).expect("width");
    assert!((noisy - 5.0).abs() / 5.0 < 0.1, "{noisy}");
    let uneven: (Vec<f64>, Vec<f32>) = {
        let (rt, signal) = gaussian_trace(5.0, 0.25, false);
        rt.into_iter()
            .zip(signal)
            .enumerate()
            .filter(|(i, _)| i % 3 != 1)
            .map(|(_, p)| p)
            .unzip()
    };
    let got = elution_sd_s(&uneven.0, &uneven.1, 0.5).expect("width");
    assert!((got - 5.0).abs() / 5.0 < 0.05, "{got}");
}

#[test]
fn a_trace_without_a_peak_has_no_elution_width() {
    let rt: Vec<f64> = (0..200).map(|i| i as f64 * 0.5).collect();
    assert_eq!(elution_sd_s(&rt, &vec![0.0; 200], 0.5), None);
    assert_eq!(elution_sd_s(&rt, &vec![100.0; 200], 0.5), None);
    assert_eq!(elution_sd_s(&rt[..4], &[1.0, 2.0, 1.0, 0.0], 0.5), None);
    assert_eq!(elution_sd_s(&rt, &[1.0, 2.0], 0.5), None);
    assert_eq!(elution_sd_s(&rt, &vec![1.0; 200], 0.0), None);
    assert_eq!(elution_sd_s(&[], &[], 0.5), None);
}

#[test]
fn the_amplitude_of_exact_areas_is_their_sum_and_ignores_one_inflated_isotope() {
    let areas: Vec<f64> = ENVELOPE.iter().map(|e| 5e6 * e).collect();
    assert!((envelope_amplitude(&areas, &ENVELOPE) - 5e6).abs() < 1e-3);
    let mut inflated = areas.clone();
    inflated[1] *= 4.0;
    assert!((envelope_amplitude(&inflated, &ENVELOPE) - 5e6).abs() < 1e-3);
    assert_eq!(envelope_amplitude(&[0.0; 7], &ENVELOPE), 0.0);
    assert_eq!(envelope_amplitude(&areas, &[0.0; 7]), 0.0);
}
