use super::*;
use crate::decoys::{Hypothesis, PEPTIDE_MASS_SLOPE};
use crate::ms1_index::Ms1Scan;
use model::mass::ISOTOPE;

type Feature = (f64, u8, [f32; 4], f64, f64, f64);

const ENV: [f32; 4] = [0.5, 0.3, 0.15, 0.05];
const MZ: f64 = 650.3271;
const AMP: f64 = 1.0e6;
const NOISE: f32 = 5.0e3;

struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        model::hash::splitmix64_next(&mut self.0)
    }

    fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    fn normal(&mut self) -> f64 {
        let (u, v) = (self.uniform(), self.uniform());
        (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
    }
}

const DENSE_CENTROIDS: usize = 200;

fn synth_scans(features: &[Feature], noise: f32, scan_dt: f64, seed: u64) -> Vec<Ms1Scan> {
    let noise = noise as f64;
    let end = features
        .iter()
        .map(|f| f.3 + 10.0 * f.4)
        .fold(60.0, f64::max)
        + 60.0;
    let mut rng = Rng(seed);
    let n = (end / scan_dt) as usize + 1;
    (0..n)
        .map(|k| {
            let t = k as f64 * scan_dt;
            let mut peaks = Vec::with_capacity(features.len() * 5 + DENSE_CENTROIDS);
            for _ in 0..DENSE_CENTROIDS {
                let mz = 300.0 + 1200.0 * rng.uniform();
                let inten = noise * (1.0 - rng.uniform().ln());
                if inten > 0.0 {
                    peaks.push((mz, inten as f32));
                }
            }
            for &(mz, z, env, apex, width, amp) in features {
                let g = (-(t - apex) * (t - apex) / (2.0 * width * width)).exp();
                for iso in -1i32..4 {
                    let signal = if iso >= 0 {
                        amp * env[iso as usize] as f64 * g
                    } else {
                        0.0
                    };
                    let inten = 5.0 * noise + signal + noise * rng.normal();
                    if inten > 0.0 {
                        peaks.push((mz + iso as f64 * ISOTOPE / z as f64, inten as f32));
                    }
                }
            }
            Ms1Scan {
                rt: t,
                peaks,
                ..Ms1Scan::default()
            }
        })
        .collect()
}

fn synth_run(features: &[Feature], noise: f32, scan_dt: f64) -> Ms1RunIndex {
    Ms1RunIndex::new(synth_scans(features, noise, scan_dt, 1))
}

fn probe(mz: f64, t_hat: f64, s: f64, w_prior: f64) -> Probe {
    Probe {
        kind: Hypothesis::Target,
        mz,
        t_hat,
        s,
        w_prior,
    }
}

fn params() -> EvidenceParams {
    EvidenceParams {
        tol_ppm: 10.0,
        n_isotopes: 4,
    }
}

fn assert_finite(e: &Evidence) {
    for (name, v) in [
        ("amplitude", e.amplitude),
        ("amplitude_se", e.amplitude_se),
        ("quantity", e.quantity),
        ("quantity_se", e.quantity_se),
        ("llr", e.llr),
        ("rt_z", e.rt_z),
        ("width_ratio", e.width_ratio),
        ("envelope_cos", e.envelope_cos),
        ("coherence", e.coherence),
        ("m1_ratio", e.m1_ratio),
        ("mass_err_ppm", e.mass_err_ppm),
        ("t0", e.t0),
        ("w", e.w),
    ] {
        assert!(v.is_finite(), "{name} = {v} is not finite");
    }
}

#[test]
fn isotope_heights_and_tailing_describe_the_peak_shape() {
    let apex = 300.0;
    let p = probe(MZ, apex, 6.0, 5.0);
    let clean = synth_run(&[(MZ, 2, ENV, apex, 5.0, AMP)], NOISE, 1.0);
    let e = fit(&clean, &p, 2, &ENV, &params()).unwrap();
    let sum: f64 = e.iso_heights.iter().sum();
    for (h, &pi) in e.iso_heights.iter().zip(&ENV) {
        assert!((h / sum - pi as f64).abs() < 0.02, "{:?}", e.iso_heights);
    }
    assert!(e.tailing.abs() < 0.05, "symmetric peak: {}", e.tailing);
    let tailed = synth_run(
        &[
            (MZ, 2, ENV, apex, 5.0, AMP),
            (MZ, 2, ENV, apex + 9.0, 8.0, 0.35 * AMP),
        ],
        NOISE,
        1.0,
    );
    let t = fit(&tailed, &p, 2, &ENV, &params()).unwrap();
    assert!(
        t.tailing > e.tailing + 0.1,
        "{} vs {}",
        t.tailing,
        e.tailing
    );
    let rt = [0.0, 1.0, 2.0];
    assert_eq!(
        tailing(&rt, &[1.0, 5.0, 1.0], &[0.5, 1.0, 0.5], 0.0, 1.0, 1.0),
        0.0
    );
}

#[test]
fn recovers_amplitude_apex_and_width() {
    let dt = 1.0;
    let apex = 300.4;
    let run = synth_run(&[(MZ, 2, ENV, apex, 5.0, AMP)], NOISE, dt);
    let e = fit(&run, &probe(MZ, 303.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
    eprintln!(
        "recover: A/A_true-1 = {:.4}, t0-apex = {:.3} s, width_ratio = {}, llr = {:.0}, \
         envelope_cos = {:.4}, coherence = {:.4}, m1_ratio = {:.4}, used = {:#06b}",
        e.amplitude / AMP - 1.0,
        e.t0 - apex,
        e.width_ratio,
        e.llr,
        e.envelope_cos,
        e.coherence,
        e.m1_ratio,
        e.isotopes_used
    );
    assert_finite(&e);
    assert!((e.amplitude / AMP - 1.0).abs() < 0.05);
    assert!((e.t0 - apex).abs() < dt);
    assert!((0.7..=1.4).contains(&e.width_ratio));
    assert_eq!(e.isotopes_used, 0b1111);
    assert!(e.envelope_cos > 0.99 && e.coherence > 0.9);
}

const NULL_LLR_BOUND: f64 = 10.8;

#[test]
fn width_is_recovered_between_grid_points() {
    for (k, w_true) in [3.3, 5.0, 7.1, 10.5].into_iter().enumerate() {
        for noise in [0.0, NOISE] {
            let run = Ms1RunIndex::new(synth_scans(
                &[(MZ, 2, ENV, 300.0, w_true, AMP)],
                noise,
                0.5,
                20 + k as u64,
            ));
            let e = fit(&run, &probe(MZ, 300.0, 6.0, 6.0), 2, &ENV, &params()).unwrap();
            eprintln!(
                "width: true {w_true}, fitted {:.3} (noise {noise}), A/A_true = {:.4}",
                e.w,
                e.amplitude / AMP
            );
            assert!((e.w / w_true - 1.0).abs() < 0.03, "{w_true}: {}", e.w);
            assert!((e.amplitude / AMP - 1.0).abs() < 0.03);
        }
    }
}

#[test]
fn golden_section_finds_the_maximum_and_never_loses_the_start() {
    let f = |x: f64| -(x - 0.3) * (x - 0.3);
    let x = golden_max(f, -1.0, 1.0, -1.0, f(-1.0));
    assert!((x - 0.3).abs() < 1e-3, "{x}");
    assert_eq!(golden_max(f, -1.0, 1.0, 0.3, 1.0), 0.3);
    assert_eq!(golden_max(f, 0.5, 0.5, 0.5, f(0.5)), 0.5);
    assert_eq!(width_bracket(&WIDTH_GRID, 1), (0.7, 1.4));
    let (lo, hi) = width_bracket(&WIDTH_GRID, 0);
    assert!((lo - 0.49).abs() < 1e-12 && hi == 1.0);
    let (lo, hi) = width_bracket(&WIDTH_GRID, 2);
    assert!(lo == 1.0 && (hi - 1.96).abs() < 1e-12);
    assert_eq!(width_bracket(&[1.0], 0), (1.0, 1.0));
}

#[test]
fn noise_estimator_scales_successive_differences() {
    let mut sort = Vec::new();
    let y: Vec<f32> = (0..21)
        .map(|k| if k % 2 == 0 { 100.0 } else { 110.0 })
        .collect();
    let expected = 1.4826 * 10.0 / std::f64::consts::SQRT_2;
    assert!((noise_sigma(&y, 1.0, &mut sort) - expected).abs() < 1e-9);
    assert_eq!(noise_sigma(&y, 50.0, &mut sort), 50.0);
    let y = [0.0f32, 1.0, 3.0, 6.0, 10.0];
    let expected = 1.4826 * 2.5 / std::f64::consts::SQRT_2;
    assert!((noise_sigma(&y, 0.0, &mut sort) - expected).abs() < 1e-9);
}

#[test]
fn fixture_noise_is_set_by_the_differences_not_the_floor() {
    let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, 0.0)], NOISE, 1.0);
    let scans = run.scan_window(250.0, 350.0);
    let limit = detection_limit(&run, scans.clone(), &mut Vec::new());
    let y = run.xic(MZ, MZ * 10e-6, scans);
    let sigma = noise_sigma(&y, 0.0, &mut Vec::new());
    eprintln!(
        "fixture: sigma = {sigma:.0}, 0.5 L = {:.0}, noise = {NOISE}",
        0.5 * limit
    );
    assert!(sigma > 0.5 * limit);
    assert!((sigma / NOISE as f64 - 1.0).abs() < 0.3);
}

#[test]
fn amplitude_unbiased_across_window_regimes() {
    let w = 5.0;
    let run = synth_run(&[(MZ, 2, ENV, 300.0, w, AMP)], NOISE, 1.0);
    for ratio in [0.25, 0.5, 1.0, 2.0] {
        let e = fit(&run, &probe(MZ, 300.0, ratio * w, w), 2, &ENV, &params()).unwrap();
        eprintln!(
            "regime s/w_prior = {ratio}: A/A_true = {:.4}, t0 = {}, w = {}",
            e.amplitude / AMP,
            e.t0,
            e.w
        );
        assert!(
            (e.amplitude / AMP - 1.0).abs() < 0.05,
            "s/w_prior = {ratio}"
        );
    }
}

#[test]
fn apex_candidates_are_spaced_by_the_narrowest_width() {
    let rt = |dt: f64| (0..100).map(|k| k as f64 * dt).collect::<Vec<_>>();
    let a = Apexes::of(&rt(1.0), 10..60, 3.5);
    assert_eq!(a.stride, 1);
    assert_eq!(a.iter().collect::<Vec<_>>(), (10..60).collect::<Vec<_>>());
    let a = Apexes::of(&rt(0.2), 10..60, 3.5);
    assert_eq!(a.stride, 8);
    let c: Vec<usize> = a.iter().collect();
    assert_eq!(c, [10, 18, 26, 34, 42, 50, 58, 59]);
    assert_eq!(
        a.around(26).collect::<Vec<_>>(),
        (19..26).chain(27..34).collect::<Vec<_>>()
    );
    assert_eq!(
        a.around(10).collect::<Vec<_>>(),
        (11..18).collect::<Vec<_>>()
    );
    assert_eq!(Apexes::of(&rt(0.2), 5..6, 3.5).stride, 1);
    assert_eq!(Apexes::of(&rt(0.2), 5..7, 3.5).stride, 1);
    assert_eq!(Apexes::of(&[0.0; 100], 0..50, 3.5).stride, 1);
}

#[test]
fn dense_scans_fit_the_same_peak() {
    let w = 6.0;
    let run = synth_run(&[(MZ, 2, ENV, 300.3, w, AMP)], NOISE, 0.3);
    let p = probe(MZ, 300.0, 8.0, w);
    let e = fit(&run, &p, 2, &ENV, &params()).unwrap();
    eprintln!(
        "dense: t0 = {}, w = {}, A/A_true = {:.4}",
        e.t0,
        e.w,
        e.amplitude / AMP
    );
    assert!((e.t0 - 300.3).abs() < 0.5, "t0 = {}", e.t0);
    assert!((e.w / w - 1.0).abs() < 0.1, "w = {}", e.w);
    assert!((e.amplitude / AMP - 1.0).abs() < 0.05);
}

#[test]
fn mass_error_is_the_weighted_centroid_offset() {
    let run = synth_run(&[(MZ * (1.0 + 3e-6), 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0);
    let e = fit(&run, &probe(MZ, 300.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
    eprintln!("mass error: {:.3} ppm (true 3)", e.mass_err_ppm);
    assert!((e.mass_err_ppm - 3.0).abs() < 0.2);
}

#[test]
fn absent_precursor_has_small_llr() {
    let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, 0.0)], NOISE, 1.0);
    let e = fit(&run, &probe(MZ, 300.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
    eprintln!(
        "absent: llr = {:.2}, amplitude = {:.0} ({:.4} of AMP, {:.2} SE)",
        e.llr,
        e.amplitude,
        e.amplitude / AMP,
        e.amplitude / e.amplitude_se
    );
    assert_finite(&e);
    assert!(e.llr < NULL_LLR_BOUND);
    assert!(e.amplitude < 0.05 * AMP);
}

#[test]
fn interfered_isotope_is_excluded_and_quantity_unbiased() {
    let w = 5.0;
    let iso2 = MZ + 2.0 * ISOTOPE / 2.0;
    let interferer = (
        iso2 * (1.0 + 3e-6),
        2,
        [1.0, 0.0, 0.0, 0.0],
        300.0,
        w,
        5.0 * AMP * ENV[2] as f64,
    );
    let run = synth_run(&[(MZ, 2, ENV, 300.0, w, AMP), interferer], NOISE, 1.0);
    let e = fit(&run, &probe(MZ, 300.0, 6.0, w), 2, &ENV, &params()).unwrap();
    let truth = AMP * w * (2.0 * std::f64::consts::PI).sqrt();
    eprintln!(
        "interfered: used = {:#06b}, quantity/true-1 = {:.4}, w = {}",
        e.isotopes_used,
        e.quantity / truth - 1.0,
        e.w
    );
    assert_finite(&e);
    assert_eq!(e.isotopes_used & 0b0100, 0);
    assert_eq!(e.isotopes_used & 0b0001, 1);
    assert!((e.quantity / truth - 1.0).abs() < 0.1);
}

#[test]
fn apex_is_refitted_after_an_isotope_is_excluded() {
    let (w, dt) = (5.0, 1.0);
    let iso1 = MZ + ISOTOPE / 2.0;
    let truth = AMP * w * (2.0 * std::f64::consts::PI).sqrt();
    for strength in [20.0, 100.0] {
        let interferer = (
            iso1 * (1.0 + 3e-6),
            2,
            [1.0, 0.0, 0.0, 0.0],
            114.0,
            w,
            strength * AMP * ENV[1] as f64,
        );
        let run = synth_run(&[(MZ, 2, ENV, 100.0, w, AMP), interferer], NOISE, dt);
        let e = fit(&run, &probe(MZ, 100.0, 6.0, w), 2, &ENV, &params()).unwrap();
        eprintln!(
            "M+1 interferer at 114 s ({strength}x): used = {:#06b}, t0 = {}, w = {}, \
             quantity/true = {:.4}",
            e.isotopes_used,
            e.t0,
            e.w,
            e.quantity / truth
        );
        assert_finite(&e);
        assert_eq!(e.isotopes_used & 0b0010, 0, "M+1 is excluded");
        assert!((e.t0 - 100.0).abs() <= dt, "t0 = {}", e.t0);
        assert!((e.quantity / truth - 1.0).abs() < 0.1);
    }
}

#[test]
fn rt_prior_prefers_peak_near_prediction() {
    let dt = 1.0;
    let run = synth_run(
        &[(MZ, 2, ENV, 300.0, 1.0, AMP), (MZ, 2, ENV, 305.0, 1.0, AMP)],
        0.0,
        dt,
    );
    let near = fit(&run, &probe(MZ, 300.0, 3.0, 1.0), 2, &ENV, &params()).unwrap();
    let far = fit(&run, &probe(MZ, 305.0, 3.0, 1.0), 2, &ENV, &params()).unwrap();
    eprintln!(
        "rt prior: t0 = {} (t_hat 300), {} (t_hat 305)",
        near.t0, far.t0
    );
    assert!((near.t0 - 300.0).abs() < dt);
    assert!((far.t0 - 305.0).abs() < dt);
}

#[test]
fn m1_ratio_flags_non_monoisotopic_hypothesis() {
    let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0);
    let mono = fit(&run, &probe(MZ, 300.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
    let plus1 = fit(
        &run,
        &probe(MZ + ISOTOPE / 2.0, 300.0, 6.0, 5.0),
        2,
        &ENV,
        &params(),
    )
    .unwrap();
    eprintln!(
        "m1_ratio: mono = {:.4}, +1 = {:.4}",
        mono.m1_ratio, plus1.m1_ratio
    );
    assert_finite(&mono);
    assert_finite(&plus1);
    assert!(mono.m1_ratio < 0.1);
    assert!(plus1.m1_ratio > 0.5);
}

#[test]
fn decoy_probe_in_empty_region_scores_below_target() {
    let decoy_mz = MZ + 7.0 * PEPTIDE_MASS_SLOPE / 2.0;
    let run = synth_run(
        &[
            (MZ, 2, ENV, 300.0, 5.0, AMP),
            (decoy_mz, 2, ENV, 300.0, 5.0, 0.0),
        ],
        NOISE,
        1.0,
    );
    let t = fit(&run, &probe(MZ, 300.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
    let mut dp = probe(decoy_mz, 300.0, 6.0, 5.0);
    dp.kind = Hypothesis::MassShiftDecoy(7);
    let d = fit(&run, &dp, 2, &ENV, &params()).unwrap();
    eprintln!("decoy: llr target = {:.0}, llr decoy = {:.2}", t.llr, d.llr);
    assert_finite(&d);
    assert!(t.llr > 20.0 * d.llr.max(1.0));
}

#[test]
fn extract_uses_donor_faims_group() {
    let at = |cv: f32, offset: f64, amp: f64, seed: u64| {
        synth_scans(&[(MZ, 2, ENV, 300.0, 5.0, amp)], NOISE, 2.0, seed)
            .into_iter()
            .map(move |mut s| {
                s.rt += offset;
                s.faims_cv = Some(cv);
                s
            })
    };
    let scans: Vec<Ms1Scan> = at(-45.0, 0.0, AMP, 1)
        .chain(at(-60.0, 1.0, 0.0, 2))
        .collect();
    let groups = Ms1RunIndex::new(scans).split_by_faims();
    let group = |cv: f32| &groups.iter().find(|g| g.0 == cv).unwrap().1;
    let pr = probe(MZ, 300.0, 6.0, 5.0);
    let on = fit(group(-45.0), &pr, 2, &ENV, &params()).unwrap();
    let off = fit(group(-60.0), &pr, 2, &ENV, &params()).unwrap();
    eprintln!(
        "faims: -45 llr = {:.0}, A/A_true = {:.4}; -60 llr = {:.2}",
        on.llr,
        on.amplitude / AMP,
        off.llr
    );
    assert!((on.amplitude / AMP - 1.0).abs() < 0.1);
    assert!(off.llr < NULL_LLR_BOUND);
}

#[test]
fn short_window_is_none_and_empty_signal_is_zero() {
    let scans = |n: usize, peaks: Vec<(f64, f32)>| {
        Ms1RunIndex::new(
            (0..n)
                .map(|k| Ms1Scan {
                    rt: k as f64,
                    peaks: peaks.clone(),
                    ..Ms1Scan::default()
                })
                .collect(),
        )
    };
    let pr = probe(MZ, 10.5, 0.25, 0.25);
    assert_eq!(fit(&scans(30, vec![]), &pr, 2, &ENV, &params()), None);
    let wide = probe(MZ, 10.0, 3.0, 2.0);
    let e = fit(&scans(30, vec![(400.0, 1e5)]), &wide, 2, &ENV, &params()).unwrap();
    assert_finite(&e);
    assert_eq!(e.amplitude, 0.0);
    assert_eq!(e.llr, 0.0);
    assert_eq!(e.coherence, 0.0);
    assert!(e.n_points > 0 && e.isotopes_used & 1 == 1);
}

#[test]
fn degenerate_probe_and_envelope_stay_finite() {
    let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0);
    let wide = probe(MZ, 300.0, 6.0, 5.0);
    assert_eq!(
        fit(&run, &probe(MZ, 300.0, 0.0, 0.0), 2, &ENV, &params()),
        None
    );
    for (pr, env) in [
        (probe(MZ, 300.0, 0.0, 5.0), ENV),
        (probe(MZ, 300.0, f64::NAN, 5.0), ENV),
        (wide, [0.0; 4]),
        (wide, [1.0, 0.0, 0.0, 0.0]),
        (wide, [f32::NAN, 1.0, 0.0, 0.0]),
    ] {
        let e = fit(&run, &pr, 2, &env, &params()).unwrap();
        assert_finite(&e);
    }
    let one = EvidenceParams {
        tol_ppm: 10.0,
        n_isotopes: 1,
    };
    let e = fit(&run, &wide, 2, &ENV, &one).unwrap();
    assert_finite(&e);
    assert_eq!(e.isotopes_used, 1);
}

#[test]
fn fit_is_deterministic() {
    let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0);
    let pr = probe(MZ, 301.0, 6.0, 5.0);
    let a = fit(&run, &pr, 2, &ENV, &params());
    assert!(a.is_some());
    assert_eq!(a, fit(&run, &pr, 2, &ENV, &params()));
}

#[test]
fn missing_isotope_trace_does_not_collapse_amplitude() {
    let scans = synth_scans(&[(MZ, 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0, 5);
    let m3 = MZ + 3.0 * ISOTOPE / 2.0;
    let without_m3: Vec<Ms1Scan> = scans
        .iter()
        .cloned()
        .map(|mut s| {
            s.peaks.retain(|p| (p.0 - m3).abs() > 0.01);
            s
        })
        .collect();
    let pr = probe(MZ, 300.0, 6.0, 5.0);
    let intact = fit(&Ms1RunIndex::new(scans), &pr, 2, &ENV, &params()).unwrap();
    let missing = fit(&Ms1RunIndex::new(without_m3), &pr, 2, &ENV, &params()).unwrap();
    eprintln!(
        "missing M+3: A/A_true = {:.4} (intact {:.4}), llr = {:.0} (intact {:.0}), used = {:#06b}",
        missing.amplitude / AMP,
        intact.amplitude / AMP,
        missing.llr,
        intact.llr,
        missing.isotopes_used
    );
    assert_finite(&missing);
    assert!((missing.amplitude / AMP - 1.0).abs() < 0.15);
    assert_eq!(missing.isotopes_used & 0b0010, 0b0010, "isotope 1 kept");
    assert_eq!(missing.isotopes_used & 0b1000, 0, "missing M+3 dropped");
    assert!(missing.llr > 0.5 * intact.llr && missing.llr < 2.0 * intact.llr);
}

fn every(range: std::ops::Range<usize>) -> Apexes {
    Apexes { range, stride: 1 }
}

#[test]
fn cached_grid_profiles_are_the_direct_profiles() {
    let rt: Vec<f64> = (0..400)
        .map(|k| 100.0 + 0.7 * k as f64 + 0.01 * (k % 3) as f64)
        .collect();
    let prior = Prior {
        t_hat: 240.0,
        s: 10.0,
        w_prior: 2.0,
    };
    let mut cache = GridCache::default();
    let g = cache
        .get(&rt, &prior, &WIDTH_GRID, &every(150..250))
        .unwrap();
    assert_eq!(g.cells.len(), 100 * WIDTH_GRID.len());
    let mut cells = g.cells.iter();
    let mut outside = 0;
    for k in 150..250 {
        for &f in &WIDTH_GRID {
            let cell = cells.next().unwrap();
            let (mut want_sg, mut want_sgg) = (0.0, 0.0);
            for (j, &t) in rt.iter().enumerate() {
                let want = gaussian::gaussian(t, rt[k], f * prior.w_prior);
                if cell.scans.contains(&j) {
                    let got = g.g[cell.at + j - cell.scans.start];
                    assert_eq!(got.to_bits(), want.to_bits());
                    want_sg += want;
                    want_sgg += want * want;
                } else {
                    assert!(want < 3e-18, "scan {j} of candidate {k}: {want}");
                    outside += 1;
                }
            }
            assert_eq!(
                (cell.sg.to_bits(), cell.sgg.to_bits()),
                (want_sg.to_bits(), want_sgg.to_bits())
            );
        }
    }
    assert!(cells.next().is_none());
    assert!(outside > 0, "the test exercises the support");
    let apexes = Apexes {
        range: 150..250,
        stride: 4,
    };
    let mut cache = GridCache::default();
    let g = cache.get(&rt, &prior, &WIDTH_GRID, &apexes).unwrap();
    assert_eq!(g.cells.len(), 26 * WIDTH_GRID.len());
    let mut cells = g.cells.iter();
    for k in apexes.iter() {
        for &f in &WIDTH_GRID {
            let cell = cells.next().unwrap();
            let j = cell.scans.start + (cell.scans.len() / 2);
            let want = gaussian::gaussian(rt[j], rt[k], f * prior.w_prior);
            assert_eq!(
                g.g[cell.at + j - cell.scans.start].to_bits(),
                want.to_bits()
            );
        }
    }
}

#[test]
fn grid_cache_hits_only_on_the_same_window_prior_and_widths() {
    let rt: Vec<f64> = (0..30).map(|k| k as f64).collect();
    let prior = Prior {
        t_hat: 15.0,
        s: 2.0,
        w_prior: 3.0,
    };
    let mut cache = GridCache::default();
    let first = cache
        .get(&rt, &prior, &WIDTH_GRID, &every(9..22))
        .unwrap()
        .g
        .clone();
    assert_eq!(cache.entries.len(), 1);
    assert_eq!(
        cache
            .get(&rt, &prior, &WIDTH_GRID, &every(9..22))
            .unwrap()
            .g,
        first
    );
    assert_eq!(cache.entries.len(), 1, "a hit computes nothing new");
    let wider = Prior {
        w_prior: 3.5,
        ..prior
    };
    assert_ne!(
        cache
            .get(&rt, &wider, &WIDTH_GRID, &every(9..22))
            .unwrap()
            .g,
        first
    );
    let mut shifted = rt.clone();
    shifted[29] += 0.5;
    cache.get(&shifted, &prior, &WIDTH_GRID, &every(9..22));
    cache.get(&rt, &prior, &FIVE_POINT_GRID, &every(9..22));
    assert_eq!(cache.entries.len(), GRID_CACHE_ENTRIES);
    assert!(
        cache.entries.iter().all(|e| !e.matches(
            &rt,
            [
                prior.t_hat.to_bits(),
                prior.s.to_bits(),
                prior.w_prior.to_bits()
            ],
            &WIDTH_GRID
        )),
        "the oldest entry was replaced"
    );
}

#[test]
fn fits_sharing_a_workspace_equal_fresh_fits() {
    let interferer = (
        MZ + 2.0 * ISOTOPE / 2.0,
        1,
        [1.0, 0.0, 0.0, 0.0],
        330.0,
        4.0,
        AMP,
    );
    let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, AMP), interferer], NOISE, 0.5);
    let mut scratch = FitScratch::default();
    for t_hat in [290.0, 300.0, 300.0, 318.0] {
        for mz in [MZ, MZ + 7.0 * PEPTIDE_MASS_SLOPE / 2.0, MZ] {
            let pr = probe(mz, t_hat, 12.0, 5.0);
            let shared = fit_with(&run, &pr, 2, &ENV, &params(), &mut scratch);
            let fresh = fit(&run, &pr, 2, &ENV, &params());
            assert_eq!(format!("{shared:?}"), format!("{fresh:?}"));
        }
    }
}

fn assert_fit_many_matches(run: &Ms1RunIndex, jobs: &[FitJob], p: &EvidenceParams) -> usize {
    let many = fit_many(run, jobs, p);
    assert_eq!(many.len(), jobs.len());
    let mut scratch = FitScratch::default();
    let mut fitted = 0;
    for (k, (job, got)) in jobs.iter().zip(&many).enumerate() {
        let want = fit_on_grid(
            run,
            &job.probe,
            job.charge,
            job.envelope,
            p,
            job.widths,
            &mut scratch,
        );
        let same =
            |a: &Option<Evidence>, b: &Option<Evidence>| format!("{a:?}") == format!("{b:?}");
        assert!(same(got, &want), "job {k}: {job:?}\n{got:?}\n{want:?}");
        if job.widths == &WIDTH_GRID[..] {
            let reference = fit(run, &job.probe, job.charge, job.envelope, p);
            assert!(same(got, &reference), "job {k}: {job:?}");
        }
        fitted += usize::from(got.is_some());
    }
    fitted
}

const FIVE_POINT_GRID: [f64; 5] = [0.5, 0.7, 1.0, 1.4, 2.0];

#[test]
fn fit_many_equals_fit() {
    let mut rng = Rng(11);
    let features: Vec<Feature> = (0..120)
        .map(|k| {
            (
                400.0 + 800.0 * rng.uniform(),
                1 + (k % 4) as u8,
                ENV,
                30.0 + 540.0 * rng.uniform(),
                2.0 + 6.0 * rng.uniform(),
                AMP * (0.01 + rng.uniform()),
            )
        })
        .collect();
    let run = synth_run(&features, NOISE, 1.0);
    let end = run.rt_range().unwrap().1;
    let envelopes = [ENV, [0.7, 0.2, 0.08, 0.02], [1.0, 0.0, 0.0, 0.0]];
    let grids: [&[f64]; 2] = [&WIDTH_GRID, &FIVE_POINT_GRID];
    let probes: Vec<(Probe, u8)> = (0..2000)
        .map(|k| {
            let f = &features[k % features.len()];
            let charge = if k % 5 == 0 {
                1 + (rng.next_u64() % 4) as u8
            } else {
                f.1
            };
            let mut pr = probe(
                f.0,
                f.3 + 20.0 * (rng.uniform() - 0.5),
                1.0 + 10.0 * rng.uniform(),
                1.0 + 8.0 * rng.uniform(),
            );
            match k % 4 {
                0 => {}
                1 => {
                    let n = [-11i8, -7, 7, 11][(k / 4) % 4];
                    pr.kind = Hypothesis::MassShiftDecoy(n);
                    pr.mz += n as f64 * PEPTIDE_MASS_SLOPE / charge as f64;
                }
                2 => {
                    pr.kind = Hypothesis::RtShiftDecoy;
                    pr.t_hat = end * rng.uniform();
                }
                _ => pr.mz = 300.0 + 1200.0 * rng.uniform(),
            }
            if k % 50 == 0 {
                pr.t_hat = [-5.0, 0.0, end, end + 5.0][(k / 50) % 4];
            }
            if k % 97 == 0 {
                pr.s = 0.2;
                pr.w_prior = 0.2;
            }
            (pr, charge)
        })
        .collect();
    let jobs: Vec<FitJob> = probes
        .iter()
        .enumerate()
        .map(|(k, &(probe, charge))| FitJob {
            probe,
            charge,
            envelope: &envelopes[k % envelopes.len()],
            widths: grids[(k / 3) % 2],
            mobility: None,
        })
        .collect();
    let narrow = EvidenceParams {
        tol_ppm: 25.0,
        n_isotopes: 2,
    };
    for p in [params(), narrow] {
        let fitted = assert_fit_many_matches(&run, &jobs, &p);
        eprintln!("fit_many == fit: {fitted} of {} fitted", jobs.len());
        assert!(fitted > 1900 && fitted < jobs.len());
    }
}

#[test]
fn fit_many_equals_fit_with_overlapping_windows_at_run_edges() {
    let mut scans = synth_scans(&[(MZ, 2, ENV, 20.0, 3.0, AMP)], NOISE, 1.0, 4);
    scans.truncate(40);
    for s in scans.iter_mut() {
        s.peaks.push((MZ * (1.0 - 4e-6), 7.0e5));
        s.peaks.push((MZ * (1.0 + 4e-6), 7.0e5));
    }
    scans[25].peaks.push((f64::NAN, 1.0e6));
    let run = Ms1RunIndex::new(scans);
    let end = run.rt_range().unwrap().1;
    let mut probes = Vec::new();
    for &t_hat in &[
        -30.0,
        -3.0,
        0.0,
        1.5,
        20.0,
        20.0,
        20.0,
        end - 1.0,
        end,
        end + 3.0,
        1e9,
    ] {
        for &mz in &[
            MZ,
            MZ,
            MZ * (1.0 + 2e-6),
            MZ + ISOTOPE / 2.0,
            MZ - ISOTOPE / 2.0,
        ] {
            for &(s, w) in &[(6.0, 5.0), (0.5, 0.5), (0.25, 0.25), (20.0, 2.0)] {
                probes.push(probe(mz, t_hat, s, w));
            }
        }
    }
    probes.push(probe(f64::NAN, 20.0, 6.0, 5.0));
    probes.push(probe(MZ, f64::NAN, 6.0, 5.0));
    probes.push(probe(f64::INFINITY, 20.0, 6.0, 5.0));
    let jobs: Vec<FitJob> = probes
        .iter()
        .enumerate()
        .map(|(k, &probe)| FitJob {
            probe,
            charge: [2, 2, 1, 3][k % 4],
            envelope: &ENV,
            widths: if k % 2 == 0 {
                &WIDTH_GRID[..]
            } else {
                &FIVE_POINT_GRID[..]
            },
            mobility: None,
        })
        .collect();
    for n_isotopes in 1..=MAX_ISOTOPES {
        let p = EvidenceParams {
            tol_ppm: 10.0,
            n_isotopes,
        };
        let fitted = assert_fit_many_matches(&run, &jobs, &p);
        assert!(fitted > 0 && fitted < jobs.len());
    }
    assert!(fit_many(&run, &[], &params()).is_empty());
}

#[test]
#[ignore = "timing"]
fn timing_fit_many_50k_probes_on_3000_scan_run() {
    use rayon::prelude::*;

    let mut rng = Rng(13);
    let features: Vec<Feature> = (0..600)
        .map(|k| {
            (
                400.0 + 800.0 * rng.uniform(),
                2 + (k % 2) as u8,
                ENV,
                60.0 + 2800.0 * rng.uniform(),
                5.0,
                AMP * (0.01 + rng.uniform()),
            )
        })
        .collect();
    let mut scans = synth_scans(&features, NOISE, 1.0, 3);
    scans.truncate(3000);
    let run = Ms1RunIndex::new(scans);
    let jobs: Vec<FitJob> = (0..50_000)
        .map(|k| {
            let f = &features[k % features.len()];
            let mz = if k % 2 == 0 {
                f.0
            } else {
                400.0 + 800.0 * rng.uniform()
            };
            FitJob {
                probe: probe(mz, f.3 + 4.0 * (rng.uniform() - 0.5), 6.0, 5.0),
                charge: f.1,
                envelope: &ENV,
                widths: &WIDTH_GRID,
                mobility: None,
            }
        })
        .collect();
    let p = params();
    let per_probe = || -> Vec<Option<Evidence>> {
        jobs.par_iter()
            .map_init(FitScratch::default, |scratch, j| {
                fit_with(&run, &j.probe, j.charge, j.envelope, &p, scratch)
            })
            .collect()
    };
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let cpu = || -> f64 {
        std::fs::read_to_string("/proc/self/stat")
            .ok()
            .and_then(|stat| {
                let fields: Vec<&str> = stat.rsplit(')').next()?.split_whitespace().collect();
                let ticks = fields.get(11)?.parse::<f64>().ok()?
                    + fields.get(12)?.parse::<f64>().ok()?;
                Some(ticks / 100.0)
            })
            .unwrap_or(f64::NAN)
    };
    let time = |label: &str, f: &dyn Fn() -> Vec<Option<Evidence>>| {
        let (start, cpu0) = (std::time::Instant::now(), cpu());
        let out = f();
        eprintln!(
            "timing {label}: {:.3} s wall, {:.2} s CPU for {} probes on {} scans ({} centroids)",
            start.elapsed().as_secs_f64(),
            cpu() - cpu0,
            jobs.len(),
            run.len(),
            run.total_peaks()
        );
        out
    };
    let threads = rayon::current_num_threads();
    let a1 = time("per-probe fit, 1 thread", &|| one.install(per_probe));
    let b1 = time("fit_many, 1 thread", &|| {
        one.install(|| fit_many(&run, &jobs, &p))
    });
    let a = time(&format!("per-probe fit, {threads} threads"), &per_probe);
    let b = time(&format!("fit_many, {threads} threads"), &|| {
        fit_many(&run, &jobs, &p)
    });
    assert_eq!(a1, b1);
    assert_eq!(a, b);
    assert_eq!(a, a1);
}

#[test]
fn timing_10k_probes_on_3000_scan_run() {
    let mut rng = Rng(7);
    let features: Vec<Feature> = (0..300)
        .map(|_| {
            (
                400.0 + 800.0 * rng.uniform(),
                2,
                ENV,
                60.0 + 2800.0 * rng.uniform(),
                5.0,
                AMP * (0.01 + rng.uniform()),
            )
        })
        .chain([(900.0, 2, ENV, 2900.0, 5.0, 0.0)])
        .collect();
    let mut scans = synth_scans(&features, NOISE, 1.0, 3);
    scans.truncate(3000);
    let run = Ms1RunIndex::new(scans);
    let probes: Vec<Probe> = (0..10_000)
        .map(|k| {
            let f = &features[k % features.len()];
            let mz = if k % 2 == 0 {
                f.0
            } else {
                400.0 + 800.0 * rng.uniform()
            };
            probe(mz, f.3 + 4.0 * (rng.uniform() - 0.5), 6.0, 5.0)
        })
        .collect();
    let start = std::time::Instant::now();
    let mut fitted = 0usize;
    for pr in &probes {
        if fit(&run, pr, 2, &ENV, &params()).is_some() {
            fitted += 1;
        }
    }
    eprintln!(
        "timing: {} probes ({} fitted) on {} scans ({} centroids) in {:.3} s",
        probes.len(),
        fitted,
        run.len(),
        run.total_peaks(),
        start.elapsed().as_secs_f64()
    );
    assert_eq!(fitted, probes.len());
}

fn coeluting_by_mobility() -> (Vec<Ms1Scan>, Vec<Ms1Scan>) {
    let target: Feature = (MZ, 2, ENV, 300.0, 6.0, AMP);
    let other: Feature = (MZ + 1e-4, 2, ENV, 300.0, 6.0, 5.0 * AMP);
    let a = synth_scans(&[target], NOISE, 1.0, 11);
    let b = synth_scans(&[other], NOISE, 1.0, 12);
    assert_eq!(a.len(), b.len());
    let resolved = a
        .iter()
        .zip(&b)
        .map(|(x, y)| {
            let mut rows: Vec<(f64, f32, f32)> = x
                .peaks
                .iter()
                .map(|p| (p.0, p.1, 0.90))
                .chain(y.peaks.iter().map(|p| (p.0, p.1, 1.05)))
                .collect();
            rows.sort_by(|p, q| p.0.total_cmp(&q.0));
            Ms1Scan {
                rt: x.rt,
                peaks: rows.iter().map(|r| (r.0, r.1)).collect(),
                mobility: rows.iter().map(|r| r.2).collect(),
                ..Ms1Scan::default()
            }
        })
        .collect();
    (a, resolved)
}

#[test]
fn a_mobility_gate_keeps_a_coeluting_interferer_out_of_the_traces() {
    let (alone, resolved) = coeluting_by_mobility();
    let alone = Ms1RunIndex::new(alone);
    let run = Ms1RunIndex::new(resolved);
    let pr = probe(MZ, 300.0, 6.0, 6.0);
    let p = params();
    let gate = MobilityGate { lo: 0.88, hi: 0.92 };
    let mut scratch = FitScratch::default();
    let gated = fit_on_grid_gated(
        &run,
        &pr,
        2,
        &ENV,
        &p,
        &WIDTH_GRID,
        Some(gate),
        &mut scratch,
    )
    .unwrap();
    let ungated = fit_on_grid(&run, &pr, 2, &ENV, &p, &WIDTH_GRID, &mut scratch).unwrap();
    let reference = fit_on_grid(&alone, &pr, 2, &ENV, &p, &WIDTH_GRID, &mut scratch).unwrap();
    assert!(
        (gated.amplitude / reference.amplitude - 1.0).abs() < 0.03,
        "{gated:?} {reference:?}"
    );
    assert!((gated.amplitude / AMP - 1.0).abs() < 0.05, "{gated:?}");
    assert!(ungated.amplitude > 5.0 * AMP, "{ungated:?}");
    let jobs = [
        FitJob {
            probe: pr,
            charge: 2,
            envelope: &ENV,
            widths: &WIDTH_GRID,
            mobility: Some(gate),
        },
        FitJob {
            probe: pr,
            charge: 2,
            envelope: &ENV,
            widths: &WIDTH_GRID,
            mobility: None,
        },
    ];
    let swept = fit_many(&run, &jobs, &p);
    assert_eq!(swept[0].as_ref(), Some(&gated));
    assert_eq!(swept[1].as_ref(), Some(&ungated));
}

#[test]
fn a_mobility_gate_is_ignored_on_a_run_without_mobility() {
    let (alone, _) = coeluting_by_mobility();
    let run = Ms1RunIndex::new(alone);
    let pr = probe(MZ, 300.0, 6.0, 6.0);
    let p = params();
    let gate = MobilityGate { lo: 0.0, hi: 0.1 };
    let mut scratch = FitScratch::default();
    let gated = fit_on_grid_gated(
        &run,
        &pr,
        2,
        &ENV,
        &p,
        &WIDTH_GRID,
        Some(gate),
        &mut scratch,
    );
    let plain = fit_on_grid(&run, &pr, 2, &ENV, &p, &WIDTH_GRID, &mut scratch);
    assert_eq!(gated, plain);
    let job = FitJob {
        probe: pr,
        charge: 2,
        envelope: &ENV,
        widths: &WIDTH_GRID,
        mobility: Some(gate),
    };
    assert_eq!(fit_many(&run, &[job], &p)[0], plain);
}

#[test]
fn trace_area_matches_the_template_on_a_gaussian_peak() {
    let w = 5.0;
    let truth = AMP * w * (2.0 * std::f64::consts::PI).sqrt();
    for (k, dt) in [0.5, 1.0, 2.0].into_iter().enumerate() {
        let run = Ms1RunIndex::new(synth_scans(
            &[(MZ, 2, ENV, 300.0, w, AMP)],
            NOISE,
            dt,
            40 + k as u64,
        ));
        let e = fit(&run, &probe(MZ, 301.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
        eprintln!(
            "gaussian dt {dt}: trace/template = {:.4}, trace/truth = {:.4}, se/trace = {:.5}",
            e.trace_quantity / e.quantity,
            e.trace_quantity / truth,
            e.trace_quantity_se / e.trace_quantity
        );
        assert!((e.trace_quantity / e.quantity - 1.0).abs() < 0.02);
        assert!((e.trace_quantity / truth - 1.0).abs() < 0.03);
        assert!(e.trace_quantity_se > 0.0 && e.trace_quantity_se < 0.02 * e.trace_quantity);
    }
}

#[test]
fn trace_area_follows_a_tailing_peak() {
    let apex = 300.0;
    let parts = [(apex, 5.0, AMP), (apex + 9.0, 8.0, 0.35 * AMP)];
    let run = synth_run(&parts.map(|(t, w, a)| (MZ, 2, ENV, t, w, a)), NOISE, 1.0);
    let e = fit(&run, &probe(MZ, apex, 6.0, 5.0), 2, &ENV, &params()).unwrap();
    let step = 1e-3;
    let truth: f64 = (0..600_000)
        .map(|k| {
            let t = k as f64 * step;
            parts
                .iter()
                .map(|&(c, w, a)| a * (-(t - c) * (t - c) / (2.0 * w * w)).exp())
                .fold(0.0, f64::max)
                * step
        })
        .sum();
    eprintln!(
        "tailing: template/truth = {:.4}, trace/truth = {:.4} (t0 {:.1}, w {:.2})",
        e.quantity / truth,
        e.trace_quantity / truth,
        e.t0,
        e.w
    );
    assert!((e.quantity / truth - 1.0).abs() > 0.1);
    assert!((e.trace_quantity / truth - 1.0).abs() < 0.08);
}

#[test]
fn trace_area_of_an_absent_analyte_averages_zero() {
    let p = probe(MZ, 300.0, 6.0, 5.0);
    let (mut z_sum, mut z2_sum, mut template_sum) = (0.0, 0.0, 0.0);
    let n = 60;
    for seed in 0..n {
        let run = Ms1RunIndex::new(synth_scans(
            &[(MZ, 2, ENV, 300.0, 5.0, 0.0)],
            NOISE,
            1.0,
            100 + seed,
        ));
        let e = fit(&run, &p, 2, &ENV, &params()).unwrap();
        assert!(e.trace_quantity.is_finite() && e.trace_quantity_se > 0.0);
        let z = e.trace_quantity / e.trace_quantity_se;
        z_sum += z;
        z2_sum += z * z;
        template_sum += e.quantity / e.trace_quantity_se;
    }
    let mean = z_sum / n as f64;
    let sd = (z2_sum / n as f64 - mean * mean).sqrt();
    let template = template_sum / n as f64;
    eprintln!(
        "absent: trace z mean {mean:.3}, sd {sd:.3}; template in trace SEs {template:.3}"
    );
    assert!(mean.abs() < 0.5, "{mean}");
    assert!((0.5..2.0).contains(&sd), "{sd}");
    assert!(template > mean);
}

#[test]
fn trace_area_without_scans_in_range_is_zero() {
    let rt = [0.0, 10.0, 20.0];
    let y: [f32; 3] = [1.0, 5.0, 1.0];
    let traces: [&[f32]; 1] = [&y];
    let pi = [1.0, 0.0, 0.0, 0.0];
    let sigma = [1.0; MAX_ISOTOPES];
    let g = [0.5, 1.0, 0.5];
    let mut sort = Vec::new();
    assert_eq!(
        trace_area(
            &rt,
            &traces,
            &pi,
            &sigma,
            0b1,
            4.0,
            (9.0, 11.0),
            &g,
            &mut sort
        ),
        (0.0, 0.0)
    );
    let (q, se) = trace_area(
        &rt,
        &traces,
        &pi,
        &sigma,
        0b1,
        4.0,
        (10.0, 20.0),
        &g,
        &mut sort,
    );
    let base = 7.0 / 3.0 - 4.0 * 2.0 / 3.0;
    assert!((q - 10.0 * (0.5 * (5.0 - base) + 0.5 * (1.0 - base))).abs() < 1e-9);
    assert!((se - (50.0 + 100.0 / 3.0f64).sqrt()).abs() < 1e-9, "{se}");
}

#[test]
fn trace_baseline_comes_from_the_lower_side_when_the_sides_disagree() {
    let rt: Vec<f64> = (0..25).map(f64::from).collect();
    let mut y = vec![1.0f32; 25];
    for v in &mut y[10..15] {
        *v = 101.0;
    }
    for v in &mut y[15..] {
        *v = 50.0;
    }
    let traces: [&[f32]; 1] = [&y];
    let (pi, sigma) = ([1.0, 0.0, 0.0, 0.0], [1.0; MAX_ISOTOPES]);
    let g = vec![0.0; 25];
    let mut sort = Vec::new();
    let (q, _) = trace_area(
        &rt,
        &traces,
        &pi,
        &sigma,
        0b1,
        0.0,
        (10.0, 14.0),
        &g,
        &mut sort,
    );
    assert!((q - 400.0).abs() < 1e-9, "{q}");
    y[15..].fill(1.0);
    y[15] = 1.5;
    let traces: [&[f32]; 1] = [&y];
    let (q, _) = trace_area(
        &rt,
        &traces,
        &pi,
        &sigma,
        0b1,
        0.0,
        (10.0, 14.0),
        &g,
        &mut sort,
    );
    assert!((q - 400.0).abs() < 1e-9, "{q}");
}

#[test]
fn trace_baseline_pools_agreeing_sides_with_different_medians() {
    let rt: Vec<f64> = (0..25).map(f64::from).collect();
    let mut y = vec![10.0f32; 25];
    for v in &mut y[10..15] {
        *v = 101.0;
    }
    for v in &mut y[15..] {
        *v = 11.0;
    }
    let traces: [&[f32]; 1] = [&y];
    let (pi, sigma) = ([1.0, 0.0, 0.0, 0.0], [1.0; MAX_ISOTOPES]);
    let g = vec![0.0; 25];
    let mut sort = Vec::new();
    let (q, _) = trace_area(
        &rt,
        &traces,
        &pi,
        &sigma,
        0b1,
        0.0,
        (10.0, 14.0),
        &g,
        &mut sort,
    );
    assert!((q - 362.0).abs() < 1e-9, "{q}");
}

#[test]
fn median_of_an_unsorted_slice() {
    assert_eq!(median_f32_as_f64_mut(&mut []).unwrap_or(0.0), 0.0);
    assert_eq!(median_f32_as_f64_mut(&mut [3.0, 1.0, 2.0]), Some(2.0));
    assert_eq!(median_f32_as_f64_mut(&mut [4.0, 1.0, 3.0, 2.0]), Some(2.5));
}
