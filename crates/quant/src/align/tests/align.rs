use super::*;

struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        model::hash::splitmix64_next(&mut self.0)
    }

    fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64 + 1.0) / (1u64 << 53) as f64
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

fn synth_truth(
    n_p: usize,
    warps: &[fn(f64) -> f64],
    noise: f64,
    seed: u64,
) -> (Vec<f64>, Vec<Observation>) {
    let mut rng = Rng(seed);
    let truth: Vec<f64> = (0..n_p).map(|_| 600.0 + 3000.0 * rng.uniform()).collect();
    let mut obs = Vec::with_capacity(n_p * warps.len());
    for (p, &t) in truth.iter().enumerate() {
        for (r, w) in warps.iter().enumerate() {
            obs.push(Observation {
                precursor: p as u32,
                run: r as u16,
                rt: w(t) + noise * rng.normal(),
            });
        }
    }
    (truth, obs)
}

fn synth(n_p: usize, warps: &[fn(f64) -> f64], noise: f64, seed: u64) -> Vec<Observation> {
    synth_truth(n_p, warps, noise, seed).1
}

fn median(v: &[f64]) -> f64 {
    let mut v = v.to_vec();
    median_in_place(&mut v)
}

fn three_warps() -> [fn(f64) -> f64; 3] {
    [
        |t| t,
        |t| 1.1 * t + 30.0,
        |t| t + 60.0 * (t / 3600.0).powi(2),
    ]
}

#[test]
fn recovers_monotone_warps() {
    let warps = three_warps();
    let obs = synth(2000, &warps, 5.0, 1);
    let a = align(&obs, 2000, 3, &AlignParams::default());
    let errs: Vec<f64> = obs
        .iter()
        .map(|o| (a.predict(o.precursor, o.run).unwrap().0 - o.rt).abs())
        .collect();
    let m = median(&errs);
    let held: Vec<f64> = obs
        .iter()
        .map(|o| (a.predict_excluding(o.precursor, o.run).unwrap().0 - o.rt).abs())
        .collect();
    let mh = median(&held);
    eprintln!("recovers_monotone_warps: median |t_hat - rt| = {m:.3}, held out = {mh:.3}");
    let sds: Vec<f64> = a.residual_sd[0]
        .as_ref()
        .unwrap()
        .iter()
        .map(|k| k.1)
        .collect();
    eprintln!(
        "recovers_monotone_warps: run 0 median sigma_r = {:.3} (noise 5)",
        median(&sds)
    );
    assert!(m < 7.5, "{m}");
    assert!(mh < 7.5, "{mh}");
}

#[test]
fn robust_to_ten_percent_outliers() {
    let warps = three_warps();
    let (truth, mut obs) = synth_truth(2000, &warps, 5.0, 2);
    let mut rng = Rng(99);
    let mut is_outlier = vec![false; obs.len()];
    for (i, o) in obs.iter_mut().enumerate() {
        if rng.uniform() <= 0.10 {
            o.rt = 4000.0 * rng.uniform();
            is_outlier[i] = true;
        }
    }
    let a = align(&obs, 2000, 3, &AlignParams::default());
    let errs: Vec<f64> = obs
        .iter()
        .map(|o| (a.predict(o.precursor, o.run).unwrap().0 - o.rt).abs())
        .collect();
    let m = median(&errs);
    let errs_truth: Vec<f64> = obs
        .iter()
        .zip(&is_outlier)
        .filter(|(_, &out)| !out)
        .map(|(o, _)| {
            let t = truth[o.precursor as usize];
            (a.predict(o.precursor, o.run).unwrap().0 - warps[o.run as usize](t)).abs()
        })
        .collect();
    let mt = median(&errs_truth);
    eprintln!(
        "robust_to_ten_percent_outliers: median |t_hat - rt| = {m:.3}, \
         median |t_hat - truth| (inliers) = {mt:.3}"
    );
    assert!(m < 7.5, "{m}");
    assert!(mt < 5.0, "{mt}");
}

#[test]
fn maps_are_monotone_and_invertible() {
    let warps = three_warps();
    let mut obs = synth(2000, &warps, 5.0, 3);
    let mut rng = Rng(7);
    for o in obs.iter_mut() {
        if rng.uniform() <= 0.10 {
            o.rt = 4000.0 * rng.uniform();
        }
    }
    let a = align(&obs, 2000, 3, &AlignParams::default());
    let mut n_aligned = 0;
    for m in a.maps.iter().flatten() {
        n_aligned += 1;
        let mut prev = f64::NEG_INFINITY;
        for i in 0..=5000 {
            let t = -500.0 + i as f64;
            let y = m.apply(t);
            assert!(y.is_finite());
            assert!(y >= prev, "apply decreases at {t}: {prev} -> {y}");
            prev = y;
            let back = m.invert(y);
            assert!((back - t).abs() < 1e-6, "invert(apply({t})) = {back}");
        }
    }
    assert_eq!(n_aligned, 3);
}

fn bits(v: Option<(f64, f64)>) -> Option<(u64, u64)> {
    v.map(|(a, b)| (a.to_bits(), b.to_bits()))
}

#[test]
fn rebuilt_alignment_predicts_bit_for_bit() {
    let warps = three_warps();
    let mut obs = synth(600, &warps, 5.0, 8);
    obs.extend((0..10u32).map(|p| Observation {
        precursor: p,
        run: 3,
        rt: 50.0 + p as f64,
    }));
    obs.retain(|o| !(o.precursor >= 590 && o.run > 0));
    let a = align(&obs, 610, 4, &AlignParams::default());
    assert!(a.maps[3].is_none());
    let runs = a.runs();
    assert!(runs[3].observations.is_empty());
    let b = Alignment::from_runs(runs.clone(), 610).unwrap();
    assert_eq!(b.runs(), runs);
    assert_eq!(a.n_runs_seen, b.n_runs_seen);
    assert_eq!(a.anchors, b.anchors);
    for p in 0..610u32 {
        let i = p as usize;
        assert_eq!(a.consensus[i].to_bits(), b.consensus[i].to_bits(), "{p}");
        assert_eq!(a.tau[i].to_bits(), b.tau[i].to_bits(), "{p}");
        for r in 0..4u16 {
            assert_eq!(bits(a.predict(p, r)), bits(b.predict(p, r)), "{p} {r}");
            assert_eq!(
                bits(a.predict_excluding(p, r)),
                bits(b.predict_excluding(p, r)),
                "{p} {r}"
            );
        }
    }
}

#[test]
fn added_observations_are_predicted_without_moving_the_fit() {
    let warps = three_warps();
    let (truth, mut obs) = synth_truth(620, &warps, 5.0, 9);
    let extra: Vec<Observation> = obs
        .iter()
        .filter(|o| o.precursor >= 600 && o.run < 2)
        .map(|o| Observation {
            rt: o.rt + 400.0,
            ..*o
        })
        .collect();
    obs.retain(|o| o.precursor < 600);
    let fitted = align(&obs, 620, 3, &AlignParams::default());
    let mut a = fitted.clone();
    a.add_observations(&extra).unwrap();
    for r in 0..3 {
        let (x, y) = (fitted.maps[r].as_ref(), a.maps[r].as_ref());
        assert_eq!(x.map(|m| m.knots()), y.map(|m| m.knots()));
    }
    assert_eq!(fitted.anchors, a.anchors);
    for p in 0..600u32 {
        for r in 0..3u16 {
            assert_eq!(bits(fitted.predict(p, r)), bits(a.predict(p, r)));
        }
    }
    let map = |r: usize| a.maps[r].as_ref().unwrap();
    for p in 600..620u32 {
        assert!(fitted.predict(p, 2).is_none());
        assert_eq!(a.n_runs_seen[p as usize], 2);
        let rt = |r: u16| {
            extra
                .iter()
                .find(|o| o.precursor == p && o.run == r)
                .unwrap()
                .rt
        };
        let tc = 0.5 * (map(0).invert(rt(0)) + map(1).invert(rt(1)));
        let (t, s) = a.predict(p, 2).unwrap();
        assert!((t - map(2).apply(tc)).abs() < 1e-6, "{p}");
        assert!(t > warps[2](truth[p as usize]) + 200.0, "{p}");
        assert!(s.is_finite() && s > 0.0);
        assert!(a.predict_excluding(p, 0).is_some());
    }
    let b = Alignment::from_runs(a.runs(), 620).unwrap();
    assert_eq!(bits(a.predict(605, 2)), bits(b.predict(605, 2)));
    assert!(a.add_observations(&extra[..1]).is_err());
}

#[test]
fn from_runs_rejects_inconsistent_state() {
    let a = align(
        &synth(300, &three_warps(), 5.0, 2),
        300,
        3,
        &AlignParams::default(),
    );
    let mut runs = a.runs();
    let first = runs[0].observations[0];
    runs[0].observations.push(first);
    let e = Alignment::from_runs(runs, 300).unwrap_err();
    assert!(e.contains("observed twice"), "{e}");
    let e = Alignment::from_runs(a.runs(), 10).unwrap_err();
    assert!(e.contains("outside"), "{e}");
    let mut runs = a.runs();
    runs[1].map = None;
    assert!(Alignment::from_runs(runs, 300).is_err());
    assert!(RunMap::from_parts(vec![], vec![], 1.0, 1.0).is_err());
    assert!(RunMap::from_parts(vec![1.0], vec![2.0], 0.0, 1.0).is_err());
    assert!(RunMap::from_parts(vec![1.0], vec![f64::NAN], 1.0, 1.0).is_err());
}

#[test]
fn run_with_few_anchors_is_unaligned() {
    let warps = three_warps();
    let mut obs = synth(2000, &warps, 5.0, 4);
    for p in 0..20u32 {
        obs.push(Observation {
            precursor: p,
            run: 3,
            rt: 1000.0 + 10.0 * p as f64,
        });
    }
    let a = align(&obs, 2000, 4, &AlignParams::default());
    assert!(a.maps[3].is_none());
    assert!(a.residual_sd[3].is_none());
    assert_eq!(a.anchors[3], 20);
    for p in 0..2000u32 {
        assert!(a.predict(p, 3).is_none());
        assert!(a.predict_excluding(p, 3).is_none());
    }
    assert_eq!(a.n_runs_seen[0], 3);
    let errs: Vec<f64> = obs
        .iter()
        .filter(|o| o.run < 3)
        .map(|o| (a.predict(o.precursor, o.run).unwrap().0 - o.rt).abs())
        .collect();
    let m = median(&errs);
    assert!(m < 7.5, "{m}");
}

#[test]
fn run_with_narrow_rt_range_is_unaligned() {
    let warps = three_warps();
    let (truth, mut obs) = synth_truth(2000, &warps, 5.0, 5);
    let mut added = 0;
    for (p, &t) in truth.iter().enumerate() {
        if t < 800.0 && added < 200 {
            obs.push(Observation {
                precursor: p as u32,
                run: 3,
                rt: t,
            });
            added += 1;
        }
    }
    assert!(added >= 100);
    let a = align(&obs, 2000, 4, &AlignParams::default());
    assert!(a.maps[3].is_none());
    assert!(a.coverage[3] < 0.30, "{}", a.coverage[3]);
    assert!(a.maps[..3].iter().all(Option::is_some));
}

#[test]
fn single_run_precursor_gets_consensus_without_anchoring() {
    let warps = three_warps();
    let mut obs = synth(2000, &warps, 5.0, 6);
    let base = align(&obs, 2001, 3, &AlignParams::default());
    obs.push(Observation {
        precursor: 2000,
        run: 1,
        rt: 2000.0,
    });
    let a = align(&obs, 2001, 3, &AlignParams::default());
    assert!(base.consensus[2000].is_nan());
    assert!(base.predict(2000, 0).is_none());
    assert!(a.consensus[2000].is_finite());
    assert_eq!(a.n_runs_seen[2000], 1);
    assert_eq!(a.maps, base.maps);
    assert_eq!(a.residual_sd, base.residual_sd);
    for p in 0..2000 {
        assert_eq!(a.consensus[p].to_bits(), base.consensus[p].to_bits());
    }
    let (t1, _) = a.predict(2000, 1).unwrap();
    assert!((t1 - 2000.0).abs() < 1e-6, "{t1}");
    assert!(a.predict(2000, 0).is_some());
    assert!(a.predict_excluding(2000, 1).is_none());
}

#[test]
fn transfer_sd_shrinks_with_more_donors() {
    let warps: [fn(f64) -> f64; 7] = [
        |t| t,
        |t| 1.05 * t + 20.0,
        |t| 0.95 * t - 10.0,
        |t| t + 40.0,
        |t| t + 60.0 * (t / 3600.0).powi(2),
        |t| 1.02 * t,
        |t| t - 30.0,
    ];
    let mut obs = synth(1500, &warps, 5.0, 8);
    for (r, w) in warps.iter().enumerate().take(3).skip(1) {
        obs.push(Observation {
            precursor: 1500,
            run: r as u16,
            rt: w(2000.0),
        });
    }
    for (r, w) in warps.iter().enumerate().skip(1) {
        obs.push(Observation {
            precursor: 1501,
            run: r as u16,
            rt: w(2000.0),
        });
    }
    let a = align(&obs, 1502, 7, &AlignParams::default());
    assert_eq!(a.n_runs_seen[1500], 2);
    assert_eq!(a.n_runs_seen[1501], 6);
    let (t2, s2) = a.predict(1500, 0).unwrap();
    let (t6, s6) = a.predict(1501, 0).unwrap();
    eprintln!(
        "transfer_sd_shrinks_with_more_donors: s(2 donors) = {s2:.3}, s(6 donors) = {s6:.3}"
    );
    assert!(
        (t2 - 2000.0).abs() < 5.0 && (t6 - 2000.0).abs() < 5.0,
        "{t2} {t6}"
    );
    assert!(s2 > s6, "{s2} <= {s6}");
}

fn noisy_singleton_donor() -> (Alignment, f64) {
    let mut rng = Rng(21);
    let mut obs = Vec::new();
    for p in 0..300u32 {
        let t = 600.0 + 3000.0 * rng.uniform();
        for r in 0..4u16 {
            let jitter = if r == 0 {
                [-40.0, -20.0, 20.0, 40.0][(rng.next_u64() % 4) as usize]
            } else {
                0.0
            };
            obs.push(Observation {
                precursor: p,
                run: r,
                rt: t + jitter,
            });
        }
    }
    for (precursor, run, rt) in [(300, 0, 2047.0), (301, 0, 2047.0), (301, 1, 2000.0)] {
        obs.push(Observation { precursor, run, rt });
    }
    let a = align(&obs, 302, 4, &AlignParams::default());
    let sds: Vec<f64> = a.residual_sd[0]
        .as_ref()
        .unwrap()
        .iter()
        .map(|k| k.1)
        .collect();
    let sd0 = median(&sds);
    (a, sd0)
}

#[test]
fn transfer_window_includes_the_donor_uncertainty() {
    use crate::decoys::{Hypothesis, Probe};
    use crate::evidence::{fit, EvidenceParams, Prior};
    use crate::ms1_index::{Ms1RunIndex, Ms1Scan};
    use model::mass::ISOTOPE;

    let (a, sd0) = noisy_singleton_donor();
    let (t_hat, s) = a.predict(300, 1).unwrap();
    eprintln!(
        "transfer_window_includes_the_donor_uncertainty: donor sigma = {sd0:.2}, \
         t_hat = {t_hat:.2}, s = {s:.3}"
    );
    let probe = Probe {
        kind: Hypothesis::Target,
        mz: 600.0,
        t_hat,
        s,
        w_prior: 5.0,
    };
    let prior = Prior::of(&probe);

    let env = [0.55f32, 0.30, 0.11, 0.04];
    let scans: Vec<Ms1Scan> = (1700..2300)
        .map(|k| {
            let t = k as f64;
            let g = (-(t - 2000.0).powi(2) / (2.0 * 25.0)).exp();
            let peaks = (0..4)
                .map(|i| {
                    (
                        600.0 + i as f64 * ISOTOPE / 2.0,
                        (100.0 + 1e6 * env[i] as f64 * g) as f32,
                    )
                })
                .collect();
            Ms1Scan {
                rt: t,
                peaks,
                ..Ms1Scan::default()
            }
        })
        .collect();
    let run = Ms1RunIndex::new(scans);
    let e = fit(&run, &probe, 2, &env, &EvidenceParams::default()).unwrap();
    eprintln!(
        "transfer_window_includes_the_donor_uncertainty: quantity = {:.0}, t0 = {:.2}",
        e.quantity, e.t0
    );
    assert!(e.quantity > 0.0, "{e:?}");
    assert!((e.t0 - 2000.0).abs() <= 1.0, "{e:?}");
    assert!((2000.0 - t_hat).abs() <= 3.0 * s, "{t_hat} {s}");
    assert!(t_hat - prior.half_window() < 2000.0);
    assert!(sd0 > 20.0, "{sd0}");
    assert!(s >= 0.8 * sd0, "s = {s}, donor sigma = {sd0}");
}

#[test]
fn excluded_prediction_carries_the_donor_uncertainty() {
    let (a, sd0) = noisy_singleton_donor();
    let (_, s) = a.predict_excluding(0, 1).unwrap();
    eprintln!(
        "excluded_prediction_carries_the_donor_uncertainty: s = {s:.3}, \
         donor sigma = {sd0:.2}"
    );
    assert!(s > 0.2 * sd0 && s < 0.6 * sd0, "{s}");
    let (_, s) = a.predict_excluding(301, 1).unwrap();
    assert!(s >= 0.8 * sd0, "{s}");
    assert_eq!(a.predict_excluding(300, 1), a.predict(300, 1));
}

fn two_pairs() -> Alignment {
    let mut rng = Rng(31);
    let mut obs = Vec::new();
    for p in 0..200u32 {
        let t = 600.0 + 3000.0 * rng.uniform();
        let (runs, shift) = if p < 100 {
            (0..2u16, 0.0)
        } else {
            (2..4u16, 120.0)
        };
        for r in runs {
            obs.push(Observation {
                precursor: p,
                run: r,
                rt: t + shift + 0.5 * rng.normal(),
            });
        }
    }
    obs.push(Observation {
        precursor: 200,
        run: 0,
        rt: 1500.0,
    });
    obs.push(Observation {
        precursor: 201,
        run: 2,
        rt: 1620.0,
    });
    align(&obs, 202, 4, &AlignParams::default())
}

#[test]
fn no_transfer_between_disconnected_runs() {
    let a = two_pairs();
    assert!(a.maps.iter().all(Option::is_some));
    assert_eq!(a.component(0), a.component(1));
    assert_eq!(a.component(2), a.component(3));
    assert_ne!(a.component(0), a.component(2));
    if let Some((t, s)) = a.predict(200, 2) {
        panic!("precursor 200 transferred across components: {t:.2} +- {s:.3}");
    }
    for (p, r) in [(200u32, 3u16), (201, 0), (201, 1), (5, 2), (150, 0)] {
        assert!(a.predict(p, r).is_none(), "{p} into {r}");
        assert!(a.predict_excluding(p, r).is_none(), "{p} into {r}");
        assert!(a.blocked_by_component(p, r, None), "{p} into {r}");
    }
    let (t, _) = a.predict(200, 1).unwrap();
    assert!((t - 1500.0).abs() < 3.0, "{t}");
    let (t, _) = a.predict(201, 3).unwrap();
    assert!((t - 1620.0).abs() < 3.0, "{t}");
    assert!(a.predict_excluding(5, 1).is_some());
    assert!(a.predict_excluding(150, 3).is_some());
    assert!(!a.blocked_by_component(200, 1, None));
    assert!(!a.blocked_by_component(5, 1, Some(1)));
}

#[test]
fn precursor_seen_in_two_components_uses_the_destination_component() {
    let mut rng = Rng(31);
    let mut obs = Vec::new();
    for p in 0..1200u32 {
        let t = 600.0 + 3000.0 * rng.uniform();
        let (runs, shift) = if p < 600 {
            (0..2u16, 0.0)
        } else {
            (2..4u16, 120.0)
        };
        for r in runs {
            obs.push(Observation {
                precursor: p,
                run: r,
                rt: t + shift + 0.5 * rng.normal(),
            });
        }
    }
    obs.push(Observation {
        precursor: 1200,
        run: 0,
        rt: 1500.0,
    });
    obs.push(Observation {
        precursor: 1200,
        run: 2,
        rt: 1620.0,
    });
    let a = align(&obs, 1201, 4, &AlignParams::default());
    assert_eq!(a.n_components(), 2);
    let (t1, _) = a.predict(1200, 1).unwrap();
    let (t3, _) = a.predict(1200, 3).unwrap();
    eprintln!(
        "precursor_seen_in_two_components_uses_the_destination_component: \
         run 1 {t1:.2} (1500), run 3 {t3:.2} (1620)"
    );
    assert!((t1 - 1500.0).abs() < 3.0, "{t1}");
    assert!((t3 - 1620.0).abs() < 3.0, "{t3}");
    assert!(a.predict_excluding(1200, 2).is_none());
    assert!(a.blocked_by_component(1200, 2, Some(2)));
}

#[test]
fn each_run_gets_its_own_residual_width_floor() {
    let mut rng = Rng(7);
    let mut obs = Vec::new();
    for p in 0..300u32 {
        let t = 600.0 + 3000.0 * rng.uniform();
        for r in 0..3u16 {
            obs.push(Observation {
                precursor: p,
                run: r,
                rt: t + 0.2 * rng.normal(),
            });
        }
    }
    let ap = AlignParams {
        run_min_spacing_s: vec![0.0, 1.5],
        ..AlignParams::default()
    };
    assert_eq!(
        (ap.floor_s(0), ap.floor_s(1), ap.floor_s(2)),
        (0.5, 1.5, 0.5)
    );
    let a = align(&obs, 300, 3, &ap);
    let widths = |r: usize| a.residual_sd[r].as_ref().unwrap().iter().map(|k| k.1);
    assert!(
        widths(1).all(|w| w >= 1.5),
        "run 1 is floored at its own spacing"
    );
    assert!(
        widths(0).all(|w| (0.5..1.5).contains(&w)),
        "run 0 keeps the default floor"
    );
    assert!(widths(2).all(|w| (0.5..1.5).contains(&w)));
}

#[test]
fn run_sharing_too_few_precursors_with_any_run_is_unaligned() {
    let mut rng = Rng(41);
    let mut obs = Vec::new();
    for p in 0..200u32 {
        let t = 600.0 + 3000.0 * rng.uniform();
        let mut runs: Vec<u16> = if p < 100 { vec![0, 1] } else { vec![2, 3] };
        if p % 100 < 12 {
            runs.push(4);
        }
        for r in runs {
            obs.push(Observation {
                precursor: p,
                run: r,
                rt: t + 0.5 * rng.normal(),
            });
        }
    }
    let ap = AlignParams {
        min_anchors: 20,
        ..AlignParams::default()
    };
    let a = align(&obs, 200, 5, &ap);
    assert_eq!(a.anchors[4], 24);
    assert!(a.maps[4].is_none() && a.component(4).is_none());
    assert!(a.maps[..4].iter().all(Option::is_some));
    assert_eq!(a.n_components(), 2);
    assert!(a.predict(0, 4).is_none());
    assert_eq!(a.n_runs_seen[0], 2);
}

#[test]
fn predict_excluding_ignores_own_observation() {
    let warps = three_warps();
    let (truth, mut obs) = synth_truth(2000, &warps, 1.0, 9);
    let t = 2000.0;
    obs.push(Observation {
        precursor: 2000,
        run: 0,
        rt: warps[0](t) + 300.0,
    });
    obs.push(Observation {
        precursor: 2000,
        run: 1,
        rt: warps[1](t),
    });
    let a = align(&obs, 2001, 3, &AlignParams::default());
    let (own, _) = a.predict(2000, 0).unwrap();
    let (excl, s_excl) = a.predict_excluding(2000, 0).unwrap();
    eprintln!(
        "predict_excluding_ignores_own_observation: predict = {own:.3}, \
         predict_excluding = {excl:.3}, truth = {:.3}",
        warps[0](t)
    );
    assert!((own - warps[0](t)).abs() > 100.0, "{own}");
    assert!((excl - warps[0](t)).abs() < 3.0, "{excl}");
    let (_, s_other) = a.predict(1, 0).unwrap();
    assert!(s_excl < 3.0 * s_other.max(1.0), "{s_excl}");
    let p = 5u32;
    let (e, _) = a.predict_excluding(p, 2).unwrap();
    assert!((e - warps[2](truth[p as usize])).abs() < 5.0, "{e}");
    assert!(a.predict_excluding(2000, 2).is_some());
    assert_eq!(a.predict_excluding(2000, 2), a.predict(2000, 2));
}

#[test]
fn mean_map_is_the_identity_after_the_gauge_fix() {
    let warps = three_warps();
    let obs = synth(2000, &warps, 5.0, 10);
    let a = align(&obs, 2000, 3, &AlignParams::default());
    let maps: Vec<&RunMap> = a.maps.iter().flatten().collect();
    assert_eq!(maps.len(), 3);
    let dev: Vec<f64> = (0..=30)
        .map(|i| {
            let t = 700.0 + 100.0 * i as f64;
            maps.iter().map(|m| m.apply(t)).sum::<f64>() / maps.len() as f64 - t
        })
        .collect();
    let mean = dev.iter().sum::<f64>() / dev.len() as f64;
    let worst = dev.iter().fold(0.0f64, |w, d| w.max(d.abs()));
    eprintln!(
        "mean_map_is_the_identity_after_the_gauge_fix: mean (g(t) - t) = {mean:.3}, \
         max |g(t) - t| = {worst:.3}"
    );
    assert!(mean.abs() < 0.5, "{mean}");
    assert!(worst < 4.0, "{worst}");
}

#[test]
fn sigma_tracks_true_noise() {
    let three = three_warps();
    let six: [fn(f64) -> f64; 6] = [
        |t| t,
        |t| 1.1 * t + 30.0,
        |t| t + 60.0 * (t / 3600.0).powi(2),
        |t| 0.95 * t - 10.0,
        |t| t + 40.0,
        |t| 1.02 * t,
    ];
    for (name, warps) in [("3 runs", &three[..]), ("6 runs", &six[..])] {
        let obs = synth(2000, warps, 5.0, 12);
        let a = align(&obs, 2000, warps.len(), &AlignParams::default());
        for (r, knots) in a.residual_sd.iter().enumerate() {
            let sds: Vec<f64> = knots.as_ref().unwrap().iter().map(|k| k.1).collect();
            let m = median(&sds);
            eprintln!("sigma_tracks_true_noise: {name}, run {r}: median sigma_r = {m:.3}");
            assert!((m - 5.0).abs() < 1.25, "{name} run {r}: {m}");
        }
    }
}

#[test]
fn output_is_independent_of_input_order() {
    let warps = three_warps();
    let obs = synth(2000, &warps, 5.0, 11);
    let a = align(&obs, 2000, 3, &AlignParams::default());
    let mut rev = obs.clone();
    rev.reverse();
    let b = align(&rev, 2000, 3, &AlignParams::default());
    assert_eq!(a.maps, b.maps);
    assert_eq!(a.residual_sd, b.residual_sd);
    assert_eq!(a.n_runs_seen, b.n_runs_seen);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&a.consensus), bits(&b.consensus));
    assert_eq!(bits(&a.tau), bits(&b.tau));
}

#[test]
fn degenerate_inputs_do_not_panic() {
    let p = AlignParams::default();
    let a = align(&[], 0, 0, &p);
    assert!(a.maps.is_empty() && a.consensus.is_empty());
    let a = align(&[], 3, 2, &p);
    assert!(a.maps.iter().all(Option::is_none));
    assert!(a.predict(0, 0).is_none() && a.predict(9, 9).is_none());
    let obs: Vec<Observation> = (0..100)
        .map(|i| Observation {
            precursor: i,
            run: 0,
            rt: 600.0 + i as f64,
        })
        .collect();
    let a = align(&obs, 100, 1, &p);
    assert!(a.maps[0].is_none());
    assert!(a.consensus.iter().all(|c| c.is_nan()));
    let mut obs: Vec<Observation> = (0..200u32)
        .flat_map(|i| {
            (0..2u16).map(move |r| Observation {
                precursor: i,
                run: r,
                rt: 1000.0,
            })
        })
        .collect();
    obs.push(Observation {
        precursor: 200,
        run: 0,
        rt: f64::NAN,
    });
    let a = align(&obs, 201, 2, &p);
    assert!(a.consensus[200].is_nan());
}
