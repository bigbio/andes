//! Tests of labelled-pair quantification.

use super::classify::{
    channel_qvalues, classify_channel, decoy_min_sep, decoy_shifts, learn_rt_shifts, log2_ratio,
    mass_factor, weighted_qvalues, ChannelStatus, ClassifyParams,
};
use super::fit::{fit_group, nn_active_set, nn_enumerate, ChannelEnvelope, FitParams, GroupProbe};
use super::spec::{LabelSpec, Site};
use crate::ms1::Ms1RunIndex;
use model::mass::ISOTOPE;

use super::*;
use crate::ms1::Ms1Scan;

#[test]
fn presets_parse() {
    let s = LabelSpec::parse("silac2").unwrap();
    assert_eq!(s.len(), 2);
    assert!(s.has_unlabelled_channel());
    assert_eq!(s.channels[1].name, "Lys8+Arg10");
    assert_eq!(s.channels[1].qpx_label.as_deref(), Some("SILAC heavy"));
    let s = LabelSpec::parse("silac3").unwrap();
    assert_eq!(s.channels[1].qpx_label.as_deref(), Some("SILAC medium"));
    assert!(s.has_deuterium());
    let s = LabelSpec::parse("[][Lys4][Lys6][Lys8]").unwrap();
    assert!(s.channels.iter().all(|c| c.qpx_label.is_none()));
    let s = LabelSpec::parse("dimethyl3").unwrap();
    assert_eq!(s.len(), 3);
    assert!(!s.has_unlabelled_channel());
    assert_eq!(s.channels[2].qpx_label.as_deref(), Some("DIMETHYL8"));
    assert_eq!(s.sites().len(), 2);
    let s = LabelSpec::parse("[][K+8.014199:6C2N,R+10.008269:6C4N]").unwrap();
    assert_eq!(s.channels[1].sites.len(), 2);
    assert_eq!(s.channels[1].sites[0].heavy_atoms, 8);
    assert_eq!(s.channels[1].sites[0].carbon13, 6);
    assert!(!s.channels[1].sites[0].chemical);
    assert!(s.channels[1].qpx_label.is_some());
    for l in &s.channels[1].sites {
        assert!(!l.name.contains('+'), "{}", l.name);
        assert_eq!(s.channel_of_label(l.site, &l.name), Some(1));
    }
}

#[test]
fn bad_specs_are_refused() {
    assert!(LabelSpec::parse("[Lys8]").is_err());
    assert!(LabelSpec::parse("[Lys8][]").is_err());
    assert!(LabelSpec::parse("[][Lys8][Arg10]").is_err());
    assert!(LabelSpec::parse("[][Lys8,Arg10][Lys4,Arg6]").is_err());
    assert!(LabelSpec::parse("[][Foo]").is_err());
    assert!(LabelSpec::parse("[][K+8.014199:8C]").is_err());
    assert!(LabelSpec::parse("[Nterm+28.0313][Nterm+32.056407:4D]").is_err());
    assert!(LabelSpec::parse("[][K+8.014199:6X]").is_err());
    assert!(LabelSpec::parse("[][K+8.014199:6C6C]").is_err());
}

#[test]
fn custom_labels_carry_their_composition() {
    let custom =
        LabelSpec::parse("[Nterm+28.0313][Nterm+32.056407:chem4D][Nterm+36.07567:chem2C6D]")
            .unwrap();
    let preset = LabelSpec::parse("[Dimethyl0][Dimethyl4][Dimethyl8]").unwrap();
    for (c, p) in custom.channels.iter().zip(&preset.channels) {
        let (a, b) = (&c.sites[0], p.label_at(Site::NTerm).unwrap());
        assert_eq!(
            (
                a.deuterium,
                a.heavy_atoms,
                a.carbon13,
                a.natural_c,
                a.natural_h,
                a.chemical
            ),
            (
                b.deuterium,
                b.heavy_atoms,
                b.carbon13,
                b.natural_c,
                b.natural_h,
                b.chemical
            ),
            "{}",
            c.name
        );
    }
    let s = LabelSpec::parse("[][K+8.014199]").unwrap();
    let l = &s.channels[1].sites[0];
    assert!(l.chemical);
    assert_eq!((l.deuterium, l.heavy_atoms), (0, 0));
}

#[test]
fn channel_delta_counts_sites() {
    let s = LabelSpec::parse("silac2").unwrap();
    let sites = s.sites();
    let k = sites
        .iter()
        .position(|x| *x == Site::Residue(b'K'))
        .unwrap();
    let mut counts = vec![0; sites.len()];
    counts[k] = 2;
    assert!((s.channel_delta(1, &sites, &counts) - 16.028398).abs() < 1e-6);
    assert_eq!(s.channel_delta(0, &sites, &counts), 0.0);
}

#[test]
fn deuterium_enrichment_gives_an_m_minus_1() {
    let light = [0.6, 0.3, 0.08, 0.02];
    let e = ChannelEnvelope::new(&light, 8, 0.96, 0, 1.0);
    let r = e.at(-1) / e.at(0);
    assert!(r > 0.2 && r < 0.6, "{r}");
    let total: f64 = e.0.iter().sum();
    assert!((total - 1.0).abs() < 1e-12);
    let plain = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    assert_eq!(plain.at(-1), 0.0);
    assert!((plain.at(0) - 0.6).abs() < 1e-12);
}

pub(crate) fn synthetic_run(
    species: &[(f64, ChannelEnvelope, f64, f64)],
    w: f64,
    noise: f32,
) -> Ms1RunIndex {
    let scans = (0..240)
        .map(|i| {
            let rt = i as f64;
            let mut peaks: Vec<(f64, f32)> = Vec::new();
            for &(mz, env, amp, apex) in species {
                let g = (-(rt - apex).powi(2) / (2.0 * w * w)).exp();
                for off in LOWEST_OFFSET..=HIGHEST_OFFSET {
                    let v = amp * env.at(off) * g;
                    if v > noise as f64 {
                        peaks.push((mz + off as f64 * ISOTOPE / 2.0, v as f32));
                    }
                }
            }
            for k in 0..50 {
                peaks.push((300.0 + k as f64 * 17.3, noise));
            }
            peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mut merged: Vec<(f64, f32)> = Vec::new();
            for p in peaks {
                match merged.last_mut() {
                    Some(l) if (p.0 - l.0).abs() < l.0 * 8e-6 => l.1 += p.1,
                    _ => merged.push(p),
                }
            }
            Ms1Scan {
                rt,
                peaks: merged,
                id: format!("scan={i}"),
                ..Default::default()
            }
        })
        .collect();
    Ms1RunIndex::new(scans)
}

fn probe(mzs: &[f64], envs: &[ChannelEnvelope], t_hat: f64) -> GroupProbe {
    GroupProbe {
        charge: 2,
        channel_mz: mzs.to_vec(),
        envelopes: envs.to_vec(),
        t_hat,
        s: 3.0,
        w_prior: 4.0,
        rt_shift: vec![0.0; mzs.len()],
        free: vec![true; mzs.len()],
        mobility: None,
        pinned: false,
    }
}

#[test]
fn active_set_matches_the_subset_enumeration() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut rnd = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut settled = 0;
    for trial in 0..400 {
        let n = 1 + trial % MAX_COMPONENTS;
        let rows = n + 3;
        let mut x: Vec<f64> = (0..rows * n).map(|_| rnd()).collect();
        if trial % 4 == 0 && n >= 2 {
            for r in 0..rows {
                x[r * n + 1] = x[r * n] * (1.0 + 1e-4 * rnd());
            }
        }
        let mut h = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                h[i * n + j] = (0..rows).map(|r| x[r * n + i] * x[r * n + j]).sum();
            }
        }
        let b: Vec<f64> = (0..n).map(|_| 2.0 * rnd() - 0.5).collect();
        let free: Vec<bool> = (0..n).map(|c| trial % 7 != 3 || c != 0).collect();
        let (a_e, llr_e, _) = nn_enumerate(&h, &b, n, &free);
        let Some((a_s, llr_s, passive)) = nn_active_set(&h, &b, n, &free) else {
            continue;
        };
        settled += 1;
        assert!(
            (llr_e - llr_s).abs() <= 1e-9 * llr_e.abs().max(1.0),
            "trial {trial} n {n}: {llr_e} vs {llr_s}"
        );
        let scale = b.iter().chain(&h).fold(0.0f64, |m, &x| m.max(x.abs()));
        for c in 0..n {
            let g = b[c] - (0..n).map(|k| h[c * n + k] * a_s[k]).sum::<f64>();
            if passive[c] {
                assert!(
                    g.abs() <= 1e-7 * scale,
                    "trial {trial}: passive {c} gradient {g}"
                );
            } else if free[c] {
                assert!(g <= 1e-7 * scale, "trial {trial}: free {c} gradient {g}");
            }
        }
        for c in 0..n {
            assert!(a_s[c] >= 0.0 && (free[c] || a_s[c] == 0.0));
            assert_eq!(passive[c], a_s[c] > 0.0);
        }
        if trial % 4 != 0 {
            for c in 0..n {
                assert!((a_e[c] - a_s[c]).abs() <= 1e-6 * (1.0 + a_e[c].abs()));
            }
        }
    }
    assert_eq!(settled, 400, "active set settled on {settled} of 400");
}

#[test]
fn fit_takes_at_most_max_components() {
    let light = [0.55, 0.3, 0.11, 0.04];
    let env = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    let mzs: Vec<f64> = (0..=MAX_COMPONENTS)
        .map(|i| 500.0 + 7.0 * i as f64)
        .collect();
    let species: Vec<(f64, ChannelEnvelope, f64, f64)> =
        mzs.iter().map(|&m| (m, env, 1e6, 120.0)).collect();
    let run = synthetic_run(&species, 4.0, 50.0);
    let envs = vec![env; mzs.len()];
    let fit = fit_group(
        &run,
        &probe(&mzs[..MAX_COMPONENTS], &envs[..MAX_COMPONENTS], 121.0),
        &FitParams::default(),
    )
    .expect("a fit of MAX_COMPONENTS components");
    for c in 1..MAX_COMPONENTS {
        let (r, _) = log2_ratio(&fit.channels[c], &fit.channels[0]).unwrap();
        assert!(r.abs() < 0.05, "component {c}: {r}");
    }
    assert!(fit_group(&run, &probe(&mzs, &envs, 121.0), &FitParams::default()).is_none());
}

#[test]
fn pinned_probe_keeps_the_profile() {
    let light = [0.55, 0.3, 0.11, 0.04];
    let env = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    let (ml, mo) = (600.3, 612.8);
    let run = synthetic_run(&[(ml, env, 1e6, 100.0), (mo, env, 2e7, 108.0)], 4.0, 50.0);
    let mut p = probe(&[ml, mo], &[env, env], 100.0);
    let free = fit_group(&run, &p, &FitParams::default()).unwrap();
    assert!((free.t0 - 108.0).abs() < 1.5, "free fit moves: {}", free.t0);
    p.pinned = true;
    p.w_prior = 4.0;
    let pinned = fit_group(&run, &p, &FitParams::default()).unwrap();
    assert_eq!(pinned.t0, 100.0);
    assert_eq!(pinned.w, 4.0);
    assert!(
        pinned.channels[1].llr < 0.2 * free.channels[1].llr,
        "{} vs {}",
        pinned.channels[1].llr,
        free.channels[1].llr
    );
    let (r, _) = log2_ratio(&pinned.channels[0], &free.channels[0]).unwrap();
    assert!(r > 0.5, "{r}");
}

#[test]
fn duplex_ratio_is_recovered() {
    let light = [0.55, 0.3, 0.11, 0.04];
    let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    let eh = ChannelEnvelope::new(&light, 0, 1.0, 8, 0.99);
    let ml = 600.3;
    let mh = ml + 8.014199 / 2.0;
    for ratio in [1.0, 0.25, 0.05] {
        let run = synthetic_run(
            &[(ml, el, 1e6, 120.0), (mh, eh, 1e6 * ratio, 120.0)],
            4.0,
            50.0,
        );
        let fit = fit_group(
            &run,
            &probe(&[ml, mh], &[el, eh], 121.0),
            &FitParams::default(),
        )
        .unwrap();
        let (r, se) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
        assert!((r - ratio.log2()).abs() < 0.05, "ratio {ratio}: {r} ± {se}");
        assert!((fit.t0 - 120.0).abs() <= 1.0);
        assert!(fit.channels[1].envelope_cos > 0.95);
        assert!(
            fit.channels[1].m1_ratio < 0.1,
            "{}",
            fit.channels[1].m1_ratio
        );
    }
}

#[test]
fn absent_channel_is_below_detection() {
    let light = [0.55, 0.3, 0.11, 0.04];
    let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    let ml = 600.3;
    let mh = ml + 8.014199 / 2.0;
    let run = synthetic_run(&[(ml, el, 1e6, 120.0)], 4.0, 50.0);
    let fit = fit_group(
        &run,
        &probe(&[ml, mh], &[el, el], 121.0),
        &FitParams::default(),
    )
    .unwrap();
    let h = &fit.channels[1];
    assert!(h.amplitude < 1e3, "{}", h.amplitude);
    let status = classify_channel(
        h,
        0.0,
        fit.channels[0].amplitude,
        false,
        &ClassifyParams::default(),
    );
    assert_eq!(status, ChannelStatus::BelowDetection);
}

#[test]
fn overlapping_lys4_channel_is_unmixed() {
    let light = [0.25, 0.30, 0.22, 0.12, 0.07, 0.04];
    let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    let em = ChannelEnvelope::new(&light, 4, 0.96, 0, 1.0);
    let z = 2.0;
    let ml = 1200.4;
    let mm = ml + 4.025107 / z;
    for ratio in [1.0, 0.2, 0.05] {
        let run = synthetic_run(
            &[(ml, el, 1e6, 100.0), (mm, em, 1e6 * ratio, 100.0)],
            4.0,
            50.0,
        );
        let fit = fit_group(
            &run,
            &probe(&[ml, mm], &[el, em], 100.0),
            &FitParams::default(),
        )
        .unwrap();
        let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
        assert!((r - ratio.log2()).abs() < 0.05, "ratio {ratio}: {r}");
        assert!(
            fit.channels[1].m1_ratio < 0.3,
            "{}",
            fit.channels[1].m1_ratio
        );
        assert!(fit.channels[1].envelope_cos > 0.95);
    }
}

#[test]
fn deuterium_shift_is_learned_and_fitted() {
    let light = [0.55, 0.3, 0.11, 0.04];
    let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    let eh = ChannelEnvelope::new(&light, 6, 0.99, 2, 0.99);
    let ml = 700.0;
    let mh = ml + (36.07567 - 28.0313) / 2.0;
    let run = synthetic_run(&[(ml, el, 1e6, 120.0), (mh, eh, 5e5, 117.0)], 4.0, 50.0);
    let mut alone = probe(&[ml, mh], &[el, eh], 120.0);
    alone.free = vec![true, false];
    let a0 = fit_group(&run, &alone, &FitParams::default()).unwrap().t0;
    alone.free = vec![false, true];
    let a1 = fit_group(&run, &alone, &FitParams::default()).unwrap().t0;
    let shifts = learn_rt_shifts(&[(1, a1 - a0)], 2, 1);
    assert!((shifts[1] + 3.0).abs() <= 1.0, "{shifts:?}");
    let mut p = probe(&[ml, mh], &[el, eh], 120.0);
    p.rt_shift = shifts;
    let fit = fit_group(&run, &p, &FitParams::default()).unwrap();
    let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
    assert!((r + 1.0).abs() < 0.1, "{r}");
}

#[test]
fn single_trace_interferer_does_not_move_the_apex() {
    let light = [0.55, 0.3, 0.11, 0.04];
    let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    let eh = ChannelEnvelope::new(&light, 0, 1.0, 8, 0.99);
    let lone = ChannelEnvelope::new(&[1.0], 0, 1.0, 0, 1.0);
    let ml = 600.3;
    let mh = ml + 8.014199 / 2.0;
    let other = mh + 2.0 * ISOTOPE / 2.0;
    let run = synthetic_run(
        &[
            (ml, el, 1e6, 100.0),
            (mh, eh, 1e6, 100.0),
            (other, lone, 8e6, 107.0),
        ],
        4.0,
        50.0,
    );
    let fit = fit_group(
        &run,
        &probe(&[ml, mh], &[el, eh], 100.0),
        &FitParams::default(),
    )
    .unwrap();
    assert!((fit.t0 - 100.0).abs() < 1.0, "apex {}", fit.t0);
    let (r, _) = log2_ratio(&fit.channels[1], &fit.channels[0]).unwrap();
    assert!(r.abs() < 0.15, "{r}");
}

#[test]
fn unmeasured_channel_earns_no_mass_credit() {
    assert_eq!(mass_factor(None, 0.0, 2.0), 0.0);
    assert!((mass_factor(Some(0.0), 0.0, 2.0) - 1.0).abs() < 1e-12);
    assert!(mass_factor(Some(4.0), 0.0, 2.0) < 0.2);
}

#[test]
fn interferer_at_channel_is_flagged() {
    let light = [0.55, 0.3, 0.11, 0.04];
    let el = ChannelEnvelope::new(&light, 0, 1.0, 0, 1.0);
    let ml = 600.3;
    let mh = ml + 8.014199 / 2.0;
    let other = mh - ISOTOPE / 2.0;
    let run = synthetic_run(&[(ml, el, 1e6, 120.0), (other, el, 2e6, 121.0)], 4.0, 50.0);
    let fit = fit_group(
        &run,
        &probe(&[ml, mh], &[el, el], 121.0),
        &FitParams::default(),
    )
    .unwrap();
    let h = &fit.channels[1];
    let status = classify_channel(
        h,
        0.0,
        fit.channels[0].amplitude,
        false,
        &ClassifyParams::default(),
    );
    assert_ne!(status, ChannelStatus::Matched, "{h:?}");
}

#[test]
fn qvalues_rank_targets_over_decoys() {
    let rows = vec![
        (10.0, false),
        (9.0, false),
        (8.0, true),
        (7.0, false),
        (1.0, true),
    ];
    let q = channel_qvalues(&rows, 1);
    assert!((q[0] - 0.5).abs() < 1e-12);
    let q4 = channel_qvalues(&rows, 4);
    assert!((q4[0] - 0.125).abs() < 1e-12);
    let w: Vec<(f64, f64)> = rows
        .iter()
        .map(|&(s, d)| (s, if d { 0.25 } else { 0.0 }))
        .collect();
    let qw = weighted_qvalues(&w);
    for (a, b) in qw.iter().zip(&q4) {
        assert!((a - b).abs() < 1e-12);
    }
    assert!(q[3] >= q[0]);
    assert_eq!(q.len(), 5);
}

#[test]
fn decoy_shifts_avoid_every_isotope() {
    let deltas = [0.0, 8.014199];
    let min_sep = decoy_min_sep(1200.0, 10.0);
    let all = decoy_shifts(&[0, 1], &deltas, min_sep, 0, 12);
    assert!(!all.is_empty());
    for s in &all {
        for m in [0, 1] {
            for d in deltas {
                for k in -3..=6 {
                    assert!(
                        (deltas[m] + s - d - k as f64 * ISOTOPE).abs() >= min_sep,
                        "{s}"
                    );
                }
            }
        }
    }
    assert!(!all
        .iter()
        .any(|s| (s + 7.0 * crate::decoys::PEPTIDE_MASS_SLOPE).abs() < 1e-9));
    let heavy_only = decoy_shifts(&[1], &deltas, min_sep, 3, 4);
    assert_eq!(heavy_only.len(), 4);
}
