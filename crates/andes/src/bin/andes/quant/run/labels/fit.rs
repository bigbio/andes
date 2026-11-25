//! Learning the fit of a run: envelope widths and channel RT shifts from the identified groups,
//! decoy probes, and the group score.

use model::mass::PROTON;
use quant::label::{
    fit_group, learn_rt_shifts, ChannelFit, FitParams, GroupFit, GroupProbe, LabelSpec, FREE_WIDTHS,
};
use quant::Ms1RunIndex;
use rayon::prelude::*;

use super::hits::{ChannelHit, Group};

pub(super) struct Fitted {
    pub(super) gid: usize,
    pub(super) run: usize,
    pub(super) fit: GroupFit,
    pub(super) ids: Vec<Option<ChannelHit>>,
    pub(super) decoy_scores: Vec<Vec<(f64, Option<f64>, f64)>>,
    pub(super) claimed: Vec<bool>,
    pub(super) transfer_q: Option<f64>,
    pub(super) rt_shift: Vec<f64>,
    pub(super) converted: Vec<Option<f64>>,
}

pub(super) struct Scales {
    pub(super) dt: f64,
    pub(super) w: f64,
    pub(super) shifts: Vec<f64>,
}

pub(super) const MIN_WIDTH_LLR: f64 = 25.0;
pub(super) const MIN_WIDTH_FITS: usize = 20;
const MAX_WIDTH_FITS: usize = 500;
const MIN_SHIFT_PAIRS: usize = 20;
pub(super) const MIN_NORM_RATIOS: usize = 10;
pub(super) const N_DECOYS: usize = 4;
pub(super) const N_TRANSFER_DECOYS: usize = 8;

pub(super) fn probe_of(g: &Group, t_hat: f64, s: f64, w: f64, shifts: &[f64]) -> GroupProbe {
    let mut channel_mz = g.channel_mz.clone();
    let mut envelopes = g.envelopes.clone();
    let mut rt_shift = shifts.to_vec();
    for &(parent, mz) in &g.satellites {
        channel_mz.push(mz);
        envelopes.push(g.envelopes[parent]);
        rt_shift.push(shifts[parent]);
    }
    GroupProbe {
        charge: g.charge,
        free: vec![true; channel_mz.len()],
        channel_mz,
        envelopes,
        t_hat,
        s,
        w_prior: w,
        rt_shift,
        mobility: None,
        pinned: false,
    }
}

pub(super) fn fit_folded(
    ms1: &Ms1RunIndex,
    g: &Group,
    probe: &GroupProbe,
    params: &FitParams,
) -> Option<(GroupFit, Vec<Option<f64>>)> {
    let mut fit = fit_group(ms1, probe, params)?;
    let n = g.channel_mz.len();
    let mut converted = vec![None; n];
    if g.satellites.is_empty() {
        return Some((fit, converted));
    }
    let sats: Vec<ChannelFit> = fit.channels.split_off(n);
    for (&(parent, _), sat) in g.satellites.iter().zip(&sats) {
        let c = &mut fit.channels[parent];
        let total = c.amplitude + sat.amplitude;
        converted[parent] = Some(if total > 0.0 {
            sat.amplitude / total
        } else {
            0.0
        });
        c.amplitude = total;
        c.quantity += sat.quantity;
        c.trace_quantity += sat.trace_quantity;
        c.amplitude_se = c.amplitude_se.hypot(sat.amplitude_se);
        c.quantity_se = c.quantity_se.hypot(sat.quantity_se);
        c.llr += sat.llr;
    }
    Some((fit, converted))
}

pub(super) type ProfileFit = ((GroupFit, Vec<Option<f64>>), GroupProbe);

pub(super) fn fit_at_identified_profile(
    ms1: &Ms1RunIndex,
    g: &Group,
    probe: &GroupProbe,
    ids: &[Option<ChannelHit>],
    params: &FitParams,
) -> Option<ProfileFit> {
    let n = g.channel_mz.len();
    let mut shape = probe.clone();
    shape.free = (0..shape.channel_mz.len())
        .map(|k| ids[if k < n { k } else { g.satellites[k - n].0 }].is_some())
        .collect();
    let (profile, _) = fit_folded(ms1, g, &shape, params)?;
    let mut pinned = probe.clone();
    pinned.t_hat = profile.t0;
    pinned.w_prior = profile.w;
    pinned.pinned = true;
    let (mut fit, converted) = fit_folded(ms1, g, &pinned, params)?;
    fit.rt_z = profile.rt_z;
    Some(((fit, converted), pinned))
}

pub(super) fn learn_width(
    ms1: &Ms1RunIndex,
    groups: &[Group],
    anchors: &[(usize, f64)],
    tol: f64,
) -> f64 {
    let dt = ms1.median_scan_spacing().unwrap_or(1.0).max(1e-3);
    let mut order: Vec<&(usize, f64)> = anchors.iter().collect();
    order.sort_by_key(|(gid, _)| seed_of(&groups[*gid].key, 0));
    let sample: Vec<&(usize, f64)> = order.into_iter().take(MAX_WIDTH_FITS).collect();
    let mut params = FitParams {
        tol_ppm: tol,
        widths: FREE_WIDTHS,
    };
    let mut w = 6.0f64.max(dt);
    for _ in 0..16 {
        let s = w.max(dt);
        let mut fitted: Vec<f64> = sample
            .par_iter()
            .filter_map(|(gid, t)| {
                let g = &groups[*gid];
                let p = probe_of(g, *t, s, w, &vec![0.0; g.channel_mz.len()]);
                fit_folded(ms1, g, &p, &params).map(|x| x.0)
            })
            .filter(|f| f.llr >= MIN_WIDTH_LLR)
            .map(|f| f.w.ln())
            .collect();
        if fitted.len() < MIN_WIDTH_FITS {
            return (6.0f64).max(dt);
        }
        fitted.sort_by(f64::total_cmp);
        let cut = fitted.len() / 10;
        let core = &fitted[cut..fitted.len() - cut];
        let next = (core.iter().sum::<f64>() / core.len() as f64).exp().max(dt);
        let moved = (next - w).abs() / w;
        w = next;
        if params.widths == FREE_WIDTHS {
            if moved < 0.1 {
                params.widths = REFINE_WIDTHS;
            }
        } else if moved < 0.01 {
            break;
        }
    }
    w
}

const REFINE_WIDTHS: &[f64] = &[0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.25];

pub(super) fn learn_shifts(
    ms1: &Ms1RunIndex,
    groups: &[Group],
    anchors: &[(usize, f64)],
    scales: (f64, f64),
    spec: &LabelSpec,
    tol: f64,
) -> Vec<f64> {
    let n = spec.len();
    if !spec.has_deuterium() {
        return vec![0.0; n];
    }
    let (dt, w) = scales;
    let params = FitParams {
        tol_ppm: tol,
        widths: quant::label::DEFAULT_WIDTHS,
    };
    let step = (anchors.len() / 2000).max(1);
    let r = spec.reference;
    let pairs: Vec<(usize, f64)> = anchors
        .par_iter()
        .step_by(step)
        .flat_map_iter(|(gid, t)| {
            let g = &groups[*gid];
            let alone = |c: usize| {
                let mut p = probe_of(g, *t, w.max(dt), w, &vec![0.0; n]);
                p.free = (0..p.channel_mz.len())
                    .map(|k| k == c || (k >= n && g.satellites[k - n].0 == c))
                    .collect();
                fit_folded(ms1, g, &p, &params)
                    .map(|x| x.0)
                    .filter(|f| f.llr >= MIN_WIDTH_LLR)
            };
            let Some(fr) = alone(r) else {
                return Vec::new();
            };
            (0..n)
                .filter(|&c| c != r)
                .filter_map(|c| alone(c).map(|f| (c, f.t0 - fr.t0)))
                .collect::<Vec<_>>()
        })
        .collect();
    learn_rt_shifts(&pairs, n, MIN_SHIFT_PAIRS)
}

pub(super) fn decoy_probes(
    g: &Group,
    probe: &GroupProbe,
    moved: &[usize],
    tol_ppm: f64,
    seed: u64,
    n: usize,
) -> Vec<GroupProbe> {
    let z = g.charge.max(1) as f64;
    let mass = (g.channel_mz[0] - PROTON) * z;
    let min_sep = quant::label::decoy_min_sep(mass, tol_ppm);
    quant::label::decoy_shifts(moved, &g.deltas, min_sep, seed, n)
        .into_iter()
        .map(|shift| {
            let mut p = probe.clone();
            for &c in moved {
                p.channel_mz[c] += shift / z;
            }
            let n_ch = g.channel_mz.len();
            for (k, &(parent, _)) in g.satellites.iter().enumerate() {
                if moved.contains(&parent) {
                    p.channel_mz[n_ch + k] += shift / z;
                }
            }
            p
        })
        .collect()
}

pub(super) fn group_score(fit: &GroupFit) -> f64 {
    let Some(best) = fit
        .channels
        .iter()
        .max_by(|a, b| a.amplitude.total_cmp(&b.amplitude))
    else {
        return 0.0;
    };
    (2.0 * fit.llr.max(0.0)).sqrt() * best.envelope_cos.max(0.0) * best.coherence.max(0.0)
}

pub(super) fn seed_of(s: &str, run: usize) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes().chain((run as u64).to_le_bytes()) {
        h ^= b as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}
