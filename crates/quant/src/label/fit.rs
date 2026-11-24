//! The joint fit of a group's channel envelopes over its MS1 traces: non-negative least squares
//! with free or fixed widths, and the fit's quality measures.

use model::mass::ISOTOPE;
use model::stats::{median_mut, Even, MAD_TO_SD};

use super::{HIGHEST_OFFSET, LOWEST_OFFSET, MAX_CHANNELS, MAX_COMPONENTS, N_OFFSETS, READ_OFFSETS};
use crate::ms1::{MobilityGate, Ms1RunIndex};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelEnvelope(pub [f64; N_OFFSETS]);

impl ChannelEnvelope {
    pub fn new(light: &[f64], deuterium: u32, p_2h: f64, heavy: u32, p_heavy: f64) -> Self {
        let lower = |n: u32, p: f64| -> [f64; 3] {
            if n == 0 || p >= 1.0 {
                return [1.0, 0.0, 0.0];
            }
            let q = (1.0 - p).clamp(0.0, 1.0);
            let nf = n as f64;
            let p0 = p.powf(nf);
            let p1 = nf * q * p.powf(nf - 1.0);
            let p2 = if n >= 2 {
                nf * (nf - 1.0) / 2.0 * q * q * p.powf(nf - 2.0)
            } else {
                0.0
            };
            [p0, p1, p2]
        };
        let a = lower(deuterium, p_2h);
        let b = lower(heavy, p_heavy);
        let mut shift = [0.0; 3];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                if i + j < 3 {
                    shift[i + j] += x * y;
                }
            }
        }
        let mut out = [0.0; N_OFFSETS];
        for (k, &v) in light.iter().take((HIGHEST_OFFSET + 1) as usize).enumerate() {
            for (j, &s) in shift.iter().enumerate() {
                let idx = k as i32 - j as i32 - LOWEST_OFFSET;
                out[idx as usize] += v * s;
            }
        }
        let total: f64 = out.iter().sum();
        if total > 0.0 {
            for v in &mut out {
                *v /= total;
            }
        }
        ChannelEnvelope(out)
    }

    pub fn at(&self, offset: i32) -> f64 {
        if (LOWEST_OFFSET..=HIGHEST_OFFSET).contains(&offset) {
            self.0[(offset - LOWEST_OFFSET) as usize]
        } else {
            0.0
        }
    }
}

#[derive(Debug, Clone)]
pub struct GroupProbe {
    pub charge: u8,
    pub channel_mz: Vec<f64>,
    pub envelopes: Vec<ChannelEnvelope>,
    pub t_hat: f64,
    pub s: f64,
    pub w_prior: f64,
    pub rt_shift: Vec<f64>,
    pub free: Vec<bool>,
    pub mobility: Option<MobilityGate>,
    pub pinned: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct FitParams {
    pub tol_ppm: f64,
    pub widths: &'static [f64],
}

pub const DEFAULT_WIDTHS: &[f64] = &[0.7, 1.0, 1.4];
pub const FREE_WIDTHS: &[f64] = &[0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0];

impl Default for FitParams {
    fn default() -> Self {
        Self {
            tol_ppm: 10.0,
            widths: DEFAULT_WIDTHS,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChannelFit {
    pub amplitude: f64,
    pub amplitude_se: f64,
    pub quantity: f64,
    pub quantity_se: f64,
    pub llr: f64,
    pub envelope_cos: f64,
    pub coherence: f64,
    pub m1_ratio: f64,
    pub mass_err_ppm: Option<f64>,
    pub apex_offset: f64,
    pub trace_quantity: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GroupFit {
    pub t0: f64,
    pub w: f64,
    pub llr: f64,
    pub rt_z: f64,
    pub n_points: u32,
    pub channels: Vec<ChannelFit>,
}

pub(super) struct Traces {
    mz: Vec<f64>,
    mix: Vec<f64>,
    own: Vec<[usize; READ_OFFSETS.len()]>,
}

fn build_traces(probe: &GroupProbe, tol_ppm: f64) -> Traces {
    let z = probe.charge.max(1) as f64;
    let n_ch = probe.channel_mz.len();
    let mut mz: Vec<f64> = Vec::new();
    let tol = |m: f64| m * tol_ppm * 1e-6;
    let mut own = vec![[0usize; READ_OFFSETS.len()]; n_ch];
    for (c, &mc) in probe.channel_mz.iter().enumerate() {
        for (k, &off) in READ_OFFSETS.iter().enumerate() {
            let m = mc + off as f64 * ISOTOPE / z;
            let t = match mz.iter().position(|&x| (x - m).abs() <= tol(m)) {
                Some(t) => t,
                None => {
                    mz.push(m);
                    mz.len() - 1
                }
            };
            own[c][k] = t;
        }
    }
    let n_t = mz.len();
    let mut mix = vec![0.0; n_t * n_ch];
    for (t, &mt) in mz.iter().enumerate() {
        for (c, &mc) in probe.channel_mz.iter().enumerate() {
            for off in LOWEST_OFFSET..=HIGHEST_OFFSET {
                let m = mc + off as f64 * ISOTOPE / z;
                if (m - mt).abs() <= tol(mt) {
                    mix[t * n_ch + c] += probe.envelopes[c].at(off);
                }
            }
        }
    }
    Traces { mz, mix, own }
}

fn nn_quadratic(h: &[f64], b: &[f64], n: usize, free: &[bool]) -> (Vec<f64>, f64, Vec<bool>) {
    if n > MAX_CHANNELS {
        if let Some(r) = nn_active_set(h, b, n, free) {
            return r;
        }
    }
    nn_enumerate(h, b, n, free)
}

pub(super) fn nn_active_set(
    h: &[f64],
    b: &[f64],
    n: usize,
    free: &[bool],
) -> Option<(Vec<f64>, f64, Vec<bool>)> {
    let scale = b.iter().fold(0.0f64, |m, &x| m.max(x.abs()));
    if scale == 0.0 {
        return Some((vec![0.0; n], 0.0, vec![false; n]));
    }
    let tol = 1e-12 * scale;
    let mut a = vec![0.0; n];
    let mut passive = vec![false; n];
    let mut held = vec![false; n];
    let solve_on = |p: &[bool]| -> Option<Vec<f64>> {
        let idx: Vec<usize> = (0..n).filter(|&c| p[c]).collect();
        let k = idx.len();
        let mut m = vec![0.0; k * k];
        let mut v = vec![0.0; k];
        for (i, &ci) in idx.iter().enumerate() {
            v[i] = b[ci];
            for (j, &cj) in idx.iter().enumerate() {
                m[i * k + j] = h[ci * n + cj];
            }
        }
        let x = solve(&mut m, &mut v, k)?;
        let mut z = vec![0.0; n];
        for (&c, &xc) in idx.iter().zip(&x) {
            z[c] = xc;
        }
        Some(z)
    };
    for _ in 0..3 * n + 3 {
        let grad = |c: usize| b[c] - (0..n).map(|j| h[c * n + j] * a[j]).sum::<f64>();
        let Some(j) = (0..n)
            .filter(|&c| free[c] && !passive[c] && !held[c])
            .map(|c| (c, grad(c)))
            .filter(|&(_, g)| g > tol)
            .max_by(|x, y| x.1.total_cmp(&y.1))
            .map(|(c, _)| c)
        else {
            let llr = (0..n).map(|c| a[c] * b[c]).sum::<f64>() / 2.0;
            return Some((a, llr, passive));
        };
        passive[j] = true;
        let mut inner = 0;
        loop {
            inner += 1;
            if inner > 3 * n + 3 {
                return None;
            }
            let z = solve_on(&passive)?;
            if (0..n).filter(|&c| passive[c]).all(|c| z[c] > 0.0) {
                a = z;
                held.iter_mut().for_each(|h| *h = false);
                break;
            }
            if inner == 1 && z[j] <= 0.0 {
                passive[j] = false;
                held[j] = true;
                break;
            }
            let alpha = (0..n)
                .filter(|&c| passive[c] && z[c] <= 0.0)
                .map(|c| a[c] / (a[c] - z[c]))
                .fold(f64::INFINITY, f64::min);
            for c in 0..n {
                if passive[c] {
                    a[c] += alpha * (z[c] - a[c]);
                }
            }
            let amax = a.iter().fold(0.0f64, |m, &x| m.max(x.abs()));
            for c in 0..n {
                if passive[c] && a[c] <= 1e-12 * amax {
                    a[c] = 0.0;
                    passive[c] = false;
                }
            }
        }
    }
    None
}

pub(super) fn nn_enumerate(
    h: &[f64],
    b: &[f64],
    n: usize,
    free: &[bool],
) -> (Vec<f64>, f64, Vec<bool>) {
    let mut best = (vec![0.0; n], 0.0, vec![false; n]);
    for mask in 1u32..(1 << n) {
        if (0..n).any(|c| mask & (1 << c) != 0 && !free[c]) {
            continue;
        }
        let idx: Vec<usize> = (0..n).filter(|c| mask & (1 << c) != 0).collect();
        let k = idx.len();
        let mut m = vec![0.0; k * k];
        let mut v = vec![0.0; k];
        for (i, &ci) in idx.iter().enumerate() {
            v[i] = b[ci];
            for (j, &cj) in idx.iter().enumerate() {
                m[i * k + j] = h[ci * n + cj];
            }
        }
        let Some(x) = solve(&mut m, &mut v, k) else {
            continue;
        };
        if x.iter().any(|&a| a < 0.0) {
            continue;
        }
        let llr: f64 = idx.iter().zip(&x).map(|(&c, &a)| a * b[c]).sum::<f64>() / 2.0;
        if llr > best.1 {
            let mut a = vec![0.0; n];
            let mut p = vec![false; n];
            for (&c, &xa) in idx.iter().zip(&x) {
                a[c] = xa;
                p[c] = true;
            }
            best = (a, llr, p);
        }
    }
    best
}

fn solve(m: &mut [f64], v: &mut [f64], k: usize) -> Option<Vec<f64>> {
    for col in 0..k {
        let piv =
            (col..k).max_by(|&a, &b| m[a * k + col].abs().total_cmp(&m[b * k + col].abs()))?;
        if m[piv * k + col].abs() < 1e-300 {
            return None;
        }
        if piv != col {
            for j in 0..k {
                m.swap(piv * k + j, col * k + j);
            }
            v.swap(piv, col);
        }
        for r in col + 1..k {
            let f = m[r * k + col] / m[col * k + col];
            for j in col..k {
                m[r * k + j] -= f * m[col * k + j];
            }
            v[r] -= f * v[col];
        }
    }
    let mut x = vec![0.0; k];
    for r in (0..k).rev() {
        let mut s = v[r];
        for j in r + 1..k {
            s -= m[r * k + j] * x[j];
        }
        x[r] = s / m[r * k + r];
        if !x[r].is_finite() {
            return None;
        }
    }
    Some(x)
}

fn invert(h: &[f64], k: usize) -> Option<Vec<f64>> {
    let mut inv = vec![0.0; k * k];
    for col in 0..k {
        let mut m = h.to_vec();
        let mut e = vec![0.0; k];
        e[col] = 1.0;
        let x = solve(&mut m, &mut e, k)?;
        for r in 0..k {
            inv[r * k + col] = x[r];
        }
    }
    Some(inv)
}

fn robust_sigma(y: &[f64], floor: f64) -> f64 {
    if y.len() < 2 {
        return floor.max(1.0);
    }
    let mut d: Vec<f64> = y.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    let med = median_mut(&mut d, Even::Upper).unwrap_or(0.0);
    (MAD_TO_SD * med / std::f64::consts::SQRT_2)
        .max(floor)
        .max(1e-9)
}

pub fn fit_group(ms1: &Ms1RunIndex, probe: &GroupProbe, params: &FitParams) -> Option<GroupFit> {
    let n_ch = probe.channel_mz.len();
    if n_ch == 0 || n_ch > MAX_COMPONENTS {
        return None;
    }
    let s = probe.s.max(1e-3);
    let w_prior = probe.w_prior.max(1e-3);
    let widths: &[f64] = if probe.pinned { &[1.0] } else { params.widths };
    let max_shift = probe.rt_shift.iter().fold(0.0f64, |a, &d| a.max(d.abs()));
    let w_max = w_prior * widths.iter().fold(1.0f64, |a, &b| a.max(b));
    let half = (4.0 * s + 4.0 * w_max).max(8.0 * w_max) + max_shift;
    let window = ms1.scan_window(probe.t_hat - half, probe.t_hat + half);
    if window.len() < 3 {
        return None;
    }
    let rt: Vec<f64> = window.clone().map(|i| ms1.rt(i)).collect();
    let traces = build_traces(probe, params.tol_ppm);
    let n_t = traces.mz.len();
    let read = |m: f64| -> Vec<f64> {
        let tol = m * params.tol_ppm * 1e-6;
        window
            .clone()
            .map(|i| ms1.slice_in(i, m - tol, m + tol, probe.mobility).1 as f64)
            .collect()
    };
    let y: Vec<Vec<f64>> = traces.mz.iter().map(|&m| read(m)).collect();
    let mut limits: Vec<f64> = window
        .clone()
        .map(|i| ms1.detection_limit(i) as f64)
        .collect();
    let limit = median_mut(&mut limits, Even::Upper).unwrap_or(f64::NAN);
    let away: Vec<bool> = rt
        .iter()
        .map(|&t| (t - probe.t_hat).abs() > 3.0 * s.max(w_max) + max_shift)
        .collect();
    let n_away = away.iter().filter(|&&a| a).count();
    let sigma: Vec<f64> = y
        .iter()
        .map(|t| {
            if n_away >= 10 {
                let base: Vec<f64> = t
                    .iter()
                    .zip(&away)
                    .filter(|(_, &a)| a)
                    .map(|(v, _)| *v)
                    .collect();
                robust_sigma(&base, 0.5 * limit)
            } else {
                robust_sigma(t, 0.5 * limit)
            }
        })
        .collect();
    let w_inv: Vec<f64> = sigma.iter().map(|s| 1.0 / (s * s)).collect();
    let reach = 3.0 * s + 3.0 * w_max + max_shift;
    let j0 = rt.partition_point(|&t| t < probe.t_hat - reach);
    let j1 = rt.partition_point(|&t| t <= probe.t_hat + reach);
    let (rt, y, window) = if j1 >= j0 + 5 {
        (
            rt[j0..j1].to_vec(),
            y.iter().map(|t| t[j0..j1].to_vec()).collect::<Vec<_>>(),
            window.start + j0..window.start + j1,
        )
    } else {
        (rt, y, window)
    };
    let n_j = rt.len();
    let ybar: Vec<f64> = y
        .iter()
        .map(|t| t.iter().sum::<f64>() / n_j as f64)
        .collect();
    let yc: Vec<Vec<f64>> = y
        .iter()
        .zip(&ybar)
        .map(|(t, m)| t.iter().map(|v| v - m).collect())
        .collect();

    let profiles = |t0: f64, w: f64| -> Vec<Vec<f64>> {
        (0..n_ch)
            .map(|c| {
                let tc = t0 + probe.rt_shift.get(c).copied().unwrap_or(0.0);
                let g: Vec<f64> = rt
                    .iter()
                    .map(|&t| (-(t - tc) * (t - tc) / (2.0 * w * w)).exp())
                    .collect();
                let gbar = g.iter().sum::<f64>() / n_j as f64;
                g.into_iter().map(|v| v - gbar).collect()
            })
            .collect()
    };
    let normal =
        |gc: &[Vec<f64>], w_inv: &[f64]| -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>, Vec<f64>) {
            let mut sgg = vec![0.0; n_ch * n_ch];
            for a in 0..n_ch {
                for b in a..n_ch {
                    let v: f64 = gc[a].iter().zip(&gc[b]).map(|(x, y)| x * y).sum();
                    sgg[a * n_ch + b] = v;
                    sgg[b * n_ch + a] = v;
                }
            }
            let proj: Vec<Vec<f64>> = yc
                .iter()
                .map(|t| {
                    (0..n_ch)
                        .map(|c| t.iter().zip(&gc[c]).map(|(x, y)| x * y).sum())
                        .collect()
                })
                .collect();
            let mut h = vec![0.0; n_ch * n_ch];
            let mut b = vec![0.0; n_ch];
            for t in 0..n_t {
                let mt = &traces.mix[t * n_ch..(t + 1) * n_ch];
                for a in 0..n_ch {
                    if mt[a] == 0.0 {
                        continue;
                    }
                    b[a] += w_inv[t] * mt[a] * proj[t][a];
                    for c in 0..n_ch {
                        h[a * n_ch + c] += w_inv[t] * mt[a] * mt[c] * sgg[a * n_ch + c];
                    }
                }
            }
            (h, b, proj, sgg)
        };

    let lo = probe.t_hat - 3.0 * s;
    let hi = probe.t_hat + 3.0 * s;
    let mut t0s: Vec<f64> = if probe.pinned {
        vec![probe.t_hat]
    } else {
        rt.iter().copied().filter(|&t| t >= lo && t <= hi).collect()
    };
    if t0s.is_empty() {
        let k = rt
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1 - probe.t_hat)
                    .abs()
                    .total_cmp(&(b.1 - probe.t_hat).abs())
            })
            .map(|x| x.0)?;
        t0s.push(rt[k]);
    }
    let search_apex = |wt: &[f64]| -> Option<(f64, f64)> {
        let mut best: Option<(f64, f64, f64)> = None;
        for &t0 in &t0s {
            for &wm in widths {
                let w = wm * w_prior;
                let gc = profiles(t0, w);
                if gc
                    .iter()
                    .all(|g| g.iter().map(|v| v * v).sum::<f64>() < 1e-9)
                {
                    continue;
                }
                let (h, b, _, _) = normal(&gc, wt);
                let (_, llr, _) = nn_quadratic(&h, &b, n_ch, &probe.free);
                let obj = llr - (t0 - probe.t_hat).powi(2) / (2.0 * s * s);
                if best.is_none_or(|(o, _, _)| obj > o) {
                    best = Some((obj, t0, w));
                }
            }
        }
        let (best_obj, t0_grid, w) = best?;
        let objective = |t0: f64| -> f64 {
            let gc = profiles(t0, w);
            let (h, b, _, _) = normal(&gc, wt);
            let (_, llr, _) = nn_quadratic(&h, &b, n_ch, &probe.free);
            llr - (t0 - probe.t_hat).powi(2) / (2.0 * s * s)
        };
        let mut t0 = t0_grid;
        if let Some(k) = t0s.iter().position(|&x| x == t0_grid) {
            if k > 0 && k + 1 < t0s.len() {
                let (x0, x1, x2) = (t0s[k - 1], t0_grid, t0s[k + 1]);
                let (y0, y1, y2) = (objective(x0), best_obj, objective(x2));
                let den = (x0 - x1) * (x0 - x2) * (x1 - x2);
                if den.abs() > 0.0 {
                    let a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den;
                    let b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den;
                    if a < 0.0 {
                        let v = (-b / (2.0 * a)).clamp(x0, x2);
                        if v.is_finite() && objective(v) > best_obj {
                            t0 = v;
                        }
                    }
                }
            }
        }
        Some((t0, w))
    };
    let (mut t0, mut w) = search_apex(&w_inv)?;
    let mut w_fit = w_inv.clone();
    let mut gc = profiles(t0, w);
    let (mut h, mut b, mut proj, mut sgg) = normal(&gc, &w_fit);
    let (mut a, mut llr, mut passive) = nn_quadratic(&h, &b, n_ch, &probe.free);
    for round in 0..3 {
        let mut reweighted = false;
        for _ in 0..2 {
            let mut changed = false;
            for t in 0..n_t {
                let mt = &traces.mix[t * n_ch..(t + 1) * n_ch];
                if mt.iter().all(|&m| m == 0.0) {
                    continue;
                }
                let mut ss = 0.0;
                let mut n = 0usize;
                for (j, &tj) in rt.iter().enumerate() {
                    if (tj - t0).abs() > 2.0 * w {
                        continue;
                    }
                    let model: f64 = (0..n_ch).map(|c| mt[c] * a[c] * gc[c][j]).sum();
                    let r = yc[t][j] - model;
                    ss += r * r;
                    n += 1;
                }
                if n < 3 {
                    continue;
                }
                let chi = ss / n as f64 * w_inv[t];
                let target = if chi > 9.0 {
                    w_inv[t] * 9.0 / chi
                } else {
                    w_inv[t]
                };
                if (target - w_fit[t]).abs() > 1e-12 * w_inv[t] {
                    w_fit[t] = target;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
            reweighted = true;
            (h, b, proj, sgg) = normal(&gc, &w_fit);
            (a, llr, passive) = nn_quadratic(&h, &b, n_ch, &probe.free);
        }
        if !reweighted || round == 2 {
            break;
        }
        let Some((t1, w1)) = search_apex(&w_fit) else {
            break;
        };
        if t1 == t0 && w1 == w {
            break;
        }
        (t0, w) = (t1, w1);
        w_fit.clone_from(&w_inv);
        gc = profiles(t0, w);
        (h, b, proj, sgg) = normal(&gc, &w_fit);
        (a, llr, passive) = nn_quadratic(&h, &b, n_ch, &probe.free);
    }
    let area = w * (2.0 * std::f64::consts::PI).sqrt();

    let pidx: Vec<usize> = (0..n_ch).filter(|&c| passive[c]).collect();
    let mut se = vec![0.0; n_ch];
    let k = pidx.len();
    let mut block = vec![0.0; k * k];
    for (i, &ci) in pidx.iter().enumerate() {
        for (j, &cj) in pidx.iter().enumerate() {
            block[i * k + j] = h[ci * n_ch + cj];
        }
    }
    let inv = if k > 0 { invert(&block, k) } else { None };
    for c in 0..n_ch {
        se[c] = match (inv.as_ref(), pidx.iter().position(|&p| p == c)) {
            (Some(inv), Some(i)) => inv[i * k + i].max(0.0).sqrt(),
            _ => 1.0 / h[c * n_ch + c].max(1e-300).sqrt(),
        };
    }

    let n_points = rt.iter().filter(|&&t| (t - t0).abs() <= 2.0 * w).count() as u32;
    let channels = (0..n_ch)
        .map(|c| {
            let mut free = probe.free.clone();
            free[c] = false;
            let (_, llr0, _) = nn_quadratic(&h, &b, n_ch, &free);
            let llr_c = (llr - llr0).max(0.0);
            let s_cc = sgg[c * n_ch + c].max(1e-300);
            let height = |t: usize| -> f64 {
                let mut p = proj[t][c];
                for c2 in 0..n_ch {
                    if c2 != c {
                        p -= traces.mix[t * n_ch + c2] * a[c2] * sgg[c * n_ch + c2];
                    }
                }
                p / s_cc
            };
            let mut obs = Vec::new();
            let mut exp = Vec::new();
            for (k, &off) in READ_OFFSETS.iter().enumerate() {
                if off < 0 {
                    continue;
                }
                let t = traces.own[c][k];
                obs.push(height(t));
                exp.push(traces.mix[t * n_ch + c]);
            }
            let envelope_cos = cosine(&obs, &exp);
            let t_ch = t0 + probe.rt_shift.get(c).copied().unwrap_or(0.0);
            let mono = traces.own[c][1];
            let mut xs = Vec::new();
            let mut gs = Vec::new();
            for (j, &t) in rt.iter().enumerate() {
                if (t - t_ch).abs() > 2.0 * w {
                    continue;
                }
                let mut v = yc[mono][j];
                for c2 in 0..n_ch {
                    if c2 != c {
                        v -= traces.mix[mono * n_ch + c2] * a[c2] * gc[c2][j];
                    }
                }
                xs.push(v);
                gs.push(gc[c][j]);
            }
            let coherence = if xs.len() >= 3 {
                pearson(&xs, &gs)
            } else {
                0.0
            };
            let m1 = traces.own[c][0];
            let mut r = proj[m1][c];
            for c2 in 0..n_ch {
                r -= traces.mix[m1 * n_ch + c2] * a[c2] * sgg[c * n_ch + c2];
            }
            let pi0 = probe.envelopes[c].at(0).max(1e-9);
            let m1_ratio = (r / s_cc).max(0.0) / (pi0 * a[c].max(se[c]).max(1e-300));
            let tc = t0 + probe.rt_shift.get(c).copied().unwrap_or(0.0);
            let (mut sy, mut syt) = (0.0, 0.0);
            for k in [1usize, 2] {
                let tr = traces.own[c][k];
                for (j, &t) in rt.iter().enumerate() {
                    if (t - tc).abs() > 2.0 * w {
                        continue;
                    }
                    let mut v = y[tr][j] - ybar[tr];
                    for c2 in 0..n_ch {
                        if c2 != c {
                            v -= traces.mix[tr * n_ch + c2] * a[c2] * gc[c2][j];
                        }
                    }
                    if v > 0.0 {
                        sy += v;
                        syt += v * (t - tc);
                    }
                }
            }
            let apex_offset = if sy > 0.0 { syt / sy } else { 0.0 };
            let trace_quantity = {
                let (lo, hi) = (tc - 3.0 * w, tc + 4.0 * w);
                let mut area = 0.0;
                let mut share = 0.0;
                for k in 1..READ_OFFSETS.len() {
                    let tr = traces.own[c][k];
                    let pi = traces.mix[tr * n_ch + c];
                    if pi <= 0.0 {
                        continue;
                    }
                    share += pi;
                    let v: Vec<f64> = (0..rt.len())
                        .map(|j| {
                            let mut x = yc[tr][j];
                            for c2 in 0..n_ch {
                                if c2 != c {
                                    x -= traces.mix[tr * n_ch + c2] * a[c2] * gc[c2][j];
                                }
                            }
                            x
                        })
                        .collect();
                    let out: Vec<f64> = (0..rt.len())
                        .filter(|&j| rt[j] < lo || rt[j] > hi)
                        .map(|j| v[j])
                        .collect();
                    let b = if out.is_empty() {
                        0.0
                    } else {
                        out.iter().sum::<f64>() / out.len() as f64
                    };
                    for j in 1..rt.len() {
                        if rt[j - 1] >= lo && rt[j] <= hi {
                            area += 0.5 * (v[j - 1] + v[j] - 2.0 * b) * (rt[j] - rt[j - 1]);
                        }
                    }
                }
                if share > 0.0 {
                    area / share
                } else {
                    0.0
                }
            };
            let mz = probe.channel_mz[c];
            let tol = mz * params.tol_ppm * 1e-6;
            let (mut sw, mut swe) = (0.0, 0.0);
            let t_ch = t0 + probe.rt_shift.get(c).copied().unwrap_or(0.0);
            for (j, i) in window.clone().enumerate() {
                if (rt[j] - t_ch).abs() > w {
                    continue;
                }
                if let Some((m, inten)) = ms1.apex_centroid(i, mz - tol, mz + tol) {
                    sw += inten as f64;
                    swe += inten as f64 * (m - mz) / mz * 1e6;
                }
            }
            ChannelFit {
                amplitude: a[c],
                amplitude_se: se[c],
                quantity: a[c] * area,
                quantity_se: se[c] * area,
                llr: llr_c,
                envelope_cos,
                coherence,
                m1_ratio,
                mass_err_ppm: (sw > 0.0).then(|| swe / sw),
                apex_offset,
                trace_quantity,
            }
        })
        .collect();
    Some(GroupFit {
        t0,
        w,
        llr,
        rt_z: (t0 - probe.t_hat) / s,
        n_points,
        channels,
    })
}

fn cosine(a: &[f64], b: &[f64]) -> f64 {
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let na: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na <= 0.0 || nb <= 0.0 {
        0.0
    } else {
        (dot / (na * nb)).clamp(-1.0, 1.0)
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma) * (x - ma);
        sbb += (y - mb) * (y - mb);
    }
    if saa <= 0.0 || sbb <= 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}
