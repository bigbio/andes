use model::mass::ISOTOPE;
use model::stats::{median_f32_as_f64_mut, MAD_TO_SD};

use crate::decoys::Probe;
use crate::gaussian;
use crate::ms1_index::{slice_value, MobilityGate, Ms1RunIndex};

pub const MAX_ISOTOPES: usize = 4;
pub const MIN_SCANS: usize = 5;
pub const WIDTH_GRID: [f64; 3] = [0.7, 1.0, 1.4];
const DEFAULT_DETECTION_LIMIT: f64 = 1.0;
const NOISE_FLOOR_OF_LIMIT: f64 = 0.5;
const MIN_CENTRED_SGG: f64 = 1e-9;
const EXCLUSION_Z: f64 = 3.0;
const MIN_SCALE_S: f64 = 1e-3;
const TRACE_WIDTHS_BEFORE: f64 = 3.0;
const TRACE_WIDTHS_AFTER: f64 = 4.0;
const MIN_BASELINE_SCANS: usize = 3;
const BASELINE_SIDE_Z: f64 = 3.0;

#[derive(Debug, Clone, PartialEq)]
pub struct EvidenceParams {
    pub tol_ppm: f64,
    pub n_isotopes: usize,
}

impl Default for EvidenceParams {
    fn default() -> Self {
        Self {
            tol_ppm: 10.0,
            n_isotopes: MAX_ISOTOPES,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Evidence {
    pub amplitude: f64,
    pub amplitude_se: f64,
    pub quantity: f64,
    pub quantity_se: f64,
    pub llr: f64,
    pub rt_z: f64,
    pub width_ratio: f64,
    pub envelope_cos: f64,
    pub coherence: f64,
    pub m1_ratio: f64,
    pub mass_err_ppm: f64,
    pub n_points: u32,
    pub isotopes_used: u8,
    pub t0: f64,
    pub w: f64,
    pub iso_heights: [f64; MAX_ISOTOPES],
    pub tailing: f64,
    pub trace_quantity: f64,
    pub trace_quantity_se: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Prior {
    pub t_hat: f64,
    pub s: f64,
    pub w_prior: f64,
}

impl Prior {
    pub fn of(probe: &Probe) -> Self {
        Self {
            t_hat: probe.t_hat,
            s: probe.s.max(MIN_SCALE_S),
            w_prior: probe.w_prior.max(MIN_SCALE_S),
        }
    }

    pub fn half_window(&self) -> f64 {
        4.0 * self.s + 4.0 * self.w_prior
    }
}

#[derive(Debug, Default)]
pub struct FitScratch {
    rt: Vec<f64>,
    y: [Vec<f32>; MAX_ISOTOPES],
    m1: Vec<f32>,
    mz0: Vec<f64>,
    limits: Vec<f32>,
    work: Workspace,
}

#[derive(Debug, Default)]
pub struct Workspace {
    sort: Vec<f32>,
    z: Vec<f64>,
    g: Vec<f64>,
    grids: GridCache,
}

const GRID_CACHE_ENTRIES: usize = 3;
const GRID_CACHE_MAX_VALUES: usize = 1 << 21;
const PROFILE_REACH: f64 = 9.0;

fn support(rt: &[f64], finite: bool, t0: f64, w: f64) -> std::ops::Range<usize> {
    let reach = PROFILE_REACH * w;
    if !(finite && t0.is_finite() && reach.is_finite()) {
        return 0..rt.len();
    }
    let lo = rt.partition_point(|&t| t < t0 - reach);
    let hi = rt.partition_point(|&t| t <= t0 + reach);
    lo..hi.max(lo)
}

const APEX_SPACING_OF_WIDTH: f64 = 0.5;

#[derive(Debug, Clone)]
struct Apexes {
    range: std::ops::Range<usize>,
    stride: usize,
}

impl Apexes {
    fn of(rt: &[f64], range: std::ops::Range<usize>, w_min: f64) -> Self {
        let n = range.len();
        let mut stride = 1;
        if n > 2 {
            let spacing = (rt[range.end - 1] - rt[range.start]) / (n - 1) as f64;
            let k = APEX_SPACING_OF_WIDTH * w_min / spacing;
            if k.is_finite() && k >= 2.0 {
                stride = (k as usize).min(n - 1);
            }
        }
        Self { range, stride }
    }

    fn iter(&self) -> impl Iterator<Item = usize> + Clone + '_ {
        let last = self.range.end - 1;
        let tail = ((last - self.range.start) % self.stride != 0).then_some(last);
        self.range.clone().step_by(self.stride).chain(tail)
    }

    fn around(&self, k: usize) -> impl Iterator<Item = usize> + '_ {
        let lo = k.saturating_sub(self.stride - 1).max(self.range.start);
        let hi = (k + self.stride).min(self.range.end);
        (lo..hi).filter(move |&j| j != k)
    }
}

#[derive(Debug, Clone)]
struct Cell {
    scans: std::ops::Range<usize>,
    at: usize,
    sg: f64,
    sgg: f64,
}

#[derive(Debug, Default)]
struct GridProfiles {
    rt: Vec<u64>,
    prior: [u64; 3],
    widths: Vec<u64>,
    g: Vec<f64>,
    cells: Vec<Cell>,
}

impl GridProfiles {
    fn matches(&self, rt: &[f64], prior: [u64; 3], widths: &[f64]) -> bool {
        self.prior == prior
            && self.rt.len() == rt.len()
            && self.widths.len() == widths.len()
            && self.rt.iter().zip(rt).all(|(a, b)| *a == b.to_bits())
            && self
                .widths
                .iter()
                .zip(widths)
                .all(|(a, b)| *a == b.to_bits())
    }
}

#[derive(Debug, Default)]
struct GridCache {
    entries: Vec<GridProfiles>,
}

impl GridCache {
    fn get(
        &mut self,
        rt: &[f64],
        prior: &Prior,
        widths: &[f64],
        candidates: &Apexes,
    ) -> Option<&GridProfiles> {
        let key = [
            prior.t_hat.to_bits(),
            prior.s.to_bits(),
            prior.w_prior.to_bits(),
        ];
        let at = match self.entries.iter().position(|e| e.matches(rt, key, widths)) {
            Some(i) => i,
            None => {
                let finite = rt.iter().all(|t| t.is_finite());
                let mut cells = Vec::with_capacity(candidates.iter().count() * widths.len());
                let mut values = 0;
                for k in candidates.iter() {
                    for &f in widths {
                        let scans = support(rt, finite, rt[k], f * prior.w_prior);
                        values += scans.len();
                        cells.push(Cell {
                            scans,
                            at: 0,
                            sg: 0.0,
                            sgg: 0.0,
                        });
                    }
                }
                if values > GRID_CACHE_MAX_VALUES {
                    return None;
                }
                let mut e = if self.entries.len() < GRID_CACHE_ENTRIES {
                    GridProfiles::default()
                } else {
                    self.entries.remove(0)
                };
                e.rt.clear();
                e.rt.extend(rt.iter().map(|t| t.to_bits()));
                e.prior = key;
                e.widths.clear();
                e.widths.extend(widths.iter().map(|f| f.to_bits()));
                e.g.clear();
                let mut c = cells.iter_mut();
                for k in candidates.iter() {
                    let t0 = rt[k];
                    for &f in widths {
                        let w = f * prior.w_prior;
                        let cell = c.next().expect("one cell per candidate and width");
                        cell.at = e.g.len();
                        gaussian::extend(&mut e.g, &rt[cell.scans.clone()], t0, w);
                        let (mut sg, mut sgg) = (0.0, 0.0);
                        for &g in &e.g[cell.at..] {
                            sg += g;
                            sgg += g * g;
                        }
                        (cell.sg, cell.sgg) = (sg, sgg);
                    }
                }
                e.cells = cells;
                self.entries.push(e);
                self.entries.len() - 1
            }
        };
        let e = self.entries.remove(at);
        self.entries.push(e);
        self.entries.last()
    }
}

pub fn fit(
    ms1: &Ms1RunIndex,
    probe: &Probe,
    charge: u8,
    envelope: &[f32; 4],
    p: &EvidenceParams,
) -> Option<Evidence> {
    fit_with(ms1, probe, charge, envelope, p, &mut FitScratch::default())
}

pub fn fit_with(
    ms1: &Ms1RunIndex,
    probe: &Probe,
    charge: u8,
    envelope: &[f32; 4],
    p: &EvidenceParams,
    scratch: &mut FitScratch,
) -> Option<Evidence> {
    fit_on_grid(ms1, probe, charge, envelope, p, &WIDTH_GRID, scratch)
}

pub fn fit_on_grid(
    ms1: &Ms1RunIndex,
    probe: &Probe,
    charge: u8,
    envelope: &[f32; 4],
    p: &EvidenceParams,
    widths: &[f64],
    scratch: &mut FitScratch,
) -> Option<Evidence> {
    fit_on_grid_gated(ms1, probe, charge, envelope, p, widths, None, scratch)
}

#[allow(clippy::too_many_arguments)]
pub fn fit_on_grid_gated(
    ms1: &Ms1RunIndex,
    probe: &Probe,
    charge: u8,
    envelope: &[f32; 4],
    p: &EvidenceParams,
    widths: &[f64],
    gate: Option<MobilityGate>,
    scratch: &mut FitScratch,
) -> Option<Evidence> {
    debug_assert!(charge > 0, "charge must be positive");
    debug_assert!(!widths.is_empty(), "width grid must not be empty");
    let prior = Prior::of(probe);
    let half = prior.half_window();
    let scans = ms1.scan_window(prior.t_hat - half, prior.t_hat + half);
    if scans.len() < MIN_SCANS {
        return None;
    }
    let n_iso = p.n_isotopes.clamp(1, MAX_ISOTOPES);
    let detection_limit = detection_limit(ms1, scans.clone(), &mut scratch.limits);
    traces(
        ms1,
        scans,
        probe.mz,
        charge.max(1),
        n_iso,
        p.tol_ppm,
        gate,
        scratch,
    );
    let FitScratch {
        rt,
        y,
        m1,
        mz0,
        work,
        ..
    } = scratch;
    let mut ys: [&[f32]; MAX_ISOTOPES] = [&[]; MAX_ISOTOPES];
    for (slot, trace) in ys.iter_mut().zip(y.iter()).take(n_iso) {
        *slot = trace;
    }
    Some(fit_traces_on_grid(
        rt,
        &ys[..n_iso],
        m1,
        mz0,
        probe.mz,
        envelope,
        &prior,
        detection_limit,
        widths,
        work,
    ))
}

fn detection_limit(ms1: &Ms1RunIndex, scans: std::ops::Range<usize>, buf: &mut Vec<f32>) -> f64 {
    buf.clear();
    buf.extend(scans.map(|j| ms1.detection_limit(j)));
    median_f32_as_f64_mut(buf).unwrap_or(0.0)
}

#[allow(clippy::too_many_arguments)]
fn traces(
    ms1: &Ms1RunIndex,
    scans: std::ops::Range<usize>,
    mu: f64,
    charge: u8,
    n_iso: usize,
    tol_ppm: f64,
    gate: Option<MobilityGate>,
    out: &mut FitScratch,
) {
    let slices = target_slices(mu, charge, n_iso, tol_ppm);
    out.rt.clear();
    out.rt.extend(scans.clone().map(|j| ms1.rt(j)));
    out.mz0.clear();
    out.y[0].clear();
    let (lo, hi) = slices[0];
    for j in scans.clone() {
        let (mz, inten) = ms1.slice_in(j, lo, hi, gate);
        out.mz0.push(mz);
        out.y[0].push(inten);
    }
    for (i, &(lo, hi)) in slices.iter().enumerate().take(n_iso).skip(1) {
        out.y[i].clear();
        out.y[i].extend(scans.clone().map(|j| ms1.slice_in(j, lo, hi, gate).1));
    }
    let (lo, hi) = slices[n_iso];
    out.m1.clear();
    out.m1
        .extend(scans.map(|j| ms1.slice_in(j, lo, hi, gate).1));
}

const MAX_TRACES: usize = MAX_ISOTOPES + 1;

fn target_slices(mu: f64, charge: u8, n_iso: usize, tol_ppm: f64) -> [(f64, f64); MAX_TRACES] {
    let step = ISOTOPE / charge as f64;
    let slice = |target: f64| {
        let tol = target * tol_ppm * 1e-6;
        (target - tol, target + tol)
    };
    let mut out = [(0.0, 0.0); MAX_TRACES];
    out[0] = slice(mu);
    for (i, o) in out.iter_mut().enumerate().take(n_iso).skip(1) {
        *o = slice(mu + i as f64 * step);
    }
    out[n_iso] = slice(mu - step);
    out
}

#[derive(Debug, Clone, Copy)]
pub struct FitJob<'a> {
    pub probe: Probe,
    pub charge: u8,
    pub envelope: &'a [f32; 4],
    pub widths: &'a [f64],
    pub mobility: Option<MobilityGate>,
}

pub fn fit_many(ms1: &Ms1RunIndex, jobs: &[FitJob], p: &EvidenceParams) -> Vec<Option<Evidence>> {
    use rayon::prelude::*;

    let n_iso = p.n_isotopes.clamp(1, MAX_ISOTOPES);
    let mut planned = Vec::new();
    let mut fallback = Vec::new();
    for (k, job) in jobs.iter().enumerate() {
        debug_assert!(job.charge > 0, "charge must be positive");
        debug_assert!(!job.widths.is_empty(), "width grid must not be empty");
        let prior = Prior::of(&job.probe);
        let half = prior.half_window();
        let scans = ms1.scan_window(prior.t_hat - half, prior.t_hat + half);
        if scans.len() < MIN_SCANS {
            continue;
        }
        let slices = target_slices(job.probe.mz, job.charge.max(1), n_iso, p.tol_ppm);
        if slices[..=n_iso]
            .iter()
            .all(|&(lo, hi)| lo.is_finite() && hi.is_finite())
        {
            planned.push(Planned {
                job: k,
                prior,
                scans,
                slices,
            });
        } else {
            fallback.push(k);
        }
    }
    planned.sort_by_key(|w| (w.scans.start, w.job));
    let ordered: Vec<bool> = (0..ms1.len()).map(|j| mz_ascending(ms1.peaks(j))).collect();
    let chunk = planned
        .len()
        .div_ceil(SWEEP_CHUNKS_PER_THREAD * rayon::current_num_threads())
        .max(MIN_SWEEP_CHUNK);
    let swept: Vec<Vec<(usize, Evidence)>> = planned
        .par_chunks(chunk)
        .map(|c| sweep(ms1, jobs, c, n_iso, &ordered))
        .collect();
    let refits: Vec<(usize, Option<Evidence>)> = fallback
        .par_iter()
        .map_init(FitScratch::default, |scratch, &k| {
            let job = &jobs[k];
            let e = fit_on_grid_gated(
                ms1,
                &job.probe,
                job.charge,
                job.envelope,
                p,
                job.widths,
                job.mobility,
                scratch,
            );
            (k, e)
        })
        .collect();
    let mut out: Vec<Option<Evidence>> = vec![None; jobs.len()];
    for (k, e) in swept.into_iter().flatten() {
        out[k] = Some(e);
    }
    for (k, e) in refits {
        out[k] = e;
    }
    out
}

const MIN_SWEEP_CHUNK: usize = 64;
const SWEEP_CHUNKS_PER_THREAD: usize = 32;

struct Planned {
    job: usize,
    prior: Prior,
    scans: std::ops::Range<usize>,
    slices: [(f64, f64); MAX_TRACES],
}

fn mz_ascending(peaks: &[(f64, f32)]) -> bool {
    peaks.iter().all(|p| !p.0.is_nan()) && peaks.windows(2).all(|w| w[0].0 <= w[1].0)
}

#[derive(Debug, Clone, Copy)]
struct Slice {
    lo: f64,
    hi: f64,
    slot: usize,
    trace: usize,
    gate: Option<MobilityGate>,
}

#[derive(Debug, Default)]
struct Buffers {
    planned: usize,
    y: [Vec<f32>; MAX_ISOTOPES],
    m1: Vec<f32>,
    mz0: Vec<f64>,
}

fn sweep(
    ms1: &Ms1RunIndex,
    jobs: &[FitJob],
    planned: &[Planned],
    n_iso: usize,
    ordered: &[bool],
) -> Vec<(usize, Evidence)> {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;

    let mut out = Vec::with_capacity(planned.len());
    let (Some(first), Some(last)) = (
        planned.first().map(|w| w.scans.start),
        planned.iter().map(|w| w.scans.end).max(),
    ) else {
        return out;
    };
    let mut slots: Vec<Buffers> = Vec::new();
    let mut free: Vec<usize> = Vec::new();
    let mut closing: BinaryHeap<Reverse<(usize, usize)>> = BinaryHeap::new();
    let mut active: Vec<Slice> = Vec::new();
    let mut incoming: Vec<Slice> = Vec::new();
    let mut merged: Vec<Slice> = Vec::new();
    let mut closed: Vec<bool> = Vec::new();
    let mut work = Workspace::default();
    let mut limits = Vec::new();
    let mut rt = Vec::new();
    let mut next = 0;
    for (j, &ascending) in ordered.iter().enumerate().take(last).skip(first) {
        incoming.clear();
        while next < planned.len() && planned[next].scans.start == j {
            let slot = free.pop().unwrap_or_else(|| {
                slots.push(Buffers::default());
                closed.push(false);
                slots.len() - 1
            });
            let b = &mut slots[slot];
            b.planned = next;
            for v in b.y.iter_mut() {
                v.clear();
            }
            b.m1.clear();
            b.mz0.clear();
            closed[slot] = false;
            let w = &planned[next];
            let gate = jobs[w.job].mobility;
            for (trace, &(lo, hi)) in w.slices[..=n_iso].iter().enumerate() {
                incoming.push(Slice {
                    lo,
                    hi,
                    slot,
                    trace,
                    gate,
                });
            }
            closing.push(Reverse((w.scans.end, slot)));
            next += 1;
        }
        if !incoming.is_empty() {
            incoming.sort_unstable_by(|a, b| a.lo.total_cmp(&b.lo));
            merged.clear();
            let (mut a, mut b) = (0, 0);
            while a < active.len() || b < incoming.len() {
                let take_active = b == incoming.len()
                    || (a < active.len() && active[a].lo.total_cmp(&incoming[b].lo).is_le());
                if take_active {
                    merged.push(active[a]);
                    a += 1;
                } else {
                    merged.push(incoming[b]);
                    b += 1;
                }
            }
            std::mem::swap(&mut active, &mut merged);
        }

        let peaks = ms1.peaks(j);
        let mobility = ms1.mobility(j);
        let mut start = 0;
        for s in &active {
            if ascending {
                while start < peaks.len() && peaks[start].0 < s.lo {
                    start += 1;
                }
            } else {
                start = peaks.partition_point(|&(mz, _)| mz < s.lo);
            }
            let (best_mz, best) = slice_value(peaks, mobility, start, s.hi, s.gate);
            let b = &mut slots[s.slot];
            if s.trace == 0 {
                b.mz0.push(best_mz);
                b.y[0].push(best);
            } else if s.trace < n_iso {
                b.y[s.trace].push(best);
            } else {
                b.m1.push(best);
            }
        }

        let mut any_closed = false;
        while let Some(&Reverse((end, slot))) = closing.peek() {
            if end > j + 1 {
                break;
            }
            closing.pop();
            any_closed = true;
            closed[slot] = true;
            let b = &slots[slot];
            let w = &planned[b.planned];
            let job = &jobs[w.job];
            let detection_limit = detection_limit(ms1, w.scans.clone(), &mut limits);
            rt.clear();
            rt.extend(w.scans.clone().map(|j| ms1.rt(j)));
            let mut ys: [&[f32]; MAX_ISOTOPES] = [&[]; MAX_ISOTOPES];
            for (y, trace) in ys.iter_mut().zip(b.y.iter()).take(n_iso) {
                *y = trace;
            }
            let e = fit_traces_on_grid(
                &rt,
                &ys[..n_iso],
                &b.m1,
                &b.mz0,
                job.probe.mz,
                job.envelope,
                &w.prior,
                detection_limit,
                job.widths,
                &mut work,
            );
            out.push((w.job, e));
            free.push(slot);
        }
        if any_closed {
            active.retain(|s| !closed[s.slot]);
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
pub fn fit_traces(
    rt: &[f64],
    traces: &[&[f32]],
    m1: &[f32],
    mz0: &[f64],
    mu: f64,
    envelope: &[f32; 4],
    prior: &Prior,
    detection_limit: f64,
    work: &mut Workspace,
) -> Evidence {
    fit_traces_on_grid(
        rt,
        traces,
        m1,
        mz0,
        mu,
        envelope,
        prior,
        detection_limit,
        &WIDTH_GRID,
        work,
    )
}

#[allow(clippy::too_many_arguments)]
fn fit_traces_on_grid(
    rt: &[f64],
    traces: &[&[f32]],
    m1: &[f32],
    mz0: &[f64],
    mu: f64,
    envelope: &[f32; 4],
    prior: &Prior,
    detection_limit: f64,
    widths: &[f64],
    work: &mut Workspace,
) -> Evidence {
    debug_assert!(!rt.is_empty() && !traces.is_empty() && traces.len() <= MAX_ISOTOPES);
    let n_iso = traces.len().min(MAX_ISOTOPES);
    let traces = &traces[..n_iso];
    let pi = normalised(envelope);
    let limit = if detection_limit.is_finite() && detection_limit > 0.0 {
        detection_limit
    } else {
        DEFAULT_DETECTION_LIMIT
    };
    let floor = NOISE_FLOOR_OF_LIMIT * limit;
    let mut ybar = [0.0f64; MAX_ISOTOPES];
    let mut sigma = [1.0f64; MAX_ISOTOPES];
    let mut mask = 0u8;
    for i in 0..n_iso {
        ybar[i] = mean(traces[i]);
        sigma[i] = noise_sigma(traces[i], floor, &mut work.sort);
        if pi[i] > 0.0 {
            mask |= 1 << i;
        }
    }
    let all = mask;
    let fit_at = |used: u8, work: &mut Workspace| {
        let (t0, w) = best_alignment(rt, traces, &pi, &ybar, &sigma, used, prior, widths, work);
        let (sgg, proj) = profile(rt, traces, &ybar, t0, w, &mut work.g);
        (t0, w, sgg, proj)
    };
    let (mut t0, mut w, mut sgg, mut proj) = fit_at(all, work);
    let mut mask = exclude_interfered(&proj, sgg, &pi, &sigma, all);
    if mask != all {
        (t0, w, sgg, proj) = fit_at(mask, work);
        let again = exclude_interfered(&proj, sgg, &pi, &sigma, all);
        if again != all && again != mask {
            mask = again;
            (t0, w, sgg, proj) = fit_at(mask, work);
        }
    }
    let g = &work.g;
    let m1_proj = project(m1, mean(m1), g);
    features(
        rt,
        traces,
        mz0,
        mu,
        prior,
        &pi,
        &sigma,
        &proj,
        sgg,
        m1_proj,
        mask,
        t0,
        w,
        g,
        &mut work.sort,
    )
}

fn profile(
    rt: &[f64],
    traces: &[&[f32]],
    ybar: &[f64; MAX_ISOTOPES],
    t0: f64,
    w: f64,
    g: &mut Vec<f64>,
) -> (f64, [f64; MAX_ISOTOPES]) {
    g.clear();
    gaussian::extend(g, rt, t0, w);
    let gbar = g.iter().sum::<f64>() / g.len() as f64;
    let sgg = g
        .iter()
        .map(|v| (v - gbar) * (v - gbar))
        .sum::<f64>()
        .max(MIN_CENTRED_SGG);
    let mut proj = [0.0f64; MAX_ISOTOPES];
    for (p, (trace, &yb)) in proj.iter_mut().zip(traces.iter().zip(ybar)) {
        *p = project(trace, yb, g);
    }
    (sgg, proj)
}

fn normalised(envelope: &[f32; 4]) -> [f64; MAX_ISOTOPES] {
    let e = envelope.map(|v| v as f64);
    let sum: f64 = e.iter().sum();
    if !sum.is_finite() || sum <= 0.0 || e.iter().any(|&v| v < 0.0) || e[0] <= 0.0 {
        return [1.0, 0.0, 0.0, 0.0];
    }
    e.map(|v| v / sum)
}

fn project(y: &[f32], ybar: f64, g: &[f64]) -> f64 {
    y.iter()
        .zip(g)
        .map(|(&v, &gj)| (v as f64 - ybar) * gj)
        .sum()
}

fn mean(y: &[f32]) -> f64 {
    if y.is_empty() {
        0.0
    } else {
        y.iter().map(|&v| v as f64).sum::<f64>() / y.len() as f64
    }
}

fn noise_sigma(y: &[f32], floor: f64, sort: &mut Vec<f32>) -> f64 {
    sort.clear();
    sort.extend(y.windows(2).map(|w| (w[1] - w[0]).abs()));
    let mad = median_f32_as_f64_mut(sort).unwrap_or(0.0);
    (MAD_TO_SD * mad / std::f64::consts::SQRT_2).max(floor)
}

#[allow(clippy::too_many_arguments)]
fn best_alignment(
    rt: &[f64],
    traces: &[&[f32]],
    pi: &[f64; MAX_ISOTOPES],
    ybar: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    mask: u8,
    prior: &Prior,
    widths: &[f64],
    work: &mut Workspace,
) -> (f64, f64) {
    let z = &mut work.z;
    z.clear();
    z.resize(rt.len(), 0.0);
    let mut d = 0.0;
    for (i, trace) in traces.iter().enumerate() {
        if mask & (1 << i) == 0 {
            continue;
        }
        let c = pi[i] / (sigma[i] * sigma[i]);
        d += pi[i] * c;
        for (zj, &v) in z.iter_mut().zip(trace.iter()) {
            *zj += c * (v as f64 - ybar[i]);
        }
    }
    let (t_hat, s) = (prior.t_hat, prior.s);
    let lo = rt.partition_point(|&t| t < t_hat - 3.0 * s);
    let hi = rt.partition_point(|&t| t <= t_hat + 3.0 * s);
    let range = if lo < hi {
        lo..hi
    } else {
        let k = nearest(rt, t_hat);
        k..k + 1
    };
    let apexes = Apexes::of(rt, range, widths[0] * prior.w_prior);
    let n = rt.len() as f64;
    let llr_of = |num: f64, sg: f64, sgg: f64| {
        let centred = sgg - sg * sg / n;
        if centred < MIN_CENTRED_SGG {
            return f64::NEG_INFINITY;
        }
        let den = d * centred;
        let a = (num / den).max(0.0);
        a * num - 0.5 * a * a * den
    };
    let finite = rt.iter().all(|t| t.is_finite());
    let llr = |t0: f64, w: f64| {
        let on = support(rt, finite, t0, w);
        let (num, sg, sgg) = gaussian::sums(&rt[on.clone()], &z[on], t0, w);
        llr_of(num, sg, sgg)
    };
    let penalty_at = |k: usize| (rt[k] - t_hat) * (rt[k] - t_hat) / (2.0 * s * s);
    let mut best = (f64::NEG_INFINITY, apexes.range.start, usize::MAX);
    if let Some(grid) = work.grids.get(rt, prior, widths, &apexes) {
        let mut cells = grid.cells.iter();
        for k in apexes.iter() {
            let penalty = penalty_at(k);
            for i in 0..widths.len() {
                let cell = cells.next().expect("one cell per candidate and width");
                let g = &grid.g[cell.at..cell.at + cell.scans.len()];
                let mut num = 0.0;
                for (&zj, &gj) in z[cell.scans.clone()].iter().zip(g) {
                    num += zj * gj;
                }
                let objective = llr_of(num, cell.sg, cell.sgg) - penalty;
                if objective > best.0 {
                    best = (objective, k, i);
                }
            }
        }
    } else {
        for k in apexes.iter() {
            let penalty = penalty_at(k);
            for (i, &f) in widths.iter().enumerate() {
                let objective = llr(rt[k], f * prior.w_prior) - penalty;
                if objective > best.0 {
                    best = (objective, k, i);
                }
            }
        }
    }
    if apexes.stride > 1 && best.2 != usize::MAX {
        let f = widths[best.2];
        for k in apexes.around(best.1) {
            let objective = llr(rt[k], f * prior.w_prior) - penalty_at(k);
            if objective > best.0 {
                best = (objective, k, best.2);
            }
        }
    }
    let (value, k, i) = best;
    let t0 = rt[k];
    if i == usize::MAX {
        return (t0, prior.w_prior);
    }
    let (lo, hi) = width_bracket(widths, i);
    let w = golden_max(
        |log_f| llr(t0, log_f.exp() * prior.w_prior),
        lo.ln(),
        hi.ln(),
        widths[i].ln(),
        value + (t0 - t_hat) * (t0 - t_hat) / (2.0 * s * s),
    )
    .exp()
        * prior.w_prior;
    (t0, w)
}

const WIDTH_REFINE_STEPS: usize = 16;

fn width_bracket(widths: &[f64], i: usize) -> (f64, f64) {
    let f = widths[i];
    match widths.len() {
        0 | 1 => (f, f),
        n => {
            let lo = if i == 0 {
                f * f / widths[1]
            } else {
                widths[i - 1]
            };
            let hi = if i + 1 == n {
                f * f / widths[n - 2]
            } else {
                widths[i + 1]
            };
            (lo, hi)
        }
    }
}

fn golden_max<F: Fn(f64) -> f64>(f: F, mut lo: f64, mut hi: f64, x0: f64, f0: f64) -> f64 {
    const INV_PHI: f64 = 0.618_033_988_749_894_9;
    let (mut best_x, mut best_f) = (x0, f0);
    if lo.is_nan() || hi.is_nan() || hi <= lo {
        return best_x;
    }
    let mut a = hi - INV_PHI * (hi - lo);
    let mut b = lo + INV_PHI * (hi - lo);
    let (mut fa, mut fb) = (f(a), f(b));
    for _ in 0..WIDTH_REFINE_STEPS {
        for (x, v) in [(a, fa), (b, fb)] {
            if v > best_f {
                (best_x, best_f) = (x, v);
            }
        }
        if fa >= fb {
            hi = b;
            b = a;
            fb = fa;
            a = hi - INV_PHI * (hi - lo);
            fa = f(a);
        } else {
            lo = a;
            a = b;
            fa = fb;
            b = lo + INV_PHI * (hi - lo);
            fb = f(b);
        }
    }
    for (x, v) in [(a, fa), (b, fb)] {
        if v > best_f {
            (best_x, best_f) = (x, v);
        }
    }
    best_x
}

fn nearest(rt: &[f64], t: f64) -> usize {
    let i = rt.partition_point(|&x| x < t);
    if i == 0 {
        0
    } else if i == rt.len() || t - rt[i - 1] <= rt[i] - t {
        i - 1
    } else {
        i
    }
}

fn num_den(
    proj: &[f64; MAX_ISOTOPES],
    sgg: f64,
    pi: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    mask: u8,
) -> (f64, f64) {
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..MAX_ISOTOPES {
        if mask & (1 << i) != 0 {
            let inv = 1.0 / (sigma[i] * sigma[i]);
            num += pi[i] * proj[i] * inv;
            den += pi[i] * pi[i] * sgg * inv;
        }
    }
    (num, den)
}

fn exclude_interfered(
    proj: &[f64; MAX_ISOTOPES],
    sgg: f64,
    pi: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    mask: u8,
) -> u8 {
    let mut worst: Option<(usize, f64)> = None;
    for i in 1..MAX_ISOTOPES {
        if mask & (1 << i) == 0 {
            continue;
        }
        let (num, den) = num_den(proj, sgg, pi, sigma, mask & !(1 << i));
        let a_without = (num / den).max(0.0);
        let se_without2 = 1.0 / den;
        let a_i = proj[i] / (pi[i] * sgg);
        let se_i = sigma[i] / (pi[i] * sgg.sqrt());
        let score = (a_i - a_without).abs() / (se_i * se_i + se_without2).sqrt();
        if worst.is_none_or(|(_, s)| score > s) {
            worst = Some((i, score));
        }
    }
    match worst {
        Some((i, score)) if score > EXCLUSION_Z => mask & !(1 << i),
        _ => mask,
    }
}

#[allow(clippy::too_many_arguments)]
fn features(
    rt: &[f64],
    traces: &[&[f32]],
    mz0: &[f64],
    mu: f64,
    prior: &Prior,
    pi: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    proj: &[f64; MAX_ISOTOPES],
    sgg: f64,
    m1_proj: f64,
    mask: u8,
    t0: f64,
    w: f64,
    g: &[f64],
    sort: &mut Vec<f32>,
) -> Evidence {
    let (num, den) = num_den(proj, sgg, pi, sigma, mask);
    let amplitude = (num / den).max(0.0);
    let amplitude_se = 1.0 / den.sqrt();
    let llr = if amplitude > 0.0 {
        amplitude * num - 0.5 * amplitude * amplitude * den
    } else {
        0.0
    };
    let area = w * (2.0 * std::f64::consts::PI).sqrt();

    let (mut dot, mut hh, mut pp) = (0.0, 0.0, 0.0);
    for i in 0..MAX_ISOTOPES {
        if mask & (1 << i) != 0 {
            let h = proj[i] / sgg;
            dot += h * pi[i];
            hh += h * h;
            pp += pi[i] * pi[i];
        }
    }
    let envelope_cos = if hh > 0.0 && pp > 0.0 {
        dot / (hh * pp).sqrt()
    } else {
        0.0
    };

    let lo = rt.partition_point(|&t| t < t0 - 2.0 * w);
    let hi = rt.partition_point(|&t| t <= t0 + 2.0 * w);
    let n_points = hi - lo;
    let mut first_two = (0..traces.len()).filter(|&i| mask & (1 << i) != 0 && pi[i] > 0.0);
    let coherence = match (first_two.next(), first_two.next()) {
        (Some(a), Some(b)) if n_points >= 3 => {
            pearson(&traces[a][lo..hi], &g[lo..hi]).min(pearson(&traces[b][lo..hi], &g[lo..hi]))
        }
        _ => 0.0,
    };

    let m1_ratio = m1_proj.max(0.0) / sgg / (pi[0] * amplitude.max(amplitude_se));

    let lo1 = rt.partition_point(|&t| t < t0 - w);
    let hi1 = rt.partition_point(|&t| t <= t0 + w);
    let (mut sw, mut swe) = (0.0, 0.0);
    for (&mz, &v) in mz0[lo1..hi1].iter().zip(&traces[0][lo1..hi1]) {
        if mz > 0.0 && v > 0.0 && mu > 0.0 {
            sw += v as f64;
            swe += v as f64 * (mz - mu) / mu * 1e6;
        }
    }
    let mass_err_ppm = if sw > 0.0 { swe / sw } else { 0.0 };

    let mut iso_heights = [0.0f64; MAX_ISOTOPES];
    for (h, p) in iso_heights.iter_mut().zip(proj).take(traces.len()) {
        *h = p / sgg;
    }
    let tailing = tailing(rt, traces[0], g, amplitude * pi[0], t0, w);
    let reach = w.max(prior.w_prior);
    let range = (
        t0 - TRACE_WIDTHS_BEFORE * reach,
        t0 + TRACE_WIDTHS_AFTER * reach,
    );
    let (trace_quantity, trace_quantity_se) =
        trace_area(rt, traces, pi, sigma, mask, amplitude, range, g, sort);

    Evidence {
        amplitude,
        amplitude_se,
        quantity: amplitude * area,
        quantity_se: amplitude_se * area,
        llr,
        rt_z: (t0 - prior.t_hat) / prior.s,
        width_ratio: w / prior.w_prior,
        envelope_cos,
        coherence,
        m1_ratio,
        mass_err_ppm,
        n_points: n_points as u32,
        isotopes_used: mask,
        t0,
        w,
        iso_heights,
        tailing,
        trace_quantity,
        trace_quantity_se,
    }
}

#[allow(clippy::too_many_arguments)]
fn trace_area(
    rt: &[f64],
    traces: &[&[f32]],
    pi: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    mask: u8,
    a: f64,
    range: (f64, f64),
    g: &[f64],
    sort: &mut Vec<f32>,
) -> (f64, f64) {
    let lo = rt.partition_point(|&t| t < range.0);
    let hi = rt.partition_point(|&t| t <= range.1).max(lo);
    if hi - lo < 2 {
        return (0.0, 0.0);
    }
    let dt = |j: usize| {
        let before = if j > lo { rt[j] - rt[j - 1] } else { 0.0 };
        let after = if j + 1 < hi { rt[j + 1] - rt[j] } else { 0.0 };
        0.5 * (before + after)
    };
    let span = rt[hi - 1] - rt[lo];
    let sum_dt2: f64 = (lo..hi).map(|j| dt(j) * dt(j)).sum();
    let (n_left, n_right) = (lo, rt.len() - hi);
    let gbar = if g.is_empty() {
        0.0
    } else {
        g.iter().sum::<f64>() / g.len() as f64
    };
    let (mut num, mut den, mut var) = (0.0, 0.0, 0.0);
    for (i, trace) in traces.iter().enumerate() {
        if mask & (1 << i) == 0 {
            continue;
        }
        let s2 = sigma[i] * sigma[i];
        let median_var = |n: usize| std::f64::consts::FRAC_PI_2 * s2 / n as f64;
        let mut side = |part: &[f32]| {
            (part.len() >= MIN_BASELINE_SCANS).then(|| {
                sort.clear();
                sort.extend_from_slice(part);
                (
                    median_f32_as_f64_mut(sort).unwrap_or(0.0),
                    median_var(part.len()),
                )
            })
        };
        let (left, right) = (side(&trace[..lo]), side(&trace[hi..]));
        let (base, base_var) = match (left, right) {
            (Some(l), Some(r)) if (l.0 - r.0).abs() > BASELINE_SIDE_Z * (l.1 + r.1).sqrt() => {
                if l.0 <= r.0 {
                    l
                } else {
                    r
                }
            }
            (Some(_), Some(_)) => {
                sort.clear();
                sort.extend_from_slice(&trace[..lo]);
                sort.extend_from_slice(&trace[hi..]);
                (
                    median_f32_as_f64_mut(sort).unwrap_or(0.0),
                    median_var(n_left + n_right),
                )
            }
            (Some(one), None) | (None, Some(one)) => one,
            (None, None) => (mean(trace) - a * pi[i] * gbar, s2 / trace.len() as f64),
        };
        let area: f64 = (lo..hi).map(|j| dt(j) * (trace[j] as f64 - base)).sum();
        let c = pi[i] / s2;
        num += c * area;
        den += c * pi[i];
        var += c * c * (s2 * sum_dt2 + span * span * base_var);
    }
    if den > 0.0 {
        (num / den, var.sqrt() / den)
    } else {
        (0.0, 0.0)
    }
}

fn tailing(rt: &[f64], y: &[f32], g: &[f64], h0: f64, t0: f64, w: f64) -> f64 {
    if h0.is_nan() || h0 <= 0.0 || g.is_empty() {
        return 0.0;
    }
    let gbar = g.iter().sum::<f64>() / g.len() as f64;
    let b0 = mean(y) - h0 * gbar;
    let lo = rt.partition_point(|&t| t < t0 - 3.0 * w);
    let hi = rt.partition_point(|&t| t <= t0 + 3.0 * w);
    let (mut num, mut den) = (0.0, 0.0);
    for ((&t, &v), &gj) in rt[lo..hi].iter().zip(&y[lo..hi]).zip(&g[lo..hi]) {
        let fitted = h0 * gj;
        let r = v as f64 - b0 - fitted;
        if t > t0 {
            num += r;
        } else if t < t0 {
            num -= r;
        }
        den += fitted;
    }
    if den > 0.0 {
        num / den
    } else {
        0.0
    }
}

fn pearson(y: &[f32], g: &[f64]) -> f64 {
    let n = y.len() as f64;
    let my = y.iter().map(|&v| v as f64).sum::<f64>() / n;
    let mg = g.iter().sum::<f64>() / n;
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (&v, &gj) in y.iter().zip(g) {
        let (dy, dg) = (v as f64 - my, gj - mg);
        sxy += dy * dg;
        sxx += dy * dy;
        syy += dg * dg;
    }
    if sxx > 0.0 && syy > 0.0 {
        (sxy / (sxx * syy).sqrt()).clamp(-1.0, 1.0)
    } else {
        0.0
    }
}

#[cfg(test)]
#[path = "evidence/tests/evidence.rs"]
mod tests;
