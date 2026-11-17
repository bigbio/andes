use model::mass::ISOTOPE;

use crate::decoys::Probe;
use crate::gaussian;
use crate::ms1_index::{slice_value, MobilityGate, Ms1RunIndex};

pub const MAX_ISOTOPES: usize = 4;
pub const MIN_SCANS: usize = 5;
pub const WIDTH_GRID: [f64; 3] = [0.7, 1.0, 1.4];
const MAD_SCALE: f64 = 1.4826;
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
    buf.sort_unstable_by(f32::total_cmp);
    median_sorted(buf)
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
    sort.sort_unstable_by(f32::total_cmp);
    let mad = median_sorted(sort);
    (MAD_SCALE * mad / std::f64::consts::SQRT_2).max(floor)
}

fn median_sorted(v: &[f32]) -> f64 {
    match v.len() {
        0 => 0.0,
        n if n % 2 == 1 => v[n / 2] as f64,
        n => 0.5 * (v[n / 2 - 1] as f64 + v[n / 2] as f64),
    }
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
                (median_unsorted(sort), median_var(part.len()))
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
                (median_unsorted(sort), median_var(n_left + n_right))
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

fn median_unsorted(v: &mut [f32]) -> f64 {
    let n = v.len();
    if n == 0 {
        return 0.0;
    }
    let (left, mid, _) = v.select_nth_unstable_by(n / 2, f32::total_cmp);
    let upper = *mid as f64;
    if n % 2 == 1 {
        upper
    } else {
        let lower = left.iter().copied().fold(f32::NEG_INFINITY, f32::max) as f64;
        0.5 * (lower + upper)
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
mod tests {
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
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = self.0;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
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
        assert_eq!(median_unsorted(&mut []), 0.0);
        assert_eq!(median_unsorted(&mut [3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median_unsorted(&mut [4.0, 1.0, 3.0, 2.0]), 2.5);
    }
}
