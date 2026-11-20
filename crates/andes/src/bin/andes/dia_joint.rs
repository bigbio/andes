//! DIA joint model: a generative likelihood of a nominee's fragment traces around its
//! calibrated RT, marginalised over the elution offset, the peak shape and the amplitude.
//!
//! Per run the model learns, from confident targets and a sample of decoys: three elution
//! shapes (width terciles), the detection probability by expected log intensity, the intensity
//! residual, the mass error by intensity, the chance rate of a fragment by ion rank and m/z,
//! the intensity distribution of chance peaks and the amplitude prior. A nominee's score is
//! log P(traces | peptide) / P(traces | chance), summed over the offsets within the RT prior.

use model::stats::MAD_TO_SD;
use rayon::prelude::*;

use crate::dia_pass2::{Nominee, Window, CYCLES, HALF_CYCLES, TOP_IONS};

/// Run-level inputs of the joint model.
pub(crate) struct RunInput<'a> {
    pub(crate) wins: &'a [&'a Window<'a>],
    pub(crate) by_win: &'a [Vec<usize>],
    pub(crate) noms: &'a [&'a Nominee],
    pub(crate) rt_cal: &'a [f64],
    pub(crate) shift: f64,
    pub(crate) rt_sigma_min: f64,
}

pub(crate) const COLUMNS: [&str; 10] = [
    "dia_jf_llr",
    "dia_jf_llr_apex",
    "dia_jf_llr_max",
    "dia_jf_t_off",
    "dia_jf_log_a",
    "dia_jf_det",
    "dia_jf_int",
    "dia_jf_mass",
    "dia_jf_pcorr",
    "dia_jf_cor",
];

const SEARCH_PPM: f64 = 15.0;
const PRESENT_PPM: f64 = 7.0;
const MIN_PRESENT: usize = 6;
const DECOY_SAMPLE_EVERY: usize = 5;
const PROFILE_HALF: usize = 8;
const PROFILE: usize = 2 * PROFILE_HALF + 1;
const FIT_HALF: isize = 5;
const COELUTE_HALF: usize = 3;
const FAR_CYCLES: usize = 8;
const OFFSETS: usize = 28;
const SHAPES: usize = 3;
/// Width of the expected-log10-intensity bins of the detection, residual and bias tables.
const LM_STEP: f64 = 0.25;
const LM_BINS: usize = 31;
/// Log10 amplitudes of the marginalised peak height, `AMP_STEP` apart.
const AMP_GRID: [f64; 10] = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0];
const AMP_STEP: f64 = 0.5;
const MZ_EDGES: [f64; 7] = [300.0, 400.0, 500.0, 650.0, 800.0, 1000.0, 1300.0];
const MZ_BINS: usize = 8;
/// Bins of the log10 intensity of chance peaks, `INT_STEP` wide from 0 up to `INT_MAX_LOG`.
const INT_BINS: usize = 40;
const INT_STEP: f64 = 0.25;
const INT_MAX_LOG: f64 = 10.0;
/// Bins of the mass-error spread by log10 intensity, `DEV_STEP` wide from `DEV_LOG_MIN`.
const DEV_BINS: usize = 10;
const DEV_LOG_MIN: f64 = 2.0;
const DEV_STEP: f64 = 0.5;
/// Ions of a bin below which its mass-error spread takes the smallest measured one.
const DEV_MIN_IONS: usize = 200;
/// Floor of a bin's mass-error spread, and its value when no bin is measured (ppm).
const DEV_SD_MIN_PPM: f64 = 0.5;
const DEV_SD_FALLBACK_PPM: f64 = 2.0;
/// Share of a present fragment's mass error from a uniform outlier density.
const MASS_EPS: f64 = 0.02;
/// Width (ppm) of that uniform density: the presence window either side of the ion.
const MASS_UNIFORM_PPM: f64 = 2.0 * PRESENT_PPM;
const HALF_LOG_2PI: f64 = 0.918_938_533_204_672_7;
const SQRT_2PI: f64 = 2.506_628_274_631_000_5;
/// Pseudo-counts of the chance rate of a fragment (hits, trials), and its bounds per run and
/// per (ion rank, m/z bin) cell.
const P0_PRIOR_HITS: f64 = 0.5;
const P0_PRIOR_TRIALS: f64 = 50.0;
const P0_RUN_MIN: f64 = 1e-3;
const P0_RUN_MAX: f64 = 0.3;
const P0_CELL_MIN: f64 = 1e-4;
const P0_CELL_MAX: f64 = 0.5;
/// Trials of an expected-intensity bin below which its detection probability is the chance rate.
const DET_MIN_TRIALS: f64 = 30.0;
/// Upper bound of a detection probability.
const DET_MAX: f64 = 0.995;
/// Detections of a bin below which its intensity residual takes the default spread and no bias.
const RES_MIN_HITS: f64 = 20.0;
const RES_VAR_MIN: f64 = 0.01;
const RES_SD_DEFAULT: f64 = 0.5;
/// Percentiles of the central width that split the confident profiles into the three shapes.
const SHAPE_SPLIT_PCT: [f64; 2] = [33.3, 66.7];
/// Floors of a relative predicted intensity and of a shape's relative height before their log.
const REL_INT_MIN: f64 = 1e-4;
const SHAPE_MIN: f64 = 1e-6;
/// Added to denominators and maxima that may be zero.
const DIV_EPS: f64 = 1e-9;
/// Floors of the cycle time (min) and of the RT prior's spread in cycles.
const CYCLE_MIN_FLOOR: f64 = 1e-6;
const PRIOR_SD_MIN_CYCLES: f64 = 1e-3;

/// One nominee's traces: intensity and ppm error of the nearest peak per cycle and ion.
struct Traces {
    v: Vec<[f32; TOP_IONS]>,
    d: Vec<[f32; TOP_IONS]>,
    valid_ion: [bool; TOP_IONS],
    logr: [f64; TOP_IONS],
    mzb: [usize; TOP_IONS],
    wn: [f64; TOP_IONS],
}

impl Traces {
    fn present(&self, k: usize, j: usize) -> bool {
        self.valid_ion[j] && self.v[k][j] > 0.0 && f64::from(self.d[k][j]).abs() <= PRESENT_PPM
    }
}

fn traces(w: &Window, n: &Nominee, rt_cal: f64, shift: f64) -> Traces {
    let c0 = w.rt.partition_point(|&t| t < rt_cal) as isize - HALF_CYCLES as isize;
    let mut v = vec![[0f32; TOP_IONS]; CYCLES];
    let mut d = vec![[f32::NAN; TOP_IONS]; CYCLES];
    let mut valid_ion = [false; TOP_IONS];
    let mut logr = [0f64; TOP_IONS];
    let mut mzb = [0usize; TOP_IONS];
    let mut wn = [0f64; TOP_IONS];
    let wmax = n.ions.iter().map(|i| f64::from(i.1)).fold(0.0, f64::max).max(DIV_EPS);
    for (j, &(mz, p)) in n.ions.iter().enumerate().take(TOP_IONS) {
        valid_ion[j] = mz > 0.0;
        wn[j] = f64::from(p) / wmax;
        logr[j] = wn[j].max(REL_INT_MIN).log10();
        mzb[j] = MZ_EDGES.iter().filter(|&&e| mz >= e).count();
    }
    for k in 0..CYCLES {
        let c = c0 + k as isize;
        if c < 0 || c as usize >= w.rt.len() {
            continue;
        }
        let peaks = w.scans[c as usize];
        for (j, &ok) in valid_ion.iter().enumerate() {
            if !ok {
                continue;
            }
            let target = n.ions[j].0 * shift;
            let i = peaks.partition_point(|p| p.0 < target);
            let near = [i.checked_sub(1), (i < peaks.len()).then_some(i)]
                .into_iter()
                .flatten()
                .min_by(|&a, &b| {
                    (peaks[a].0 - target)
                        .abs()
                        .total_cmp(&(peaks[b].0 - target).abs())
                });
            if let Some(p) = near {
                let ppm = (peaks[p].0 - target) / target * 1e6;
                if ppm.abs() <= SEARCH_PPM {
                    v[k][j] = peaks[p].1;
                    d[k][j] = ppm as f32;
                }
            }
        }
    }
    Traces {
        v,
        d,
        valid_ion,
        logr,
        mzb,
        wn,
    }
}

/// Apex cycle (weighted presence), co-elution with the other fragments, and the summed
/// profile over +-8 cycles of the apex.
struct Summary {
    cs: usize,
    npres: usize,
    cor: f64,
    prof: [f64; PROFILE],
}

fn clampc(c: isize) -> usize {
    c.clamp(0, CYCLES as isize - 1) as usize
}

fn summarise(t: &Traces) -> Summary {
    let vp = |k: usize, j: usize| if t.present(k, j) { f64::from(t.v[k][j]) } else { 0.0 };
    let mut cs = 0usize;
    let mut best = f64::NEG_INFINITY;
    for k in 0..CYCLES {
        let s: f64 = (0..TOP_IONS).map(|j| vp(k, j).sqrt() * t.wn[j].sqrt()).sum();
        if s > best {
            best = s;
            cs = k;
        }
    }
    let span: Vec<usize> = (-(COELUTE_HALF as isize)..=COELUTE_HALF as isize)
        .map(|d| clampc(cs as isize + d))
        .collect();
    let mut npres = 0usize;
    let mut cor_sum = 0.0;
    for j in (0..TOP_IONS).filter(|&j| t.present(cs, j)) {
        npres += 1;
        let a: Vec<f64> = span.iter().map(|&k| vp(k, j)).collect();
        let b: Vec<f64> = span
            .iter()
            .map(|&k| (0..TOP_IONS).filter(|&i| i != j).map(|i| vp(k, i)).sum())
            .collect();
        cor_sum += pearson(&a, &b);
    }
    let mut prof = [0f64; PROFILE];
    for (p, slot) in prof.iter_mut().enumerate() {
        let k = clampc(cs as isize + p as isize - PROFILE_HALF as isize);
        *slot = (0..TOP_IONS).map(|j| vp(k, j)).sum();
    }
    Summary {
        cs,
        npres,
        cor: cor_sum / npres.max(1) as f64,
        prof,
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut s, mut sa, mut sb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        s += (x - ma) * (y - mb);
        sa += (x - ma) * (x - ma);
        sb += (y - mb) * (y - mb);
    }
    s / ((sa * sb).sqrt() + DIV_EPS)
}

fn normalised(p: &[f64; PROFILE]) -> [f64; PROFILE] {
    let m = p.iter().copied().fold(0.0, f64::max).max(DIV_EPS);
    let mut o = [0.0; PROFILE];
    for (a, b) in o.iter_mut().zip(p) {
        *a = b / m;
    }
    o
}

fn central_width(p: &[f64; PROFILE]) -> f64 {
    let lo = (PROFILE_HALF as isize - FIT_HALF) as usize;
    let hi = (PROFILE_HALF as isize + FIT_HALF) as usize;
    p[lo..=hi].iter().sum()
}

fn percentile(v: &mut [f64], q: f64) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(f64::total_cmp);
    let pos = q / 100.0 * (v.len() - 1) as f64;
    let (i, f) = (pos.floor() as usize, pos - pos.floor());
    v[i] + f * (v[(i + 1).min(v.len() - 1)] - v[i])
}

fn median(v: &mut [f64]) -> f64 {
    percentile(v, 50.0)
}

/// Per-run fitted model.
struct Model {
    log_shape: [[f64; PROFILE]; SHAPES],
    width_edges: [f64; 2],
    e_prof: [f64; PROFILE],
    lpd: [f64; LM_BINS],
    l1mpd: [f64; LM_BINS],
    ls: [f64; LM_BINS],
    bias: [f64; LM_BINS],
    sd_dev: [f64; DEV_BINS],
    p0: [[f64; MZ_BINS]; TOP_IONS],
    f0: [f64; INT_BINS],
    log_pa: [f64; AMP_GRID.len()],
}

fn lm_bin(lm: f64) -> usize {
    if !lm.is_finite() || lm < 0.0 {
        return 0;
    }
    ((lm / LM_STEP) as usize).min(LM_BINS - 1)
}

/// Bin of a chance peak's log10 intensity.
fn int_bin(l: f64) -> usize {
    ((l / INT_STEP) as usize).min(INT_BINS - 1)
}

/// Mass-error bin of a present peak's log10 intensity.
fn dev_bin(l: f64) -> usize {
    (((l - DEV_LOG_MIN) / DEV_STEP).max(0.0) as usize).min(DEV_BINS - 1)
}

fn l10(v: f32) -> f64 {
    f64::from(v).max(1.0).log10()
}

/// Fit sufficient statistics, merged over windows.
#[derive(Clone)]
struct Stats {
    n: [f64; LM_BINS],
    k: [f64; LM_BINS],
    r1: [f64; LM_BINS],
    r2: [f64; LM_BINS],
    prof: [f64; PROFILE],
    nprof: f64,
    la: Vec<f64>,
    dev: Vec<Vec<f64>>,
    p0n: [[f64; MZ_BINS]; TOP_IONS],
    p0k: [[f64; MZ_BINS]; TOP_IONS],
    f0: [f64; INT_BINS],
}

impl Default for Stats {
    fn default() -> Self {
        Self {
            n: [0.0; LM_BINS],
            k: [0.0; LM_BINS],
            r1: [0.0; LM_BINS],
            r2: [0.0; LM_BINS],
            prof: [0.0; PROFILE],
            nprof: 0.0,
            la: Vec::new(),
            dev: vec![Vec::new(); DEV_BINS],
            p0n: [[0.0; MZ_BINS]; TOP_IONS],
            p0k: [[0.0; MZ_BINS]; TOP_IONS],
            f0: [0.0; INT_BINS],
        }
    }
}

impl Stats {
    fn add(mut self, o: Stats) -> Stats {
        for b in 0..LM_BINS {
            self.n[b] += o.n[b];
            self.k[b] += o.k[b];
            self.r1[b] += o.r1[b];
            self.r2[b] += o.r2[b];
        }
        for p in 0..PROFILE {
            self.prof[p] += o.prof[p];
        }
        self.nprof += o.nprof;
        self.la.extend(o.la);
        for (a, b) in self.dev.iter_mut().zip(o.dev) {
            a.extend(b);
        }
        for j in 0..TOP_IONS {
            for b in 0..MZ_BINS {
                self.p0n[j][b] += o.p0n[j][b];
                self.p0k[j][b] += o.p0k[j][b];
            }
        }
        for b in 0..INT_BINS {
            self.f0[b] += o.f0[b];
        }
        self
    }
}

fn shapes(confident: &[[f64; PROFILE]]) -> ([[f64; PROFILE]; SHAPES], [f64; 2]) {
    let mut widths: Vec<f64> = confident.iter().map(central_width).collect();
    let edges = [
        percentile(&mut widths, SHAPE_SPLIT_PCT[0]),
        percentile(&mut widths, SHAPE_SPLIT_PCT[1]),
    ];
    let mut sum = [[0f64; PROFILE]; SHAPES];
    let mut cnt = [0f64; SHAPES];
    for p in confident {
        let g = edges.iter().filter(|&&e| central_width(p) >= e).count();
        for (a, b) in sum[g].iter_mut().zip(p) {
            *a += b;
        }
        cnt[g] += 1.0;
    }
    let mut log_shape = [[0f64; PROFILE]; SHAPES];
    for g in 0..SHAPES {
        let mean: Vec<f64> = sum[g].iter().map(|s| s / cnt[g].max(1.0)).collect();
        let m = mean.iter().copied().fold(0.0, f64::max).max(DIV_EPS);
        for (ls, &mp) in log_shape[g].iter_mut().zip(&mean) {
            *ls = (mp / m).max(SHAPE_MIN).log10();
        }
    }
    (log_shape, edges)
}

fn shape_of(edges: &[f64; 2], prof: &[f64; PROFILE]) -> usize {
    let w = central_width(&normalised(prof));
    edges.iter().filter(|&&e| w >= e).count()
}

fn fit_row(
    st: &mut Stats,
    t: &Traces,
    s: &Summary,
    confident: bool,
    decoy: bool,
    log_shape: &[[f64; PROFILE]; SHAPES],
    edges: &[f64; 2],
) {
    if confident && s.npres >= MIN_PRESENT {
        let mut diffs: Vec<f64> = (0..TOP_IONS)
            .filter(|&j| t.present(s.cs, j))
            .map(|j| l10(t.v[s.cs][j]) - t.logr[j])
            .collect();
        let la = median(&mut diffs);
        let g = shape_of(edges, &s.prof);
        for d in -FIT_HALF..=FIT_HALF {
            let k = clampc(s.cs as isize + d);
            for j in (0..TOP_IONS).filter(|&j| t.valid_ion[j]) {
                let lm = la + t.logr[j] + log_shape[g][(PROFILE_HALF as isize + d) as usize];
                let b = lm_bin(lm);
                st.n[b] += 1.0;
                if t.present(k, j) {
                    let r = l10(t.v[k][j]) - lm;
                    st.k[b] += 1.0;
                    st.r1[b] += r;
                    st.r2[b] += r * r;
                }
            }
        }
        let p = normalised(&s.prof);
        for (a, b) in st.prof.iter_mut().zip(p) {
            *a += b;
        }
        st.nprof += 1.0;
        st.la.push(la);
        for j in (0..TOP_IONS).filter(|&j| t.present(s.cs, j)) {
            let lv = l10(t.v[s.cs][j]);
            let q = dev_bin(lv);
            st.dev[q].push(f64::from(t.d[s.cs][j]));
        }
    }
    if decoy {
        for k in (0..CYCLES).filter(|&k| k.abs_diff(s.cs) > FAR_CYCLES) {
            for j in (0..TOP_IONS).filter(|&j| t.valid_ion[j]) {
                st.p0n[j][t.mzb[j]] += 1.0;
                if t.present(k, j) {
                    st.p0k[j][t.mzb[j]] += 1.0;
                    let lv = l10(t.v[k][j]);
                    if lv <= INT_MAX_LOG {
                        st.f0[int_bin(lv)] += 1.0;
                    }
                }
            }
        }
    }
}

fn model(st: Stats, log_shape: [[f64; PROFILE]; SHAPES], edges: [f64; 2]) -> Model {
    let p0k: f64 = st.p0k.iter().flatten().sum();
    let p0n: f64 = st.p0n.iter().flatten().sum();
    let p0m = ((p0k + P0_PRIOR_HITS) / (p0n + P0_PRIOR_TRIALS)).clamp(P0_RUN_MIN, P0_RUN_MAX);
    let mut pdet = [0f64; LM_BINS];
    let mut run = 0f64;
    for (b, pd) in pdet.iter_mut().enumerate() {
        let v = if st.n[b] < DET_MIN_TRIALS {
            p0m
        } else {
            ((st.k[b] + 0.5) / (st.n[b] + 1.0)).clamp(p0m, DET_MAX)
        };
        run = run.max(v);
        *pd = run;
    }
    let mut m = Model {
        log_shape,
        width_edges: edges,
        e_prof: [0.0; PROFILE],
        lpd: [0.0; LM_BINS],
        l1mpd: [0.0; LM_BINS],
        ls: [0.0; LM_BINS],
        bias: [0.0; LM_BINS],
        sd_dev: [0.0; DEV_BINS],
        p0: [[0.0; MZ_BINS]; TOP_IONS],
        f0: [0.0; INT_BINS],
        log_pa: [0.0; AMP_GRID.len()],
    };
    for (b, &pd) in pdet.iter().enumerate() {
        let cnt = st.k[b].max(1.0);
        let (bias, sres) = if st.k[b] > RES_MIN_HITS {
            let bias = st.r1[b] / cnt;
            (bias, (st.r2[b] / cnt - bias * bias).max(RES_VAR_MIN).sqrt())
        } else {
            (0.0, RES_SD_DEFAULT)
        };
        m.lpd[b] = pd.ln();
        m.l1mpd[b] = (1.0 - pd).ln();
        m.ls[b] = sres.ln();
        m.bias[b] = bias;
    }
    let pm = st.prof.iter().map(|p| p / st.nprof.max(1.0)).fold(0.0, f64::max).max(DIV_EPS);
    for p in 0..PROFILE {
        m.e_prof[p] = st.prof[p] / st.nprof.max(1.0) / pm;
    }
    let mut sd: Vec<Option<f64>> = Vec::with_capacity(DEV_BINS);
    for mut v in st.dev {
        if v.len() > DEV_MIN_IONS {
            let med = median(&mut v);
            let mut dev: Vec<f64> = v.iter().map(|x| (x - med).abs()).collect();
            sd.push(Some((MAD_TO_SD * median(&mut dev)).max(DEV_SD_MIN_PPM)));
        } else {
            sd.push(None);
        }
    }
    let floor = sd.iter().flatten().copied().fold(f64::INFINITY, f64::min);
    let floor = if floor.is_finite() { floor } else { DEV_SD_FALLBACK_PPM };
    for (q, s) in sd.into_iter().enumerate() {
        m.sd_dev[q] = s.unwrap_or(floor);
    }
    for j in 0..TOP_IONS {
        for b in 0..MZ_BINS {
            m.p0[j][b] = ((st.p0k[j][b] + P0_PRIOR_HITS) / (st.p0n[j][b] + P0_PRIOR_TRIALS))
                .clamp(P0_CELL_MIN, P0_CELL_MAX);
        }
    }
    let ftot: f64 = st.f0.iter().map(|c| c + 1.0).sum();
    for b in 0..INT_BINS {
        m.f0[b] = (st.f0[b] + 1.0) / ftot / INT_STEP;
    }
    let mut pa = [1f64; AMP_GRID.len()];
    let lo = AMP_GRID[0] - 0.5 * AMP_STEP;
    let hi = AMP_GRID[AMP_GRID.len() - 1] + 0.5 * AMP_STEP;
    for &la in st.la.iter().filter(|&&la| la >= lo && la <= hi) {
        pa[(((la - lo) / AMP_STEP) as usize).min(AMP_GRID.len() - 1)] += 1.0;
    }
    let tot: f64 = pa.iter().sum();
    for (lp, &c) in m.log_pa.iter_mut().zip(&pa) {
        *lp = (c / tot).ln();
    }
    m
}

/// The ten joint-model columns of one nominee.
fn score(m: &Model, t: &Traces, s: &Summary, log_prior: &[f64]) -> [f32; COLUMNS.len()] {
    let mut base = vec![[0f64; TOP_IONS]; CYCLES];
    let mut ly = vec![[0f64; TOP_IONS]; CYCLES];
    for k in 0..CYCLES {
        for j in (0..TOP_IONS).filter(|&j| t.valid_ion[j]) {
            let lp0 = m.p0[j][t.mzb[j]];
            if t.present(k, j) {
                let l = l10(t.v[k][j]);
                ly[k][j] = l;
                let cint = -m.f0[int_bin(l)].ln();
                let q = dev_bin(l);
                let sdd = m.sd_dev[q];
                let dn = f64::from(t.d[k][j]);
                let cmass = ((1.0 - MASS_EPS) * (-0.5 * (dn / sdd).powi(2)).exp()
                    / (sdd * SQRT_2PI)
                    + MASS_EPS / MASS_UNIFORM_PPM)
                    .ln()
                    + MASS_UNIFORM_PPM.ln();
                base[k][j] = -lp0.ln() + cint + cmass - HALF_LOG_2PI;
            } else {
                base[k][j] = -(1.0 - lp0).ln();
            }
        }
    }
    let ln_shapes = (SHAPES as f64).ln();
    let eval = |wd: isize| -> (f64, f64, usize, usize, usize) {
        let mut best = f64::NEG_INFINITY;
        let (mut bt, mut ba, mut bs) = (0usize, 0usize, 0usize);
        let mut acc = f64::NEG_INFINITY;
        for (ti, lpi) in log_prior.iter().enumerate() {
            let c0 = HALF_CYCLES as isize - OFFSETS as isize + ti as isize;
            for sh in 0..SHAPES {
                for (a, &ga) in AMP_GRID.iter().enumerate() {
                    let mut tot = 0.0;
                    for d in -wd..=wd {
                        let c = c0 + d;
                        if c < 0 || c >= CYCLES as isize {
                            continue;
                        }
                        let c = c as usize;
                        let lef = m.log_shape[sh][(PROFILE_HALF as isize + d) as usize];
                        for j in (0..TOP_IONS).filter(|&j| t.valid_ion[j]) {
                            let lm = ga + t.logr[j] + lef;
                            let b = lm_bin(lm);
                            if t.present(c, j) {
                                let z = (ly[c][j] - lm - m.bias[b]) / m.ls[b].exp();
                                tot += base[c][j] + m.lpd[b] - m.ls[b] - 0.5 * z * z;
                            } else {
                                tot += base[c][j] + m.l1mpd[b];
                            }
                        }
                    }
                    let s = tot + lpi + m.log_pa[a] - ln_shapes;
                    if s > best {
                        best = s;
                        bt = ti;
                        ba = a;
                        bs = sh;
                    }
                    let hi = s.max(acc);
                    acc = hi + ((acc - hi).exp() + (s - hi).exp()).ln();
                }
            }
        }
        (acc, best, bt, ba, bs)
    };
    let (llr, llr_max, bt, ba, bs) = eval(FIT_HALF);
    let (llr_apex, _, _, _, _) = eval(0);
    let c0 = HALF_CYCLES as isize - OFFSETS as isize + bt as isize;
    let (mut det, mut int, mut mass) = (0.0, 0.0, 0.0);
    for d in -FIT_HALF..=FIT_HALF {
        let c = c0 + d;
        if c < 0 || c >= CYCLES as isize {
            continue;
        }
        let c = c as usize;
        let lef = m.log_shape[bs][(PROFILE_HALF as isize + d) as usize];
        for j in (0..TOP_IONS).filter(|&j| t.valid_ion[j]) {
            let lm = AMP_GRID[ba] + t.logr[j] + lef;
            let b = lm_bin(lm);
            let lp0 = m.p0[j][t.mzb[j]];
            if t.present(c, j) {
                let z = (ly[c][j] - lm - m.bias[b]) / m.ls[b].exp();
                let l = ly[c][j];
                let cint = -m.f0[int_bin(l)].ln();
                det += -lp0.ln() + m.lpd[b];
                int += cint - m.ls[b] - 0.5 * z * z - HALF_LOG_2PI;
                mass += base[c][j] - (-lp0.ln() + cint - HALF_LOG_2PI);
            } else {
                det += -(1.0 - lp0).ln() + m.l1mpd[b];
            }
        }
    }
    let pcorr = pearson(&s.prof, &m.e_prof);
    [
        llr as f32,
        llr_apex as f32,
        llr_max as f32,
        (bt as isize - OFFSETS as isize) as f32,
        AMP_GRID[ba] as f32,
        det as f32,
        int as f32,
        mass as f32,
        pcorr as f32,
        s.cor as f32,
    ]
}

/// RT prior over the elution offsets: Gaussian with the run's RT calibration spread.
fn log_prior(sigma_min: f64, cycle_min: f64) -> Vec<f64> {
    let sg = (sigma_min / cycle_min.max(CYCLE_MIN_FLOOR)).max(PRIOR_SD_MIN_CYCLES);
    let lp: Vec<f64> = (-(OFFSETS as isize)..=OFFSETS as isize)
        .map(|o| -0.5 * (o as f64 / sg).powi(2))
        .collect();
    let hi = lp.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let lse = hi + lp.iter().map(|v| (v - hi).exp()).sum::<f64>().ln();
    lp.into_iter().map(|v| v - lse).collect()
}

/// Joint-model columns of one run's nominees, grouped by window (`by_win[w]` = indices into
/// `noms` and `rt_cal`). Returns one row of `COLUMNS` per nominee index.
pub(crate) fn run(inp: &RunInput) -> Vec<(usize, [f32; COLUMNS.len()])> {
    let (wins, by_win, noms, rt_cal, shift, rt_sigma_min) = (
        inp.wins,
        inp.by_win,
        inp.noms,
        inp.rt_cal,
        inp.shift,
        inp.rt_sigma_min,
    );
    let confident = |n: &Nominee| crate::dia_pass2::confident(n);
    let mut decoy_rank = vec![usize::MAX; noms.len()];
    let mut nd = 0usize;
    for (k, n) in noms.iter().enumerate() {
        if crate::dia_pass2::null_member(n) {
            decoy_rank[k] = nd;
            nd += 1;
        }
    }
    let sampled = |k: usize| decoy_rank[k] != usize::MAX && decoy_rank[k] % DECOY_SAMPLE_EVERY == 0;
    let profiles: Vec<[f64; PROFILE]> = by_win
        .par_iter()
        .enumerate()
        .flat_map_iter(|(w, ks)| {
            ks.iter()
                .filter(|&&k| confident(noms[k]))
                .filter_map(|&k| {
                    let t = traces(wins[w], noms[k], rt_cal[k], shift);
                    let s = summarise(&t);
                    (s.npres >= MIN_PRESENT).then(|| normalised(&s.prof))
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let (log_shape, edges) = shapes(&profiles);
    let per_win: Vec<Stats> = by_win
        .par_iter()
        .enumerate()
        .map(|(w, ks)| {
            let mut st = Stats::default();
            for &k in ks.iter().filter(|&&k| confident(noms[k]) || sampled(k)) {
                let t = traces(wins[w], noms[k], rt_cal[k], shift);
                let s = summarise(&t);
                fit_row(&mut st, &t, &s, confident(noms[k]), noms[k].decoy, &log_shape, &edges);
            }
            st
        })
        .collect();
    let stats = per_win.into_iter().fold(Stats::default(), Stats::add);
    let n_conf = stats.la.len();
    let m = model(stats, log_shape, edges);
    log::info!(
        target: "dia",
        "joint model: {n_conf} confident targets with >= {MIN_PRESENT} fragments, shape width edges \
         {:.2} / {:.2}, RT prior sigma {rt_sigma_min:.3} min",
        m.width_edges[0],
        m.width_edges[1]
    );
    by_win
        .par_iter()
        .enumerate()
        .flat_map_iter(|(w, ks)| {
            let dt = {
                let mut d: Vec<f64> = wins[w].rt.windows(2).map(|p| p[1] - p[0]).collect();
                median(&mut d)
            };
            let prior = log_prior(rt_sigma_min, dt);
            ks.iter()
                .map(|&k| {
                    let t = traces(wins[w], noms[k], rt_cal[k], shift);
                    let s = summarise(&t);
                    (k, score(&m, &t, &s, &prior))
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nominee(decoy: bool, apex_rt: f64, ions: Vec<(f64, f32)>) -> Nominee {
        Nominee {
            run: 0,
            charge: 2,
            decoy,
            prec_mz: 500.0,
            pred_rt: apex_rt,
            rt_missing: false,
            raw_rt: f64::NAN,
            apex_rt,
            dia_llr: if decoy { 0.0 } else { 20.0 },
            mobility: f64::NAN,
            ions,
            seq: "PEPTIDEK".into(),
            peptide: "-.PEPTIDEK.-".into(),
            proteins: "P1".into(),
            spec_id: "s".into(),
            xfer: false,
        }
    }

    /// Windows of 160 cycles with co-eluting targets of non-integer intensities and decoys
    /// without signal: (scans per window, RT per window, nominees, window of each nominee).
    #[allow(clippy::type_complexity)]
    fn study() -> (Vec<Vec<Vec<(f64, f32)>>>, Vec<f64>, Vec<Nominee>, Vec<Vec<usize>>) {
        let n_win = 9;
        let n_cyc = 160;
        let rt: Vec<f64> = (0..n_cyc).map(|c| c as f64 * 0.021).collect();
        let mut scans = vec![vec![Vec::new(); n_cyc]; n_win];
        let mut noms = Vec::new();
        let mut by_win = vec![Vec::new(); n_win];
        for w in 0..n_win {
            for i in 0..24 {
                let decoy = i % 4 == 3;
                let apex = 30 + (i * 37 + w * 11) % 100;
                let ions: Vec<(f64, f32)> = (0..TOP_IONS)
                    .map(|j| {
                        let mz = 200.0 + 3.17 * i as f64 + 61.3 * j as f64 + 0.37 * w as f64;
                        (mz, 1.0 / (1.0 + 0.71 * j as f32))
                    })
                    .collect();
                if !decoy {
                    for c in apex.saturating_sub(9)..(apex + 10).min(n_cyc) {
                        let d = c as f64 - apex as f64;
                        let h = 1.0e5 * (1.0 + 0.13 * i as f64) * (-0.5 * (d / 2.7).powi(2)).exp();
                        for (j, &(mz, p)) in ions.iter().enumerate() {
                            let ppm = 1.3 * ((i + j + c) % 5) as f64 - 2.6;
                            scans[w][c].push((mz * (1.0 + ppm * 1e-6), (h * f64::from(p)) as f32 + 17.3));
                        }
                    }
                }
                by_win[w].push(noms.len());
                noms.push(nominee(decoy, rt[apex], ions));
            }
            for (c, sc) in scans[w].iter_mut().enumerate() {
                for k in 0..40 {
                    sc.push((150.0 + 23.71 * k as f64 + 0.011 * c as f64, 50.0 + 7.3 * k as f32));
                }
                sc.sort_by(|a, b| a.0.total_cmp(&b.0));
            }
        }
        (scans, rt, noms, by_win)
    }

    #[test]
    fn joint_columns_do_not_depend_on_the_thread_count() {
        let (scans, rt, noms, by_win) = study();
        let wins: Vec<Window> = scans
            .iter()
            .map(|w| Window::for_test(rt.clone(), w.iter().map(|s| s.as_slice()).collect()))
            .collect();
        let win_refs: Vec<&Window> = wins.iter().collect();
        let refs: Vec<&Nominee> = noms.iter().collect();
        let rt_cal: Vec<f64> = noms.iter().map(|n| n.pred_rt).collect();
        let inp = RunInput {
            wins: &win_refs,
            by_win: &by_win,
            noms: &refs,
            rt_cal: &rt_cal,
            shift: 1.0,
            rt_sigma_min: 0.2,
        };
        let bits = |threads: usize| -> Vec<(usize, Vec<u32>)> {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            let mut out = pool.install(|| run(&inp));
            out.sort_by_key(|r| r.0);
            out.into_iter()
                .map(|(k, v)| (k, v.iter().map(|x| x.to_bits()).collect()))
                .collect()
        };
        let one = bits(1);
        assert_eq!(one.len(), noms.len());
        assert!(one.iter().any(|(_, v)| f32::from_bits(v[0]) > 0.0), "targets score");
        for threads in [2, 5, 8] {
            assert_eq!(bits(threads), one, "{threads} threads");
        }
    }

    #[test]
    fn log_prior_is_normalised_and_centred() {
        let lp = log_prior(0.27, 0.02);
        let total: f64 = lp.iter().map(|v| v.exp()).sum();
        assert!((total - 1.0).abs() < 1e-9);
        let best = lp
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap()
            .0;
        assert_eq!(best, OFFSETS);
    }

    #[test]
    fn central_width_counts_the_middle_cycles() {
        let mut p = [0f64; PROFILE];
        p[PROFILE_HALF] = 1.0;
        p[0] = 1.0;
        assert_eq!(central_width(&p), 1.0);
    }

    #[test]
    fn lm_bins_clip_at_both_ends() {
        assert_eq!(lm_bin(-3.0), 0);
        assert_eq!(lm_bin(f64::NEG_INFINITY), 0);
        assert_eq!(lm_bin(0.3), 1);
        assert_eq!(lm_bin(99.0), LM_BINS - 1);
    }
}
