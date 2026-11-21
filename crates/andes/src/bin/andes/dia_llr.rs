//! Per-scan fragment evidence of the DIA candidate pool: a likelihood ratio of the predicted
//! b/y ions' presence against chance, a same-apex term over the neighbouring cycles, and two
//! explaining-away variants over the rows that share a scan.
//!
//! For each pool row (one candidate of one DIA scan) every predicted b/y ion of charge at most the
//! precursor's, inside the scan's peak range, is looked up in the scan (`x` = 1 when a peak lies
//! within the fragment tolerance). Then
//!
//! `LLR = Σ x·ln(q/p0) + (1 − x)·ln((1 − q)/(1 − p0))`
//!
//! with `q = P(present | true)` by the ion's predicted-intensity rank and the precursor's MS1
//! intensity quartile, learned from anchors (best row of a scan at RawScore q ≤ 1 %, targets), and
//! `p0 = P(present | chance)` by ion length, m/z and charge, learned from the decoy rows. Both are
//! learned per run on one half of the peptidoforms and applied to the other half, so no row
//! contributes to the tables that score it.
//!
//! The same-apex term is `−ln P(X ≥ m)`, `X ~ Binomial(k, π0)`, where `k` is the number of ions
//! present in the scan or the previous or next cycle of its isolation window and `m` the number
//! whose most intense of the three is the row's own scan; `π0` is the share measured on the decoy
//! rows.
//!
//! Explaining away, within one scan: `greedy` takes rows in descending LLR, and once a row is
//! taken the peaks it matched no longer count for the rows after it (their ions become
//! uninformative); a row's value is its LLR when taken. `soft` shares every matched peak's
//! intensity among the rows that match it, in proportion to each row's predicted relative
//! intensity times its estimated abundance (two rounds), then scores fractional presence with the
//! same tables and the spectral angle against the shared intensities.

use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

use model::scan::Ms1Scan;
use model::stats::MAD_TO_SD;
use model::Spectrum;
use onnx_infer::predictors::form::{Mass, H2O, PROTON};
use onnx_infer::predictors::spectral::FragmentTolerance;
use quant::Ms1RunIndex;
use rayon::prelude::*;

use crate::predicted::{anchor_rows, read_pin, write_pin, PredictedMs2, NONE};

/// PIN columns this module writes, in order.
pub(crate) const COLUMNS: [&str; 7] = [
    "dia_llr",
    "dia_apex_llr",
    "dia_llr_greedy",
    "dia_llr_soft",
    "dia_angle_soft",
    "dia_ion_frac",
    "dia_ppm_rms",
];

/// Ions whose ppm error calibrates the run's fragment m/z shift and spread.
const CALIB_RANK: u16 = 6;
const CALIB_MIN_IONS: usize = 500;
/// The RMS ppm error reported for a row with no top ion present.
const TOL_FLOOR_PPM: f64 = 5.0;

/// A run's fragment m/z calibration: the median ppm error (observed − theoretical) of the
/// confident rows' top ions and the tolerance around the shifted position (the search tolerance).
#[derive(Clone, Copy)]
struct Calib {
    shift_ppm: f64,
    tol: FragmentTolerance,
}

impl Calib {
    fn centre(&self, mz: f64) -> f64 {
        mz * (1.0 + self.shift_ppm * 1e-6)
    }
}

/// m/z tolerance (ppm) of the precursor's monoisotopic peak in the MS1 scan nearest the pool row,
/// whose intensity gives the row's MS1 quartile (a coarse abundance bin, not an evidence term).
const MS1_PPM: f64 = 10.0;
const RANK_BINS: usize = 12;
const MS1_BINS: usize = 5;
const LEN_BINS: usize = 4;
const MZ_BINS: usize = 20;
const Z_BINS: usize = 2;
const SOFT_ROUNDS: usize = 2;
const P_MIN: f64 = 1e-4;
/// Anchors of a (rank, MS1 quartile) cell below which its presence rate pools the rank's quartiles.
const Q_MIN_ANCHORS: f64 = 20.0;
/// Width (Th) of the fragment m/z bins of the chance table.
const MZ_BIN_TH: f64 = 100.0;

fn rank_bin(rank: usize, pred: f32) -> usize {
    if pred <= 0.0 {
        return 11;
    }
    match rank {
        1..=6 => rank - 1,
        7 | 8 => 6,
        9 | 10 => 7,
        11..=15 => 8,
        16..=20 => 9,
        _ => 10,
    }
}

fn len_bin(len: u8) -> usize {
    (len.max(1) as usize).min(LEN_BINS) - 1
}

fn mz_bin(mz: f64) -> usize {
    ((mz / MZ_BIN_TH) as usize).min(MZ_BINS - 1)
}

#[derive(Clone, Copy)]
struct IonObs {
    len: u8,
    z: u8,
    mz: f64,
    pred: f32,
    rank: u16,
    prev: f32,
    apex: f32,
    next: f32,
    peak: u32,
    ppm: f32,
}

impl IonObs {
    fn present(&self) -> bool {
        self.peak != NONE
    }
    fn any(&self) -> bool {
        self.prev > 0.0 || self.apex > 0.0 || self.next > 0.0
    }
    fn apex_here(&self) -> bool {
        self.apex > 0.0 && self.apex >= self.prev && self.apex >= self.next
    }
}

fn best_peak_ix(peaks: &[(f64, f32)], mz: f64, tol: FragmentTolerance) -> Option<usize> {
    let w = match tol {
        FragmentTolerance::Da(d) => d,
        FragmentTolerance::Ppm(p) => mz * p * 1e-6,
    };
    let lo = peaks.partition_point(|p| p.0 < mz - w);
    let mut best: Option<usize> = None;
    for (j, p) in peaks[lo..].iter().enumerate() {
        if p.0 > mz + w {
            break;
        }
        if best.is_none_or(|b| peaks[b].1 < p.1) {
            best = Some(lo + j);
        }
    }
    best
}

#[allow(clippy::too_many_arguments)]
/// The row's predicted b/y ions inside the scan's peak range, with their presence in the scan
/// (`here`) and the previous and next cycle of its window.
fn row_ions(
    mass: &Mass,
    charge: u8,
    pred: &[f32],
    layout: &[onnx_infer::predictors::spectral::Ion],
    here: &[(f64, f32)],
    prev: Option<&[(f64, f32)]>,
    next: Option<&[(f64, f32)]>,
    cal: Calib,
) -> Vec<IonObs> {
    let tol = cal.tol;
    let (Some(lo), Some(hi)) = (here.first().map(|p| p.0), here.last().map(|p| p.0)) else {
        return Vec::new();
    };
    let n = mass.res.len();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut prefix = mass.nterm;
    let mut out = Vec::new();
    let int_at = |s: Option<&[(f64, f32)]>, mz: f64| {
        s.and_then(|s| best_peak_ix(s, cal.centre(mz), tol).map(|j| s[j].1))
            .unwrap_or(0.0)
    };
    for i in 0..n.saturating_sub(1) {
        prefix += mass.res[i];
        for (k, &(is_b, z)) in layout.iter().enumerate() {
            if z > charge {
                continue;
            }
            let neutral = if is_b { prefix } else { total - prefix + H2O };
            let mz = (neutral + z as f64 * PROTON) / z as f64;
            if mz < lo || mz > hi {
                continue;
            }
            let peak = best_peak_ix(here, cal.centre(mz), tol);
            out.push(IonObs {
                len: if is_b { i + 1 } else { n - i - 1 } as u8,
                z,
                mz,
                pred: pred.get(i * layout.len() + k).copied().unwrap_or(0.0),
                rank: 0,
                prev: int_at(prev, mz),
                apex: peak.map_or(0.0, |j| here[j].1),
                next: int_at(next, mz),
                peak: peak.map_or(NONE, |j| j as u32),
                ppm: peak.map_or(0.0, |j| ((here[j].0 - mz) / mz * 1e6) as f32),
            });
        }
    }
    let mut order: Vec<usize> = (0..out.len()).collect();
    order.sort_by(|&a, &b| out[b].pred.total_cmp(&out[a].pred).then(a.cmp(&b)));
    for (r, &i) in order.iter().enumerate() {
        out[i].rank = (r + 1).min(u16::MAX as usize) as u16;
    }
    out
}

#[derive(Clone)]
struct Counts {
    q_n: [[f64; MS1_BINS]; RANK_BINS],
    q_hit: [[f64; MS1_BINS]; RANK_BINS],
    p_n: [[[f64; Z_BINS]; MZ_BINS]; LEN_BINS],
    p_hit: [[[f64; Z_BINS]; MZ_BINS]; LEN_BINS],
    apex_k: f64,
    apex_m: f64,
}

impl Default for Counts {
    fn default() -> Self {
        Self {
            q_n: [[0.0; MS1_BINS]; RANK_BINS],
            q_hit: [[0.0; MS1_BINS]; RANK_BINS],
            p_n: [[[0.0; Z_BINS]; MZ_BINS]; LEN_BINS],
            p_hit: [[[0.0; Z_BINS]; MZ_BINS]; LEN_BINS],
            apex_k: 0.0,
            apex_m: 0.0,
        }
    }
}

#[allow(clippy::needless_range_loop)]
impl Counts {
    fn add(&mut self, o: &Counts) {
        for r in 0..RANK_BINS {
            for b in 0..MS1_BINS {
                self.q_n[r][b] += o.q_n[r][b];
                self.q_hit[r][b] += o.q_hit[r][b];
            }
        }
        for l in 0..LEN_BINS {
            for m in 0..MZ_BINS {
                for z in 0..Z_BINS {
                    self.p_n[l][m][z] += o.p_n[l][m][z];
                    self.p_hit[l][m][z] += o.p_hit[l][m][z];
                }
            }
        }
        self.apex_k += o.apex_k;
        self.apex_m += o.apex_m;
    }
}

/// Log-likelihood terms of one (run, fold): `pos`/`neg` per (rank bin, MS1 bin, length bin,
/// m/z bin, charge bin) are built on the fly from `q` and `p0`.
struct Tables {
    q: [[f64; MS1_BINS]; RANK_BINS],
    p0: [[[f64; Z_BINS]; MZ_BINS]; LEN_BINS],
    pi0: f64,
}

#[allow(clippy::needless_range_loop)]
impl Tables {
    fn from_counts(c: &Counts) -> Self {
        let rate = |hit: f64, n: f64| ((hit + 1.0) / (n + 2.0)).clamp(P_MIN, 1.0 - P_MIN);
        let mut q = [[0.0; MS1_BINS]; RANK_BINS];
        for r in 0..RANK_BINS {
            for b in 0..MS1_BINS {
                let (mut h, mut n) = (c.q_hit[r][b], c.q_n[r][b]);
                if n < Q_MIN_ANCHORS {
                    h = c.q_hit[r].iter().sum();
                    n = c.q_n[r].iter().sum();
                }
                q[r][b] = rate(h, n);
            }
        }
        let mut p0 = [[[0.0; Z_BINS]; MZ_BINS]; LEN_BINS];
        for l in 0..LEN_BINS {
            for m in 0..MZ_BINS {
                for z in 0..Z_BINS {
                    p0[l][m][z] = rate(c.p_hit[l][m][z], c.p_n[l][m][z]);
                }
            }
        }
        let pi0 = rate(c.apex_m, c.apex_k);
        Self { q, p0, pi0 }
    }

    fn terms(&self, ion: &IonObs, ms1: usize) -> (f64, f64) {
        let q = self.q[rank_bin(ion.rank as usize, ion.pred)][ms1];
        let p = self.p0[len_bin(ion.len)][mz_bin(ion.mz)][(ion.z.max(1) as usize - 1).min(1)];
        ((q / p).ln(), ((1.0 - q) / (1.0 - p)).ln())
    }
}

/// `−ln P(X ≥ m)` for `X ~ Binomial(k, p)`.
fn binom_tail_nlog(k: usize, m: usize, p: f64) -> f64 {
    if m == 0 {
        return 0.0;
    }
    let lp = p.ln();
    let lq = (1.0 - p).ln();
    let mut log_c = 0.0f64;
    let mut terms = Vec::with_capacity(k + 1);
    for i in 0..=k {
        if i > 0 {
            log_c += ((k - i + 1) as f64).ln() - (i as f64).ln();
        }
        if i >= m {
            terms.push(log_c + i as f64 * lp + (k - i) as f64 * lq);
        }
    }
    let mx = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !mx.is_finite() {
        return 0.0;
    }
    let s: f64 = terms.iter().map(|t| (t - mx).exp()).sum();
    -(mx + s.ln()).min(0.0)
}

struct RowOut {
    row: u32,
    vals: [f64; COLUMNS.len()],
}

/// Adds the `COLUMNS` to every row of the PIN (the full pool, before re-ranking).
#[allow(clippy::too_many_arguments)]
/// The pool table with its `dia_llr` per row and each run's learned fragment m/z shift.
pub(crate) struct PoolScores {
    pub(crate) table: crate::predicted::PinTable,
    pub(crate) llr: Vec<f64>,
    pub(crate) shift_ppm: Vec<f64>,
}

pub(crate) fn add_to_pin(
    pin: &Path,
    spectra: &[Spectrum],
    run_ends: &[usize],
    ms1: &[Vec<Ms1Scan>],
    predicted: &PredictedMs2,
    tol0: FragmentTolerance,
) -> Result<PoolScores, String> {
    let t0 = Instant::now();
    let cols: Vec<String> = COLUMNS.iter().map(|c| c.to_string()).collect();
    let table = read_pin(pin, spectra, &cols)?;
    let n = table.rows.len();
    let run_of = |i: usize| run_ends.partition_point(|&e| e <= i);
    let n_runs = run_ends.len().max(1);
    let (slot, windows) = crate::dia::windows(spectra, run_ends);
    let ms1: Vec<Option<Ms1RunIndex>> = ms1
        .iter()
        .map(|s| (!s.is_empty()).then(|| Ms1RunIndex::new(s.clone())))
        .collect();

    let mut by_spec: Vec<Vec<u32>> = vec![Vec::new(); spectra.len()];
    for (i, r) in table.rows.iter().enumerate() {
        if r.spec != NONE && r.form != NONE && r.charge > 0 {
            by_spec[r.spec as usize].push(i as u32);
        }
    }
    let mut anchor = vec![false; n];
    for i in anchor_rows(&table.rows, n_runs, &|s| run_of(s as usize))
        .into_iter()
        .flatten()
    {
        anchor[i] = true;
    }
    let fold = |form: u32| ((form as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 63) as usize;

    let ms1_int: Vec<f32> = table
        .rows
        .par_iter()
        .map(|r| {
            if r.spec == NONE || r.form == NONE || r.charge == 0 {
                return 0.0;
            }
            let Some(Some(idx)) = ms1.get(run_of(r.spec as usize)) else {
                return 0.0;
            };
            let mass = &table.forms[r.form as usize].1;
            let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
            let mono = (total + H2O) / r.charge as f64 + PROTON;
            let rt = spectra[r.spec as usize].rt_seconds.unwrap_or(0.0);
            idx.nearest_scan(rt)
                .and_then(|i| {
                    let peaks = &idx.scan(i).peaks;
                    let tol = mono * MS1_PPM * 1e-6;
                    let lo = peaks.partition_point(|p| p.0 < mono - tol);
                    peaks[lo..]
                        .iter()
                        .take_while(|p| p.0 <= mono + tol)
                        .map(|p| p.1)
                        .reduce(f32::max)
                })
                .unwrap_or(0.0)
        })
        .collect();
    let mut cuts: Vec<Vec<f32>> = vec![Vec::new(); n_runs];
    for (i, r) in table.rows.iter().enumerate() {
        if anchor[i] && ms1_int[i] > 0.0 {
            cuts[run_of(r.spec as usize).min(n_runs - 1)].push(ms1_int[i]);
        }
    }
    let cuts: Vec<[f32; 3]> = cuts
        .into_iter()
        .map(|mut v| {
            v.sort_by(f32::total_cmp);
            let at = |f: f64| v.get((v.len() as f64 * f) as usize).copied().unwrap_or(0.0);
            [at(0.25), at(0.5), at(0.75)]
        })
        .collect();
    let ms1_bin = |i: usize, run: usize| -> usize {
        let v = ms1_int[i];
        if v <= 0.0 {
            0
        } else {
            1 + cuts[run.min(n_runs - 1)].iter().filter(|&&c| v > c).count()
        }
    };

    let layout = &predicted.ions;
    let spec_rows = |s: usize, cals: &[Calib]| -> Vec<(u32, Vec<IonObs>)> {
        let cal = cals[run_of(s).min(n_runs - 1)];
        let (w, p) = slot[s];
        let win = &windows[w as usize];
        let p = p as usize;
        let nb = |k: Option<usize>| k.map(|k| spectra[win[k] as usize].peaks.as_slice());
        let prev = nb(p.checked_sub(1));
        let next = nb((p + 1 < win.len()).then_some(p + 1));
        let here = spectra[s].peaks.as_slice();
        by_spec[s]
            .iter()
            .map(|&i| {
                let r = &table.rows[i as usize];
                let (form, mass) = &table.forms[r.form as usize];
                let ions = predicted
                    .by_form
                    .get(&(form.clone(), r.charge))
                    .map(|pred| row_ions(mass, r.charge, pred, layout, here, prev, next, cal))
                    .unwrap_or_default();
                (i, ions)
            })
            .collect()
    };

    let first = vec![
        Calib {
            shift_ppm: 0.0,
            tol: tol0
        };
        n_runs
    ];
    let errs: Vec<Vec<f32>> = (0..spectra.len())
        .into_par_iter()
        .filter(|&s| !by_spec[s].is_empty())
        .fold(
            || vec![Vec::new(); n_runs],
            |mut acc, s| {
                let run = run_of(s).min(n_runs - 1);
                for (i, ions) in spec_rows(s, &first) {
                    if anchor[i as usize] {
                        acc[run].extend(
                            ions.iter()
                                .filter(|o| o.present() && o.rank <= CALIB_RANK && o.pred > 0.0)
                                .map(|o| o.ppm),
                        );
                    }
                }
                acc
            },
        )
        .reduce(
            || vec![Vec::new(); n_runs],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(b) {
                    x.extend(y);
                }
                a
            },
        );
    let cals: Vec<Calib> = errs
        .into_iter()
        .enumerate()
        .map(|(run, mut e)| {
            if e.len() < CALIB_MIN_IONS {
                return first[run];
            }
            e.sort_by(f32::total_cmp);
            let med = f64::from(e[e.len() / 2]);
            let mut dev: Vec<f32> = e.iter().map(|&x| (f64::from(x) - med).abs() as f32).collect();
            dev.sort_by(f32::total_cmp);
            let sigma = MAD_TO_SD * f64::from(dev[dev.len() / 2]);
            log::info!(
                target: "dia",
                "run {run}: fragment shift {med:.2} ppm, spread {sigma:.2} ppm from {} ions",
                e.len()
            );
            Calib { shift_ppm: med, tol: tol0 }
        })
        .collect();

    let counts: Vec<Counts> = (0..spectra.len())
        .into_par_iter()
        .filter(|&s| !by_spec[s].is_empty())
        .fold(
            || vec![Counts::default(); n_runs * 2],
            |mut acc, s| {
                let run = run_of(s).min(n_runs - 1);
                for (i, ions) in spec_rows(s, &cals) {
                    let i = i as usize;
                    let r = &table.rows[i];
                    let c = &mut acc[run * 2 + fold(r.form)];
                    if anchor[i] {
                        let b = ms1_bin(i, run);
                        for o in &ions {
                            let rb = rank_bin(o.rank as usize, o.pred);
                            c.q_n[rb][b] += 1.0;
                            c.q_hit[rb][b] += f64::from(u8::from(o.present()));
                        }
                    }
                    if r.decoy {
                        for o in &ions {
                            let (l, m, z) =
                                (len_bin(o.len), mz_bin(o.mz), (o.z as usize - 1).min(1));
                            c.p_n[l][m][z] += 1.0;
                            c.p_hit[l][m][z] += f64::from(u8::from(o.present()));
                            if o.any() {
                                c.apex_k += 1.0;
                                c.apex_m += f64::from(u8::from(o.apex_here()));
                            }
                        }
                    }
                }
                acc
            },
        )
        .reduce(
            || vec![Counts::default(); n_runs * 2],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(&b) {
                    x.add(y);
                }
                a
            },
        );
    let tables: Vec<Tables> = (0..n_runs * 2)
        .map(|k| Tables::from_counts(&counts[k ^ 1]))
        .collect();

    let per_spec: Vec<Vec<RowOut>> = (0..spectra.len())
        .into_par_iter()
        .filter(|&s| !by_spec[s].is_empty())
        .map(|s| {
            let run = run_of(s).min(n_runs - 1);
            let rows = spec_rows(s, &cals);
            let mut outs = Vec::with_capacity(rows.len());
            let mut terms: Vec<Vec<(f64, f64)>> = Vec::with_capacity(rows.len());
            for (i, ions) in &rows {
                let i = *i as usize;
                let r = &table.rows[i];
                let t = &tables[run * 2 + fold(r.form)];
                let b = ms1_bin(i, run);
                let tv: Vec<(f64, f64)> = ions.iter().map(|o| t.terms(o, b)).collect();
                let llr: f64 = ions
                    .iter()
                    .zip(&tv)
                    .map(|(o, (pos, neg))| if o.present() { *pos } else { *neg })
                    .sum();
                let top: Vec<&IonObs> = ions
                    .iter()
                    .filter(|o| o.rank <= CALIB_RANK && o.pred > 0.0)
                    .collect();
                let frac = if top.is_empty() {
                    0.0
                } else {
                    top.iter().filter(|o| o.any()).count() as f64 / top.len() as f64
                };
                let hit: Vec<f64> = top
                    .iter()
                    .filter(|o| o.present())
                    .map(|o| f64::from(o.ppm) - cals[run].shift_ppm)
                    .collect();
                let rms = if hit.is_empty() {
                    TOL_FLOOR_PPM
                } else {
                    (hit.iter().map(|x| x * x).sum::<f64>() / hit.len() as f64).sqrt()
                };
                let k = ions.iter().filter(|o| o.any()).count();
                let m = ions.iter().filter(|o| o.any() && o.apex_here()).count();
                let apex = binom_tail_nlog(k, m, t.pi0);
                outs.push(RowOut {
                    row: i as u32,
                    vals: [llr, apex, llr, llr, 0.0, frac, rms],
                });
                terms.push(tv);
            }
            explain_greedy(&rows, &terms, &mut outs);
            explain_soft(&rows, &terms, &spectra[s].peaks, &mut outs);
            outs
        })
        .collect();

    let mut values = vec![vec![0.0; n]; COLUMNS.len()];
    for outs in &per_spec {
        for o in outs {
            for (c, v) in values.iter_mut().zip(o.vals) {
                c[o.row as usize] = v;
            }
        }
    }
    crate::memlimit::mark_stage("pool scores: columns computed");
    write_pin(pin, &table, &cols, &values)?;
    let shift_ppm: Vec<f64> = cals.iter().map(|c| c.shift_ppm).collect();
    let llr = values.swap_remove(0);
    let report: Vec<String> = tables
        .iter()
        .enumerate()
        .map(|(k, t)| {
            format!(
                "run {} fold {}: q(rank 1, mid MS1) {:.3}, q(rank 21+) {:.3}, pi0 {:.3}",
                k / 2,
                k % 2,
                t.q[0][2],
                t.q[10][2],
                t.pi0
            )
        })
        .collect();
    log::info!(
        target: "dia",
        "pool scores: {} columns on {n} PIN rows; {} [{:.2}s]",
        COLUMNS.len(),
        report.join("; "),
        t0.elapsed().as_secs_f64()
    );
    Ok(PoolScores {
        table,
        llr,
        shift_ppm,
    })
}

fn explain_greedy(rows: &[(u32, Vec<IonObs>)], terms: &[Vec<(f64, f64)>], outs: &mut [RowOut]) {
    let mut cur: Vec<f64> = outs.iter().map(|o| o.vals[0]).collect();
    let mut taken = vec![false; rows.len()];
    let mut explained: std::collections::HashSet<u32> = std::collections::HashSet::new();
    for _ in 0..rows.len() {
        let Some(best) = (0..rows.len())
            .filter(|&j| !taken[j])
            .max_by(|&a, &b| cur[a].total_cmp(&cur[b]).then(b.cmp(&a)))
        else {
            break;
        };
        taken[best] = true;
        outs[best].vals[2] = cur[best];
        let new: Vec<u32> = rows[best]
            .1
            .iter()
            .filter(|o| o.present() && !explained.contains(&o.peak))
            .map(|o| o.peak)
            .collect();
        if new.is_empty() {
            continue;
        }
        explained.extend(new.iter().copied());
        for j in (0..rows.len()).filter(|&j| !taken[j]) {
            for (o, (pos, _)) in rows[j].1.iter().zip(&terms[j]) {
                if o.present() && new.contains(&o.peak) {
                    cur[j] -= pos;
                }
            }
        }
    }
}

fn explain_soft(
    rows: &[(u32, Vec<IonObs>)],
    terms: &[Vec<(f64, f64)>],
    peaks: &[(f64, f32)],
    outs: &mut [RowOut],
) {
    let norm: Vec<f64> = rows
        .iter()
        .map(|(_, ions)| ions.iter().map(|o| o.pred.max(0.0) as f64).sum::<f64>())
        .collect();
    let rel = |j: usize, o: &IonObs| {
        if norm[j] > 0.0 {
            o.pred.max(0.0) as f64 / norm[j]
        } else {
            0.0
        }
    };
    let mut abundance: Vec<f64> = rows
        .iter()
        .enumerate()
        .map(|(j, (_, ions))| {
            let (mut obs, mut pred) = (0.0, 0.0);
            for o in ions.iter().filter(|o| o.present()) {
                obs += peaks[o.peak as usize].1 as f64;
                pred += rel(j, o);
            }
            if pred > 0.0 {
                obs / pred
            } else {
                0.0
            }
        })
        .collect();
    let mut frac: Vec<Vec<f64>> = rows.iter().map(|(_, v)| vec![0.0; v.len()]).collect();
    for _ in 0..SOFT_ROUNDS {
        let mut demand: HashMap<u32, f64> = HashMap::new();
        for (j, (_, ions)) in rows.iter().enumerate() {
            for o in ions.iter().filter(|o| o.present()) {
                *demand.entry(o.peak).or_default() += abundance[j] * rel(j, o);
            }
        }
        for (j, (_, ions)) in rows.iter().enumerate() {
            let (mut obs, mut pred) = (0.0, 0.0);
            for (x, o) in frac[j].iter_mut().zip(ions) {
                if !o.present() {
                    continue;
                }
                let d = demand[&o.peak];
                let w = abundance[j] * rel(j, o);
                *x = if d > 0.0 { w / d } else { 0.0 };
                obs += *x * peaks[o.peak as usize].1 as f64;
                pred += rel(j, o);
            }
            abundance[j] = if pred > 0.0 { obs / pred } else { 0.0 };
        }
    }
    for (j, (_, ions)) in rows.iter().enumerate() {
        let llr: f64 = frac[j]
            .iter()
            .zip(&terms[j])
            .map(|(x, (pos, neg))| x * pos + (1.0 - x) * neg)
            .sum();
        let (mut dot, mut pp, mut oo) = (0.0, 0.0, 0.0);
        for (x, o) in frac[j].iter().zip(ions) {
            let p = rel(j, o);
            let v = if o.present() {
                x * peaks[o.peak as usize].1 as f64
            } else {
                0.0
            };
            dot += p * v;
            pp += p * p;
            oo += v * v;
        }
        let angle = if pp > 0.0 && oo > 0.0 {
            1.0 - 2.0 * (dot / (pp * oo).sqrt()).clamp(-1.0, 1.0).acos() / std::f64::consts::PI
        } else {
            0.0
        };
        outs[j].vals[3] = llr;
        outs[j].vals[4] = angle;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use onnx_infer::predictors::form::parse_peptide;

    #[test]
    fn binomial_tail_is_zero_for_no_hits_and_grows_with_agreement() {
        assert_eq!(binom_tail_nlog(10, 0, 0.5), 0.0);
        let a = binom_tail_nlog(10, 5, 0.5);
        let b = binom_tail_nlog(10, 10, 0.5);
        assert!(a < b);
        assert!((b - 10.0 * 2f64.ln()).abs() < 1e-9);
        assert!((binom_tail_nlog(1, 1, 0.25) - 4f64.ln()).abs() < 1e-9);
    }

    #[test]
    fn rows_find_their_ions_and_rank_them_by_prediction() {
        let (_, mass) = parse_peptide("K.PEPTIDEK.R").unwrap();
        let layout = [(true, 1u8), (false, 1u8)];
        let n = mass.res.len();
        let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
        let mut prefix = mass.nterm;
        let mut here = Vec::new();
        let mut pred = Vec::new();
        for i in 0..n - 1 {
            prefix += mass.res[i];
            here.push((total - prefix + H2O + PROTON, 100.0 + i as f32));
            pred.extend([0.1f32, 1.0 + i as f32]);
        }
        here.sort_by(|a, b| a.0.total_cmp(&b.0));
        let ions = row_ions(
            &mass,
            2,
            &pred,
            &layout,
            &here,
            None,
            Some(&here),
            Calib {
                shift_ppm: 0.0,
                tol: FragmentTolerance::Ppm(7.0),
            },
        );
        let ys: Vec<&IonObs> = ions.iter().filter(|o| o.pred >= 1.0).collect();
        assert_eq!(ys.len(), n - 1);
        assert!(ys.iter().all(|o| o.present() && o.apex_here()));
        let top = ions.iter().find(|o| o.rank == 1).unwrap();
        assert!(top.pred >= 1.0 && top.len == 1);
    }

    fn obs(peak: u32, pred: f32) -> IonObs {
        IonObs {
            len: 5,
            z: 1,
            mz: 500.0,
            pred,
            rank: 1,
            prev: 0.0,
            apex: if peak == NONE { 0.0 } else { 1.0 },
            next: 0.0,
            peak,
            ppm: 0.0,
        }
    }

    fn out(row: u32, llr: f64) -> RowOut {
        RowOut {
            row,
            vals: [llr, 0.0, llr, llr, 0.0, 0.0, 0.0],
        }
    }

    #[test]
    fn greedy_takes_shared_peaks_from_the_weaker_row() {
        let rows = vec![
            (0u32, vec![obs(0, 1.0), obs(1, 1.0)]),
            (1u32, vec![obs(0, 1.0), obs(NONE, 1.0)]),
        ];
        let terms = vec![vec![(2.0, -1.0); 2], vec![(2.0, -1.0); 2]];
        let mut outs = vec![out(0, 4.0), out(1, 1.0)];
        explain_greedy(&rows, &terms, &mut outs);
        assert_eq!(outs[0].vals[2], 4.0);
        assert_eq!(outs[1].vals[2], -1.0);
    }

    #[test]
    fn soft_sharing_splits_a_peak_by_predicted_abundance() {
        let peaks = vec![(500.0, 100.0f32), (600.0, 100.0)];
        let rows = vec![
            (0u32, vec![obs(0, 1.0), obs(1, 1.0)]),
            (1u32, vec![obs(0, 1.0), obs(NONE, 1.0)]),
        ];
        let terms = vec![vec![(2.0, -1.0); 2], vec![(2.0, -1.0); 2]];
        let mut outs = vec![out(0, 4.0), out(1, 1.0)];
        explain_soft(&rows, &terms, &peaks, &mut outs);
        assert!(outs[0].vals[3] > outs[1].vals[3]);
        assert!(outs[0].vals[3] <= 4.0 && outs[1].vals[3] < 1.0);
        assert!(outs[0].vals[4] > 0.8);
    }
}
