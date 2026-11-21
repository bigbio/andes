//! DIA pass 2: precursor-centric evidence for the nominated precursors of each run.
//!
//! A nominee is the best pool row of a (peptidoform, charge) in a run by `dia_llr`, targets and
//! decoys alike. Each nominee is traced over the cycles of its isolation window around its
//! calibrated predicted RT; the presence of its top predicted fragments per cycle is scored
//! against chance with per-run rates, and the best cycle gives the pass-2 columns.
//!
//! Each row also gets its quantities in the run, over the elution window: the apex ± the run's
//! learned peak width, the median half-maximum width (cycles) of the summed fragment traces of
//! the run's confident targets. The apex is the local maximum of the row's summed fragment
//! trace reached from its best cycle. The fragment quantity is the summed intensity over that
//! window of the [`QUANT_TOP_FRAGMENTS`] most intense predicted fragments, each above its
//! background: the median of its trace over the cycles beyond one more learned width on each
//! side of the window, out to [`QUANT_BACKGROUND_REACH_HALF_WIDTHS`] learned widths from the
//! apex. The same fragments are summed in every run, whatever their per-run interference, so a
//! ratio between runs compares like with like; the MS1 quantity is the summed intensity of the
//! monoisotopic and M+1 MS1 traces over the same retention-time span. A row whose window does
//! not lie inside its traced cycles (the peak is cut by the run's first or last spectrum or by
//! the trace) or that has no signal has no quantity; a quantity is never 0.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use model::mass::ISOTOPE;
use model::stats::MAD_TO_SD;
use model::Spectrum;
use onnx_infer::predictors::form::{H2O, PROTON};
use quant::Ms1RunIndex;
use rayon::prelude::*;

use crate::predicted::{PinTable, PredictedMs2, NONE};

pub(crate) const TOP_IONS: usize = 12;
const NOMINATE_MIN_LLR: f64 = -15.0;
const CONFIDENT_LLR: f64 = 15.0;
/// `dia_llr` of a nominee transferred from another run of a study (no pool row in this run).
pub(crate) const TRANSFER_DIA_LLR: f64 = -50.0;
/// Name of the pass-2 column that marks transferred nominees (1) against the run's own (0).
pub(crate) const XFER_COLUMN: &str = "xfer";
pub(crate) const HALF_CYCLES: usize = 36;
pub(crate) const CYCLES: usize = 2 * HALF_CYCLES + 1;
/// m/z tolerance (ppm) of a predicted fragment's peak in a cycle of its trace, around the run's
/// shifted fragment position; tighter than the pool search's tolerance because the shift is
/// learned first.
const FRAG_TRACE_PPM: f64 = 7.0;
/// m/z tolerance (ppm) of the precursor's isotope peaks in the MS1 scans around the fragment apex.
const MS1_TRACE_PPM: f64 = 7.0;
const MS1_HALF: usize = 5;
/// Isotope peaks (monoisotopic first) looked up for `p2_ms1_iso`.
const MS1_ISOTOPES: u8 = 4;
/// Half-width (Th) of a window's isolation range when the spectrum gives no isolation offsets.
const WINDOW_REACH_MZ: f64 = 1.6;
const RIDGE_LAMBDA: f64 = 1e-3;
/// Pivot magnitude below which the RT ridge system counts as singular.
const PIVOT_MIN: f64 = 1e-12;
/// Ions, most intense first, that `p2_n6` and `p2_any6` count.
const TOP_SUBSET: usize = 6;
/// Floor of a per-cycle log-likelihood ratio where a cycle has no value: the second-best cycle of
/// `p2_gap` and `p2_llr2`, and a neighbour of the best cycle outside the trace or the run in
/// `p2_lsum3`.
const LLR_FLOOR: f64 = -50.0;
/// Bounds of the per-run presence rates at the apex (`q`) and by chance (`p0`).
const Q_MIN: f64 = 1e-3;
const P0_MIN: f64 = 1e-4;
const RATE_MAX: f64 = 0.999;
/// Added to the spread of the per-cycle log-likelihood ratios in `p2_z`.
const Z_SD_EPS: f64 = 1e-3;
/// Added to denominators of correlations and fractions.
const DIV_EPS: f64 = 1e-9;
const COELUTE_HALF: usize = 3;
const COELUTE_MIN_CORR: f64 = 0.6;
const FAR_CYCLES: usize = 3;
const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// RT prior spread when a run has too few confident targets to measure it.
const RT_SIGMA_FALLBACK_MIN: f64 = 0.5;
const RT_TERMS: usize = 4 + 20;
/// Fragments summed into the fragment quantity: the most intense predicted.
const QUANT_TOP_FRAGMENTS: usize = 6;
/// Reach (learned half-widths from the apex) of the cycles that give a fragment's background.
const QUANT_BACKGROUND_REACH_HALF_WIDTHS: usize = 6;
/// Background cycles a fragment's background needs; with fewer it is 0.
const QUANT_BACKGROUND_MIN_CYCLES: usize = 4;
/// Smallest half-width (cycles) of the quantity window.
const QUANT_MIN_HALF_CYCLES: usize = 2;
/// Half-width (cycles) of the quantity window when a run has too few confident targets to
/// learn its peak width.
const QUANT_HALF_CYCLES_FALLBACK: usize = 3;
/// Confident targets with a measurable peak a run needs to learn its peak width.
const QUANT_WIDTH_MIN_TARGETS: usize = 50;

pub(crate) const COLUMNS: [&str; 27] = [
    "p2_llr",
    "p2_llr2",
    "p2_z",
    "p2_gap",
    "p2_n6",
    "p2_nall",
    "p2_any6",
    "p2_expl",
    "p2_angle",
    "p2_cor",
    "p2_cor6",
    "p2_apex",
    "p2_lsum3",
    "p2_logI",
    "p2_width",
    "p2_rterr",
    "p2_rtstar",
    "p2_ms1_log",
    "p2_ms1_iso",
    "p2_ms1_off",
    "p2_ms1_corr",
    "p2_shared",
    "p2_lone",
    "p2_lone_gap",
    "p2_dia_llr",
    "p2_rt_pred_err",
    "p2_ok",
];

/// One nominated precursor of one run.
#[derive(Clone, Debug)]
pub(crate) struct Nominee {
    pub(crate) run: u32,
    pub(crate) charge: u8,
    pub(crate) decoy: bool,
    pub(crate) prec_mz: f64,
    /// Predicted RT (min) of the peptidoform.
    pub(crate) pred_rt: f64,
    /// The PIN row had no andes-rt prediction (`andesrt_missing`): `pred_rt` is an imputed value,
    /// so the row stays out of the RT ridge and its spread.
    pub(crate) rt_missing: bool,
    /// The raw (uncalibrated) andes-rt prediction of the peptidoform; NaN without one. Stored
    /// for study mode, where a transferred nominee is calibrated with the receiving run's map.
    pub(crate) raw_rt: f64,
    /// RT (min) of the scan whose pool row nominated it.
    pub(crate) apex_rt: f64,
    pub(crate) dia_llr: f64,
    /// Centre (scan index) of the nominating scan's mobility window; NaN without one.
    pub(crate) mobility: f64,
    /// Top predicted fragments, (m/z, predicted intensity), most intense first.
    pub(crate) ions: Vec<(f64, f32)>,
    /// Residue letters without modifications.
    pub(crate) seq: Box<str>,
    pub(crate) peptide: Box<str>,
    pub(crate) proteins: Box<str>,
    pub(crate) spec_id: Box<str>,
    /// Nominated in another run of the study and scored here: never in a fit set.
    pub(crate) xfer: bool,
}

/// A confident target of the run: the RT ridge, the presence rates at the apex and the joint
/// model's fits use these. Transferred nominees never count.
pub(crate) fn confident(n: &Nominee) -> bool {
    !n.xfer && !n.decoy && n.dia_llr > CONFIDENT_LLR
}

/// A decoy of the run's chance model (presence rates, joint-model decoy sample); transferred
/// decoys never count.
pub(crate) fn null_member(n: &Nominee) -> bool {
    n.decoy && !n.xfer
}

/// Pass-2 columns of the nominees, in `COLUMNS` order; rows outside every isolation window
/// carry `p2_ok` = 0 and are left out downstream.
pub(crate) struct Pass2 {
    pub(crate) cols: Vec<Vec<f32>>,
    pub(crate) joint: Vec<Vec<f32>>,
    pub(crate) quant: Vec<RowQuant>,
}

/// Quantities of one row in its run (module doc).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct RowQuant {
    pub(crate) frag: Option<f64>,
    pub(crate) n_frag: u8,
    pub(crate) ms1: Option<f64>,
}

impl Pass2 {
    pub(crate) fn col(&self, name: &str) -> &[f32] {
        let i = COLUMNS.iter().position(|c| *c == name).unwrap_or(0);
        &self.cols[i]
    }
}

/// The best pool row of every (run, peptidoform, charge) with `llr` above the nomination
/// threshold. `llr` is the pool `dia_llr` of each PIN row of `table`.
pub(crate) fn nominate(
    pin: &Path,
    spectra: &[Spectrum],
    run_ends: &[usize],
    table: &PinTable,
    llr: &[f64],
    predicted: &PredictedMs2,
) -> Result<Vec<Nominee>, String> {
    let run_of = |s: usize| run_ends.partition_point(|&e| e <= s);
    if !predicted.rt_raw.is_empty() && predicted.rt_raw.len() != table.forms.len() {
        return Err(format!(
            "{}: {} peptidoforms with an andes-rt prediction, {} in the pool table",
            pin.display(),
            predicted.rt_raw.len(),
            table.forms.len()
        ));
    }
    let mut best: HashMap<(u32, u8, usize), usize> = HashMap::new();
    for (i, r) in table.rows.iter().enumerate() {
        if r.spec == NONE || r.form == NONE || r.charge == 0 || llr[i] <= NOMINATE_MIN_LLR {
            continue;
        }
        let e = best
            .entry((r.form, r.charge, run_of(r.spec as usize)))
            .or_insert(i);
        if llr[i] > llr[*e] {
            *e = i;
        }
    }
    let mut want: Vec<usize> = best.into_values().collect();
    want.sort_unstable();
    let fields = pin_fields(pin, &want)?;
    let layout = &predicted.ions;
    let out = want
        .par_iter()
        .zip(fields.into_par_iter())
        .map(|(&i, (pred_rt, rt_missing, peptide, proteins, spec_id))| {
            let r = &table.rows[i];
            let (form, mass) = &table.forms[r.form as usize];
            let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
            let mut ions: Vec<(f64, f32)> = Vec::new();
            if let Some(pred) = predicted.get(r.form, r.charge) {
                let mut prefix = mass.nterm;
                for site in 0..mass.res.len().saturating_sub(1) {
                    prefix += mass.res[site];
                    for (k, &(is_b, z)) in layout.iter().enumerate() {
                        if z > r.charge {
                            continue;
                        }
                        let neutral = if is_b { prefix } else { total - prefix + H2O };
                        let mz = (neutral + f64::from(z) * PROTON) / f64::from(z);
                        let p = pred.get(site * layout.len() + k).copied().unwrap_or(0.0);
                        ions.push((mz, p.max(0.0)));
                    }
                }
                ions.sort_by(|a, b| b.1.total_cmp(&a.1));
                ions.truncate(TOP_IONS);
            }
            let sp = &spectra[r.spec as usize];
            Nominee {
                run: run_of(r.spec as usize) as u32,
                charge: r.charge,
                decoy: r.decoy,
                prec_mz: (total + H2O) / f64::from(r.charge) + PROTON,
                pred_rt,
                rt_missing,
                raw_rt: predicted
                    .rt_raw
                    .get(r.form as usize)
                    .map_or(f64::NAN, |&v| f64::from(v)),
                apex_rt: sp.rt_seconds.unwrap_or(0.0) / 60.0,
                dia_llr: llr[i],
                mobility: mobility_centre(sp.mobility_window),
                ions,
                seq: form.seq.as_str().into(),
                peptide: peptide.into(),
                proteins: proteins.into(),
                spec_id: spec_id.into(),
                xfer: false,
            }
        })
        .filter(|n| !n.ions.is_empty())
        .collect();
    Ok(out)
}

/// Predicted RT, its missing flag, Peptide, Proteins and SpecId of a PIN row.
type PinFields = (f64, bool, String, String, String);

/// Predicted RT, its missing flag, Peptide, Proteins and SpecId of the wanted PIN rows (sorted
/// row indices).
fn pin_fields(
    pin: &Path,
    want: &[usize],
) -> Result<Vec<PinFields>, String> {
    let f = std::fs::File::open(pin).map_err(|e| format!("{}: {e}", pin.display()))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("{}: empty PIN", pin.display()))?
        .map_err(|e| format!("{}: {e}", pin.display()))?;
    let h: Vec<&str> = header.split('\t').collect();
    let at = |n: &str| {
        h.iter()
            .position(|c| *c == n)
            .ok_or_else(|| format!("{}: PIN has no {n} column", pin.display()))
    };
    let (c_id, c_rt, c_miss, c_pep, c_prot) = (
        at("SpecId")?,
        at("andesrt_predicted_rt")?,
        at("andesrt_missing")?,
        at("Peptide")?,
        at("Proteins")?,
    );
    let mut out = Vec::with_capacity(want.len());
    let mut w = 0usize;
    let mut idx = 0usize;
    for line in lines {
        if w == want.len() {
            break;
        }
        let line = line.map_err(|e| format!("{}: {e}", pin.display()))?;
        if line.is_empty() {
            continue;
        }
        if want[w] == idx {
            let x: Vec<&str> = line.split('\t').collect();
            if x.len() <= c_prot.max(c_pep).max(c_rt).max(c_miss) {
                return Err(format!("{}: PIN row {} is short", pin.display(), idx + 2));
            }
            out.push((
                x[c_rt].trim().parse().unwrap_or(f64::NAN),
                x[c_miss].trim().parse::<f64>() != Ok(0.0),
                x[c_pep].to_string(),
                x[c_prot..].join(";"),
                x[c_id].to_string(),
            ));
            w += 1;
        }
        idx += 1;
    }
    if out.len() != want.len() {
        return Err(format!(
            "{}: {} of {} nominee rows found",
            pin.display(),
            out.len(),
            want.len()
        ));
    }
    Ok(out)
}

/// Ridge fit of the observed apex RT on predicted RT (cubic) and residue counts over the
/// confident targets of one run; returns the calibrated RT of every nominee given and the
/// robust spread (min) of the fit residuals.
fn calibrate_rt(noms: &[&Nominee]) -> (Vec<f64>, f64) {
    let design = |n: &Nominee| -> [f64; RT_TERMS] {
        let mut x = [0.0; RT_TERMS];
        let p = n.pred_rt;
        x[0] = 1.0;
        x[1] = p;
        x[2] = p * p;
        x[3] = p * p * p;
        for b in n.seq.bytes() {
            if let Some(k) = AA.iter().position(|&a| a == b) {
                x[4 + k] += 1.0;
            }
        }
        x
    };
    let mut a = [[0.0f64; RT_TERMS]; RT_TERMS];
    let mut v = [0.0f64; RT_TERMS];
    let mut n_fit = 0usize;
    for n in noms
        .iter()
        .filter(|n| rt_fit_member(n))
    {
        let x = design(n);
        for i in 0..RT_TERMS {
            v[i] += x[i] * n.apex_rt;
            for j in 0..RT_TERMS {
                a[i][j] += x[i] * x[j];
            }
        }
        n_fit += 1;
    }
    let beta = if n_fit > RT_TERMS {
        for (i, row) in a.iter_mut().enumerate() {
            row[i] += RIDGE_LAMBDA;
        }
        solve(a, v)
    } else {
        None
    };
    let cal: Vec<f64> = noms
        .iter()
        .map(|n| match beta {
            Some(b) => design(n).iter().zip(b.iter()).map(|(x, c)| x * c).sum(),
            None => n.pred_rt,
        })
        .collect();
    let mut res: Vec<f64> = noms
        .iter()
        .zip(&cal)
        .filter(|(n, _)| rt_fit_member(n))
        .map(|(n, c)| n.apex_rt - c)
        .collect();
    let sigma = robust_sigma(&mut res).unwrap_or(RT_SIGMA_FALLBACK_MIN);
    (cal, sigma)
}

/// `MAD_TO_SD` x median absolute deviation; `None` for fewer than 50 values.
fn robust_sigma(v: &mut [f64]) -> Option<f64> {
    if v.len() < 50 {
        return None;
    }
    v.sort_by(f64::total_cmp);
    let med = v[v.len() / 2];
    let mut dev: Vec<f64> = v.iter().map(|x| (x - med).abs()).collect();
    dev.sort_by(f64::total_cmp);
    Some(MAD_TO_SD * dev[dev.len() / 2])
}

/// A confident target with an andes-rt prediction: a row of the RT ridge.
fn rt_fit_member(n: &Nominee) -> bool {
    confident(n) && !n.rt_missing && n.pred_rt.is_finite()
}

/// Gaussian elimination with partial pivoting; `None` when singular.
fn solve(mut a: [[f64; RT_TERMS]; RT_TERMS], mut b: [f64; RT_TERMS]) -> Option<[f64; RT_TERMS]> {
    for c in 0..RT_TERMS {
        let p = (c..RT_TERMS).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[p][c].abs() < PIVOT_MIN {
            return None;
        }
        a.swap(c, p);
        b.swap(c, p);
        for r in c + 1..RT_TERMS {
            let f = a[r][c] / a[c][c];
            for k in c..RT_TERMS {
                a[r][k] -= f * a[c][k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = [0.0; RT_TERMS];
    for r in (0..RT_TERMS).rev() {
        let s: f64 = (r + 1..RT_TERMS).map(|k| a[r][k] * x[k]).sum();
        x[r] = (b[r] - s) / a[r][r];
    }
    Some(x)
}

/// One isolation window of one run: its cycles in RT order.
pub(crate) struct Window<'a> {
    centre: f64,
    /// Isolation range (Th), lower and upper bound.
    iso: (f64, f64),
    mobility: Option<(u32, u32)>,
    pub(crate) rt: Vec<f64>,
    pub(crate) scans: Vec<&'a [(f64, f32)]>,
}

impl Window<'_> {
    /// Most intense peak within the tolerance of `mz` in cycle `c`: (intensity, peak id).
    fn peak(&self, c: usize, mz: f64) -> Option<(f32, u32)> {
        let peaks = self.scans[c];
        let tol = mz * FRAG_TRACE_PPM * 1e-6;
        let lo = peaks.partition_point(|p| p.0 < mz - tol);
        peaks[lo..]
            .iter()
            .enumerate()
            .take_while(|(_, p)| p.0 <= mz + tol)
            .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
            .map(|(k, p)| (p.1, (lo + k) as u32))
    }
}

#[cfg(test)]
impl<'a> Window<'a> {
    /// A window without isolation or mobility bounds, for tests of the per-window models.
    pub(crate) fn for_test(rt: Vec<f64>, scans: Vec<&'a [(f64, f32)]>) -> Self {
        Window {
            centre: 0.0,
            iso: (0.0, f64::INFINITY),
            mobility: None,
            rt,
            scans,
        }
    }
}

/// Intensities of a nominee's ions over its cycles, and the cycle range used.
struct Trace {
    c0: isize,
    valid: [bool; CYCLES],
    x: Vec<[f32; TOP_IONS]>,
    peak: Vec<[u32; TOP_IONS]>,
}

fn trace(w: &Window, n: &Nominee, rt_cal: f64, shift: f64) -> Trace {
    let c_first = w.rt.partition_point(|&t| t < rt_cal) as isize;
    let c0 = c_first - HALF_CYCLES as isize;
    let mut valid = [false; CYCLES];
    let mut x = vec![[0f32; TOP_IONS]; CYCLES];
    let mut peak = vec![[u32::MAX; TOP_IONS]; CYCLES];
    for (k, v) in valid.iter_mut().enumerate() {
        let c = c0 + k as isize;
        if c < 0 || c as usize >= w.rt.len() {
            continue;
        }
        *v = true;
        for (j, &(mz, _)) in n.ions.iter().enumerate() {
            if let Some((i, id)) = w.peak(c as usize, mz * shift) {
                x[k][j] = i;
                peak[k][j] = id;
            }
        }
    }
    Trace {
        c0,
        valid,
        x,
        peak,
    }
}

impl Trace {
    fn rt(&self, w: &Window, k: usize) -> f64 {
        let c = (self.c0 + k as isize).clamp(0, w.rt.len() as isize - 1);
        w.rt[c as usize]
    }
}

/// Per-run presence rates by ion rank: at the apex of confident targets (`q`) and over every
/// traced cycle of decoys (`p0`).
#[derive(Default, Clone)]
struct Rates {
    qn: [f64; TOP_IONS],
    qh: [f64; TOP_IONS],
    pn: [f64; TOP_IONS],
    ph: [f64; TOP_IONS],
}

impl Rates {
    fn add(&mut self, o: &Rates) {
        for j in 0..TOP_IONS {
            self.qn[j] += o.qn[j];
            self.qh[j] += o.qh[j];
            self.pn[j] += o.pn[j];
            self.ph[j] += o.ph[j];
        }
    }

    /// Per-ion log terms (present, absent).
    fn terms(&self) -> [(f64, f64); TOP_IONS] {
        let mut t = [(0.0, 0.0); TOP_IONS];
        for (j, tj) in t.iter_mut().enumerate() {
            let q = ((self.qh[j] + 1.0) / (self.qn[j] + 2.0)).clamp(Q_MIN, RATE_MAX);
            let p0 = ((self.ph[j] + 1.0) / (self.pn[j] + 2.0)).clamp(P0_MIN, RATE_MAX);
            *tj = ((q / p0).ln(), ((1.0 - q) / (1.0 - p0)).ln());
        }
        t
    }
}

fn collect_rates(w: &Window, n: &Nominee, t: &Trace) -> Rates {
    let mut r = Rates::default();
    let nj = n.ions.len();
    if confident(n) {
        let apex = (0..CYCLES)
            .filter(|&k| t.valid[k])
            .min_by(|&a, &b| {
                (t.rt(w, a) - n.apex_rt)
                    .abs()
                    .total_cmp(&(t.rt(w, b) - n.apex_rt).abs())
            });
        if let Some(a) = apex {
            for j in 0..nj {
                r.qn[j] += 1.0;
                r.qh[j] += f64::from(u8::from(t.x[a][j] > 0.0));
            }
        }
    }
    if null_member(n) {
        for k in (0..CYCLES).filter(|&k| t.valid[k]) {
            for j in 0..nj {
                r.pn[j] += 1.0;
                r.ph[j] += f64::from(u8::from(t.x[k][j] > 0.0));
            }
        }
    }
    r
}

/// A peak claimed by a nominee at its best cycle: (cycle spectrum, peak id) -> claimants.
type PeakKey = (u32, u32);

struct RowFeatures {
    vals: [f32; COLUMNS.len()],
    tr_full: [f32; CYCLES],
    tt: [f64; CYCLES],
    claims: Vec<(PeakKey, f64)>,
    /// Best cycle (index into the traced cycles).
    apex: usize,
    quant: RowQuant,
    /// Traced cycles of the quantity window; `None` without a contained peak.
    quant_window: Option<Range<usize>>,
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

fn row_features(
    w: &Window,
    win_scans: &[u32],
    n: &Nominee,
    t: &Trace,
    terms: &[(f64, f64); TOP_IONS],
    rt_cal: f64,
) -> RowFeatures {
    let nj = n.ions.len();
    let wt: Vec<f64> = n.ions.iter().map(|i| f64::from(i.1)).collect();
    let pres = |k: usize, j: usize| t.x[k][j] > 0.0;
    let mut lc = [f64::NEG_INFINITY; CYCLES];
    for (k, l) in lc.iter_mut().enumerate() {
        if t.valid[k] {
            *l = (0..nj)
                .map(|j| if pres(k, j) { terms[j].0 } else { terms[j].1 })
                .sum();
        }
    }
    let mut cs = 0usize;
    for k in 1..CYCLES {
        if lc[k] > lc[cs] {
            cs = k;
        }
    }
    let best = lc[cs];
    let vk: Vec<usize> = (0..CYCLES).filter(|&k| t.valid[k]).collect();
    let vm = vk.len().max(1) as f64;
    let mean = vk.iter().map(|&k| lc[k]).sum::<f64>() / vm;
    let sd = (vk.iter().map(|&k| (lc[k] - mean).powi(2)).sum::<f64>() / vm).sqrt() + Z_SD_EPS;
    let second = vk
        .iter()
        .filter(|&&k| k.abs_diff(cs) > FAR_CYCLES)
        .map(|&k| lc[k])
        .reduce(f64::max)
        .unwrap_or(LLR_FLOOR);
    let lo1 = cs.saturating_sub(1);
    let hi1 = (cs + 1).min(CYCLES - 1);
    let beside = |k: usize| if k != cs && t.valid[k] { lc[k] } else { LLR_FLOOR };
    let n6 = nj.min(TOP_SUBSET);
    let n_top6 = (0..n6).filter(|&j| pres(cs, j)).count();
    let n_all = (0..nj).filter(|&j| pres(cs, j)).count();
    let any6 = (0..n6)
        .filter(|&j| pres(lo1, j) || pres(cs, j) || pres(hi1, j))
        .count() as f64
        / n6.max(1) as f64;
    let wsum: f64 = wt.iter().sum();
    let expl = (0..nj).filter(|&j| pres(cs, j)).map(|j| wt[j]).sum::<f64>() / wsum.max(DIV_EPS);
    let xc: Vec<f64> = (0..nj).map(|j| f64::from(t.x[cs][j])).collect();
    let num: f64 = xc.iter().zip(&wt).map(|(a, b)| a * b).sum();
    let den = (xc.iter().map(|a| a * a).sum::<f64>() * wt.iter().map(|b| b * b).sum::<f64>())
        .sqrt()
        + DIV_EPS;
    let ang = num / den;
    let span: Vec<usize> = (0..=2 * COELUTE_HALF)
        .map(|d| (cs + d).saturating_sub(COELUTE_HALF).min(CYCLES - 1))
        .collect();
    let tr: Vec<f64> = span
        .iter()
        .map(|&k| (0..nj).map(|j| f64::from(t.x[k][j])).sum())
        .collect();
    let mut cor_sum = 0.0;
    let mut cor6 = 0usize;
    let mut apex_here = 0usize;
    for j in (0..nj).filter(|&j| pres(cs, j)) {
        let ion: Vec<f64> = span.iter().map(|&k| f64::from(t.x[k][j])).collect();
        let c = pearson(&ion, &tr);
        cor_sum += c;
        cor6 += usize::from(c > COELUTE_MIN_CORR);
        apex_here += usize::from(t.x[cs][j] >= t.x[lo1][j] && t.x[cs][j] >= t.x[hi1][j]);
    }
    let np = n_all.max(1) as f64;
    let tmax = tr.iter().copied().fold(0.0, f64::max);
    let width = tr.iter().filter(|&&v| v > 0.5 * tmax).count();
    let mut tt = [0f64; CYCLES];
    let mut tr_full = [0f32; CYCLES];
    for k in 0..CYCLES {
        tt[k] = t.rt(w, k);
        tr_full[k] = (0..nj).map(|j| t.x[k][j]).sum();
    }
    let claims = if t.valid[cs] {
        let spec = win_scans[(t.c0 + cs as isize) as usize];
        (0..nj)
            .filter(|&j| pres(cs, j))
            .map(|j| ((spec, t.peak[cs][j]), terms[j].0))
            .collect()
    } else {
        Vec::new()
    };
    let mut vals = [0f32; COLUMNS.len()];
    let set = |vals: &mut [f32; COLUMNS.len()], name: &str, v: f64| {
        if let Some(i) = COLUMNS.iter().position(|c| *c == name) {
            vals[i] = v as f32;
        }
    };
    set(&mut vals, "p2_llr", best);
    set(&mut vals, "p2_llr2", second);
    set(&mut vals, "p2_z", (best - mean) / sd);
    set(&mut vals, "p2_gap", best - second.max(LLR_FLOOR));
    set(&mut vals, "p2_n6", n_top6 as f64);
    set(&mut vals, "p2_nall", n_all as f64);
    set(&mut vals, "p2_any6", any6);
    set(&mut vals, "p2_expl", expl);
    set(&mut vals, "p2_angle", ang);
    set(&mut vals, "p2_cor", cor_sum / np);
    set(&mut vals, "p2_cor6", cor6 as f64 / np);
    set(&mut vals, "p2_apex", apex_here as f64 / np);
    set(&mut vals, "p2_lsum3", beside(lo1) + best + beside(hi1));
    set(&mut vals, "p2_logI", (tr.iter().sum::<f64>() + 1.0).log10());
    set(&mut vals, "p2_width", width as f64);
    set(&mut vals, "p2_rterr", (tt[cs] - rt_cal).abs());
    set(&mut vals, "p2_rtstar", tt[cs]);
    set(&mut vals, "p2_dia_llr", n.dia_llr);
    set(
        &mut vals,
        "p2_rt_pred_err",
        if n.xfer { 0.0 } else { (rt_cal - n.apex_rt).abs() },
    );
    set(&mut vals, "p2_ok", 1.0);
    RowFeatures {
        vals,
        tr_full,
        tt,
        claims,
        apex: cs,
        quant: RowQuant::default(),
        quant_window: None,
    }
}

/// Half-maximum width (cycles) of the summed fragment trace of a confident target around its
/// observed apex: the contiguous traced cycles at or above half the largest summed intensity
/// within `COELUTE_HALF` cycles of the cycle nearest `n.apex_rt`. `None` without signal there.
fn half_max_width(w: &Window, n: &Nominee, t: &Trace) -> Option<usize> {
    let near = (0..CYCLES).filter(|&k| t.valid[k]).min_by(|&a, &b| {
        (t.rt(w, a) - n.apex_rt)
            .abs()
            .total_cmp(&(t.rt(w, b) - n.apex_rt).abs())
    })?;
    let nj = n.ions.len();
    let sum = |k: usize| -> f64 { (0..nj).map(|j| f64::from(t.x[k][j])).sum() };
    let lo = near.saturating_sub(COELUTE_HALF);
    let hi = (near + COELUTE_HALF).min(CYCLES - 1);
    let top = (lo..=hi)
        .filter(|&k| t.valid[k])
        .max_by(|&a, &b| sum(a).total_cmp(&sum(b)).then(b.cmp(&a)))?;
    let half_max = 0.5 * sum(top);
    if half_max.is_nan() || half_max <= 0.0 {
        return None;
    }
    let (mut l, mut r) = (top, top);
    while l > 0 && t.valid[l - 1] && sum(l - 1) >= half_max {
        l -= 1;
    }
    while r + 1 < CYCLES && t.valid[r + 1] && sum(r + 1) >= half_max {
        r += 1;
    }
    Some(r - l + 1)
}

/// Half-width (cycles) of a run's quantity window from the half-maximum widths of its
/// confident targets: their median, at least [`QUANT_MIN_HALF_CYCLES`] and inside the traced
/// cycles; the fallback with fewer than [`QUANT_WIDTH_MIN_TARGETS`] widths.
fn quant_half_cycles(widths: &mut [usize]) -> usize {
    if widths.len() < QUANT_WIDTH_MIN_TARGETS {
        return QUANT_HALF_CYCLES_FALLBACK;
    }
    let m = widths.len() / 2;
    let (_, med, _) = widths.select_nth_unstable(m);
    (*med).clamp(QUANT_MIN_HALF_CYCLES, HALF_CYCLES - 1)
}

/// The traced cycles of the quantity window of a row whose best cycle is `best`: from `best`,
/// step to the largest `consensus` within ± `half` valid cycles until it stays (the apex), then
/// the apex ± `half`. `None` when `best` is not a valid cycle or the window leaves the valid
/// cycles (the peak is not contained).
fn quant_window(
    valid: &[bool; CYCLES],
    consensus: &[f32],
    best: usize,
    half: usize,
) -> Option<Range<usize>> {
    let ok = |k: usize| k < CYCLES && valid[k];
    if !ok(best) {
        return None;
    }
    let value = |k: usize| consensus.get(k).copied().unwrap_or(0.0);
    let mut apex = best;
    for _ in 0..CYCLES {
        let next = (apex.saturating_sub(half)..=(apex + half).min(CYCLES - 1))
            .filter(|&k| ok(k))
            .max_by(|&a, &b| value(a).total_cmp(&value(b)).then(b.cmp(&a)))
            .unwrap_or(apex);
        if value(next) <= value(apex) {
            break;
        }
        apex = next;
    }
    let lo = apex.checked_sub(half)?;
    let hi = apex + half;
    (ok(lo) && ok(hi)).then_some(lo..hi + 1)
}

/// Fragment quantity over `window` (the apex ± `half`): the summed intensity there of the first
/// [`QUANT_TOP_FRAGMENTS`] of `traces` (most intense predicted first), each less its background
/// and floored at 0 per cycle, and how many of them have signal above background. The background
/// of a trace is its median over the `valid` cycles at most
/// [`QUANT_BACKGROUND_REACH_HALF_WIDTHS`] × `half` from the apex and more than `half` outside the
/// window; 0 with fewer than [`QUANT_BACKGROUND_MIN_CYCLES`] such cycles. `None` (and 0) without signal above background.
pub(crate) fn fragment_quantity(
    traces: &[Vec<f32>],
    valid: &[bool],
    window: Range<usize>,
    half: usize,
) -> (Option<f64>, u8) {
    if window.is_empty() {
        return (None, 0);
    }
    let apex = (window.start + window.end - 1) / 2;
    let reach = QUANT_BACKGROUND_REACH_HALF_WIDTHS * half;
    let background_cycles: Vec<usize> = (apex.saturating_sub(reach)..=apex + reach)
        .filter(|&k| valid.get(k).copied().unwrap_or(false))
        .filter(|&k| k + half < window.start || k >= window.end + half)
        .collect();
    let mut total = 0f64;
    let mut n_frag = 0u8;
    for tr in traces.iter().take(QUANT_TOP_FRAGMENTS) {
        let Some(x) = tr.get(window.clone()) else {
            continue;
        };
        let background = if background_cycles.len() < QUANT_BACKGROUND_MIN_CYCLES {
            0.0
        } else {
            let mut b: Vec<f64> = background_cycles
                .iter()
                .filter_map(|&k| tr.get(k).map(|&v| f64::from(v)))
                .collect();
            model::stats::median_mut(&mut b, model::stats::Even::Mean).unwrap_or(0.0)
        };
        let above: f64 = x.iter().map(|&v| (f64::from(v) - background).max(0.0)).sum();
        if above > 0.0 && above.is_finite() {
            total += above;
            n_frag += 1;
        }
    }
    if n_frag == 0 || !total.is_finite() || total <= 0.0 {
        return (None, 0);
    }
    (Some(total), n_frag)
}

/// MS1 quantity over `window`: the summed intensity of the monoisotopic and M+1 traces;
/// `None` without monoisotopic signal there.
pub(crate) fn ms1_quantity(mono: &[f32], m1: &[f32], window: Range<usize>) -> Option<f64> {
    let end = window.end.min(mono.len()).min(m1.len());
    let w = window.start.min(end)..end;
    let sum = |v: &[f32]| -> f64 { v[w.clone()].iter().map(|&x| f64::from(x)).sum() };
    let (s_mono, s_m1) = (sum(mono), sum(m1));
    (s_mono > 0.0 && (s_mono + s_m1).is_finite()).then_some(s_mono + s_m1)
}

/// Fragment quantity of a traced row with its window (module doc).
fn row_quantity(t: &Trace, nj: usize, best: usize, half: usize) -> (RowQuant, Option<Range<usize>>) {
    let consensus: Vec<f32> = (0..CYCLES)
        .map(|k| (0..nj).map(|j| t.x[k][j]).sum())
        .collect();
    let Some(win) = quant_window(&t.valid, &consensus, best, half) else {
        return (RowQuant::default(), None);
    };
    let traces: Vec<Vec<f32>> = (0..nj)
        .map(|j| (0..CYCLES).map(|k| t.x[k][j]).collect())
        .collect();
    let (frag, n_frag) = fragment_quantity(&traces, &t.valid, win.clone(), half);
    (
        RowQuant {
            frag,
            n_frag,
            ms1: None,
        },
        Some(win),
    )
}

/// MS1 quantity of a row over the retention times `rt_lo_min..=rt_hi_min` (minutes).
fn ms1_row_quantity(ms1: &Ms1RunIndex, n: &Nominee, rt_lo_min: f64, rt_hi_min: f64) -> Option<f64> {
    let len = ms1.len();
    let first_at = |rt_s: f64, inclusive: bool| {
        let (mut lo, mut hi) = (0usize, len);
        while lo < hi {
            let mid = (lo + hi) / 2;
            let before = if inclusive {
                ms1.rt(mid) < rt_s
            } else {
                ms1.rt(mid) <= rt_s
            };
            if before {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    };
    let (s0, s1) = (first_at(rt_lo_min * 60.0, true), first_at(rt_hi_min * 60.0, false));
    if s0 >= s1 {
        return None;
    }
    let z = f64::from(n.charge.max(1));
    let at = |s: usize, mz: f64| -> f32 {
        let peaks = &ms1.scan(s).peaks;
        let tol = mz * MS1_TRACE_PPM * 1e-6;
        let lo = peaks.partition_point(|p| p.0 < mz - tol);
        peaks[lo..]
            .iter()
            .take_while(|p| p.0 <= mz + tol)
            .map(|p| p.1)
            .fold(0.0, f32::max)
    };
    let mono: Vec<f32> = (s0..s1).map(|s| at(s, n.prec_mz)).collect();
    let m1: Vec<f32> = (s0..s1).map(|s| at(s, n.prec_mz + ISOTOPE / z)).collect();
    ms1_quantity(&mono, &m1, 0..mono.len())
}

/// MS1 columns at the fragment apex: (log intensity, isotopes seen, apex offset, trace corr).
fn ms1_columns(ms1: &Ms1RunIndex, n: &Nominee, rt_star: f64, tt: &[f64], tr: &[f32]) -> [f64; 4] {
    let len = ms1.len();
    if len == 0 {
        return [0.0; 4];
    }
    let rt_s = rt_star * 60.0;
    let (mut lo, mut hi) = (0usize, len);
    while lo < hi {
        let mid = (lo + hi) / 2;
        if ms1.rt(mid) < rt_s {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let c1 = lo;
    let scans: Vec<usize> = (0..=2 * MS1_HALF)
        .map(|d| (c1 + d).saturating_sub(MS1_HALF).min(len - 1))
        .collect();
    let z = f64::from(n.charge.max(1));
    let at = |s: usize, mz: f64| -> f64 {
        let peaks = &ms1.scan(s).peaks;
        let tol = mz * MS1_TRACE_PPM * 1e-6;
        let lo = peaks.partition_point(|p| p.0 < mz - tol);
        peaks[lo..]
            .iter()
            .take_while(|p| p.0 <= mz + tol)
            .map(|p| f64::from(p.1))
            .fold(0.0, f64::max)
    };
    let mono: Vec<f64> = scans.iter().map(|&s| at(s, n.prec_mz)).collect();
    let mid = MS1_HALF - 1..=MS1_HALF + 1;
    let log_i = (mono[mid.clone()].iter().copied().fold(0.0, f64::max) + 1.0).log10();
    let iso = (0..MS1_ISOTOPES)
        .filter(|&k| {
            let mz = n.prec_mz + f64::from(k) * ISOTOPE / z;
            mid.clone().map(|d| at(scans[d], mz)).sum::<f64>() > 0.0
        })
        .count();
    let mut apex = 0usize;
    for k in 1..mono.len() {
        if mono[k] > mono[apex] {
            apex = k;
        }
    }
    let off = (apex as f64 - MS1_HALF as f64).abs();
    let frag: Vec<f64> = scans
        .iter()
        .map(|&s| {
            let t = ms1.rt(s) / 60.0;
            let k = tt.partition_point(|&x| x < t);
            if k == 0 {
                f64::from(tr[0])
            } else if k >= tt.len() {
                f64::from(tr[tt.len() - 1])
            } else {
                let (t0, t1) = (tt[k - 1], tt[k]);
                let f = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
                f64::from(tr[k - 1]) + f * (f64::from(tr[k]) - f64::from(tr[k - 1]))
            }
        })
        .collect();
    [log_i, iso as f64, off, pearson(&frag, &mono)]
}

/// Centre (scan index) of a mobility window; NaN without one.
pub(crate) fn mobility_centre(w: Option<(u32, u32)>) -> f64 {
    w.map_or(f64::NAN, |(lo, hi)| 0.5 * (f64::from(lo) + f64::from(hi)))
}

/// Where a window sits.
#[derive(Clone, Copy, Debug)]
struct Place {
    /// Isolation target (Th).
    centre: f64,
    /// Isolation range (Th), lower and upper bound.
    iso: (f64, f64),
    mobility: Option<(u32, u32)>,
    /// Index of the window's first spectrum.
    first: u32,
}

/// Isolation range of a window's spectrum: target minus its lower offset to target plus its
/// upper offset, or [`WINDOW_REACH_MZ`] either side of the target when the offsets are missing.
fn isolation_range(s: &Spectrum) -> (f64, f64) {
    let t = s.isolation_target_mz.unwrap_or(s.precursor_mz);
    let off = |o: Option<f64>| o.filter(|v| v.is_finite() && *v >= 0.0);
    match (off(s.isolation_lower_offset), off(s.isolation_upper_offset)) {
        (Some(lo), Some(hi)) if lo + hi > 0.0 => (t - lo, t + hi),
        _ => (t - WINDOW_REACH_MZ, t + WINDOW_REACH_MZ),
    }
}

/// The window of a nominee at precursor `mz` and `mobility` (NaN without one) among `places`,
/// one [`Place`] per window: of the windows whose isolation range holds `mz` (within the
/// precursor m/z accuracy `MS1_TRACE_PPM`, as the search widens a window by its precursor
/// tolerance; adjacent windows may leave sub-mTh gaps), the first by
/// (mobility outside the window, m/z distance to the centre, distance to the window's mobility
/// centre, first spectrum index). Without mobility on either side the mobility keys are equal,
/// so the nearest centre wins.
fn window_of(places: &[Place], mz: f64, mobility: f64) -> Option<usize> {
    let key = |p: &Place| {
        let (outside, dmob) = match p.mobility {
            Some((lo, hi)) if mobility.is_finite() => (
                mobility < f64::from(lo) || mobility > f64::from(hi),
                (mobility_centre(p.mobility) - mobility).abs(),
            ),
            _ => (false, 0.0),
        };
        (outside, (p.centre - mz).abs(), dmob, p.first)
    };
    places
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            let tol = mz * MS1_TRACE_PPM * 1e-6;
            mz >= p.iso.0 - tol && mz <= p.iso.1 + tol
        })
        .min_by(|(_, a), (_, b)| {
            let (ka, kb) = (key(a), key(b));
            ka.0.cmp(&kb.0)
                .then(ka.1.total_cmp(&kb.1))
                .then(ka.2.total_cmp(&kb.2))
                .then(ka.3.cmp(&kb.3))
        })
        .map(|(w, _)| w)
}

/// Pass-2 columns of every nominee. `ms1` is indexed by run; `shift_ppm` is each run's learned
/// fragment m/z shift.
pub(crate) fn extract(
    noms: &[Nominee],
    spectra: &[Spectrum],
    run_ends: &[usize],
    ms1: &[Option<Ms1RunIndex>],
    shift_ppm: &[f64],
) -> Pass2 {
    let n = noms.len();
    let mut cols = vec![vec![0f32; n]; COLUMNS.len()];
    let mut joint = vec![vec![0f32; n]; crate::dia_joint::COLUMNS.len()];
    let mut quant = vec![RowQuant::default(); n];
    let (_, windows) = crate::dia::windows(spectra, run_ends);
    let run_of = |s: usize| run_ends.partition_point(|&e| e <= s);
    let n_runs = run_ends.len().max(1);
    for run in 0..n_runs {
        let wins: Vec<(&Vec<u32>, Window)> = windows
            .iter()
            .filter(|v| !v.is_empty() && run_of(v[0] as usize) == run)
            .map(|v| {
                let s0 = &spectra[v[0] as usize];
                (
                    v,
                    Window {
                        centre: s0.isolation_target_mz.unwrap_or(s0.precursor_mz),
                        iso: isolation_range(s0),
                        mobility: s0.mobility_window,
                        rt: v
                            .iter()
                            .map(|&i| spectra[i as usize].rt_seconds.unwrap_or(0.0) / 60.0)
                            .collect(),
                        scans: v.iter().map(|&i| spectra[i as usize].peaks.as_slice()).collect(),
                    },
                )
            })
            .collect();
        let idx: Vec<usize> = (0..n).filter(|&i| noms[i].run as usize == run).collect();
        let refs: Vec<&Nominee> = idx.iter().map(|&i| &noms[i]).collect();
        let (rt_cal, rt_sigma) = calibrate_rt(&refs);
        let mut by_win: Vec<Vec<usize>> = vec![Vec::new(); wins.len()];
        let places: Vec<Place> = wins
            .iter()
            .map(|(v, w)| Place {
                centre: w.centre,
                iso: w.iso,
                mobility: w.mobility,
                first: v[0],
            })
            .collect();
        for (k, &i) in idx.iter().enumerate() {
            if let Some(w) = window_of(&places, noms[i].prec_mz, noms[i].mobility) {
                by_win[w].push(k);
            }
        }
        let shift = 1.0 + shift_ppm.get(run).copied().unwrap_or(0.0) * 1e-6;
        let per_win_rates: Vec<(Rates, Vec<usize>)> = by_win
            .par_iter()
            .enumerate()
            .map(|(w, ks)| {
                let mut r = Rates::default();
                let mut widths = Vec::new();
                for &k in ks.iter().filter(|&&k| confident(refs[k]) || null_member(refs[k])) {
                    let t = trace(&wins[w].1, refs[k], rt_cal[k], shift);
                    r.add(&collect_rates(&wins[w].1, refs[k], &t));
                    if confident(refs[k]) {
                        widths.extend(half_max_width(&wins[w].1, refs[k], &t));
                    }
                }
                (r, widths)
            })
            .collect();
        let mut widths: Vec<usize> = Vec::new();
        let rates = per_win_rates.iter().fold(Rates::default(), |mut a, (b, wb)| {
            a.add(b);
            widths.extend(wb);
            a
        });
        let terms = rates.terms();
        let n_widths = widths.len();
        let half = quant_half_cycles(&mut widths);
        if n_widths >= QUANT_WIDTH_MIN_TARGETS {
            log::info!(
                target: "dia",
                "pass 2 quantities, run {run}: window apex ± {half} cycles (median half-maximum \
                 width of {n_widths} confident targets)"
            );
        } else {
            log::info!(
                target: "dia",
                "pass 2 quantities, run {run}: window apex ± {half} cycles (fallback: {n_widths} \
                 confident targets with a peak, {QUANT_WIDTH_MIN_TARGETS} needed)"
            );
        }
        let per_win: Vec<Vec<(usize, RowFeatures)>> = by_win
            .par_iter()
            .enumerate()
            .map(|(w, ks)| {
                let (scans, win) = &wins[w];
                let mut out: Vec<(usize, RowFeatures)> = ks
                    .iter()
                    .map(|&k| {
                        let t = trace(win, refs[k], rt_cal[k], shift);
                        let mut f = row_features(win, scans, refs[k], &t, &terms, rt_cal[k]);
                        (f.quant, f.quant_window) = row_quantity(&t, refs[k].ions.len(), f.apex, half);
                        (k, f)
                    })
                    .collect();
                share_peaks(&mut out);
                if let Some(Some(m1)) = ms1.get(run) {
                    for (k, f) in out.iter_mut() {
                        let rt_star = f.vals[col_ix("p2_rtstar")] as f64;
                        let m = ms1_columns(m1, refs[*k], rt_star, &f.tt, &f.tr_full);
                        for (name, v) in ["p2_ms1_log", "p2_ms1_iso", "p2_ms1_off", "p2_ms1_corr"]
                            .iter()
                            .zip(m)
                        {
                            f.vals[col_ix(name)] = v as f32;
                        }
                        if let Some(win) = &f.quant_window {
                            let (lo, hi) = (f.tt[win.start], f.tt[win.end - 1]);
                            f.quant.ms1 = ms1_row_quantity(m1, refs[*k], lo, hi);
                        }
                    }
                }
                out
            })
            .collect();
        for (k, f) in per_win.into_iter().flatten() {
            let row = idx[k];
            for (c, v) in cols.iter_mut().zip(f.vals) {
                c[row] = v;
            }
            quant[row] = f.quant;
        }
        let win_refs: Vec<&Window> = wins.iter().map(|(_, w)| w).collect();
        let jf = crate::dia_joint::run(&crate::dia_joint::RunInput {
            wins: &win_refs,
            by_win: &by_win,
            noms: &refs,
            rt_cal: &rt_cal,
            shift,
            rt_sigma_min: rt_sigma,
        });
        for (k, vals) in jf {
            for (c, v) in joint.iter_mut().zip(vals) {
                c[idx[k]] = v;
            }
        }
    }
    Pass2 { cols, joint, quant }
}

/// Nominates, extracts and writes `<pin>.dia_pass2.parquet` beside the PIN. With
/// `nominees_out` (study mode), the nominees are also stored for the study's candidate
/// transfer (`dia_xfer`). The pool table and the kept predictions are only needed to
/// nominate, so they are taken by value and freed before the extraction.
pub(crate) fn run(
    pin: &Path,
    spectra: &[Spectrum],
    run_ends: &[usize],
    ms1: &[Option<Ms1RunIndex>],
    pool: crate::dia_llr::PoolScores,
    predicted: PredictedMs2,
    nominees_out: Option<crate::dia_xfer::NomineesOut>,
) -> Result<output::dia_pass2_table::Pass2Table, String> {
    let t0 = std::time::Instant::now();
    let noms = nominate(pin, spectra, run_ends, &pool.table, &pool.llr, &predicted)?;
    let shift_ppm = pool.shift_ppm.clone();
    let rt_cal = predicted.rt_cal.first().cloned().flatten();
    drop(pool);
    drop(predicted);
    crate::memlimit::release_free_memory();
    log::info!(
        target: "dia",
        "pass 2: nominated {} precursors [{:.2}s]",
        noms.len(),
        t0.elapsed().as_secs_f64()
    );
    let shift = shift_ppm.first().copied().unwrap_or(0.0);
    let rt_cal = rt_cal.as_ref();
    let (stored, p2) = std::thread::scope(|scope| {
        let stored = nominees_out.map(|out| {
            let noms = &noms;
            scope.spawn(move || {
                let t = std::time::Instant::now();
                crate::dia_xfer::write_nominees_table(&out, noms, shift, rt_cal)
                    .map(|()| t.elapsed().as_secs_f64())
            })
        });
        let p2 = extract(&noms, spectra, run_ends, ms1, &shift_ppm);
        let stored = stored.map(|h| {
            h.join()
                .unwrap_or_else(|_| Err("writing the nominee table panicked".to_string()))
        });
        (stored, p2)
    });
    if let Some(secs) = stored.transpose()? {
        log::info!(
            target: "dia",
            "pass 2: stored the nominees for the study, beside the extraction [{secs:.2}s]"
        );
    }
    let table = to_table(&noms, &p2);
    let path = pin.with_extension("dia_pass2.parquet");
    output::dia_pass2_table::write_pass2(&path, &table)?;
    let targets = noms.iter().filter(|n| !n.decoy).count();
    log::info!(
        target: "dia",
        "pass 2: {} nominees ({targets} targets), {} inside a window, {} [{:.2}s]",
        noms.len(),
        table.len(),
        path.display(),
        t0.elapsed().as_secs_f64()
    );
    Ok(table)
}

/// The rows inside an isolation window as the versioned pass-2 table.
pub(crate) fn to_table(noms: &[Nominee], p2: &Pass2) -> output::dia_pass2_table::Pass2Table {
    let ok = p2.col("p2_ok");
    let keep: Vec<usize> = (0..noms.len()).filter(|&i| ok[i] > 0.0).collect();
    let pick = |f: &dyn Fn(&Nominee) -> String| keep.iter().map(|&i| f(&noms[i])).collect();
    output::dia_pass2_table::Pass2Table {
        spec_id: pick(&|n| n.spec_id.to_string()),
        run: keep.iter().map(|&i| noms[i].run).collect(),
        decoy: keep.iter().map(|&i| noms[i].decoy).collect(),
        charge: keep.iter().map(|&i| noms[i].charge).collect(),
        peptide: pick(&|n| n.peptide.to_string()),
        proteins: pick(&|n| n.proteins.to_string()),
        seq: pick(&|n| n.seq.to_string()),
        prec_mz: keep.iter().map(|&i| noms[i].prec_mz).collect(),
        features: COLUMNS
            .iter()
            .zip(&p2.cols)
            .chain(crate::dia_joint::COLUMNS.iter().zip(&p2.joint))
            .map(|(n, v)| (n.to_string(), keep.iter().map(|&i| v[i]).collect()))
            .chain(noms.iter().any(|n| n.xfer).then(|| {
                (
                    XFER_COLUMN.to_string(),
                    keep.iter().map(|&i| f32::from(u8::from(noms[i].xfer))).collect(),
                )
            }))
            .collect(),
        quant_frag: keep.iter().map(|&i| p2.quant[i].frag).collect(),
        quant_ms1: keep.iter().map(|&i| p2.quant[i].ms1).collect(),
        quant_n_frag: keep.iter().map(|&i| p2.quant[i].n_frag).collect(),
    }
}

fn col_ix(name: &str) -> usize {
    COLUMNS.iter().position(|c| *c == name).unwrap_or(0)
}

/// Peaks claimed at their best cycle by several nominees of a window go to the claimant with
/// the highest `p2_llr`; the others count them as shared and lose their log term.
fn share_peaks(rows: &mut [(usize, RowFeatures)]) {
    let mut owner: HashMap<PeakKey, (f64, usize)> = HashMap::new();
    for (r, (_, f)) in rows.iter().enumerate() {
        let llr = f.vals[col_ix("p2_llr")] as f64;
        for (key, _) in &f.claims {
            let e = owner.entry(*key).or_insert((llr, r));
            if llr > e.0 {
                *e = (llr, r);
            }
        }
    }
    for (r, (_, f)) in rows.iter_mut().enumerate() {
        let (mut n_shared, mut lost) = (0usize, 0f64);
        for (key, term) in &f.claims {
            if owner.get(key).is_some_and(|o| o.1 != r) {
                n_shared += 1;
                lost += term;
            }
        }
        let llr = f.vals[col_ix("p2_llr")] as f64;
        f.vals[col_ix("p2_shared")] = (n_shared as f64 / f.claims.len().max(1) as f64) as f32;
        f.vals[col_ix("p2_lone")] = (llr - lost) as f32;
        f.vals[col_ix("p2_lone_gap")] = (-lost) as f32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nominee(decoy: bool, llr: f64, ions: Vec<(f64, f32)>) -> Nominee {
        Nominee {
            run: 0,
            charge: 2,
            decoy,
            prec_mz: 500.0,
            pred_rt: 10.0,
            rt_missing: false,
            raw_rt: f64::NAN,
            apex_rt: 10.0,
            dia_llr: llr,
            mobility: f64::NAN,
            ions,
            seq: "PEPTIDEK".into(),
            peptide: "-.PEPTIDEK.-".into(),
            proteins: "P1".into(),
            spec_id: "s_1_1".into(),
            xfer: false,
        }
    }

    fn window(n_cycles: usize, ions: &[f64], on: std::ops::Range<usize>) -> (Vec<Vec<(f64, f32)>>, Vec<f64>) {
        let mut scans = Vec::new();
        let mut rt = Vec::new();
        for c in 0..n_cycles {
            let mut p: Vec<(f64, f32)> = vec![(150.0, 10.0), (999.0, 10.0)];
            if on.contains(&c) {
                let h = 1.0 + (5.0 - (c as f64 - 33.0).abs()).max(0.0);
                for (j, &mz) in ions.iter().enumerate() {
                    p.push((mz, (1000.0 * h * (1.0 + j as f64)) as f32));
                }
            }
            p.sort_by(|a, b| a.0.total_cmp(&b.0));
            scans.push(p);
            rt.push(c as f64 * 0.02);
        }
        (scans, rt)
    }

    fn place(centre: f64, half: f64, mobility: Option<(u32, u32)>, first: u32) -> Place {
        Place {
            centre,
            iso: (centre - half, centre + half),
            mobility,
            first,
        }
    }

    #[test]
    fn window_of_picks_by_mz_then_mobility_and_breaks_ties_by_first_spectrum() {
        let slices = [
            place(500.0, 1.6, Some((600, 700)), 40),
            place(500.0, 1.6, Some((650, 750)), 10),
            place(500.0, 1.6, Some((700, 800)), 20),
            place(520.0, 1.6, Some((650, 750)), 5),
        ];
        assert_eq!(window_of(&slices, 500.2, 700.0), Some(1));
        assert_eq!(window_of(&slices, 500.2, 750.0), Some(2));
        assert_eq!(window_of(&slices, 500.2, 640.0), Some(0));
        assert_eq!(window_of(&slices, 500.2, 900.0), Some(2), "nearest slice when none holds it");
        assert_eq!(window_of(&slices, 500.2, f64::NAN), Some(1), "earliest of equal slices");
        let fenced = [place(500.0, 1.6, Some((600, 640)), 1), place(501.0, 1.6, Some((650, 700)), 2)];
        assert_eq!(window_of(&fenced, 500.2, 660.0), Some(1), "the window that holds it");
        assert_eq!(window_of(&fenced, 503.0, 660.0), None, "outside both isolation ranges");
        let astral = [place(512.0, 1.6, None, 30), place(500.0, 1.6, None, 7), place(503.0, 1.6, None, 3)];
        assert_eq!(window_of(&astral, 501.0, f64::NAN), Some(1));
        assert_eq!(window_of(&astral, 501.5, f64::NAN), Some(2), "a tie goes to the first spectrum");
        assert_eq!(window_of(&astral, 520.0, f64::NAN), None);
    }

    #[test]
    fn a_precursor_inside_a_wide_isolation_window_is_placed_in_it() {
        let wide = [place(500.0, 10.0, None, 1), place(520.0, 10.0, None, 2)];
        assert_eq!(window_of(&wide, 503.0, f64::NAN), Some(0), "3 Th from the centre of a 20 Th window");
        assert_eq!(window_of(&wide, 509.9, f64::NAN), Some(0));
        assert_eq!(window_of(&wide, 512.5, f64::NAN), Some(1));
        assert_eq!(window_of(&wide, 531.0, f64::NAN), None, "above every window");
        let pasef = [
            place(612.5, 12.5, Some((300, 420)), 4),
            place(612.5, 12.5, Some((420, 540)), 9),
            place(637.5, 12.5, Some((300, 420)), 2),
        ];
        assert_eq!(window_of(&pasef, 621.0, 500.0), Some(1), "the m/z window, then the mobility slice");
        assert_eq!(window_of(&pasef, 630.0, 350.0), Some(2));
        let asym = [Place { centre: 700.0, iso: (699.0, 712.0), mobility: None, first: 0 }];
        assert_eq!(window_of(&asym, 711.0, f64::NAN), Some(0), "asymmetric offsets");
        assert_eq!(window_of(&asym, 698.5, f64::NAN), None);
        let gap = [place(381.42325, 1.0, None, 1), place(383.42416, 1.0, None, 2)];
        assert_eq!(window_of(&gap, 382.4237, f64::NAN), Some(0), "inside the 0.9 mTh gap");
    }

    #[test]
    fn isolation_range_uses_the_offsets_and_falls_back_to_the_reach() {
        let mut s = Spectrum {
            isolation_target_mz: Some(600.0),
            isolation_lower_offset: Some(12.5),
            isolation_upper_offset: Some(12.5),
            ..Default::default()
        };
        assert_eq!(isolation_range(&s), (587.5, 612.5));
        s.isolation_upper_offset = None;
        assert_eq!(isolation_range(&s), (600.0 - WINDOW_REACH_MZ, 600.0 + WINDOW_REACH_MZ));
    }

    fn edge_features(on: std::ops::Range<usize>, rt_cal_cycle: usize, n_cycles: usize) -> RowFeatures {
        let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
        let (scans, rt) = window(n_cycles, &mzs, on);
        let w = Window {
            centre: 500.0,
            iso: (499.0, 501.0),
            mobility: None,
            rt: rt.clone(),
            scans: scans.iter().map(|s| s.as_slice()).collect(),
        };
        let ids: Vec<u32> = (0..n_cycles as u32).collect();
        let n = nominee(false, 20.0, mzs.iter().map(|&m| (m, 1.0)).collect());
        let rt_cal = rt[rt_cal_cycle];
        let t = trace(&w, &n, rt_cal, 1.0);
        let terms = [(2.0, -1.0); TOP_IONS];
        row_features(&w, &ids, &n, &t, &terms, rt_cal)
    }

    #[test]
    fn lsum3_stays_bounded_when_the_best_cycle_is_the_first_of_the_run() {
        let f = edge_features(0..1, 0, 120);
        let v = |c: &str| f.vals[col_ix(c)];
        assert_eq!(v("p2_llr"), 24.0);
        assert_eq!(v("p2_lsum3"), (LLR_FLOOR + 24.0 - 12.0) as f32);
        let f = edge_features(119..120, 119, 120);
        assert_eq!(f.vals[col_ix("p2_lsum3")], (LLR_FLOOR + 24.0 - 12.0) as f32, "last cycle of the run");
    }

    #[test]
    fn lsum3_counts_the_best_cycle_once_at_the_trace_edge() {
        let f = edge_features(0..120, 60, 200);
        let v = |c: &str| f.vals[col_ix(c)];
        assert_eq!(v("p2_llr"), 24.0);
        assert_eq!(v("p2_lsum3"), (LLR_FLOOR + 24.0 + 24.0) as f32, "best at trace cycle 0");
        let short = edge_features(0..3, 1, 3);
        let v = |c: &str| short.vals[col_ix(c)];
        assert_eq!(v("p2_llr2"), LLR_FLOOR as f32, "no far cycle in a three-cycle run");
        assert_eq!(v("p2_lsum3"), (LLR_FLOOR + 24.0 + 24.0) as f32);
    }

    #[test]
    fn rows_without_an_rt_prediction_stay_out_of_the_ridge() {
        let mut noms: Vec<Nominee> = (0..200)
            .map(|i| {
                let mut n = nominee(false, 20.0, vec![(300.0, 1.0)]);
                n.pred_rt = i as f64 * 0.1;
                n.apex_rt = 1.0 + 0.9 * n.pred_rt;
                n
            })
            .collect();
        let clean: Vec<&Nominee> = noms.iter().collect();
        let (cal_clean, sigma_clean) = calibrate_rt(&clean);
        noms.extend((0..200).map(|i| {
            let mut n = nominee(false, 20.0, vec![(300.0, 1.0)]);
            n.pred_rt = 10.0;
            n.apex_rt = i as f64 * 0.3;
            n.rt_missing = true;
            n
        }));
        let all: Vec<&Nominee> = noms.iter().collect();
        let (cal_all, sigma_all) = calibrate_rt(&all);
        assert_eq!(&cal_all[..200], &cal_clean[..]);
        assert_eq!(sigma_all, sigma_clean, "and out of the RT prior's spread");
    }

    #[test]
    fn pass2_llr_peaks_at_the_coeluting_cycle() {
        let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
        let (scans, rt) = window(120, &mzs, 32..36);
        let w = Window {
            centre: 500.0,
            iso: (499.0, 501.0),
            mobility: None,
            rt: rt.clone(),
            scans: scans.iter().map(|s| s.as_slice()).collect(),
        };
        let ids: Vec<u32> = (0..120).collect();
        let n = nominee(false, 20.0, mzs.iter().map(|&m| (m, 1.0)).collect());
        let rt_cal = rt[33];
        let t = trace(&w, &n, rt_cal, 1.0);
        let mut terms = [(2.0, -1.0); TOP_IONS];
        terms[11] = (1.0, -0.5);
        let f = row_features(&w, &ids, &n, &t, &terms, rt_cal);
        let v = |c: &str| f.vals[col_ix(c)];
        assert!(v("p2_rtstar") >= rt[30] as f32 && v("p2_rtstar") <= rt[36] as f32);
        assert_eq!(v("p2_n6"), 6.0);
        assert_eq!(v("p2_nall"), 12.0);
        assert!(v("p2_llr") > v("p2_llr2"));
        assert!(v("p2_cor") > 0.9);
    }

    #[test]
    fn pass2_shared_peak_goes_to_the_higher_llr() {
        let mk = |llr: f32, key: PeakKey| {
            let mut vals = [0f32; COLUMNS.len()];
            vals[col_ix("p2_llr")] = llr;
            RowFeatures {
                vals,
                tr_full: [0.0; CYCLES],
                tt: [0.0; CYCLES],
                claims: vec![(key, 2.0)],
                apex: 0,
                quant: RowQuant::default(),
                quant_window: None,
            }
        };
        let mut rows = vec![(0, mk(10.0, (5, 7))), (1, mk(4.0, (5, 7)))];
        share_peaks(&mut rows);
        assert_eq!(rows[0].1.vals[col_ix("p2_shared")], 0.0);
        assert_eq!(rows[1].1.vals[col_ix("p2_shared")], 1.0);
        assert_eq!(rows[1].1.vals[col_ix("p2_lone")], 2.0);
    }

    #[test]
    fn transferred_rows_stay_out_of_the_rt_ridge_and_the_rates() {
        let mut noms: Vec<Nominee> = (0..200)
            .map(|i| {
                let mut n = nominee(false, 20.0, vec![(300.0, 1.0)]);
                n.pred_rt = i as f64 * 0.1;
                n.apex_rt = 1.0 + 0.9 * n.pred_rt;
                n
            })
            .collect();
        let own: Vec<&Nominee> = noms.iter().collect();
        let (cal_own, sigma_own) = calibrate_rt(&own);
        noms.extend((0..200).map(|i| {
            let mut n = nominee(i % 2 == 0, 20.0, vec![(300.0, 1.0)]);
            n.pred_rt = i as f64 * 0.1;
            n.apex_rt = 50.0 - n.pred_rt;
            n.xfer = true;
            n
        }));
        let all: Vec<&Nominee> = noms.iter().collect();
        let (cal_all, sigma_all) = calibrate_rt(&all);
        assert_eq!(&cal_all[..200], &cal_own[..]);
        assert_eq!(sigma_all, sigma_own);

        let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
        let (scans, rt) = window(120, &mzs, 32..36);
        let w = Window {
            centre: 500.0,
            iso: (499.0, 501.0),
            mobility: None,
            rt: rt.clone(),
            scans: scans.iter().map(|s| s.as_slice()).collect(),
        };
        let ions: Vec<(f64, f32)> = mzs.iter().map(|&m| (m, 1.0)).collect();
        for decoy in [false, true] {
            let mut n = nominee(decoy, 20.0, ions.clone());
            n.apex_rt = rt[33];
            let t = trace(&w, &n, rt[33], 1.0);
            let r = collect_rates(&w, &n, &t);
            assert!(r.qn.iter().chain(&r.pn).sum::<f64>() > 0.0);
            n.xfer = true;
            let r = collect_rates(&w, &n, &t);
            assert_eq!(r.qn.iter().chain(&r.pn).sum::<f64>(), 0.0);
            assert!(!confident(&n) && !null_member(&n));
        }
    }

    #[test]
    fn xfer_column_marks_transferred_rows() {
        let mut noms = vec![
            nominee(false, 20.0, vec![(300.0, 1.0)]),
            nominee(true, 20.0, vec![(300.0, 1.0)]),
            nominee(false, TRANSFER_DIA_LLR, vec![(300.0, 1.0)]),
        ];
        let p2 = Pass2 {
            cols: (0..COLUMNS.len())
                .map(|c| if COLUMNS[c] == "p2_ok" { vec![1.0; 3] } else { vec![0.0; 3] })
                .collect(),
            joint: vec![vec![0.0; 3]; crate::dia_joint::COLUMNS.len()],
            quant: vec![
                RowQuant::default(),
                RowQuant {
                    frag: Some(5.0),
                    n_frag: 2,
                    ms1: Some(7.0),
                },
                RowQuant::default(),
            ],
        };
        assert!(to_table(&noms, &p2).feature(XFER_COLUMN).is_none());
        noms[2].xfer = true;
        let t = to_table(&noms, &p2);
        assert_eq!(t.feature(XFER_COLUMN).unwrap(), &[0.0, 0.0, 1.0]);
        assert_eq!(t.features.last().unwrap().0, XFER_COLUMN);
        assert_eq!(t.quant_frag, vec![None, Some(5.0), None]);
        assert_eq!(t.quant_ms1, vec![None, Some(7.0), None]);
        assert_eq!(t.quant_n_frag, vec![0, 2, 0]);
        assert!(t.features.iter().all(|(n, _)| !n.starts_with("quant")));
    }

    /// Gaussian elution profile over `cycles` cycles, apex at `apex`, sd `sd` cycles.
    fn gauss(cycles: usize, apex: f64, sd: f64) -> Vec<f32> {
        (0..cycles)
            .map(|k| (-((k as f64 - apex).powi(2)) / (2.0 * sd * sd)).exp() as f32)
            .collect()
    }

    fn area(t: &[f32], w: Range<usize>) -> f64 {
        t[w].iter().map(|&v| f64::from(v)).sum()
    }

    #[test]
    fn fragment_quantity_sums_the_six_most_intense_predicted_fragments() {
        let g = gauss(41, 20.0, 2.0);
        let valid = [true; 41];
        let mut traces: Vec<Vec<f32>> = (0..6)
            .map(|j| g.iter().map(|v| v * 1000.0 * (j + 1) as f32).collect())
            .collect();
        traces.push(g.iter().map(|v| v * 1e6).collect());
        let w = 17..24;
        let (q, n) = fragment_quantity(&traces, &valid, w.clone(), 3);
        assert_eq!(n, 6);
        let want: f64 = traces[..6].iter().map(|t| area(t, w.clone())).sum();
        let got = q.unwrap();
        assert!((got - want).abs() < 1e-6 * want, "{got} vs {want}");
    }

    #[test]
    fn fragment_quantity_subtracts_each_fragment_background() {
        let g = gauss(41, 20.0, 2.0);
        let valid = [true; 41];
        let floor = [50.0f32, 0.0, 400.0];
        let traces: Vec<Vec<f32>> = [300.0f32, 900.0, 150.0]
            .iter()
            .zip(floor)
            .map(|(a, b)| g.iter().map(|v| v * a + b).collect())
            .collect();
        let w = 17..24;
        let (q, n) = fragment_quantity(&traces, &valid, w.clone(), 3);
        assert_eq!(n, 3);
        let want: f64 = [300.0f64, 900.0, 150.0]
            .iter()
            .map(|a| a * area(&g, w.clone()))
            .sum();
        assert!((q.unwrap() - want).abs() < 1e-3 * want, "{q:?} vs {want}");
        let mut few = [false; 41];
        few[14..27].fill(true);
        let (q, _) = fragment_quantity(&traces, &few, w.clone(), 3);
        let raw: f64 = traces.iter().map(|t| area(t, w.clone())).sum();
        assert!((q.unwrap() - raw).abs() < 1e-6 * raw, "without background cycles nothing is subtracted");
    }

    #[test]
    fn quantity_missing_without_signal() {
        let valid = [true; 41];
        let zero = vec![vec![0f32; 41]; 6];
        assert_eq!(fragment_quantity(&zero, &valid, 17..24, 3), (None, 0));
        let flat = vec![vec![70f32; 41]; 6];
        assert_eq!(fragment_quantity(&flat, &valid, 17..24, 3), (None, 0), "background only");
        let g = gauss(41, 20.0, 2.0);
        let one = vec![g.clone()];
        assert_eq!(fragment_quantity(&one, &valid, 20..20, 3), (None, 0), "empty window");
        assert_eq!(fragment_quantity(&[], &valid, 17..24, 3), (None, 0), "no fragments");
        let g = gauss(21, 10.0, 2.0);
        assert_eq!(ms1_quantity(&[0.0; 9], &[0.0; 9], 0..9), None);
        assert_eq!(ms1_quantity(&[0.0; 9], &[5.0; 9], 0..9), None, "M+1 alone is not a quantity");
        let mono: Vec<f32> = g.iter().map(|v| v * 1e4).collect();
        let m1: Vec<f32> = g.iter().map(|v| v * 5e3).collect();
        let q = ms1_quantity(&mono, &m1, 6..15).unwrap();
        let want = area(&mono, 6..15) + area(&m1, 6..15);
        assert!((q - want).abs() < 1e-6 * want);
        assert_eq!(ms1_quantity(&mono, &m1, 30..40), None, "window past the trace");
    }

    #[test]
    fn transferred_row_at_window_edge_is_missing() {
        let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
        let rt: Vec<f64> = (0..40).map(|c| c as f64 * 0.02).collect();
        let scans: Vec<Vec<(f64, f32)>> = (0..40)
            .map(|c| {
                let h = (c as usize).saturating_sub(29) as f32 * 1000.0;
                mzs.iter().filter(|_| h > 0.0).map(|&m| (m, h)).collect()
            })
            .collect();
        let w = Window {
            centre: 500.0,
            iso: (499.0, 501.0),
            mobility: None,
            rt: rt.clone(),
            scans: scans.iter().map(|s| s.as_slice()).collect(),
        };
        let ions: Vec<(f64, f32)> = mzs.iter().map(|&m| (m, 1.0)).collect();
        let mut n = nominee(false, TRANSFER_DIA_LLR, ions);
        n.xfer = true;
        let terms = [(2.0, -1.0); TOP_IONS];
        let ids: Vec<u32> = (0..40).collect();
        for rt_cal in [rt[39], rt[39] + 100.0] {
            let t = trace(&w, &n, rt_cal, 1.0);
            let f = row_features(&w, &ids, &n, &t, &terms, rt_cal);
            let (q, win) = row_quantity(&t, n.ions.len(), f.apex, 3);
            assert_eq!(q, RowQuant::default(), "rt_cal {rt_cal}");
            assert!(win.is_none());
        }
        let (scans, rt) = window(120, &mzs, 28..39);
        let w = Window {
            centre: 500.0,
            iso: (499.0, 501.0),
            mobility: None,
            rt: rt.clone(),
            scans: scans.iter().map(|s| s.as_slice()).collect(),
        };
        let ids: Vec<u32> = (0..120).collect();
        let t = trace(&w, &n, rt[33], 1.0);
        let f = row_features(&w, &ids, &n, &t, &terms, rt[33]);
        let (q, win) = row_quantity(&t, n.ions.len(), f.apex, 3);
        let win = win.unwrap();
        assert_eq!(win.len(), 7);
        assert_eq!(t.rt(&w, win.start + 3), rt[33], "the window is centred on the trace maximum");
        assert_eq!(q.n_frag, 6);
        assert!(q.frag.is_some_and(|v| v > 0.0));
    }

    #[test]
    fn run_peak_width_is_the_median_half_maximum_width() {
        let mut few = vec![9usize; QUANT_WIDTH_MIN_TARGETS - 1];
        assert_eq!(quant_half_cycles(&mut few), QUANT_HALF_CYCLES_FALLBACK);
        let mut w: Vec<usize> = (0..101).map(|i| 1 + i % 9).collect();
        assert_eq!(quant_half_cycles(&mut w), 5);
        let mut narrow = vec![1usize; 60];
        assert_eq!(quant_half_cycles(&mut narrow), QUANT_MIN_HALF_CYCLES);
        let mut wide = vec![500usize; 60];
        assert_eq!(quant_half_cycles(&mut wide), HALF_CYCLES - 1);

        let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
        let (scans, rt) = window(120, &mzs, 28..39);
        let win = Window {
            centre: 500.0,
            iso: (499.0, 501.0),
            mobility: None,
            rt: rt.clone(),
            scans: scans.iter().map(|s| s.as_slice()).collect(),
        };
        let mut n = nominee(false, 20.0, mzs.iter().map(|&m| (m, 1.0)).collect());
        n.apex_rt = rt[33];
        let t = trace(&win, &n, rt[33], 1.0);
        assert_eq!(half_max_width(&win, &n, &t), Some(7));
    }

    #[test]
    fn ridge_recovers_a_linear_rt_map() {
        let noms: Vec<Nominee> = (0..200)
            .map(|i| {
                let mut n = nominee(false, 20.0, vec![(300.0, 1.0)]);
                n.pred_rt = i as f64 * 0.1;
                n.apex_rt = 1.0 + 0.9 * n.pred_rt;
                n
            })
            .collect();
        let refs: Vec<&Nominee> = noms.iter().collect();
        let (cal, _) = calibrate_rt(&refs);
        for (n, c) in noms.iter().zip(&cal) {
            assert!((n.apex_rt - c).abs() < 0.05, "{} vs {}", n.apex_rt, c);
        }
    }
}
