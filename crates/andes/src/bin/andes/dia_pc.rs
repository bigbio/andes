//! `--dia` Phase B: peptide-centric DIA scoring, after the spectrum-centric search.
//!
//! The spectrum-centric search (Phase A) asks each scan which peptide it holds; a
//! weak precursor loses that ranking to the other forms of a crowded window.
//! This pass asks the opposite question, as DIA-NN and OpenSWATH do: for each
//! library precursor, is there a chromatographic peak group where its fragments
//! co-elute at its predicted retention time?
//!
//! 1. **Library.** Every target and decoy peptidoform of the database at every
//!    searched charge whose m/z falls in an isolation window of the run.
//! 2. **Retention time.** andes-rt predicts each form; the run's calibration
//!    is fitted on the Phase A PIN (the best row of each spectrum that is a
//!    target at q <= 1 %), and its residual width sets the extraction window.
//! 3. **Fragments.** andes-ms2 predicts b/y 1+/2+ intensities; the
//!    [`N_FRAGMENTS`] most intense ions (length >= 3, m/z in the scan range)
//!    are the precursor's assay.
//! 4. **Extraction.** Each fragment is traced, at [`FRAG_PPM`], over the scans
//!    of the precursor's isolation window within the RT window; the
//!    monoisotopic precursor is traced in MS1.
//! 5. **Peak groups.** Candidate apexes are local maxima of the
//!    predicted-intensity-weighted sum of the smoothed fragment traces; the
//!    [`PEAK_GROUPS`] highest become rows of the precursor in the PIN, sharing its
//!    `ScanNr`, so target-decoy competition (`-Y`) keeps one per precursor.
//! 6. **Scores** per peak group ([`COLUMNS`]): fragment co-elution
//!    (correlation with the summed trace, how many correlate), spectral angle
//!    of the apex and of the integrated intensities against the prediction,
//!    matched fragments at the apex and their mass error, MS1 co-elution,
//!    error and isotope fit, the retention-time error in units of the
//!    calibration width, and intensities.
//!
//! The PIN `<stem>.dia.pin` is precursor-level: one target-decoy competition
//! per precursor (peptidoform, charge), which is the level DIA results are
//! reported at.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use model::isotope::averagine_isotope_envelope;
use model::mass::ISOTOPE;
use model::scan::Ms1Scan;
use model::Spectrum;
use onnx_infer::predictors::calibration::RtCalibration;
use onnx_infer::predictors::form::{parse_peptide, Form, Mass, H2O, PROTON};
use onnx_infer::predictors::{NativeMs2, NativeRt, Predictor};
use quant::Ms1RunIndex;
use rayon::prelude::*;
use search::candidate_gen::Candidate;

use crate::predicted::{read_pin, ANDES_MS2_WEIGHTS, ANDES_RT_WEIGHTS, NONE};

/// Fragments in a precursor's assay.
const N_FRAGMENTS: usize = 12;
/// Fewest assay fragments a precursor needs.
const MIN_FRAGMENTS: usize = 4;
/// Shortest b/y ion (residues).
const MIN_ION_LEN: usize = 3;
/// Fragment and MS1 extraction tolerance (ppm) before calibration; after it,
/// [`CAL_SIGMAS`] robust widths around the measured offset, at least
/// [`MIN_TOL_PPM`].
const FRAG_PPM: f64 = 10.0;
const MS1_PPM: f64 = 10.0;
const CAL_SIGMAS: f64 = 3.0;
const MIN_TOL_PPM: f64 = 3.0;
/// Precursors (target, in Phase A at q <= 1 %) a mass calibration is learned on.
const CAL_PRECURSORS: usize = 2000;
/// Peak groups kept per precursor, of the [`CANDIDATE_APEXES`] most intense
/// local maxima scored.
const PEAK_GROUPS: usize = 3;
const CANDIDATE_APEXES: usize = 8;
/// RT window half-width in calibration widths, and its floor (minutes).
/// ponytail: one width for the whole gradient; a width per RT bin if the
/// calibration residuals turn out to vary along it.
const RT_WINDOW_WIDTHS: f64 = 4.0;
const RT_WINDOW_FLOOR_MIN: f64 = 0.5;
/// Scans either side of an apex that a peak group's correlations read. Astral
/// peaks are ~3 s at half height on a 1.65 s cycle: ±3 scans hold the peak
/// and its shoulders.
const PEAK_HALF_SCANS: usize = 3;

/// The columns of the Phase B PIN, in order.
pub(crate) const COLUMNS: [&str; 23] = [
    "pc_corr_mean",
    "pc_corr_n06",
    "pc_corr_top3",
    "pc_apex_angle",
    "pc_area_angle",
    "pc_apex_matched",
    "pc_apex_matched_frac",
    "pc_apex_ppm",
    "pc_log_apex",
    "pc_log_area",
    "pc_ms1_corr",
    "pc_ms1_found",
    "pc_ms1_abs_ppm",
    "pc_ms1_iso_cos",
    "pc_ms1_log",
    "pc_rt_err",
    "pc_peak_rank",
    "pc_n_fragments",
    "pc_apex_coincide",
    "pc_apex_top3",
    "pc_apex_vs_bg",
    "pc_xcorr_shift",
    "pc_apex_pred_dot",
];

/// Mass calibration of a run: signed offsets (observed - theoretical, ppm)
/// and extraction tolerances (ppm) of fragments and of the MS1 precursor.
#[derive(Clone, Copy, Debug)]
struct MassCal {
    frag_offset: f64,
    frag_tol: f64,
    ms1_offset: f64,
    ms1_tol: f64,
}

impl Default for MassCal {
    fn default() -> Self {
        Self {
            frag_offset: 0.0,
            frag_tol: FRAG_PPM,
            ms1_offset: 0.0,
            ms1_tol: MS1_PPM,
        }
    }
}

/// Median and robust width (1.4826 MAD) of `v`; `None` with under 50 values.
fn median_sigma(mut v: Vec<f64>) -> Option<(f64, f64)> {
    if v.len() < 50 {
        return None;
    }
    v.sort_by(f64::total_cmp);
    let m = v[v.len() / 2];
    let mut d: Vec<f64> = v.iter().map(|x| (x - m).abs()).collect();
    d.sort_by(f64::total_cmp);
    Some((m, 1.4826 * d[d.len() / 2]))
}

/// A precursor of the library.
struct Precursor {
    /// Index into the candidate list (peptide text, decoy flag, protein).
    cand: u32,
    form: u32,
    charge: u8,
    mz: f64,
}

/// One isolation window of a run: its scans in RT order.
struct Window {
    lo: f64,
    hi: f64,
    scans: Vec<u32>,
}

/// Inputs of one run.
pub(crate) struct Run {
    /// `spectra` index range of the run.
    pub range: std::ops::Range<usize>,
    pub ms1: Vec<Ms1Scan>,
    pub stem: String,
}

/// Builds the Phase B PIN `out` for `runs` from the Phase A PIN `pin`
/// (RT anchors) and the database `cands` (all target and decoy peptidoforms).
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_pin(
    out: &Path,
    pin: &Path,
    all_spectra: &[Spectrum],
    runs: Vec<Run>,
    cands: &[Candidate],
    accession: &(dyn Fn(&Candidate) -> String + Sync),
    charges: std::ops::RangeInclusive<u8>,
    ms2_conditions: (Option<f32>, &str, &str),
) -> Result<(), String> {
    let t0 = Instant::now();
    let rt_model = NativeRt::from_bytes(ANDES_RT_WEIGHTS).map_err(|e| e.to_string())?;
    let (nce, act, inst) = ms2_conditions;
    let ms2_model =
        NativeMs2::from_bytes(ANDES_MS2_WEIGHTS, nce, act, inst).map_err(|e| e.to_string())?;

    // Parsed forms, one per distinct peptide text.
    let mut forms: Vec<(Form, Mass)> = Vec::new();
    let mut form_of: HashMap<String, u32> = HashMap::new();
    let mut cand_form: Vec<u32> = Vec::with_capacity(cands.len());
    for c in cands {
        let text = c.peptide.to_string();
        let f = *form_of
            .entry(text.clone())
            .or_insert_with(|| match parse_peptide(&text) {
                Some(fm) => {
                    forms.push(fm);
                    (forms.len() - 1) as u32
                }
                None => NONE,
            });
        cand_form.push(f);
    }

    // Predicted RT per form.
    let form_refs: Vec<&Form> = forms.iter().map(|f| &f.0).collect();
    let pred_rt: Vec<f32> = form_refs
        .par_chunks(4096)
        .flat_map_iter(|chunk| {
            rt_model
                .rt(&mut [], chunk)
                .unwrap_or_else(|_| vec![f32::NAN; chunk.len()])
        })
        .collect();

    // RT anchors from the Phase A PIN: best row per spectrum, target, q <= 1 %.
    let table = read_pin(pin, all_spectra, &[])?;
    let n_runs = runs.len();
    let run_of = |spec: u32| {
        runs.iter()
            .position(|r| r.range.contains(&(spec as usize)))
            .unwrap_or(0)
    };
    let anchors = crate::predicted::anchor_rows(&table.rows, n_runs, &run_of);
    let mut out_rows: Vec<String> = Vec::new();
    let mut n_prec_total = 0usize;
    for (ri, run) in runs.into_iter().enumerate() {
        let pts: Vec<(f64, f64)> = anchors[ri]
            .iter()
            .filter_map(|&i| {
                let r = &table.rows[i];
                let f = &table.forms.get(r.form as usize)?.0;
                let p = rt_model.rt(&mut [], &[f]).ok()?[0] as f64;
                let o = all_spectra[r.spec as usize].rt_seconds? / 60.0;
                Some((p, o))
            })
            .collect();
        let cal = RtCalibration::fit(&pts).ok_or("DIA Phase B: no RT calibration")?;
        // Width out of sample: fit on half the anchors, measure on the other
        // half (an isotonic fit's own residuals are far too small).
        let half = |k: usize| -> Vec<(f64, f64)> {
            pts.iter()
                .enumerate()
                .filter(|(i, _)| i % 2 == k)
                .map(|(_, p)| *p)
                .collect()
        };
        let mut res: Vec<f64> = Vec::new();
        for k in 0..2 {
            if let Some(c) = RtCalibration::fit(&half(k)) {
                res.extend(half(1 - k).iter().map(|(p, o)| (o - c.apply(*p)).abs()));
            }
        }
        if res.is_empty() {
            return Err("DIA Phase B: too few RT anchors".into());
        }
        res.sort_by(f64::total_cmp);
        let width = 1.4826 * res[res.len() / 2];
        let half_win = (RT_WINDOW_WIDTHS * width).max(RT_WINDOW_FLOOR_MIN);
        eprintln!(
            "DIA Phase B [{}]: RT calibration on {} anchors, width {:.3} min, window ±{:.2} min",
            run.stem,
            pts.len(),
            width,
            half_win
        );

        let windows = isolation_windows(all_spectra, run.range.clone());
        // Library: (cand, charge) whose m/z falls in a window.
        let mut precs: Vec<Precursor> = Vec::new();
        for (ci, c) in cands.iter().enumerate() {
            let f = cand_form[ci];
            if f == NONE {
                continue;
            }
            for z in charges.clone() {
                let mz = c.peptide.mass() / z as f64 + PROTON;
                if window_of(&windows, mz).is_some() {
                    precs.push(Precursor {
                        cand: ci as u32,
                        form: f,
                        charge: z,
                        mz,
                    });
                }
            }
        }
        // One precursor per (form, charge, decoy): the database repeats a
        // peptide once per protein.
        let mut seen: HashMap<(u32, u8, bool), ()> = HashMap::new();
        precs.retain(|p| {
            seen.insert((p.form, p.charge, cands[p.cand as usize].is_decoy), ())
                .is_none()
        });
        n_prec_total += precs.len();

        // Predicted fragments per precursor.
        let items: Vec<(&Form, u8)> = precs
            .iter()
            .map(|p| (&forms[p.form as usize].0, p.charge))
            .collect();
        let pred_ms2: Vec<Vec<f32>> = items
            .par_chunks(1024)
            .flat_map_iter(|chunk| {
                ms2_model
                    .ms2(&mut [], chunk)
                    .unwrap_or_else(|_| vec![Vec::new(); chunk.len()])
            })
            .collect();

        let ms1 = (!run.ms1.is_empty()).then(|| Ms1RunIndex::new(run.ms1));
        // Precursor -> its extraction inputs (window scans, assay, RT).
        let prepare =
            |p: &Precursor, pred: &[f32]| -> Option<(Vec<&Spectrum>, Vec<(f64, f64)>, f64)> {
                let form = &forms[p.form as usize];
                let rt_pred = cal.apply(pred_rt[p.form as usize] as f64);
                let w = window_of(&windows, p.mz)?;
                let assay = assay(&form.1, p.charge, pred, ms2_model.ions());
                if assay.len() < MIN_FRAGMENTS || !rt_pred.is_finite() {
                    return None;
                }
                let win = &windows[w];
                let lo = win
                    .scans
                    .partition_point(|&s| rt_min(all_spectra, s) < rt_pred - half_win);
                let hi = win
                    .scans
                    .partition_point(|&s| rt_min(all_spectra, s) <= rt_pred + half_win);
                if hi <= lo + 2 {
                    return None;
                }
                let scans = win.scans[lo..hi]
                    .iter()
                    .map(|&s| &all_spectra[s as usize])
                    .collect();
                Some((scans, assay, rt_pred))
            };

        // Mass calibration: the apex errors of the best peak group of target
        // precursors Phase A identified, extracted at the default tolerance.
        let confident: std::collections::HashSet<(&Form, u8)> = anchors[ri]
            .iter()
            .filter_map(|&i| {
                let r = &table.rows[i];
                Some((&table.forms.get(r.form as usize)?.0, r.charge))
            })
            .collect();
        let cal_idx: Vec<usize> = precs
            .iter()
            .enumerate()
            .filter(|(_, p)| {
                !cands[p.cand as usize].is_decoy
                    && confident.contains(&(&forms[p.form as usize].0, p.charge))
            })
            .map(|(i, _)| i)
            .take(CAL_PRECURSORS)
            .collect();
        let errs: Vec<(Vec<f64>, Option<f64>)> = cal_idx
            .par_iter()
            .filter_map(|&i| {
                let p = &precs[i];
                let (scans, assay, _) = prepare(p, &pred_ms2[i])?;
                Some(apex_errors(&scans, &assay, p.mz, ms1.as_ref()))
            })
            .collect();
        let frag = median_sigma(errs.iter().flat_map(|e| e.0.iter().copied()).collect());
        let ms1e = median_sigma(errs.iter().filter_map(|e| e.1).collect());
        let mcal = MassCal {
            frag_offset: frag.map_or(0.0, |f| f.0),
            frag_tol: frag.map_or(FRAG_PPM, |f| {
                (CAL_SIGMAS * f.1).clamp(MIN_TOL_PPM, FRAG_PPM)
            }),
            ms1_offset: ms1e.map_or(0.0, |f| f.0),
            ms1_tol: ms1e.map_or(MS1_PPM, |f| (CAL_SIGMAS * f.1).clamp(MIN_TOL_PPM, MS1_PPM)),
        };
        eprintln!(
            "DIA Phase B [{}]: mass calibration on {} precursors: fragments {:+.2} ppm, tol {:.2} ppm; MS1 {:+.2} ppm, tol {:.2} ppm",
            run.stem,
            errs.len(),
            mcal.frag_offset,
            mcal.frag_tol,
            mcal.ms1_offset,
            mcal.ms1_tol
        );

        let rows: Vec<Vec<String>> = precs
            .par_iter()
            .zip(pred_ms2.par_iter())
            .map(|(p, pred)| {
                let form = &forms[p.form as usize];
                let Some((scans, assay, rt_pred)) = prepare(p, pred) else {
                    return Vec::new();
                };
                let groups = score_precursor(
                    &scans,
                    &assay,
                    p.mz,
                    p.charge,
                    &form.1,
                    rt_pred,
                    width,
                    ms1.as_ref(),
                    &mcal,
                );
                let c = &cands[p.cand as usize];
                let label = if c.is_decoy { -1 } else { 1 };
                let id = format!("{}_{}_{}", c.peptide, p.charge, label);
                // ScanNr: one competition group per precursor (form, charge),
                // so `-Y` keeps its best peak group. A decoy is its own form,
                // so targets and decoys are counted, not paired.
                let scan_nr = p.form as u64 * 8 + p.charge as u64;
                let onehot: String = charges
                    .clone()
                    .map(|z| if z == p.charge { "\t1" } else { "\t0" })
                    .collect();
                groups
                    .into_iter()
                    .enumerate()
                    .map(|(k, f)| {
                        let mut s = format!(
                            "{}/{id}_{k}\t{label}\t{scan_nr}\t{:.5}{onehot}",
                            run.stem,
                            c.peptide.mass(),
                        );
                        for v in f {
                            s.push('\t');
                            s.push_str(&format!("{v:.5}"));
                        }
                        s.push('\t');
                        s.push_str(&c.peptide.to_string());
                        s.push('\t');
                        s.push_str(&accession(c));
                        s
                    })
                    .collect()
            })
            .collect();
        out_rows.extend(rows.into_iter().flatten());
    }

    let mut w = std::io::BufWriter::new(
        std::fs::File::create(out).map_err(|e| format!("{}: {e}", out.display()))?,
    );
    let io = |e: std::io::Error| format!("{}: {e}", out.display());
    let mut header: Vec<String> = ["SpecId", "Label", "ScanNr", "ExpMass"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    header.extend(charges.clone().map(|z| format!("charge{z}")));
    header.extend(COLUMNS.iter().map(|s| s.to_string()));
    header.extend(["Peptide".to_string(), "Proteins".to_string()]);
    writeln!(w, "{}", header.join("\t")).map_err(io)?;
    for r in &out_rows {
        writeln!(w, "{r}").map_err(io)?;
    }
    w.flush().map_err(io)?;
    eprintln!(
        "DIA Phase B: {} library precursors, {} peak-group rows -> {} [PHASE dia_pc: {:.2}s]",
        n_prec_total,
        out_rows.len(),
        out.display(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn rt_min(spectra: &[Spectrum], s: u32) -> f64 {
    spectra[s as usize].rt_seconds.unwrap_or(0.0) / 60.0
}

/// The run's isolation windows (by recorded target m/z), sorted by lower edge.
fn isolation_windows(spectra: &[Spectrum], range: std::ops::Range<usize>) -> Vec<Window> {
    let mut by: HashMap<i64, Window> = HashMap::new();
    for i in range {
        let s = &spectra[i];
        let (lo, hi) = search::precursor_matching::isolation_window_mz(s, 1.0);
        let t = s.isolation_target_mz.unwrap_or(s.precursor_mz);
        by.entry((t * 100.0).round() as i64)
            .or_insert_with(|| Window {
                lo,
                hi,
                scans: Vec::new(),
            })
            .scans
            .push(i as u32);
    }
    let mut v: Vec<Window> = by.into_values().collect();
    for w in &mut v {
        w.scans.sort_by(|&a, &b| {
            spectra[a as usize]
                .rt_seconds
                .unwrap_or(0.0)
                .total_cmp(&spectra[b as usize].rt_seconds.unwrap_or(0.0))
        });
    }
    v.sort_by(|a, b| a.lo.total_cmp(&b.lo));
    v
}

/// The window whose range holds `mz` and whose centre is nearest it (windows
/// overlap at their edges).
fn window_of(windows: &[Window], mz: f64) -> Option<usize> {
    let end = windows.partition_point(|w| w.lo <= mz);
    windows[..end]
        .iter()
        .enumerate()
        .rev()
        .take(4)
        .filter(|(_, w)| mz <= w.hi)
        .min_by(|a, b| {
            let c = |w: &Window| ((w.lo + w.hi) / 2.0 - mz).abs();
            c(a.1).total_cmp(&c(b.1))
        })
        .map(|(i, _)| i)
}

/// The precursor's assay: (m/z, predicted relative intensity) of its most
/// intense predicted b/y ions of length >= [`MIN_ION_LEN`], intensity-sorted.
fn assay(mass: &Mass, charge: u8, pred: &[f32], ions: &[(bool, u8)]) -> Vec<(f64, f64)> {
    let n = mass.res.len();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut out: Vec<(f64, f64)> = Vec::new();
    let mut prefix = mass.nterm;
    for i in 0..n.saturating_sub(1) {
        prefix += mass.res[i];
        for (k, &(is_b, z)) in ions.iter().enumerate() {
            let len = if is_b { i + 1 } else { n - i - 1 };
            if z > charge || len < MIN_ION_LEN {
                continue;
            }
            let v = pred.get(i * ions.len() + k).copied().unwrap_or(0.0) as f64;
            if v <= 0.0 {
                continue;
            }
            let neutral = if is_b { prefix } else { total - prefix + H2O };
            let mz = (neutral + z as f64 * PROTON) / z as f64;
            if mz < 150.0 {
                continue;
            }
            out.push((mz, v));
        }
    }
    out.sort_by(|a, b| b.1.total_cmp(&a.1));
    out.truncate(N_FRAGMENTS);
    let top = out.first().map_or(1.0, |f| f.1);
    for f in &mut out {
        f.1 /= top;
    }
    out
}

/// The most intense peak within `ppm` of `mz` in m/z-sorted `peaks`.
fn best_peak(peaks: &[(f64, f32)], mz: f64, ppm: f64) -> Option<(f64, f32)> {
    let tol = mz * ppm * 1e-6;
    let lo = peaks.partition_point(|p| p.0 < mz - tol);
    peaks[lo..]
        .iter()
        .take_while(|p| p.0 <= mz + tol)
        .copied()
        .max_by(|a, b| a.1.total_cmp(&b.1))
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if a.len() < 3 {
        return 0.0;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
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

/// Normalised spectral contrast angle of sqrt intensities, in [0, 1].
fn angle(pred: &[f64], obs: &[f64]) -> f64 {
    let (mut d, mut a, mut b) = (0.0, 0.0, 0.0);
    for (p, o) in pred.iter().zip(obs) {
        let (p, o) = (p.max(0.0).sqrt(), o.max(0.0).sqrt());
        d += p * o;
        a += p * p;
        b += o * o;
    }
    if a <= 0.0 || b <= 0.0 {
        return 0.0;
    }
    let c = (d / (a * b).sqrt()).clamp(0.0, 1.0);
    1.0 - 2.0 * c.acos() / std::f64::consts::PI
}

/// Three-point moving average (peaks are ~2 cycles wide at half height).
fn smooth(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    (0..n)
        .map(|i| {
            let (a, b) = (i.saturating_sub(1), (i + 2).min(n));
            x[a..b].iter().sum::<f64>() / (b - a) as f64
        })
        .collect()
}

/// Signed apex errors (ppm, uncalibrated) of the precursor's best peak group by
/// fragment co-elution: each matched fragment's, and the MS1 monoisotope's.
fn apex_errors(
    scans: &[&Spectrum],
    assay: &[(f64, f64)],
    mz: f64,
    ms1: Option<&Ms1RunIndex>,
) -> (Vec<f64>, Option<f64>) {
    let n = scans.len();
    let summed: Vec<f64> = (0..n)
        .map(|s| {
            assay
                .iter()
                .map(|&(fmz, rel)| {
                    best_peak(&scans[s].peaks, fmz, FRAG_PPM)
                        .map_or(0.0, |p| p.1 as f64 * rel.sqrt())
                })
                .sum()
        })
        .collect();
    let sm = smooth(&summed);
    let Some(a) = (0..n).max_by(|&x, &y| sm[x].total_cmp(&sm[y])) else {
        return (Vec::new(), None);
    };
    let frag = assay
        .iter()
        .filter_map(|&(fmz, _)| {
            best_peak(&scans[a].peaks, fmz, FRAG_PPM).map(|p| (p.0 - fmz) / fmz * 1e6)
        })
        .collect();
    let m = ms1.and_then(|ms1| {
        let i = ms1.nearest_scan(scans[a].rt_seconds?)?;
        best_peak(&ms1.scan(i).peaks, mz, MS1_PPM).map(|p| (p.0 - mz) / mz * 1e6)
    });
    (frag, m)
}

/// The peak groups of one precursor on its window's `scans` (RT order), as
/// [`COLUMNS`] values, best candidate apex first.
#[allow(clippy::too_many_arguments)]
fn score_precursor(
    scans: &[&Spectrum],
    assay: &[(f64, f64)],
    mz: f64,
    charge: u8,
    mass: &Mass,
    rt_pred: f64,
    rt_width: f64,
    ms1: Option<&Ms1RunIndex>,
    cal: &MassCal,
) -> Vec<[f64; COLUMNS.len()]> {
    let n = scans.len();
    let nf = assay.len();
    // traces[f][s], and the mass error of each point.
    let mut traces = vec![vec![0.0f64; n]; nf];
    let mut errs = vec![vec![f64::NAN; n]; nf];
    for (s, sp) in scans.iter().enumerate() {
        for (f, &(fmz, _)) in assay.iter().enumerate() {
            let c = fmz * (1.0 + cal.frag_offset * 1e-6);
            if let Some((m, i)) = best_peak(&sp.peaks, c, cal.frag_tol) {
                traces[f][s] = i as f64;
                errs[f][s] = (m - fmz) / fmz * 1e6 - cal.frag_offset;
            }
        }
    }
    let sm: Vec<Vec<f64>> = traces.iter().map(|t| smooth(t)).collect();
    let summed: Vec<f64> = (0..n)
        .map(|s| (0..nf).map(|f| sm[f][s] * assay[f].1.sqrt()).sum::<f64>())
        .collect();
    // Candidate apexes: local maxima of the summed trace, highest first.
    let mut apexes: Vec<usize> = (0..n)
        .filter(|&s| {
            summed[s] > 0.0
                && (s == 0 || summed[s] >= summed[s - 1])
                && (s + 1 == n || summed[s] > summed[s + 1])
        })
        .collect();
    apexes.sort_by(|&a, &b| summed[b].total_cmp(&summed[a]));
    // Score the most intense local maxima; keep the best by co-elution and
    // spectral agreement, not by intensity (an interferer is often brighter).
    apexes.truncate(CANDIDATE_APEXES);

    let pred: Vec<f64> = assay.iter().map(|f| f.1).collect();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let theo_iso = averagine_isotope_envelope(total + H2O, 3);
    let mut groups: Vec<[f64; COLUMNS.len()]> = apexes
        .into_iter()
        .enumerate()
        .map(|(rank, a)| {
            let (lo, hi) = (
                a.saturating_sub(PEAK_HALF_SCANS),
                (a + PEAK_HALF_SCANS + 1).min(n),
            );
            let seg_sum: Vec<f64> = (lo..hi).map(|s| summed[s]).collect();
            // Each raw fragment trace against the sum of the others (leave one
            // out, so a dominant fragment does not correlate with itself).
            let raw_sum: Vec<f64> = (lo..hi)
                .map(|s| (0..nf).map(|f| traces[f][s] * assay[f].1.sqrt()).sum())
                .collect();
            let corrs: Vec<f64> = (0..nf)
                .map(|f| {
                    let w = assay[f].1.sqrt();
                    let own: Vec<f64> = traces[f][lo..hi].to_vec();
                    let rest: Vec<f64> = raw_sum.iter().zip(&own).map(|(s, x)| s - w * x).collect();
                    pearson(&own, &rest)
                })
                .collect();
            let mut sorted = corrs.clone();
            sorted.sort_by(|x, y| y.total_cmp(x));
            let corr_mean = corrs.iter().sum::<f64>() / nf as f64;
            let n06 = corrs.iter().filter(|&&c| c >= 0.6).count() as f64;
            let top3 = sorted.iter().take(3).sum::<f64>() / 3.0;
            let apex_obs: Vec<f64> = (0..nf).map(|f| traces[f][a]).collect();
            // Sparse-peak evidence (weak precursors are 1-3 scans wide):
            // fragments whose own maximum in the segment lies within one scan
            // of the apex, and whether the three most intense predicted ones
            // are present at the apex.
            let coincide = (0..nf)
                .filter(|&f| {
                    let seg = &traces[f][lo..hi];
                    let (m, v) = seg.iter().enumerate().fold((0usize, 0.0f64), |b, (i, &x)| {
                        if x > b.1 {
                            (i, x)
                        } else {
                            b
                        }
                    });
                    v > 0.0 && (lo + m).abs_diff(a) <= 1
                })
                .count() as f64;
            let top3_at_apex = (0..nf.min(3)).filter(|&f| apex_obs[f] > 0.0).count() as f64;
            // Matched fragments at the apex against the window's median scan:
            // a chance bump does not beat the background by much.
            let mut per_scan: Vec<usize> = (0..n)
                .map(|s| (0..nf).filter(|&f| traces[f][s] > 0.0).count())
                .collect();
            let at_apex = per_scan[a] as f64;
            per_scan.sort_unstable();
            let bg = per_scan[n / 2] as f64;
            let vs_bg = at_apex - bg;
            // Mean |lag| (scans) at which each fragment's trace best matches the
            // summed trace of the others.
            let mut lags = 0.0f64;
            let mut nl = 0usize;
            for f in 0..nf {
                let mut best = (0i64, f64::MIN);
                for lag in -2i64..=2 {
                    let mut dot = 0.0;
                    for s in lo..hi {
                        let t = s as i64 + lag;
                        if t >= 0 && (t as usize) < n {
                            dot +=
                                traces[f][t as usize] * (summed[s] - sm[f][s] * assay[f].1.sqrt());
                        }
                    }
                    if dot > best.1 {
                        best = (lag, dot);
                    }
                }
                if best.1 > 0.0 {
                    lags += best.0.unsigned_abs() as f64;
                    nl += 1;
                }
            }
            let xshift = if nl > 0 { lags / nl as f64 } else { 2.0 };
            // Predicted-intensity-weighted share of fragments present at the apex.
            let pred_dot = (0..nf)
                .filter(|&f| apex_obs[f] > 0.0)
                .map(|f| assay[f].1)
                .sum::<f64>()
                / assay.iter().map(|f| f.1).sum::<f64>().max(1e-9);
            let area_obs: Vec<f64> = (0..nf).map(|f| traces[f][lo..hi].iter().sum()).collect();
            let matched = apex_obs.iter().filter(|&&v| v > 0.0).count() as f64;
            let ppms: Vec<f64> = (0..nf)
                .filter_map(|f| {
                    let e = errs[f][a];
                    e.is_finite().then_some(e.abs())
                })
                .collect();
            let apex_ppm = if ppms.is_empty() {
                cal.frag_tol
            } else {
                ppms.iter().sum::<f64>() / ppms.len() as f64
            };
            let rt_apex = scans[a].rt_seconds.unwrap_or(0.0) / 60.0;
            let rt_err = ((rt_apex - rt_pred) / rt_width.max(1e-3)).abs();
            // MS1 at the precursor's mono m/z across the segment.
            let (mut ms1_corr, mut found, mut abs_ppm, mut iso_cos, mut ms1_log) =
                (0.0, 0.0, cal.ms1_tol, 0.0, 0.0);
            let mz_c = mz * (1.0 + cal.ms1_offset * 1e-6);
            if let Some(ms1) = ms1 {
                let tr: Vec<f64> = (lo..hi)
                    .map(|s| {
                        ms1.nearest_scan(scans[s].rt_seconds.unwrap_or(0.0))
                            .and_then(|i| best_peak(&ms1.scan(i).peaks, mz_c, cal.ms1_tol))
                            .map_or(0.0, |p| p.1 as f64)
                    })
                    .collect();
                ms1_corr = pearson(&tr, &seg_sum);
                if let Some(i) = ms1.nearest_scan(scans[a].rt_seconds.unwrap_or(0.0)) {
                    let peaks = &ms1.scan(i).peaks;
                    if let Some((m, v)) = best_peak(peaks, mz_c, cal.ms1_tol) {
                        found = 1.0;
                        abs_ppm = ((m - mz) / mz * 1e6 - cal.ms1_offset).abs();
                        ms1_log = (v as f64).max(1.0).log10();
                        let obs: Vec<f64> = (0..3)
                            .map(|k| {
                                best_peak(
                                    peaks,
                                    mz_c + k as f64 * ISOTOPE / charge as f64,
                                    cal.ms1_tol,
                                )
                                .map_or(0.0, |p| p.1 as f64)
                            })
                            .collect();
                        let d: f64 = obs.iter().zip(&theo_iso).map(|(x, y)| x * y).sum();
                        let (na, nb) = (
                            obs.iter().map(|v| v * v).sum::<f64>().sqrt(),
                            theo_iso.iter().map(|v| v * v).sum::<f64>().sqrt(),
                        );
                        iso_cos = if na > 0.0 && nb > 0.0 {
                            d / (na * nb)
                        } else {
                            0.0
                        };
                    }
                }
            }
            [
                corr_mean,
                n06,
                top3,
                angle(&pred, &apex_obs),
                angle(&pred, &area_obs),
                matched,
                matched / nf as f64,
                apex_ppm,
                apex_obs.iter().sum::<f64>().max(1.0).log10(),
                area_obs.iter().sum::<f64>().max(1.0).log10(),
                ms1_corr,
                found,
                abs_ppm,
                iso_cos,
                ms1_log,
                rt_err,
                rank as f64,
                nf as f64,
                coincide,
                top3_at_apex,
                vs_bg,
                xshift,
                pred_dot,
            ]
        })
        .collect();
    // Peak-group choice: co-elution x spectral agreement, with the sparse-peak
    // evidence for weak precursors.
    let key = |g: &[f64; COLUMNS.len()]| (g[0].max(0.0) + 0.1 * g[18]) * g[4] * (0.2 + g[22]);
    groups.sort_by(|a, b| key(b).total_cmp(&key(a)));
    groups.truncate(PEAK_GROUPS);
    for (k, g) in groups.iter_mut().enumerate() {
        g[16] = k as f64;
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A precursor whose predicted fragments co-elute, in predicted ratios, at
    /// the predicted RT scores its peak group high; a second, unrelated bump
    /// of one fragment ranks below it.
    #[test]
    fn coeluting_assay_scores_its_peak_group_first() {
        let (_, mass) = parse_peptide("K.PEPTIDEKLR.A").unwrap();
        let assay: Vec<(f64, f64)> = vec![
            (400.0, 1.0),
            (500.0, 0.6),
            (600.0, 0.4),
            (700.0, 0.3),
            (800.0, 0.2),
        ];
        let shape = [
            0.0, 0.05, 0.2, 0.6, 1.0, 0.7, 0.3, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let spectra: Vec<Spectrum> = shape
            .iter()
            .enumerate()
            .map(|(k, &h)| {
                let mut peaks: Vec<(f64, f32)> = assay
                    .iter()
                    .map(|&(mz, rel)| (mz, (1e5 * h * rel) as f32 + 1.0))
                    .collect();
                // Interference on one fragment only, late in the window.
                if k == 12 {
                    peaks[0].1 = 5e4;
                }
                Spectrum {
                    rt_seconds: Some(60.0 * (10.0 + k as f64 * 0.03)),
                    peaks,
                    ..Default::default()
                }
            })
            .collect();
        let refs: Vec<&Spectrum> = spectra.iter().collect();
        let g = score_precursor(
            &refs,
            &assay,
            500.0,
            2,
            &mass,
            10.12,
            0.05,
            None,
            &MassCal::default(),
        );
        let f = g[0];
        assert_eq!(f[16], 0.0, "best peak group first");
        assert!(f[0] > 0.9, "corr mean {}", f[0]);
        assert!(f[3] > 0.95, "apex angle {}", f[3]);
        assert_eq!(f[5], 5.0);
        assert!(f[15] < 1.0, "rt err {}", f[15]);
        if let Some(second) = g.get(1) {
            assert!(second[0] < f[0] && second[3] < f[3]);
        }
    }

    #[test]
    fn window_of_picks_the_window_holding_the_mz() {
        let w = |lo: f64, hi: f64| Window {
            lo,
            hi,
            scans: Vec::new(),
        };
        let ws = vec![w(399.0, 401.0), w(401.0, 403.0), w(403.0, 405.0)];
        assert_eq!(window_of(&ws, 400.2), Some(0));
        assert_eq!(window_of(&ws, 402.0), Some(1));
        assert_eq!(window_of(&ws, 406.0), None);
        assert_eq!(window_of(&ws, 398.0), None);
    }
}
