use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use rayon::prelude::*;

use model::stats::{median_mut, Even};

use crate::acquisition::Acquisition;
use model::{InstrumentType, Spectrum};
use onnx_infer::predictors::calibration::RtCalibration;
use onnx_infer::predictors::form::{Form, Mass};
use onnx_infer::predictors::spectral::{
    ms2_features, ms2_features_with, FragmentTolerance, MS2_FEATURES,
};
use onnx_infer::predictors::{
    DeepLc, Engine, NativeMs2, NativeRt, PeptDeep, Predictor, DEEPLC_FILE, PEPTDEEP_FILES,
};
use onnx_infer::{BackendKind, BackendOptions, OnnxSession};
use search::tdc::{qvalues, ScoredLabel};

use crate::cli::{Cli, PredictFeature, PredictorBackend};
pub(crate) use crate::dia::pin::read_pin;

static ANDES_RT_WEIGHTS: &[u8] =
    include_bytes!("../../../../../../resources/predictors/andes_rt.artn");

static ANDES_MS2_WEIGHTS: &[u8] =
    include_bytes!("../../../../../../resources/predictors/andes_ms2.ams2");

pub(crate) const ANCHOR_Q: f64 = 0.01;

/// The bundled andes-rt model as the per-form RT predictor of the DIA vote cut: the raw
/// prediction of a peptidoform's own sequence, NaN when the model cannot encode it.
pub(crate) fn andes_rt_form_model(
) -> Result<std::sync::Arc<search::dia::FormRtModel>, String> {
    use onnx_infer::rt::{RtModel, RtPeptide, Scratch};
    thread_local! {
        static RT_SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::default());
    }
    let rt_model =
        RtModel::from_bytes(ANDES_RT_WEIGHTS).map_err(|e| format!("andes-rt weights: {e}"))?;
    Ok(std::sync::Arc::new(move |p: &model::peptide::Peptide| {
        RtPeptide::from_andes(&p.to_string())
            .and_then(|f| RT_SCRATCH.with(|s| rt_model.predict_with(&f, &mut s.borrow_mut())))
            .unwrap_or(f32::NAN)
    }))
}

struct Source {
    pred: Box<dyn Predictor>,
    sessions: Vec<OnnxSession>,
    rt: bool,
    ms2: bool,
}

/// Predicted fragment intensities of the first MS2 source, kept for the DIA pool scores.
/// [`PredictedMs2::get`] gives, for a (peptidoform id of the PIN table, charge) pair, one value
/// per cleavage site and ion of `ions`, site-major.
pub(crate) struct PredictedMs2 {
    pub(crate) ions: Vec<onnx_infer::predictors::spectral::Ion>,
    pair: HashMap<(u32, u8), u32>,
    values: Vec<Vec<f32>>,
    /// The andes-rt source's raw (uncalibrated) prediction of every peptidoform of the PIN
    /// table, by form id (NaN where it has none); empty without that source.
    pub(crate) rt_raw: Vec<f32>,
    /// Each run's andes-rt calibration from raw prediction to minutes (`andesrt_predicted_rt`);
    /// `None` for a run with too few anchors.
    pub(crate) rt_cal: Vec<Option<RtCalibration>>,
}

impl PredictedMs2 {
    /// The prediction of peptidoform `form` (its id in the PIN table) at `charge`.
    pub(crate) fn get(&self, form: u32, charge: u8) -> Option<&[f32]> {
        self.pair
            .get(&(form, charge))
            .map(|&k| self.values[k as usize].as_slice())
    }

    /// Frees the predictions of every (form, charge) pair not in `keep`; [`Self::get`] then
    /// answers `None` for them and is unchanged for the kept pairs.
    pub(crate) fn keep_only(&mut self, keep: impl Iterator<Item = (u32, u8)>) {
        let mut wanted = vec![false; self.values.len()];
        for p in keep {
            if let Some(&k) = self.pair.get(&p) {
                wanted[k as usize] = true;
            }
        }
        self.pair.retain(|_, k| wanted[*k as usize]);
        for (v, w) in self.values.iter_mut().zip(wanted) {
            if !w {
                *v = Vec::new();
            }
        }
    }
}

pub(crate) struct Predictors {
    sources: Vec<Source>,
    tol: FragmentTolerance,
    threads: usize,
    batch: usize,
    dia: bool,
}

fn requested(list: &[PredictFeature]) -> Result<Option<Vec<PredictFeature>>, String> {
    let mut v: Vec<PredictFeature> = list.to_vec();
    v.sort_unstable();
    v.dedup();
    if v.contains(&PredictFeature::None) {
        if v.len() > 1 {
            return Err("--predict none cannot be combined with other entries".into());
        }
        return Ok(None);
    }
    Ok((!v.is_empty()).then_some(v))
}

fn peptdeep_instrument(i: InstrumentType) -> &'static str {
    match i {
        InstrumentType::QExactive => "QE",
        InstrumentType::TimsTOF => "timsTOF",
        InstrumentType::TOF => "SciexTOF",
        _ => "Lumos",
    }
}

impl Predictors {
    pub(crate) fn from_cli(
        cli: &Cli,
        acquisition: &Acquisition,
        model_instrument: InstrumentType,
    ) -> Result<Option<Self>, String> {
        let instrument = acquisition.predictor_instrument(model_instrument);
        let Some(want) = requested(&cli.predict)? else {
            return Ok(None);
        };
        let has = |f: PredictFeature| want.contains(&f);
        let engine = Engine {
            kind: match cli.predictor_backend {
                PredictorBackend::Tract => BackendKind::Tract,
                PredictorBackend::Ort => BackendKind::Ort,
            },
            opts: BackendOptions::default(),
        };
        let third_party = has(PredictFeature::PeptDeepMs2)
            || has(PredictFeature::PeptDeepRt)
            || has(PredictFeature::DeepLcRt);
        let dir: Option<PathBuf> = cli.predictor_models.clone();
        if third_party && dir.is_none() {
            return Err(format!(
                "--predict {} loads exported ONNX models: pass --predictor-models DIR, a directory holding {} (AlphaPeptDeep) and/or {} \
                 (DeepLC), as written by the onnx-infer export scripts",
                want.iter()
                    .filter_map(|f| match f {
                        PredictFeature::PeptDeepMs2 => Some("peptdeep-ms2"),
                        PredictFeature::PeptDeepRt => Some("peptdeep-rt"),
                        PredictFeature::DeepLcRt => Some("deeplc-rt"),
                        _ => None,
                    })
                    .collect::<Vec<_>>()
                    .join(","),
                PEPTDEEP_FILES.join(", "),
                DEEPLC_FILE
            ));
        }
        let mut sources = Vec::new();
        if has(PredictFeature::Rt) {
            sources.push(Source {
                pred: Box::new(NativeRt::from_bytes(ANDES_RT_WEIGHTS).map_err(|e| e.to_string())?),
                sessions: Vec::new(),
                rt: true,
                ms2: false,
            });
        }
        if has(PredictFeature::Ms2) {
            let nce = cli
                .predictor_nce
                .map(f64::from)
                .or_else(|| acquisition.predictor_nce())
                .map(|v| v as f32);
            let activation = acquisition.ms2_activation();
            let inst = cli
                .predictor_instrument
                .clone()
                .unwrap_or_else(|| acquisition.ms2_instrument().to_string());
            log::info!(
                target: "predict",
                "andes-ms2 (activation {activation}, instrument {}, NCE {})",
                if inst.is_empty() { "unknown" } else { &inst },
                nce.map_or("unknown".to_string(), |v| v.to_string())
            );
            sources.push(Source {
                pred: Box::new(
                    NativeMs2::from_bytes(ANDES_MS2_WEIGHTS, nce, activation, &inst)
                        .map_err(|e| e.to_string())?,
                ),
                sessions: Vec::new(),
                rt: false,
                ms2: true,
            });
        }
        if has(PredictFeature::PeptDeepMs2) || has(PredictFeature::PeptDeepRt) {
            let dir = dir.as_deref().expect("checked above");
            let nce = cli
                .predictor_nce
                .or_else(|| acquisition.predictor_nce().map(|v| v as f32))
                .unwrap_or(match instrument {
                    InstrumentType::QExactive => 28.0,
                    _ => 30.0,
                });
            let inst = cli
                .predictor_instrument
                .clone()
                .unwrap_or_else(|| peptdeep_instrument(instrument).to_string());
            let p = PeptDeep::new(
                dir,
                nce,
                &inst,
                engine.clone(),
                has(PredictFeature::PeptDeepRt),
                has(PredictFeature::PeptDeepMs2),
            )
            .map_err(|e| format!("--predict: {e}"))?;
            log::info!(
                target: "predict",
                "AlphaPeptDeep from {} (instrument {inst}, NCE {nce}); modifications: {}",
                dir.display(),
                {
                    let mut m: Vec<&str> = p.known_mods().collect();
                    m.sort_unstable();
                    m.join(", ")
                }
            );
            sources.push(Source {
                pred: Box::new(p),
                sessions: Vec::new(),
                rt: has(PredictFeature::PeptDeepRt),
                ms2: has(PredictFeature::PeptDeepMs2),
            });
        }
        if has(PredictFeature::DeepLcRt) {
            let dir = dir.as_deref().expect("checked above");
            sources.push(Source {
                pred: Box::new(
                    DeepLc::new(dir, engine.clone()).map_err(|e| format!("--predict: {e}"))?,
                ),
                sessions: Vec::new(),
                rt: true,
                ms2: false,
            });
        }
        for s in &mut sources {
            s.sessions = s
                .pred
                .sessions()
                .map_err(|e| format!("--predict {}: {e}", s.pred.prefix()))?;
        }
        let high_res = acquisition.high_res;
        let tol = match (cli.fragment_tol_ppm, cli.fragment_tol_da) {
            (Some(p), _) => FragmentTolerance::Ppm(p),
            (None, Some(d)) => FragmentTolerance::Da(d),
            (None, None) if high_res => FragmentTolerance::Ppm(20.0),
            (None, None) => FragmentTolerance::Da(0.5),
        };
        Ok(Some(Self {
            sources,
            tol,
            threads: cli.threads.max(1),
            batch: cli.predictor_batch.max(1),
            dia: cli.dia,
        }))
    }

    pub(crate) fn tolerance(&self) -> FragmentTolerance {
        self.tol
    }

    pub(crate) fn needs_peaks(&self) -> bool {
        self.sources.iter().any(|s| s.ms2)
    }

    pub(crate) fn describe(&self) -> String {
        let blocks: Vec<String> = self
            .sources
            .iter()
            .map(|s| {
                let kinds: Vec<&str> = [(s.ms2, "MS2"), (s.rt, "RT")]
                    .iter()
                    .filter(|(on, _)| *on)
                    .map(|(_, k)| *k)
                    .collect();
                format!("{} {}", s.pred.prefix(), kinds.join("+"))
            })
            .collect();
        format!("{} (fragment tolerance {:?})", blocks.join(", "), self.tol)
    }

    fn columns(&self) -> Vec<String> {
        let mut cols = Vec::new();
        for s in &self.sources {
            let p = s.pred.prefix();
            if s.ms2 {
                cols.extend(MS2_FEATURES.iter().map(|f| format!("{p}_{f}")));
                if self.dia {
                    cols.extend(
                        crate::dia::EVIDENCE_COLUMNS
                            .iter()
                            .map(|f| format!("{p}_{f}")),
                    );
                }
            }
            if s.rt {
                cols.push(format!("{p}_rt_diff"));
                cols.push(format!("{p}_predicted_rt"));
            }
        }
        for s in &self.sources {
            cols.push(format!("{}_missing", s.pred.prefix()));
        }
        cols
    }

    /// Adds the predicted-feature columns to the PIN; with `keep_ms2`, also returns the first
    /// MS2 source's predictions and the PIN table they index, which describes the rewritten
    /// PIN (same rows, the new columns before `Peptide`).
    pub(crate) fn add_to_pin(
        &self,
        pin: &Path,
        spectra: &[Spectrum],
        run_ends: &[usize],
        keep_ms2: bool,
    ) -> Result<Option<(PredictedMs2, PinTable)>, String> {
        let t0 = Instant::now();
        let c0 = process_cpu_seconds();
        let new_cols = self.columns();
        let mut table = read_pin(pin, spectra, &new_cols)?;
        crate::memory::mark_stage("predictions: PIN read");
        let n = table.rows.len();
        let run_of = |spec: u32| run_ends.partition_point(|&e| e <= spec as usize);
        let n_runs = run_ends.len().max(1);

        let obs: Vec<f64> = table
            .rows
            .iter()
            .map(|r| {
                (r.spec != NONE)
                    .then(|| spectra[r.spec as usize].rt_seconds)
                    .flatten()
                    .filter(|t| t.is_finite())
                    .map_or(f64::NAN, |t| t / 60.0)
            })
            .collect();
        let anchors = anchor_rows(&table.rows, n_runs, &run_of);

        let mut blocks: Vec<(Vec<Vec<f64>>, Vec<bool>)> = Vec::new();
        let mut report: Vec<String> = Vec::new();
        let mut kept: Option<PredictedMs2> = None;
        let mut kept_rt: Option<(Vec<f32>, Vec<Option<RtCalibration>>)> = None;
        for src in &self.sources {
            let pred = src.pred.as_ref();
            let mut miss = vec![false; n];
            let mut cols: Vec<Vec<f64>> = Vec::new();
            if src.ms2 {
                let keep = (keep_ms2 && kept.is_none()).then_some(&mut kept);
                let (ms2_cols, msg) =
                    self.ms2_block(src, &table, spectra, run_ends, &mut miss, keep)?;
                cols.extend(ms2_cols);
                report.push(msg);
            }
            if src.rt {
                let raw = self.predict_rt(src, &table.forms)?;
                let mut diff = vec![f64::NAN; n];
                let mut cal_rt = vec![f64::NAN; n];
                let mut per_run = Vec::new();
                let mut run_cals: Vec<Option<RtCalibration>> = Vec::with_capacity(anchors.len());
                for (run, rows) in anchors.iter().enumerate() {
                    let anchor_pts: Vec<(usize, f64, f64)> = rows
                        .iter()
                        .filter_map(|&i| {
                            let f = table.rows[i].form;
                            (f != NONE).then(|| (i, raw[f as usize] as f64, obs[i]))
                        })
                        .filter(|(_, p, o)| p.is_finite() && o.is_finite())
                        .collect();
                    let pts: Vec<(f64, f64)> =
                        anchor_pts.iter().map(|&(_, p, o)| (p, o)).collect();
                    let cal = RtCalibration::fit(&pts);
                    let held_out = cross_fit_calibrations(&table.rows, &anchor_pts);
                    let resid = cal.as_ref().map(|c| {
                        let mut r: Vec<f64> =
                            pts.iter().map(|(p, o)| (o - c.apply(*p)).abs()).collect();
                        r.sort_by(f64::total_cmp);
                        r.get(r.len() / 2).copied().unwrap_or(f64::NAN)
                    });
                    per_run.push(format!(
                        "file {run}: {} anchors{}",
                        pts.len(),
                        resid.map_or(" (no calibration)".to_string(), |m| format!(
                            ", median |residual| {m:.2} min"
                        ))
                    ));
                    run_cals.push(cal.clone());
                    if let Some(c) = cal {
                        let fold_of: HashMap<usize, usize> = anchor_pts
                            .iter()
                            .map(|&(i, _, _)| (i, rt_fold(&table.rows[i])))
                            .collect();
                        cal_rt
                            .par_iter_mut()
                            .zip(diff.par_iter_mut())
                            .enumerate()
                            .for_each(|(i, (cal_i, diff_i))| {
                                let r = &table.rows[i];
                                if r.spec == NONE || run_of(r.spec) != run || r.form == NONE {
                                    return;
                                }
                                let p = raw[r.form as usize] as f64;
                                if p.is_finite() && obs[i].is_finite() {
                                    let own = fold_of.get(&i).and_then(|&k| held_out[k].as_ref());
                                    *cal_i = own.unwrap_or(&c).apply(p);
                                    *diff_i = (obs[i] - *cal_i).abs();
                                }
                            });
                    }
                }
                for (m, d) in miss.iter_mut().zip(&diff) {
                    *m |= !d.is_finite();
                }
                report.push(format!("{} RT: {}", pred.prefix(), per_run.join("; ")));
                if keep_ms2 && kept_rt.is_none() && pred.prefix() == ANDES_RT_PREFIX {
                    kept_rt = Some((raw, run_cals));
                }
                cols.push(diff);
                cols.push(cal_rt);
            }
            let n_miss = miss.iter().filter(|m| **m).count();
            report.push(format!("{} missing rows {n_miss} of {n}", pred.prefix()));
            if n > 0 && n_miss * 5 > n {
                log::warn!(
                    target: "predict",
                    "{} could not score {n_miss} of {n} PIN rows ({:.0}%); \
                     a modification the model does not encode (for AlphaPeptDeep: one missing \
                     from peptdeep_mods.json) or missing retention times are the usual causes",
                    pred.prefix(),
                    100.0 * n_miss as f64 / n as f64
                );
            }
            for group in cols.chunks_mut(IMPUTE_COLUMNS_AT_ONCE) {
                group
                    .par_iter_mut()
                    .for_each(|col| impute_median(col, &miss));
            }
            blocks.push((cols, miss));
        }

        for (src, (_, flags)) in self.sources.iter().zip(&blocks) {
            if n > 0 && flags.iter().all(|&m| m) {
                log::warn!(
                    target: "predict",
                    "{} scored no PIN row; its columns are constant",
                    src.pred.prefix()
                );
            }
        }
        let mut values: Vec<Vec<f64>> = Vec::new();
        let mut flag_values: Vec<Vec<bool>> = Vec::new();
        for (cols, flags) in blocks {
            values.extend(cols);
            flag_values.push(flags);
        }
        crate::memory::mark_stage("predictions: columns computed");

        write_pin_with_flags(pin, &table, &new_cols, &values, &flag_values)?;
        drop((values, flag_values));
        table.peptide_col += new_cols.len();
        for line in &report {
            log::info!(target: "predict", "{line}");
        }
        log::info!(
            target: "predict",
            "added {} columns to {} PIN rows ({} peptidoforms) in {:.2} s wall, {:.1} s CPU",
            new_cols.len(),
            n,
            table.forms.len(),
            t0.elapsed().as_secs_f64(),
            process_cpu_seconds() - c0
        );
        if let (Some(k), Some((raw, cals))) = (kept.as_mut(), kept_rt) {
            k.rt_raw = raw;
            k.rt_cal = cals;
        }
        Ok(kept.map(|k| (k, table)))
    }

    fn predict_rt(&self, src: &Source, forms: &[(Form, Mass)]) -> Result<Vec<f32>, String> {
        let pred = src.pred.as_ref();
        let ix: Vec<usize> = (0..forms.len())
            .filter(|&i| pred.supports_rt(&forms[i].0))
            .collect();
        let jobs = batches(ix.iter().map(|&i| pred.batch_key(&forms[i].0)), self.batch);
        let out = Mutex::new(vec![f32::NAN; forms.len()]);
        let t = Instant::now();
        run_parallel(&src.sessions, self.threads, &jobs, |s, job| {
            let fs: Vec<&Form> = job.iter().map(|&j| &forms[ix[j]].0).collect();
            let v = pred
                .rt(s, &fs)
                .map_err(|e| format!("{} RT: {e}", pred.prefix()))?;
            let mut g = out.lock().expect("RT output lock");
            for (&j, x) in job.iter().zip(v) {
                g[ix[j]] = if x.is_finite() { x } else { f32::NAN };
            }
            Ok(())
        })?;
        log::info!(
            target: "predict",
            "{} RT for {} of {} peptidoforms in {:.1}s",
            pred.prefix(),
            ix.len(),
            forms.len(),
            t.elapsed().as_secs_f64()
        );
        Ok(out.into_inner().expect("RT output lock"))
    }

    fn ms2_block(
        &self,
        src: &Source,
        table: &PinTable,
        spectra: &[Spectrum],
        run_ends: &[usize],
        miss: &mut [bool],
        keep: Option<&mut Option<PredictedMs2>>,
    ) -> Result<(Vec<Vec<f64>>, String), String> {
        let pred = src.pred.as_ref();
        let forms = &table.forms;
        let mut fz_ix: HashMap<(u32, u8), u32> = HashMap::new();
        let mut fz: Vec<(u32, u8)> = Vec::new();
        let mut no_peaks = 0usize;
        for r in &table.rows {
            if r.form == NONE || r.spec == NONE || r.charge == 0 {
                continue;
            }
            if spectra[r.spec as usize].peaks.is_empty() {
                no_peaks += 1;
                continue;
            }
            if pred.supports_ms2(&forms[r.form as usize].0, r.charge) {
                fz_ix.entry((r.form, r.charge)).or_insert_with(|| {
                    fz.push((r.form, r.charge));
                    (fz.len() - 1) as u32
                });
            }
        }
        let mut row_start = vec![0u32; fz.len() + 1];
        let fz_of_row = prediction_of_rows(&table.rows, spectra, &fz_ix);
        let fz_ix = keep.is_some().then_some(fz_ix);
        for &k in &fz_of_row {
            if k != NONE {
                row_start[k as usize + 1] += 1;
            }
        }
        for k in 0..fz.len() {
            row_start[k + 1] += row_start[k];
        }
        let mut fill = row_start.clone();
        let mut rows_of_fz = vec![0u32; row_start[fz.len()] as usize];
        for (i, &k) in fz_of_row.iter().enumerate() {
            if k != NONE {
                rows_of_fz[fill[k as usize] as usize] = i as u32;
                fill[k as usize] += 1;
            }
        }
        drop(fill);
        for (m, &k) in miss.iter_mut().zip(&fz_of_row) {
            *m |= k == NONE;
        }
        drop(fz_of_row);
        let jobs = batches(
            fz.iter()
                .map(|(f, _)| pred.batch_key(&forms[*f as usize].0)),
            self.batch,
        );
        let tol = self.tol;
        let dia = self.dia.then(|| crate::dia::windows(spectra, run_ends));
        let width = MS2_FEATURES.len()
            + if dia.is_some() {
                crate::dia::EVIDENCE_COLUMNS.len()
            } else {
                0
            };
        let row_features = |prediction: &[f32], r: &PinRow| -> Vec<f64> {
            let spec = &spectra[r.spec as usize];
            let (_, mass) = &forms[r.form as usize];
            let Some((slot, windows)) = &dia else {
                return ms2_features(prediction, pred.ions(), mass, r.charge, &spec.peaks, tol)
                    .to_vec();
            };
            let (scans, center) = crate::dia::neighbours(spectra, slot, windows, r.spec as usize);
            let (frags, traces) =
                crate::dia::fragment_evidence(prediction, pred.ions(), mass, r.charge, &scans, tol);
            let ns = scans.len();
            let mut v = ms2_features_with(prediction, pred.ions(), mass, r.charge, |k, _| {
                traces[k * ns + center]
            })
            .to_vec();
            v.extend(crate::dia::evidence_of(&frags, &traces, ns, center));
            v
        };
        let out = Mutex::new(vec![vec![f64::NAN; table.rows.len()]; width]);
        let kept_pred: Option<Mutex<Vec<Vec<f32>>>> =
            keep.is_some().then(|| Mutex::new(vec![Vec::new(); fz.len()]));
        let t = Instant::now();
        run_parallel(&src.sessions, self.threads, &jobs, |s, job| {
            let items: Vec<(&Form, u8)> = job
                .iter()
                .map(|&i| (&forms[fz[i].0 as usize].0, fz[i].1))
                .collect();
            let v = pred
                .ms2(s, &items)
                .map_err(|e| format!("{} MS2: {e}", pred.prefix()))?;
            let mut done: Vec<(u32, Vec<f64>)> = Vec::new();
            for (&i, x) in job.iter().zip(&v) {
                for &row in &rows_of_fz[row_start[i] as usize..row_start[i + 1] as usize] {
                    done.push((row, row_features(x, &table.rows[row as usize])));
                }
            }
            let mut g = out.lock().expect("MS2 output lock");
            for (row, f) in done {
                for (c, x) in g.iter_mut().zip(f) {
                    c[row as usize] = x;
                }
            }
            drop(g);
            if let Some(kp) = &kept_pred {
                let mut g = kp.lock().expect("MS2 prediction lock");
                for (&i, x) in job.iter().zip(v) {
                    g[i] = x;
                }
            }
            Ok(())
        })?;
        let cols = out.into_inner().expect("MS2 output lock");
        let predict_s = t.elapsed().as_secs_f64();
        if no_peaks > 0 {
            log::warn!(
                target: "predict",
                "{no_peaks} PIN rows have a spectrum without peaks; \
                 they are marked missing"
            );
        }
        let msg = format!(
            "{} MS2 for {} (peptidoform, charge) in {} batches, {predict_s:.1}s",
            pred.prefix(),
            fz.len(),
            jobs.len()
        );
        if let (Some(slot), Some(kp), Some(pair)) = (keep, kept_pred, fz_ix) {
            *slot = Some(PredictedMs2 {
                ions: pred.ions().to_vec(),
                pair,
                values: kp.into_inner().expect("MS2 prediction lock"),
                rt_raw: Vec::new(),
                rt_cal: Vec::new(),
            });
        }
        Ok((cols, msg))
    }
}

pub(crate) const NONE: u32 = u32::MAX;
/// Feature columns imputed at once; each holds a copy of its finite values while it runs.
const IMPUTE_COLUMNS_AT_ONCE: usize = 4;
/// Column prefix of the andes-rt source, whose calibrated prediction is `andesrt_predicted_rt`.
const ANDES_RT_PREFIX: &str = "andesrt";

/// The index of each row's (peptidoform, charge) prediction, or `NONE` for a row that gets no
/// MS2 features: one whose spectrum was not located or has no peaks, or whose pair has no
/// prediction. Another row with the same pair does not lend such a row its prediction.
fn prediction_of_rows(
    rows: &[PinRow],
    spectra: &[Spectrum],
    fz_ix: &HashMap<(u32, u8), u32>,
) -> Vec<u32> {
    rows.iter()
        .map(|r| {
            if r.spec == NONE || spectra[r.spec as usize].peaks.is_empty() {
                return NONE;
            }
            fz_ix.get(&(r.form, r.charge)).copied().unwrap_or(NONE)
        })
        .collect()
}

pub(crate) struct PinRow {
    pub(crate) spec: u32,
    pub(crate) form: u32,
    pub(crate) charge: u8,
    pub(crate) decoy: bool,
    pub(crate) score: f64,
}

pub(crate) struct PinTable {
    pub(crate) rows: Vec<PinRow>,
    pub(crate) forms: Vec<(Form, Mass)>,
    pub(super) peptide_col: usize,
}

pub(super) fn spectrum_bases(spectra: &[Spectrum]) -> HashMap<String, u32> {
    let mut m: HashMap<String, u32> = HashMap::with_capacity(spectra.len());
    for (i, s) in spectra.iter().enumerate() {
        let scan = s.scan.unwrap_or(0);
        let base = if s.title.is_empty() {
            format!("scan={scan}_{scan}")
        } else {
            format!("{}_{scan}", s.title)
        };
        m.entry(base).and_modify(|v| *v = NONE).or_insert(i as u32);
    }
    m
}

pub(super) fn locate(
    bases: &HashMap<String, u32>,
    spectra: &[Spectrum],
    spec_id: &str,
    scan: i32,
) -> u32 {
    let mut s = spec_id;
    for _ in 0..2 {
        let Some((head, _)) = s.rsplit_once('_') else {
            return NONE;
        };
        if let Some(&i) = bases.get(head) {
            if i != NONE && spectra[i as usize].scan.unwrap_or(0) == scan {
                return i;
            }
        }
        s = head;
    }
    NONE
}

pub(super) fn row_charge(
    hot: Option<u8>,
    edges: (u8, u8),
    exp_mass: Option<f64>,
    precursor_mz: Option<f64>,
) -> u8 {
    let Some(z) = hot else {
        return 0;
    };
    if z != edges.0 && z != edges.1 {
        return z;
    }
    let (Some(m), Some(mz)) = (exp_mass, precursor_mz) else {
        return z;
    };
    let per_charge = mz - model::PROTON;
    if !(m.is_finite() && per_charge > 0.0) {
        return z;
    }
    let whole = (m / per_charge).round();
    if !(1.0..=f64::from(u8::MAX)).contains(&whole) || (m - whole * per_charge).abs() > 0.01 {
        return z;
    }
    let derived = whole as u8;
    if (z == edges.1 && derived > z) || (z == edges.0 && derived < z) {
        derived
    } else {
        z
    }
}

/// Folds of the RT calibration's cross-fit; a peptidoform's anchors share a fold.
const RT_FOLDS: usize = 3;

fn rt_fold(row: &PinRow) -> usize {
    ((u64::from(row.form).wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 40) % RT_FOLDS as u64) as usize
}

/// The run's RT calibration fitted without each fold's anchors: an anchor's `cal_rt` comes
/// from the calibration that did not see it, so the rescorer's RT error of an anchor is not
/// shrunk by its own fit. `None` for a fold whose other anchors are too few to fit; that fold
/// takes the calibration of all anchors.
fn cross_fit_calibrations(
    rows: &[PinRow],
    anchors: &[(usize, f64, f64)],
) -> Vec<Option<RtCalibration>> {
    (0..RT_FOLDS)
        .map(|k| {
            let train: Vec<(f64, f64)> = anchors
                .iter()
                .filter(|&&(i, _, _)| rt_fold(&rows[i]) != k)
                .map(|&(_, p, o)| (p, o))
                .collect();
            RtCalibration::fit(&train)
        })
        .collect()
}

/// Fails when the PIN's header already has one of `new_cols` (read from the header line only).
pub(crate) fn check_new_columns(pin: &Path, new_cols: &[String]) -> Result<(), String> {
    let f = File::open(pin).map_err(|e| format!("{}: {e}", pin.display()))?;
    let header = BufReader::new(f)
        .lines()
        .next()
        .ok_or_else(|| format!("{}: empty PIN", pin.display()))?
        .map_err(|e| format!("{}: {e}", pin.display()))?;
    let h: Vec<&str> = header.split('\t').collect();
    match new_cols.iter().find(|c| h.contains(&c.as_str())) {
        Some(c) => Err(format!("{}: PIN already has a {c} column", pin.display())),
        None => Ok(()),
    }
}

pub(crate) fn anchor_rows(
    rows: &[PinRow],
    n_runs: usize,
    run_of: &dyn Fn(u32) -> usize,
) -> Vec<Vec<usize>> {
    let mut best: HashMap<u32, usize> = HashMap::new();
    for (i, r) in rows.iter().enumerate() {
        if r.spec == NONE {
            continue;
        }
        best.entry(r.spec)
            .and_modify(|b| {
                if rows[*b].score < r.score {
                    *b = i;
                }
            })
            .or_insert(i);
    }
    let mut per_run: Vec<Vec<usize>> = vec![Vec::new(); n_runs];
    let mut best: Vec<(u32, usize)> = best.into_iter().collect();
    best.sort_unstable();
    for (spec, i) in best {
        per_run[run_of(spec).min(n_runs - 1)].push(i);
    }
    per_run
        .into_iter()
        .map(|ix| {
            let items: Vec<ScoredLabel> = ix
                .iter()
                .map(|&i| ScoredLabel {
                    score: rows[i].score as f32,
                    is_decoy: rows[i].decoy,
                })
                .collect();
            let q = qvalues(&items);
            ix.into_iter()
                .zip(q)
                .filter(|(i, q)| !rows[*i].decoy && *q <= ANCHOR_Q)
                .map(|(i, _)| i)
                .collect()
        })
        .collect()
}

fn impute_median(col: &mut [f64], miss: &[bool]) {
    let mut ok: Vec<f64> = col
        .iter()
        .zip(miss)
        .filter(|(v, m)| v.is_finite() && !**m)
        .map(|(v, _)| *v)
        .collect();
    let med = median_mut(&mut ok, Even::Mean).unwrap_or(0.0);
    for (v, m) in col.iter_mut().zip(miss) {
        if !v.is_finite() || *m {
            *v = med;
        }
    }
}

pub(crate) fn write_pin(
    pin: &Path,
    table: &PinTable,
    cols: &[String],
    values: &[Vec<f64>],
) -> Result<(), String> {
    write_pin_with_flags(pin, table, cols, values, &[])
}

/// [`write_pin`] with 0/1 flag columns after the `values` columns; `cols` names both, values
/// first.
pub(crate) fn write_pin_with_flags(
    pin: &Path,
    table: &PinTable,
    cols: &[String],
    values: &[Vec<f64>],
    flags: &[Vec<bool>],
) -> Result<(), String> {
    let tmp = pin.with_file_name(format!(
        "{}.predict.tmp",
        pin.file_name()
            .map_or_else(|| "pin".into(), |n| n.to_string_lossy())
    ));
    let io = |e: std::io::Error| format!("{}: {e}", tmp.display());
    let w = std::cell::RefCell::new(BufWriter::new(File::create(&tmp).map_err(io)?));
    let k = table.peptide_col;
    let n = table.rows.len();
    let mut rows = 0usize;
    crate::dia::blocks::read_blocks(
        pin,
        |header| {
            let h: Vec<&str> = header.split('\t').collect();
            let k = k.min(h.len());
            let (a, b) = (h[..k].join("\t"), h[k..].join("\t"));
            writeln!(w.borrow_mut(), "{a}\t{}\t{b}", cols.join("\t")).map_err(io)
        },
        |first, lines| {
            if first + lines.len() > n {
                return Err(format!(
                    "{}: changed while adding predicted features ({} rows, expected {n})",
                    pin.display(),
                    first + lines.len()
                ));
            }
            let out: Vec<Result<String, String>> = lines
                .par_chunks(crate::dia::blocks::LINES_PER_TASK)
                .enumerate()
                .map(|(c, chunk)| {
                    let mut buf = String::new();
                    for (j, line) in chunk.iter().enumerate() {
                        let row = first + c * crate::dia::blocks::LINES_PER_TASK + j;
                        let (head, tail) = split_before_field(line, k).ok_or_else(|| {
                            format!(
                                "{}: changed while adding predicted features (row {} has at \
                                 most {k} columns)",
                                pin.display(),
                                row + 1
                            )
                        })?;
                        buf.push_str(head);
                        for v in values {
                            buf.push('\t');
                            push_value(&mut buf, v[row]);
                        }
                        for f in flags {
                            buf.push('\t');
                            push_value(&mut buf, f64::from(u8::from(f[row])));
                        }
                        buf.push('\t');
                        buf.push_str(tail);
                        buf.push('\n');
                    }
                    Ok(buf)
                })
                .collect();
            let mut w = w.borrow_mut();
            for part in out {
                w.write_all(part?.as_bytes()).map_err(io)?;
            }
            rows = first + lines.len();
            Ok(())
        },
    )?;
    if rows != n {
        return Err(format!(
            "{}: changed while adding predicted features ({} rows, expected {})",
            pin.display(),
            rows,
            n
        ));
    }
    let mut w = w.into_inner();
    w.flush().map_err(io)?;
    drop(w);
    std::fs::rename(&tmp, pin).map_err(|e| format!("{} -> {}: {e}", tmp.display(), pin.display()))
}

fn split_before_field(line: &str, k: usize) -> Option<(&str, &str)> {
    if k == 0 {
        return Some(("", line));
    }
    let (i, _) = line.match_indices('\t').nth(k - 1)?;
    Some((&line[..i], &line[i + 1..]))
}

fn push_value(buf: &mut String, v: f64) {
    use std::fmt::Write as _;
    if v.fract() == 0.0 && v.abs() < 1e15 {
        let _ = write!(buf, "{}", v as i64);
    } else if v != 0.0 && (v.abs() < 1e-4 || v.abs() >= 1e6) {
        let _ = write!(buf, "{v:.5e}");
    } else {
        let digits = (5 - v.abs().log10().floor() as i32).clamp(0, 10) as usize;
        let _ = write!(buf, "{v:.digits$}");
    }
}

fn batches(keys: impl Iterator<Item = usize>, batch: usize) -> Vec<Vec<usize>> {
    let mut by: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, k) in keys.enumerate() {
        by.entry(k).or_default().push(i);
    }
    let mut groups: Vec<_> = by.into_iter().collect();
    groups.sort_by_key(|(k, _)| *k);
    groups
        .into_iter()
        .flat_map(|(_, ix)| ix.chunks(batch).map(<[usize]>::to_vec).collect::<Vec<_>>())
        .collect()
}

fn run_parallel<J: Sync, F>(
    master: &[OnnxSession],
    threads: usize,
    jobs: &[J],
    work: F,
) -> Result<(), String>
where
    F: Fn(&mut [OnnxSession], &J) -> Result<(), String> + Sync,
{
    if jobs.is_empty() {
        return Ok(());
    }
    let next = AtomicUsize::new(0);
    let failed: Mutex<Option<String>> = Mutex::new(None);
    let workers = threads.min(jobs.len()).max(1);
    let handles: Vec<Vec<OnnxSession>> = (0..workers)
        .map(|_| {
            master
                .iter()
                .map(OnnxSession::share)
                .collect::<Result<_, _>>()
        })
        .collect::<Result<_, _>>()
        .map_err(|e| format!("sharing a predictor session: {e}"))?;
    std::thread::scope(|sc| {
        for mut sessions in handles {
            let (next, work, failed) = (&next, &work, &failed);
            sc.spawn(move || loop {
                let j = next.fetch_add(1, Ordering::Relaxed);
                if j >= jobs.len() {
                    break;
                }
                if let Err(e) = work(&mut sessions, &jobs[j]) {
                    failed.lock().expect("error lock").get_or_insert(e);
                    next.store(jobs.len(), Ordering::Relaxed);
                    break;
                }
            });
        }
    });
    match failed.into_inner().expect("error lock") {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

fn process_cpu_seconds() -> f64 {
    std::fs::read_to_string("/proc/self/stat")
        .ok()
        .and_then(|s| {
            let rest = s.rsplit_once(')')?.1;
            let f: Vec<&str> = rest.split_whitespace().collect();
            let t = f.get(11)?.parse::<f64>().ok()? + f.get(12)?.parse::<f64>().ok()?;
            Some(t / 100.0)
        })
        .unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spectrum(title: &str, scan: i32, rt_min: f64) -> Spectrum {
        Spectrum {
            title: title.to_string(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: Some(rt_min * 60.0),
            scan: Some(scan),
            peaks: vec![(100.0, 1.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        }
    }

    #[test]
    fn requested_dedups_and_rejects_none_with_others() {
        use PredictFeature::*;
        assert_eq!(requested(&[None]).unwrap(), Option::None);
        assert_eq!(requested(&[]).unwrap(), Option::None);
        assert_eq!(requested(&[Ms2, Rt, Rt]).unwrap(), Some(vec![Rt, Ms2]));
        assert!(requested(&[None, Rt]).is_err());
    }

    #[test]
    fn rows_locate_their_spectrum_by_spec_id_and_scan() {
        let spectra = vec![
            spectrum("run.10.10.2", 10, 1.0),
            spectrum("", 11, 2.0),
            spectrum("a_b", 12, 3.0),
        ];
        let b = spectrum_bases(&spectra);
        assert_eq!(locate(&b, &spectra, "run.10.10.2_10_1", 10), 0);
        assert_eq!(locate(&b, &spectra, "run.10.10.2_10_2_1", 10), 0);
        assert_eq!(locate(&b, &spectra, "scan=11_11_1", 11), 1);
        assert_eq!(locate(&b, &spectra, "a_b_12_1_0", 12), 2);
        assert_eq!(locate(&b, &spectra, "run.10.10.2_10_1", 9), NONE);
        assert_eq!(locate(&b, &spectra, "other_5_1", 5), NONE);
    }

    #[test]
    fn unlocated_and_peakless_rows_get_no_prediction_of_a_located_row() {
        let mut empty = spectrum("run.11.11.2", 11, 2.0);
        empty.peaks.clear();
        let spectra = vec![spectrum("run.10.10.2", 10, 1.0), empty];
        let row = |spec: u32, form: u32| PinRow {
            spec,
            form,
            charge: 2,
            decoy: false,
            score: 0.0,
        };
        let rows = [row(0, 7), row(NONE, 7), row(1, 7), row(0, 8)];
        let fz_ix = HashMap::from([((7u32, 2u8), 0u32)]);
        assert_eq!(
            prediction_of_rows(&rows, &spectra, &fz_ix),
            [0, NONE, NONE, NONE]
        );
    }

    #[test]
    fn row_charge_recovers_a_charge_beyond_the_columns() {
        let mz = 500.0;
        let exp = |z: f64| Some(z * (mz - model::PROTON));
        assert_eq!(row_charge(Some(3), (2, 4), exp(3.0), Some(mz)), 3);
        assert_eq!(row_charge(Some(4), (2, 4), exp(6.0), Some(mz)), 6);
        assert_eq!(row_charge(Some(2), (2, 4), exp(1.0), Some(mz)), 1);
        assert_eq!(row_charge(Some(4), (2, 4), exp(4.0), Some(mz)), 4);
        assert_eq!(row_charge(Some(4), (2, 4), Some(3456.789), Some(mz)), 4);
        assert_eq!(row_charge(Some(4), (2, 4), None, Some(mz)), 4);
        assert_eq!(row_charge(Some(4), (2, 4), exp(6.0), None), 4);
        assert_eq!(row_charge(None, (2, 4), exp(6.0), Some(mz)), 0);
    }

    #[test]
    fn anchors_are_confident_best_rows_per_file() {
        let mut rows = Vec::new();
        for s in 0..40u32 {
            rows.push(PinRow {
                spec: s,
                form: 0,
                charge: 2,
                decoy: s == 39,
                score: 100.0 - s as f64,
            });
            rows.push(PinRow {
                spec: s,
                form: 0,
                charge: 2,
                decoy: false,
                score: 0.0,
            });
        }
        let a = anchor_rows(&rows, 1, &|_| 0);
        assert_eq!(a.len(), 1);
        assert_eq!(a[0].len(), 39);
        assert!(a[0].iter().all(|&i| i % 2 == 0 && !rows[i].decoy));
    }

    #[test]
    fn an_anchor_is_calibrated_without_its_own_retention_time() {
        let rows: Vec<PinRow> = (0..90u32)
            .map(|f| PinRow {
                spec: f,
                form: f,
                charge: 2,
                decoy: false,
                score: 1.0,
            })
            .collect();
        let anchors: Vec<(usize, f64, f64)> = (0..90)
            .map(|i| (i, i as f64, if i == 45 { 90.0 } else { i as f64 }))
            .collect();
        let held = cross_fit_calibrations(&rows, &anchors);
        for (k, cal) in held.iter().enumerate() {
            let outside: Vec<(f64, f64)> = anchors
                .iter()
                .filter(|&&(i, _, _)| rt_fold(&rows[i]) != k)
                .map(|&(_, p, o)| (p, o))
                .collect();
            assert_eq!(cal, &RtCalibration::fit(&outside));
        }
        let all: Vec<(f64, f64)> = anchors.iter().map(|&(_, p, o)| (p, o)).collect();
        let full = RtCalibration::fit(&all).unwrap();
        let own = held[rt_fold(&rows[45])].as_ref().unwrap();
        assert!((own.apply(45.0) - 45.0).abs() < 1e-9);
        assert!(full.apply(45.0) > 50.0, "{}", full.apply(45.0));
    }

    #[test]
    fn impute_uses_the_median_of_scored_rows() {
        let mut col = vec![1.0, f64::NAN, 3.0, 100.0, 2.0];
        let miss = vec![false, true, false, true, false];
        impute_median(&mut col, &miss);
        assert_eq!(col, vec![1.0, 2.0, 3.0, 2.0, 2.0]);
    }

    #[test]
    fn split_before_field_matches_split_and_join() {
        for line in ["a\tb\tc\td", "a\t\t\tb", "\ta\tb", "a\tb\t"] {
            let x: Vec<&str> = line.split('\t').collect();
            for k in 0..x.len() {
                assert_eq!(
                    split_before_field(line, k),
                    Some((x[..k].join("\t").as_str(), x[k..].join("\t").as_str())),
                    "{line:?} k {k}"
                );
            }
            assert_eq!(split_before_field(line, x.len()), None);
        }
    }

    #[test]
    fn values_format_compactly() {
        let mut s = String::new();
        for v in [
            0.0,
            1.0,
            0.5,
            0.123456789,
            12.3456789,
            1234.56789,
            1e-7,
            2.5e7,
            -0.25,
        ] {
            push_value(&mut s, v);
            s.push(' ');
        }
        assert_eq!(
            s.trim_end(),
            "0 1 0.500000 0.123457 12.3457 1234.57 1.00000e-7 25000000 -0.250000"
        );
    }

    #[test]
    fn the_table_of_a_pin_describes_it_after_columns_are_added() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("x.pin");
        let mut text = String::from(
            "SpecId\tLabel\tScanNr\tRawScore\tcharge2\tcharge3\tPeptide\tProteins\n",
        );
        let mut spectra = Vec::new();
        let peptides = ["K.AEFVEVTK.L", "K.YLYEIAR.R", "K.AEFVEVTK.L", "R.HPEYAVSVLLR.L"];
        for (i, p) in peptides.iter().enumerate() {
            let scan = i as i32 + 1;
            let title = format!("t.{scan}");
            spectra.push(spectrum(&title, scan, 10.0 + i as f64));
            let (z2, z3) = if i % 2 == 0 { (1, 0) } else { (0, 1) };
            text.push_str(&format!(
                "{title}_{scan}_1\t1\t{scan}\t{}\t{z2}\t{z3}\t{p}\tP1\n",
                5 - i
            ));
        }
        std::fs::write(&pin, &text).unwrap();
        let cols = vec!["a_x".to_string(), "a_y".to_string()];
        let mut table = read_pin(&pin, &spectra, &cols).unwrap();
        write_pin(&pin, &table, &cols, &[vec![1.0; 4], vec![2.0; 4]]).unwrap();
        table.peptide_col += cols.len();
        let fresh = read_pin(&pin, &spectra, &["b".to_string()]).unwrap();
        assert_eq!(table.peptide_col, fresh.peptide_col);
        assert_eq!(table.rows.len(), fresh.rows.len());
        for (a, b) in table.rows.iter().zip(&fresh.rows) {
            assert_eq!(
                (a.spec, a.form, a.charge, a.decoy, a.score.to_bits()),
                (b.spec, b.form, b.charge, b.decoy, b.score.to_bits())
            );
        }
        let forms =
            |t: &PinTable| -> Vec<Form> { t.forms.iter().map(|f| f.0.clone()).collect() };
        assert_eq!(forms(&table), forms(&fresh));
        assert!(check_new_columns(&pin, &cols).is_err());
        assert!(check_new_columns(&pin, &["b".to_string()]).is_ok());
    }

    #[test]
    fn pin_gets_native_rt_columns_before_peptide() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("x.pin");
        let peptides = [
            "K.LGEYGFQNALIVR.Y",
            "K.AEFVEVTK.L",
            "R.HPEYAVSVLLR.L",
            "K.YLYEIAR.R",
            "K.SLHTLFGDELCK.V",
        ];
        let mut text =
            String::from("SpecId\tLabel\tScanNr\tRawScore\tcharge2\tPeptide\tProteins\n");
        let mut spectra = Vec::new();
        for i in 0..60 {
            let p = peptides[i % peptides.len()];
            let scan = i as i32 + 1;
            let title = format!("t.{scan}");
            spectra.push(spectrum(&title, scan, 10.0 + i as f64));
            let decoy = i >= 55;
            text.push_str(&format!(
                "{title}_{scan}_1\t{}\t{scan}\t{}\t1\t{p}\tP1\n",
                if decoy { -1 } else { 1 },
                200 - i
            ));
        }
        std::fs::write(&pin, &text).unwrap();
        let preds = Predictors {
            sources: vec![Source {
                pred: Box::new(NativeRt::from_bytes(ANDES_RT_WEIGHTS).unwrap()),
                sessions: Vec::new(),
                rt: true,
                ms2: false,
            }],
            tol: FragmentTolerance::Ppm(20.0),
            threads: 2,
            batch: 7,
            dia: false,
        };
        preds
            .add_to_pin(&pin, &spectra, &[spectra.len()], false)
            .unwrap();
        let out = std::fs::read_to_string(&pin).unwrap();
        let mut lines = out.lines();
        assert_eq!(
            lines.next().unwrap(),
            "SpecId\tLabel\tScanNr\tRawScore\tcharge2\tandesrt_rt_diff\tandesrt_predicted_rt\t\
             andesrt_missing\tPeptide\tProteins"
        );
        let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
        assert_eq!(rows.len(), 60);
        for r in &rows {
            assert_eq!(r.len(), 10);
            assert_eq!(r[7], "0");
            let d: f64 = r[5].parse().unwrap();
            assert!(d.is_finite() && d >= 0.0);
        }
        assert!(!dir.path().join("x.pin.predict.tmp").exists());
        assert!(preds
            .add_to_pin(&pin, &spectra, &[spectra.len()], false)
            .is_err());

        std::fs::write(&pin, &text).unwrap();
        for s in &mut spectra {
            s.rt_seconds = None;
        }
        preds
            .add_to_pin(&pin, &spectra, &[spectra.len()], false)
            .unwrap();
        let out = std::fs::read_to_string(&pin).unwrap();
        let mut lines = out.lines();
        assert_eq!(
            lines.next().unwrap(),
            "SpecId\tLabel\tScanNr\tRawScore\tcharge2\tandesrt_rt_diff\tandesrt_predicted_rt\t\
             andesrt_missing\tPeptide\tProteins"
        );
        let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
        assert_eq!(rows.len(), 60);
        for r in &rows {
            assert_eq!((r[5], r[6], r[7]), ("0", "0", "1"));
        }
    }
}
