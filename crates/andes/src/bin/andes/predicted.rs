//! Predicted-spectrum and retention-time features added to the PIN (`--predict`).
//!
//! After the PIN is written, each requested predictor scores every row
//! (targets and decoys alike) and the PIN is rewritten with the new columns
//! inserted before `Peptide`. Percolator then weighs them with the search
//! features; andes computes no FDR from them.
//!
//! - `rt`: the bundled andes-rt model ([`ANDES_RT_WEIGHTS`]).
//! - `ms2`, `peptdeep-rt`: AlphaPeptDeep graphs exported to ONNX.
//! - `deeplc-rt`: the DeepLC graph exported to ONNX.
//!
//! The third-party graphs are loaded from `--predictor-models`; they are not
//! part of andes.
//!
//! Columns per predictor (`<p>` = `andesrt`, `pdeep`, `deeplc`):
//! - RT: `<p>_rt_diff` = |observed − calibrated predicted RT| (minutes) and
//!   `<p>_predicted_rt`. Predicted RT is mapped to minutes per input file by
//!   an isotonic fit on anchor PSMs: per spectrum the row with the highest
//!   `RawScore`, kept when it is a target at target-decoy q ≤ [`ANCHOR_Q`]
//!   within its file. These anchors only place the calibration curve.
//! - MS2: `pdeep_<f>` for each spectral-match feature `f` in
//!   [`MS2_FEATURES`], comparing the predicted b/y (z ≤ 2) intensities with
//!   the most intense observed peak within the fragment tolerance.
//! - `<p>_missing`: 1 when the predictor could not score the row (an
//!   unsupported modification or length, no observed RT, no calibration);
//!   that row's other columns of the predictor hold the median over scored rows.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use model::{InstrumentType, Spectrum};
use onnx_infer::predictors::calibration::RtCalibration;
use onnx_infer::predictors::form::{parse_peptide, Form, Mass};
use onnx_infer::predictors::spectral::{ms2_features, FragmentTolerance, MS2_FEATURES};
use onnx_infer::predictors::{
    DeepLc, Engine, NativeRt, PeptDeep, Predictor, DEEPLC_FILE, PEPTDEEP_FILES,
};
use onnx_infer::{BackendKind, BackendOptions, OnnxSession};
use rayon::prelude::*;
use search::tdc::{qvalues, ScoredLabel};

use crate::cli::{Cli, PredictFeature, PredictorBackend};

/// andes-rt v1 weights (`ARTN0001`): a 77,794-parameter retention-time CNN
/// trained on retention times from public reanalysed datasets (see
/// `crates/onnx-infer/scripts/student/NOTES.md`).
static ANDES_RT_WEIGHTS: &[u8] =
    include_bytes!("../../../../../resources/predictors/andes_rt.artn");

/// Target-decoy q-value at or below which a file's best-scoring targets anchor
/// the RT calibration.
pub(crate) const ANCHOR_Q: f64 = 0.01;

/// One predictor and the PIN blocks it fills.
struct Source {
    pred: Box<dyn Predictor>,
    sessions: Vec<OnnxSession>,
    rt: bool,
    ms2: bool,
}

/// The predictors `--predict` asked for, loaded.
pub(crate) struct Predictors {
    sources: Vec<Source>,
    tol: FragmentTolerance,
    threads: usize,
    batch: usize,
}

/// The requested blocks, deduplicated; `None` when nothing is requested.
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

/// The AlphaPeptDeep instrument name for a resolved instrument.
fn peptdeep_instrument(i: InstrumentType) -> &'static str {
    match i {
        InstrumentType::QExactive => "QE",
        InstrumentType::TimsTOF => "timsTOF",
        InstrumentType::TOF => "SciexTOF",
        _ => "Lumos",
    }
}

impl Predictors {
    /// Load the predictors `cli.predict` names. `None` when it names none.
    /// A requested third-party model whose files are missing, or an ONNX
    /// engine this build lacks, is an error here, before the search runs.
    pub(crate) fn from_cli(cli: &Cli, instrument: InstrumentType) -> Result<Option<Self>, String> {
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
        let third_party = has(PredictFeature::Ms2)
            || has(PredictFeature::PeptDeepRt)
            || has(PredictFeature::DeepLcRt);
        let dir: Option<PathBuf> = cli.predictor_models.clone();
        if third_party && dir.is_none() {
            return Err(format!(
                "--predict {} loads exported ONNX models: pass --predictor-models DIR, a directory holding {} (AlphaPeptDeep) and/or {} \
                 (DeepLC), as written by the onnx-infer export scripts",
                want.iter()
                    .filter(|f| **f != PredictFeature::Rt)
                    .map(|f| match f {
                        PredictFeature::Ms2 => "ms2",
                        PredictFeature::PeptDeepRt => "peptdeep-rt",
                        _ => "deeplc-rt",
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
        if has(PredictFeature::Ms2) || has(PredictFeature::PeptDeepRt) {
            let dir = dir.as_deref().expect("checked above");
            let nce = cli.predictor_nce.unwrap_or(match instrument {
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
                has(PredictFeature::Ms2),
            )
            .map_err(|e| format!("--predict: {e}"))?;
            eprintln!(
                "Predict: AlphaPeptDeep from {} (instrument {inst}, NCE {nce}); modifications: {}",
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
                ms2: has(PredictFeature::Ms2),
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
        let high_res = instrument.is_high_resolution();
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
        }))
    }

    /// Whether the MS2 features need the spectra's peaks after the search.
    pub(crate) fn needs_peaks(&self) -> bool {
        self.sources.iter().any(|s| s.ms2)
    }

    /// One line naming the loaded blocks.
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

    /// The new PIN columns, in order.
    fn columns(&self) -> Vec<String> {
        let mut cols = Vec::new();
        for s in &self.sources {
            let p = s.pred.prefix();
            if s.ms2 {
                cols.extend(MS2_FEATURES.iter().map(|f| format!("{p}_{f}")));
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

    /// Add the predicted columns to the PIN at `pin` (rewritten in place).
    /// `spectra` are the searched spectra (with peaks when
    /// [`Self::needs_peaks`]) and `run_ends` the exclusive end index of each
    /// input file in `spectra`.
    pub(crate) fn add_to_pin(
        &self,
        pin: &Path,
        spectra: &[Spectrum],
        run_ends: &[usize],
    ) -> Result<(), String> {
        let t0 = Instant::now();
        let c0 = process_cpu_seconds();
        let new_cols = self.columns();
        let table = read_pin(pin, spectra, &new_cols)?;
        let n = table.rows.len();
        let run_of = |spec: u32| run_ends.partition_point(|&e| e <= spec as usize);
        let n_runs = run_ends.len().max(1);

        // Observed RT (minutes) and run per row.
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

        // Per predictor: its feature columns and its missing flags.
        let mut blocks: Vec<(Vec<Vec<f64>>, Vec<f64>)> = Vec::new();
        let mut report: Vec<String> = Vec::new();
        for src in &self.sources {
            let pred = src.pred.as_ref();
            let mut miss = vec![false; n];
            let mut cols: Vec<Vec<f64>> = Vec::new();
            if src.ms2 {
                let (ms2_cols, msg) = self.ms2_block(src, &table, spectra, &mut miss)?;
                cols.extend(ms2_cols);
                report.push(msg);
            }
            if src.rt {
                let raw = self.predict_rt(src, &table.forms)?;
                let mut diff = vec![f64::NAN; n];
                let mut cal_rt = vec![f64::NAN; n];
                let mut per_run = Vec::new();
                for (run, rows) in anchors.iter().enumerate() {
                    let pts: Vec<(f64, f64)> = rows
                        .iter()
                        .filter_map(|&i| {
                            let f = table.rows[i].form;
                            (f != NONE).then(|| (raw[f as usize] as f64, obs[i]))
                        })
                        .filter(|(p, o)| p.is_finite() && o.is_finite())
                        .collect();
                    let cal = RtCalibration::fit(&pts);
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
                    if let Some(c) = cal {
                        for (i, r) in table.rows.iter().enumerate() {
                            if r.spec == NONE || run_of(r.spec) != run || r.form == NONE {
                                continue;
                            }
                            let p = raw[r.form as usize] as f64;
                            if p.is_finite() && obs[i].is_finite() {
                                cal_rt[i] = c.apply(p);
                                diff[i] = (obs[i] - cal_rt[i]).abs();
                            }
                        }
                    }
                }
                for (m, d) in miss.iter_mut().zip(&diff) {
                    *m |= !d.is_finite();
                }
                report.push(format!("{} RT: {}", pred.prefix(), per_run.join("; ")));
                cols.push(diff);
                cols.push(cal_rt);
            }
            let n_miss = miss.iter().filter(|m| **m).count();
            report.push(format!("{} missing rows {n_miss} of {n}", pred.prefix()));
            if n > 0 && n_miss * 5 > n {
                eprintln!(
                    "WARN: --predict: {} could not score {n_miss} of {n} PIN rows ({:.0}%); \
                     a modification the model does not encode (for AlphaPeptDeep: one missing \
                     from peptdeep_mods.json) or missing retention times are the usual causes",
                    pred.prefix(),
                    100.0 * n_miss as f64 / n as f64
                );
            }
            let flags: Vec<f64> = miss.iter().map(|&m| f64::from(u8::from(m))).collect();
            // A missing row takes the column median over the predictor's scored rows.
            for col in &mut cols {
                impute_median(col, &flags);
            }
            blocks.push((cols, flags));
        }

        // Every configured predictor adds its columns, also when it scored no
        // row (e.g. no retention times in the input; its columns are then 0 and
        // its missing flag 1). The columns follow the --predict configuration
        // only, so the PINs of separately searched runs share one header.
        for (src, (_, flags)) in self.sources.iter().zip(&blocks) {
            if n > 0 && flags.iter().all(|&m| m != 0.0) {
                eprintln!(
                    "WARN: --predict: {} scored no PIN row; its columns are constant",
                    src.pred.prefix()
                );
            }
        }
        let mut values: Vec<Vec<f64>> = Vec::new();
        let mut flag_values: Vec<Vec<f64>> = Vec::new();
        for (cols, flags) in blocks {
            values.extend(cols);
            flag_values.push(flags);
        }
        values.extend(flag_values);

        write_pin(pin, &table, &new_cols, &values)?;
        for line in &report {
            eprintln!("Predict: {line}");
        }
        eprintln!(
            "Predict: added {} columns to {} PIN rows ({} peptidoforms) [PHASE predict: wall {:.2}s, cpu {:.1}s]",
            new_cols.len(),
            n,
            table.forms.len(),
            t0.elapsed().as_secs_f64(),
            process_cpu_seconds() - c0
        );
        Ok(())
    }

    /// RT prediction per distinct form (NaN when unsupported).
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
        eprintln!(
            "Predict: {} RT for {} of {} peptidoforms in {:.1}s",
            pred.prefix(),
            ix.len(),
            forms.len(),
            t.elapsed().as_secs_f64()
        );
        Ok(out.into_inner().expect("RT output lock"))
    }

    /// The MS2 feature columns; marks rows the model could not score.
    fn ms2_block(
        &self,
        src: &Source,
        table: &PinTable,
        spectra: &[Spectrum],
        miss: &mut [bool],
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
        let jobs = batches(
            fz.iter()
                .map(|(f, _)| pred.batch_key(&forms[*f as usize].0)),
            self.batch,
        );
        let out = Mutex::new(vec![Vec::new(); fz.len()]);
        let t = Instant::now();
        run_parallel(&src.sessions, self.threads, &jobs, |s, job| {
            let items: Vec<(&Form, u8)> = job
                .iter()
                .map(|&i| (&forms[fz[i].0 as usize].0, fz[i].1))
                .collect();
            let v = pred
                .ms2(s, &items)
                .map_err(|e| format!("{} MS2: {e}", pred.prefix()))?;
            let mut g = out.lock().expect("MS2 output lock");
            for (&i, x) in job.iter().zip(v) {
                g[i] = x;
            }
            Ok(())
        })?;
        let ms2_pred = out.into_inner().expect("MS2 output lock");
        let predict_s = t.elapsed().as_secs_f64();
        let tol = self.tol;
        let feats: Vec<Option<[f64; 8]>> = table
            .rows
            .par_iter()
            .map(|r| {
                let k = fz_ix.get(&(r.form, r.charge))?;
                let spec = &spectra[r.spec as usize];
                let (_, mass) = &forms[r.form as usize];
                Some(ms2_features(
                    &ms2_pred[*k as usize],
                    pred.ions(),
                    mass,
                    r.charge,
                    &spec.peaks,
                    tol,
                ))
            })
            .collect();
        let mut cols = vec![vec![f64::NAN; table.rows.len()]; MS2_FEATURES.len()];
        for (i, f) in feats.iter().enumerate() {
            match f {
                Some(v) => {
                    for (c, x) in cols.iter_mut().zip(v) {
                        c[i] = *x;
                    }
                }
                None => miss[i] = true,
            }
        }
        if no_peaks > 0 {
            eprintln!(
                "WARN: --predict ms2: {no_peaks} PIN rows have a spectrum without peaks; \
                 they are marked missing"
            );
        }
        let msg = format!(
            "{} MS2 for {} (peptidoform, charge) in {} batches, {predict_s:.1}s",
            pred.prefix(),
            fz.len(),
            jobs.len()
        );
        Ok((cols, msg))
    }
}

/// Row index sentinel: no spectrum / no parsed peptidoform.
const NONE: u32 = u32::MAX;

/// The fields of one PIN row the predictors need.
struct PinRow {
    spec: u32,
    form: u32,
    charge: u8,
    decoy: bool,
    score: f64,
}

/// Pass 1 over the PIN.
struct PinTable {
    rows: Vec<PinRow>,
    forms: Vec<(Form, Mass)>,
    /// Column index of `Peptide` (new columns go before it).
    peptide_col: usize,
}

/// `spec_id` base (`<title or scan=N>_<scan>`) → spectrum index; bases that
/// occur twice map to [`NONE`].
fn spectrum_bases(spectra: &[Spectrum]) -> HashMap<String, u32> {
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

/// The spectrum of a PIN row: its SpecId is the spectrum base followed by
/// `_<rank>` or `_<rank>_<row>`; the spectrum's scan must equal `ScanNr`.
fn locate(bases: &HashMap<String, u32>, spectra: &[Spectrum], spec_id: &str, scan: i32) -> u32 {
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

fn read_pin(pin: &Path, spectra: &[Spectrum], new_cols: &[String]) -> Result<PinTable, String> {
    let f = File::open(pin).map_err(|e| format!("{}: {e}", pin.display()))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("{}: empty PIN", pin.display()))?
        .map_err(|e| format!("{}: {e}", pin.display()))?;
    let h: Vec<&str> = header.split('\t').collect();
    let col = |name: &str| {
        h.iter()
            .position(|c| *c == name)
            .ok_or_else(|| format!("{}: PIN has no {name} column", pin.display()))
    };
    let (c_id, c_label, c_scan, c_pep) = (
        col("SpecId")?,
        col("Label")?,
        col("ScanNr")?,
        col("Peptide")?,
    );
    let c_score = col("RawScore").or_else(|_| col("RankScore"))?;
    for c in new_cols {
        if h.contains(&c.as_str()) {
            return Err(format!("{}: PIN already has a {c} column", pin.display()));
        }
    }
    let charge_cols: Vec<(usize, u8)> = h
        .iter()
        .enumerate()
        .filter_map(|(i, c)| Some((i, c.strip_prefix("charge")?.parse().ok()?)))
        .collect();
    let need = charge_cols
        .iter()
        .map(|(i, _)| *i)
        .chain([c_id, c_label, c_scan, c_pep, c_score])
        .max()
        .unwrap_or(0);
    let bases = spectrum_bases(spectra);
    let mut forms: Vec<(Form, Mass)> = Vec::new();
    let mut form_ix: HashMap<Form, u32> = HashMap::new();
    let mut by_peptide: HashMap<String, u32> = HashMap::new();
    let mut rows = Vec::new();
    let mut unlocated = 0usize;
    for (lineno, line) in lines.enumerate() {
        let line = line.map_err(|e| format!("{}: {e}", pin.display()))?;
        if line.is_empty() {
            continue;
        }
        let x: Vec<&str> = line.split('\t').collect();
        if x.len() <= need {
            return Err(format!(
                "{}: line {} has {} columns, expected more than {need}",
                pin.display(),
                lineno + 2,
                x.len()
            ));
        }
        let scan: i32 = x[c_scan].trim().parse().map_err(|_| {
            format!(
                "{}: line {} ScanNr {:?} is not an integer",
                pin.display(),
                lineno + 2,
                x[c_scan]
            )
        })?;
        let spec = locate(&bases, spectra, x[c_id], scan);
        unlocated += usize::from(spec == NONE);
        let form = match by_peptide.get(x[c_pep]) {
            Some(&f) => f,
            None => {
                let f = parse_peptide(x[c_pep]).map_or(NONE, |(f, m)| {
                    *form_ix.entry(f.clone()).or_insert_with(|| {
                        forms.push((f, m));
                        (forms.len() - 1) as u32
                    })
                });
                by_peptide.insert(x[c_pep].to_string(), f);
                f
            }
        };
        rows.push(PinRow {
            spec,
            form,
            charge: charge_cols
                .iter()
                .find(|(i, _)| x[*i] == "1")
                .map_or(0, |(_, z)| *z),
            decoy: x[c_label].trim() == "-1",
            score: x[c_score].trim().parse().unwrap_or(f64::NEG_INFINITY),
        });
    }
    if unlocated > 0 {
        eprintln!(
            "WARN: --predict: {unlocated} PIN rows match no searched spectrum; they are marked missing"
        );
    }
    Ok(PinTable {
        rows,
        forms,
        peptide_col: c_pep,
    })
}

/// Per input file, the rows that anchor the RT calibration: the best-scoring
/// row of each spectrum, when it is a target at q ≤ [`ANCHOR_Q`] among the
/// file's best rows.
fn anchor_rows(rows: &[PinRow], n_runs: usize, run_of: &dyn Fn(u32) -> usize) -> Vec<Vec<usize>> {
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

/// Replace non-finite values and rows flagged in `miss` by the median of the rest.
fn impute_median(col: &mut [f64], miss: &[f64]) {
    let mut ok: Vec<f64> = col
        .iter()
        .zip(miss)
        .filter(|(v, m)| v.is_finite() && **m == 0.0)
        .map(|(v, _)| *v)
        .collect();
    ok.sort_by(f64::total_cmp);
    let med = if ok.is_empty() {
        0.0
    } else if ok.len() % 2 == 1 {
        ok[ok.len() / 2]
    } else {
        0.5 * (ok[ok.len() / 2 - 1] + ok[ok.len() / 2])
    };
    for (v, m) in col.iter_mut().zip(miss) {
        if !v.is_finite() || *m != 0.0 {
            *v = med;
        }
    }
}

/// Pass 2: copy the PIN with the new columns inserted before `Peptide`, then
/// replace the original.
fn write_pin(
    pin: &Path,
    table: &PinTable,
    cols: &[String],
    values: &[Vec<f64>],
) -> Result<(), String> {
    let tmp = pin.with_file_name(format!(
        "{}.predict.tmp",
        pin.file_name()
            .map_or_else(|| "pin".into(), |n| n.to_string_lossy())
    ));
    let io = |e: std::io::Error| format!("{}: {e}", tmp.display());
    let f = File::open(pin).map_err(|e| format!("{}: {e}", pin.display()))?;
    let mut lines = BufReader::new(f).lines();
    let mut w = BufWriter::new(File::create(&tmp).map_err(io)?);
    let k = table.peptide_col;
    let header = lines
        .next()
        .ok_or_else(|| format!("{}: empty PIN", pin.display()))?
        .map_err(|e| format!("{}: {e}", pin.display()))?;
    let h: Vec<&str> = header.split('\t').collect();
    writeln!(
        w,
        "{}\t{}\t{}",
        h[..k].join("\t"),
        cols.join("\t"),
        h[k..].join("\t")
    )
    .map_err(io)?;
    let mut row = 0usize;
    let mut buf = String::new();
    for line in lines {
        let line = line.map_err(|e| format!("{}: {e}", pin.display()))?;
        if line.is_empty() {
            continue;
        }
        let x: Vec<&str> = line.split('\t').collect();
        buf.clear();
        buf.push_str(&x[..k].join("\t"));
        for v in values {
            buf.push('\t');
            push_value(&mut buf, v[row]);
        }
        buf.push('\t');
        buf.push_str(&x[k..].join("\t"));
        writeln!(w, "{buf}").map_err(io)?;
        row += 1;
    }
    if row != table.rows.len() {
        return Err(format!(
            "{}: changed while adding predicted features ({} rows, expected {})",
            pin.display(),
            row,
            table.rows.len()
        ));
    }
    w.flush().map_err(io)?;
    drop(w);
    std::fs::rename(&tmp, pin).map_err(|e| format!("{} -> {}: {e}", tmp.display(), pin.display()))
}

/// Integers as integers, other values with six significant digits.
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

/// Group item indices by `key` and chunk each group into batches of `batch`.
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

/// Run `jobs` on `threads` workers, each with its own handles on `master`.
/// The first error stops the remaining jobs and is returned.
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

/// Process CPU time (user + system) in seconds; 0 where `/proc` is absent.
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
        // ScanNr must agree with the spectrum.
        assert_eq!(locate(&b, &spectra, "run.10.10.2_10_1", 9), NONE);
        assert_eq!(locate(&b, &spectra, "other_5_1", 5), NONE);
    }

    #[test]
    fn anchors_are_confident_best_rows_per_file() {
        // 40 spectra in one file: best rows are targets with descending scores,
        // a decoy at the bottom; each spectrum also has a lower-scoring row.
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
    fn impute_uses_the_median_of_scored_rows() {
        let mut col = vec![1.0, f64::NAN, 3.0, 100.0, 2.0];
        let miss = vec![0.0, 1.0, 0.0, 1.0, 0.0];
        impute_median(&mut col, &miss);
        assert_eq!(col, vec![1.0, 2.0, 3.0, 2.0, 2.0]);
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
        };
        preds.add_to_pin(&pin, &spectra, &[spectra.len()]).unwrap();
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
        // A second pass refuses to add the same columns again.
        assert!(preds.add_to_pin(&pin, &spectra, &[spectra.len()]).is_err());

        // Without retention times no row is scored: the header is the same,
        // the columns are 0 and every row is flagged missing.
        std::fs::write(&pin, &text).unwrap();
        for s in &mut spectra {
            s.rt_seconds = None;
        }
        preds.add_to_pin(&pin, &spectra, &[spectra.len()]).unwrap();
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
