//! Match-between-runs stages shared by the stage subcommands and the
//! all-in-one run.
//!
//! Extract: for one run, every library precursor is measured on the run's
//! cached MS1 by the template fit of [`quant::evidence`], whether the run
//! identified it or not, and the rows go to `<stem>.extract.parquet`.
//!
//! Per library entry and run:
//!
//! * identified in the run (a donor run: a target PSM at PSM q-value within
//!   the library's donor threshold): an `Identified` target fit with the
//!   prior centred on the best donor PSM's RT and `s = max(w_prior, 6 s)`,
//!   plus a `PseudoTransfer` target fit with the alignment prior computed
//!   without this run's own observation;
//! * not identified: a `Transfer` target and its decoy probes, with the
//!   alignment prior;
//! * an unaligned run gets its `Identified` rows only (no pseudo-transfers,
//!   no transfers); pairs without a prediction, and identified pairs whose
//!   best PSM RT is not finite, are skipped and counted.
//!
//! An entry is fitted on the MS1 series of its FAIMS compensation voltage.
//! An entry without a voltage, or any entry in a run without FAIMS, is fitted
//! on the whole run; an entry whose voltage the FAIMS run lacks gets no rows
//! there and is counted (never fitted on another voltage's series or on the
//! mixed series). Only one copy of the run's MS1 is held at a time: the whole run
//! serves the width-prior fits of the whole-run entries; the scans are then
//! moved into one series per voltage for the voltage entries' width-prior
//! fits and rows; the split is dropped and the cache read again for the
//! whole-run entries' rows. A run with at most one voltage is never split.
//!
//! Width prior: the identified pairs are first fitted with the width free on
//! `{0.5, 0.7, 1, 1.4, 2} * w0` (`w0 = 6 s`, `s = 6 s`); `w_prior(t)` is the
//! fitted width through 20 equal-count bins of fitted apex RT (median per bin,
//! then a running median of 3 over the bins), interpolated linearly between
//! the bin medians and held constant beyond them. With fewer than 20 usable
//! fits, `w_prior = w0` everywhere.
//!
//! Finalize: the transfer targets and decoys of every run go to the
//! cross-validated transfer FDR of [`crate::transfer_fdr`] (features
//! [`TRANSFER_FEATURES`], built by [`transfer_features`]); rows without
//! evidence or with a non-finite feature take no part and get no q-value. Each
//! pseudo-transfer is scored by the model of its precursor's fold, which never
//! saw it, and placed against that fold's transfer scores. `<base>.mbr.tsv`
//! gets one row per identified or transferred (precursor, run) target, and
//! `<base>.mbr_validation.tsv` the pseudo-transfer acceptance rate, the
//! wrong-peak rate (`|t0_pseudo - t0_identified| > w_identified`) and the
//! accepted transfer counts of the pooled and per-family estimates. The
//! acceptance rate is optimistic for real transfers: the alignment map was
//! fitted with the pseudo-transfer's own run observation as an anchor, and
//! only the consensus RT leaves it out.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use output::mbr_out::{
    write_mbr_tsv, write_mbr_validation, Accepted, EvidenceType, MbrRow, MbrValidation,
    VALIDATION_Q,
};
use quant::decoys::{probes, Hypothesis, Probe};
use quant::evidence::{fit_many, Evidence, EvidenceParams, FitJob, WIDTH_GRID};
use quant::library::{Library, LibraryEntry};
use quant::ms1_index::Ms1RunIndex;
use rayon::prelude::*;

use crate::quant_run::ms1_cache::Ms1Cache;
use crate::transfer_fdr::{transfer_fdr, FdrCounts, FdrRow, TRANSFER_FEATURES};

pub(crate) use output::extract_table::{ExtractRole as Role, ExtractRow};

/// Default elution half-width `w0`, seconds.
const DEFAULT_W0_S: f64 = 6.0;
/// Floor of the RT uncertainty of an identified pair (apex vs trigger
/// offset), seconds.
const IDENTIFIED_MIN_S: f64 = 6.0;
/// Width grid of the width-prior fits, as multiples of `w0`.
const WIDTH_PRIOR_GRID: [f64; 5] = [0.5, 0.7, 1.0, 1.4, 2.0];
/// Equal-count RT bins of the width prior; also the fewest usable fits.
const WIDTH_PRIOR_BINS: usize = 20;
/// Smallest log-likelihood ratio of an identified fit that enters the width
/// prior; fits below it are not used.
const WIDTH_FIT_MIN_LLR: f64 = 25.0;

/// What one run's extraction skipped, for the caller to log.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct ExtractCounts {
    /// The run is unaligned: only its identified rows were produced.
    pub unaligned: bool,
    /// Entries whose FAIMS voltage the run lacks (no rows).
    pub foreign_cv: usize,
    /// Identified pairs whose best PSM RT is not finite (no identified and
    /// no pseudo-transfer row).
    pub nonfinite_best_rt: usize,
    /// Identified pairs without a prediction from the other runs (no
    /// pseudo-transfer row).
    pub no_pseudo_prediction: usize,
    /// Unidentified pairs without a prediction (no transfer rows).
    pub no_prediction: usize,
    /// Identified fits that set the run's width prior.
    pub width_fits: usize,
    /// Fewer than the needed identified fits: the width prior is the default.
    pub width_prior_fallback: bool,
}

impl ExtractCounts {
    fn add(&mut self, o: &ExtractCounts) {
        self.foreign_cv += o.foreign_cv;
        self.nonfinite_best_rt += o.nonfinite_best_rt;
        self.no_pseudo_prediction += o.no_pseudo_prediction;
        self.no_prediction += o.no_prediction;
    }
}

/// The rows of one run, in library order then probe order, and the counts.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RunExtract {
    pub rows: Vec<ExtractRow>,
    pub counts: ExtractCounts,
}

/// One in-memory copy of the run's MS1 (the whole run or its split). In tests
/// it counts the live copies of the calling thread.
struct Held;

impl Held {
    fn new() -> Self {
        #[cfg(test)]
        tracking::acquire();
        Held
    }
}

impl Drop for Held {
    fn drop(&mut self) {
        #[cfg(test)]
        tracking::release();
    }
}

/// Per-thread counters of MS1 copies and cache reads, for tests.
#[cfg(test)]
mod tracking {
    use std::cell::Cell;

    thread_local! {
        static LIVE: Cell<usize> = const { Cell::new(0) };
        static PEAK: Cell<usize> = const { Cell::new(0) };
        static LOADS: Cell<usize> = const { Cell::new(0) };
    }

    pub(super) fn acquire() {
        let live = LIVE.get() + 1;
        LIVE.set(live);
        PEAK.set(PEAK.get().max(live));
    }

    pub(super) fn release() {
        LIVE.set(LIVE.get() - 1);
    }

    pub(super) fn loaded() {
        LOADS.set(LOADS.get() + 1);
    }

    pub(super) fn reset() {
        PEAK.set(LIVE.get());
        LOADS.set(0);
    }

    /// `(peak live copies, cache reads)` since the last reset.
    pub(super) fn stats() -> (usize, usize) {
        (PEAK.get(), LOADS.get())
    }
}

/// Where an entry is fitted in one run.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Group {
    /// The whole run.
    Whole,
    /// The series at `voltages[k]` (the run's split order).
    Series(usize),
    /// Not at all: the FAIMS run lacks the entry's voltage.
    Foreign,
}

/// Per entry, where it is fitted in a run with FAIMS `voltages`. A run
/// without FAIMS, and an entry without a voltage, use the whole run; a run
/// with one voltage is fitted whole, and so is an entry at that voltage.
fn entry_groups(voltages: &[f32], entries: &[LibraryEntry]) -> Vec<Group> {
    entries
        .iter()
        .map(|e| match e.faims_cv {
            None => Group::Whole,
            Some(_) if voltages.is_empty() => Group::Whole,
            Some(cv) => match voltages.iter().position(|&v| v == cv) {
                None => Group::Foreign,
                Some(_) if voltages.len() == 1 => Group::Whole,
                Some(k) => Group::Series(k),
            },
        })
        .collect()
}

/// `w_prior(t)` as `(RT, width)` knots; empty means [`DEFAULT_W0_S`].
#[derive(Debug, Clone, PartialEq)]
struct WidthPrior {
    knots: Vec<(f64, f64)>,
}

impl WidthPrior {
    /// From fitted `(apex RT, width)` pairs; the default with fewer than
    /// [`WIDTH_PRIOR_BINS`] of them.
    fn from_fits(mut fits: Vec<(f64, f64)>) -> Self {
        let n = fits.len();
        if n < WIDTH_PRIOR_BINS {
            return Self { knots: Vec::new() };
        }
        fits.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        let mut ts = Vec::with_capacity(WIDTH_PRIOR_BINS);
        let mut ws = Vec::with_capacity(WIDTH_PRIOR_BINS);
        for k in 0..WIDTH_PRIOR_BINS {
            let bin = &fits[k * n / WIDTH_PRIOR_BINS..(k + 1) * n / WIDTH_PRIOR_BINS];
            ts.push(median(bin.iter().map(|f| f.0).collect()));
            ws.push(median(bin.iter().map(|f| f.1).collect()));
        }
        let smoothed: Vec<f64> = (0..ws.len())
            .map(|k| {
                if k == 0 || k + 1 == ws.len() {
                    ws[k]
                } else {
                    median(ws[k - 1..=k + 1].to_vec())
                }
            })
            .collect();
        Self {
            knots: ts.into_iter().zip(smoothed).collect(),
        }
    }

    fn is_default(&self) -> bool {
        self.knots.is_empty()
    }

    /// `w_prior` at RT `t`.
    fn at(&self, t: f64) -> f64 {
        let k = &self.knots;
        let Some((first, last)) = k.first().zip(k.last()) else {
            return DEFAULT_W0_S;
        };
        if t.is_nan() || t <= first.0 {
            return first.1;
        }
        if t >= last.0 {
            return last.1;
        }
        let i = k.partition_point(|x| x.0 <= t);
        let ((x0, y0), (x1, y1)) = (k[i - 1], k[i]);
        if x1 > x0 {
            y0 + (y1 - y0) * (t - x0) / (x1 - x0)
        } else {
            y0
        }
    }
}

/// Median of a non-empty vector (mean of the middle two for an even length).
fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// The RT of `entry`'s best PSM in `run`, when the run identified it.
fn identified_rt(entry: &LibraryEntry, run: u16) -> Option<f64> {
    entry
        .best_rt
        .iter()
        .find(|(r, _)| *r == run)
        .map(|&(_, rt)| rt)
}

/// Fits `jobs[k].1` on `series[jobs[k].0]`, in job order; one [`fit_many`]
/// sweep per series.
fn fit_on_series(
    series: &[&Ms1RunIndex],
    jobs: &[(usize, FitJob)],
    ep: &EvidenceParams,
) -> Vec<Option<Evidence>> {
    let mut out: Vec<Option<Evidence>> = vec![None; jobs.len()];
    for (g, ms1) in series.iter().enumerate() {
        let at: Vec<usize> = (0..jobs.len()).filter(|&k| jobs[k].0 == g).collect();
        if at.is_empty() {
            continue;
        }
        let group: Vec<FitJob> = at.iter().map(|&k| jobs[k].1).collect();
        for (k, e) in at.into_iter().zip(fit_many(ms1, &group, ep)) {
            out[k] = e;
        }
    }
    out
}

/// The usable width-free fits `(t0, w)` of the identified pairs among the
/// entries `ids`, each fitted on `series[series_of(id)]`.
fn width_fits<F>(
    ids: &[usize],
    run: u16,
    lib: &Library,
    ep: &EvidenceParams,
    series: &[&Ms1RunIndex],
    series_of: F,
) -> Vec<(f64, f64)>
where
    F: Fn(usize) -> usize,
{
    let jobs: Vec<(usize, FitJob)> = ids
        .iter()
        .filter_map(|&i| {
            let e = &lib.entries[i];
            let rt = identified_rt(e, run).filter(|rt| rt.is_finite())?;
            let meta = &lib.meta[e.precursor as usize];
            let probe = Probe {
                kind: Hypothesis::Target,
                mz: meta.mz,
                t_hat: rt,
                s: DEFAULT_W0_S.max(IDENTIFIED_MIN_S),
                w_prior: DEFAULT_W0_S,
            };
            let job = FitJob {
                probe,
                charge: meta.charge.max(1),
                envelope: &meta.envelope,
                widths: &WIDTH_PRIOR_GRID,
            };
            Some((series_of(i), job))
        })
        .collect();
    fit_on_series(series, &jobs, ep)
        .into_iter()
        .flatten()
        .filter(|ev| ev.amplitude > 0.0 && ev.llr >= WIDTH_FIT_MIN_LLR)
        .map(|ev| (ev.t0, ev.w))
        .collect()
}

type EntryRows = (Vec<ExtractRow>, ExtractCounts);

/// Fills `slots[id]` with the rows of each entry in `ids`, fitted on
/// `series[series_of(id)]`.
#[allow(clippy::too_many_arguments)]
fn rows_into<F>(
    slots: &mut [Option<EntryRows>],
    ids: &[usize],
    run: u16,
    aligned: bool,
    lib: &Library,
    wp: &WidthPrior,
    ep: &EvidenceParams,
    series: &[&Ms1RunIndex],
    series_of: F,
) where
    F: Fn(usize) -> usize + Sync,
{
    let planned: Vec<(Vec<(ExtractRow, Probe)>, ExtractCounts)> = ids
        .par_iter()
        .map(|&i| entry_rows(&lib.entries[i], run, aligned, lib, series[series_of(i)], wp))
        .collect();
    let mut jobs = Vec::new();
    for (&i, (rows, _)) in ids.iter().zip(&planned) {
        let meta = &lib.meta[lib.entries[i].precursor as usize];
        for (_, probe) in rows {
            let job = FitJob {
                probe: *probe,
                charge: meta.charge.max(1),
                envelope: &meta.envelope,
                widths: &WIDTH_GRID,
            };
            jobs.push((series_of(i), job));
        }
    }
    let mut evidence = fit_on_series(series, &jobs, ep).into_iter();
    for (&i, (rows, counts)) in ids.iter().zip(planned) {
        let rows = rows
            .into_iter()
            .map(|(mut row, _)| {
                row.evidence = evidence.next().flatten();
                row
            })
            .collect();
        slots[i] = Some((rows, counts));
    }
}

/// The rows of one library entry in `run`, each with the probe to fit on
/// `ms1` and no evidence yet. An unaligned run gets the identified row only.
fn entry_rows(
    e: &LibraryEntry,
    run: u16,
    run_aligned: bool,
    lib: &Library,
    ms1: &Ms1RunIndex,
    wp: &WidthPrior,
) -> (Vec<(ExtractRow, Probe)>, ExtractCounts) {
    let mut counts = ExtractCounts::default();
    let mut rows = Vec::new();
    let meta = &lib.meta[e.precursor as usize];
    let charge = meta.charge.max(1);
    let aligned = |r: u16| {
        lib.alignment
            .maps
            .get(r as usize)
            .is_some_and(|m| m.is_some())
    };
    let donors = e.donors.iter().filter(|&&r| aligned(r)).count() as u16;
    let other_donors = e.donors.iter().filter(|&&r| r != run && aligned(r)).count() as u16;
    let row = |kind: Hypothesis, role: Role, probe: &Probe, n_donor_runs: u16| {
        let row = ExtractRow {
            precursor: e.precursor,
            run,
            kind,
            role,
            evidence: None,
            n_donor_runs,
            library_q: e.q,
        };
        (row, *probe)
    };
    let target = |t_hat: f64, s: f64, w_prior: f64| Probe {
        kind: Hypothesis::Target,
        mz: meta.mz,
        t_hat,
        s,
        w_prior,
    };
    match identified_rt(e, run) {
        Some(rt) if !rt.is_finite() => counts.nonfinite_best_rt += 1,
        Some(rt) => {
            let w = wp.at(rt);
            let probe = target(rt, w.max(IDENTIFIED_MIN_S), w);
            rows.push(row(Hypothesis::Target, Role::Identified, &probe, donors));
            if !run_aligned {
                return (rows, counts);
            }
            match lib.alignment.predict_excluding(e.precursor, run) {
                Some((t_hat, s)) => {
                    let probe = target(t_hat, s, wp.at(t_hat));
                    rows.push(row(
                        Hypothesis::Target,
                        Role::PseudoTransfer,
                        &probe,
                        other_donors,
                    ));
                }
                None => counts.no_pseudo_prediction += 1,
            }
        }
        None if !run_aligned => {}
        None => match lib.alignment.predict(e.precursor, run) {
            Some((t_hat, s)) => {
                let run_rt = ms1.rt_range().unwrap_or((f64::NAN, f64::NAN));
                for probe in probes(
                    e.precursor,
                    run,
                    meta.mz,
                    charge,
                    t_hat,
                    s,
                    wp.at(t_hat),
                    run_rt,
                ) {
                    rows.push(row(probe.kind, Role::Transfer, &probe, other_donors));
                }
            }
            None => counts.no_prediction += 1,
        },
    }
    (rows, counts)
}

/// Measures every library precursor in run `run` from its MS1 cache.
///
/// An unaligned run yields its identified rows only and `counts.unaligned`.
/// Errors: `run` outside the library's runs, an entry whose precursor has no
/// metadata, or an unreadable cache (the message names the file).
pub(crate) fn extract_run(
    run: u16,
    cache: &Path,
    lib: &Library,
    ep: &EvidenceParams,
) -> Result<RunExtract, String> {
    let maps = &lib.alignment.maps;
    let Some(map) = maps.get(run as usize) else {
        return Err(format!(
            "run {run} is outside the library's {} runs",
            maps.len()
        ));
    };
    if let Some(e) = lib
        .entries
        .iter()
        .find(|e| e.precursor as usize >= lib.meta.len())
    {
        return Err(format!(
            "library entry {} has no precursor metadata ({} precursors)",
            e.precursor,
            lib.meta.len()
        ));
    }
    let aligned = map.is_some();
    let load = || {
        #[cfg(test)]
        tracking::loaded();
        Ms1Cache::load_from(cache).map_err(|e| e.to_string())
    };
    let n = lib.entries.len();
    let ms1 = load()?;
    let held = Held::new();
    let groups = entry_groups(&ms1.faims_voltages(), &lib.entries);
    let whole_ids: Vec<usize> = (0..n).filter(|&i| groups[i] == Group::Whole).collect();
    let cv_ids: Vec<usize> = (0..n)
        .filter(|&i| matches!(groups[i], Group::Series(_)))
        .collect();
    let foreign = groups.iter().filter(|&&g| g == Group::Foreign).count();
    let mut slots: Vec<Option<EntryRows>> = (0..n).map(|_| None).collect();
    let finish = |slots, width_fits: usize, wp: &WidthPrior| {
        let mut out = assemble(slots, width_fits, wp.is_default());
        out.counts.unaligned = !aligned;
        out.counts.foreign_cv = foreign;
        out
    };

    let mut fits = width_fits(&whole_ids, run, lib, ep, &[&ms1], |_| 0);
    let width_fit_count;
    if cv_ids.is_empty() {
        width_fit_count = fits.len();
        let wp = WidthPrior::from_fits(fits);
        rows_into(
            &mut slots,
            &whole_ids,
            run,
            aligned,
            lib,
            &wp,
            ep,
            &[&ms1],
            |_| 0,
        );
        drop(ms1);
        drop(held);
        return Ok(finish(slots, width_fit_count, &wp));
    }

    drop(held);
    let split = ms1.split_by_faims();
    let held = Held::new();
    let series: Vec<&Ms1RunIndex> = split.iter().map(|(_, s)| s).collect();
    let series_of = |i: usize| match groups[i] {
        Group::Series(k) => k,
        g => unreachable!("entry {i} is fitted on {g:?}"),
    };
    fits.extend(width_fits(&cv_ids, run, lib, ep, &series, series_of));
    width_fit_count = fits.len();
    let wp = WidthPrior::from_fits(fits);
    rows_into(
        &mut slots, &cv_ids, run, aligned, lib, &wp, ep, &series, series_of,
    );
    drop(series);
    drop(split);
    drop(held);

    if !whole_ids.is_empty() {
        let ms1 = load()?;
        let _held = Held::new();
        rows_into(
            &mut slots,
            &whole_ids,
            run,
            aligned,
            lib,
            &wp,
            ep,
            &[&ms1],
            |_| 0,
        );
    }
    Ok(finish(slots, width_fit_count, &wp))
}

/// Floor of `amplitude_se` in the amplitude signal-to-noise feature.
const MIN_AMPLITUDE_SE: f64 = 1e-12;

/// The transfer FDR features of a row, in [`TRANSFER_FEATURES`] order;
/// `None` without evidence or with a non-finite value.
pub(crate) fn transfer_features(row: &ExtractRow) -> Option<Vec<f32>> {
    let e = row.evidence.as_ref()?;
    let snr = e.amplitude / e.amplitude_se.max(MIN_AMPLITUDE_SE);
    let x = [
        e.llr,
        e.rt_z,
        e.rt_z.abs(),
        e.width_ratio,
        e.envelope_cos,
        e.coherence,
        e.m1_ratio,
        e.mass_err_ppm.abs(),
        e.n_points as f64,
        e.isotopes_used.count_ones() as f64,
        (snr + 1.0).log10(),
        row.n_donor_runs as f64,
        row.library_q,
    ];
    debug_assert_eq!(x.len(), TRANSFER_FEATURES.len());
    let x: Vec<f32> = x.iter().map(|&v| v as f32).collect();
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// `<base>.mbr.tsv` and `<base>.mbr_validation.tsv`.
pub(crate) fn mbr_paths(base: &Path) -> (PathBuf, PathBuf) {
    let with = |suffix: &str| {
        let mut s = base.as_os_str().to_owned();
        s.push(suffix);
        PathBuf::from(s)
    };
    (with(".mbr.tsv"), with(".mbr_validation.tsv"))
}

/// What one finalize did, for the caller to log.
#[derive(Debug, Clone, Default)]
pub(crate) struct FinalizeSummary {
    pub validation: MbrValidation,
    /// Rows of `<base>.mbr.tsv`.
    pub rows_written: usize,
    /// Transfer targets and decoys without evidence or with a non-finite
    /// feature (no q-value).
    pub no_features: usize,
    pub fdr: FdrCounts,
    /// The rows of `<base>.mbr.tsv`, in file order.
    pub features: Vec<MbrFeature>,
}

/// One row of `<base>.mbr.tsv` by precursor id and run index.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MbrFeature {
    pub precursor: u32,
    pub run: u16,
    pub evidence_type: EvidenceType,
    pub quantity: Option<f64>,
    pub quantity_se: Option<f64>,
    pub transfer_q: Option<f64>,
    pub apex_rt_s: Option<f64>,
    pub library_q: f64,
    /// Intensity-weighted isotope-0 mass error of the fit, ppm.
    pub mass_err_ppm: Option<f64>,
}

/// Transfer FDR, pseudo-transfer validation and the MBR outputs over the
/// extract rows of every run (`runs[i]` holds rows of any runs; a row's `run`
/// indexes `stems`, its `precursor` indexes `lib.meta`).
///
/// Errors: `mbr_fdr` outside `(0, 1]`, a row whose run or precursor is out of
/// range, or an unwritable output (the message names the file).
pub(crate) fn finalize(
    stems: &[String],
    lib: &Library,
    runs: &[Vec<ExtractRow>],
    mbr_fdr: f64,
    seed: u64,
    mbr_tsv: &Path,
    validation_tsv: &Path,
) -> Result<FinalizeSummary, String> {
    if !(mbr_fdr > 0.0 && mbr_fdr <= 1.0) {
        return Err(format!(
            "transfer q-value threshold {mbr_fdr} is outside (0, 1]"
        ));
    }
    let rows: Vec<&ExtractRow> = runs.iter().flatten().collect();
    if let Some(r) = rows
        .iter()
        .find(|r| r.run as usize >= stems.len() || r.precursor as usize >= lib.meta.len())
    {
        return Err(format!(
            "extract row of precursor {} in run {} is outside the {} precursors and {} runs",
            r.precursor,
            r.run,
            lib.meta.len(),
            stems.len()
        ));
    }

    let mut summary = FinalizeSummary::default();
    let mut fdr_rows = Vec::new();
    // Per row: the index of its FdrRow, for transfer targets and decoys.
    let mut fdr_index: Vec<Option<usize>> = vec![None; rows.len()];
    let mut identified: HashMap<(u32, u16), &Evidence> = HashMap::new();
    for (i, r) in rows.iter().enumerate() {
        match (r.role, r.kind) {
            (Role::Identified, Hypothesis::Target) => {
                if let Some(e) = &r.evidence {
                    let previous = identified.insert((r.precursor, r.run), e);
                    debug_assert!(
                        previous.is_none(),
                        "two identified rows of precursor {} in run {}",
                        r.precursor,
                        r.run
                    );
                }
            }
            (Role::Transfer, kind) => {
                let Some(x) = transfer_features(r) else {
                    summary.no_features += 1;
                    continue;
                };
                let family = match kind {
                    Hypothesis::Target => 0,
                    Hypothesis::MassShiftDecoy(_) => 1,
                    Hypothesis::RtShiftDecoy => 2,
                };
                fdr_index[i] = Some(fdr_rows.len());
                fdr_rows.push(FdrRow {
                    precursor: r.precursor,
                    is_decoy: family != 0,
                    family,
                    x,
                });
            }
            _ => {}
        }
    }
    let fdr = transfer_fdr(&fdr_rows, seed);
    let accepted = |q: Option<f64>| q.is_some_and(|q| q <= mbr_fdr);
    // The validation thresholds, VALIDATION_Q then mbr_fdr, and their counts.
    let thresholds = [VALIDATION_Q, mbr_fdr];
    let mut at = [Accepted::default(), Accepted::default()];

    let mut v = MbrValidation {
        q_threshold: mbr_fdr,
        ..MbrValidation::default()
    };
    let mut out: Vec<(u32, u16, MbrRow<'_>, Option<f64>)> = Vec::new();
    for (i, r) in rows.iter().enumerate() {
        if r.kind != Hypothesis::Target {
            continue;
        }
        let e = r.evidence.as_ref();
        let (evidence_type, transfer_q) = match r.role {
            Role::Identified => (EvidenceType::Identified, None),
            Role::Transfer => {
                let q = fdr_index[i].and_then(|k| fdr.result.q[k]);
                v.transfers += 1;
                v.transfers_scored += usize::from(q.is_some());
                for (a, &t) in at.iter_mut().zip(&thresholds) {
                    let ok = |q: Option<f64>| usize::from(q.is_some_and(|q| q <= t));
                    a.pooled += ok(q);
                    if let Some(k) = fdr_index[i] {
                        a.mass_only += ok(fdr.result.q_mass_only[k]);
                        a.rt_only += ok(fdr.result.q_rt_only[k]);
                    }
                }
                let t = if accepted(q) {
                    EvidenceType::TransferMs1
                } else {
                    EvidenceType::BelowThreshold
                };
                (t, q)
            }
            Role::PseudoTransfer => {
                let q = transfer_features(r).and_then(|x| fdr.pseudo_q(r.precursor, &x));
                v.pseudo_transfers += 1;
                v.pseudo_scored += usize::from(q.is_some());
                let wrong = match (e, identified.get(&(r.precursor, r.run))) {
                    (Some(p), Some(id)) => Some((p.t0 - id.t0).abs() > id.w),
                    _ => None,
                };
                if let Some(w) = wrong {
                    v.wrong_peak_pairs += 1;
                    v.wrong_peak += usize::from(w);
                }
                for (a, &t) in at.iter_mut().zip(&thresholds) {
                    if !q.is_some_and(|q| q <= t) {
                        continue;
                    }
                    a.pseudo += 1;
                    if let Some(w) = wrong {
                        a.pairs += 1;
                        a.wrong_peak += usize::from(w);
                    }
                }
                continue;
            }
        };
        out.push((
            r.precursor,
            r.run,
            MbrRow {
                precursor_key: &lib.meta[r.precursor as usize].key,
                run: &stems[r.run as usize],
                evidence_type,
                quantity: e.map(|e| e.quantity),
                quantity_se: e.map(|e| e.quantity_se),
                transfer_q,
                apex_rt_s: e.map(|e| e.t0),
                library_q: r.library_q,
            },
            e.map(|e| e.mass_err_ppm),
        ));
    }
    let [at_q01, at_threshold] = at;
    v.at_q01 = at_q01;
    v.at_threshold = at_threshold;
    out.sort_by_key(|&(p, run, ..)| (p, run));
    summary.features = out
        .iter()
        .map(|(p, run, r, mass_err)| MbrFeature {
            precursor: *p,
            run: *run,
            evidence_type: r.evidence_type,
            quantity: r.quantity,
            quantity_se: r.quantity_se,
            transfer_q: r.transfer_q,
            apex_rt_s: r.apex_rt_s,
            library_q: r.library_q,
            mass_err_ppm: *mass_err,
        })
        .collect();
    let out: Vec<MbrRow<'_>> = out.into_iter().map(|(_, _, row, _)| row).collect();
    write_mbr_tsv(mbr_tsv, &out)?;
    write_mbr_validation(validation_tsv, &v)?;

    eprintln!(
        "[mbr] transfer FDR (internal quant-quality score) features: {}",
        TRANSFER_FEATURES.join(",")
    );
    for f in &fdr.counts.failed_folds {
        eprintln!("[mbr] transfer FDR fold without a model: {f}");
    }
    eprintln!(
        "[mbr] transfers {}: {} with a q-value; {} transfer rows without features, {} unscored, {} scored without a q-value",
        v.transfers,
        v.transfers_scored,
        summary.no_features,
        fdr.counts.unscored,
        fdr.counts.no_q
    );
    eprintln!(
        "[mbr] pseudo-transfers {}: wrong peak {} of {}",
        v.pseudo_transfers, v.wrong_peak, v.wrong_peak_pairs
    );
    let shown: &[(f64, &Accepted)] = if mbr_fdr == VALIDATION_Q {
        &[(VALIDATION_Q, &v.at_q01)]
    } else {
        &[(VALIDATION_Q, &v.at_q01), (mbr_fdr, &v.at_threshold)]
    };
    for &(t, a) in shown {
        eprintln!(
            "[mbr] at q <= {t}: {} transfers (mass-shift decoys only {}, RT-shift decoys only {}); {} pseudo-transfers (acceptance rate {}, optimistic for transfers: the alignment kept their runs' anchors), {} of {} accepted at the wrong peak",
            a.pooled,
            a.mass_only,
            a.rt_only,
            a.pseudo,
            v.pseudo_acceptance(a).map_or("NA".into(), |s| format!("{s:.4}")),
            a.wrong_peak,
            a.pairs
        );
    }
    summary.rows_written = out.len();
    summary.validation = v;
    summary.fdr = fdr.counts;
    Ok(summary)
}

/// Concatenates the per-entry rows in library order and sums the counts.
fn assemble(slots: Vec<Option<EntryRows>>, width_fits: usize, fallback: bool) -> RunExtract {
    let mut counts = ExtractCounts {
        width_fits,
        width_prior_fallback: fallback,
        ..ExtractCounts::default()
    };
    let mut rows = Vec::new();
    for (r, c) in slots.into_iter().flatten() {
        rows.extend(r);
        counts.add(&c);
    }
    RunExtract { rows, counts }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use model::mass::ISOTOPE;
    use model::scan::Ms1Scan;
    use output::extract_table::{read_extract_table, write_extract_table};
    use quant::align::AlignParams;
    use quant::decoys::Hypothesis;
    use quant::library::{build_library, PrecursorMeta, PsmIn};
    use quant::ms1_index::Ms1RunIndex;

    use super::*;
    use crate::quant_run::ms1_cache::Ms1Cache;

    const ENVELOPE: [f32; 4] = [0.55, 0.30, 0.11, 0.04];
    const CHARGE: u8 = 2;
    const WIDTH: f64 = 4.2;
    const AMPLITUDE: f64 = 1e6;
    /// Precursors 0..BOTH are identified in runs 0 and 1.
    const BOTH: u32 = 30;
    /// Identified in run 0 only.
    const ONLY_RUN0: u32 = 30;
    /// Identified in run 1 only.
    const ONLY_RUN1: u32 = 31;
    const N_PREC: u32 = 32;
    const N_RUNS: usize = 3;
    const FAIMS_CVS: [f32; 2] = [-45.0, -60.0];

    fn mz_of(p: u32) -> f64 {
        450.0 + 13.7 * p as f64
    }

    fn apex(p: u32, run: u16) -> f64 {
        200.0 + 30.0 * p as f64 + 4.0 * run as f64
    }

    struct Lcg(u64);
    impl Lcg {
        fn uniform(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// MS1 scans every second with dense random centroids and every
    /// precursor's isotope pattern at its apex in `run`. With `faims`, scans
    /// alternate between two voltages and precursor `p` elutes only in the
    /// scans at `FAIMS_CVS[p % 2]`.
    /// Precursors in `silent` give no signal in `run`.
    fn run_scans(run: u16, faims: bool, silent: &dyn Fn(u32) -> bool) -> Vec<Ms1Scan> {
        let mut rng = Lcg(17 + run as u64);
        (0..1300)
            .map(|k| {
                let t = k as f64;
                let cv = faims.then_some(FAIMS_CVS[k % 2]);
                let mut peaks = Vec::new();
                for _ in 0..40 {
                    peaks.push((
                        300.0 + 1200.0 * rng.uniform(),
                        (500.0 + 1000.0 * rng.uniform()) as f32,
                    ));
                }
                for p in 0..N_PREC {
                    let elutes = !silent(p) && cv.is_none_or(|cv| cv == FAIMS_CVS[p as usize % 2]);
                    let g = (-(t - apex(p, run)).powi(2) / (2.0 * WIDTH * WIDTH)).exp();
                    for (i, &e) in ENVELOPE.iter().enumerate() {
                        let signal = if elutes {
                            AMPLITUDE * e as f64 * g
                        } else {
                            0.0
                        };
                        let inten = 200.0 + signal + 50.0 * (rng.uniform() - 0.5);
                        peaks.push((mz_of(p) + i as f64 * ISOTOPE / CHARGE as f64, inten as f32));
                    }
                }
                Ms1Scan {
                    rt: t,
                    peaks,
                    faims_cv: cv,
                    ..Ms1Scan::default()
                }
            })
            .collect()
    }

    struct Fixture {
        _dir: tempfile::TempDir,
        lib: Library,
        caches: Vec<PathBuf>,
    }

    /// Three runs: precursors 0..30 identified in runs 0 and 1, 30 only in
    /// run 0, 31 only in run 1; run 2 identifies precursor 0 only and is
    /// unaligned. Run 1 uses FAIMS when `faims`.
    fn fixture(faims: bool) -> Fixture {
        fixture_with(faims, &|_| false, &|_| false)
    }

    /// The PSM q-value threshold of the donor PSMs.
    const DONOR_Q: f64 = 0.01;

    /// [`fixture`] where the precursors in `silent_in_run1` give no MS1
    /// signal in run 1 (they stay identified there), and the run-1 PSMs of
    /// the precursors in `weak_in_run1` have a PSM q-value above
    /// [`DONOR_Q`].
    fn fixture_with(
        faims: bool,
        silent_in_run1: &dyn Fn(u32) -> bool,
        weak_in_run1: &dyn Fn(u32) -> bool,
    ) -> Fixture {
        let dir = tempfile::tempdir().unwrap();
        let meta: Vec<PrecursorMeta> = (0..N_PREC)
            .map(|p| PrecursorMeta {
                key: format!("PEP{p}/2"),
                mz: mz_of(p),
                charge: CHARGE,
                envelope: ENVELOPE,
            })
            .collect();
        let mut psms = Vec::new();
        let mut psm = |p: u32, run: u16| {
            psms.push(PsmIn {
                run,
                precursor: p,
                rt: apex(p, run) + 2.0,
                score: 10.0,
                q: if run == 1 && weak_in_run1(p) {
                    0.2
                } else {
                    0.001
                },
                pep: 0.0,
                is_decoy: false,
                faims_cv: (faims && run == 1).then_some(FAIMS_CVS[p as usize % 2]),
            })
        };
        for p in 0..BOTH {
            psm(p, 0);
            psm(p, 1);
        }
        psm(ONLY_RUN0, 0);
        psm(ONLY_RUN1, 1);
        psm(0, 2);
        let ap = AlignParams {
            min_anchors: 10,
            ..AlignParams::default()
        };
        let lib = build_library(&psms, meta, N_RUNS, 1.0, DONOR_Q, None, &ap);
        let caches = (0..N_RUNS as u16)
            .map(|r| {
                let path = dir.path().join(format!("run{r}.ms1cache"));
                let silent = |p: u32| r == 1 && silent_in_run1(p);
                let index = Ms1RunIndex::new(run_scans(r, faims && r == 1, &silent));
                Ms1Cache::store_at(&index, &path).unwrap();
                path
            })
            .collect();
        Fixture {
            _dir: dir,
            lib,
            caches,
        }
    }

    fn params() -> EvidenceParams {
        EvidenceParams::default()
    }

    fn rows_of(rows: &[ExtractRow], p: u32) -> Vec<&ExtractRow> {
        rows.iter().filter(|r| r.precursor == p).collect()
    }

    #[test]
    fn identified_pairs_get_identified_and_pseudo_rows() {
        let fx = fixture(false);
        assert!(fx.lib.alignment.maps[1].is_some());
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        for p in 0..BOTH {
            let rows = rows_of(&out.rows, p);
            assert_eq!(rows.len(), 2, "precursor {p}: {rows:?}");
            let (id, pseudo) = (rows[0], rows[1]);
            assert_eq!((id.role, id.kind), (Role::Identified, Hypothesis::Target));
            assert_eq!(
                (pseudo.role, pseudo.kind),
                (Role::PseudoTransfer, Hypothesis::Target)
            );
            // Precursor 0's donor run 2 is unaligned and not counted.
            assert_eq!((id.n_donor_runs, pseudo.n_donor_runs), (2, 1));
            assert_eq!((id.run, pseudo.run), (1, 1));
            for r in [id, pseudo] {
                let e = r.evidence.as_ref().unwrap();
                assert!(e.llr > 100.0, "precursor {p}: {e:?}");
                assert!((e.t0 - apex(p, 1)).abs() <= 1.0, "precursor {p}: {e:?}");
                assert!((e.amplitude / AMPLITUDE - 1.0).abs() < 0.1, "{e:?}");
            }
        }
        // Seen only in run 1: no other run predicts it.
        let rows = rows_of(&out.rows, ONLY_RUN1);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].role, Role::Identified);
        assert_eq!(out.counts.no_pseudo_prediction, 1);
        assert!(!out.counts.unaligned);
    }

    #[test]
    fn missing_pairs_get_target_and_decoys() {
        let fx = fixture(false);
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        let rows = rows_of(&out.rows, ONLY_RUN0);
        assert_eq!(rows.len(), 3, "{rows:?}");
        assert!(rows
            .iter()
            .all(|r| r.role == Role::Transfer && r.n_donor_runs == 1));
        assert_eq!(rows[0].kind, Hypothesis::Target);
        assert!(matches!(rows[1].kind, Hypothesis::MassShiftDecoy(_)));
        assert_eq!(rows[2].kind, Hypothesis::RtShiftDecoy);
        let target = rows[0].evidence.as_ref().unwrap();
        assert!((target.t0 - apex(ONLY_RUN0, 1)).abs() <= 1.0, "{target:?}");
        assert!(target.llr > 100.0, "{target:?}");
        for d in &rows[1..] {
            let e = d.evidence.as_ref().unwrap();
            assert!(e.llr < target.llr / 10.0, "{d:?}");
        }
        // The run-1-only precursor is a transfer into run 0.
        let out0 = extract_run(0, &fx.caches[0], &fx.lib, &params()).unwrap();
        let rows = rows_of(&out0.rows, ONLY_RUN1);
        assert_eq!(rows.len(), 3);
        assert!(rows.iter().all(|r| r.role == Role::Transfer));
        let lq = fx
            .lib
            .entries
            .iter()
            .find(|e| e.precursor == ONLY_RUN1)
            .unwrap()
            .q;
        assert!(rows.iter().all(|r| r.library_q == lq));
    }

    #[test]
    fn unaligned_run_keeps_its_identified_rows() {
        let fx = fixture(false);
        assert!(fx.lib.alignment.maps[2].is_none());
        let out = extract_run(2, &fx.caches[2], &fx.lib, &params()).unwrap();
        assert!(out.counts.unaligned);
        // Run 2 identifies precursor 0 only: its identified row, and no
        // pseudo-transfer or transfer rows.
        assert_eq!(out.rows.len(), 1, "{:?}", out.rows);
        let r = &out.rows[0];
        assert_eq!(
            (r.precursor, r.run, r.role, r.kind),
            (0, 2, Role::Identified, Hypothesis::Target)
        );
        let e = r.evidence.as_ref().unwrap();
        assert!((e.t0 - apex(0, 2)).abs() <= 1.0, "{e:?}");
        assert!((e.amplitude / AMPLITUDE - 1.0).abs() < 0.1, "{e:?}");
        assert_eq!(
            (out.counts.no_prediction, out.counts.no_pseudo_prediction),
            (0, 0)
        );
        assert!(extract_run(N_RUNS as u16, &fx.caches[2], &fx.lib, &params()).is_err());

        // The identified row reaches the mbr.tsv.
        let runs: Vec<Vec<ExtractRow>> = (0..N_RUNS)
            .map(|r| {
                extract_run(r as u16, &fx.caches[r], &fx.lib, &params())
                    .unwrap()
                    .rows
            })
            .collect();
        let stems: Vec<String> = ["r0", "r1", "r2"].map(String::from).to_vec();
        let (mbr, val) = mbr_paths(&fx._dir.path().join("x"));
        finalize(&stems, &fx.lib, &runs, 0.01, 7, &mbr, &val).unwrap();
        let text = std::fs::read_to_string(&mbr).unwrap();
        let in_r2: Vec<&str> = text
            .lines()
            .filter(|l| l.split('\t').nth(1) == Some("r2"))
            .collect();
        assert_eq!(in_r2.len(), 1, "{text}");
        assert!(
            in_r2[0].starts_with("PEP0/2\tr2\tidentified\t"),
            "{in_r2:?}"
        );
        assert!(!in_r2[0].contains("\tNA\tNA\tNA"), "{in_r2:?}");
    }

    #[test]
    fn psms_above_the_donor_threshold_make_a_transfer() {
        // Precursor 7's only PSM in run 1 is above the PSM q-value threshold.
        let fx = fixture_with(false, &|_| false, &|p| p == 7);
        let e = fx.lib.entries.iter().find(|e| e.precursor == 7).unwrap();
        assert_eq!(e.donors, vec![0]);
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        let rows = rows_of(&out.rows, 7);
        assert_eq!(rows.len(), 3, "{rows:?}");
        assert!(rows.iter().all(|r| r.role == Role::Transfer), "{rows:?}");
        assert_eq!(rows[0].kind, Hypothesis::Target);
        let target = rows[0].evidence.as_ref().unwrap();
        assert!((target.t0 - apex(7, 1)).abs() <= 1.0, "{target:?}");
        // The other precursors are still identified in run 1.
        let rows = rows_of(&out.rows, 8);
        assert_eq!(rows[0].role, Role::Identified);
    }

    #[test]
    fn foreign_voltage_entries_are_skipped_and_counted() {
        let mut fx = fixture(true);
        // Run 1 is acquired at FAIMS_CVS only; run 0 has no FAIMS.
        for e in fx.lib.entries.iter_mut() {
            if e.precursor == ONLY_RUN0 || e.precursor == 3 {
                e.faims_cv = Some(-30.0);
            }
        }
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        assert!(rows_of(&out.rows, ONLY_RUN0).is_empty());
        assert!(rows_of(&out.rows, 3).is_empty());
        assert_eq!(out.counts.foreign_cv, 2);
        assert!(!rows_of(&out.rows, 4).is_empty());
        // A run without FAIMS fits them on the whole run.
        let out0 = extract_run(0, &fx.caches[0], &fx.lib, &params()).unwrap();
        assert_eq!(out0.counts.foreign_cv, 0);
        assert_eq!(rows_of(&out0.rows, ONLY_RUN0)[0].role, Role::Identified);
        assert_eq!(rows_of(&out0.rows, 3)[0].role, Role::Identified);
    }

    #[test]
    fn output_order_is_deterministic() {
        let fx = fixture(false);
        let a = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        let b = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        assert_eq!(a, b);
        // Library order: each entry's rows are contiguous and in entry order.
        let order: Vec<u32> = fx.lib.entries.iter().map(|e| e.precursor).collect();
        let mut seen: Vec<u32> = a.rows.iter().map(|r| r.precursor).collect();
        seen.dedup();
        let expected: Vec<u32> = order.into_iter().filter(|p| seen.contains(p)).collect();
        assert_eq!(seen, expected);
        assert_eq!(seen.len(), N_PREC as usize);
    }

    #[test]
    fn nonfinite_best_rt_is_skipped_and_counted() {
        let mut fx = fixture(false);
        let e = fx
            .lib
            .entries
            .iter_mut()
            .find(|e| e.precursor == 5)
            .unwrap();
        for (r, rt) in e.best_rt.iter_mut() {
            if *r == 1 {
                *rt = f64::NAN;
            }
        }
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        assert!(rows_of(&out.rows, 5).is_empty());
        assert_eq!(out.counts.nonfinite_best_rt, 1);
    }

    #[test]
    fn width_prior_follows_identified_fits() {
        let fx = fixture(false);
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        assert!(out.counts.width_fits >= 20, "{:?}", out.counts);
        assert!(!out.counts.width_prior_fallback);
        for r in out.rows.iter().filter(|r| r.role == Role::Identified) {
            let e = r.evidence.as_ref().unwrap();
            let w_prior = e.w / e.width_ratio;
            assert!(
                (w_prior / WIDTH - 1.0).abs() < 0.1,
                "w_prior {w_prior}: {e:?}"
            );
        }
    }

    #[test]
    fn faims_entry_is_fitted_on_its_voltage_series() {
        let fx = fixture(true);
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        for r in out.rows.iter().filter(|r| r.role == Role::Identified) {
            let e = r.evidence.as_ref().unwrap();
            assert!(
                (e.amplitude / AMPLITUDE - 1.0).abs() < 0.1,
                "precursor {}: {e:?}",
                r.precursor
            );
        }
    }

    #[test]
    fn rows_round_trip_through_the_extract_table() {
        let fx = fixture(false);
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        let path = fx._dir.path().join("run1.extract.parquet");
        write_extract_table(&path, &out.rows).unwrap();
        assert_eq!(read_extract_table(&path).unwrap(), out.rows);
    }

    #[test]
    fn faims_run_with_whole_run_entry_holds_one_ms1_copy() {
        let fx = fixture(true);
        let whole = fx
            .lib
            .entries
            .iter()
            .find(|e| e.precursor == ONLY_RUN0)
            .unwrap();
        assert_eq!(whole.faims_cv, None);
        tracking::reset();
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        let (peak, loads) = tracking::stats();
        assert_eq!(peak, 1, "at most one MS1 copy is live");
        assert_eq!(loads, 2, "the whole run is read again for its entries");
        // The whole-run entry is still measured, in library order.
        assert_eq!(rows_of(&out.rows, ONLY_RUN0).len(), 3);
        let mut seen: Vec<u32> = out.rows.iter().map(|r| r.precursor).collect();
        seen.dedup();
        let order: Vec<u32> = fx.lib.entries.iter().map(|e| e.precursor).collect();
        assert_eq!(seen, order);

        // Without FAIMS the cache is read once.
        let fx = fixture(false);
        tracking::reset();
        extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        assert_eq!(tracking::stats(), (1, 1));
    }

    #[test]
    fn width_prior_ignores_identified_pairs_without_signal() {
        // About 30% of run 1's identified pairs have no MS1 signal.
        let silent = |p: u32| p % 10 < 3;
        let fx = fixture_with(false, &silent, &|_| false);
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        assert!(out.counts.width_fits >= 20, "{:?}", out.counts);
        assert!(!out.counts.width_prior_fallback);
        for r in out.rows.iter().filter(|r| r.role == Role::Identified) {
            let e = r.evidence.as_ref().unwrap();
            let w_prior = e.w / e.width_ratio;
            assert!(
                (w_prior / WIDTH - 1.0).abs() < 0.1,
                "precursor {}: w_prior {w_prior}",
                r.precursor
            );
        }
    }

    #[test]
    fn faims_transfer_decoys_use_the_target_voltage_series() {
        let mut fx = fixture(true);
        // Precursor ONLY_RUN0 elutes at FAIMS_CVS[0] in run 1, where it is a
        // transfer.
        let e = fx
            .lib
            .entries
            .iter_mut()
            .find(|e| e.precursor == ONLY_RUN0)
            .unwrap();
        e.faims_cv = Some(FAIMS_CVS[ONLY_RUN0 as usize % 2]);
        let out = extract_run(1, &fx.caches[1], &fx.lib, &params()).unwrap();
        let rows = rows_of(&out.rows, ONLY_RUN0);
        assert_eq!(rows.len(), 3);
        let target = rows[0].evidence.as_ref().unwrap();
        assert!(
            (target.amplitude / AMPLITUDE - 1.0).abs() < 0.1,
            "{target:?}"
        );
        // The voltage series has one scan every 2 s; the whole run one every
        // second, which would put about 4w + 1 scans in [t0 - 2w, t0 + 2w].
        for r in &rows {
            let e = r.evidence.as_ref().unwrap();
            assert!(
                (e.n_points as f64) <= 2.0 * e.w + 1.5,
                "{:?}: {e:?}",
                r.kind
            );
        }
    }
}

#[cfg(test)]
mod finalize_tests {
    use quant::align::AlignParams;
    use quant::library::{build_library, PrecursorMeta};

    use super::*;

    const N: u32 = 900;
    /// Every `WRONG_EVERY`-th pseudo-transfer sits 60 s from its identified apex.
    const WRONG_EVERY: u32 = 10;

    struct Rng(u64);
    impl Rng {
        fn unit(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        }
        fn normal(&mut self) -> f64 {
            let (u1, u2) = (self.unit(), self.unit());
            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        }
    }

    /// Evidence of a real co-eluting pattern at `t0` (`present`) or of noise.
    fn evidence(g: &mut Rng, present: bool, t0: f64) -> Evidence {
        let (llr, cos, coh) = if present {
            (60.0 + 10.0 * g.normal(), 0.97 + 0.01 * g.normal(), 0.9)
        } else {
            (
                3.0 * g.normal().abs(),
                0.6 + 0.2 * g.normal(),
                0.2 * g.normal(),
            )
        };
        let amplitude = if present { 1e6 } else { 2e3 * g.unit() };
        Evidence {
            amplitude,
            amplitude_se: 1e4,
            quantity: amplitude * 10.0,
            quantity_se: 1e5,
            llr,
            rt_z: if present {
                0.3 * g.normal()
            } else {
                2.0 * g.normal()
            },
            width_ratio: 1.0 + 0.1 * g.normal(),
            envelope_cos: cos,
            coherence: coh,
            m1_ratio: 0.05 * g.unit(),
            mass_err_ppm: if present {
                g.normal()
            } else {
                5.0 * g.normal()
            },
            n_points: 12,
            isotopes_used: 0b1111,
            t0,
            w: 4.0,
        }
    }

    fn row(p: u32, run: u16, kind: Hypothesis, role: Role, e: Option<Evidence>) -> ExtractRow {
        ExtractRow {
            precursor: p,
            run,
            kind,
            role,
            evidence: e,
            n_donor_runs: 1,
            library_q: 0.001,
        }
    }

    /// Run 0 identifies every precursor (with a pseudo-transfer, every
    /// `WRONG_EVERY`-th at the wrong peak); in run 1 even precursors are
    /// present and odd ones absent, each with a mass-shift and an RT-shift
    /// decoy. Precursor 1 has no evidence in run 1.
    fn extract() -> Vec<Vec<ExtractRow>> {
        let mut g = Rng(5);
        let mut run0 = Vec::new();
        let mut run1 = Vec::new();
        for p in 0..N {
            let apex = 100.0 + p as f64;
            run0.push(row(
                p,
                0,
                Hypothesis::Target,
                Role::Identified,
                Some(evidence(&mut g, true, apex)),
            ));
            let pseudo_t0 = if p % WRONG_EVERY == 0 {
                apex + 60.0
            } else {
                apex + 0.5
            };
            run0.push(row(
                p,
                0,
                Hypothesis::Target,
                Role::PseudoTransfer,
                Some(evidence(&mut g, true, pseudo_t0)),
            ));
            let target = (p != 1).then(|| evidence(&mut g, p % 2 == 0, apex));
            run1.push(row(p, 1, Hypothesis::Target, Role::Transfer, target));
            run1.push(row(
                p,
                1,
                Hypothesis::MassShiftDecoy(7),
                Role::Transfer,
                Some(evidence(&mut g, false, apex)),
            ));
            run1.push(row(
                p,
                1,
                Hypothesis::RtShiftDecoy,
                Role::Transfer,
                Some(evidence(&mut g, false, apex + 300.0)),
            ));
        }
        vec![run0, run1]
    }

    fn library() -> Library {
        let meta = (0..N)
            .map(|p| PrecursorMeta {
                key: format!("PEP{p}/2"),
                mz: 500.0 + p as f64,
                charge: 2,
                envelope: [0.55, 0.30, 0.11, 0.04],
            })
            .collect();
        build_library(&[], meta, 2, 0.01, 0.01, None, &AlignParams::default())
    }

    fn stems() -> Vec<String> {
        vec!["a".into(), "b".into()]
    }

    fn table(path: &Path) -> Vec<Vec<String>> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .map(|l| l.split('\t').map(String::from).collect())
            .collect()
    }

    #[test]
    fn finalize_writes_one_row_per_target_pair_with_evidence_types() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let s = finalize(&stems(), &library(), &extract(), 0.01, 7, &mbr, &val).unwrap();
        let t = table(&mbr);
        assert_eq!(t[0][2], "evidence_type");
        let body = &t[1..];
        assert_eq!(body.len(), 2 * N as usize);
        assert_eq!(s.rows_written, 2 * N as usize);
        let (mut true_acc, mut false_acc) = (0, 0);
        for (k, r) in body.iter().enumerate() {
            let p = (k / 2) as u32;
            assert_eq!(r[0], format!("PEP{p}/2"));
            if r[1] == "a" {
                assert_eq!((r[2].as_str(), r[5].as_str()), ("identified", "NA"));
                continue;
            }
            assert_eq!(r[1], "b");
            assert!(r[2] == "transfer_ms1" || r[2] == "below_threshold", "{r:?}");
            if p == 1 {
                assert_eq!(&r[2..7], ["below_threshold", "NA", "NA", "NA", "NA"]);
            } else {
                assert_ne!(r[3], "NA", "a below-threshold row keeps its quantity");
            }
            if r[2] == "transfer_ms1" {
                if p % 2 == 0 {
                    true_acc += 1;
                } else {
                    false_acc += 1;
                }
            }
        }
        eprintln!("finalize: {true_acc} present and {false_acc} absent transfers accepted");
        assert!(
            true_acc >= 400,
            "{true_acc} of 450 present transfers accepted"
        );
        assert!(false_acc <= 5, "{false_acc} absent transfers accepted");
        assert_eq!(s.no_features, 1);
    }

    #[test]
    fn finalize_reports_pseudo_transfer_validation() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let s = finalize(&stems(), &library(), &extract(), 0.01, 7, &mbr, &val).unwrap();
        let v = &s.validation;
        eprintln!("validation: {v:?}");
        assert_eq!(v.pseudo_transfers, N as usize);
        assert_eq!(v.wrong_peak_pairs, N as usize);
        assert_eq!(v.wrong_peak, (N / WRONG_EVERY) as usize);
        assert!(v.pseudo_acceptance(&v.at_q01).unwrap() > 0.9, "{v:?}");
        assert_eq!(v.at_q01, v.at_threshold);
        assert_eq!(v.transfers, N as usize);
        assert_eq!(v.transfers_scored, N as usize - 1);
        assert!(v.at_q01.mass_only > 0 && v.at_q01.rt_only > 0);
        let t = table(&val);
        assert_eq!(t.len(), 2);
        let get = |name: &str| &t[1][t[0].iter().position(|h| h == name).unwrap()];
        assert_eq!(get("wrong_peak_rate"), "0.1");
        assert_eq!(get("transfers"), &N.to_string());
        assert_eq!(get("accepted_pooled_q01"), &v.at_q01.pooled.to_string());

        // A looser threshold keeps the q <= 0.01 counts and adds its own.
        let loose = finalize(&stems(), &library(), &extract(), 0.5, 7, &mbr, &val).unwrap();
        let w = &loose.validation;
        assert_eq!(w.q_threshold, 0.5);
        assert_eq!(w.at_q01, v.at_q01);
        assert!(w.at_threshold.pooled > w.at_q01.pooled, "{w:?}");
        let t = table(&val);
        let get = |name: &str| &t[1][t[0].iter().position(|h| h == name).unwrap()];
        assert_eq!(get("accepted_pooled_q01"), &v.at_q01.pooled.to_string());
        assert_eq!(
            get("accepted_pooled_at_threshold"),
            &w.at_threshold.pooled.to_string()
        );
    }

    #[test]
    fn finalize_is_deterministic_for_a_seed() {
        let dir = tempfile::tempdir().unwrap();
        let runs = extract();
        let (m1, v1) = mbr_paths(&dir.path().join("one"));
        let (m2, v2) = mbr_paths(&dir.path().join("two"));
        finalize(&stems(), &library(), &runs, 0.01, 3, &m1, &v1).unwrap();
        finalize(&stems(), &library(), &runs, 0.01, 3, &m2, &v2).unwrap();
        assert_eq!(std::fs::read(&m1).unwrap(), std::fs::read(&m2).unwrap());
        assert_eq!(std::fs::read(&v1).unwrap(), std::fs::read(&v2).unwrap());
    }

    #[test]
    fn finalize_rejects_out_of_range_rows_and_thresholds() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let mut runs = extract();
        runs[1].push(row(0, 2, Hypothesis::Target, Role::Transfer, None));
        let e = finalize(&stems(), &library(), &runs, 0.01, 1, &mbr, &val).unwrap_err();
        assert!(e.contains("run 2"), "{e}");
        runs[1].pop();
        runs[1].push(row(N, 1, Hypothesis::Target, Role::Transfer, None));
        let e = finalize(&stems(), &library(), &runs, 0.01, 1, &mbr, &val).unwrap_err();
        assert!(e.contains(&format!("precursor {N}")), "{e}");
        for bad in [0.0, -0.1, 1.5, f64::NAN] {
            assert!(finalize(&stems(), &library(), &extract(), bad, 1, &mbr, &val).is_err());
        }
        assert!(!mbr.exists());
    }

    #[test]
    fn transfer_features_follow_the_feature_names() {
        let mut g = Rng(1);
        let mut r = row(
            0,
            1,
            Hypothesis::Target,
            Role::Transfer,
            Some(evidence(&mut g, true, 50.0)),
        );
        r.n_donor_runs = 3;
        let e = r.evidence.clone().unwrap();
        let x = transfer_features(&r).unwrap();
        assert_eq!(x.len(), TRANSFER_FEATURES.len());
        let at = |name: &str| x[TRANSFER_FEATURES.iter().position(|n| *n == name).unwrap()];
        assert_eq!(at("abs_rt_z"), e.rt_z.abs() as f32);
        assert_eq!(at("n_isotopes_used"), 4.0);
        assert_eq!(at("n_donor_runs"), 3.0);
        assert_eq!(
            at("log10_amplitude_snr"),
            (1e6f64 / 1e4 + 1.0).log10() as f32
        );
        r.evidence.as_mut().unwrap().amplitude_se = 0.0;
        assert!(transfer_features(&r).unwrap().iter().all(|v| v.is_finite()));
        r.evidence.as_mut().unwrap().llr = f64::NAN;
        assert!(transfer_features(&r).is_none());
        r.evidence = None;
        assert!(transfer_features(&r).is_none());
    }
}
