//! Match-between-runs stages shared by the stage subcommands and the
//! all-in-one run.
//!
//! Extract: for one run, every library precursor is measured on the run's
//! MS1 table (`<stem>.ms1.parquet`) by the template fit of
//! [`quant::evidence`], whether the run identified it or not, and the rows go
//! to `<stem>.extract.parquet`.
//!
//! Per library entry and run:
//!
//! * identified in the run (a donor run: a target PSM at PSM q-value within
//!   the library's donor threshold): an `Identified` target fit with the
//!   prior centred on the best donor PSM's RT and `s = max(w_prior(t), dt)`
//!   (the run's own elution width there, floored at its MS1 scan spacing),
//!   with its decoy probes ([`quant::decoys::probes`], same prior), plus a
//!   `PseudoTransfer` target fit with the alignment prior computed without
//!   this run's own observation;
//! * not identified: a `Transfer` target and its decoy probes, with the
//!   alignment prior, when the [`TransferSelection`] keeps the pair (it
//!   decides from the library alone, before any fit, and skips the target
//!   and its decoy probes together);
//! * an unaligned run gets its `Identified` rows (target and decoys) only (no
//!   pseudo-transfers, no transfers); pairs without a prediction, and identified pairs whose
//!   best PSM RT is not finite, are skipped and counted.
//!
//! An entry is fitted on the MS1 series of its FAIMS compensation voltage.
//! An entry without a voltage, or any entry in a run without FAIMS, is fitted
//! on the whole run; an entry whose voltage the FAIMS run lacks gets no rows
//! there and is counted (never fitted on another voltage's series or on the
//! mixed series). Only one copy of the run's MS1 is held at a time: the scans
//! of a FAIMS run are moved into one series per voltage for the voltage
//! entries' width scale, width-prior fits and rows; the split is dropped and
//! the MS1 table read again for the whole-run entries' rows, which use the
//! same width prior. A run with at most one voltage is never split.
//!
//! Run scales: every elution-time scale of the fits is derived from the run's
//! own MS1 before any row is fitted (the extract log and the
//! `andes.extract.run_scales` tag record them):
//!
//! * `dt`, the median spacing of the run's MS1 scans of one FAIMS voltage:
//!   the floor of every RT prior width `s` and of the width scale;
//! * `w_run`, the width scale: up to 500 identified pairs, spread evenly
//!   over the library (identity decoys left out), are fitted with the width
//!   free on `{0.5, 0.7, 1, 1.4, 2} * w` with `s = max(w, dt)`, and `w` becomes
//!   the median fitted width of the usable fits (positive amplitude, LLR at
//!   least 25), floored at `dt`. The rounds start at `w = 6 s` and stop when a
//!   round moves `w` by less than 10 % (at most 6 rounds), so the result is
//!   the run's own width, wherever it starts. With fewer than 20 usable fits
//!   in a round, `w_run = max(6 s, dt)`.
//!
//! Width prior: the identified pairs (of the voltage series in a FAIMS run)
//! are fitted with the width free on `{0.5, 0.7, 1, 1.4, 2} * w_run` and
//! `s = max(w_run, dt)`; `w_prior(t)` is the fitted width through 20
//! equal-count bins of fitted apex RT (median per bin, then a running median
//! of 3 over the bins), interpolated linearly between the bin medians and held
//! constant beyond them. With fewer than 20 usable fits, `w_prior = w_run`
//! everywhere.
//!
//! Finalize: the transfer targets and decoys of every run go to the
//! cross-validated transfer FDR of [`crate::transfer_fdr`] (features
//! [`TRANSFER_FEATURES`], built by [`transfer_features`]); rows without
//! evidence or with a non-finite feature take no part and get no q-value.
//! The identified targets and their decoys go, separately, through the same
//! machinery with their own models (features [`IDENTIFIED_FEATURES`], built
//! by [`identified_features`]): their prior is centred on the PSM's RT, so
//! their feature distribution is not the transfers', and the transfer
//! q-values stay what they are without them. Each identified target gets the
//! quality q-value `quant_q`; a transfer's `quant_q` is its `transfer_q`.
//! With a quality threshold, an identified row above it (or without a
//! q-value) is `identified_below_threshold`. Each
//! pseudo-transfer is scored by the model of its precursor's fold, which never
//! saw it, and placed against that fold's transfer scores. `<base>.mbr.parquet`
//! gets one row per identified or transferred (precursor, run) target, and
//! `<base>.mbr_validation.parquet` the pseudo-transfer acceptance rate, the
//! wrong-peak rate (`|t0_pseudo - t0_identified| > w_identified`) and the
//! accepted transfer counts of the pooled and per-family estimates. The
//! acceptance rate is optimistic for real transfers: the alignment map was
//! fitted with the pseudo-transfer's own run observation as an anchor, and
//! only the consensus RT leaves it out.
//!
//! Reported quantity: one elution width per precursor. Its consensus width
//! `w_c` is the median, over its identified rows in aligned runs with a fit
//! of positive amplitude, of `w_r / f_r'(f_r^-1(t0_r))`: the fitted width
//! carried onto the consensus RT scale by the run map's local slope. A row
//! in aligned run `r` with apex `t0` reports `A * w * sqrt(2 pi)` and
//! `amplitude_se * w * sqrt(2 pi)` with `w = w_c * f_r'(f_r^-1(t0))`. A row
//! of a precursor without such a fit, and every row of an unaligned run,
//! reports the extract's own `quantity` and `quantity_se`. The extract
//! tables keep the per-run values.
//!
//! Normalised quantity: `quantity_normalized` is the reported quantity with
//! its run's offset at the row's apex RT removed ([`quant::normalize`]: a
//! level per run and an RT-local term). The offsets are fitted on the
//! identified targets with `quant_q <= 0.01` and a positive quantity; a run
//! without enough of them reports null. `quantity` stays the raw measurement.

#[cfg(test)]
use std::collections::HashMap;
use std::path::{Path, PathBuf};

use output::ms1_table::read_ms1_table;
use quant::align::RunMap;
use quant::decoys::{probes, Hypothesis, Probe, PEPTIDE_MASS_SLOPE};
use quant::evidence::{
    fit_many, fit_on_grid, Evidence, EvidenceParams, FitJob, FitScratch, WIDTH_GRID,
};
use quant::library::{Library, LibraryEntry};
use quant::ms1_index::Ms1RunIndex;
use quant::ms2::{elution_window, LibrarySpectrum, Ms2RunIndex};
use rayon::prelude::*;

use crate::transfer_fdr::{IDENTIFIED_FEATURES, TRANSFER_FEATURES};

pub(crate) use output::extract_table::{ExtractRole as Role, ExtractRow};

/// Elution width, seconds, at which the run's width-scale rounds start, and
/// the width scale of a run with too few usable identified fits.
const SEED_WIDTH_S: f64 = 6.0;
/// Width grid of the width-scale and width-prior fits, as multiples of the
/// current width scale.
const WIDTH_PRIOR_GRID: [f64; 5] = [0.5, 0.7, 1.0, 1.4, 2.0];
/// Most identified pairs fitted per round of the width-scale estimate.
const SCALE_PAIRS: usize = 500;
/// Most rounds of the width-scale estimate.
const SCALE_MAX_ROUNDS: usize = 6;
/// The width-scale estimate stops once a round changes it by less than this
/// ratio.
const SCALE_CONVERGED: f64 = 1.1;
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
    /// Identified pairs whose other donor runs are all in other alignment
    /// components (no pseudo-transfer row; not in `no_pseudo_prediction`).
    pub pseudo_other_component: usize,
    /// Unidentified pairs whose donor runs are all in other alignment
    /// components (no transfer rows; not in `no_prediction`).
    pub other_component: usize,
    /// Unidentified pairs with a prediction that the transfer selection
    /// skips (no transfer rows).
    pub not_selected: usize,
    /// Identified fits that set the run's width prior.
    pub width_fits: usize,
    /// Fewer than the needed identified fits: the width prior is the run's
    /// width scale everywhere.
    pub width_prior_fallback: bool,
    /// The run's MS1 scan spacing and width scale.
    pub scales: RunScales,
}

/// The elution-time scales of one run, derived from its MS1 before any row is
/// fitted (see the module docs).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct RunScales {
    /// Median spacing of the MS1 scans of one FAIMS voltage, seconds; `None`
    /// with fewer than two scans of any voltage.
    pub spacing_s: Option<f64>,
    /// Width scale `w_run`, seconds.
    pub width_s: f64,
    /// Usable fits of the last width-scale round.
    pub width_fits: usize,
    /// Width-scale rounds run.
    pub rounds: usize,
    /// Too few usable fits: `width_s` is [`SEED_WIDTH_S`] (floored at the
    /// spacing).
    pub fallback: bool,
}

impl RunScales {
    /// Floor of the RT prior width `s` and of the width scale, seconds: the
    /// scan spacing (0 when unknown).
    fn floor_s(&self) -> f64 {
        self.spacing_s.unwrap_or(0.0)
    }

    /// The `andes.extract.run_scales` tag value.
    pub(crate) fn describe(&self) -> String {
        let spacing = self
            .spacing_s
            .map_or_else(|| "unknown".to_string(), |s| format!("{s:.4}"));
        format!(
            "ms1_spacing_s={spacing};width_scale_s={:.4};width_scale_fits={};width_scale_rounds={};width_scale_fallback={}",
            self.width_s, self.width_fits, self.rounds, self.fallback
        )
    }
}

impl ExtractCounts {
    fn add(&mut self, o: &ExtractCounts) {
        self.foreign_cv += o.foreign_cv;
        self.nonfinite_best_rt += o.nonfinite_best_rt;
        self.no_pseudo_prediction += o.no_pseudo_prediction;
        self.no_prediction += o.no_prediction;
        self.pseudo_other_component += o.pseudo_other_component;
        self.other_component += o.other_component;
        self.not_selected += o.not_selected;
    }
}

/// Which unidentified (precursor, run) pairs of an aligned run extract
/// measures as transfers. It decides from the library alone, before any fit,
/// and a skipped pair gets no rows at all: neither its target nor its decoy
/// probes. The transfer FDR therefore sees the kept pairs exactly as without
/// a selection (each with its full set of decoys), and an identity-decoy
/// entry is selected by the same rule as a target entry, from its own donor
/// runs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TransferSelection {
    /// Skip the transfers of an entry whose best donor PSM is above
    /// [`TRANSFER_DONOR_Q`]: finalize never transfers such an entry.
    pub donor_gate: bool,
    /// Skip a pair with fewer other aligned donor runs than this fraction of
    /// the library's aligned runs other than the acceptor, rounded up; 0
    /// keeps every pair.
    pub min_donor_fraction: f64,
}

/// Default of [`TransferSelection::min_donor_fraction`]: every pair whose
/// entry finalize can transfer is measured.
pub(crate) const DEFAULT_MIN_DONOR_FRACTION: f64 = 0.0;

impl TransferSelection {
    /// Every pair with a prediction is measured.
    pub(crate) const OFF: TransferSelection = TransferSelection {
        donor_gate: false,
        min_donor_fraction: 0.0,
    };

    pub(crate) fn is_off(&self) -> bool {
        !self.donor_gate && self.min_donor_fraction <= 0.0
    }

    /// The fewest other aligned donor runs a pair needs in a library of
    /// `aligned_runs` aligned runs.
    fn min_donors(&self, aligned_runs: usize) -> u16 {
        let others = aligned_runs.saturating_sub(1) as f64;
        (self.min_donor_fraction * others)
            .ceil()
            .clamp(0.0, u16::MAX as f64) as u16
    }

    /// One line for the extract table metadata and the log.
    pub(crate) fn describe(&self, aligned_runs: usize) -> String {
        format!(
            "donor_gate={} (best donor PSM q <= {TRANSFER_DONOR_Q}), min_donor_fraction={} ({} of {} other aligned runs)",
            self.donor_gate,
            self.min_donor_fraction,
            self.min_donors(aligned_runs),
            aligned_runs.saturating_sub(1)
        )
    }
}

impl Default for TransferSelection {
    fn default() -> Self {
        TransferSelection {
            donor_gate: true,
            min_donor_fraction: DEFAULT_MIN_DONOR_FRACTION,
        }
    }
}

/// [`TransferSelection`] resolved for one library.
#[derive(Debug, Clone, Copy)]
struct PairGate {
    donor_gate: bool,
    min_donors: u16,
}

impl PairGate {
    fn new(sel: &TransferSelection, lib: &Library) -> Self {
        let aligned = lib.alignment.maps.iter().filter(|m| m.is_some()).count();
        PairGate {
            donor_gate: sel.donor_gate,
            min_donors: sel.min_donors(aligned),
        }
    }

    /// True when the transfer of `e` with `other_donors` other aligned donor
    /// runs is measured.
    fn keeps(&self, e: &LibraryEntry, other_donors: u16) -> bool {
        let donor_ok = !self.donor_gate || e.donor_psm.q <= TRANSFER_DONOR_Q;
        donor_ok && other_donors >= self.min_donors
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

/// Per-thread counters of MS1 copies and MS1 table reads, for tests.
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

    /// `(peak live copies, MS1 table reads)` since the last reset.
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

/// `w_prior(t)` as `(RT, width)` knots, or the run's width scale when there
/// are none; with the floor of an identified pair's RT prior width.
#[derive(Debug, Clone, PartialEq)]
struct WidthPrior {
    knots: Vec<(f64, f64)>,
    /// `w_prior` without knots, seconds.
    default_s: f64,
    /// Floor of an identified pair's `s`, seconds.
    floor_s: f64,
}

impl WidthPrior {
    /// From fitted `(apex RT, width)` pairs; the run's width scale everywhere
    /// with fewer than [`WIDTH_PRIOR_BINS`] of them.
    fn from_fits(mut fits: Vec<(f64, f64)>, scales: &RunScales) -> Self {
        let n = fits.len();
        let default_s = scales.width_s;
        let floor_s = scales.floor_s();
        if n < WIDTH_PRIOR_BINS {
            return Self {
                knots: Vec::new(),
                default_s,
                floor_s,
            };
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
            default_s,
            floor_s,
        }
    }

    fn is_default(&self) -> bool {
        self.knots.is_empty()
    }

    /// `w_prior` at RT `t`.
    fn at(&self, t: f64) -> f64 {
        let k = &self.knots;
        let Some((first, last)) = k.first().zip(k.last()) else {
            return self.default_s;
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

impl WidthPrior {
    /// RT prior width `s` of an identified pair whose PSM is at `t`: the
    /// elution width there, floored at the scan spacing.
    fn identified_s(&self, t: f64) -> f64 {
        self.at(t).max(self.floor_s)
    }
}

/// Median of a non-empty vector (mean of the middle two for an even length).
pub(crate) fn median(mut v: Vec<f64>) -> f64 {
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

/// The entries among `ids` that the run identified at a finite RT, identity
/// decoys left out, with that RT.
fn identified_targets(ids: &[usize], run: u16, lib: &Library) -> Vec<(usize, f64)> {
    ids.iter()
        .filter_map(|&i| {
            let e = &lib.entries[i];
            if e.identity_decoy {
                return None;
            }
            let rt = identified_rt(e, run).filter(|rt| rt.is_finite())?;
            Some((i, rt))
        })
        .collect()
}

/// The width-free fit job of an identified pair: prior centred on its PSM RT
/// `rt` with width `s`, widths [`WIDTH_PRIOR_GRID`] times `scale`.
fn width_job(lib: &Library, i: usize, rt: f64, scale: f64, s: f64) -> FitJob<'_> {
    let meta = &lib.meta[lib.entries[i].precursor as usize];
    FitJob {
        probe: Probe {
            kind: Hypothesis::Target,
            mz: meta.mz,
            t_hat: rt,
            s,
            w_prior: scale,
        },
        charge: meta.charge.max(1),
        envelope: &meta.envelope,
        widths: &WIDTH_PRIOR_GRID,
    }
}

/// Whether a width-free fit enters the width scale and the width prior.
fn usable_width_fit(ev: &Evidence) -> bool {
    ev.amplitude > 0.0 && ev.llr >= WIDTH_FIT_MIN_LLR
}

/// The run's width scale (see the module docs) from the identified pairs
/// among the entries `ids`, each fitted on `series[series_of(id)]`, with
/// `spacing_s` the run's MS1 scan spacing. Each round fits at most
/// [`SCALE_PAIRS`] pairs, every k-th of the identified ones in library
/// order, one at a time ([`fit_on_grid`] gives the result [`fit_many`]
/// would, without its sweep over the whole run).
fn width_scale<F>(
    ids: &[usize],
    run: u16,
    lib: &Library,
    ep: &EvidenceParams,
    series: &[&Ms1RunIndex],
    series_of: F,
    spacing_s: Option<f64>,
) -> RunScales
where
    F: Fn(usize) -> usize + Sync,
{
    let floor = spacing_s.unwrap_or(0.0);
    let all = identified_targets(ids, run, lib);
    let stride = all.len().div_ceil(SCALE_PAIRS).max(1);
    let pilot: Vec<(usize, f64)> = all.into_iter().step_by(stride).collect();
    let mut out = RunScales {
        spacing_s,
        width_s: SEED_WIDTH_S.max(floor),
        ..RunScales::default()
    };
    while out.rounds < SCALE_MAX_ROUNDS {
        let w = out.width_s;
        let widths: Vec<f64> = pilot
            .par_iter()
            .map_init(FitScratch::default, |scratch, &(i, rt)| {
                let job = width_job(lib, i, rt, w, w.max(floor));
                fit_on_grid(
                    series[series_of(i)],
                    &job.probe,
                    job.charge,
                    job.envelope,
                    ep,
                    job.widths,
                    scratch,
                )
            })
            .flatten()
            .filter(usable_width_fit)
            .map(|ev| ev.w)
            .collect();
        out.rounds += 1;
        out.width_fits = widths.len();
        if widths.len() < WIDTH_PRIOR_BINS {
            out.width_s = SEED_WIDTH_S.max(floor);
            out.fallback = true;
            return out;
        }
        let next = median(widths).max(floor);
        out.width_s = next;
        if (next / w).max(w / next) < SCALE_CONVERGED {
            break;
        }
    }
    out
}

/// The usable width-free fits `(t0, w)` of the identified pairs among the
/// entries `ids` (identity decoys left out), each fitted on
/// `series[series_of(id)]` with widths [`WIDTH_PRIOR_GRID`] times the run's
/// width scale and `s = max(w_run, dt)`.
fn width_fits<F>(
    ids: &[usize],
    run: u16,
    lib: &Library,
    ep: &EvidenceParams,
    series: &[&Ms1RunIndex],
    series_of: F,
    scales: &RunScales,
) -> Vec<(f64, f64)>
where
    F: Fn(usize) -> usize,
{
    let (w, s) = (scales.width_s, scales.width_s.max(scales.floor_s()));
    let jobs: Vec<(usize, FitJob)> = identified_targets(ids, run, lib)
        .into_iter()
        .map(|(i, rt)| (series_of(i), width_job(lib, i, rt, w, s)))
        .collect();
    fit_on_series(series, &jobs, ep)
        .into_iter()
        .flatten()
        .filter(usable_width_fit)
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
    gate: &PairGate,
    wp: &WidthPrior,
    ep: &EvidenceParams,
    series: &[&Ms1RunIndex],
    series_of: F,
) where
    F: Fn(usize) -> usize + Sync,
{
    let planned: Vec<(Vec<(ExtractRow, Probe)>, ExtractCounts)> = ids
        .par_iter()
        .map(|&i| {
            let e = &lib.entries[i];
            entry_rows(e, run, aligned, lib, gate, series[series_of(i)], wp)
        })
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
/// `ms1` and no evidence yet. An unaligned run gets the identified row only;
/// a transfer that `gate` skips gets no rows.
fn entry_rows(
    e: &LibraryEntry,
    run: u16,
    run_aligned: bool,
    lib: &Library,
    gate: &PairGate,
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
            ms2: None,
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
    let run_rt = || ms1.rt_range().unwrap_or((f64::NAN, f64::NAN));
    match identified_rt(e, run) {
        Some(rt) if !rt.is_finite() => counts.nonfinite_best_rt += 1,
        Some(rt) => {
            let w = wp.at(rt);
            let s = wp.identified_s(rt);
            for probe in probes(e.precursor, run, meta.mz, charge, rt, s, w, run_rt()) {
                rows.push(row(probe.kind, Role::Identified, &probe, donors));
            }
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
                None if lib
                    .alignment
                    .blocked_by_component(e.precursor, run, Some(run)) =>
                {
                    counts.pseudo_other_component += 1
                }
                None => counts.no_pseudo_prediction += 1,
            }
        }
        None if !run_aligned => {}
        None => match lib.alignment.predict(e.precursor, run) {
            Some(_) if !gate.keeps(e, other_donors) => counts.not_selected += 1,
            Some((t_hat, s)) => {
                for probe in probes(
                    e.precursor,
                    run,
                    meta.mz,
                    charge,
                    t_hat,
                    s,
                    wp.at(t_hat),
                    run_rt(),
                ) {
                    rows.push(row(probe.kind, Role::Transfer, &probe, other_donors));
                }
            }
            None if lib.alignment.blocked_by_component(e.precursor, run, None) => {
                counts.other_component += 1
            }
            None => counts.no_prediction += 1,
        },
    }
    (rows, counts)
}

/// Measures every library precursor in run `run` from its MS1 table: its
/// identified pairs, and the transfers that `selection` keeps.
///
/// An unaligned run yields its identified rows only and `counts.unaligned`.
/// Errors: `run` outside the library's runs, an entry whose precursor has no
/// metadata, or an unreadable MS1 table (the message names the file).
pub(crate) fn extract_run(
    run: u16,
    ms1_table: &Path,
    lib: &Library,
    ep: &EvidenceParams,
    selection: &TransferSelection,
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
    let gate = PairGate::new(selection, lib);
    let load = || {
        #[cfg(test)]
        tracking::loaded();
        read_ms1_table(ms1_table)
    };
    let n = lib.entries.len();
    let ms1 = load()?;
    let held = Held::new();
    let spacing = ms1.median_scan_spacing();
    let groups = entry_groups(&ms1.faims_voltages(), &lib.entries);
    let whole_ids: Vec<usize> = (0..n).filter(|&i| groups[i] == Group::Whole).collect();
    let cv_ids: Vec<usize> = (0..n)
        .filter(|&i| matches!(groups[i], Group::Series(_)))
        .collect();
    let foreign = groups.iter().filter(|&&g| g == Group::Foreign).count();
    let mut slots: Vec<Option<EntryRows>> = (0..n).map(|_| None).collect();
    let finish = |slots, width_fits: usize, wp: &WidthPrior, scales: RunScales| {
        let mut out = assemble(slots, width_fits, wp.is_default());
        out.counts.unaligned = !aligned;
        out.counts.foreign_cv = foreign;
        out.counts.scales = scales;
        out
    };

    if cv_ids.is_empty() {
        let scales = width_scale(&whole_ids, run, lib, ep, &[&ms1], |_| 0, spacing);
        let fits = width_fits(&whole_ids, run, lib, ep, &[&ms1], |_| 0, &scales);
        let width_fit_count = fits.len();
        let wp = WidthPrior::from_fits(fits, &scales);
        rows_into(
            &mut slots,
            &whole_ids,
            run,
            aligned,
            lib,
            &gate,
            &wp,
            ep,
            &[&ms1],
            |_| 0,
        );
        drop(ms1);
        drop(held);
        return Ok(finish(slots, width_fit_count, &wp, scales));
    }

    drop(held);
    let split = ms1.split_by_faims();
    let held = Held::new();
    let series: Vec<&Ms1RunIndex> = split.iter().map(|(_, s)| s).collect();
    let series_of = |i: usize| match groups[i] {
        Group::Series(k) => k,
        g => unreachable!("entry {i} is fitted on {g:?}"),
    };
    let scales = width_scale(&cv_ids, run, lib, ep, &series, series_of, spacing);
    let fits = width_fits(&cv_ids, run, lib, ep, &series, series_of, &scales);
    let width_fit_count = fits.len();
    let wp = WidthPrior::from_fits(fits, &scales);
    rows_into(
        &mut slots, &cv_ids, run, aligned, lib, &gate, &wp, ep, &series, series_of,
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
            &gate,
            &wp,
            ep,
            &[&ms1],
            |_| 0,
        );
    }
    Ok(finish(slots, width_fit_count, &wp, scales))
}

/// MS2 evidence of one run's transfer and pseudo-transfer rows, by role and
/// decoy family, for the extract log.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct Ms2Counts {
    /// Per role (0 transfer, 1 pseudo-transfer) and family (0 target, 1 mass
    /// shift, 2 RT shift): rows with MS2 evidence, of them with at least one
    /// scan in the window, and of them supported
    /// ([`quant::ms2::Ms2Evidence::supports`]).
    pub rows: [[usize; 3]; 2],
    pub with_scan: [[usize; 3]; 2],
    pub supported: [[usize; 3]; 2],
    /// Transfer and pseudo-transfer rows with a fit whose entry has no
    /// library spectrum (no MS2 evidence).
    pub no_spectrum: usize,
}

/// The m/z of a probe of precursor `meta` with hypothesis `kind`, as
/// [`quant::decoys::probes`] builds it.
fn probe_mz(meta: &quant::library::PrecursorMeta, kind: Hypothesis) -> f64 {
    match kind {
        Hypothesis::MassShiftDecoy(n) => {
            meta.mz + n as f64 * PEPTIDE_MASS_SLOPE / meta.charge.max(1) as f64
        }
        Hypothesis::Target | Hypothesis::RtShiftDecoy => meta.mz,
    }
}

/// Sets the MS2 evidence of every transfer and pseudo-transfer row with a
/// fit: the best scan of run `ms2` for the entry's library spectrum
/// (`spectra`, per precursor id) among the scans isolating the row's probe
/// m/z within its fitted elution window, at the entry's FAIMS voltage. A
/// decoy row is scored with its target's spectrum on its own m/z and window.
/// A pseudo-transfer leaves out the scans of its run's donor PSMs
/// (`donor_scans`, per precursor id), whose identifications make it a
/// pseudo-transfer; a transfer has none in its run. Rows of entries without
/// a spectrum, identified rows and rows without a fit keep `None`.
pub(crate) fn attach_ms2(
    rows: &mut [ExtractRow],
    lib: &Library,
    spectra: &[Option<LibrarySpectrum>],
    donor_scans: &[Vec<(u16, f64)>],
    ms2: &Ms2RunIndex,
    tol: model::Tolerance,
) -> Ms2Counts {
    let mut faims: Vec<Option<f32>> = vec![None; lib.meta.len()];
    for e in &lib.entries {
        faims[e.precursor as usize] = e.faims_cv;
    }
    rows.par_iter_mut().for_each(|r| {
        r.ms2 = None;
        if r.role == Role::Identified {
            return;
        }
        let (Some(e), Some(Some(lib_spec))) = (&r.evidence, spectra.get(r.precursor as usize))
        else {
            return;
        };
        let Some(window) = elution_window(e.t0, e.w) else {
            return;
        };
        let p = r.precursor as usize;
        let exclude: Vec<f64> = if r.role == Role::PseudoTransfer {
            donor_scans[p]
                .iter()
                .filter(|(run, _)| *run == r.run)
                .map(|&(_, t)| t)
                .collect()
        } else {
            Vec::new()
        };
        let mz = probe_mz(&lib.meta[p], r.kind);
        r.ms2 = Some(ms2.best_match(lib_spec, mz, window, faims[p], &exclude, tol));
    });
    let mut c = Ms2Counts::default();
    for r in rows.iter() {
        let role = match r.role {
            Role::Identified => continue,
            Role::Transfer => 0,
            Role::PseudoTransfer => 1,
        };
        let Some(m) = &r.ms2 else {
            if r.evidence.is_some() {
                c.no_spectrum += 1;
            }
            continue;
        };
        let f = family(r.kind) as usize;
        c.rows[role][f] += 1;
        c.with_scan[role][f] += usize::from(m.n_scans > 0);
        c.supported[role][f] += usize::from(m.supports());
    }
    c
}

/// `d run RT / d consensus RT` of run `run`'s map at run RT `rt`; `None` for
/// an unaligned run or a slope that is not finite and positive.
pub(crate) fn local_slope(maps: &[Option<RunMap>], run: u16, rt: f64) -> Option<f64> {
    let map = maps.get(run as usize)?.as_ref()?;
    let slope = map.slope_at(map.invert(rt));
    (slope.is_finite() && slope > 0.0).then_some(slope)
}

/// Per precursor, the consensus width `w_c` of the module docs, from the
/// identified target rows among `rows`.
#[cfg(test)]
pub(crate) fn consensus_widths<'a>(
    rows: impl IntoIterator<Item = &'a ExtractRow>,
    maps: &[Option<RunMap>],
) -> HashMap<u32, f64> {
    let mut widths: HashMap<u32, Vec<f64>> = HashMap::new();
    for r in rows {
        if let Some(w) = consensus_scale_width(r, maps) {
            widths.entry(r.precursor).or_default().push(w);
        }
    }
    widths.into_iter().map(|(p, w)| (p, median(w))).collect()
}

/// `(quantity, quantity_se)` of evidence `e` in run `run`: at the consensus
/// width `w_c` carried into the run, or the extract's own values without a
/// consensus width or in an unaligned run.
pub(crate) fn reported_quantity(
    e: &Evidence,
    run: u16,
    w_c: Option<f64>,
    maps: &[Option<RunMap>],
) -> (f64, f64) {
    match w_c.zip(local_slope(maps, run, e.t0)) {
        Some((w_c, slope)) => {
            let area = w_c * slope * (2.0 * std::f64::consts::PI).sqrt();
            (e.amplitude * area, e.amplitude_se * area)
        }
        None => (e.quantity, e.quantity_se),
    }
}

/// Floor of `amplitude_se` in the amplitude signal-to-noise feature.
const MIN_AMPLITUDE_SE: f64 = 1e-12;

/// Shape features of a row without a donor reference: no width or tailing
/// difference, and an isotope cosine of -1 (below any real cosine).
pub(crate) const NO_SHAPE: [f64; 3] = [0.0, 0.0, -1.0];

/// The peak shapes of every precursor's identified fits, the donor reference
/// of the transfer shape features: per precursor, `(run, width on the
/// consensus scale, tailing, isotope profile)` of its identified targets in
/// aligned runs with a positive amplitude.
#[cfg(test)]
#[derive(Debug, Clone, Default)]
pub(crate) struct DonorShapes {
    fits: HashMap<u32, Vec<DonorFit>>,
}

/// One identified fit: run, width on the consensus scale,
/// tailing and isotope profile.
pub(crate) type DonorFit = (u16, f64, f64, [f64; 4]);

/// Isotope heights clamped at 0 and scaled to sum 1; `None` when all are 0
/// or not finite.
pub(crate) fn iso_profile(h: &[f64; 4]) -> Option<[f64; 4]> {
    let c = h.map(|v| if v.is_finite() { v.max(0.0) } else { 0.0 });
    let sum: f64 = c.iter().sum();
    (sum > 0.0).then(|| c.map(|v| v / sum))
}

/// `a . b / (|a| |b|)`; -1 when either is zero.
pub(crate) fn cosine(a: &[f64; 4], b: &[f64; 4]) -> f64 {
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let na = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na > 0.0 && nb > 0.0 {
        dot / (na * nb)
    } else {
        -1.0
    }
}

/// The donor fit of an identified target row: run, width on the consensus
/// scale, tailing and isotope profile; `None` for any other row, or without
/// a positive amplitude, a finite positive width, a finite tailing, the
/// run's map or an isotope profile.
pub(crate) fn donor_fit(r: &ExtractRow, maps: &[Option<RunMap>]) -> Option<DonorFit> {
    if (r.role, r.kind) != (Role::Identified, Hypothesis::Target) {
        return None;
    }
    let e = r.evidence.as_ref()?;
    if !(e.amplitude > 0.0 && e.w.is_finite() && e.w > 0.0 && e.tailing.is_finite()) {
        return None;
    }
    let slope = local_slope(maps, r.run, e.t0)?;
    let iso = iso_profile(&e.iso_heights)?;
    Some((r.run, e.w / slope, e.tailing, iso))
}

/// The width of an identified target row on the consensus scale, the input
/// of the consensus width `w_c`; `None` for any other row, or without a
/// positive amplitude, a finite positive width or the run's map.
pub(crate) fn consensus_scale_width(r: &ExtractRow, maps: &[Option<RunMap>]) -> Option<f64> {
    if (r.role, r.kind) != (Role::Identified, Hypothesis::Target) {
        return None;
    }
    let e = r.evidence.as_ref()?;
    if !(e.amplitude > 0.0 && e.w.is_finite() && e.w > 0.0) {
        return None;
    }
    local_slope(maps, r.run, e.t0).map(|slope| e.w / slope)
}

#[cfg(test)]
impl DonorShapes {
    /// From the identified target rows among `rows`.
    pub(crate) fn from_rows<'a>(
        rows: impl IntoIterator<Item = &'a ExtractRow>,
        maps: &[Option<RunMap>],
    ) -> Self {
        let mut fits: HashMap<u32, Vec<DonorFit>> = HashMap::new();
        for r in rows {
            if let Some(f) = donor_fit(r, maps) {
                fits.entry(r.precursor).or_default().push(f);
            }
        }
        Self { fits }
    }
}

/// The donor reference of the transfer shape features.
pub(crate) trait ShapeReference {
    /// `[log width ratio, tailing difference, isotope cosine]` of `row`'s fit
    /// against the median of its precursor's identified fits in the other
    /// runs (the row's own run left out, so a pseudo-transfer is not
    /// compared with the identification it shadows): the width on the
    /// consensus scale, the tailing and the isotope profile. [`NO_SHAPE`]
    /// without a reference, a positive amplitude or the run's map.
    fn features(&self, row: &ExtractRow, maps: &[Option<RunMap>]) -> [f64; 3];
}

#[cfg(test)]
impl ShapeReference for DonorShapes {
    fn features(&self, row: &ExtractRow, maps: &[Option<RunMap>]) -> [f64; 3] {
        let Some(e) = row.evidence.as_ref().filter(|e| e.amplitude > 0.0) else {
            return NO_SHAPE;
        };
        let refs: Vec<&DonorFit> = self
            .fits
            .get(&row.precursor)
            .map(|v| v.iter().filter(|f| f.0 != row.run).collect())
            .unwrap_or_default();
        let mine = local_slope(maps, row.run, e.t0).zip(iso_profile(&e.iso_heights));
        let Some((slope, iso)) = mine.filter(|_| !refs.is_empty()) else {
            return NO_SHAPE;
        };
        let w_ref = median(refs.iter().map(|f| f.1).collect());
        let tail_ref = median(refs.iter().map(|f| f.2).collect());
        let mut iso_ref = [0.0; 4];
        for (i, v) in iso_ref.iter_mut().enumerate() {
            *v = median(refs.iter().map(|f| f.3[i]).collect());
        }
        [
            ((e.w / slope) / w_ref).ln(),
            e.tailing - tail_ref,
            cosine(&iso, &iso_ref),
        ]
    }
}

/// The PSM q-value a library entry's best donor PSM needs for the entry to
/// be transferred. Stricter than the donor threshold of the identified rows:
/// a single marginal donor PSM is the commonest false donor, and its
/// transfers sit on a real recurring feature that the transfer decoys
/// cannot tell from a true one.
pub(crate) const TRANSFER_DONOR_Q: f64 = 0.001;

/// The transfer FDR features of a row, in [`TRANSFER_FEATURES`] order, with
/// `neighbour_rate` the share of the run's most similar runs that identified
/// the precursor ([`crate::quant_finalize::EntryIndex::neighbour_rate`]);
/// `None` without evidence or with a non-finite value. The last three
/// compare the row's peak shape with its precursor's identified fits in the
/// other runs ([`DonorShapes::features`]): width, tailing and isotope
/// ratios, none of them a quantity. The MS2 evidence is
/// not a feature: MS2 support shows that the peptide was isolated, not that
/// the quantified MS1 peak is its peak, and the decoys cannot imitate a
/// transfer of the right peptide on the wrong peak (with MS2 features the
/// accepted wrong-peak rate rose from 1.0 % to 3.6 % on ProteoBench). It only
/// labels accepted transfers `transfer_ms2`.
pub(crate) fn transfer_features(
    row: &ExtractRow,
    shapes: &impl ShapeReference,
    maps: &[Option<RunMap>],
    neighbour_rate: f64,
) -> Option<Vec<f32>> {
    let e = row.evidence.as_ref()?;
    let snr = e.amplitude / e.amplitude_se.max(MIN_AMPLITUDE_SE);
    let [shape_width, shape_tailing, shape_iso] = shapes.features(row, maps);
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
        neighbour_rate,
        shape_width,
        shape_tailing,
        shape_iso,
    ];
    debug_assert_eq!(x.len(), TRANSFER_FEATURES.len());
    let x: Vec<f32> = x.iter().map(|&v| v as f32).collect();
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// The identified-row quality features of a row, in [`IDENTIFIED_FEATURES`]
/// order: the fit evidence only (no donor count, library q-value or other
/// run's measurement); `None` without evidence or with a non-finite value.
pub(crate) fn identified_features(row: &ExtractRow) -> Option<Vec<f32>> {
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
    ];
    debug_assert_eq!(x.len(), IDENTIFIED_FEATURES.len());
    let x: Vec<f32> = x.iter().map(|&v| v as f32).collect();
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// The decoy family of a hypothesis: 0 target, 1 mass shift, 2 RT shift.
pub(crate) fn family(kind: Hypothesis) -> u8 {
    match kind {
        Hypothesis::Target => 0,
        Hypothesis::MassShiftDecoy(_) => 1,
        Hypothesis::RtShiftDecoy => 2,
    }
}

/// `<base>.mbr.parquet` and `<base>.mbr_validation.parquet`.
pub(crate) fn mbr_paths(base: &Path) -> (PathBuf, PathBuf) {
    let with = |suffix: &str| {
        let mut s = base.as_os_str().to_owned();
        s.push(suffix);
        PathBuf::from(s)
    };
    (with(".mbr.parquet"), with(".mbr_validation.parquet"))
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
    use output::ms1_table::write_ms1_table;
    use quant::align::AlignParams;
    use quant::decoys::Hypothesis;
    use quant::library::{build_library, PrecursorMeta, PsmIn};
    use quant::ms1_index::Ms1RunIndex;

    use super::*;
    use crate::quant_finalize::finalize;

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
        run_scans_shaped(run, faims, silent, WIDTH, 1.0, 0.0)
    }

    /// Seconds between a precursor's apex and its echo in
    /// [`run_scans_shaped`].
    const ECHO_AFTER: f64 = 19.0;

    /// [`run_scans`] with elution width `width`, one scan every `step`
    /// seconds, and, when `echo > 0`, a second peak of every precursor's
    /// pattern [`ECHO_AFTER`] seconds after its apex at `echo` times its
    /// amplitude (another species at the same m/z).
    fn run_scans_shaped(
        run: u16,
        faims: bool,
        silent: &dyn Fn(u32) -> bool,
        width: f64,
        step: f64,
        echo: f64,
    ) -> Vec<Ms1Scan> {
        let mut rng = Lcg(17 + run as u64);
        (0..(1300.0 / step) as usize)
            .map(|k| {
                let t = k as f64 * step;
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
                    let profile = |c: f64| (-(t - c).powi(2) / (2.0 * width * width)).exp();
                    let g = profile(apex(p, run)) + echo * profile(apex(p, run) + ECHO_AFTER);
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
        ms1_tables: Vec<PathBuf>,
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
        fixture_decoys(faims, silent_in_run1, weak_in_run1, &|_| false)
    }

    /// [`fixture_with`] where the precursors in `decoy` are decoy
    /// precursors: their PSMs are decoy PSMs, so they are identity decoys.
    fn fixture_decoys(
        faims: bool,
        silent_in_run1: &dyn Fn(u32) -> bool,
        weak_in_run1: &dyn Fn(u32) -> bool,
        decoy: &dyn Fn(u32) -> bool,
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
                is_decoy: decoy(p),
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
        let ms1_tables = (0..N_RUNS as u16)
            .map(|r| {
                let path = dir.path().join(format!("run{r}.ms1.parquet"));
                let silent = |p: u32| r == 1 && silent_in_run1(p);
                let index = Ms1RunIndex::new(run_scans(r, faims && r == 1, &silent));
                write_ms1_table(&path, &format!("run{r}"), &index).unwrap();
                path
            })
            .collect();
        Fixture {
            _dir: dir,
            lib,
            ms1_tables,
        }
    }

    fn params() -> EvidenceParams {
        EvidenceParams::default()
    }

    /// [`extract_run`] measuring every transfer pair.
    fn extract_all(
        run: u16,
        ms1_table: &Path,
        lib: &Library,
        ep: &EvidenceParams,
    ) -> Result<RunExtract, String> {
        extract_run(run, ms1_table, lib, ep, &TransferSelection::OFF)
    }

    fn rows_of(rows: &[ExtractRow], p: u32) -> Vec<&ExtractRow> {
        rows.iter().filter(|r| r.precursor == p).collect()
    }

    fn is_identified_target(r: &ExtractRow) -> bool {
        (r.role, r.kind) == (Role::Identified, Hypothesis::Target)
    }

    #[test]
    fn identified_pairs_get_identified_and_pseudo_rows() {
        let fx = fixture(false);
        assert!(fx.lib.alignment.maps[1].is_some());
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        for p in 0..BOTH {
            let rows = rows_of(&out.rows, p);
            assert_eq!(rows.len(), 4, "precursor {p}: {rows:?}");
            let (id, pseudo) = (rows[0], rows[3]);
            assert_eq!((id.role, id.kind), (Role::Identified, Hypothesis::Target));
            // The identified target's decoys follow it, with its role.
            assert_eq!(rows[1].role, Role::Identified);
            assert!(matches!(rows[1].kind, Hypothesis::MassShiftDecoy(_)));
            assert_eq!(
                (rows[2].role, rows[2].kind),
                (Role::Identified, Hypothesis::RtShiftDecoy)
            );
            for d in &rows[1..3] {
                assert_eq!(
                    (d.n_donor_runs, d.library_q),
                    (id.n_donor_runs, id.library_q)
                );
                let e = d.evidence.as_ref().unwrap();
                let target_llr = id.evidence.as_ref().unwrap().llr;
                assert!(e.llr < target_llr / 10.0, "precursor {p}: {d:?}");
            }
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
        assert_eq!(rows.len(), 3);
        assert!(rows.iter().all(|r| r.role == Role::Identified));
        assert_eq!(out.counts.no_pseudo_prediction, 1);
        assert!(!out.counts.unaligned);
    }

    #[test]
    fn identified_decoys_are_the_transfer_probes_at_the_identified_prior() {
        let fx = fixture(false);
        let ms1 = Ms1RunIndex::new(run_scans(1, false, &|_| false));
        let scales = RunScales {
            spacing_s: ms1.median_scan_spacing(),
            width_s: WIDTH,
            ..RunScales::default()
        };
        let wp = WidthPrior::from_fits(Vec::new(), &scales);
        let entry = fx.lib.entries.iter().find(|e| e.precursor == 4).unwrap();
        let off = PairGate::new(&TransferSelection::OFF, &fx.lib);
        let (rows, _) = entry_rows(entry, 1, true, &fx.lib, &off, &ms1, &wp);
        let id: Vec<&(ExtractRow, Probe)> = rows
            .iter()
            .filter(|(r, _)| r.role == Role::Identified)
            .collect();
        assert_eq!(id.len(), 3, "{rows:?}");
        let rt = identified_rt(entry, 1).unwrap();
        let w = wp.at(rt);
        // Exactly the probes a transfer gets, with the identified prior.
        let want = probes(
            4,
            1,
            mz_of(4),
            CHARGE,
            rt,
            w.max(1.0),
            w,
            ms1.rt_range().unwrap(),
        );
        let got: Vec<Probe> = id.iter().map(|(_, p)| *p).collect();
        assert_eq!(got, want);
        let target = got[0];
        assert_eq!(target.kind, Hypothesis::Target);
        assert_eq!((target.t_hat, target.mz), (rt, mz_of(4)));
        for (row, probe) in &id {
            assert_eq!(row.kind, probe.kind);
            // Same search breadth for target and decoys.
            assert_eq!((probe.s, probe.w_prior), (target.s, target.w_prior));
        }
        let mass = got[1];
        let Hypothesis::MassShiftDecoy(n) = mass.kind else {
            panic!("{mass:?}")
        };
        assert_eq!(mass.t_hat, target.t_hat);
        let shift = n as f64 * quant::decoys::PEPTIDE_MASS_SLOPE / CHARGE as f64;
        assert!((mass.mz - target.mz - shift).abs() < 1e-9);
        let rt_decoy = got[2];
        assert_eq!(rt_decoy.mz, target.mz);
        assert!((rt_decoy.t_hat - rt).abs() >= 8.0 * target.s + 4.0 * target.w_prior);
    }

    #[test]
    fn identity_decoys_are_measured_as_targets() {
        // Precursor ONLY_RUN0 as a decoy precursor: its decoy PSM in run 0
        // makes it an identity decoy on the same real MS1 feature.
        let target = fixture(false);
        let decoy = fixture_decoys(false, &|_| false, &|_| false, &|p| p == ONLY_RUN0);
        let e = decoy
            .lib
            .entries
            .iter()
            .find(|e| e.precursor == ONLY_RUN0)
            .unwrap();
        assert!(e.identity_decoy);
        assert_eq!(e.donors, vec![0]);
        // The decoy moves every entry gate q-value; nothing else differs.
        let same = |rows: Vec<&ExtractRow>| -> Vec<ExtractRow> {
            rows.into_iter()
                .map(|r| ExtractRow {
                    library_q: 0.0,
                    ..r.clone()
                })
                .collect()
        };
        for run in [0u16, 1] {
            let a = extract_all(
                run,
                &target.ms1_tables[run as usize],
                &target.lib,
                &params(),
            )
            .unwrap();
            let b =
                extract_all(run, &decoy.ms1_tables[run as usize], &decoy.lib, &params()).unwrap();
            let rows = rows_of(&b.rows, ONLY_RUN0);
            // Identified in run 0 (target and decoy probes), a transfer with
            // its decoys in run 1: the target entry's rows, fitted alike.
            assert_eq!(rows.len(), 3, "run {run}: {rows:?}");
            let role = if run == 0 {
                Role::Identified
            } else {
                Role::Transfer
            };
            assert!(rows.iter().all(|r| r.role == role));
            if run == 1 {
                let ev = rows[0].evidence.as_ref().unwrap();
                assert!(ev.llr > 100.0, "{ev:?}");
                assert_eq!(same(rows_of(&a.rows, ONLY_RUN0)), same(rows));
                assert_eq!(same(a.rows.iter().collect()), same(b.rows.iter().collect()));
            } else {
                // The width prior of run 0 leaves the identity decoy out.
                assert_eq!(a.counts.width_fits, b.counts.width_fits + 1);
            }
        }
    }

    #[test]
    fn missing_pairs_get_target_and_decoys() {
        let fx = fixture(false);
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
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
        let out0 = extract_all(0, &fx.ms1_tables[0], &fx.lib, &params()).unwrap();
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

    /// `lib` with the entry of precursor `p` changed by `f`.
    fn with_entry(lib: &Library, p: u32, f: impl Fn(&mut LibraryEntry)) -> Library {
        let mut lib = lib.clone();
        f(lib.entries.iter_mut().find(|e| e.precursor == p).unwrap());
        lib
    }

    /// The rows of `all` without the transfer rows of the pairs in `skipped`.
    fn without(all: &[ExtractRow], skipped: &[(u32, u16)]) -> Vec<ExtractRow> {
        all.iter()
            .filter(|r| !(r.role == Role::Transfer && skipped.contains(&(r.precursor, r.run))))
            .cloned()
            .collect()
    }

    const SELECT: TransferSelection = TransferSelection {
        donor_gate: true,
        min_donor_fraction: 0.5,
    };

    #[test]
    fn selection_skips_transfers_finalize_cannot_accept() {
        let fx = fixture(false);
        // The best donor PSM of precursor ONLY_RUN0 misses the transfer donor
        // threshold: its transfer into run 1 (target and decoys) is not
        // measured; everything else is measured as without a selection.
        let lib = with_entry(&fx.lib, ONLY_RUN0, |e| e.donor_psm.q = 0.002);
        let all = extract_all(1, &fx.ms1_tables[1], &lib, &params()).unwrap();
        let sel = extract_run(1, &fx.ms1_tables[1], &lib, &params(), &SELECT).unwrap();
        assert_eq!(rows_of(&all.rows, ONLY_RUN0).len(), 3);
        assert!(rows_of(&sel.rows, ONLY_RUN0).is_empty());
        assert_eq!(sel.rows, without(&all.rows, &[(ONLY_RUN0, 1)]));
        assert_eq!((all.counts.not_selected, sel.counts.not_selected), (0, 1));
        // Without the donor gate the pair is measured.
        let no_gate = TransferSelection {
            donor_gate: false,
            ..SELECT
        };
        let kept = extract_run(1, &fx.ms1_tables[1], &lib, &params(), &no_gate).unwrap();
        assert_eq!(kept.rows, all.rows);
    }

    #[test]
    fn selection_skips_transfers_with_too_few_donor_runs() {
        let fx = fixture(false);
        // Precursor ONLY_RUN1 keeps its alignment observation in run 1 (so
        // run 0 still has a prediction) but has no donor run left.
        let lib = with_entry(&fx.lib, ONLY_RUN1, |e| e.donors.clear());
        let all = extract_all(0, &fx.ms1_tables[0], &lib, &params()).unwrap();
        assert_eq!(rows_of(&all.rows, ONLY_RUN1).len(), 3);
        let sel = extract_run(0, &fx.ms1_tables[0], &lib, &params(), &SELECT).unwrap();
        assert_eq!(sel.rows, without(&all.rows, &[(ONLY_RUN1, 0)]));
        assert_eq!(sel.counts.not_selected, 1);
        // The identified rows and pseudo-transfers of every precursor stay.
        let roles =
            |rows: &[ExtractRow], role: Role| rows.iter().filter(|r| r.role == role).count();
        for role in [Role::Identified, Role::PseudoTransfer] {
            assert_eq!(roles(&sel.rows, role), roles(&all.rows, role));
        }
        // The fixture's own transfers have one donor run each: kept.
        let sel = extract_run(0, &fx.ms1_tables[0], &fx.lib, &params(), &SELECT).unwrap();
        assert_eq!(
            sel,
            extract_all(0, &fx.ms1_tables[0], &fx.lib, &params()).unwrap()
        );
    }

    #[test]
    fn min_donors_scale_with_the_aligned_runs() {
        let sel = |f: f64| TransferSelection {
            donor_gate: false,
            min_donor_fraction: f,
        };
        assert_eq!(sel(0.0).min_donors(500), 0);
        assert_eq!(sel(0.004).min_donors(6), 1);
        assert_eq!(sel(0.004).min_donors(250), 1);
        assert_eq!(sel(0.004).min_donors(500), 2);
        assert_eq!(sel(0.004).min_donors(1000), 4);
        assert_eq!(sel(1.0).min_donors(1), 0);
        assert!(TransferSelection::OFF.is_off());
        assert!(!TransferSelection::default().is_off());
    }

    #[test]
    fn transfers_stay_inside_an_alignment_component() {
        // Runs 0/1 share precursors 0..60, runs 2/3 share 60..120; precursor
        // 120 is identified in run 0 only, 121 in runs 0 and 2.
        let dir = tempfile::tempdir().unwrap();
        let n = 122u32;
        let meta: Vec<PrecursorMeta> = (0..n)
            .map(|p| PrecursorMeta {
                key: format!("PEP{p}/2"),
                mz: 400.0 + 7.3 * p as f64,
                charge: CHARGE,
                envelope: ENVELOPE,
            })
            .collect();
        let mut psms = Vec::new();
        let mut psm = |p: u32, run: u16| {
            psms.push(PsmIn {
                run,
                precursor: p,
                rt: 100.0 + 8.0 * p as f64 + if run >= 2 { 120.0 } else { 0.0 },
                score: 10.0,
                q: 0.001,
                pep: 0.0,
                is_decoy: false,
                faims_cv: None,
            })
        };
        for p in 0..120 {
            let pair = if p < 60 { [0, 1] } else { [2, 3] };
            pair.into_iter().for_each(|r| psm(p, r));
        }
        psm(120, 0);
        psm(121, 0);
        psm(121, 2);
        let ap = AlignParams {
            min_anchors: 10,
            ..AlignParams::default()
        };
        let lib = build_library(&psms, meta, 4, 1.0, DONOR_Q, None, &ap);
        assert_eq!(lib.alignment.n_components(), 2);
        let path = dir.path().join("run3.ms1.parquet");
        let index = Ms1RunIndex::new(run_scans(3, false, &|_| false));
        write_ms1_table(&path, "run3", &index).unwrap();
        let out = extract_all(3, &path, &lib, &params()).unwrap();
        // 0..60 and 120 come from the other component: no rows, counted.
        assert_eq!(out.counts.other_component, 61);
        assert_eq!(out.counts.no_prediction, 0);
        for p in (0..60).chain([120]) {
            assert!(rows_of(&out.rows, p).is_empty(), "{p}");
        }
        // 121 is transferred from run 2, its donor in this component.
        let rows = rows_of(&out.rows, 121);
        assert!(!rows.is_empty());
        assert!(rows.iter().all(|r| r.role == Role::Transfer), "{rows:?}");
        assert_eq!(rows[0].kind, Hypothesis::Target);
        // Run 2 identifies 121, but its only other donor is in the other
        // component: no pseudo-transfer.
        let path = dir.path().join("run2.ms1.parquet");
        let index = Ms1RunIndex::new(run_scans(2, false, &|_| false));
        write_ms1_table(&path, "run2", &index).unwrap();
        let out = extract_all(2, &path, &lib, &params()).unwrap();
        assert_eq!(out.counts.pseudo_other_component, 1);
        let rows = rows_of(&out.rows, 121);
        assert!(rows.iter().all(|r| r.role == Role::Identified), "{rows:?}");
        assert_eq!(rows.iter().filter(|r| is_identified_target(r)).count(), 1);
    }

    #[test]
    fn unaligned_run_keeps_its_identified_rows() {
        let fx = fixture(false);
        assert!(fx.lib.alignment.maps[2].is_none());
        let out = extract_all(2, &fx.ms1_tables[2], &fx.lib, &params()).unwrap();
        assert!(out.counts.unaligned);
        // Run 2 identifies precursor 0 only: its identified row and decoys,
        // and no pseudo-transfer or transfer rows.
        assert_eq!(out.rows.len(), 3, "{:?}", out.rows);
        assert!(out.rows.iter().all(|r| r.role == Role::Identified));
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
        assert!(extract_all(N_RUNS as u16, &fx.ms1_tables[2], &fx.lib, &params()).is_err());

        // The identified row reaches the mbr.parquet.
        let runs: Vec<Vec<ExtractRow>> = (0..N_RUNS)
            .map(|r| {
                extract_all(r as u16, &fx.ms1_tables[r], &fx.lib, &params())
                    .unwrap()
                    .rows
            })
            .collect();
        let stems: Vec<String> = ["r0", "r1", "r2"].map(String::from).to_vec();
        let (mbr, val) = mbr_paths(&fx._dir.path().join("x"));
        finalize(&stems, &fx.lib, &runs, 0.01, None, 7, &mbr, &val).unwrap();
        let rows = output::mbr_out::read_mbr_parquet(&mbr).unwrap();
        let in_r2: Vec<_> = rows.iter().filter(|r| r.run == "r2").collect();
        assert_eq!(in_r2.len(), 1, "{rows:?}");
        assert_eq!(
            (
                in_r2[0].precursor_key.as_str(),
                in_r2[0].evidence_type.as_str()
            ),
            ("PEP0/2", "identified"),
            "{in_r2:?}"
        );
        // Its quantity and standard error are reported.
        assert!(
            in_r2[0].quantity.is_some() && in_r2[0].quantity_se.is_some(),
            "{in_r2:?}"
        );
    }

    #[test]
    fn psms_above_the_donor_threshold_make_a_transfer() {
        // Precursor 7's only PSM in run 1 is above the PSM q-value threshold.
        let fx = fixture_with(false, &|_| false, &|p| p == 7);
        let e = fx.lib.entries.iter().find(|e| e.precursor == 7).unwrap();
        assert_eq!(e.donors, vec![0]);
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
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
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        assert!(rows_of(&out.rows, ONLY_RUN0).is_empty());
        assert!(rows_of(&out.rows, 3).is_empty());
        assert_eq!(out.counts.foreign_cv, 2);
        assert!(!rows_of(&out.rows, 4).is_empty());
        // A run without FAIMS fits them on the whole run.
        let out0 = extract_all(0, &fx.ms1_tables[0], &fx.lib, &params()).unwrap();
        assert_eq!(out0.counts.foreign_cv, 0);
        assert_eq!(rows_of(&out0.rows, ONLY_RUN0)[0].role, Role::Identified);
        assert_eq!(rows_of(&out0.rows, 3)[0].role, Role::Identified);
    }

    #[test]
    fn output_order_is_deterministic() {
        let fx = fixture(false);
        let a = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        let b = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
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
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        assert!(rows_of(&out.rows, 5).is_empty());
        assert_eq!(out.counts.nonfinite_best_rt, 1);
    }

    #[test]
    fn width_prior_follows_identified_fits() {
        let fx = fixture(false);
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        assert!(out.counts.width_fits >= 20, "{:?}", out.counts);
        assert!(!out.counts.width_prior_fallback);
        for r in out.rows.iter().filter(|r| is_identified_target(r)) {
            let e = r.evidence.as_ref().unwrap();
            let w_prior = e.w / e.width_ratio;
            assert!(
                (w_prior / WIDTH - 1.0).abs() < 0.1,
                "w_prior {w_prior}: {e:?}"
            );
        }
    }

    #[test]
    fn run_scales_come_from_the_run() {
        let fx = fixture(false);
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        let sc = out.counts.scales;
        assert_eq!(sc.spacing_s, Some(1.0));
        assert!(!sc.fallback, "{sc:?}");
        assert!(sc.width_fits >= 20, "{sc:?}");
        assert!((sc.width_s / WIDTH - 1.0).abs() < 0.1, "{sc:?}");
        assert!(
            sc.describe().contains("ms1_spacing_s=1.0000;"),
            "{}",
            sc.describe()
        );
    }

    /// Run 1 of [`fixture`] with narrow peaks sampled every 0.25 s, and a
    /// same-m/z peak five times larger [`ECHO_AFTER`] seconds after each
    /// apex, written next to the fixture's tables.
    fn narrow_run(fx: &Fixture, width: f64) -> PathBuf {
        let path = fx._dir.path().join("narrow.ms1.parquet");
        let index = Ms1RunIndex::new(run_scans_shaped(1, false, &|_| false, width, 0.25, 5.0));
        write_ms1_table(&path, "run1", &index).unwrap();
        path
    }

    #[test]
    fn narrow_peaks_set_a_narrow_width_scale() {
        // The rounds start at 6 s; the run's peaks are 1 s wide.
        let fx = fixture(false);
        let path = narrow_run(&fx, 1.0);
        let out = extract_all(1, &path, &fx.lib, &params()).unwrap();
        let sc = out.counts.scales;
        assert_eq!(sc.spacing_s, Some(0.25));
        assert!(!sc.fallback && sc.rounds > 1, "{sc:?}");
        assert!((sc.width_s - 1.0).abs() < 0.1, "{sc:?}");
        for r in out.rows.iter().filter(|r| is_identified_target(r)) {
            let e = r.evidence.as_ref().unwrap();
            let w_prior = e.w / e.width_ratio;
            assert!(
                (w_prior - 1.0).abs() < 0.15,
                "precursor {}: {e:?}",
                r.precursor
            );
        }
    }

    #[test]
    fn identified_fit_stays_on_its_own_peak() {
        // The PSM is 2 s after a 1 s wide apex; a peak five times larger
        // elutes 17 s after the PSM, inside a fixed 18 s search. The
        // identified search spans three run widths around the PSM, so the
        // fit stays on the identified peak.
        let fx = fixture(false);
        let path = narrow_run(&fx, 1.0);
        let out = extract_all(1, &path, &fx.lib, &params()).unwrap();
        let mut n = 0;
        for r in out.rows.iter().filter(|r| is_identified_target(r)) {
            let e = r.evidence.as_ref().unwrap();
            assert!(
                (e.t0 - apex(r.precursor, 1)).abs() <= 0.5,
                "precursor {}: {e:?}",
                r.precursor
            );
            n += 1;
        }
        assert!(n >= 30, "{n}");
    }

    #[test]
    fn too_few_identified_fits_fall_back_to_the_seed_width() {
        // Every identified pair of run 1 is silent: no usable fit.
        let fx = fixture_with(false, &|_| true, &|_| false);
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        let sc = out.counts.scales;
        assert!(sc.fallback, "{sc:?}");
        assert_eq!((sc.rounds, sc.width_s), (1, SEED_WIDTH_S));
        assert!(out.counts.width_prior_fallback);
    }

    #[test]
    fn faims_entry_is_fitted_on_its_voltage_series() {
        let fx = fixture(true);
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        for r in out.rows.iter().filter(|r| is_identified_target(r)) {
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
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
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
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        let (peak, loads) = tracking::stats();
        assert_eq!(peak, 1, "at most one MS1 copy is live");
        assert_eq!(loads, 2, "the whole run is read again for its entries");
        // The whole-run entry is still measured, in library order.
        assert_eq!(rows_of(&out.rows, ONLY_RUN0).len(), 3);
        let mut seen: Vec<u32> = out.rows.iter().map(|r| r.precursor).collect();
        seen.dedup();
        let order: Vec<u32> = fx.lib.entries.iter().map(|e| e.precursor).collect();
        assert_eq!(seen, order);

        // Without FAIMS the MS1 table is read once.
        let fx = fixture(false);
        tracking::reset();
        extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        assert_eq!(tracking::stats(), (1, 1));
    }

    #[test]
    fn width_prior_ignores_identified_pairs_without_signal() {
        // About 30% of run 1's identified pairs have no MS1 signal.
        let silent = |p: u32| p % 10 < 3;
        let fx = fixture_with(false, &silent, &|_| false);
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        assert!(out.counts.width_fits >= 20, "{:?}", out.counts);
        assert!(!out.counts.width_prior_fallback);
        for r in out.rows.iter().filter(|r| is_identified_target(r)) {
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
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
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
    use output::mbr_out::{
        write_mbr_validation, Accepted, EvidenceType, MbrParquetWriter, MbrRow, MbrValidation,
        VALIDATION_Q,
    };
    use quant::align::AlignParams;
    use quant::library::{build_library, PrecursorMeta, PsmIn};
    use quant::normalize::{fit as fit_normalization, NormCell, NormParams};

    use super::*;
    use crate::quant_finalize::{
        finalize, finalize_stream, for_each_scored, table_features, EntryIndex, ExtractRuns,
        FinalizeLimits, FinalizeParams, FinalizeSummary, MbrFeature, Table, SCORE_CHUNK,
    };
    use crate::transfer_fdr::{
        row_unit, train_folds, Estimate, FdrRow, FoldScores, QValues, TrainingSample,
    };

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
            iso_heights: [0.55, 0.3, 0.11, 0.04].map(|v| v * amplitude),
            tailing: 0.0,
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
            ms2: None,
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

    /// A header row and one row of cells per table row, numbers formatted
    /// with `{}` and a null as `NA`, for `<base>.mbr.parquet` or
    /// `<base>.mbr_validation.parquet` by file name.
    fn table(path: &Path) -> Vec<Vec<String>> {
        use output::mbr_out::{read_mbr_parquet, read_mbr_validation, ValidationCell, MBR_COLUMNS};
        let num = |v: Option<f64>| v.map_or_else(|| "NA".to_string(), |v| v.to_string());
        if path.to_string_lossy().ends_with(".mbr_validation.parquet") {
            let cells = read_mbr_validation(path).unwrap();
            let header = cells.iter().map(|(n, _)| n.clone()).collect();
            let values = cells
                .iter()
                .map(|(_, c)| match *c {
                    ValidationCell::Count(n) => n.to_string(),
                    ValidationCell::Value(v) => num(v),
                })
                .collect();
            return vec![header, values];
        }
        let mut t = vec![MBR_COLUMNS.iter().map(|c| c.to_string()).collect()];
        for r in read_mbr_parquet(path).unwrap() {
            t.push(vec![
                r.precursor_key,
                r.run,
                r.evidence_type,
                num(r.quantity),
                num(r.quantity_se),
                num(r.quantity_normalized),
                num(r.transfer_q),
                num(r.quant_q),
                num(r.apex_rt_s),
                num(Some(r.library_q)),
            ]);
        }
        t
    }

    /// [`library`] with every precursor an entry identified in run 0; those
    /// in `identity` are identity decoys. Precursor `p`'s best donor PSM has
    /// q-value `donor_q(p)`.
    fn library_with(identity: &dyn Fn(u32) -> bool, donor_q: &dyn Fn(u32) -> f64) -> Library {
        let mut lib = library();
        lib.entries = (0..N)
            .map(|p| LibraryEntry {
                precursor: p,
                identity_decoy: identity(p),
                q: 0.001,
                donor_psm: quant::library::DonorPsm {
                    score: 10.0,
                    q: donor_q(p),
                    pep: 1e-4,
                },
                consensus_rt: f64::NAN,
                faims_cv: None,
                donors: vec![0],
                best_rt: vec![(0, 100.0 + p as f64)],
            })
            .collect();
        lib
    }

    /// Every 50th precursor from 2 on: 18 identity decoys, all on a real
    /// feature in run 1 (even precursors are present there), as the
    /// transfer of a false donor onto the feature it was misidentified on.
    fn on_real_features(p: u32) -> bool {
        p % 50 == 2
    }

    fn confident(_: u32) -> f64 {
        1e-4
    }

    #[test]
    fn identity_decoys_on_real_features_count_in_the_transfer_q() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let lib = library_with(&on_real_features, &confident);
        let s = finalize(&stems(), &lib, &extract(), 0.2, None, 7, &mbr, &val).unwrap();
        let v = &s.validation;
        eprintln!(
            "q 0.01 pooled {} wrong-peak only {}; identity decoys accepted at 0.2: {} of {}",
            v.at_q01.pooled,
            v.at_q01.wrong_peak_only,
            v.at_threshold.decoys[2],
            v.decoy_transfers[2]
        );
        assert_eq!((v.target_entries, v.identity_entries), (882, 18));
        assert_eq!(v.decoy_transfers[2], 18);
        // On real features they look like true transfers and pass ...
        assert!(v.at_threshold.decoys[2] >= 12, "{v:?}");
        // ... so the transfer q-value counts about 4 % false donors that
        // the wrong-peak decoys alone do not see.
        assert!(v.at_q01.wrong_peak_only >= 380, "{v:?}");
        assert!(2 * v.at_q01.pooled < v.at_q01.wrong_peak_only, "{v:?}");
        assert!(v.at_threshold.pooled >= 380, "{v:?}");
        // Without identity decoys the wrong-peak estimate is all there is.
        let plain = library_with(&|_| false, &confident);
        let s = finalize(&stems(), &plain, &extract(), 0.2, None, 7, &mbr, &val).unwrap();
        let v = &s.validation;
        assert_eq!(v.decoy_transfers[2], 0);
        assert_eq!(v.at_q01.pooled, v.at_q01.wrong_peak_only);
    }

    #[test]
    fn identity_decoy_rows_never_reach_the_outputs() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let lib = library_with(&on_real_features, &confident);
        let s = finalize(&stems(), &lib, &extract(), 0.5, None, 7, &mbr, &val).unwrap();
        let t = table(&mbr);
        let keys: Vec<&str> = t[1..].iter().map(|r| r[0].as_str()).collect();
        assert_eq!(keys.len(), 2 * (N as usize - 18));
        for p in (0..N).filter(|&p| on_real_features(p)) {
            let key = format!("PEP{p}/2");
            assert!(!keys.contains(&key.as_str()), "{key}");
        }
        assert_eq!(s.rows_written, keys.len());
        assert!(s.features.iter().all(|f| !on_real_features(f.precursor)));
        // The identified quality model and the normalisation never see them
        // either: every reported identified row is a target's.
        assert_eq!(s.validation.identified, N as usize - 18);
    }

    #[test]
    fn entries_with_a_marginal_best_donor_psm_are_not_transferred() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        // Precursors above 400 have a best donor PSM at q 0.005.
        let lib = library_with(&|_| false, &|p| {
            if p > 400 {
                0.005
            } else {
                TRANSFER_DONOR_Q
            }
        });
        let s = finalize(&stems(), &lib, &extract(), 0.01, None, 7, &mbr, &val).unwrap();
        assert_eq!(s.validation.donor_gated, N as usize - 401);
        let mut accepted_low = 0;
        for f in s.features.iter().filter(|f| f.run == 1) {
            if f.precursor > 400 {
                assert_eq!(f.evidence_type, EvidenceType::BelowThreshold);
                assert_eq!(f.transfer_q, None);
            } else {
                accepted_low += usize::from(f.evidence_type != EvidenceType::BelowThreshold);
            }
        }
        assert!(accepted_low >= 150, "{accepted_low}");
        // Their identified rows are reported as before.
        assert_eq!(s.validation.identified, N as usize);
    }

    #[test]
    fn finalize_writes_one_row_per_target_pair_with_evidence_types() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let s = finalize(&stems(), &library(), &extract(), 0.01, None, 7, &mbr, &val).unwrap();
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
        let s = finalize(&stems(), &library(), &extract(), 0.01, None, 7, &mbr, &val).unwrap();
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
        let loose = finalize(&stems(), &library(), &extract(), 0.5, None, 7, &mbr, &val).unwrap();
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
        finalize(&stems(), &library(), &runs, 0.01, None, 3, &m1, &v1).unwrap();
        finalize(&stems(), &library(), &runs, 0.01, None, 3, &m2, &v2).unwrap();
        assert_eq!(std::fs::read(&m1).unwrap(), std::fs::read(&m2).unwrap());
        assert_eq!(std::fs::read(&v1).unwrap(), std::fs::read(&v2).unwrap());
    }

    /// [`library`]'s precursors as entries, each identified in run 0.
    fn library_with_entries() -> Library {
        let psms: Vec<PsmIn> = (0..N)
            .map(|p| PsmIn {
                run: 0,
                precursor: p,
                rt: 100.0 + p as f64,
                score: 10.0 + p as f64,
                q: 0.0001,
                pep: 0.0001,
                is_decoy: false,
                faims_cv: None,
            })
            .collect();
        build_library(
            &psms,
            library().meta,
            2,
            1.0,
            0.01,
            None,
            &AlignParams::default(),
        )
    }

    /// The rows of every run as one shard, for [`reference`].
    struct One<'a>(&'a [ExtractRow]);

    impl One<'_> {
        fn for_each(
            &self,
            mut f: impl FnMut(&[ExtractRow]) -> Result<(), String>,
        ) -> Result<(), String> {
            f(self.0)
        }
    }

    /// Finalize over every row in memory, one pass per stage: the reference
    /// of the map-reduce [`finalize_stream`].
    #[allow(clippy::too_many_arguments)]
    fn reference(
        stems: &[String],
        lib: &Library,
        runs: &[Vec<ExtractRow>],
        mbr_fdr: f64,
        quant_feature_fdr: Option<f64>,
        seed: u64,
        mbr_path: &Path,
        validation_path: &Path,
    ) -> Result<FinalizeSummary, String> {
        let mut all_features = Vec::new();
        let on_features = &mut |f: Vec<MbrFeature>| -> Result<(), String> {
            all_features.extend(f);
            Ok(())
        };
        if !(mbr_fdr > 0.0 && mbr_fdr <= 1.0) {
            return Err(format!(
                "transfer q-value threshold {mbr_fdr} is outside (0, 1]"
            ));
        }
        if let Some(t) = quant_feature_fdr.filter(|t| !(*t > 0.0 && *t <= 1.0)) {
            return Err(format!(
                "identified quality q-value threshold {t} is outside (0, 1]"
            ));
        }
        let maps = &lib.alignment.maps;
        let idx = EntryIndex::new(lib, stems.len());
        let all: Vec<ExtractRow> = runs.concat();
        let shards = One(&all);
        let mut summary = FinalizeSummary::default();

        // Pass 1: the training rows of both classifiers. Spilled shards hold a
        // precursor's rows of every run, so the source order (run, then library
        // order, then the order within an entry) is restored by key.
        type Key = (u16, u32, u64);
        let mut t_sample: TrainingSample<Key> = TrainingSample::new(usize::MAX, seed);
        let mut i_sample: TrainingSample<Key> = TrainingSample::new(usize::MAX, seed);
        let mut seq = 0u64;
        shards.for_each(|rows| {
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
            let shapes = DonorShapes::from_rows(rows.iter(), maps);
            for (c, chunk) in rows.chunks(SCORE_CHUNK).enumerate() {
                let feats = table_features(chunk, &idx, &shapes, maps, |t: Table| t != Table::Pseudo);
                for (i, t, fam, x) in feats {
                    seq += 1;
                    let r = &rows[c * SCORE_CHUNK + i];
                    let Some(x) = x else {
                        match t {
                            Table::Identified => summary.identified_no_features += 1,
                            _ => summary.no_features += 1,
                        }
                        continue;
                    };
                    // Source order (run, then row).
                    let key = (0u16, 0u32, seq);
                    let row = FdrRow {
                        precursor: r.precursor,
                        is_decoy: fam != 0,
                        family: fam,
                        x,
                    };
                    match t {
                        Table::Identified => i_sample.push(key, row_unit(r.precursor, r.run), row),
                        _ => t_sample.push(key, row_unit(r.precursor, r.run), row),
                    }
                }
            }
            Ok(())
        })?;
        let t_models = train_folds(&t_sample.into_items(), seed);
        let i_models = train_folds(&i_sample.into_items(), seed);

        // Pass 2: every row's held-out score, into the q-value step functions.
        let mut t_scores = FoldScores::new();
        let mut i_scores = FoldScores::new();
        shards.for_each(|rows| {
            let shapes = DonorShapes::from_rows(rows.iter(), maps);
            let want = |t: Table| t != Table::Pseudo;
            for_each_scored(
                rows,
                &idx,
                &shapes,
                maps,
                &t_models,
                &i_models,
                want,
                |_, t, fam, s| match t {
                    Table::Identified => i_scores.add(fam, s),
                    _ => t_scores.add(fam, s),
                },
            );
            Ok(())
        })?;
        let fdr = t_scores.finish(&t_models.failed_folds);
        let id_fdr = i_scores.finish(&i_models.failed_folds);

        // Pass 3: the reference cells of the run normalisation.
        let mut cells: TrainingSample<u64, NormCell> = TrainingSample::new(usize::MAX, seed);
        let mut cell_seq = 0u64;
        let empty_shapes = DonorShapes::default();
        shards.for_each(|rows| {
            let w_c = consensus_widths(rows.iter(), maps);
            let want = |t: Table| t == Table::Identified;
            for_each_scored(
                rows,
                &idx,
                &empty_shapes,
                maps,
                &t_models,
                &i_models,
                want,
                |i, _, fam, s| {
                    let r = &rows[i];
                    if r.kind != Hypothesis::Target
                        || !id_fdr
                            .q(Estimate::Pooled, fam, s)
                            .is_some_and(|q| q <= VALIDATION_Q)
                    {
                        return;
                    }
                    let Some(e) = r.evidence.as_ref() else { return };
                    let q = reported_quantity(e, r.run, w_c.get(&r.precursor).copied(), maps).0;
                    if q > 0.0 && q.is_finite() {
                        cell_seq += 1;
                        cells.push(
                            cell_seq,
                            u64::from(r.precursor),
                            NormCell {
                                precursor: r.precursor,
                                run: r.run,
                                log2_quantity: q.log2(),
                                rt: e.t0,
                            },
                        );
                    }
                },
            );
            Ok(())
        })?;
        let norm = fit_normalization(&cells.into_items(), stems.len(), &NormParams::default());

        // Pass 4: the outputs.
        let accepted = |q: Option<f64>| q.is_some_and(|q| q <= mbr_fdr);
        // The validation thresholds, VALIDATION_Q then mbr_fdr, and their counts.
        let thresholds = [VALIDATION_Q, mbr_fdr];
        let mut at = [Accepted::default(), Accepted::default()];
        let identity_entries = lib.entries.iter().filter(|e| e.identity_decoy).count();
        let mut v = MbrValidation {
            q_threshold: mbr_fdr,
            quant_feature_threshold: quant_feature_fdr,
            target_entries: lib.entries.len() - identity_entries,
            identity_entries,
            ..MbrValidation::default()
        };
        // Decoy transfer rows with a q-value, and their acceptance, by family:
        // mass shift, RT shift, identity decoy.
        let count_decoy =
            |v: &mut MbrValidation, at: &mut [Accepted; 2], q: Option<f64>, f: usize| {
                let Some(q) = q else { return };
                v.decoy_transfers[f] += 1;
                for (a, &t) in at.iter_mut().zip(&thresholds) {
                    a.decoys[f] += usize::from(q <= t);
                }
            };
        let mut writer = MbrParquetWriter::create(mbr_path)?;
        shards.for_each(|rows| {
            let shapes = DonorShapes::from_rows(rows.iter(), maps);
            let w_c = consensus_widths(rows.iter(), maps);
            let mut identified: HashMap<(u32, u16), &Evidence> = HashMap::new();
            for r in rows {
                if let (Role::Identified, Hypothesis::Target, Some(e)) =
                    (r.role, r.kind, &r.evidence)
                {
                    let previous = identified.insert((r.precursor, r.run), e);
                    debug_assert!(
                        previous.is_none(),
                        "two identified rows of precursor {} in run {}",
                        r.precursor,
                        r.run
                    );
                }
            }
            // Per row with features: its family and held-out score.
            let mut scored: Vec<Option<(u8, Option<f32>)>> = vec![None; rows.len()];
            for_each_scored(
                rows,
                &idx,
                &shapes,
                maps,
                &t_models,
                &i_models,
                |_: Table| true,
                |i, _, fam, s| {
                    scored[i] = Some((fam, s));
                },
            );
            let q_of =
                |q: &QValues, e: Estimate, i: usize| scored[i].and_then(|(fam, s)| q.q(e, fam, s));
            // A fit whose amplitude is clamped to 0 found no signal: its
            // quantity is missing, not zero, in the MBR table and the
            // feature parquet.
            let reported = |r: &ExtractRow| {
                r.evidence
                    .as_ref()
                    .filter(|e| e.amplitude > 0.0)
                    .map(|e| reported_quantity(e, r.run, w_c.get(&r.precursor).copied(), maps))
            };
            let mut out: Vec<(u32, u16, MbrRow<'_>, Option<f64>)> = Vec::new();
            for (i, r) in rows.iter().enumerate() {
                if idx.is_identity(r.precursor) {
                    // Never reported: only the identity-decoy family's counts.
                    if (r.role, r.kind) == (Role::Transfer, Hypothesis::Target) {
                        count_decoy(&mut v, &mut at, q_of(&fdr, Estimate::Pooled, i), 2);
                    }
                    continue;
                }
                if r.role == Role::Transfer
                    && r.kind == Hypothesis::Target
                    && idx.gated(r.precursor)
                {
                    v.donor_gated += 1;
                }
                if r.kind != Hypothesis::Target {
                    if r.role == Role::Transfer {
                        let f = family(r.kind) as usize;
                        count_decoy(&mut v, &mut at, q_of(&fdr, Estimate::Pooled, i), f - 1);
                        if let Some(m) = &r.ms2 {
                            v.ms2_rows[f] += 1;
                            v.ms2_supported[f] += usize::from(m.supports());
                        }
                    }
                    continue;
                }
                let e = r.evidence.as_ref();
                let quantity = reported(r);
                let quantity_normalized = quantity
                    .zip(e)
                    .and_then(|((q, _), e)| norm.normalize(r.run, r.precursor, e.t0, q));
                let (evidence_type, transfer_q, quant_q) = match r.role {
                    Role::Identified => {
                        let q = q_of(&id_fdr, Estimate::Pooled, i);
                        v.identified += 1;
                        v.identified_scored += usize::from(q.is_some());
                        v.identified_q01 += usize::from(q.is_some_and(|q| q <= VALIDATION_Q));
                        let t = match quant_feature_fdr {
                            Some(t) if !q.is_some_and(|q| q <= t) => {
                                v.identified_below_threshold += 1;
                                EvidenceType::IdentifiedBelowThreshold
                            }
                            _ => EvidenceType::Identified,
                        };
                        (t, None, q)
                    }
                    Role::Transfer => {
                        let q = q_of(&fdr, Estimate::Pooled, i);
                        v.transfers += 1;
                        v.transfers_scored += usize::from(q.is_some());
                        let (mass, rt, wrong) = (
                            q_of(&fdr, Estimate::MassOnly, i),
                            q_of(&fdr, Estimate::RtOnly, i),
                            q_of(&fdr, Estimate::WrongPeak, i),
                        );
                        for (a, &t) in at.iter_mut().zip(&thresholds) {
                            let ok = |q: Option<f64>| usize::from(q.is_some_and(|q| q <= t));
                            a.pooled += ok(q);
                            a.mass_only += ok(mass);
                            a.rt_only += ok(rt);
                            a.wrong_peak_only += ok(wrong);
                        }
                        let supported = r.ms2.is_some_and(|m| m.supports());
                        if let Some(m) = &r.ms2 {
                            v.ms2_rows[0] += 1;
                            v.ms2_supported[0] += usize::from(m.supports());
                        }
                        let t = match (accepted(q), supported) {
                            (true, true) => {
                                v.transfer_ms2 += 1;
                                EvidenceType::TransferMs2
                            }
                            (true, false) => EvidenceType::TransferMs1,
                            (false, _) => EvidenceType::BelowThreshold,
                        };
                        (t, q, q)
                    }
                    Role::PseudoTransfer => {
                        let q = scored[i].and_then(|(_, s)| fdr.pseudo_q(s));
                        v.pseudo_transfers += 1;
                        if let Some(m) = &r.ms2 {
                            v.ms2_pseudo_rows += 1;
                            v.ms2_pseudo_supported += usize::from(m.supports());
                        }
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
                        quantity: quantity.map(|q| q.0),
                        quantity_se: quantity.map(|q| q.1),
                        quantity_normalized,
                        transfer_q,
                        quant_q,
                        apex_rt_s: e.map(|e| e.t0),
                        library_q: r.library_q,
                    },
                    e.map(|e| e.mass_err_ppm),
                ));
            }
            out.sort_by_key(|&(p, run, ..)| (p, run));
            let mut features = Vec::with_capacity(out.len());
            for (p, run, r, mass_err) in &out {
                writer.push(r)?;
                features.push(MbrFeature {
                    precursor: *p,
                    run: *run,
                    evidence_type: r.evidence_type,
                    quantity: r.quantity,
                    quantity_se: r.quantity_se,
                    quantity_normalized: r.quantity_normalized,
                    transfer_q: r.transfer_q,
                    quant_q: r.quant_q,
                    apex_rt_s: r.apex_rt_s,
                    library_q: r.library_q,
                    mass_err_ppm: *mass_err,
                });
            }
            summary.rows_written += out.len();
            on_features(features)
        })?;
        let [at_q01, at_threshold] = at;
        v.at_q01 = at_q01;
        v.at_threshold = at_threshold;
        writer.finish()?;
        write_mbr_validation(validation_path, &v)?;
        summary.validation = v;
        summary.fdr = fdr.counts;
        summary.identified_fdr = id_fdr.counts;
        summary.features = all_features;
        Ok(summary)
    }

    fn unlimited() -> FinalizeLimits {
        FinalizeLimits {
            max_train_rows: usize::MAX,
            transfer_models: 1,
            max_norm_cells: usize::MAX,
            max_donor_fits: usize::MAX,
            max_held_scores: usize::MAX,
            temp_dir: None,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn stream(
        lib: &Library,
        runs: &[Vec<ExtractRow>],
        limits: FinalizeLimits,
        batch_rows: usize,
        base: &Path,
        mbr_fdr: f64,
        quant_feature_fdr: Option<f64>,
    ) -> (FinalizeSummary, Vec<MbrFeature>) {
        let (mbr, val) = mbr_paths(base);
        let params = FinalizeParams {
            mbr_fdr,
            quant_feature_fdr,
            seed: 3,
            limits,
            mbr_path: &mbr,
            validation_path: &val,
        };
        let source = ExtractRuns { runs, batch_rows };
        let mut features = Vec::new();
        let s = finalize_stream(&stems(), lib, &source, &params, &mut |f| {
            features.extend(f);
            Ok(())
        })
        .unwrap();
        (s, features)
    }

    #[test]
    fn map_reduce_matches_the_in_memory_reference() {
        let dir = tempfile::tempdir().unwrap();
        let spill_dir = tempfile::tempdir().unwrap();
        let lib = library_with_entries();
        assert_eq!(lib.entries.len(), N as usize);
        for (name, runs) in [
            ("plain", extract()),
            ("decoys", extract_with_identified_decoys()),
        ] {
            for (mbr_fdr, quant_feature_fdr) in [(0.01, None), (0.05, Some(0.01))] {
                let (m1, v1) = mbr_paths(&dir.path().join(format!("{name}-{mbr_fdr}-ref")));
                let want = reference(
                    &stems(),
                    &lib,
                    &runs,
                    mbr_fdr,
                    quant_feature_fdr,
                    3,
                    &m1,
                    &v1,
                )
                .unwrap();
                assert!(want.rows_written > 0 && want.validation.transfers_scored > 0);
                let arms = [
                    ("one", usize::MAX, unlimited()),
                    ("batches", 7, unlimited()),
                    (
                        "groups",
                        97,
                        FinalizeLimits {
                            max_donor_fits: 100,
                            ..unlimited()
                        },
                    ),
                    (
                        "spill",
                        1000,
                        FinalizeLimits {
                            max_held_scores: 50,
                            temp_dir: Some(spill_dir.path().to_path_buf()),
                            ..unlimited()
                        },
                    ),
                ];
                for (arm, batch_rows, limits) in arms {
                    let base = dir.path().join(format!("{name}-{mbr_fdr}-{arm}"));
                    let (s, features) = stream(
                        &lib,
                        &runs,
                        limits,
                        batch_rows,
                        &base,
                        mbr_fdr,
                        quant_feature_fdr,
                    );
                    let (m2, v2) = mbr_paths(&base);
                    let what = format!("{name} {mbr_fdr} {arm}");
                    assert_eq!(
                        std::fs::read(&m1).unwrap(),
                        std::fs::read(&m2).unwrap(),
                        "{what}"
                    );
                    assert_eq!(
                        std::fs::read(&v1).unwrap(),
                        std::fs::read(&v2).unwrap(),
                        "{what}"
                    );
                    assert_eq!(features, want.features, "{what}");
                    assert_eq!(s.fdr, want.fdr, "{what}");
                    assert_eq!(s.identified_fdr, want.identified_fdr, "{what}");
                    assert_eq!(s.validation, want.validation, "{what}");
                    assert_eq!(
                        (s.no_features, s.identified_no_features, s.rows_written),
                        (
                            want.no_features,
                            want.identified_no_features,
                            want.rows_written
                        ),
                        "{what}"
                    );
                }
            }
        }
        for d in [dir.path(), spill_dir.path()] {
            let left: Vec<_> = std::fs::read_dir(d)
                .unwrap()
                .filter_map(|e| e.ok())
                .filter(|e| {
                    e.file_name()
                        .to_string_lossy()
                        .starts_with(".andes-finalize-")
                })
                .collect();
            assert!(left.is_empty(), "temporary tables are removed");
        }
    }

    #[test]
    fn sampled_training_reports_every_row() {
        let dir = tempfile::tempdir().unwrap();
        let runs = extract();
        let lib = library_with_entries();
        let (m1, v1) = mbr_paths(&dir.path().join("one"));
        let one = finalize(&stems(), &lib, &runs, 0.01, None, 3, &m1, &v1).unwrap();
        let limits = FinalizeLimits {
            max_train_rows: 600,
            max_norm_cells: 100,
            ..unlimited()
        };
        let (s, features) = stream(
            &lib,
            &runs,
            limits,
            1000,
            &dir.path().join("sampled"),
            0.01,
            None,
        );
        assert_eq!(s.rows_written, one.rows_written);
        assert_eq!(s.validation.transfers, one.validation.transfers);
        assert_eq!(s.validation.identified, one.validation.identified);
        let keys = |f: &[MbrFeature]| f.iter().map(|f| (f.precursor, f.run)).collect::<Vec<_>>();
        assert_eq!(keys(&features), keys(&one.features));
        // Sampled models still separate the present targets from the decoys.
        assert!(s.validation.transfers_scored > 0);
        let accepted = features
            .iter()
            .filter(|f| f.evidence_type == EvidenceType::TransferMs1)
            .count();
        assert!(accepted > 0, "{accepted}");
    }

    #[test]
    fn bagged_transfer_classifier_reports_every_row_and_is_reproducible() {
        let dir = tempfile::tempdir().unwrap();
        let runs = extract_with_identified_decoys();
        let lib = library_with_entries();
        let run = |name: &str, members: usize, batch_rows: usize| {
            let limits = FinalizeLimits {
                transfer_models: members,
                max_train_rows: 900,
                ..unlimited()
            };
            let base = dir.path().join(name);
            let (s, features) = stream(&lib, &runs, limits, batch_rows, &base, 0.05, None);
            let (m, v) = mbr_paths(&base);
            (
                s,
                features,
                std::fs::read(m).unwrap(),
                std::fs::read(v).unwrap(),
            )
        };
        let (one, one_f, _, _) = run("one", 1, 1000);
        let (a, a_f, a_m, a_v) = run("a", 3, 1000);
        let (_, _, b_m, b_v) = run("b", 3, 7);
        // Same seed, other batches: the same bytes.
        assert_eq!((a_m, a_v), (b_m, b_v));
        assert_eq!(a.rows_written, one.rows_written);
        assert_eq!(a.validation.transfers, one.validation.transfers);
        assert_eq!(
            a.validation.transfers_scored,
            one.validation.transfers_scored
        );
        let keys = |f: &[MbrFeature]| f.iter().map(|f| (f.precursor, f.run)).collect::<Vec<_>>();
        assert_eq!(keys(&a_f), keys(&one_f));
        assert!(a.fdr.failed_folds.is_empty(), "{:?}", a.fdr.failed_folds);
        let accepted = a_f
            .iter()
            .filter(|f| f.evidence_type == EvidenceType::TransferMs1)
            .count();
        assert!(accepted > 0, "{accepted}");
        // The members change the scores, so some q-value differs.
        assert!(a_f
            .iter()
            .zip(&one_f)
            .any(|(x, y)| x.transfer_q != y.transfer_q));
    }

    #[test]
    fn finalize_rejects_rows_of_a_precursor_apart() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let mut runs = extract();
        let first = runs[1].remove(0);
        runs[1].push(first);
        let e = finalize(&stems(), &library(), &runs, 0.01, None, 1, &mbr, &val).unwrap_err();
        assert!(e.contains("not contiguous"), "{e}");
    }

    #[test]
    fn finalize_rejects_out_of_range_rows_and_thresholds() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let mut runs = extract();
        runs[1].push(row(0, 2, Hypothesis::Target, Role::Transfer, None));
        let e = finalize(&stems(), &library(), &runs, 0.01, None, 1, &mbr, &val).unwrap_err();
        assert!(e.contains("run 2"), "{e}");
        runs[1].pop();
        runs[1].push(row(N, 1, Hypothesis::Target, Role::Transfer, None));
        let e = finalize(&stems(), &library(), &runs, 0.01, None, 1, &mbr, &val).unwrap_err();
        assert!(e.contains(&format!("precursor {N}")), "{e}");
        for bad in [0.0, -0.1, 1.5, f64::NAN] {
            assert!(finalize(&stems(), &library(), &extract(), bad, None, 1, &mbr, &val).is_err());
        }
        assert!(!mbr.exists());
    }

    /// Run 0 on the consensus scale, run 1 stretched by 1.25 and run 2
    /// unaligned.
    fn stretched_maps() -> Vec<Option<RunMap>> {
        let map = |slope: f64| {
            RunMap::from_parts(vec![0.0, 3000.0], vec![0.0, 3000.0 * slope], slope, slope).unwrap()
        };
        vec![Some(map(1.0)), Some(map(1.25)), None]
    }

    /// A precursor seen in runs 0 and 1 (and 2, unaligned) with the same
    /// abundance: in run `r` the apex is `A0 / slope_r` and the true width
    /// `w0 * slope_r`, and the fitted width is off by `jitter[r]`, as a
    /// coarse width grid makes it, while the amplitude is right.
    fn jittered_rows(p: u32, jitter: [f64; 3]) -> Vec<ExtractRow> {
        let (a0, w0) = (1e6, 5.0);
        let slopes = [1.0, 1.25, 1.0];
        let mut g = Rng(p as u64 + 1);
        (0..3u16)
            .flat_map(|run| {
                let slope = slopes[run as usize];
                let w = w0 * slope * jitter[run as usize];
                let amplitude = a0 / slope;
                let mut e = evidence(&mut g, true, 500.0 * slope);
                e.amplitude = amplitude;
                e.amplitude_se = 0.01 * amplitude;
                e.w = w;
                e.quantity = amplitude * w * (2.0 * std::f64::consts::PI).sqrt();
                e.quantity_se = 0.01 * e.quantity;
                let identified = row(
                    p,
                    run,
                    Hypothesis::Target,
                    Role::Identified,
                    Some(e.clone()),
                );
                let pseudo = row(p, run, Hypothesis::Target, Role::PseudoTransfer, Some(e));
                [identified, pseudo]
            })
            .collect()
    }

    #[test]
    fn consensus_width_gives_one_quantity_across_runs() {
        let maps = stretched_maps();
        let truth = 1e6 * 5.0 * (2.0 * std::f64::consts::PI).sqrt();
        for jitter in [[0.7, 1.4, 1.0], [1.4, 0.7, 1.4], [1.0, 1.0, 0.7]] {
            let rows = jittered_rows(3, jitter);
            let w_c = consensus_widths(&rows, &maps);
            let per_run: Vec<(f64, f64)> = rows
                .iter()
                .filter(|r| r.role == Role::Identified)
                .map(|r| {
                    let e = r.evidence.as_ref().unwrap();
                    reported_quantity(e, r.run, w_c.get(&3).copied(), &maps)
                })
                .collect();
            let raw: Vec<f64> = rows
                .iter()
                .filter(|r| r.role == Role::Identified)
                .map(|r| r.evidence.as_ref().unwrap().quantity)
                .collect();
            eprintln!("jitter {jitter:?}: raw {raw:?}, consensus {per_run:?}");
            // Runs 0 and 1 report one quantity: the jitter is gone.
            let (q0, q1) = (per_run[0].0, per_run[1].0);
            assert!((q0 / q1 - 1.0).abs() < 1e-12, "{per_run:?}");
            // w_c is the median of the two aligned runs' widths.
            let want = truth * 0.5 * (jitter[0] + jitter[1]);
            assert!((q0 / want - 1.0).abs() < 1e-12, "{q0} vs {want}");
            for &(q, se) in &per_run[..2] {
                assert!((se / q - 0.01).abs() < 1e-12);
            }
            // The unaligned run keeps the extract's own values.
            let e2 = rows[4].evidence.as_ref().unwrap();
            assert_eq!(per_run[2], (e2.quantity, e2.quantity_se));
        }
        // Raw quantities with jitter 0.7 vs 1.4 differ twofold; the
        // consensus removes it.
        let rows = jittered_rows(3, [0.7, 1.4, 1.0]);
        let raw = |run: usize| rows[2 * run].evidence.as_ref().unwrap().quantity;
        assert!((raw(1) / raw(0) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn consensus_width_skips_unfitted_and_unaligned_rows() {
        let maps = stretched_maps();
        let mut rows = jittered_rows(4, [1.0, 1.0, 1.0]);
        // Only run 0's identified fit has signal; run 2 is unaligned.
        rows[2].evidence.as_mut().unwrap().amplitude = 0.0;
        rows[4].evidence.as_mut().unwrap().w = 50.0;
        let w_c = consensus_widths(&rows, &maps);
        assert!((w_c[&4] - 5.0).abs() < 1e-12, "{w_c:?}");
        // A precursor without an identified fit keeps its own quantity.
        let mut rows = jittered_rows(5, [1.0, 1.0, 1.0]);
        for r in rows.iter_mut().filter(|r| r.role == Role::Identified) {
            r.evidence = None;
        }
        assert!(consensus_widths(&rows, &maps).is_empty());
        let e = rows[1].evidence.as_ref().unwrap();
        assert_eq!(
            reported_quantity(e, 1, None, &maps),
            (e.quantity, e.quantity_se)
        );
    }

    #[test]
    fn finalize_reports_the_consensus_quantity() {
        // Precursor 0 is identified in runs 0 and 1 with widths 4 and 8; with
        // unit-slope maps each row reports its A * 6 * sqrt(2 pi).
        let mut lib = library();
        let unit = RunMap::from_parts(vec![0.0, 3000.0], vec![0.0, 3000.0], 1.0, 1.0).unwrap();
        lib.alignment.maps = vec![Some(unit.clone()), Some(unit)];
        let mut runs = extract();
        let mut g = Rng(9);
        for (run, w) in [(0u16, 4.0), (1, 8.0)] {
            let mut e = evidence(&mut g, true, 100.0);
            e.amplitude = 4e6 / w;
            e.amplitude_se = 1e3;
            e.w = w;
            e.quantity = e.amplitude * w * (2.0 * std::f64::consts::PI).sqrt();
            let rows = &mut runs[run as usize];
            rows.retain(|r| !(r.precursor == 0 && r.role != Role::PseudoTransfer));
            // A precursor's rows of a run are contiguous.
            let at = rows
                .iter()
                .position(|r| r.precursor != 0)
                .unwrap_or(rows.len());
            rows.insert(
                at,
                row(0, run, Hypothesis::Target, Role::Identified, Some(e)),
            );
        }
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let s = finalize(&stems(), &lib, &runs, 0.01, None, 7, &mbr, &val).unwrap();
        let mut p0 = s.features.iter().filter(|f| f.precursor == 0);
        let sqrt_2pi = (2.0 * std::f64::consts::PI).sqrt();
        let a = p0.next().unwrap();
        let b = p0.next().unwrap();
        assert_eq!((a.run, b.run), (0, 1));
        assert!((a.quantity.unwrap() / (1e6 * 6.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
        assert!((b.quantity.unwrap() / (5e5 * 6.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
        assert!((a.quantity_se.unwrap() / (1e3 * 6.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
        // The extract rows keep their per-run quantities.
        let id0 = runs[0]
            .iter()
            .find(|r| r.precursor == 0 && r.role == Role::Identified);
        assert!(
            (id0.unwrap().evidence.as_ref().unwrap().quantity / (4e6 * sqrt_2pi) - 1.0).abs()
                < 1e-12
        );
    }

    /// Evidence of another precursor's feature caught in the window: strong,
    /// but with the wrong isotope pattern and a poor elution profile.
    fn foreign(g: &mut Rng, t0: f64) -> Evidence {
        let mut e = evidence(g, true, t0);
        e.llr = 50.0 + 15.0 * g.normal();
        e.envelope_cos = 0.7 + 0.1 * g.normal();
        e.coherence = 0.4 + 0.15 * g.normal();
        e.rt_z = 2.0 * g.normal();
        e.m1_ratio = 0.5 + 0.2 * g.unit();
        e.isotopes_used = 0b0111;
        e.mass_err_ppm = 5.0 * g.normal();
        e
    }

    /// Identified precursor `p` in run 0 is clean, at the wrong peak (no
    /// signal of its own: `p % 10 == 0`) or interfered (`p % 10 == 5`).
    fn identified_kind(p: u32) -> &'static str {
        match p % 10 {
            0 => "wrong",
            5 => "interfered",
            _ => "clean",
        }
    }

    /// [`extract`] with run 0's identified targets as [`identified_kind`]
    /// says, each followed by a mass-shift and an RT-shift decoy; a third of
    /// the decoys catch another precursor's feature.
    fn extract_with_identified_decoys() -> Vec<Vec<ExtractRow>> {
        let mut runs = extract();
        let mut g = Rng(11);
        let mut run0 = Vec::new();
        for r in runs[0].drain(..) {
            if r.role != Role::Identified {
                run0.push(r);
                continue;
            }
            let (p, apex) = (r.precursor, r.evidence.as_ref().unwrap().t0);
            let target = match identified_kind(p) {
                "wrong" => evidence(&mut g, false, apex),
                "interfered" => foreign(&mut g, apex),
                _ => evidence(&mut g, true, apex),
            };
            run0.push(row(
                p,
                0,
                Hypothesis::Target,
                Role::Identified,
                Some(target),
            ));
            for (k, kind) in [Hypothesis::MassShiftDecoy(-7), Hypothesis::RtShiftDecoy]
                .into_iter()
                .enumerate()
            {
                let t0 = apex + 300.0 * k as f64;
                let e = if (p + k as u32) % 3 == 0 {
                    foreign(&mut g, t0)
                } else {
                    evidence(&mut g, false, t0)
                };
                run0.push(row(p, 0, kind, Role::Identified, Some(e)));
            }
        }
        runs[0] = run0;
        runs
    }

    #[test]
    fn identified_quality_q_separates_clean_from_wrong_and_interfered_peaks() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let runs = extract_with_identified_decoys();
        let s = finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
        let mut q: HashMap<&str, Vec<f64>> = HashMap::new();
        for f in s.features.iter().filter(|f| f.run == 0) {
            assert_eq!(f.evidence_type, EvidenceType::Identified);
            assert_eq!(f.transfer_q, None);
            q.entry(identified_kind(f.precursor))
                .or_default()
                .push(f.quant_q.unwrap());
        }
        let share = |kind: &str, ok: &dyn Fn(f64) -> bool| {
            let v = &q[kind];
            v.iter().filter(|&&x| ok(x)).count() as f64 / v.len() as f64
        };
        let clean = share("clean", &|x| x <= 0.01);
        let wrong = share("wrong", &|x| x > 0.05);
        let interfered = share("interfered", &|x| x > 0.05);
        eprintln!("clean q<=0.01 {clean}, wrong q>0.05 {wrong}, interfered q>0.05 {interfered}");
        assert!(clean >= 0.95, "{clean}");
        assert!(wrong >= 0.95, "{wrong}");
        assert!(interfered >= 0.8, "{interfered}");
        let v = &s.validation;
        assert_eq!(
            (v.identified, v.identified_scored),
            (N as usize, N as usize)
        );
        assert_eq!(v.identified_below_threshold, 0);
        assert_eq!(s.identified_no_features, 0);
        // A transfer's quant_q is its transfer_q.
        for f in s.features.iter().filter(|f| f.run == 1) {
            assert_eq!(f.quant_q, f.transfer_q);
        }
    }

    /// [`extract_with_identified_decoys`] with run 1 identifying the even
    /// precursors as run 0 does (same targets and decoys) at `gain` times run
    /// 0's quantity; the odd precursors stay transfers into run 1.
    fn extract_two_identified_runs(gain: f64) -> Vec<Vec<ExtractRow>> {
        let mut runs = extract_with_identified_decoys();
        let copies: Vec<ExtractRow> = runs[0]
            .iter()
            .filter(|r| r.precursor % 2 == 0 && r.role == Role::Identified)
            .map(|r| {
                let mut r = r.clone();
                r.run = 1;
                if let Some(e) = r.evidence.as_mut() {
                    e.amplitude *= gain;
                    e.amplitude_se *= gain;
                    e.quantity *= gain;
                    e.quantity_se *= gain;
                }
                r
            })
            .collect();
        let run1: Vec<ExtractRow> = runs[1]
            .drain(..)
            .filter(|r| r.precursor % 2 == 1)
            .chain(copies)
            .collect();
        runs[1] = run1;
        runs
    }

    #[test]
    fn finalize_reports_run_normalized_quantities() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let runs = extract_two_identified_runs(2.0);
        let s = finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
        let get = |p: u32, run: u16| {
            s.features
                .iter()
                .find(|f| f.precursor == p && f.run == run)
                .unwrap()
        };
        let mut clean = 0;
        for p in (0..N).filter(|p| p % 2 == 0 && identified_kind(*p) == "clean") {
            let (a, b) = (get(p, 0), get(p, 1));
            // Raw run 1 reads twice run 0; normalised, the two runs agree and
            // sit halfway (log2 offsets -0.5 and +0.5).
            assert!((b.quantity.unwrap() / a.quantity.unwrap() - 2.0).abs() < 1e-9);
            let (na, nb) = (
                a.quantity_normalized.unwrap(),
                b.quantity_normalized.unwrap(),
            );
            assert!((nb / na - 1.0).abs() < 1e-9, "{p}: {na} {nb}");
            assert!((na / (a.quantity.unwrap() * 2f64.sqrt()) - 1.0).abs() < 1e-9);
            clean += 1;
        }
        assert!(clean > 300, "{clean}");
        // Every other row with a quantity is scaled by its run's offset:
        // transfers into run 1 and identified rows that are not references.
        for f in s.features.iter().filter(|f| f.quantity.is_some()) {
            let want = f.quantity.unwrap()
                * if f.run == 0 {
                    2f64.sqrt()
                } else {
                    0.5f64.sqrt()
                };
            let got = f.quantity_normalized.unwrap();
            assert!((got / want - 1.0).abs() < 1e-9, "{f:?}");
        }
        // The column is in the table, after quantity_se.
        let t = table(&mbr);
        let col = |name: &str| t[0].iter().position(|h| h == name).unwrap();
        assert_eq!(col("quantity_normalized"), col("quantity_se") + 1);
    }

    #[test]
    fn finalize_leaves_quantities_unnormalized_without_two_reference_runs() {
        // Only run 0 has identified rows: no run is normalised.
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let runs = extract_with_identified_decoys();
        let s = finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
        assert!(s.features.iter().all(|f| f.quantity_normalized.is_none()));
        let t = table(&mbr);
        let col = t[0]
            .iter()
            .position(|h| h == "quantity_normalized")
            .unwrap();
        assert!(t[1..].iter().all(|r| r[col] == "NA"));
    }

    #[test]
    fn identified_quality_leaves_the_transfer_q_values_unchanged() {
        let dir = tempfile::tempdir().unwrap();
        let (m1, v1) = mbr_paths(&dir.path().join("plain"));
        let (m2, v2) = mbr_paths(&dir.path().join("decoys"));
        let plain = finalize(&stems(), &library(), &extract(), 0.01, None, 7, &m1, &v1).unwrap();
        let runs = extract_with_identified_decoys();
        let with = finalize(&stems(), &library(), &runs, 0.01, None, 7, &m2, &v2).unwrap();
        let transfers = |s: &FinalizeSummary| -> Vec<(u32, Option<f64>, EvidenceType)> {
            s.features
                .iter()
                .filter(|f| f.run == 1)
                .map(|f| (f.precursor, f.transfer_q, f.evidence_type))
                .collect()
        };
        assert_eq!(transfers(&plain), transfers(&with));
        assert_eq!(plain.validation.at_q01, with.validation.at_q01);
        // Without identified decoys there is no identified quality q-value.
        assert!(plain
            .features
            .iter()
            .filter(|f| f.run == 0)
            .all(|f| f.quant_q.is_none() && f.evidence_type == EvidenceType::Identified));
    }

    #[test]
    fn quant_feature_fdr_flags_identified_rows_above_it() {
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let runs = extract_with_identified_decoys();
        let s = finalize(&stems(), &library(), &runs, 0.01, Some(0.01), 7, &mbr, &val).unwrap();
        let (mut flagged, mut kept) = (0, 0);
        for f in s.features.iter().filter(|f| f.run == 0) {
            let below = f.quant_q.is_none_or(|q| q > 0.01);
            let want = if below {
                EvidenceType::IdentifiedBelowThreshold
            } else {
                EvidenceType::Identified
            };
            assert_eq!(f.evidence_type, want, "{f:?}");
            // A flagged row keeps its quantity.
            assert!(f.quantity.is_some());
            if below {
                flagged += 1;
            } else {
                kept += 1;
            }
        }
        assert_eq!(s.validation.identified_below_threshold, flagged);
        assert!(
            flagged >= 160 && kept >= 680,
            "{flagged} flagged, {kept} kept"
        );
        let t = table(&mbr);
        assert_eq!(t[0][7], "quant_q");
        let flagged_rows = t
            .iter()
            .filter(|r| r[2] == "identified_below_threshold")
            .count();
        assert_eq!(flagged_rows, flagged);
        let t = table(&val);
        let get = |name: &str| &t[1][t[0].iter().position(|h| h == name).unwrap()];
        assert_eq!(get("quant_feature_threshold"), "0.01");
        assert_eq!(get("identified_below_threshold"), &flagged.to_string());
        for bad in [0.0, 1.5, f64::NAN] {
            assert!(finalize(&stems(), &library(), &runs, 0.01, Some(bad), 7, &mbr, &val).is_err());
        }
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
        let x = transfer_features(&r, &DonorShapes::default(), &[], 0.0).unwrap();
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
        assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0)
            .unwrap()
            .iter()
            .all(|v| v.is_finite()));
        r.evidence.as_mut().unwrap().llr = f64::NAN;
        assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0).is_none());
        r.evidence = None;
        assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0).is_none());
    }

    use quant::ms2::Ms2Evidence;

    fn supporting() -> Ms2Evidence {
        Ms2Evidence {
            n_scans: 2,
            best: quant::ms2::ScanMatch {
                n_matched: 9,
                angle: 0.9,
                explained: 0.95,
            },
        }
    }

    fn empty_ms2() -> Ms2Evidence {
        Ms2Evidence {
            n_scans: 1,
            best: quant::ms2::ScanMatch::default(),
        }
    }

    /// Identity maps for `n` runs (unit slope everywhere).
    fn identity_maps(n: usize) -> Vec<Option<RunMap>> {
        (0..n)
            .map(|_| {
                Some(RunMap::from_parts(vec![0.0, 1000.0], vec![0.0, 1000.0], 1.0, 1.0).unwrap())
            })
            .collect()
    }

    /// A fit of width `w`, tailing `tail` and isotope heights `iso`.
    fn shaped(w: f64, tail: f64, iso: [f64; 4]) -> Evidence {
        let mut e = evidence(&mut Rng(9), true, 100.0);
        e.w = w;
        e.tailing = tail;
        e.iso_heights = iso;
        e
    }

    const DONOR_ISO: [f64; 4] = [5.5e5, 3.0e5, 1.1e5, 0.4e5];

    #[test]
    fn shape_features_compare_with_the_donor_fits_of_other_runs() {
        let maps = identity_maps(3);
        let id = |run: u16, w: f64, tail: f64, iso| {
            row(
                7,
                run,
                Hypothesis::Target,
                Role::Identified,
                Some(shaped(w, tail, iso)),
            )
        };
        let donors = [id(0, 4.0, 0.10, DONOR_ISO), id(2, 4.0, 0.10, DONOR_ISO)];
        let shapes = DonorShapes::from_rows(&donors, &maps);
        // A transfer at the donors' shape agrees on every count.
        let like = row(
            7,
            1,
            Hypothesis::Target,
            Role::Transfer,
            Some(shaped(4.0, 0.10, DONOR_ISO.map(|v| 2.0 * v))),
        );
        let f = shapes.features(&like, &maps);
        assert!(
            f[0].abs() < 1e-12 && f[1].abs() < 1e-12 && (f[2] - 1.0).abs() < 1e-12,
            "{f:?}"
        );
        // A decoy of that precursor on another peak (wider, other isotope
        // ratios) is compared with the same donors and disagrees.
        let other = row(
            7,
            1,
            Hypothesis::MassShiftDecoy(7),
            Role::Transfer,
            Some(shaped(8.0, -0.3, [2.0e5, 4.0e5, 3.0e5, 1.0e5])),
        );
        let g = shapes.features(&other, &maps);
        assert!((g[0] - 2f64.ln()).abs() < 1e-12, "{g:?}");
        assert!((g[1] + 0.4).abs() < 1e-12, "{g:?}");
        assert!(g[2] < f[2] - 0.1, "{g:?}");
        // The same numbers end the transfer feature vector.
        let x = transfer_features(&other, &shapes, &maps, 0.0).unwrap();
        let n = TRANSFER_FEATURES.len();
        assert_eq!(x.len(), n);
        assert_eq!(
            &TRANSFER_FEATURES[n - 3..],
            ["shape_log_width", "shape_tailing_diff", "shape_iso_cos"]
        );
        assert_eq!(&x[n - 3..], &g.map(|v| v as f32));
    }

    #[test]
    fn shape_features_leave_out_the_rows_own_run() {
        let maps = identity_maps(2);
        // Run 1 identifies the precursor with an odd shape; run 0 is the
        // only other donor.
        let donors = [
            row(
                3,
                0,
                Hypothesis::Target,
                Role::Identified,
                Some(shaped(4.0, 0.0, DONOR_ISO)),
            ),
            row(
                3,
                1,
                Hypothesis::Target,
                Role::Identified,
                Some(shaped(9.0, 0.5, [1.0, 1.0, 1.0, 1.0])),
            ),
        ];
        let shapes = DonorShapes::from_rows(&donors, &maps);
        let pseudo = row(
            3,
            1,
            Hypothesis::Target,
            Role::PseudoTransfer,
            Some(shaped(9.0, 0.5, [1.0, 1.0, 1.0, 1.0])),
        );
        let f = shapes.features(&pseudo, &maps);
        assert!((f[0] - (9.0f64 / 4.0).ln()).abs() < 1e-12, "{f:?}");
        assert!((f[1] - 0.5).abs() < 1e-12, "{f:?}");
        assert!(f[2] < 0.99, "{f:?}");
        // Without another run, or without a map or a fit, there is no
        // reference.
        let only = DonorShapes::from_rows(&donors[1..], &maps);
        assert_eq!(only.features(&pseudo, &maps), NO_SHAPE);
        assert_eq!(shapes.features(&pseudo, &[None, None]), NO_SHAPE);
        let mut empty = pseudo.clone();
        empty.evidence.as_mut().unwrap().amplitude = 0.0;
        assert_eq!(shapes.features(&empty, &maps), NO_SHAPE);
        // Decoy and identified decoy rows never become references.
        let decoy = row(
            3,
            0,
            Hypothesis::RtShiftDecoy,
            Role::Identified,
            Some(shaped(1.0, 0.0, DONOR_ISO)),
        );
        let with_decoy = DonorShapes::from_rows(donors.iter().chain([&decoy]), &maps);
        assert_eq!(with_decoy.features(&pseudo, &maps), f);
    }

    #[test]
    fn transfer_features_leave_out_the_ms2_evidence() {
        let mut g = Rng(3);
        let mut r = row(
            0,
            1,
            Hypothesis::Target,
            Role::Transfer,
            Some(evidence(&mut g, true, 100.0)),
        );
        let x = transfer_features(&r, &DonorShapes::default(), &[], 0.0).unwrap();
        r.ms2 = Some(supporting());
        assert_eq!(
            transfer_features(&r, &DonorShapes::default(), &[], 0.0).unwrap(),
            x
        );
        assert!(TRANSFER_FEATURES.iter().all(|f| !f.starts_with("ms2")));
    }

    #[test]
    fn supported_accepted_transfers_are_transfer_ms2() {
        let mut runs = extract();
        for r in runs[1].iter_mut() {
            r.ms2 = Some(if r.kind == Hypothesis::Target && r.precursor % 4 == 0 {
                supporting()
            } else {
                empty_ms2()
            });
        }
        let dir = tempfile::tempdir().unwrap();
        let (mbr, val) = mbr_paths(&dir.path().join("x"));
        let s = finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
        let t = table(&mbr);
        let ms2_rows: Vec<&Vec<String>> =
            t[1..].iter().filter(|r| r[2] == "transfer_ms2").collect();
        assert!(!ms2_rows.is_empty());
        for r in &ms2_rows {
            let p: u32 = r[0]
                .trim_start_matches("PEP")
                .trim_end_matches("/2")
                .parse()
                .unwrap();
            assert_eq!(p % 4, 0, "{r:?}");
        }
        let v = &s.validation;
        assert_eq!(v.transfer_ms2, ms2_rows.len());
        assert_eq!(v.ms2_rows[0], N as usize);
        assert_eq!(v.ms2_supported[0], (N as usize).div_ceil(4));
        assert_eq!((v.ms2_supported[1], v.ms2_supported[2]), (0, 0));
        assert!(v.ms2_rows[1] > 0 && v.ms2_rows[2] > 0);
        let vt = table(&val);
        let get = |name: &str| &vt[1][vt[0].iter().position(|h| h == name).unwrap()];
        assert_eq!(get("transfer_ms2"), &ms2_rows.len().to_string());
    }

    /// Run 1: every transfer target and every decoy sits on a clean MS1 peak
    /// (decoys on other peptides' peaks); only the present targets (even
    /// precursors) have an MS2 scan that matches their library spectrum.
    fn ambiguous_ms1(with_ms2: bool) -> Vec<Vec<ExtractRow>> {
        let mut runs = extract();
        let mut g = Rng(11);
        let mut run1 = Vec::new();
        for p in 0..N {
            let apex = 100.0 + p as f64;
            for kind in [
                Hypothesis::Target,
                Hypothesis::MassShiftDecoy(7),
                Hypothesis::RtShiftDecoy,
            ] {
                let mut r = row(
                    p,
                    1,
                    kind,
                    Role::Transfer,
                    Some(evidence(&mut g, true, apex)),
                );
                if with_ms2 {
                    r.ms2 = Some(if kind == Hypothesis::Target && p % 2 == 0 {
                        supporting()
                    } else {
                        empty_ms2()
                    });
                }
                run1.push(r);
            }
        }
        runs[1] = run1;
        runs
    }

    #[test]
    fn ms2_evidence_does_not_change_the_transfer_q_values() {
        let accepted = |with_ms2: bool| {
            let dir = tempfile::tempdir().unwrap();
            let (mbr, val) = mbr_paths(&dir.path().join("x"));
            let runs = ambiguous_ms1(with_ms2);
            finalize(&stems(), &library(), &runs, 0.01, None, 7, &mbr, &val).unwrap();
            let (mut present, mut absent) = (0usize, 0usize);
            for r in &table(&mbr)[1..] {
                if r[1] != "b" || !r[2].starts_with("transfer_ms") {
                    continue;
                }
                let p: u32 = r[0]
                    .trim_start_matches("PEP")
                    .trim_end_matches("/2")
                    .parse()
                    .unwrap();
                if p % 2 == 0 {
                    present += 1;
                } else {
                    absent += 1;
                }
            }
            (present, absent)
        };
        // Decoys on clean peaks make targets and decoys indistinguishable in
        // MS1; MS2 support of the present targets does not change which
        // transfers are accepted (it only labels them).
        let (p0, a0) = accepted(false);
        let (p1, a1) = accepted(true);
        eprintln!(
            "MS1 only: {p0} present, {a0} absent accepted; with MS2: {p1} present, {a1} absent"
        );
        assert_eq!((p0, a0), (p1, a1));
        assert!(p0 + a0 < 50, "{p0} + {a0}");
    }

    #[test]
    fn attach_ms2_scores_each_probe_on_its_own_window() {
        use quant::ms2::{fragment_ladder, library_spectrum, residue_masses, Ms2Scan};
        let lib = library();
        let tol = model::Tolerance::Ppm(20.0);
        let ladder = fragment_ladder(
            &residue_masses("LGEYGFQNALIVR", &std::collections::HashMap::new()).unwrap(),
            2,
        );
        let donor: Vec<(f32, f32)> = ladder
            .iter()
            .enumerate()
            .map(|(k, f)| (f.mz as f32, 100.0 + k as f32))
            .collect();
        let spectrum = library_spectrum(&ladder, &donor, tol).unwrap();
        let mut spectra = vec![None; lib.meta.len()];
        spectra[0] = Some(spectrum);
        let mut donor_scans = vec![Vec::new(); lib.meta.len()];
        donor_scans[0] = vec![(0u16, 50.0), (1u16, 100.0)];
        let noise = vec![(321.1f32, 10.0f32), (654.3, 20.0)];
        let target_mz = lib.meta[0].mz;
        let decoy_mz = probe_mz(&lib.meta[0], Hypothesis::MassShiftDecoy(7));
        assert!((decoy_mz - target_mz - 7.0 * PEPTIDE_MASS_SLOPE / 2.0).abs() < 1e-12);
        let scan = |rt: f64, mz: f64, peaks: &Vec<(f32, f32)>| Ms2Scan {
            native_id: format!("scan={rt}"),
            spec_key: format!("k{rt}"),
            rt,
            precursor_mz: mz,
            charge: 2,
            isolation: Some((mz, 0.7, 0.7)),
            faims_cv: None,
            peaks: peaks.clone(),
        };
        let index = Ms2RunIndex::new(vec![
            scan(100.0, target_mz, &donor),
            scan(101.0, target_mz, &noise),
            scan(200.0, decoy_mz, &noise),
            scan(400.0, target_mz, &donor),
        ]);
        let ev = |t0: f64| {
            let mut g = Rng(1);
            Some(evidence(&mut g, true, t0))
        };
        let mut rows = vec![
            row(0, 1, Hypothesis::Target, Role::Transfer, ev(100.0)),
            row(
                0,
                1,
                Hypothesis::MassShiftDecoy(7),
                Role::Transfer,
                ev(200.0),
            ),
            row(0, 1, Hypothesis::RtShiftDecoy, Role::Transfer, ev(400.0)),
            row(0, 1, Hypothesis::Target, Role::PseudoTransfer, ev(100.0)),
            row(0, 1, Hypothesis::Target, Role::Identified, ev(100.0)),
            row(0, 1, Hypothesis::Target, Role::Transfer, None),
            row(1, 1, Hypothesis::Target, Role::Transfer, ev(100.0)),
        ];
        let c = attach_ms2(&mut rows, &lib, &spectra, &donor_scans, &index, tol);
        let m = |i: usize| rows[i].ms2.unwrap();
        // The target sees both scans at its m/z; the donor-like one wins.
        assert_eq!(m(0).n_scans, 2);
        assert!(m(0).supports());
        // The mass-shift decoy sees only the scan of its own m/z.
        assert_eq!(m(1).n_scans, 1);
        assert!(!m(1).supports());
        // The RT-shift decoy at the same m/z elsewhere in the run is scored
        // with the same spectrum: a scan there with that spectrum supports it.
        assert_eq!(m(2).n_scans, 1);
        assert!(m(2).supports());
        // The pseudo-transfer leaves out its run's donor scan (RT 100).
        assert_eq!(m(3).n_scans, 1);
        assert!(!m(3).supports());
        // Identified rows, rows without a fit and entries without a spectrum
        // get none.
        assert!(rows[4].ms2.is_none() && rows[5].ms2.is_none() && rows[6].ms2.is_none());
        assert_eq!(c.rows[0], [1, 1, 1]);
        assert_eq!(c.supported[0], [1, 0, 1]);
        assert_eq!((c.rows[1][0], c.supported[1][0]), (1, 0));
        assert_eq!(c.no_spectrum, 1);
    }
}
