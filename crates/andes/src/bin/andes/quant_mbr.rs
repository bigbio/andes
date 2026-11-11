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
//!   this run's own observation. The best PSM of a run is not always taken
//!   on the elution peak: an abundant peptide is fragmented again on its
//!   tail, minutes after the apex, and that PSM can score best. The
//!   identified fit, which stays within three elution widths of the PSM,
//!   then measures the tail. So when the fit at the position the other runs
//!   predict is a different peak (apexes further apart than the two fitted
//!   widths together) with more than [`CONSENSUS_PEAK_FACTOR`] times the
//!   amplitude of a positive identified fit, the identified row takes that
//!   fit ([`prefer_consensus_peak`]; counted in the extract log). An
//!   identified fit without signal, or a lower one at the predicted
//!   position, changes nothing;
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
//! quality q-value `quant_q`. An accepted transfer gets `measurement_q`, the
//! running mean, over its run's accepted transfers in ascending probability of
//! a wrong measurement, of that probability (a model trained on the
//! pseudo-transfers: see `quant_finalize`), and its `quant_q` is the larger of
//! `transfer_q` and `measurement_q`; a transfer below the threshold keeps its
//! `transfer_q`. `quant_q` is so the largest error estimate made for the row:
//! a cut on it keeps a row only when each of its estimates passes the cut,
//! and is not the error rate of the rows it keeps.
//! With a quality threshold, an identified row above it (or without a
//! q-value) is `identified_below_threshold`. Each
//! pseudo-transfer is scored by the model of its precursor's fold, which never
//! saw it, and placed against that fold's transfer scores. `<base>.mbr.parquet`
//! gets one row per identified or transferred (precursor, run) target, and
//! `<base>.mbr_validation.parquet` the pseudo-transfer acceptance rate, the
//! wrong-peak rate (`|t0_pseudo - t0_identified| > w_identified`), the
//! wrong-measurement rate (a wrong peak, or a fitted quantity off by more
//! than a factor of 2 from the identified fit's; the peak alone where either
//! fit has no positive quantity) and the
//! accepted transfer counts of the pooled and per-family estimates. The
//! acceptance rate is optimistic for real transfers: the alignment map was
//! fitted with the pseudo-transfer's own run observation as an anchor, and
//! only the consensus RT leaves it out.
//!
//! Reported quantity: one relative elution width per precursor. Its relative
//! width `k` is the median, over its identified target rows with a fit of
//! positive amplitude, of `width_ratio`: the fitted width over the width
//! prior `w_prior(t)` of the row's run at the row's expected RT. A row with
//! fitted width `w` reports `A * w' * sqrt(2 pi)` and
//! `amplitude_se * w' * sqrt(2 pi)` with
//! `w' = k * (w / width_ratio) * (width_ratio / k)^lambda_run`: the relative
//! width times the width prior of the row's own run, aligned or not, moved
//! the share `lambda_run` of the way (in log scale) to the row's own fitted
//! width `w`; a transfer row keeps `lambda_run = 0` ([`row_width_weight`]).
//! `lambda_run` in `[0, 1]` is learned per run from the data,
//! without labels ([`run_width_weights`]): the normalisation reference cells
//! of a run take the weight that makes them agree best with the same
//! precursors in the other runs, so a run whose peaks are really narrower or
//! wider than its width prior says keeps that change in its quantities, and a
//! run whose width deviations are fit noise keeps the relative width alone. A
//! row of a precursor without such a fit reports the extract's own
//! `quantity` and `quantity_se`. The extract tables keep the per-run values.
//!
//! Whole isotope pattern: the fits use the first four isotopes of a
//! precursor, so the values above are the area of those four. The reported
//! `quantity` and `quantity_se` are these divided by the precursor's
//! `isotope_share` (library table): the share of the whole pattern the four
//! isotopes hold, from the same formula as the envelope. The extract tables,
//! the fits, the model features and the normalisation offsets keep the
//! four-isotope values.
//!
//! Trace area (`quant finalize --quant-trace-fraction F`, default 0.1; `--mbr`
//! uses the default): the brightest share `F` of the precursors, ranked by the
//! median amplitude
//! signal-to-noise `A / amplitude_se` of their identified target fits,
//! report the extract's template-free `trace_quantity` and
//! `trace_quantity_se` in every run instead (a row whose trace area is not
//! positive keeps the template quantity, decided per row). The choice is per
//! precursor, so a precursor's runs are measured alike apart from such rows. `quant_q` and the transfer q-values
//! do not use the reported quantity and are unchanged.
//!
//! Normalised quantity: `quantity_normalized` is the reported quantity with
//! its run's offset at the row's apex RT removed ([`quant::normalize`]: a
//! level per run and an RT-local term). The offsets are fitted on the
//! identified targets with `quant_q <= 0.01` and a positive quantity; a run
//! without enough of them reports null. `quantity` stays the measurement
//! without the offset.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use output::ms1_table::read_ms1_table;
use quant::align::RunMap;
use quant::claims;
use quant::decoys::{probes, Hypothesis, Probe, PEPTIDE_MASS_SLOPE};
use quant::evidence::{
    fit_many, fit_on_grid_gated, Evidence, EvidenceParams, FitJob, FitScratch, WIDTH_GRID,
};
use quant::library::{Library, LibraryEntry};
use quant::mobility::{MobilityObservation, RunMobility};
use quant::ms1_index::{MobilityGate, Ms1RunIndex};
use quant::ms2::{elution_window, LibrarySpectrum, Ms2RunIndex};
use quant::normalize::NormCell;
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
    /// Identified pairs measured on the peak at the position the other runs
    /// predict instead of the one at their PSM ([`prefer_consensus_peak`]).
    pub consensus_peak: usize,
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
    /// The run's learned ion-mobility scale, on a mobility-resolved run with
    /// enough identified precursors that have a library mobility.
    pub mobility: Option<RunMobility>,
    /// The run's learned claim windows (`None`: too few same-compound pairs,
    /// so no claims).
    pub claim_windows: Option<claims::Windows>,
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
        let mut out = format!(
            "ms1_spacing_s={spacing};width_scale_s={:.4};width_scale_fits={};width_scale_rounds={};width_scale_fallback={}",
            self.width_s, self.width_fits, self.rounds, self.fallback
        );
        if let Some(m) = &self.mobility {
            out.push_str(&format!(
                ";mobility_offset={:.5};mobility_sd={:.5};mobility_pairs={}",
                m.offset, m.sigma, m.pairs
            ));
        }
        match &self.claim_windows {
            Some(w) => out.push_str(&format!(
                ";claim_ppm={:.3};claim_s={:.3};claim_sd_ppm={:.4};claim_sd_s={:.4};claim_pairs={}",
                w.ppm, w.seconds, w.sd_ppm, w.sd_s, w.pairs
            )),
            None => out.push_str(";claim_windows=none"),
        }
        out
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
        self.consensus_peak += o.consensus_peak;
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
fn width_job(
    lib: &Library,
    i: usize,
    rt: f64,
    scale: f64,
    s: f64,
    mobility: Option<MobilityGate>,
) -> FitJob<'_> {
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
        mobility,
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
/// order, one at a time ([`fit_on_grid_gated`] gives the result [`fit_many`]
/// would, without its sweep over the whole run).
#[allow(clippy::too_many_arguments)]
fn width_scale<F>(
    ids: &[usize],
    run: u16,
    lib: &Library,
    ep: &EvidenceParams,
    series: &[&Ms1RunIndex],
    series_of: F,
    spacing_s: Option<f64>,
    gates: &[Option<MobilityGate>],
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
                let job = width_job(lib, i, rt, w, w.max(floor), gates[i]);
                fit_on_grid_gated(
                    series[series_of(i)],
                    &job.probe,
                    job.charge,
                    job.envelope,
                    ep,
                    job.widths,
                    job.mobility,
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
#[allow(clippy::too_many_arguments)]
fn width_fits<F>(
    ids: &[usize],
    run: u16,
    lib: &Library,
    ep: &EvidenceParams,
    series: &[&Ms1RunIndex],
    series_of: F,
    scales: &RunScales,
    gates: &[Option<MobilityGate>],
) -> Vec<(f64, f64)>
where
    F: Fn(usize) -> usize,
{
    let (w, s) = (scales.width_s, scales.width_s.max(scales.floor_s()));
    let jobs: Vec<(usize, FitJob)> = identified_targets(ids, run, lib)
        .into_iter()
        .map(|(i, rt)| (series_of(i), width_job(lib, i, rt, w, s, gates[i])))
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
    gates: &[Option<MobilityGate>],
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
                mobility: gates[i],
            };
            jobs.push((series_of(i), job));
        }
    }
    let mut evidence = fit_on_series(series, &jobs, ep).into_iter();
    for (&i, (rows, mut counts)) in ids.iter().zip(planned) {
        let mut rows: Vec<ExtractRow> = rows
            .into_iter()
            .map(|(mut row, _)| {
                row.evidence = evidence.next().flatten();
                row
            })
            .collect();
        counts.consensus_peak += usize::from(prefer_consensus_peak(&mut rows));
        slots[i] = Some((rows, counts));
    }
}

/// Amplitude of the fit at the position the other runs predict, as a
/// multiple of the identified fit's, above which an identified pair is
/// measured there (see the module docs). 2 is the factor beyond which two
/// measurements of one cell count as different ([`crate::quant_finalize`],
/// the wrong-measurement label).
pub(crate) const CONSENSUS_PEAK_FACTOR: f64 = 2.0;

/// Gives the identified target row of one entry's `rows` the fit of its
/// pseudo-transfer row when that fit is another, larger peak: both fits have
/// a positive amplitude, the apexes lie further apart than the sum of the
/// two fitted widths, and the pseudo-transfer's amplitude is more than
/// [`CONSENSUS_PEAK_FACTOR`] times the identified one. Returns whether it
/// did. The PSM of such a pair was taken off the elution apex (mostly on the
/// tail of an abundant peptide), where the identified fit, held near the
/// PSM, reads a fraction of the peak the other runs measure. The decoy rows
/// and the pseudo-transfer row keep their fits.
fn prefer_consensus_peak(rows: &mut [ExtractRow]) -> bool {
    let at = |role: Role| {
        rows.iter()
            .position(|r| r.role == role && r.kind == Hypothesis::Target)
    };
    let (Some(id), Some(pseudo)) = (at(Role::Identified), at(Role::PseudoTransfer)) else {
        return false;
    };
    let (Some(a), Some(b)) = (&rows[id].evidence, &rows[pseudo].evidence) else {
        return false;
    };
    let other_peak = (b.t0 - a.t0).abs() > a.w + b.w;
    let larger = a.amplitude > 0.0 && b.amplitude > CONSENSUS_PEAK_FACTOR * a.amplitude;
    if !(other_peak && larger) {
        return false;
    }
    rows[id].evidence = rows[pseudo].evidence.clone();
    true
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
            claim: claims::Claim::default(),
            shared: None,
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

/// Most identified precursors the run's mobility scale is learned from.
const MOBILITY_PAIRS: usize = 2000;

/// The run's learned ion-mobility scale ([`quant::mobility::learn`]) from up
/// to [`MOBILITY_PAIRS`] of its identified target entries with a library
/// mobility, spread evenly over the library, each searched within
/// [`SEED_WIDTH_S`] of its best PSM's retention time; `None` on a run that
/// is not mobility-resolved.
fn run_mobility(
    run: u16,
    lib: &Library,
    ep: &EvidenceParams,
    ms1: &Ms1RunIndex,
) -> Option<RunMobility> {
    if !ms1.has_mobility() {
        return None;
    }
    let ids: Vec<usize> = (0..lib.entries.len())
        .filter(|&i| lib.entries[i].ion_mobility.is_some())
        .collect();
    let all = identified_targets(&ids, run, lib);
    let stride = all.len().div_ceil(MOBILITY_PAIRS).max(1);
    let obs: Vec<MobilityObservation> = all
        .into_iter()
        .step_by(stride)
        .filter_map(|(i, rt)| {
            let e = &lib.entries[i];
            Some(MobilityObservation {
                mz: lib.meta[e.precursor as usize].mz,
                rt,
                k0: e.ion_mobility?,
            })
        })
        .collect();
    quant::mobility::learn(ms1, &obs, ep.tol_ppm, SEED_WIDTH_S)
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
    let mobility = run_mobility(run, lib, ep, &ms1);
    let gates: Vec<Option<MobilityGate>> = lib
        .entries
        .iter()
        .map(|e| mobility.and_then(|m| e.ion_mobility.map(|k0| m.gate(k0))))
        .collect();
    let groups = entry_groups(&ms1.faims_voltages(), &lib.entries);
    let whole_ids: Vec<usize> = (0..n).filter(|&i| groups[i] == Group::Whole).collect();
    let cv_ids: Vec<usize> = (0..n)
        .filter(|&i| matches!(groups[i], Group::Series(_)))
        .collect();
    let foreign = groups.iter().filter(|&&g| g == Group::Foreign).count();
    let mut slots: Vec<Option<EntryRows>> = (0..n).map(|_| None).collect();
    let finish = |slots, width_fits: usize, wp: &WidthPrior, scales: RunScales| {
        let mut out = assemble(slots, width_fits, wp.is_default());
        let windows = assign_claims(&mut out.rows, lib);
        out.counts.unaligned = !aligned;
        out.counts.foreign_cv = foreign;
        out.counts.scales = scales;
        out.counts.scales.claim_windows = windows;
        out
    };

    if cv_ids.is_empty() {
        let mut scales = width_scale(&whole_ids, run, lib, ep, &[&ms1], |_| 0, spacing, &gates);
        scales.mobility = mobility;
        let fits = width_fits(&whole_ids, run, lib, ep, &[&ms1], |_| 0, &scales, &gates);
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
            &gates,
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
    let mut scales = width_scale(&cv_ids, run, lib, ep, &series, series_of, spacing, &gates);
    scales.mobility = mobility;
    let fits = width_fits(&cv_ids, run, lib, ep, &series, series_of, &scales, &gates);
    let width_fit_count = fits.len();
    let wp = WidthPrior::from_fits(fits, &scales);
    rows_into(
        &mut slots, &cv_ids, run, aligned, lib, &gate, &wp, ep, &series, series_of, &gates,
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
            &gates,
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

/// What [`attach_shared`] did, for the extract log.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct SharedCounts {
    /// Transfer rows with a claimant and a fit.
    pub claimed: usize,
    /// Of them: covered by at least one MS2 scan, and scored against a
    /// claimant with distinguishing ions on both sides.
    pub with_scan: usize,
    pub scored: usize,
}

/// Peptidoform of a `"<peptidoform>/<charge>"` library key.
fn peptidoform_of(key: &str) -> &str {
    key.rsplit_once('/').map_or(key, |(p, _)| p)
}

/// Sets the MS2 distinguishing-ion evidence ([`quant::shared`]) of every
/// transfer target row of run `run` that has a claimant and a fit: the
/// scans of `ms2` covering its feature, identified through the donor PSMs of
/// the run (`donor_scans`, per precursor id), against the fragment ladders
/// of the row's precursor and its claimants (`mods` resolves modification
/// accessions). Other rows keep `None`.
pub(crate) fn attach_shared(
    rows: &mut [ExtractRow],
    lib: &Library,
    donor_scans: &[Vec<(u16, f64)>],
    mods: &HashMap<String, f64>,
    ms2: &Ms2RunIndex,
    tol: model::Tolerance,
) -> SharedCounts {
    use quant::ms2::fragment_ladder;
    use quant::shared::{evidence, Ions, Query, ScanIdentity};
    let wanted = |r: &ExtractRow| {
        r.role == Role::Transfer
            && r.kind == Hypothesis::Target
            && r.evidence.is_some()
            && r.claim.n > 0
    };
    let mut c = SharedCounts::default();
    for r in rows.iter_mut() {
        r.shared = None;
        c.claimed += usize::from(wanted(r));
    }
    if c.claimed == 0 {
        return c;
    }
    let run = rows[0].run;
    let mut peptides: HashMap<&str, u32> = HashMap::new();
    let mut intern = |p: usize| {
        let next = peptides.len() as u32;
        *peptides
            .entry(peptidoform_of(&lib.meta[p].key))
            .or_insert(next)
    };
    let mut identity = ScanIdentity::default();
    for (p, scans) in donor_scans.iter().enumerate() {
        let mut id = None;
        for &(r, rt) in scans {
            if r == run {
                let id = *id.get_or_insert_with(|| intern(p));
                identity.add(rt, id);
            }
        }
    }
    let mut ladders: HashMap<u32, Ions> = HashMap::new();
    for r in rows.iter().filter(|r| wanted(r)) {
        for p in std::iter::once(r.precursor).chain(r.claim.claimants()) {
            if ladders.contains_key(&p) {
                continue;
            }
            let meta = &lib.meta[p as usize];
            let peptide = intern(p as usize);
            let ions = match quant::ms2::residue_masses(peptidoform_of(&meta.key), mods) {
                Ok(m) => Ions::new(peptide, &fragment_ladder(&m, meta.charge.max(1))),
                Err(_) => Ions {
                    peptide,
                    mz: Vec::new(),
                },
            };
            ladders.insert(p, ions);
        }
    }
    let mut faims: HashMap<u32, Option<f32>> = HashMap::new();
    for e in &lib.entries {
        faims.insert(e.precursor, e.faims_cv);
    }
    rows.par_iter_mut().filter(|r| wanted(r)).for_each(|r| {
        let e = r.evidence.as_ref().expect("wanted rows have a fit");
        let claimants: Vec<&Ions> = r.claim.claimants().map(|p| &ladders[&p]).collect();
        let q = Query {
            mz: lib.meta[r.precursor as usize].mz,
            t0: e.t0,
            w: e.w,
            faims: faims.get(&r.precursor).copied().flatten(),
            target: &ladders[&r.precursor],
            claimants: &claimants,
        };
        r.shared = Some(Box::new(evidence(ms2, &identity, &q, tol)));
    });
    for r in rows.iter() {
        if let Some(s) = &r.shared {
            c.with_scan += usize::from(s.n_scans > 0);
            c.scored += usize::from(s.scored());
        }
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

/// Per precursor, the relative width `k` of the module docs, from the
/// identified target rows among `rows`.
#[cfg(test)]
pub(crate) fn relative_widths<'a>(
    rows: impl IntoIterator<Item = &'a ExtractRow>,
) -> HashMap<u32, f64> {
    let mut ratios: HashMap<u32, Vec<f64>> = HashMap::new();
    for r in rows {
        if let Some(k) = relative_fit_width(r) {
            ratios.entry(r.precursor).or_default().push(k);
        }
    }
    ratios.into_iter().map(|(p, k)| (p, median(k))).collect()
}

/// The log2 deviation `log2(width_ratio / k)` of evidence `e` from the
/// precursor's relative width `k = relative_width`, which the run width
/// weight multiplies in [`reported_quantity`]; `None` when that quantity does
/// not depend on the weight (a trace area with `trace`, no relative width,
/// or no finite positive width prior `w / width_ratio`), and 0 when the
/// deviation is not finite.
pub(crate) fn width_deviation(
    e: &Evidence,
    relative_width: Option<f64>,
    trace: bool,
) -> Option<f64> {
    if trace && e.trace_quantity > 0.0 && e.trace_quantity.is_finite() {
        return None;
    }
    let w_prior = e.w / e.width_ratio;
    let k = relative_width.filter(|_| w_prior.is_finite() && w_prior > 0.0)?;
    let d = (e.width_ratio / k).log2();
    Some(if d.is_finite() { d } else { 0.0 })
}

/// `(quantity, quantity_se)` of evidence `e`: with `trace`, its trace area
/// when that is positive; else at the width
/// `k * w_prior * 2^(own_weight * d)`: the precursor's relative width
/// `k = relative_width` times the width prior of the row's run, `w_prior =
/// w / width_ratio`, moved by the run width weight `own_weight` (in `[0, 1]`,
/// [`run_width_weights`]) along the row's [`width_deviation`] `d`, so that
/// weight 0 is `k * w_prior` and weight 1 the row's own fitted width `w`; or
/// the extract's own values without a relative width or without a finite
/// positive width prior.
pub(crate) fn reported_quantity(
    e: &Evidence,
    relative_width: Option<f64>,
    trace: bool,
    own_weight: f64,
) -> (f64, f64) {
    if trace && e.trace_quantity > 0.0 && e.trace_quantity.is_finite() {
        return (e.trace_quantity, e.trace_quantity_se);
    }
    match (relative_width, width_deviation(e, relative_width, false)) {
        (Some(k), Some(d)) => {
            let width = k * (e.w / e.width_ratio) * (own_weight * d).exp2();
            let area = width * (2.0 * std::f64::consts::PI).sqrt();
            (e.amplitude * area, e.amplitude_se * area)
        }
        _ => (e.quantity, e.quantity_se),
    }
}

/// A normalisation reference cell's part in [`run_width_weights`]: its
/// precursor and run, the log2 of its [`reported_quantity`] at weight 0,
/// and its [`width_deviation`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct WidthCell {
    pub precursor: u32,
    pub run: u16,
    pub log2_quantity: f64,
    pub deviation: f64,
}

/// The candidate run width weights are `0, 1/STEPS, ..., 1`.
const WIDTH_WEIGHT_STEPS: u32 = 10;
/// Fewest cells of a precursor for its cells to have a leave-one-out mean.
const WIDTH_WEIGHT_MIN_PRECURSOR_CELLS: usize = 4;
/// Fewest such cells of a run for its weight to be learned; a run with fewer
/// keeps weight 0.
const WIDTH_WEIGHT_MIN_RUN_CELLS: usize = 11;
/// Passes over the runs.
const WIDTH_WEIGHT_PASSES: usize = 2;
/// A larger weight replaces a smaller one only when it lowers the run's
/// spread (log2) by more than this.
const WIDTH_WEIGHT_TOLERANCE: f64 = 1e-6;

/// Median absolute deviation from the median of `v` (reordered).
fn median_absolute_deviation(v: &mut [f64]) -> f64 {
    let mid = |v: &mut [f64]| {
        v.sort_unstable_by(f64::total_cmp);
        let n = v.len();
        if n % 2 == 1 {
            v[n / 2]
        } else {
            0.5 * (v[n / 2 - 1] + v[n / 2])
        }
    };
    let m = mid(v);
    for x in v.iter_mut() {
        *x = (*x - m).abs();
    }
    mid(v)
}

/// Per run (`runs` of them, by run index), the weight of [`reported_quantity`]
/// learned from `cells` (one per precursor and run), without labels. A
/// cell's residual is its log2 quantity at its run's weight minus the mean of
/// its precursor's other cells (precursors with at least
/// [`WIDTH_WEIGHT_MIN_PRECURSOR_CELLS`] cells); a run's spread is the median
/// absolute deviation of its cells' residuals, which ignores the run's own
/// offset (the normalisation's part). Run by run, with the other runs at
/// their current weights, the weight of the smallest spread among
/// `0, 0.1, ..., 1` is taken (the smallest weight on a tie), in
/// [`WIDTH_WEIGHT_PASSES`] passes; a run's own weight does not move the means
/// its residuals are taken against, so each candidate costs one pass over
/// the run's cells. A run with fewer than [`WIDTH_WEIGHT_MIN_RUN_CELLS`] such
/// cells keeps 0, the precursor's relative width alone.
pub(crate) fn run_width_weights(cells: &[WidthCell], runs: usize) -> Vec<f64> {
    let mut weights = vec![0.0; runs];
    let mut order: Vec<usize> = (0..cells.len())
        .filter(|&i| (cells[i].run as usize) < runs)
        .collect();
    order.sort_by_key(|&i| (cells[i].precursor, cells[i].run));
    // Per cell its precursor's group; the group sizes and sums of the cells'
    // log2 quantities at the current weights.
    let mut group = vec![usize::MAX; cells.len()];
    let (mut size, mut sum) = (Vec::new(), Vec::new());
    for g in order.chunk_by(|&a, &b| cells[a].precursor == cells[b].precursor) {
        if g.len() < WIDTH_WEIGHT_MIN_PRECURSOR_CELLS {
            continue;
        }
        for &i in g {
            group[i] = size.len();
        }
        size.push(g.len() as f64);
        sum.push(g.iter().map(|&i| cells[i].log2_quantity).sum::<f64>());
    }
    let mut by_run: Vec<Vec<usize>> = vec![Vec::new(); runs];
    for &i in &order {
        if group[i] != usize::MAX {
            by_run[cells[i].run as usize].push(i);
        }
    }
    let (mut base, mut res) = (Vec::new(), Vec::new());
    for _ in 0..WIDTH_WEIGHT_PASSES {
        for (run, members) in by_run.iter().enumerate() {
            if members.len() < WIDTH_WEIGHT_MIN_RUN_CELLS {
                continue;
            }
            let current = weights[run];
            // Residual at weight 0: the cell minus its precursor's other cells.
            base.clear();
            base.extend(members.iter().map(|&i| {
                let (c, g) = (&cells[i], group[i]);
                let own = c.log2_quantity + current * c.deviation;
                c.log2_quantity - (sum[g] - own) / (size[g] - 1.0)
            }));
            let mut best = (0.0, f64::INFINITY);
            for step in 0..=WIDTH_WEIGHT_STEPS {
                let w = f64::from(step) / f64::from(WIDTH_WEIGHT_STEPS);
                res.clear();
                res.extend(
                    base.iter()
                        .zip(members)
                        .map(|(b, &i)| b + w * cells[i].deviation),
                );
                let spread = median_absolute_deviation(&mut res);
                if spread < best.1 - WIDTH_WEIGHT_TOLERANCE {
                    best = (w, spread);
                }
            }
            for &i in members {
                sum[group[i]] += (best.0 - current) * cells[i].deviation;
            }
            weights[run] = best.0;
        }
    }
    weights
}

/// A one-line description of run width weights `weights`: how many runs took
/// each weight.
pub(crate) fn describe_width_weights(weights: &[f64]) -> String {
    let mut counts = [0usize; WIDTH_WEIGHT_STEPS as usize + 1];
    for &w in weights {
        counts[(w * f64::from(WIDTH_WEIGHT_STEPS)).round() as usize] += 1;
    }
    counts
        .iter()
        .enumerate()
        .filter(|(_, &n)| n > 0)
        .map(|(s, n)| format!("{:.1}: {n}", s as f64 / f64::from(WIDTH_WEIGHT_STEPS)))
        .collect::<Vec<_>>()
        .join(", ")
}

/// The width weight of row `r` at run width weights `weights`: its run's
/// weight for an identified row, 0 for any other. A transfer's fitted width
/// is not confirmed by an identification of the peak, and transfers of low
/// abundance fit narrower peaks than their precursor's relative width, so a
/// transfer keeps the relative width.
pub(crate) fn row_width_weight(r: &ExtractRow, weights: &[f64]) -> f64 {
    if r.role == Role::Identified {
        weights.get(r.run as usize).copied().unwrap_or(0.0)
    } else {
        0.0
    }
}

/// The normalisation reference cells of `cells` (each with its
/// [`width_deviation`]) at the run width weights learned from them
/// ([`run_width_weights`], `runs` runs), and the weights.
pub(crate) fn weigh_widths(
    cells: Vec<(NormCell, Option<f64>)>,
    runs: usize,
) -> (Vec<NormCell>, Vec<f64>) {
    let width_cells: Vec<WidthCell> = cells
        .iter()
        .filter_map(|(c, d)| {
            d.map(|deviation| WidthCell {
                precursor: c.precursor,
                run: c.run,
                log2_quantity: c.log2_quantity,
                deviation,
            })
        })
        .collect();
    let weights = run_width_weights(&width_cells, runs);
    let cells = cells
        .into_iter()
        .map(|(mut c, d)| {
            if let (Some(d), Some(&w)) = (d, weights.get(c.run as usize)) {
                c.log2_quantity += w * d;
            }
            c
        })
        .collect();
    (cells, weights)
}

/// `(quantity, quantity_se)` of the whole isotope pattern, from the pair
/// measured on the library envelope, whose isotopes hold `share` of it.
pub(crate) fn whole_pattern((quantity, se): (f64, f64), share: f32) -> (f64, f64) {
    let share = f64::from(share);
    (quantity / share, se / share)
}

/// `amplitude / amplitude_se` of a fit. NaN when the standard error is not
/// positive (zero, negative or NaN), so a row with such a fit gets no
/// features; 0 with an infinite standard error (a fit without information).
pub(crate) fn amplitude_snr(e: &Evidence) -> f64 {
    if e.amplitude_se > 0.0 {
        e.amplitude / e.amplitude_se
    } else {
        f64::NAN
    }
}

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

/// The `width_ratio` of an identified target row, the input of the relative
/// width `k`; `None` for any other row, or without a positive amplitude, a
/// finite positive width or a finite positive ratio.
pub(crate) fn relative_fit_width(r: &ExtractRow) -> Option<f64> {
    if (r.role, r.kind) != (Role::Identified, Hypothesis::Target) {
        return None;
    }
    let e = r.evidence.as_ref()?;
    let fitted = e.amplitude > 0.0 && e.w.is_finite() && e.w > 0.0;
    (fitted && e.width_ratio.is_finite() && e.width_ratio > 0.0).then_some(e.width_ratio)
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
    let snr = amplitude_snr(e);
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
    let snr = amplitude_snr(e);
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
/// Sets the same-run claim of every fitted transfer target of `rows`: the
/// other fitted identified or transferred targets of the run on the same MS1
/// feature (see [`claims`]).
fn assign_claims(rows: &mut [ExtractRow], lib: &Library) -> Option<claims::Windows> {
    let mut probes = Vec::new();
    let mut transfers = Vec::new();
    let mut at = Vec::new();
    let mut peptides: HashMap<&str, u32> = HashMap::new();
    for (i, r) in rows.iter().enumerate() {
        if r.kind != Hypothesis::Target || !matches!(r.role, Role::Identified | Role::Transfer) {
            continue;
        }
        let Some(e) = &r.evidence else { continue };
        let meta = &lib.meta[r.precursor as usize];
        if r.role == Role::Transfer {
            transfers.push(probes.len());
            at.push(i);
        }
        let next = peptides.len() as u32;
        let peptide = *peptides
            .entry(
                meta.key
                    .rsplit_once('/')
                    .map_or(meta.key.as_str(), |(p, _)| p),
            )
            .or_insert(next);
        probes.push(claims::Probe {
            precursor: r.precursor,
            charge: meta.charge,
            mass: meta.mz * f64::from(meta.charge.max(1)),
            t0: e.t0,
            w: e.w,
            identified: r.role == Role::Identified,
            peptide,
            mass_err_ppm: e.mass_err_ppm,
        });
    }
    let windows = claims::learn_windows(&probes);
    for (i, claim) in at
        .into_iter()
        .zip(claims::assign(&probes, &transfers, windows))
    {
        rows[i].claim = claim;
    }
    windows
}

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
        fixture_psm_rt(faims, silent_in_run1, weak_in_run1, decoy, &|p, run| {
            apex(p, run) + 2.0
        })
    }

    /// [`fixture_decoys`] with the PSM of precursor `p` in run `run` at
    /// `psm_rt(p, run)`.
    fn fixture_psm_rt(
        faims: bool,
        silent_in_run1: &dyn Fn(u32) -> bool,
        weak_in_run1: &dyn Fn(u32) -> bool,
        decoy: &dyn Fn(u32) -> bool,
        psm_rt: &dyn Fn(u32, u16) -> f64,
    ) -> Fixture {
        let dir = tempfile::tempdir().unwrap();
        let meta: Vec<PrecursorMeta> = (0..N_PREC)
            .map(|p| PrecursorMeta {
                key: format!("PEP{p}/2"),
                mz: mz_of(p),
                charge: CHARGE,
                envelope: ENVELOPE,
                isotope_share: 1.0,
            })
            .collect();
        let mut psms = Vec::new();
        let mut psm = |p: u32, run: u16| {
            psms.push(PsmIn {
                q_gate: f64::NAN,
                run,
                precursor: p,
                rt: psm_rt(p, run),
                score: 10.0,
                q: if run == 1 && weak_in_run1(p) {
                    0.2
                } else {
                    0.001
                },
                pep: 0.0,
                is_decoy: decoy(p),
                faims_cv: (faims && run == 1).then_some(FAIMS_CVS[p as usize % 2]),
                ion_mobility: None,
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
                isotope_share: 1.0,
            })
            .collect();
        let mut psms = Vec::new();
        let mut psm = |p: u32, run: u16| {
            psms.push(PsmIn {
                q_gate: f64::NAN,
                run,
                precursor: p,
                rt: 100.0 + 8.0 * p as f64 + if run >= 2 { 120.0 } else { 0.0 },
                score: 10.0,
                q: 0.001,
                pep: 0.0,
                is_decoy: false,
                faims_cv: None,
                ion_mobility: None,
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

    /// Run 1 of [`fixture`] with 1 s wide peaks sampled every 0.25 s and a
    /// same-m/z peak of `echo` times the amplitude [`ECHO_AFTER`] seconds
    /// after each apex, written next to the fixture's tables.
    fn echo_run(fx: &Fixture, echo: f64) -> PathBuf {
        let path = fx._dir.path().join("echo.ms1.parquet");
        let index = Ms1RunIndex::new(run_scans_shaped(1, false, &|_| false, 1.0, 0.25, echo));
        write_ms1_table(&path, "run1", &index).unwrap();
        path
    }

    /// Whether the run-1 PSM of precursor `p` sits on the echo in
    /// [`fixture_late_psms`]: three precursors apart from each other and
    /// from the ends of the gradient, so the run's map (30 anchors) still
    /// follows the other PSMs.
    fn late(p: u32) -> bool {
        matches!(p, 8 | 15 | 22)
    }

    /// [`fixture`] where the run-1 PSMs of the [`late`] precursors were
    /// taken [`ECHO_AFTER`] seconds after the apex (run 0's PSMs are at the
    /// apex).
    fn fixture_late_psms() -> Fixture {
        fixture_psm_rt(false, &|_| false, &|_| false, &|_| false, &|p, run| {
            let after = if run == 1 && late(p) { ECHO_AFTER } else { 0.0 };
            apex(p, run) + 2.0 + after
        })
    }

    #[test]
    fn identified_pair_with_a_psm_on_the_tail_is_measured_on_the_elution_peak() {
        // The run-1 PSMs of the late precursors sit on a peak of a twentieth
        // of the amplitude, 19 s after the apex that run 0 predicts.
        let fx = fixture_late_psms();
        let path = echo_run(&fx, 0.05);
        let out = extract_all(1, &path, &fx.lib, &params()).unwrap();
        assert_eq!(out.counts.consensus_peak, 3);
        for p in 0..BOTH {
            let rows = rows_of(&out.rows, p);
            let (id, pseudo) = (rows[0], rows[3]);
            assert!(is_identified_target(id));
            assert_eq!(pseudo.role, Role::PseudoTransfer);
            let e = id.evidence.as_ref().unwrap();
            assert!((e.t0 - apex(p, 1)).abs() <= 0.5, "precursor {p}: {e:?}");
            assert!((e.amplitude / AMPLITUDE - 1.0).abs() < 0.1, "{p}: {e:?}");
            if late(p) {
                assert_eq!(id.evidence, pseudo.evidence, "precursor {p}");
            }
            // The decoy probes keep the fits of the PSM-centred prior.
            for d in &rows[1..3] {
                assert_eq!(d.role, Role::Identified);
                assert_ne!(d.evidence, id.evidence);
            }
        }
    }

    #[test]
    fn identified_pair_keeps_a_psm_peak_that_is_not_much_lower() {
        // The same late PSMs on a peak of 0.6 times the amplitude at the
        // predicted position: less than the factor, so the PSM's peak stays.
        let fx = fixture_late_psms();
        let path = echo_run(&fx, 0.6);
        let out = extract_all(1, &path, &fx.lib, &params()).unwrap();
        assert_eq!(out.counts.consensus_peak, 0);
        for p in (0..BOTH).filter(|&p| late(p)) {
            let rows = rows_of(&out.rows, p);
            let e = rows[0].evidence.as_ref().unwrap();
            assert!(
                (e.t0 - apex(p, 1) - ECHO_AFTER).abs() <= 0.5,
                "precursor {p}: {e:?}"
            );
            assert!((e.amplitude / (0.6 * AMPLITUDE) - 1.0).abs() < 0.1, "{e:?}");
        }
    }

    #[test]
    fn consensus_peak_needs_another_larger_peak_and_a_signal_at_the_psm() {
        let fx = fixture(false);
        let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
        let base: Vec<ExtractRow> = rows_of(&out.rows, 3).into_iter().cloned().collect();
        let fit = base[0].evidence.clone().unwrap();
        // Rows whose identified fit is `id` and whose pseudo-transfer fit
        // is `pseudo`, and whether the identified row then takes the
        // pseudo-transfer's fit.
        let taken = |id: Option<&Evidence>, pseudo: Option<&Evidence>| {
            let mut rows = base.clone();
            rows[0].evidence = id.cloned();
            rows[3].evidence = pseudo.cloned();
            let before = rows.clone();
            let took = prefer_consensus_peak(&mut rows);
            if took {
                assert_eq!(rows[0].evidence, rows[3].evidence);
                assert_eq!(rows[1..], before[1..]);
            } else {
                assert_eq!(rows, before);
            }
            took
        };
        let at = |t0: f64, amplitude: f64| Evidence {
            t0,
            amplitude,
            w: 4.0,
            ..fit.clone()
        };
        let id = at(500.0, 1e4);
        // Another peak (apexes more than the two widths apart), more than
        // twice the amplitude.
        assert!(taken(Some(&id), Some(&at(440.0, 2.1e4))));
        assert!(taken(Some(&id), Some(&at(508.1, 1e6))));
        // Twice the amplitude or less.
        assert!(!taken(Some(&id), Some(&at(440.0, 2e4))));
        assert!(!taken(Some(&id), Some(&at(440.0, 5e3))));
        // The same peak: apexes within the two widths.
        assert!(!taken(Some(&id), Some(&at(508.0, 1e6))));
        assert!(!taken(Some(&id), Some(&at(500.0, 1e6))));
        // No signal at the PSM, or a missing fit on either side.
        assert!(!taken(Some(&at(500.0, 0.0)), Some(&at(440.0, 1e6))));
        assert!(!taken(None, Some(&at(440.0, 1e6))));
        assert!(!taken(Some(&id), None));
        // Without a pseudo-transfer row there is nothing to take.
        let mut rows = base[..3].to_vec();
        assert!(!prefer_consensus_peak(&mut rows));
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
        write_mbr_validation, Accepted, EvidenceType, MbrClaim, MbrParquetWriter, MbrRow,
        MbrValidation, VALIDATION_Q,
    };
    use quant::align::AlignParams;
    use quant::library::{build_library, PrecursorMeta, PsmIn};
    use quant::normalize::{fit as fit_normalization, NormParams};

    use super::*;

    /// One output row of the reference finalize: precursor, run, the table row,
    /// the fit's mass error and the row's claim.
    type OutEntry<'a> = (u32, u16, MbrRow<'a>, (Option<f64>, claims::Claim));
    use crate::quant_finalize::{
        finalize, finalize_stream, for_each_scored, mismatch, pseudo_example, running_mean_q,
        table_features, transfer_quant_q, wrong_features, wrong_label, BlockContext, EntryIndex,
        ExtractRuns, FinalizeLimits, FinalizeParams, FinalizeSummary, MbrFeature, Mismatch, Table,
        WrongCalibration, SCORE_CHUNK,
    };
    use crate::transfer_fdr::{
        row_unit, train_bagged_folds, train_folds, transfer_train_params, Estimate, FdrRow,
        FoldScores, QValues, TrainingSample, TARGET,
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
            trace_quantity: amplitude * 12.5,
            trace_quantity_se: 1.25e5,
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
            claim: claims::Claim::default(),
            shared: None,
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
                isotope_share: 1.0,
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
                num(r.measurement_q),
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
                ion_mobility: None,
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
        // The fixture's quantities agree, so a wrong measurement is a wrong
        // peak.
        assert_eq!(v.wrong_measurement, v.wrong_peak);
        assert_eq!(v.at_q01.wrong_measurement, v.at_q01.wrong_peak);
        assert!(v.pseudo_acceptance(&v.at_q01).unwrap() > 0.9, "{v:?}");
        assert_eq!(v.at_q01, v.at_threshold);
        assert_eq!(v.transfers, N as usize);
        assert_eq!(v.transfers_scored, N as usize - 1);
        assert!(v.at_q01.mass_only > 0 && v.at_q01.rt_only > 0);
        let t = table(&val);
        assert_eq!(t.len(), 2);
        let get = |name: &str| &t[1][t[0].iter().position(|h| h == name).unwrap()];
        assert_eq!(get("wrong_peak_rate"), "0.1");
        assert_eq!(get("wrong_measurement_rate"), "0.1");
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
                q_gate: f64::NAN,
                run: 0,
                precursor: p,
                rt: 100.0 + p as f64,
                score: 10.0 + p as f64,
                q: 0.0001,
                pep: 0.0001,
                is_decoy: false,
                faims_cv: None,
                ion_mobility: None,
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
        let t_models =
            train_bagged_folds(&[&t_sample.into_items()], seed, &transfer_train_params());
        let i_models = train_folds(&i_sample.into_items(), seed);
        // The wrong-peak model: the pseudo-transfers judged against the
        // identified fit of the same cell, and its calibration on held-out
        // margins.
        let mut w_sample: TrainingSample<u64> = TrainingSample::new(usize::MAX, seed);
        let mut w_seq = 0u64;
        shards.for_each(|rows| {
            let ctx = BlockContext::new(rows, maps)?;
            for r in rows {
                w_seq += 1;
                if let Some((y, x)) = pseudo_example(r, &ctx, &idx) {
                    w_sample.push(
                        w_seq,
                        row_unit(r.precursor, r.run),
                        FdrRow {
                            precursor: r.precursor,
                            is_decoy: y,
                            family: 0,
                            x,
                        },
                    );
                }
            }
            Ok(())
        })?;
        let w_models = train_folds(&w_sample.into_items(), seed);
        let mut w_held: Vec<(f32, bool)> = Vec::new();
        shards.for_each(|rows| {
            let ctx = BlockContext::new(rows, maps)?;
            let ex: Vec<(bool, Vec<f32>, u32)> = rows
                .iter()
                .filter_map(|r| pseudo_example(r, &ctx, &idx).map(|(y, x)| (y, x, r.precursor)))
                .collect();
            let precursors: Vec<u32> = ex.iter().map(|e| e.2).collect();
            let xs: Vec<&[f32]> = ex.iter().map(|e| e.1.as_slice()).collect();
            let mut margins = Vec::new();
            w_models.score(&precursors, &xs, &mut margins);
            w_held.extend(
                margins
                    .iter()
                    .zip(&ex)
                    .filter_map(|(m, e)| m.map(|m| (m, e.0))),
            );
            Ok(())
        })?;
        let wrong_cal = WrongCalibration::fit(w_held);

        // Pass 2: every row's held-out score, into the q-value step functions.
        let mut t_scores = FoldScores::with_strata(seed, idx.strata());
        let mut i_scores = FoldScores::new(seed);
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
                |i, t, fam, s| {
                    let p = rows[i].precursor;
                    match t {
                        Table::Identified => i_scores.add(p, fam, s),
                        _ => t_scores.add(p, fam, s),
                    }
                },
            );
            Ok(())
        })?;
        let fdr = t_scores.finish(&t_models.failed_folds);
        let id_fdr = i_scores.finish(&i_models.failed_folds);

        // Pass 3: the reference cells of the run normalisation.
        let mut cells: TrainingSample<u64, (NormCell, Option<f64>)> =
            TrainingSample::new(usize::MAX, seed);
        let mut cell_seq = 0u64;
        let empty_shapes = DonorShapes::default();
        shards.for_each(|rows| {
            let k = relative_widths(rows.iter());
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
                            .q(Estimate::Pooled, r.precursor, fam, s)
                            .is_some_and(|q| q <= VALIDATION_Q)
                    {
                        return;
                    }
                    let Some(e) = r.evidence.as_ref() else { return };
                    let k = k.get(&r.precursor).copied();
                    let q = reported_quantity(e, k, false, 0.0).0;
                    if q > 0.0 && q.is_finite() {
                        cell_seq += 1;
                        cells.push(
                            cell_seq,
                            u64::from(r.precursor),
                            (
                                NormCell {
                                    precursor: r.precursor,
                                    run: r.run,
                                    log2_quantity: q.log2(),
                                    rt: e.t0,
                                },
                                width_deviation(e, k, false),
                            ),
                        );
                    }
                },
            );
            Ok(())
        })?;
        let (cells, weights) = weigh_widths(cells.into_items(), stems.len());
        let norm = fit_normalization(&cells, stems.len(), &NormParams::default());

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
            let k = relative_widths(rows.iter());
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
            let q_of = |q: &QValues, e: Estimate, i: usize| {
                scored[i].and_then(|(fam, s)| q.q(e, rows[i].precursor, fam, s))
            };
            // A fit whose amplitude is clamped to 0 found no signal: its
            // quantity is missing, not zero, in the MBR table and the
            // feature parquet.
            let reported = |r: &ExtractRow| {
                r.evidence.as_ref().filter(|e| e.amplitude > 0.0).map(|e| {
                    whole_pattern(
                        reported_quantity(
                            e,
                            k.get(&r.precursor).copied(),
                            false,
                            row_width_weight(r, &weights),
                        ),
                        lib.meta[r.precursor as usize].isotope_share,
                    )
                })
            };
            let mut out: Vec<OutEntry<'_>> = Vec::new();
            // Probability that a transfer's measurement is wrong (1 without
            // features, a fold model or a calibration).
            let mut p_wrong: Vec<f64> = vec![1.0; rows.len()];
            {
                let at: Vec<usize> = (0..rows.len())
                    .filter(|&i| {
                        rows[i].role == Role::Transfer
                            && matches!(idx.table(&rows[i]), Some((Table::Transfer, TARGET)))
                    })
                    .collect();
                let xs: Vec<Option<Vec<f32>>> = at
                    .iter()
                    .map(|&i| rows[i].evidence.as_ref().and_then(wrong_features))
                    .collect();
                let ok: Vec<usize> = (0..at.len()).filter(|&k| xs[k].is_some()).collect();
                let precursors: Vec<u32> = ok.iter().map(|&k| rows[at[k]].precursor).collect();
                let x_refs: Vec<&[f32]> = ok.iter().map(|&k| xs[k].as_deref().unwrap()).collect();
                let mut margins = Vec::new();
                w_models.score(&precursors, &x_refs, &mut margins);
                for (&k, m) in ok.iter().zip(margins) {
                    if let (Some(m), Some(c)) = (m, &wrong_cal) {
                        p_wrong[at[k]] = c.at(m);
                    }
                }
            }
            let mut wrongs: Vec<(usize, f64)> = Vec::new();
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
                        let q = scored[i].and_then(|(fam, s)| fdr.pseudo_q(r.precursor, fam, s));
                        v.pseudo_transfers += 1;
                        if let Some(m) = &r.ms2 {
                            v.ms2_pseudo_rows += 1;
                            v.ms2_pseudo_supported += usize::from(m.supports());
                        }
                        v.pseudo_scored += usize::from(q.is_some());
                        let wrong = e
                            .zip(identified.get(&(r.precursor, r.run)))
                            .and_then(|(p, id)| mismatch(p, id));
                        if let Some(m) = wrong {
                            v.wrong_peak_pairs += 1;
                            v.wrong_peak += usize::from(m.apex);
                            v.wrong_measurement += usize::from(m.wrong());
                        }
                        for (a, &t) in at.iter_mut().zip(&thresholds) {
                            if !q.is_some_and(|q| q <= t) {
                                continue;
                            }
                            a.pseudo += 1;
                            if let Some(m) = wrong {
                                a.pairs += 1;
                                a.wrong_peak += usize::from(m.apex);
                                a.wrong_measurement += usize::from(m.wrong());
                            }
                        }
                        continue;
                    }
                };
                if r.role == Role::Transfer {
                    wrongs.push((out.len(), p_wrong[i]));
                }
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
                        claim: (r.role == Role::Transfer).then(|| MbrClaim {
                            claimed_by: r.claim.claimant.map(|p| lib.meta[p as usize].key.as_str()),
                            dt_s: r.claim.dt_s.map(f64::from),
                            identified: r.claim.identified,
                            co_apex: r.claim.co_apex,
                            n_co_apex: u32::from(r.claim.n_co_apex),
                            n_identified_co_apex: u32::from(r.claim.n_identified_co_apex),
                        }),
                        measurement_q: None,
                    },
                    (
                        e.map(|e| e.mass_err_ppm),
                        if r.role == Role::Transfer {
                            r.claim
                        } else {
                            claims::Claim::default()
                        },
                    ),
                ));
            }
            // quant_q of an accepted transfer: the running mean, over its
            // run's accepted transfers in ascending probability of a wrong
            // measurement, of that probability.
            let mut by_run: HashMap<u16, Vec<(f64, usize)>> = HashMap::new();
            for &(k, p) in &wrongs {
                if matches!(
                    out[k].2.evidence_type,
                    EvidenceType::TransferMs2 | EvidenceType::TransferMs1
                ) {
                    by_run.entry(out[k].1).or_default().push((p, k));
                }
            }
            for pool in by_run.values_mut() {
                running_mean_q(pool, |k, q| {
                    let r = &mut out[k].2;
                    r.measurement_q = Some(q);
                    r.quant_q = Some(transfer_quant_q(r.transfer_q, q));
                });
            }
            out.sort_by_key(|&(p, run, ..)| (p, run));
            let mut features = Vec::with_capacity(out.len());
            for (p, run, r, (mass_err, claim)) in &out {
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
                    claim: *claim,
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
            trace_fraction: 0.0,
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

    /// A precursor seen in three runs with the same abundance, eluting
    /// `K_TRUE` times as wide as its run's usual peak: run `r`'s width prior
    /// is `5 * scale_r`, the apex `A0 / scale_r`, and the fitted width is off
    /// by `jitter[r]` while the amplitude is right.
    fn jittered_rows(p: u32, jitter: [f64; 3]) -> Vec<ExtractRow> {
        let (a0, w0) = (1e6, 5.0);
        let scales = [1.0, 1.25, 2.0];
        let mut g = Rng(p as u64 + 1);
        (0..3u16)
            .flat_map(|run| {
                let scale = scales[run as usize];
                let amplitude = a0 / scale;
                let mut e = evidence(&mut g, true, 500.0 * scale);
                e.amplitude = amplitude;
                e.amplitude_se = 0.01 * amplitude;
                e.width_ratio = K_TRUE * jitter[run as usize];
                e.w = w0 * scale * e.width_ratio;
                e.quantity = amplitude * e.w * (2.0 * std::f64::consts::PI).sqrt();
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

    /// Width of the precursor of [`jittered_rows`] over its run's width prior.
    const K_TRUE: f64 = 1.2;

    #[test]
    fn relative_width_gives_one_quantity_across_runs() {
        let truth = 1e6 * 5.0 * K_TRUE * (2.0 * std::f64::consts::PI).sqrt();
        for jitter in [[0.7, 1.4, 1.0], [1.4, 0.7, 1.4], [1.0, 1.0, 0.7]] {
            let rows = jittered_rows(3, jitter);
            let k = relative_widths(&rows);
            let per_run: Vec<(f64, f64)> = rows
                .iter()
                .filter(|r| r.role == Role::Identified)
                .map(|r| {
                    let e = r.evidence.as_ref().unwrap();
                    reported_quantity(e, k.get(&3).copied(), false, 0.0)
                })
                .collect();
            // The three runs, of different peak widths, report one quantity:
            // the jitter is gone.
            for &(q, se) in &per_run {
                assert!((q / per_run[0].0 - 1.0).abs() < 1e-12, "{per_run:?}");
                assert!((se / q - 0.01).abs() < 1e-12);
            }
            // The relative width is the median of the three runs' ratios.
            let mut sorted = jitter;
            sorted.sort_by(f64::total_cmp);
            let want = truth * sorted[1];
            assert!((per_run[0].0 / want - 1.0).abs() < 1e-12, "{per_run:?}");
        }
        // Raw quantities with jitter 0.7 vs 1.4 differ twofold at equal
        // abundance; the relative width removes it.
        let rows = jittered_rows(3, [0.7, 1.4, 1.0]);
        let raw = |run: usize| rows[2 * run].evidence.as_ref().unwrap().quantity;
        assert!((raw(1) / raw(0) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn relative_width_skips_unfitted_rows() {
        let mut rows = jittered_rows(4, [1.0, 0.7, 1.4]);
        // Only run 0's identified fit counts: run 1's found no signal and
        // run 2's has no width ratio.
        rows[2].evidence.as_mut().unwrap().amplitude = 0.0;
        rows[4].evidence.as_mut().unwrap().width_ratio = f64::NAN;
        let k = relative_widths(&rows);
        assert!((k[&4] - K_TRUE).abs() < 1e-12, "{k:?}");
        // A row without a width prior keeps the extract's own quantity.
        let e = rows[4].evidence.as_ref().unwrap();
        assert_eq!(
            reported_quantity(e, Some(k[&4]), false, 0.0),
            (e.quantity, e.quantity_se)
        );
        // So does a precursor without an identified fit.
        let mut rows = jittered_rows(5, [1.0, 1.0, 1.0]);
        for r in rows.iter_mut().filter(|r| r.role == Role::Identified) {
            r.evidence = None;
        }
        assert!(relative_widths(&rows).is_empty());
        let e = rows[1].evidence.as_ref().unwrap();
        assert_eq!(
            reported_quantity(e, None, false, 0.0),
            (e.quantity, e.quantity_se)
        );
    }

    #[test]
    fn width_weight_moves_the_quantity_to_the_own_fitted_width() {
        let rows = jittered_rows(3, [0.7, 1.4, 1.0]);
        let k = relative_widths(&rows)[&3];
        let e = rows[0].evidence.as_ref().unwrap();
        let at = |w: f64| reported_quantity(e, Some(k), false, w);
        // Weight 1 is the row's own fitted area, A * w * sqrt(2 pi).
        assert!((at(1.0).0 / e.quantity - 1.0).abs() < 1e-12);
        assert!((at(1.0).1 / e.quantity_se - 1.0).abs() < 1e-12);
        // Weight 0.5 is the geometric mean of weight 0 and the own area.
        let mid = (at(0.0).0 * at(1.0).0).sqrt();
        assert!((at(0.5).0 / mid - 1.0).abs() < 1e-12);
        let d = width_deviation(e, Some(k), false).unwrap();
        assert!((d - (0.7f64 / 1.0).log2()).abs() < 1e-12, "{d}");
        // No relative width, or a trace area: the weight has nothing to move.
        assert_eq!(width_deviation(e, None, false), None);
        let mut t = e.clone();
        t.trace_quantity = 3.0;
        assert_eq!(width_deviation(&t, Some(k), true), None);
        assert_eq!(reported_quantity(&t, Some(k), true, 1.0).0, 3.0);
    }

    /// Cells of `precursors` precursors in four runs: run 0's peaks really
    /// change width (its weight-0 quantities miss the change `d`), run 1's
    /// width deviations are fit noise, runs 2 and 3 have none.
    fn width_cells(precursors: u32) -> Vec<WidthCell> {
        let mut g = Rng(11);
        let mut cells = Vec::new();
        for p in 0..precursors {
            let truth = 20.0 + 3.0 * g.normal();
            for run in 0..4u16 {
                let d = if run < 2 { 0.4 * g.normal() } else { 0.0 };
                let missed = if run == 0 { d } else { 0.0 };
                cells.push(WidthCell {
                    precursor: p,
                    run,
                    log2_quantity: truth - missed + 0.05 * g.normal() + f64::from(run),
                    deviation: d,
                });
            }
        }
        cells
    }

    #[test]
    fn run_width_weights_keep_real_width_changes_and_drop_fit_noise() {
        let cells = width_cells(300);
        // The run offsets (+run) are the normalisation's part, not the weights'.
        assert_eq!(run_width_weights(&cells, 4), vec![1.0, 0.0, 0.0, 0.0]);
        // A run without cells keeps 0.
        assert_eq!(run_width_weights(&cells, 5), vec![1.0, 0.0, 0.0, 0.0, 0.0]);
        // The order of the cells does not matter.
        let mut reversed = cells.clone();
        reversed.reverse();
        assert_eq!(run_width_weights(&reversed, 4), vec![1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn run_width_weights_need_enough_cells() {
        // Fewer than WIDTH_WEIGHT_MIN_RUN_CELLS precursors: every run keeps 0.
        let few = width_cells(WIDTH_WEIGHT_MIN_RUN_CELLS as u32 - 1);
        assert_eq!(run_width_weights(&few, 4), vec![0.0; 4]);
        // Precursors in fewer than WIDTH_WEIGHT_MIN_PRECURSOR_CELLS runs have
        // no leave-one-out mean.
        let three: Vec<WidthCell> = width_cells(300).into_iter().filter(|c| c.run < 3).collect();
        assert_eq!(run_width_weights(&three, 4), vec![0.0; 4]);
        assert_eq!(run_width_weights(&[], 2), vec![0.0; 2]);
    }

    #[test]
    fn only_identified_rows_take_the_run_width_weight() {
        let rows = jittered_rows(3, [0.7, 1.4, 1.0]);
        let weights = [0.6, 0.3, 1.0];
        for r in &rows {
            let want = if r.role == Role::Identified {
                weights[r.run as usize]
            } else {
                0.0
            };
            assert_eq!(row_width_weight(r, &weights), want);
        }
        // A run beyond the weights keeps 0.
        assert_eq!(row_width_weight(&rows[0], &[]), 0.0);
        let mut transfer = rows[0].clone();
        transfer.role = Role::Transfer;
        assert_eq!(row_width_weight(&transfer, &weights), 0.0);
    }

    #[test]
    fn weigh_widths_moves_only_cells_with_a_deviation() {
        let cells: Vec<(NormCell, Option<f64>)> = width_cells(300)
            .into_iter()
            .enumerate()
            .map(|(i, c)| {
                let cell = NormCell {
                    precursor: c.precursor,
                    run: c.run,
                    log2_quantity: c.log2_quantity,
                    rt: 100.0,
                };
                // Every tenth cell reports a trace area: no deviation.
                (cell, (i % 10 != 0).then_some(c.deviation))
            })
            .collect();
        let (out, weights) = weigh_widths(cells.clone(), 4);
        assert_eq!(weights, vec![1.0, 0.0, 0.0, 0.0]);
        for ((before, d), after) in cells.iter().zip(&out) {
            let want = before.log2_quantity + d.map_or(0.0, |d| weights[before.run as usize] * d);
            assert_eq!(after.log2_quantity, want);
            assert_eq!(
                (after.precursor, after.run, after.rt),
                (before.precursor, before.run, before.rt)
            );
        }
        assert_eq!(describe_width_weights(&weights), "0.0: 3, 1.0: 1");
    }

    #[test]
    fn trace_fraction_reports_the_trace_area_of_the_brightest_precursors() {
        use crate::quant_finalize::{ExtractFiles, ExtractSource};
        let lib = library();
        let mut runs = extract();
        // The identified signal-to-noise rises with the precursor index.
        for r in runs[0].iter_mut().filter(|r| r.role == Role::Identified) {
            let e = r.evidence.as_mut().unwrap();
            e.amplitude_se = e.amplitude / (10.0 + r.precursor as f64);
            if r.precursor == 0 {
                // A non-positive trace area keeps the template quantity.
                e.trace_quantity = 0.0;
            }
        }
        let dir = tempfile::tempdir().unwrap();
        // The same rows from memory and from extract tables on disk: finalize
        // reads the identified fits of the tables with a reader of its own.
        let files = ExtractFiles(
            runs.iter()
                .enumerate()
                .map(|(r, rows)| {
                    let path = dir.path().join(format!("{r}.extract.parquet"));
                    output::extract_table::write_extract_table(&path, rows).unwrap();
                    path
                })
                .collect(),
        );
        let memory = ExtractRuns {
            runs: &runs,
            batch_rows: usize::MAX,
        };
        let sources: [(&str, &dyn ExtractSource); 2] = [("memory", &memory), ("files", &files)];
        let target: HashMap<(u32, u16), &Evidence> = runs
            .iter()
            .flatten()
            .filter(|r| r.kind == Hypothesis::Target && r.role != Role::PseudoTransfer)
            .filter_map(|r| r.evidence.as_ref().map(|e| ((r.precursor, r.run), e)))
            .collect();
        for (name, source) in sources {
            let finalize_at = |fraction: f64| {
                let (mbr, val) = mbr_paths(&dir.path().join(format!("{name}{fraction}")));
                let params = FinalizeParams {
                    mbr_fdr: 0.01,
                    quant_feature_fdr: None,
                    seed: 7,
                    limits: FinalizeLimits::default(),
                    mbr_path: &mbr,
                    validation_path: &val,
                    trace_fraction: fraction,
                };
                let mut features = Vec::new();
                finalize_stream(&stems(), &lib, source, &params, &mut |f| {
                    features.extend(f);
                    Ok(())
                })
                .unwrap();
                features
            };
            let off = finalize_at(0.0);
            let all = finalize_at(1.0);
            let half = finalize_at(0.5);
            assert_eq!(off.len(), all.len(), "{name}");
            assert_eq!(off.len(), half.len(), "{name}");
            let mut traced = 0;
            for ((o, a), h) in off.iter().zip(&all).zip(&half) {
                let key = (o.precursor, o.run);
                assert_eq!(key, (a.precursor, a.run));
                assert_eq!(key, (h.precursor, h.run));
                // The q-values and the evidence type do not use the quantity.
                for f in [a, h] {
                    assert_eq!(
                        (f.evidence_type, f.transfer_q, f.quant_q),
                        (o.evidence_type, o.transfer_q, o.quant_q),
                        "{name}"
                    );
                }
                let Some(q) = o.quantity else {
                    assert!(a.quantity.is_none() && h.quantity.is_none());
                    continue;
                };
                let e = target[&key];
                assert_ne!(q, e.trace_quantity);
                if e.trace_quantity <= 0.0 {
                    for f in [a, h] {
                        assert_eq!(
                            (f.quantity, f.quantity_se),
                            (o.quantity, o.quantity_se),
                            "{name} {key:?}"
                        );
                    }
                    continue;
                }
                // Every precursor has an identified fit: share 1 traces them all.
                assert_eq!(
                    (a.quantity, a.quantity_se),
                    (Some(e.trace_quantity), Some(e.trace_quantity_se)),
                    "{name} {key:?}"
                );
                // Share 0.5 traces the brighter half, precursors 450 and up.
                if o.precursor >= N / 2 {
                    assert_eq!(h.quantity, Some(e.trace_quantity), "{name} {key:?}");
                    traced += 1;
                } else {
                    assert_eq!((h.quantity, h.quantity_se), (o.quantity, o.quantity_se));
                }
            }
            assert!(traced > 0, "{name}");
        }
    }

    #[test]
    fn finalize_reports_the_quantity_at_the_relative_width() {
        // Precursor 0 is identified in run 0 (width prior 5) at 0.8 of the
        // prior and in run 1 (width prior 10) at 1.2: its relative width is
        // 1, and each row reports its A * the run's width prior * sqrt(2 pi).
        let lib = library();
        let mut runs = extract();
        let mut g = Rng(9);
        for (run, amplitude, w, ratio) in [(0u16, 1e6, 4.0, 0.8), (1, 5e5, 12.0, 1.2)] {
            let mut e = evidence(&mut g, true, 100.0);
            e.amplitude = amplitude;
            e.amplitude_se = 1e3;
            e.w = w;
            e.width_ratio = ratio;
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
        assert!((a.quantity.unwrap() / (1e6 * 5.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
        assert!((b.quantity.unwrap() / (5e5 * 10.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
        assert!((a.quantity_se.unwrap() / (1e3 * 5.0 * sqrt_2pi) - 1.0).abs() < 1e-12);
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
        // A transfer below the threshold keeps its transfer_q as quant_q; an
        // accepted one has the larger of its transfer_q and its
        // measurement_q, the running mean of its run's wrong-measurement
        // probabilities. No other row has a measurement_q.
        for f in s.features.iter().filter(|f| f.run == 1) {
            if f.evidence_type == EvidenceType::BelowThreshold {
                assert_eq!(f.quant_q, f.transfer_q);
            } else {
                assert!(f.quant_q.is_some_and(|q| (0.0..=1.0).contains(&q)));
                assert!(f.quant_q >= f.transfer_q, "{f:?}");
            }
        }
        let rows = output::mbr_out::read_mbr_parquet(&mbr).unwrap();
        let accepted = |t: &str| matches!(t, "transfer_ms2" | "transfer_ms1" | "shared_feature");
        assert!(rows.iter().any(|r| accepted(r.evidence_type.as_str())));
        for r in &rows {
            if accepted(r.evidence_type.as_str()) {
                let (t, m) = (r.transfer_q.unwrap(), r.measurement_q.unwrap());
                assert!((0.0..=1.0).contains(&m), "{r:?}");
                assert_eq!(r.quant_q, Some(t.max(m)), "{r:?}");
            } else {
                assert_eq!(r.measurement_q, None, "{r:?}");
            }
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
    fn whole_pattern_divides_quantity_and_error_by_the_share() {
        assert_eq!(whole_pattern((3.0, 0.5), 1.0), (3.0, 0.5));
        assert_eq!(whole_pattern((3.0, 0.5), 0.5), (6.0, 1.0));
        let share = 0.871_108_8_f32;
        let (q, se) = whole_pattern((1e7, 2e5), share);
        assert_eq!((q, se), (1e7 / f64::from(share), 2e5 / f64::from(share)));
        assert!((q / 1e7 - 1.0 / 0.871_108_8).abs() < 1e-7);
    }

    #[test]
    fn finalize_reports_whole_pattern_quantities() {
        // The same rows against two libraries that differ only in the
        // isotope shares: 1 everywhere, and four values by precursor.
        let share = |p: u32| [1.0_f32, 0.948_998_3, 0.871_108_8, 0.5][p as usize % 4];
        let mut lib = library();
        for (p, m) in lib.meta.iter_mut().enumerate() {
            m.isotope_share = share(p as u32);
        }
        let dir = tempfile::tempdir().unwrap();
        let runs = extract_two_identified_runs(2.0);
        let run = |name: &str, lib: &Library| {
            let (mbr, val) = mbr_paths(&dir.path().join(name));
            let s = finalize(&stems(), lib, &runs, 0.01, None, 7, &mbr, &val).unwrap();
            (s.features, table(&mbr), table(&val))
        };
        let (four, four_table, four_validation) = run("four", &library());
        let (whole, whole_table, whole_validation) = run("whole", &lib);
        // The rows, their evidence types and their q-values do not use the
        // share, nor does the validation.
        assert_eq!(four_validation, whole_validation);
        assert_eq!(four.len(), whole.len());
        let (mut scaled, mut normalized) = (0, 0);
        for (f, w) in four.iter().zip(&whole) {
            assert_eq!((f.precursor, f.run), (w.precursor, w.run));
            assert_eq!(
                (f.evidence_type, f.transfer_q, f.quant_q, f.apex_rt_s),
                (w.evidence_type, w.transfer_q, w.quant_q, w.apex_rt_s)
            );
            let s = f64::from(share(f.precursor));
            assert_eq!(w.quantity, f.quantity.map(|q| q / s), "{f:?}");
            assert_eq!(w.quantity_se, f.quantity_se.map(|q| q / s), "{f:?}");
            assert_eq!(
                w.quantity_normalized.is_some(),
                f.quantity_normalized.is_some()
            );
            if let (Some(a), Some(b)) = (f.quantity_normalized, w.quantity_normalized) {
                assert!((b * s / a - 1.0).abs() < 1e-12, "{f:?}");
                normalized += 1;
            }
            if f.precursor % 4 == 0 {
                // Share 1: the same bits.
                assert_eq!(w.quantity_normalized, f.quantity_normalized);
            } else if f.quantity.is_some() {
                scaled += 1;
            }
        }
        assert!(scaled > 300 && normalized > 300, "{scaled} {normalized}");
        // The tables differ in the three quantity columns only.
        assert_eq!(four_table.len(), whole_table.len());
        let quantity_columns = ["quantity", "quantity_se", "quantity_normalized"];
        for (a, b) in four_table.iter().zip(&whole_table) {
            for ((name, a), b) in four_table[0].iter().zip(a).zip(b) {
                if !quantity_columns.contains(&name.as_str()) {
                    assert_eq!(a, b, "{name}");
                }
            }
        }
        let quantity = four_table[0].iter().position(|h| h == "quantity").unwrap();
        assert!(four_table
            .iter()
            .zip(&whole_table)
            .skip(1)
            .any(|(a, b)| a[quantity] != b[quantity]));
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
        // A fit without a positive standard error has no signal-to-noise
        // ratio and gets no features; one without information has ratio 0.
        for se in [0.0, -1.0, f64::NAN] {
            r.evidence.as_mut().unwrap().amplitude_se = se;
            assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0).is_none());
            assert!(identified_features(&r).is_none());
            assert!(amplitude_snr(r.evidence.as_ref().unwrap()).is_nan());
        }
        r.evidence.as_mut().unwrap().amplitude_se = f64::INFINITY;
        let x = transfer_features(&r, &DonorShapes::default(), &[], 0.0).unwrap();
        let snr = TRANSFER_FEATURES
            .iter()
            .position(|n| *n == "log10_amplitude_snr")
            .unwrap();
        assert_eq!(x[snr], 0.0);
        assert!(identified_features(&r).is_some());
        r.evidence.as_mut().unwrap().amplitude_se = 1e4;
        r.evidence.as_mut().unwrap().llr = f64::NAN;
        assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0).is_none());
        r.evidence = None;
        assert!(transfer_features(&r, &DonorShapes::default(), &[], 0.0).is_none());
    }

    #[test]
    fn a_wrong_measurement_is_a_wrong_peak_or_a_wrong_fitted_quantity() {
        let mut g = Rng(3);
        let id = evidence(&mut g, true, 50.0);
        // The same apex and amplitude, fitted 2.5 times as wide: the fitted
        // quantities differ by the widths.
        let mut p = id.clone();
        p.w = 2.5 * id.w;
        p.width_ratio = 2.5 * id.width_ratio;
        p.quantity = 2.5 * id.quantity;
        assert_eq!(
            mismatch(&p, &id),
            Some(Mismatch {
                apex: false,
                quantity: Some(true)
            })
        );
        assert_eq!(wrong_label(&p, &id), Some(true));
        // Within a factor of 2 the quantity agrees; the trace area is not
        // compared.
        p.quantity = 1.9 * id.quantity;
        p.trace_quantity = 3.0 * id.trace_quantity;
        assert!(!mismatch(&p, &id).unwrap().wrong());
        assert_eq!(wrong_label(&p, &id), Some(false));
        // An apex at the identified width is the same peak; beyond it is not.
        p.t0 = id.t0 + id.w;
        assert!(!mismatch(&p, &id).unwrap().apex);
        p.t0 = id.t0 + id.w + 0.1;
        let m = mismatch(&p, &id).unwrap();
        assert!(m.apex && m.quantity == Some(false) && m.wrong());
        // A fit without a quantity: the pair is judged on the apex alone and
        // has no training label.
        let no_quantity = |e: &Evidence| {
            let mut e = e.clone();
            e.amplitude = 0.0;
            e.quantity = 0.0;
            e
        };
        for (a, b) in [(no_quantity(&p), id.clone()), (p.clone(), no_quantity(&id))] {
            let m = mismatch(&a, &b).unwrap();
            assert!(m.apex && m.quantity.is_none() && m.wrong());
            assert_eq!(wrong_label(&a, &b), None);
            let mut a = a;
            a.t0 = b.t0;
            let m = mismatch(&a, &b).unwrap();
            assert!(!m.apex && m.quantity.is_none() && !m.wrong());
        }
        // No judgement without an apex in either fit, or without the
        // identified width.
        let (mut a, mut b) = (p.clone(), id.clone());
        a.t0 = f64::NAN;
        b.t0 = f64::NAN;
        assert_eq!(mismatch(&a, &id), None);
        assert_eq!(mismatch(&p, &b), None);
        let mut b = id.clone();
        b.w = f64::NAN;
        assert_eq!(mismatch(&p, &b), None);
    }

    #[test]
    fn a_pseudo_transfer_copied_into_its_identified_row_has_no_label() {
        let mut g = Rng(21);
        let e = evidence(&mut g, true, 500.0);
        // The identified row took the pseudo-transfer's fit: nothing to
        // compare, so no label and no wrong-peak count.
        assert_eq!(mismatch(&e, &e), None);
        assert_eq!(wrong_label(&e, &e), None);
        // Also when the copied fit carries a NaN feature.
        let mut n = e.clone();
        n.m1_ratio = f64::NAN;
        assert_eq!(mismatch(&n, &n.clone()), None);
        // Any other identified fit is compared as before.
        let mut id = e.clone();
        id.t0 += 3.0 * id.w;
        assert!(mismatch(&e, &id).unwrap().apex);
    }

    #[test]
    fn wrong_features_take_the_fitted_quantity_without_floors() {
        let mut g = Rng(5);
        let mut e = evidence(&mut g, true, 50.0);
        e.quantity = 2e7;
        e.quantity_se = 4e5;
        let x = wrong_features(&e).unwrap();
        assert_eq!(x.len(), 11);
        assert_eq!(x[2], e.width_ratio.ln().abs() as f32);
        assert_eq!(x[8], (1e6f64 / 1e4 + 1.0).log10() as f32);
        assert_eq!(x[9], (4e5f64 / 2e7).ln() as f32);
        assert_eq!(x[10], 2e7f64.log2() as f32);
        let unusable: [fn(&mut Evidence); 9] = [
            |e| e.quantity = 0.0,
            |e| e.quantity = -1.0,
            |e| e.quantity = f64::NAN,
            |e| e.quantity_se = 0.0,
            |e| e.quantity_se = -1.0,
            |e| e.quantity_se = f64::NAN,
            |e| e.amplitude_se = 0.0,
            |e| e.amplitude_se = f64::NAN,
            |e| e.width_ratio = 0.0,
        ];
        for (i, f) in unusable.iter().enumerate() {
            let mut e = e.clone();
            f(&mut e);
            assert!(wrong_features(&e).is_none(), "{i}");
        }
    }

    #[test]
    fn an_accepted_transfer_has_the_larger_of_its_two_q_values() {
        assert_eq!(transfer_quant_q(Some(0.004), 0.002), 0.004);
        assert_eq!(transfer_quant_q(Some(0.004), 0.02), 0.02);
        assert_eq!(transfer_quant_q(None, 0.3), 0.3);
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
            reporters: None,
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
