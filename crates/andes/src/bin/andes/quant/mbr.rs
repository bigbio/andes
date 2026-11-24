use std::collections::HashMap;
use std::path::{Path, PathBuf};

use model::stats::{median_mad_mut, median_mut, Even};
use output::ms1::read_ms1_table;
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

use crate::transfers::{IDENTIFIED_FEATURES, TRANSFER_FEATURES};

pub(crate) use output::extract::{ExtractRole as Role, ExtractRow};

const SEED_WIDTH_S: f64 = 6.0;
const WIDTH_PRIOR_GRID: [f64; 5] = [0.5, 0.7, 1.0, 1.4, 2.0];
const SCALE_PAIRS: usize = 500;
const SCALE_MAX_ROUNDS: usize = 6;
const SCALE_CONVERGED: f64 = 1.1;
const WIDTH_PRIOR_BINS: usize = 20;
const WIDTH_FIT_MIN_LLR: f64 = 25.0;

#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct ExtractCounts {
    pub unaligned: bool,
    pub foreign_cv: usize,
    pub nonfinite_best_rt: usize,
    pub no_pseudo_prediction: usize,
    pub no_prediction: usize,
    pub pseudo_other_component: usize,
    pub other_component: usize,
    pub not_selected: usize,
    pub consensus_peak: usize,
    pub consensus_peak_decoys: usize,
    pub width_fits: usize,
    pub width_prior_fallback: bool,
    pub scales: RunScales,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct RunScales {
    pub spacing_s: Option<f64>,
    pub width_s: f64,
    pub width_fits: usize,
    pub rounds: usize,
    pub fallback: bool,
    pub mobility: Option<RunMobility>,
    pub claim_windows: Option<claims::Windows>,
}

impl RunScales {
    fn floor_s(&self) -> f64 {
        self.spacing_s.unwrap_or(0.0)
    }

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
        self.consensus_peak_decoys += o.consensus_peak_decoys;
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TransferSelection {
    pub min_donor_fraction: f64,
}

pub(crate) const DEFAULT_MIN_DONOR_FRACTION: f64 = 0.0;

impl TransferSelection {
    pub(crate) const OFF: TransferSelection = TransferSelection {
        min_donor_fraction: 0.0,
    };

    pub(crate) fn is_off(&self) -> bool {
        self.min_donor_fraction <= 0.0
    }

    fn min_donors(&self, aligned_runs: usize) -> u16 {
        let others = aligned_runs.saturating_sub(1) as f64;
        (self.min_donor_fraction * others)
            .ceil()
            .clamp(0.0, u16::MAX as f64) as u16
    }

    pub(crate) fn describe(&self, aligned_runs: usize) -> String {
        format!(
            "min_donor_fraction={} ({} of {} other aligned runs)",
            self.min_donor_fraction,
            self.min_donors(aligned_runs),
            aligned_runs.saturating_sub(1)
        )
    }
}

impl Default for TransferSelection {
    fn default() -> Self {
        TransferSelection {
            min_donor_fraction: DEFAULT_MIN_DONOR_FRACTION,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct PairGate {
    min_donors: u16,
}

impl PairGate {
    fn new(sel: &TransferSelection, lib: &Library) -> Self {
        let aligned = lib.alignment.maps.iter().filter(|m| m.is_some()).count();
        PairGate {
            min_donors: sel.min_donors(aligned),
        }
    }

    fn keeps(&self, _e: &LibraryEntry, other_donors: u16) -> bool {
        other_donors >= self.min_donors
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RunExtract {
    pub rows: Vec<ExtractRow>,
    pub counts: ExtractCounts,
}

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

#[cfg(test)]
#[path = "mbr/tests/tracking.rs"]
mod tracking;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Group {
    Whole,
    Series(usize),
    Foreign,
}

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

#[derive(Debug, Clone, PartialEq)]
struct WidthPrior {
    knots: Vec<(f64, f64)>,
    default_s: f64,
    floor_s: f64,
}

impl WidthPrior {
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
    fn identified_s(&self, t: f64) -> f64 {
        self.at(t).max(self.floor_s)
    }
}

pub(crate) fn median(mut v: Vec<f64>) -> f64 {
    median_mut(&mut v, Even::Mean).unwrap_or(f64::NAN)
}

fn identified_rt(entry: &LibraryEntry, run: u16) -> Option<f64> {
    entry
        .best_rt
        .iter()
        .find(|(r, _)| *r == run)
        .map(|&(_, rt)| rt)
}

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

fn usable_width_fit(ev: &Evidence) -> bool {
    ev.amplitude > 0.0 && ev.llr >= WIDTH_FIT_MIN_LLR
}

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
    pool: &quant::decoys::RtPool,
) where
    F: Fn(usize) -> usize + Sync,
{
    let planned: Vec<PlannedRows> = ids
        .par_iter()
        .map(|&i| {
            let e = &lib.entries[i];
            entry_rows(e, run, aligned, lib, gate, series[series_of(i)], wp, pool)
        })
        .collect();
    let mut jobs = Vec::new();
    for (&i, (rows, counterparts, _)) in ids.iter().zip(&planned) {
        let meta = &lib.meta[lib.entries[i].precursor as usize];
        for probe in rows.iter().map(|(_, p)| p).chain(counterparts.iter().map(|(_, p)| p)) {
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
    for (&i, (rows, counterparts, mut counts)) in ids.iter().zip(planned) {
        let mut rows: Vec<ExtractRow> = rows
            .into_iter()
            .map(|(mut row, _)| {
                row.evidence = evidence.next().flatten();
                row
            })
            .collect();
        let counterparts: Vec<(usize, Option<Evidence>)> = counterparts
            .into_iter()
            .map(|(at, _)| (at, evidence.next().flatten()))
            .collect();
        let (target, decoys) = prefer_consensus_peak(&mut rows, &counterparts);
        counts.consensus_peak += usize::from(target);
        counts.consensus_peak_decoys += decoys;
        slots[i] = Some((rows, counts));
    }
}

pub(crate) const CONSENSUS_PEAK_FACTOR: f64 = 2.0;

fn takes_other_peak(own: &Evidence, other: &Evidence) -> bool {
    (other.t0 - own.t0).abs() > own.w + other.w
        && own.amplitude > 0.0
        && other.amplitude > CONSENSUS_PEAK_FACTOR * own.amplitude
}

fn prefer_consensus_peak(
    rows: &mut [ExtractRow],
    counterparts: &[(usize, Option<Evidence>)],
) -> (bool, usize) {
    let mut decoys = 0;
    for (i, other) in counterparts {
        if let (Some(own), Some(other)) = (&rows[*i].evidence, other) {
            if takes_other_peak(own, other) {
                rows[*i].evidence = Some(other.clone());
                decoys += 1;
            }
        }
    }
    let at = |role: Role| {
        rows.iter()
            .position(|r| r.role == role && r.kind == Hypothesis::Target)
    };
    let (Some(id), Some(pseudo)) = (at(Role::Identified), at(Role::PseudoTransfer)) else {
        return (false, decoys);
    };
    let (Some(a), Some(b)) = (&rows[id].evidence, &rows[pseudo].evidence) else {
        return (false, decoys);
    };
    if !takes_other_peak(a, b) {
        return (false, decoys);
    }
    rows[id].evidence = rows[pseudo].evidence.clone();
    (true, decoys)
}

type PlannedRows = (Vec<(ExtractRow, Probe)>, Vec<(usize, Probe)>, ExtractCounts);

#[allow(clippy::too_many_arguments)]
fn entry_rows(
    e: &LibraryEntry,
    run: u16,
    run_aligned: bool,
    lib: &Library,
    gate: &PairGate,
    ms1: &Ms1RunIndex,
    wp: &WidthPrior,
    pool: &quant::decoys::RtPool,
) -> PlannedRows {
    let mut counts = ExtractCounts::default();
    let mut rows = Vec::new();
    let mut counterparts = Vec::new();
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
                return (rows, counterparts, counts);
            }
            match lib.alignment.predict_excluding(e.precursor, run) {
                Some((t_hat, s)) => {
                    let (lo, hi) = run_rt();
                    for (k, (r, p)) in rows.iter().enumerate() {
                        if r.kind == Hypothesis::Target {
                            continue;
                        }
                        let t = p.t_hat + (t_hat - rt);
                        if t >= lo && t <= hi {
                            counterparts.push((
                                k,
                                Probe {
                                    t_hat: t,
                                    s,
                                    w_prior: wp.at(t),
                                    ..*p
                                },
                            ));
                        }
                    }
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
                let w = wp.at(t_hat);
                for probe in probes(e.precursor, run, meta.mz, charge, t_hat, s, w, run_rt())
                    .into_iter()
                    .filter(|p| p.kind != Hypothesis::RtShiftDecoy)
                {
                    rows.push(row(probe.kind, Role::Transfer, &probe, other_donors));
                }
                let mass = (meta.mz - model::mass::PROTON) * charge as f64;
                for probe in quant::decoys::random_rt_probes(
                    e.precursor,
                    run,
                    meta.mz,
                    mass,
                    t_hat,
                    s,
                    w,
                    run_rt(),
                    pool,
                    quant::decoys::RANDOM_RT_DECOYS,
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
    (rows, counterparts, counts)
}

const MOBILITY_PAIRS: usize = 2000;

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

fn rt_pool(lib: &Library, run: u16) -> quant::decoys::RtPool {
    quant::decoys::RtPool::new(
        lib.entries
            .iter()
            .filter_map(|e| {
                let meta = lib.meta.get(e.precursor as usize)?;
                let (t, _) = lib.alignment.predict(e.precursor, run)?;
                Some(((meta.mz - model::mass::PROTON) * meta.charge.max(1) as f64, t))
            })
            .collect(),
    )
}

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
    let pool = rt_pool(lib, run);
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
            &pool,
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
        &mut slots, &cv_ids, run, aligned, lib, &gate, &wp, ep, &series, series_of, &gates, &pool,
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
            &pool,
        );
    }
    Ok(finish(slots, width_fit_count, &wp, scales))
}

#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct Ms2Counts {
    pub rows: [[usize; 3]; 2],
    pub with_scan: [[usize; 3]; 2],
    pub supported: [[usize; 3]; 2],
    pub no_spectrum: usize,
}

fn probe_mz(meta: &quant::library::PrecursorMeta, kind: Hypothesis) -> f64 {
    match kind {
        Hypothesis::MassShiftDecoy(n) => {
            meta.mz + n as f64 * PEPTIDE_MASS_SLOPE / meta.charge.max(1) as f64
        }
        Hypothesis::Target | Hypothesis::RtShiftDecoy | Hypothesis::RandomRtDecoy(_) => meta.mz,
    }
}

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

#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct SharedCounts {
    pub claimed: usize,
    pub with_scan: usize,
    pub scored: usize,
}

fn peptidoform_of(key: &str) -> &str {
    key.rsplit_once('/').map_or(key, |(p, _)| p)
}

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

pub(crate) fn local_slope(maps: &[Option<RunMap>], run: u16, rt: f64) -> Option<f64> {
    let map = maps.get(run as usize)?.as_ref()?;
    let slope = map.slope_at(map.invert(rt));
    (slope.is_finite() && slope > 0.0).then_some(slope)
}

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

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct WidthCell {
    pub precursor: u32,
    pub run: u16,
    pub log2_quantity: f64,
    pub deviation: f64,
}

const WIDTH_WEIGHT_STEPS: u32 = 10;
const WIDTH_WEIGHT_MIN_PRECURSOR_CELLS: usize = 4;
const WIDTH_WEIGHT_MIN_RUN_CELLS: usize = 11;
const WIDTH_WEIGHT_PASSES: usize = 2;
const WIDTH_WEIGHT_TOLERANCE: f64 = 1e-6;

pub(crate) fn run_width_weights(cells: &[WidthCell], runs: usize) -> Vec<f64> {
    let mut weights = vec![0.0; runs];
    let mut order: Vec<usize> = (0..cells.len())
        .filter(|&i| (cells[i].run as usize) < runs)
        .collect();
    order.sort_by_key(|&i| (cells[i].precursor, cells[i].run));
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
                let spread = median_mad_mut(&mut res, Even::Mean).map_or(f64::NAN, |(_, mad)| mad);
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

pub(crate) fn row_width_weight(r: &ExtractRow, weights: &[f64]) -> f64 {
    if r.role == Role::Identified {
        weights.get(r.run as usize).copied().unwrap_or(0.0)
    } else {
        0.0
    }
}

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

pub(crate) fn whole_pattern((quantity, se): (f64, f64), share: f32) -> (f64, f64) {
    let share = f64::from(share);
    (quantity / share, se / share)
}

pub(crate) fn amplitude_snr(e: &Evidence) -> f64 {
    if e.amplitude_se > 0.0 {
        e.amplitude / e.amplitude_se
    } else {
        f64::NAN
    }
}

pub(crate) const NO_SHAPE: [f64; 3] = [0.0, 0.0, -1.0];

#[cfg(test)]
#[derive(Debug, Clone, Default)]
pub(crate) struct DonorShapes {
    fits: HashMap<u32, Vec<DonorFit>>,
}

pub(crate) type DonorFit = (u16, f64, f64, [f64; 4]);

pub(crate) fn iso_profile(h: &[f64; 4]) -> Option<[f64; 4]> {
    let c = h.map(|v| if v.is_finite() { v.max(0.0) } else { 0.0 });
    let sum: f64 = c.iter().sum();
    (sum > 0.0).then(|| c.map(|v| v / sum))
}

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

pub(crate) trait ShapeReference {
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

pub(crate) fn transfer_features(
    row: &ExtractRow,
    shapes: &impl ShapeReference,
    maps: &[Option<RunMap>],
    neighbour_rate: f64,
    donor_pep: f64,
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
        donor_pep.max(1e-12).log10(),
    ];
    debug_assert_eq!(x.len(), TRANSFER_FEATURES.len());
    let x: Vec<f32> = x.iter().map(|&v| v as f32).collect();
    x.iter().all(|v| v.is_finite()).then_some(x)
}

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

pub(crate) fn family(kind: Hypothesis) -> u8 {
    match kind {
        Hypothesis::Target => 0,
        Hypothesis::MassShiftDecoy(_) => 1,
        Hypothesis::RtShiftDecoy | Hypothesis::RandomRtDecoy(_) => 2,
    }
}

pub(crate) fn mbr_paths(base: &Path) -> (PathBuf, PathBuf) {
    let with = |suffix: &str| {
        let mut s = base.as_os_str().to_owned();
        s.push(suffix);
        PathBuf::from(s)
    };
    (with(".mbr.parquet"), with(".mbr_validation.parquet"))
}

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
#[path = "mbr/tests/mbr.rs"]
mod tests;

#[cfg(test)]
#[path = "mbr/tests/finalize.rs"]
mod finalize_tests;
