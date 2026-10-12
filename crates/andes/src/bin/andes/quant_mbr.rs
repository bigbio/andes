//! Match-between-runs stages shared by the stage subcommands and the
//! all-in-one run.
//!
//! Extract: for one run, every library precursor is measured on the run's
//! cached MS1 by the template fit of [`quant::evidence`], whether the run
//! identified it or not, and the rows go to `<stem>.extract.parquet`.
//!
//! Per library entry and run:
//!
//! * identified in the run: an `Identified` target fit with the prior centred
//!   on the best PSM's RT and `s = max(w_prior, 6 s)`, plus a `PseudoTransfer`
//!   target fit with the alignment prior computed without this run's own
//!   observation;
//! * not identified: a `Transfer` target and its decoy probes, with the
//!   alignment prior;
//! * an unaligned run gets no rows; pairs without a prediction, and identified
//!   pairs whose best PSM RT is not finite, are skipped and counted.
//!
//! An entry is fitted on the MS1 series of its FAIMS compensation voltage;
//! when the run has no series at that voltage, or the entry has none, on the
//! whole run. Only one copy of the run's MS1 is held at a time: the whole run
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

use std::path::Path;

use quant::decoys::{probes, Hypothesis, Probe};
use quant::evidence::{fit_on_grid, fit_with, EvidenceParams, FitScratch};
use quant::library::{Library, LibraryEntry};
use quant::ms1_index::Ms1RunIndex;
use rayon::prelude::*;

use crate::quant_run::ms1_cache::Ms1Cache;

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
    /// The run is unaligned: no rows were produced for it.
    pub unaligned: bool,
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

/// Per entry, the position of its FAIMS voltage among `voltages` (the
/// run's split order), or `None` for the whole run. A run with at most one
/// voltage is fitted whole.
fn entry_groups(voltages: &[f32], entries: &[LibraryEntry]) -> Vec<Option<usize>> {
    entries
        .iter()
        .map(|e| {
            if voltages.len() <= 1 {
                return None;
            }
            e.faims_cv
                .and_then(|cv| voltages.iter().position(|&v| v == cv))
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

/// The usable width-free fits `(t0, w)` of the identified pairs among the
/// entries `ids`, each fitted on `ms1_of(id)`.
fn width_fits<'a, F>(
    ids: &[usize],
    run: u16,
    lib: &Library,
    ep: &EvidenceParams,
    ms1_of: F,
) -> Vec<(f64, f64)>
where
    F: Fn(usize) -> &'a Ms1RunIndex + Sync,
{
    ids.par_iter()
        .map_init(FitScratch::default, |scratch, &i| {
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
            let ev = fit_on_grid(
                ms1_of(i),
                &probe,
                meta.charge.max(1),
                &meta.envelope,
                ep,
                &WIDTH_PRIOR_GRID,
                scratch,
            )?;
            (ev.amplitude > 0.0 && ev.llr >= WIDTH_FIT_MIN_LLR).then_some((ev.t0, ev.w))
        })
        .flatten()
        .collect()
}

type EntryRows = (Vec<ExtractRow>, ExtractCounts);

/// Fills `slots[id]` with the rows of each entry in `ids`, fitted on
/// `ms1_of(id)`.
#[allow(clippy::too_many_arguments)]
fn rows_into<'a, F>(
    slots: &mut [Option<EntryRows>],
    ids: &[usize],
    run: u16,
    lib: &Library,
    wp: &WidthPrior,
    ep: &EvidenceParams,
    ms1_of: F,
) where
    F: Fn(usize) -> &'a Ms1RunIndex + Sync,
{
    let done: Vec<EntryRows> = ids
        .par_iter()
        .map_init(FitScratch::default, |scratch, &i| {
            entry_rows(&lib.entries[i], run, lib, ms1_of(i), wp, ep, scratch)
        })
        .collect();
    for (&i, r) in ids.iter().zip(done) {
        slots[i] = Some(r);
    }
}

/// The rows of one library entry in `run`.
fn entry_rows(
    e: &LibraryEntry,
    run: u16,
    lib: &Library,
    ms1: &Ms1RunIndex,
    wp: &WidthPrior,
    ep: &EvidenceParams,
    scratch: &mut FitScratch,
) -> (Vec<ExtractRow>, ExtractCounts) {
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
    let mut row = |kind: Hypothesis, role: Role, probe: &Probe, n_donor_runs: u16| ExtractRow {
        precursor: e.precursor,
        run,
        kind,
        role,
        evidence: fit_with(ms1, probe, charge, &meta.envelope, ep, scratch),
        n_donor_runs,
        library_q: e.q,
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
/// An unaligned run yields no rows and `counts.unaligned` (its cache is not
/// read). Errors: `run` outside the library's runs, an entry whose precursor
/// has no metadata, or an unreadable cache (the message names the file).
#[allow(dead_code)]
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
    if map.is_none() {
        return Ok(RunExtract {
            rows: Vec::new(),
            counts: ExtractCounts {
                unaligned: true,
                ..ExtractCounts::default()
            },
        });
    }
    let load = || {
        #[cfg(test)]
        tracking::loaded();
        Ms1Cache::load_from(cache).map_err(|e| e.to_string())
    };
    let n = lib.entries.len();
    let ms1 = load()?;
    let held = Held::new();
    let groups = entry_groups(&ms1.faims_voltages(), &lib.entries);
    let whole_ids: Vec<usize> = (0..n).filter(|&i| groups[i].is_none()).collect();
    let cv_ids: Vec<usize> = (0..n).filter(|&i| groups[i].is_some()).collect();
    let mut slots: Vec<Option<EntryRows>> = (0..n).map(|_| None).collect();

    let mut fits = width_fits(&whole_ids, run, lib, ep, |_| &ms1);
    let width_fit_count;
    if cv_ids.is_empty() {
        width_fit_count = fits.len();
        let wp = WidthPrior::from_fits(fits);
        rows_into(&mut slots, &whole_ids, run, lib, &wp, ep, |_| &ms1);
        drop(ms1);
        drop(held);
        return Ok(assemble(slots, width_fit_count, wp.is_default()));
    }

    drop(held);
    let split = ms1.split_by_faims();
    let held = Held::new();
    let series_of = |i: usize| &split[groups[i].expect("voltage entry")].1;
    fits.extend(width_fits(&cv_ids, run, lib, ep, series_of));
    width_fit_count = fits.len();
    let wp = WidthPrior::from_fits(fits);
    rows_into(&mut slots, &cv_ids, run, lib, &wp, ep, series_of);
    drop(split);
    drop(held);

    if !whole_ids.is_empty() {
        let ms1 = load()?;
        let _held = Held::new();
        rows_into(&mut slots, &whole_ids, run, lib, &wp, ep, |_| &ms1);
    }
    Ok(assemble(slots, width_fit_count, wp.is_default()))
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
        fixture_with(faims, &|_| false)
    }

    /// [`fixture`] where the precursors in `silent_in_run1` give no MS1
    /// signal in run 1 (they stay identified there).
    fn fixture_with(faims: bool, silent_in_run1: &dyn Fn(u32) -> bool) -> Fixture {
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
        let lib = build_library(&psms, meta, N_RUNS, 1.0, &ap);
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
    fn unaligned_run_is_skipped_with_count() {
        let fx = fixture(false);
        assert!(fx.lib.alignment.maps[2].is_none());
        let out = extract_run(2, &fx.caches[2], &fx.lib, &params()).unwrap();
        assert!(out.rows.is_empty());
        assert!(out.counts.unaligned);
        assert!(extract_run(N_RUNS as u16, &fx.caches[2], &fx.lib, &params()).is_err());
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
        let fx = fixture_with(false, &silent);
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
