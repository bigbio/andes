//! `andes quant finalize`: the transfer FDR, the identified rows' quality
//! q-values, the run normalisation and the MBR outputs over the extract
//! tables of every run (see [`crate::quant_mbr`] for what each computes), as
//! a map-reduce over runs.
//!
//! Every quantity a row needs from other runs is a per-precursor summary or
//! a global model, so each run's table is read on its own, in record batches
//! cut at precursor boundaries, and the runs are processed in parallel. Each
//! pass reads every table once:
//!
//! 1. donor summaries: per precursor, the consensus width `w_c` and the
//!    sorted values around the middle of its identified fits (width on the
//!    consensus scale, tailing, isotope heights), which give the median of
//!    every fit and of every fit but one run's: the shape features leave the
//!    row's own run out, and that run's fit is in the same table. With more
//!    donor fits than [`FinalizeLimits::max_donor_fits`] the precursors are
//!    summarised in groups, one pass each;
//! 2. the training rows of both classifiers, a [`TrainingSample`] keyed by
//!    run and row, so that the kept rows enter training in table order (one
//!    sample per member of the bagged transfer classifier, see
//!    [`FinalizeLimits::transfer_models`]);
//! 3. every row's held-out score, into [`FoldScores`] merged run by run;
//!    above [`FinalizeLimits::max_held_scores`] scores in memory they are
//!    sorted and written to temporary tables, and read back in order when
//!    the q-value step functions are built. Each run's transfer scores are
//!    also written, in row order, to a temporary table;
//! 4. the reference cells of the run normalisation;
//! 5. per run, the output rows, sorted by precursor, into a temporary table,
//!    and the validation counts; the transfer rows take their scores from
//!    pass 3's table instead of being scored again;
//!
//! then the per-run tables are merged by (precursor, run) into
//! `<base>.mbr.parquet` and handed to the caller.
//!
//! Every reduction is exact and does not depend on the order in which the
//! runs finish, so the outputs are those of one pass over every row in
//! memory, byte for byte. Memory follows the library (the summaries), the
//! limits, the q-value step functions and one run's output rows per worker
//! thread, not the number of runs.

use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use arrow::array::{
    Array, ArrayRef, Float32Array, Float64Array, Float64Builder, UInt16Array, UInt16Builder,
    UInt32Array, UInt32Builder, UInt8Array, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use output::extract_table::ExtractTableReader;
use output::mbr_out::{
    write_mbr_validation, Accepted, EvidenceType, MbrClaim, MbrParquetWriter, MbrRow,
    MbrValidation, VALIDATION_Q,
};
use parquet::arrow::arrow_reader::{ParquetRecordBatchReader, ParquetRecordBatchReaderBuilder};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;
use quant::align::RunMap;
use quant::claims::Claim;
use quant::decoys::Hypothesis;
use quant::evidence::Evidence;
use quant::library::{Library, LibraryEntry};
use quant::normalize::{fit as fit_normalization, NormCell, NormParams};
use quant::shared;
use rayon::prelude::*;

use crate::quant_mbr::{
    consensus_scale_width, cosine, donor_fit, family, identified_features, iso_profile,
    local_slope, median, reported_quantity, transfer_features, DonorFit, ExtractRow, Role,
    ShapeReference, NO_SHAPE, TRANSFER_DONOR_Q,
};
use crate::transfer_fdr::{
    member_seed, pack_family, row_unit, sample_hash, stratum_of, train_bagged_folds, train_folds,
    transfer_train_params, Estimate, FdrCounts, FdrRow, FoldModels, FoldScores, QValues,
    TrainingSample, IDENTITY, MAX_STRATA, TARGET, TRANSFER_FEATURES,
};

/// What one finalize did, for the caller to log.
#[derive(Debug, Clone, Default)]
pub(crate) struct FinalizeSummary {
    pub validation: MbrValidation,
    /// Rows of `<base>.mbr.parquet`.
    pub rows_written: usize,
    /// Transfer targets and decoys without evidence or with a non-finite
    /// feature (no q-value).
    pub no_features: usize,
    pub fdr: FdrCounts,
    /// Identified targets and decoys without evidence or with a non-finite
    /// feature (no quality q-value).
    pub identified_no_features: usize,
    pub identified_fdr: FdrCounts,
    /// The rows of `<base>.mbr.parquet`, in file order.
    #[cfg(test)]
    pub features: Vec<MbrFeature>,
}

/// One row of `<base>.mbr.parquet` by precursor id and run index.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MbrFeature {
    pub precursor: u32,
    pub run: u16,
    pub evidence_type: EvidenceType,
    pub quantity: Option<f64>,
    /// Least-squares standard error of the fitted amplitude at the fixed
    /// peak shape, times the area factor: the fit error, not the cell's
    /// run-to-run uncertainty.
    pub quantity_se: Option<f64>,
    pub quantity_normalized: Option<f64>,
    pub transfer_q: Option<f64>,
    /// For an accepted transfer, the running mean of the wrong-peak
    /// probability of a model trained on pseudo-transfers (identified cells
    /// re-measured blind): it sees a wrong peak, not an absent peptide or a
    /// wrong donor identity.
    pub quant_q: Option<f64>,
    pub apex_rt_s: Option<f64>,
    pub library_q: f64,
    /// Intensity-weighted isotope-0 mass error of the fit, ppm.
    pub mass_err_ppm: Option<f64>,
    /// The same-run claims on a transfer row (default for other rows).
    pub claim: Claim,
}

/// Limits of finalize: memory follows the library and these limits rather
/// than the number of runs or rows. None of them changes the outputs except
/// the two sample caps and the transfer classifier's members.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FinalizeLimits {
    /// Training rows per classifier (transfers, identified rows), and per
    /// member of the transfer classifier. Above it the classifier is trained
    /// on the rows of a hashed sample of (precursor, run) pairs (every row
    /// still scored by its fold's model).
    pub max_train_rows: usize,
    /// Members of the bagged transfer classifier: each fold's score is the
    /// mean margin of this many models, each trained on its own sample of
    /// the training rows (the same rows when they are all kept) with its
    /// own seed. Training and scoring time grow with it; one member is the
    /// single model.
    pub transfer_models: usize,
    /// Reference cells of the run normalisation. Above it the fit uses a
    /// hashed sample of whole precursors (coarser RT-local offsets).
    pub max_norm_cells: usize,
    /// Donor fits (identified targets) summarised in one pass; above it the
    /// precursors are summarised in groups, one pass over the tables each.
    pub max_donor_fits: usize,
    /// Held-out scores kept in memory per classifier; above it they are
    /// sorted and written to temporary tables.
    pub max_held_scores: usize,
    /// Directory of the temporary tables (removed at the end); `None` for
    /// the directory of `<base>.mbr.parquet`.
    pub temp_dir: Option<PathBuf>,
}

/// Default of [`FinalizeLimits::max_train_rows`].
pub(crate) const DEFAULT_MAX_TRAIN_ROWS: usize = 3_000_000;
/// Default share of the precursors that report their trace areas (see
/// [`crate::quant_mbr`], "Trace area"). On six LFQ benchmarks 0.1 lowered
/// the gross errors of the brightest tenth of the identified cells by 17 to
/// 43 % and the replicate SD on every set, with differential-expression
/// calls unchanged; 1 made them worse on three of four sets.
pub(crate) const DEFAULT_TRACE_FRACTION: f64 = 0.1;
/// Default of [`FinalizeLimits::transfer_models`].
pub(crate) const DEFAULT_TRANSFER_MODELS: usize = 8;
/// Default of [`FinalizeLimits::max_norm_cells`].
pub(crate) const DEFAULT_MAX_NORM_CELLS: usize = 20_000_000;
/// Default of [`FinalizeLimits::max_donor_fits`].
pub(crate) const DEFAULT_MAX_DONOR_FITS: usize = 12_000_000;
/// Default of [`FinalizeLimits::max_held_scores`].
pub(crate) const DEFAULT_MAX_HELD_SCORES: usize = 64_000_000;

impl Default for FinalizeLimits {
    fn default() -> Self {
        Self {
            max_train_rows: DEFAULT_MAX_TRAIN_ROWS,
            transfer_models: DEFAULT_TRANSFER_MODELS,
            max_norm_cells: DEFAULT_MAX_NORM_CELLS,
            max_donor_fits: DEFAULT_MAX_DONOR_FITS,
            max_held_scores: DEFAULT_MAX_HELD_SCORES,
            temp_dir: None,
        }
    }
}

/// Rows whose features are built and scored together.
pub(crate) const SCORE_CHUNK: usize = 65_536;
/// Rows per record batch read from an extract table.
const READ_BATCH_ROWS: usize = 65_536;
/// Rows per row group of a temporary score table.
const SCORE_GROUP_ROWS: usize = 65_536;
/// Rows per row group of a temporary output-row table, per batch read back
/// and per merged batch: small, as the merge reads every run's table at once.
const OUT_GROUP_ROWS: usize = 4_096;
/// Temporary tables merged at once; more are merged in groups first.
const MAX_OPEN_PARTS: usize = 512;

/// The extract rows of every run, as tables ("chunks") read in batches. A
/// precursor's rows of one run are contiguous and in one chunk, as `andes
/// quant extract` writes them (one table per run, rows by library entry).
pub(crate) trait ExtractSource: Sync {
    /// Tables.
    fn chunks(&self) -> usize;
    /// The table's name, for messages.
    fn name(&self, chunk: usize) -> String;
    /// Calls `f` with the rows of table `chunk` in order, a batch at a time.
    fn for_each_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String>;

    /// [`ExtractSource::for_each_batch`] over the identified targets only,
    /// of which only the keys and the peak shape (see [`donor_fit`] and
    /// [`consensus_scale_width`]) are read.
    fn for_each_fit_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String> {
        self.for_each_batch(chunk, &mut |rows| {
            f(rows
                .into_iter()
                .filter(|r| (r.role, r.kind) == (Role::Identified, Hypothesis::Target))
                .collect())
        })
    }
}

/// Extract rows in memory, one vector per run, handed out in batches of
/// `batch_rows` rows.
#[cfg(test)]
pub(crate) struct ExtractRuns<'a> {
    pub runs: &'a [Vec<ExtractRow>],
    pub batch_rows: usize,
}

#[cfg(test)]
impl ExtractSource for ExtractRuns<'_> {
    fn chunks(&self) -> usize {
        self.runs.len()
    }

    fn name(&self, chunk: usize) -> String {
        format!("extract rows {chunk}")
    }

    fn for_each_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String> {
        self.runs[chunk]
            .chunks(self.batch_rows.max(1))
            .try_for_each(|b| f(b.to_vec()))
    }
}

/// The `<stem>.extract.parquet` files of every run.
pub(crate) struct ExtractFiles(pub Vec<PathBuf>);

impl ExtractSource for ExtractFiles {
    fn chunks(&self) -> usize {
        self.0.len()
    }

    fn name(&self, chunk: usize) -> String {
        self.0[chunk].display().to_string()
    }

    fn for_each_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String> {
        for b in ExtractTableReader::open(&self.0[chunk], READ_BATCH_ROWS)? {
            f(b?)?;
        }
        Ok(())
    }

    fn for_each_fit_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String> {
        for b in ExtractTableReader::open_identified_fits(&self.0[chunk], READ_BATCH_ROWS)? {
            f(b?)?;
        }
        Ok(())
    }
}

/// The runs and precursors extract rows may index.
#[derive(Debug, Clone, Copy)]
struct Bounds {
    runs: usize,
    precursors: usize,
}

/// Calls `f(state, chunk, first, block)` on the rows of every table of
/// `source` in blocks of whole precursors (`first` is the index of the
/// block's first row in its table), tables in parallel, each table with its
/// own `state` from `init`, handed to `done` after its last block. With
/// `fits_only`, the identified targets only
/// ([`ExtractSource::for_each_fit_batch`]).
///
/// Errors: a row outside `bounds`, a precursor whose rows are not contiguous
/// in its table, and any error of `source`, `f` or `done`.
fn for_each_block<S>(
    source: &dyn ExtractSource,
    fits_only: bool,
    bounds: Bounds,
    init: impl Fn(usize) -> S + Sync,
    f: impl Fn(&mut S, usize, u64, &[ExtractRow]) -> Result<(), String> + Sync,
    done: impl Fn(usize, S) -> Result<(), String> + Sync,
) -> Result<(), String> {
    (0..source.chunks()).into_par_iter().try_for_each(|chunk| {
        let mut state = init(chunk);
        let mut closed: HashSet<u32> = HashSet::new();
        let mut emit = |state: &mut S, rows: &[ExtractRow], first: u64| {
            if let Some(r) = rows
                .iter()
                .find(|r| r.run as usize >= bounds.runs || r.precursor as usize >= bounds.precursors)
            {
                return Err(format!(
                    "extract row of precursor {} in run {} is outside the {} precursors and {} runs",
                    r.precursor, r.run, bounds.precursors, bounds.runs
                ));
            }
            for g in rows.chunk_by(|a, b| a.precursor == b.precursor) {
                if !closed.insert(g[0].precursor) {
                    return Err(format!(
                        "{}: the rows of precursor {} are not contiguous",
                        source.name(chunk),
                        g[0].precursor
                    ));
                }
            }
            f(state, chunk, first, rows)
        };
        let mut buf: Vec<ExtractRow> = Vec::new();
        let mut first = 0u64;
        let mut each = |batch: Vec<ExtractRow>| {
            buf.extend(batch);
            let Some(last) = buf.last().map(|r| r.precursor) else {
                return Ok(());
            };
            // The last precursor may go on in the next batch.
            let cut = buf.len() - buf.iter().rev().take_while(|r| r.precursor == last).count();
            if cut > 0 {
                emit(&mut state, &buf[..cut], first)?;
                first += cut as u64;
                buf.drain(..cut);
            }
            Ok(())
        };
        if fits_only {
            source.for_each_fit_batch(chunk, &mut each)?;
        } else {
            source.for_each_batch(chunk, &mut each)?;
        }
        if !buf.is_empty() {
            emit(&mut state, &buf, first)?;
        }
        done(chunk, state)
    })
}

/// Per (precursor, run) of a block: the identified target's evidence and
/// donor fit.
pub(crate) struct BlockContext<'a> {
    identified: HashMap<(u32, u16), &'a Evidence>,
    own: HashMap<(u32, u16), DonorFit>,
}

impl<'a> BlockContext<'a> {
    /// Errors on two identified targets with evidence of one (precursor,
    /// run).
    pub(crate) fn new(rows: &'a [ExtractRow], maps: &[Option<RunMap>]) -> Result<Self, String> {
        let mut identified = HashMap::new();
        let mut own = HashMap::new();
        for r in rows {
            if let (Role::Identified, Hypothesis::Target, Some(e)) = (r.role, r.kind, &r.evidence) {
                if identified.insert((r.precursor, r.run), e).is_some() {
                    return Err(format!(
                        "two identified rows of precursor {} in run {}",
                        r.precursor, r.run
                    ));
                }
                if let Some(f) = donor_fit(r, maps) {
                    own.insert((r.precursor, r.run), f);
                }
            }
        }
        Ok(Self { identified, own })
    }
}

/// Sorted values around the middle of one value of a precursor's donor fits:
/// enough for the median of them all and of them all but one.
#[derive(Debug, Clone, Copy)]
struct MedianWindow {
    n: u32,
    /// Index of `v[0]` among the sorted values.
    first: u32,
    v: [f64; 3],
}

impl MedianWindow {
    /// From every value, sorted ascending by [`f64::total_cmp`].
    fn new(sorted: &[f64]) -> Self {
        let n = sorted.len();
        let first = (n / 2).saturating_sub(1);
        let mut v = [f64::NAN; 3];
        for (slot, x) in v.iter_mut().zip(sorted.iter().skip(first)) {
            *slot = *x;
        }
        Self {
            n: n as u32,
            first: first as u32,
            v,
        }
    }

    fn at(&self, i: usize) -> f64 {
        self.v[i - self.first as usize]
    }

    /// The median of every value, as [`median`] gives it; `None` without
    /// values.
    fn median(&self) -> Option<f64> {
        let n = self.n as usize;
        match n {
            0 => None,
            _ if n % 2 == 1 => Some(self.at(n / 2)),
            _ => Some(0.5 * (self.at(n / 2 - 1) + self.at(n / 2))),
        }
    }

    /// The median of every value but one copy of `x`, which is one of them;
    /// `None` when it is the only one. With `i` the first position of `x`
    /// in the sorted values, the `j`-th value left is the `j`-th when
    /// `j < i` (that is, when it is below `x`) and the `(j + 1)`-th
    /// otherwise.
    fn median_without(&self, x: f64) -> Option<f64> {
        let m = (self.n as usize).saturating_sub(1);
        if m == 0 {
            return None;
        }
        let left = |j: usize| {
            let v = self.at(j);
            if v.total_cmp(&x).is_lt() {
                v
            } else {
                self.at(j + 1)
            }
        };
        Some(if m % 2 == 1 {
            left(m / 2)
        } else {
            0.5 * (left(m / 2 - 1) + left(m / 2))
        })
    }

    /// [`MedianWindow::median_without`] `x` when given, else
    /// [`MedianWindow::median`].
    fn median_of(&self, x: Option<f64>) -> Option<f64> {
        match x {
            Some(x) => self.median_without(x),
            None => self.median(),
        }
    }
}

/// The donor fits of one precursor: width on the consensus scale, tailing
/// and isotope profile.
#[derive(Debug, Clone, Copy)]
struct ShapeSummary {
    width: MedianWindow,
    tailing: MedianWindow,
    iso: [MedianWindow; 4],
}

/// What finalize needs of a precursor's identified rows in every run.
#[derive(Debug, Clone, Copy, Default)]
struct PrecursorSummary {
    /// The consensus width `w_c`.
    w_c: Option<f64>,
    shape: Option<ShapeSummary>,
    /// The median amplitude signal-to-noise of the identified target fits,
    /// read only with [`FinalizeParams::trace_fraction`].
    snr: Option<f64>,
    /// Whether the precursor reports its trace areas (see
    /// [`crate::quant_mbr`], "Trace area").
    trace: bool,
}

type Summaries = HashMap<u32, PrecursorSummary>;

/// One identified row's part of its precursor's summary: precursor, width
/// on the consensus scale, donor fit and amplitude signal-to-noise.
type SummaryItem = (u32, Option<f64>, Option<DonorFit>, Option<f64>);

/// The summaries of the precursors of `items`, in any order.
fn summarise(mut items: Vec<SummaryItem>) -> Summaries {
    items.sort_unstable_by_key(|it| it.0);
    let mut out = Summaries::new();
    for g in items.chunk_by(|a, b| a.0 == b.0) {
        let widths: Vec<f64> = g.iter().filter_map(|it| it.1).collect();
        let snrs: Vec<f64> = g.iter().filter_map(|it| it.3).collect();
        let fits: Vec<&DonorFit> = g.iter().filter_map(|it| it.2.as_ref()).collect();
        let window = |value: &dyn Fn(&DonorFit) -> f64| {
            let mut v: Vec<f64> = fits.iter().map(|f| value(f)).collect();
            v.sort_by(f64::total_cmp);
            MedianWindow::new(&v)
        };
        let shape = (!fits.is_empty()).then(|| ShapeSummary {
            width: window(&|f| f.1),
            tailing: window(&|f| f.2),
            iso: [0, 1, 2, 3].map(|i| window(&|f: &DonorFit| f.3[i])),
        });
        out.insert(
            g[0].0,
            PrecursorSummary {
                w_c: (!widths.is_empty()).then(|| median(widths)),
                shape,
                snr: (!snrs.is_empty()).then(|| median(snrs)),
                trace: false,
            },
        );
    }
    out
}

/// `A / amplitude_se` of an identified target row's fit; `None` for any
/// other row, or without a positive amplitude and a finite positive
/// standard error.
fn identified_snr(r: &ExtractRow) -> Option<f64> {
    if (r.role, r.kind) != (Role::Identified, Hypothesis::Target) {
        return None;
    }
    let e = r.evidence.as_ref()?;
    (e.amplitude > 0.0 && e.amplitude_se > 0.0 && e.amplitude_se.is_finite())
        .then(|| e.amplitude / e.amplitude_se)
        .filter(|v| v.is_finite())
}

/// Mark the brightest share `fraction` (in `(0, 1]`) of the precursors with
/// a signal-to-noise, by that median, to report their trace areas: those at
/// or above the value of rank `ceil(fraction n)` from the top. Logs the cut
/// and the distribution it was taken from.
fn mark_trace_precursors(summaries: &mut Summaries, fraction: f64) {
    let mut snrs: Vec<f64> = summaries.values().filter_map(|s| s.snr).collect();
    if snrs.is_empty() {
        eprintln!("[mbr] finalize: trace area for 0 precursors (no identified fit with a signal-to-noise)");
        return;
    }
    snrs.sort_unstable_by(f64::total_cmp);
    let n = snrs.len();
    let k = ((fraction * n as f64).ceil() as usize).clamp(1, n);
    let cut = snrs[n - k];
    let mut marked = 0usize;
    for s in summaries.values_mut() {
        s.trace = s.snr.is_some_and(|v| v >= cut);
        marked += usize::from(s.trace);
    }
    let at = |q: f64| snrs[((q * (n - 1) as f64).round() as usize).min(n - 1)];
    eprintln!(
        "[mbr] finalize: trace area for {marked} of {n} precursors with an identified fit (share {fraction}): median amplitude signal-to-noise >= {cut:.1}; distribution p50 {:.1}, p90 {:.1}, p99 {:.1}",
        at(0.5),
        at(0.9),
        at(0.99)
    );
}

/// The shape reference of the rows of one block: the precursors' summaries
/// and the block's own donor fits, left out of a row's reference.
struct BlockShapes<'a> {
    summaries: &'a Summaries,
    own: &'a HashMap<(u32, u16), DonorFit>,
}

impl ShapeReference for BlockShapes<'_> {
    fn features(&self, row: &ExtractRow, maps: &[Option<RunMap>]) -> [f64; 3] {
        let Some(e) = row.evidence.as_ref().filter(|e| e.amplitude > 0.0) else {
            return NO_SHAPE;
        };
        let shape = self
            .summaries
            .get(&row.precursor)
            .and_then(|s| s.shape.as_ref());
        let own = self.own.get(&(row.precursor, row.run));
        let refs = shape
            .map_or(0, |s| s.width.n as usize)
            .saturating_sub(usize::from(own.is_some()));
        let mine = local_slope(maps, row.run, e.t0).zip(iso_profile(&e.iso_heights));
        let (Some((slope, iso)), Some(s)) = (mine, shape.filter(|_| refs > 0)) else {
            return NO_SHAPE;
        };
        let reference = |w: &MedianWindow, x: Option<f64>| {
            w.median_of(x)
                .expect("a precursor with references has a median")
        };
        let w_ref = reference(&s.width, own.map(|f| f.1));
        let tail_ref = reference(&s.tailing, own.map(|f| f.2));
        let mut iso_ref = [0.0; 4];
        for (i, v) in iso_ref.iter_mut().enumerate() {
            *v = reference(&s.iso[i], own.map(|f| f.3[i]));
        }
        [
            ((e.w / slope) / w_ref).ln(),
            e.tailing - tail_ref,
            cosine(&iso, &iso_ref),
        ]
    }
}

/// The classifier a row belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Table {
    Transfer,
    Identified,
    /// A pseudo-transfer: scored by the transfer models, outside their fit.
    Pseudo,
}

/// Most similar runs whose identifications make a precursor's presence
/// prior ([`EntryIndex::presence_prior`]), at most [`NEIGHBOUR_RUNS`].
pub(crate) const PRIOR_NEIGHBOURS: usize = 20;
/// Library entries the strata are cut from, at most (a stride over them).
const PRIOR_SAMPLE: usize = 20_000;

/// Runs compared with a transfer's run for the `neighbour_donor_rate`
/// feature: its most similar runs, at most this many.
pub(crate) const NEIGHBOUR_RUNS: usize = 50;

/// Library lookups of finalize: per precursor id, its entry; per run, its
/// most similar runs.
pub(crate) struct EntryIndex<'a> {
    entry_of: Vec<Option<&'a LibraryEntry>>,
    /// Per run, the [`NEIGHBOUR_RUNS`] other runs whose target entries'
    /// donor sets are most similar to its own (Jaccard index, ties by run
    /// index).
    neighbours: Vec<Vec<u16>>,
    /// Ascending edges of the prior strata ([`strata_edges`]).
    edges: Vec<f64>,
}

impl<'a> EntryIndex<'a> {
    /// The index of `lib`'s entries for `runs` runs.
    pub(crate) fn new(lib: &'a Library, runs: usize) -> Self {
        let mut entry_of: Vec<Option<&LibraryEntry>> = vec![None; lib.meta.len()];
        for e in &lib.entries {
            if let Some(slot) = entry_of.get_mut(e.precursor as usize) {
                *slot = Some(e);
            }
        }
        let neighbours = run_neighbours(lib, runs, NEIGHBOUR_RUNS);
        let edges = strata_edges(&entry_of, &neighbours, runs, PRIOR_NEIGHBOURS);
        Self {
            entry_of,
            neighbours,
            edges,
        }
    }

    /// Number of prior strata.
    pub(crate) fn strata(&self) -> usize {
        self.edges.len() + 1
    }

    /// The prior that the row's precursor is present in its run: the share
    /// of the run's [`PRIOR_NEIGHBOURS`] most similar runs that identified
    /// it. 0 without an entry or neighbours; the row's own run is never a
    /// neighbour. A precursor identified nowhere near this run has a prior of
    /// 0, however often it was identified in runs unlike it.
    pub(crate) fn presence_prior(&self, r: &ExtractRow) -> f64 {
        let (Some(Some(e)), Some(nn)) = (
            self.entry_of.get(r.precursor as usize),
            self.neighbours.get(r.run as usize),
        ) else {
            return 0.0;
        };
        let nn = &nn[..nn.len().min(PRIOR_NEIGHBOURS)];
        if nn.is_empty() {
            return 0.0;
        }
        let hits = nn
            .iter()
            .filter(|n| e.donors.binary_search(n).is_ok())
            .count();
        hits as f64 / nn.len() as f64
    }

    /// `family` of a row of table `t` with its prior stratum
    /// ([`pack_family`]); the stratum is 0 outside the transfer tables.
    pub(crate) fn packed(&self, r: &ExtractRow, t: Table, family: u8) -> u8 {
        match t {
            Table::Transfer | Table::Pseudo => {
                pack_family(family, stratum_of(self.presence_prior(r), &self.edges))
            }
            Table::Identified => family,
        }
    }

    /// The share of the row's run's neighbours ([`run_neighbours`]) that
    /// are donor runs of its precursor's entry: how often the runs most
    /// like this one identified the precursor. 0 without an entry or
    /// neighbours. The row's own run is never its neighbour, so a
    /// pseudo-transfer does not see its own identification.
    pub(crate) fn neighbour_rate(&self, r: &ExtractRow) -> f64 {
        let (Some(Some(e)), Some(nn)) = (
            self.entry_of.get(r.precursor as usize),
            self.neighbours.get(r.run as usize),
        ) else {
            return 0.0;
        };
        if nn.is_empty() {
            return 0.0;
        }
        let hits = nn
            .iter()
            .filter(|n| e.donors.binary_search(n).is_ok())
            .count();
        hits as f64 / nn.len() as f64
    }

    pub(crate) fn is_identity(&self, p: u32) -> bool {
        self.entry_of[p as usize].is_some_and(|e| e.identity_decoy)
    }

    /// Not transferred: the entry's best donor PSM is above TRANSFER_DONOR_Q.
    pub(crate) fn gated(&self, p: u32) -> bool {
        self.entry_of[p as usize]
            .is_some_and(|e| e.donor_psm.q.is_nan() || e.donor_psm.q > TRANSFER_DONOR_Q)
    }

    /// The classifier of `r` and its family; `None` for a row that takes no
    /// part in any.
    pub(crate) fn table(&self, r: &ExtractRow) -> Option<(Table, u8)> {
        // An identity decoy's rows take part only through its transfer
        // targets, as the identity-decoy family.
        let identity = self.is_identity(r.precursor);
        let t = match r.role {
            Role::Identified if identity => return None,
            Role::Identified => Table::Identified,
            Role::Transfer if identity && r.kind != Hypothesis::Target => return None,
            Role::Transfer if self.gated(r.precursor) => return None,
            Role::Transfer => Table::Transfer,
            Role::PseudoTransfer => Table::Pseudo,
        };
        Some((t, if identity { IDENTITY } else { family(r.kind) }))
    }
}

/// The ascending edges of the prior strata, learned from the library alone:
/// cut from the prior ([`EntryIndex::presence_prior`] with `k` neighbours)
/// of the (precursor, run) pairs without an identification, so that each of
/// at most [`MAX_STRATA`] strata holds about an equal share of them, as far
/// as the discrete prior allows. A target and its decoys share their pair,
/// so these are equal decoy counts. Taken over a stride of at most
/// [`PRIOR_SAMPLE`] entries.
pub(crate) fn strata_edges(
    entry_of: &[Option<&LibraryEntry>],
    neighbours: &[Vec<u16>],
    runs: usize,
    k: usize,
) -> Vec<f64> {
    let entries: Vec<&LibraryEntry> = entry_of
        .iter()
        .flatten()
        .copied()
        .filter(|e| !e.identity_decoy)
        .collect();
    // The neighbours every run has, at most `k`: the prior is hits / m.
    let m = neighbours.iter().map(Vec::len).min().unwrap_or(0).min(k);
    if m == 0 {
        return Vec::new();
    }
    let stride = entries.len().div_ceil(PRIOR_SAMPLE).max(1);
    // Per number of neighbours identifying the pair: how many pairs.
    let hist = entries
        .iter()
        .step_by(stride)
        .par_bridge()
        .fold(
            || vec![0u64; m + 1],
            |mut h, e| {
                for (r, nn) in neighbours.iter().enumerate().take(runs) {
                    if e.donors.binary_search(&(r as u16)).is_ok() || nn.is_empty() {
                        continue;
                    }
                    let nn = &nn[..m];
                    let hits = nn
                        .iter()
                        .filter(|n| e.donors.binary_search(n).is_ok())
                        .count();
                    h[hits] += 1;
                }
                h
            },
        )
        .reduce(
            || vec![0u64; m + 1],
            |mut a, b| {
                a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                a
            },
        );
    let edges = edges_of_hist(&hist);
    eprintln!("[mbr] presence prior strata edges {edges:?}");
    edges
}

/// The strata edges cut from `hist`, the pair counts per number of
/// neighbour hits `0..=m` (the prior is hits / m): each stratum closes once
/// it holds a share `total / MAX_STRATA`; an edge is added only while there
/// are fewer than `MAX_STRATA - 1`, so at most `MAX_STRATA - 1` edges and
/// [`MAX_STRATA`] strata.
fn edges_of_hist(hist: &[u64]) -> Vec<f64> {
    let m = hist.len().saturating_sub(1);
    let total: u64 = hist.iter().sum();
    let mut edges = Vec::new();
    if total > 0 && m > 0 {
        let share = total as f64 / MAX_STRATA as f64;
        let mut open = 0.0;
        for (i, &c) in hist.iter().enumerate().take(m) {
            open += c as f64;
            if open >= share && edges.len() + 1 < MAX_STRATA {
                edges.push(i as f64 / m as f64);
                open = 0.0;
            }
        }
        // A last stratum below one share merges into the one before it.
        if open + (hist[m] as f64) < share && !edges.is_empty() {
            edges.pop();
        }
    }
    edges
}

/// Per run of `runs`, the `k` other runs most similar to it, most similar
/// first: the Jaccard index of the sets of target entries they are donors
/// of (identity decoys are left out), ties by run index. Computed from the
/// donor runs only, the same for targets and decoys.
pub(crate) fn run_neighbours(lib: &Library, runs: usize, k: usize) -> Vec<Vec<u16>> {
    let words = lib.entries.len().div_ceil(64);
    let mut bits = vec![vec![0u64; words]; runs];
    for (i, e) in lib.entries.iter().enumerate() {
        if e.identity_decoy {
            continue;
        }
        for &r in &e.donors {
            if let Some(b) = bits.get_mut(r as usize) {
                b[i / 64] |= 1u64 << (i % 64);
            }
        }
    }
    let ones = |b: &[u64]| -> u64 { b.iter().map(|w| u64::from(w.count_ones())).sum() };
    let count: Vec<u64> = bits.iter().map(|b| ones(b)).collect();
    (0..runs)
        .into_par_iter()
        .map(|a| {
            let mut sim: Vec<(f64, u16)> = (0..runs)
                .filter(|&b| b != a)
                .map(|b| {
                    let inter: u64 = bits[a]
                        .iter()
                        .zip(&bits[b])
                        .map(|(x, y)| u64::from((x & y).count_ones()))
                        .sum();
                    let union = count[a] + count[b] - inter;
                    let j = if union == 0 {
                        0.0
                    } else {
                        inter as f64 / union as f64
                    };
                    (j, b as u16)
                })
                .collect();
            sim.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
            sim.truncate(k);
            sim.into_iter().map(|(_, b)| b).collect()
        })
        .collect()
}

/// The features of the rows of `rows` whose table passes `want`, built in
/// parallel (in row order): `(index in rows, table, family, features)`.
pub(crate) fn table_features(
    rows: &[ExtractRow],
    idx: &EntryIndex<'_>,
    shapes: &(impl ShapeReference + Sync),
    maps: &[Option<RunMap>],
    want: impl Fn(Table) -> bool + Sync,
) -> Vec<(usize, Table, u8, Option<Vec<f32>>)> {
    rows.par_iter()
        .enumerate()
        .filter_map(|(i, r)| {
            let (t, fam) = idx.table(r).filter(|(t, _)| want(*t))?;
            let x = match t {
                Table::Identified => identified_features(r),
                Table::Transfer | Table::Pseudo => {
                    transfer_features(r, shapes, maps, idx.neighbour_rate(r))
                }
            };
            Some((i, t, fam, x))
        })
        .collect()
}

/// Calls `f(index in rows, table, family, score)` for every row of `rows`
/// with features whose table passes `want`, in row order: `score` is the
/// held-out raw margin of its classifier's fold model (`None` for a fold
/// without a model). Rows are scored [`SCORE_CHUNK`] at a time.
#[allow(clippy::too_many_arguments)]
pub(crate) fn for_each_scored(
    rows: &[ExtractRow],
    idx: &EntryIndex<'_>,
    shapes: &(impl ShapeReference + Sync),
    maps: &[Option<RunMap>],
    transfer: &FoldModels,
    identified: &FoldModels,
    want: impl Fn(Table) -> bool + Sync + Copy,
    mut f: impl FnMut(usize, Table, u8, Option<f32>),
) {
    let mut out = Vec::new();
    for (c, chunk) in rows.chunks(SCORE_CHUNK).enumerate() {
        let base = c * SCORE_CHUNK;
        let feats = table_features(chunk, idx, shapes, maps, want);
        let mut scores: Vec<Option<f32>> = vec![None; feats.len()];
        for (models, is_id) in [(transfer, false), (identified, true)] {
            let at: Vec<usize> = (0..feats.len())
                .filter(|&k| feats[k].3.is_some() && ((feats[k].1 == Table::Identified) == is_id))
                .collect();
            if at.is_empty() {
                continue;
            }
            let precursors: Vec<u32> = at.iter().map(|&k| chunk[feats[k].0].precursor).collect();
            let xs: Vec<&[f32]> = at
                .iter()
                .map(|&k| feats[k].3.as_deref().expect("filtered on features"))
                .collect();
            models.score(&precursors, &xs, &mut out);
            for (&k, s) in at.iter().zip(&out) {
                scores[k] = *s;
            }
        }
        for ((i, t, fam, x), s) in feats.iter().zip(scores) {
            if x.is_some() {
                f(base + i, *t, idx.packed(&chunk[*i], *t, *fam), s);
            }
        }
    }
}

/// Metadata key and version of the temporary score tables.
const SCORES_KEY: &str = "andes.finalize_scores.schema";
/// Metadata key and version of the temporary output-row tables.
const ROWS_KEY: &str = "andes.finalize_rows.schema";
const TEMP_VERSION: &str = "1.0";
const TEMP_MAJOR: &str = "1";

/// A writer of a temporary table with schema version `key`.
fn temp_writer(
    path: &Path,
    key: &str,
    schema: Arc<Schema>,
    compression: Compression,
    group_rows: usize,
) -> Result<ArrowWriter<File>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let props = WriterProperties::builder()
        .set_compression(compression)
        .set_max_row_group_size(group_rows)
        .set_key_value_metadata(Some(vec![KeyValue::new(
            key.to_string(),
            TEMP_VERSION.to_string(),
        )]))
        .build();
    let file = File::create(path).map_err(|e| err(&e))?;
    ArrowWriter::try_new(file, schema, Some(props)).map_err(|e| err(&e))
}

/// A reader of a temporary table with schema version `key` (major 1).
fn temp_reader(
    path: &Path,
    key: &str,
    batch_rows: usize,
) -> Result<ParquetRecordBatchReader, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let file = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
    let version = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .and_then(|kv| kv.iter().find(|k| k.key == key))
        .and_then(|k| k.value.clone());
    if version.as_deref().and_then(|v| v.split('.').next()) != Some(TEMP_MAJOR) {
        return Err(err(&format!(
            "not a temporary table of schema {key} {TEMP_MAJOR}.x"
        )));
    }
    builder
        .with_batch_size(batch_rows)
        .build()
        .map_err(|e| err(&e))
}

fn column<'a, T: 'static>(b: &'a RecordBatch, i: usize, path: &Path) -> Result<&'a T, String> {
    b.column(i)
        .as_any()
        .downcast_ref::<T>()
        .ok_or_else(|| format!("{}: column {i} has an unexpected type", path.display()))
}

/// Writes scores, in the order given, as a temporary score table.
fn write_scores(path: &Path, scores: &[f32]) -> Result<(), String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let schema = Arc::new(Schema::new(vec![Field::new(
        "score",
        DataType::Float32,
        false,
    )]));
    let mut w = temp_writer(
        path,
        SCORES_KEY,
        schema.clone(),
        Compression::UNCOMPRESSED,
        SCORE_GROUP_ROWS,
    )?;
    for c in scores.chunks(SCORE_GROUP_ROWS) {
        let col: ArrayRef = Arc::new(Float32Array::from(c.to_vec()));
        let b = RecordBatch::try_new(schema.clone(), vec![col]).map_err(|e| err(&e))?;
        w.write(&b).map_err(|e| err(&e))?;
    }
    w.close().map_err(|e| err(&e))?;
    Ok(())
}

/// Whether `r` is a row whose held-out score pass 3 keeps for the later
/// passes: a transfer row (target or decoy) or an identified row.
fn is_cached(idx: &EntryIndex<'_>, r: &ExtractRow) -> bool {
    matches!(idx.table(r), Some((Table::Transfer | Table::Identified, _)))
}

/// The cached value of an identified row that has features but no score (its
/// fold has no model). A transfer row without features or without a score,
/// and an identified row without features, cache NaN.
const NO_SCORE: f32 = f32::INFINITY;

/// The temporary table of the held-out scores of the transfer and identified rows of table `chunk`.
fn cached_scores_path(dir: &Path, chunk: usize) -> PathBuf {
    dir.join(format!("transfer-scores-{chunk:06}.parquet"))
}

/// The scores of a table written by [`write_scores`], in order.
fn read_cached_scores(path: &Path) -> Result<Vec<f32>, String> {
    let error = RefCell::new(None);
    let scores: Vec<f32> = ScoreFile::new(path, &error).collect();
    match error.into_inner() {
        Some(e) => Err(e),
        None => Ok(scores),
    }
}

/// The scores of a temporary score table, in order. A read error ends the
/// stream and is kept in `error`.
struct ScoreFile<'a> {
    path: PathBuf,
    reader: Option<ParquetRecordBatchReader>,
    batch: Option<Float32Array>,
    pos: usize,
    error: &'a RefCell<Option<String>>,
}

impl<'a> ScoreFile<'a> {
    fn new(path: &Path, error: &'a RefCell<Option<String>>) -> Self {
        let reader = match temp_reader(path, SCORES_KEY, SCORE_GROUP_ROWS) {
            Ok(r) => Some(r),
            Err(e) => {
                error.borrow_mut().get_or_insert(e);
                None
            }
        };
        Self {
            path: path.to_path_buf(),
            reader,
            batch: None,
            pos: 0,
            error,
        }
    }

    fn fail(&mut self, e: String) -> Option<f32> {
        self.error.borrow_mut().get_or_insert(e);
        self.reader = None;
        None
    }
}

impl Iterator for ScoreFile<'_> {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        loop {
            if let Some(b) = &self.batch {
                if self.pos < b.len() {
                    self.pos += 1;
                    return Some(b.value(self.pos - 1));
                }
            }
            let next = self.reader.as_mut()?.next();
            match next {
                None => {
                    self.reader = None;
                    return None;
                }
                Some(Err(e)) => return self.fail(format!("{}: {e}", self.path.display())),
                Some(Ok(b)) => match column::<Float32Array>(&b, 0, &self.path) {
                    Ok(a) => {
                        self.batch = Some(a.clone());
                        self.pos = 0;
                    }
                    Err(e) => return self.fail(e),
                },
            }
        }
    }
}

/// The held-out scores of one classifier, merged run by run: in memory up
/// to a cap, then sorted and written to temporary tables.
struct ScoreStore<'d> {
    name: &'static str,
    dir: &'d Path,
    cap: usize,
    scores: Mutex<FoldScores>,
    /// Per spill, per group (fold and stratum), the table of each family.
    spills: Mutex<Vec<Vec<[PathBuf; 4]>>>,
    next: AtomicUsize,
}

impl<'d> ScoreStore<'d> {
    fn new(name: &'static str, dir: &'d Path, cap: usize, seed: u64, strata: usize) -> Self {
        Self {
            name,
            dir,
            cap,
            scores: Mutex::new(FoldScores::with_strata(seed, strata)),
            spills: Mutex::new(Vec::new()),
            next: AtomicUsize::new(0),
        }
    }

    /// Adds the scores of one run.
    fn add(&self, s: FoldScores) -> Result<(), String> {
        let held = {
            let mut all = self.scores.lock().expect("score store lock");
            all.merge(s);
            (all.held() > self.cap).then(|| all.take_held())
        };
        let Some(held) = held else { return Ok(()) };
        let k = self.next.fetch_add(1, Ordering::Relaxed);
        let mut paths = Vec::with_capacity(held.len());
        for (fold, by_family) in held.into_iter().enumerate() {
            let mut fold_paths: [PathBuf; 4] = Default::default();
            for (f, mut v) in by_family.into_iter().enumerate() {
                v.par_sort_unstable_by(f32::total_cmp);
                let p = self
                    .dir
                    .join(format!("scores-{}-{k:04}-{fold}-{f}.parquet", self.name));
                write_scores(&p, &v)?;
                fold_paths[f] = p;
            }
            paths.push(fold_paths);
        }
        self.spills.lock().expect("score store lock").push(paths);
        Ok(())
    }

    /// The q-values over every score added (the estimates other than the
    /// pooled one kept at `thresholds`), and the number of spills.
    fn finish(
        self,
        failed_folds: &[String],
        thresholds: &[f64],
    ) -> Result<(QValues, usize), String> {
        let spills = self.spills.into_inner().expect("score store lock");
        let scores = self.scores.into_inner().expect("score store lock");
        let error = RefCell::new(None);
        let q = scores.finish_with(failed_folds, Some(thresholds), |fold, f| {
            spills
                .iter()
                .map(|s| {
                    Box::new(ScoreFile::new(&s[fold][f], &error))
                        as Box<dyn Iterator<Item = f32> + '_>
                })
                .collect()
        });
        if let Some(e) = error.into_inner() {
            return Err(e);
        }
        for p in spills.iter().flatten().flatten() {
            let _ = std::fs::remove_file(p);
        }
        Ok((q, spills.len()))
    }
}

/// One output row by precursor id and run index, as the temporary tables
/// hold it.
#[derive(Debug, Clone, Copy, PartialEq)]
struct OutRow {
    precursor: u32,
    run: u16,
    evidence_type: EvidenceType,
    quantity: Option<f64>,
    quantity_se: Option<f64>,
    quantity_normalized: Option<f64>,
    transfer_q: Option<f64>,
    quant_q: Option<f64>,
    apex_rt_s: Option<f64>,
    library_q: f64,
    mass_err_ppm: Option<f64>,
    claim: Claim,
}

impl OutRow {
    fn feature(&self) -> MbrFeature {
        MbrFeature {
            precursor: self.precursor,
            run: self.run,
            evidence_type: self.evidence_type,
            quantity: self.quantity,
            quantity_se: self.quantity_se,
            quantity_normalized: self.quantity_normalized,
            transfer_q: self.transfer_q,
            quant_q: self.quant_q,
            apex_rt_s: self.apex_rt_s,
            library_q: self.library_q,
            mass_err_ppm: self.mass_err_ppm,
            claim: self.claim,
        }
    }
}

/// The evidence types by their code in the temporary tables.
const EVIDENCE_CODES: [EvidenceType; 6] = [
    EvidenceType::Identified,
    EvidenceType::IdentifiedBelowThreshold,
    EvidenceType::TransferMs2,
    EvidenceType::TransferMs1,
    EvidenceType::BelowThreshold,
    EvidenceType::SharedFeature,
];

fn evidence_code(e: EvidenceType) -> u8 {
    EVIDENCE_CODES
        .iter()
        .position(|&x| x == e)
        .expect("every evidence type has a code") as u8
}

/// The optional values of [`OutRow`], in column order after the evidence
/// type.
const OUT_VALUES: [&str; 6] = [
    "quantity",
    "quantity_se",
    "quantity_normalized",
    "transfer_q",
    "quant_q",
    "apex_rt_s",
];

/// The claim values of [`OutRow`], after the mass error: claimant
/// precursor, apex offset, claimant identified, claimant co-apex, claimants,
/// co-apex claimants, identified co-apex claimants.
const OUT_CLAIM: [&str; 7] = [
    "claimant",
    "claim_dt_s",
    "claim_identified",
    "claim_co_apex",
    "claim_n",
    "claim_n_co_apex",
    "claim_n_identified_co_apex",
];

fn out_schema() -> Arc<Schema> {
    let mut fields = vec![
        Field::new("precursor", DataType::UInt32, false),
        Field::new("run", DataType::UInt16, false),
        Field::new("evidence_type", DataType::UInt8, false),
    ];
    fields.extend(
        OUT_VALUES
            .iter()
            .map(|n| Field::new(*n, DataType::Float64, true)),
    );
    fields.push(Field::new("library_q", DataType::Float64, false));
    fields.push(Field::new("mass_err_ppm", DataType::Float64, true));
    fields.extend(
        OUT_CLAIM
            .iter()
            .map(|n| Field::new(*n, DataType::Float64, true)),
    );
    Arc::new(Schema::new(fields))
}

/// A temporary table of output rows being written, rows in the order given.
struct OutWriter {
    path: PathBuf,
    schema: Arc<Schema>,
    writer: ArrowWriter<File>,
}

impl OutWriter {
    fn create(path: &Path) -> Result<Self, String> {
        let schema = out_schema();
        let writer = temp_writer(
            path,
            ROWS_KEY,
            schema.clone(),
            Compression::ZSTD(ZstdLevel::try_new(1).map_err(|e| e.to_string())?),
            OUT_GROUP_ROWS,
        )?;
        Ok(Self {
            path: path.to_path_buf(),
            schema,
            writer,
        })
    }

    fn write(&mut self, rows: &[OutRow]) -> Result<(), String> {
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", self.path.display());
        for c in rows.chunks(OUT_GROUP_ROWS) {
            let mut precursor = UInt32Builder::with_capacity(c.len());
            let mut run = UInt16Builder::with_capacity(c.len());
            let mut ev = UInt8Builder::with_capacity(c.len());
            let mut values: Vec<Float64Builder> = (0..OUT_VALUES.len() + 2 + OUT_CLAIM.len())
                .map(|_| Float64Builder::with_capacity(c.len()))
                .collect();
            for r in c {
                precursor.append_value(r.precursor);
                run.append_value(r.run);
                ev.append_value(evidence_code(r.evidence_type));
                let v = [
                    r.quantity,
                    r.quantity_se,
                    r.quantity_normalized,
                    r.transfer_q,
                    r.quant_q,
                    r.apex_rt_s,
                    Some(r.library_q),
                    r.mass_err_ppm,
                    r.claim.claimant.map(f64::from),
                    r.claim.dt_s.map(f64::from),
                    Some(f64::from(u8::from(r.claim.identified))),
                    Some(f64::from(u8::from(r.claim.co_apex))),
                    Some(f64::from(r.claim.n)),
                    Some(f64::from(r.claim.n_co_apex)),
                    Some(f64::from(r.claim.n_identified_co_apex)),
                ];
                for (b, x) in values.iter_mut().zip(v) {
                    b.append_option(x);
                }
            }
            let mut cols: Vec<ArrayRef> = vec![
                Arc::new(precursor.finish()),
                Arc::new(run.finish()),
                Arc::new(ev.finish()),
            ];
            cols.extend(values.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
            let b = RecordBatch::try_new(self.schema.clone(), cols).map_err(|e| err(&e))?;
            self.writer.write(&b).map_err(|e| err(&e))?;
        }
        Ok(())
    }

    fn finish(self) -> Result<(), String> {
        let path = self.path;
        self.writer
            .close()
            .map(|_| ())
            .map_err(|e| format!("{}: {e}", path.display()))
    }
}

/// A temporary table of output rows, read in order.
struct OutReader {
    path: PathBuf,
    reader: ParquetRecordBatchReader,
    rows: std::vec::IntoIter<OutRow>,
}

impl OutReader {
    fn open(path: &Path) -> Result<Self, String> {
        Ok(Self {
            path: path.to_path_buf(),
            reader: temp_reader(path, ROWS_KEY, OUT_GROUP_ROWS)?,
            rows: Vec::new().into_iter(),
        })
    }

    fn next_row(&mut self) -> Result<Option<OutRow>, String> {
        loop {
            if let Some(r) = self.rows.next() {
                return Ok(Some(r));
            }
            let Some(b) = self.reader.next() else {
                return Ok(None);
            };
            let b = b.map_err(|e| format!("{}: {e}", self.path.display()))?;
            self.rows = self.decode(&b)?.into_iter();
        }
    }

    fn decode(&self, b: &RecordBatch) -> Result<Vec<OutRow>, String> {
        let path = self.path.as_path();
        let expected = OUT_VALUES.len() + 5 + OUT_CLAIM.len();
        if b.num_columns() != expected {
            return Err(format!(
                "{}: {} columns where {} were expected",
                path.display(),
                b.num_columns(),
                expected
            ));
        }
        let precursor = column::<UInt32Array>(b, 0, path)?;
        let run = column::<UInt16Array>(b, 1, path)?;
        let ev = column::<UInt8Array>(b, 2, path)?;
        let values: Vec<&Float64Array> = (3..b.num_columns())
            .map(|i| column::<Float64Array>(b, i, path))
            .collect::<Result<_, _>>()?;
        let opt = |k: usize, i: usize| (!values[k].is_null(i)).then(|| values[k].value(i));
        (0..b.num_rows())
            .map(|i| {
                let evidence_type = *EVIDENCE_CODES
                    .get(ev.value(i) as usize)
                    .ok_or_else(|| format!("{}: unknown evidence code", path.display()))?;
                Ok(OutRow {
                    precursor: precursor.value(i),
                    run: run.value(i),
                    evidence_type,
                    quantity: opt(0, i),
                    quantity_se: opt(1, i),
                    quantity_normalized: opt(2, i),
                    transfer_q: opt(3, i),
                    quant_q: opt(4, i),
                    apex_rt_s: opt(5, i),
                    library_q: values[6].value(i),
                    mass_err_ppm: opt(7, i),
                    claim: Claim {
                        claimant: opt(8, i).map(|x| x as u32),
                        dt_s: opt(9, i).map(|x| x as f32),
                        identified: opt(10, i).is_some_and(|x| x != 0.0),
                        co_apex: opt(11, i).is_some_and(|x| x != 0.0),
                        n: opt(12, i).map_or(0, |x| x as u16),
                        n_co_apex: opt(13, i).map_or(0, |x| x as u16),
                        n_identified_co_apex: opt(14, i).map_or(0, |x| x as u16),
                        ..Claim::default()
                    },
                })
            })
            .collect()
    }
}

/// Merges the temporary output tables `paths` (each in (precursor, run)
/// order, no key in two of them) and calls `sink` with the rows in
/// (precursor, run) order, in batches. More than [`MAX_OPEN_PARTS`] tables
/// are merged in groups into tables of `dir` first (named by `level`).
fn merge_out(
    paths: &[PathBuf],
    dir: &Path,
    level: usize,
    sink: &mut dyn FnMut(&[OutRow]) -> Result<(), String>,
) -> Result<(), String> {
    if paths.len() > MAX_OPEN_PARTS {
        let mut merged = Vec::new();
        for (g, group) in paths.chunks(MAX_OPEN_PARTS).enumerate() {
            let p = dir.join(format!("rows-merged-{level}-{g:05}.parquet"));
            let mut w = OutWriter::create(&p)?;
            merge_out(group, dir, level + 1, &mut |rows| w.write(rows))?;
            w.finish()?;
            for q in group {
                let _ = std::fs::remove_file(q);
            }
            merged.push(p);
        }
        return merge_out(&merged, dir, level + 1, sink);
    }
    let mut readers: Vec<OutReader> = paths
        .iter()
        .map(|p| OutReader::open(p))
        .collect::<Result<_, _>>()?;
    let mut heads: Vec<Option<OutRow>> = Vec::with_capacity(readers.len());
    let mut heap = std::collections::BinaryHeap::new();
    for (i, r) in readers.iter_mut().enumerate() {
        let head = r.next_row()?;
        if let Some(h) = &head {
            heap.push(std::cmp::Reverse((h.precursor, h.run, i)));
        }
        heads.push(head);
    }
    let mut out = Vec::with_capacity(OUT_GROUP_ROWS);
    while let Some(std::cmp::Reverse((_, _, i))) = heap.pop() {
        out.push(heads[i].take().expect("a queued table has a head row"));
        if let Some(h) = readers[i].next_row()? {
            heap.push(std::cmp::Reverse((h.precursor, h.run, i)));
            heads[i] = Some(h);
        }
        if out.len() == OUT_GROUP_ROWS {
            sink(&out)?;
            out.clear();
        }
    }
    if !out.is_empty() {
        sink(&out)?;
    }
    Ok(())
}

/// Transfer FDR, identified-row quality q-values, pseudo-transfer validation
/// and the MBR outputs over the extract rows of `runs` (one vector per run;
/// a row's `run` indexes `stems`, its `precursor` indexes `lib.meta`), held
/// in memory, without sample caps; [`FinalizeSummary::features`] gets the
/// rows of `<base>.mbr.parquet`. See [`finalize_stream`].
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn finalize(
    stems: &[String],
    lib: &Library,
    runs: &[Vec<ExtractRow>],
    mbr_fdr: f64,
    quant_feature_fdr: Option<f64>,
    seed: u64,
    mbr_path: &Path,
    validation_path: &Path,
) -> Result<FinalizeSummary, String> {
    let params = FinalizeParams {
        mbr_fdr,
        quant_feature_fdr,
        seed,
        limits: FinalizeLimits {
            max_train_rows: usize::MAX,
            transfer_models: 1,
            max_norm_cells: usize::MAX,
            max_donor_fits: usize::MAX,
            max_held_scores: usize::MAX,
            temp_dir: None,
        },
        mbr_path,
        validation_path,
        trace_fraction: 0.0,
    };
    let source = ExtractRuns {
        runs,
        batch_rows: usize::MAX,
    };
    let mut features = Vec::new();
    let mut summary = finalize_stream(stems, lib, &source, &params, &mut |f| {
        features.extend(f);
        Ok(())
    })?;
    summary.features = features;
    Ok(summary)
}

/// What [`finalize_stream`] computes and where it writes.
pub(crate) struct FinalizeParams<'a> {
    /// Transfer q-value threshold of the reported transfers.
    pub mbr_fdr: f64,
    /// Identified-row quality threshold (`identified_below_threshold`).
    pub quant_feature_fdr: Option<f64>,
    /// Seed of the folds, the models and the samples.
    pub seed: u64,
    pub limits: FinalizeLimits,
    pub mbr_path: &'a Path,
    pub validation_path: &'a Path,
    /// Share, in `[0, 1]`, of the precursors with an identified fit that
    /// report their trace areas, the brightest by median amplitude
    /// signal-to-noise (see [`crate::quant_mbr`], "Trace area"); 0 is off.
    pub trace_fraction: f64,
}

/// Number of features of the wrong-peak model.
const WRONG_FEATURES: usize = 11;
/// Wrong-peak training and calibration rows kept at most (calibration pairs
/// are thinned by whole units above it).
const WRONG_CALIBRATION_ROWS: usize = 4_000_000;

/// Features of the wrong-peak model: the fit evidence only (no donor count,
/// library value or MS2 evidence), `None` with a non-finite feature or no
/// positive quantity. The width ratio and the two standard errors are
/// floored with `f64::max` first, which also turns a NaN among them into the
/// floor rather than a rejection.
pub(crate) fn wrong_features(e: &Evidence) -> Option<Vec<f32>> {
    if e.quantity.is_nan() || e.quantity <= 0.0 {
        return None;
    }
    let x = [
        e.llr,
        e.rt_z.abs(),
        e.width_ratio.max(1e-3).ln().abs(),
        e.envelope_cos,
        e.coherence,
        e.m1_ratio,
        e.mass_err_ppm.abs(),
        e.n_points as f64,
        (e.amplitude / e.amplitude_se.max(1e-9) + 1.0).log10(),
        (e.quantity_se / e.quantity).max(1e-6).ln(),
        e.quantity.max(1e-3).log2(),
    ];
    debug_assert_eq!(x.len(), WRONG_FEATURES);
    let x: Vec<f32> = x.iter().map(|&v| v as f32).collect();
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// quant_q of each `(p_wrong, row)` in `pool`: the mean `p_wrong` of the rows
/// at or below its own `p_wrong`. Rows with equal `p_wrong` are one step and
/// share the mean at the end of their group, so the result does not depend on
/// the order of the rows. `set(row, q)` receives every row once; `pool` is
/// sorted in place.
pub(crate) fn running_mean_q(pool: &mut [(f64, usize)], mut set: impl FnMut(usize, f64)) {
    pool.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut sum = 0.0;
    let mut start = 0;
    while start < pool.len() {
        let p = pool[start].0;
        let end = start + pool[start..].partition_point(|r| r.0.total_cmp(&p).is_eq());
        sum += pool[start..end].iter().map(|r| r.0).sum::<f64>();
        let q = sum / end as f64;
        for &(_, i) in &pool[start..end] {
            set(i, q);
        }
        start = end;
    }
}

/// Whether a pseudo-transfer's measurement is wrong against the same cell's
/// identified fit: the apex is further than the identified width from the
/// identified apex, or the quantity is off by more than a factor of 2.
fn wrong_label(pseudo: &Evidence, id: &Evidence) -> Option<bool> {
    let ok = pseudo.t0.is_finite()
        && id.t0.is_finite()
        && id.w.is_finite()
        && pseudo.quantity > 0.0
        && id.quantity > 0.0;
    ok.then(|| {
        (pseudo.t0 - id.t0).abs() > id.w || (pseudo.quantity / id.quantity).log2().abs() > 1.0
    })
}

/// A non-increasing step function of the classifier margin: the isotonic
/// estimate of the probability that a measurement is wrong.
pub(crate) struct WrongCalibration {
    /// Lowest margin of each block, ascending.
    lo: Vec<f32>,
    /// Probability of each block, non-increasing.
    p: Vec<f64>,
}

impl WrongCalibration {
    /// Pool-adjacent-violators fit on `(margin, wrong)` pairs; tied margins
    /// are one point. `None` without pairs.
    pub(crate) fn fit(mut pairs: Vec<(f32, bool)>) -> Option<Self> {
        if pairs.is_empty() {
            return None;
        }
        pairs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        // Distinct margins, descending: (margin, wrong count, rows).
        let mut pts: Vec<(f32, f64, f64)> = Vec::new();
        for &(s, w) in pairs.iter().rev() {
            match pts.last_mut() {
                Some(l) if l.0 == s => {
                    l.1 += f64::from(u8::from(w));
                    l.2 += 1.0;
                }
                _ => pts.push((s, f64::from(u8::from(w)), 1.0)),
            }
        }
        // Blocks with non-decreasing mean as the margin descends.
        let mut blocks: Vec<(f32, f64, f64)> = Vec::new();
        for pt in pts {
            blocks.push(pt);
            while blocks.len() > 1 {
                let n = blocks.len();
                let (a, b) = (blocks[n - 2], blocks[n - 1]);
                if a.1 / a.2 <= b.1 / b.2 {
                    break;
                }
                blocks.truncate(n - 2);
                blocks.push((b.0, a.1 + b.1, a.2 + b.2));
            }
        }
        blocks.reverse();
        Some(Self {
            lo: blocks.iter().map(|b| b.0).collect(),
            p: blocks.iter().map(|b| b.1 / b.2).collect(),
        })
    }

    /// The probability at margin `s`: that of the block of the nearest
    /// margin at or below `s` (the highest one below every margin seen).
    pub(crate) fn at(&self, s: f32) -> f64 {
        let k = self.lo.partition_point(|&b| b <= s);
        self.p[k.saturating_sub(1)]
    }
}

/// A pseudo-transfer target's wrong-peak label and features; `None` for any
/// other row, an identity decoy's row, or a row without an identified fit of
/// the same cell.
pub(crate) fn pseudo_example(
    r: &ExtractRow,
    ctx: &BlockContext<'_>,
    idx: &EntryIndex<'_>,
) -> Option<(bool, Vec<f32>)> {
    if r.role != Role::PseudoTransfer
        || r.kind != Hypothesis::Target
        || idx.is_identity(r.precursor)
    {
        return None;
    }
    let e = r.evidence.as_ref()?;
    let id = ctx.identified.get(&(r.precursor, r.run))?;
    Some((wrong_label(e, id)?, wrong_features(e)?))
}

/// The output rows and validation counts of one run.
#[derive(Default)]
struct RunOutput {
    rows: Vec<OutRow>,
    /// Per transfer row: its index in `rows` and its probability of being a
    /// wrong measurement.
    wrong: Vec<(usize, f64)>,
    validation: MbrValidation,
}

/// Seconds since `t`, for the pass log.
fn secs(t: Instant) -> f64 {
    t.elapsed().as_secs_f64()
}

/// Transfer FDR, identified-row quality q-values, pseudo-transfer validation
/// and the MBR outputs over the extract rows of `source` (see the module
/// docs): `<base>.mbr.parquet` in `(precursor, run)` order, its rows handed
/// to `on_features` in the same order, and `<base>.mbr_validation.parquet`.
/// With `quant_feature_fdr`, identified rows whose quality q-value is above
/// it, or missing, are `identified_below_threshold`.
///
/// Errors: `mbr_fdr` or `quant_feature_fdr` outside `(0, 1]`, a row whose run
/// or precursor is out of range, a precursor whose rows of a run are not
/// contiguous, two identified rows of one (precursor, run), or an unreadable
/// or unwritable file (the message names it).
pub(crate) fn finalize_stream(
    stems: &[String],
    lib: &Library,
    source: &dyn ExtractSource,
    params: &FinalizeParams<'_>,
    on_features: &mut dyn FnMut(Vec<MbrFeature>) -> Result<(), String>,
) -> Result<FinalizeSummary, String> {
    let (mbr_fdr, quant_feature_fdr, seed) =
        (params.mbr_fdr, params.quant_feature_fdr, params.seed);
    let limits = &params.limits;
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
    let trace_fraction = params.trace_fraction;
    if !(0.0..=1.0).contains(&trace_fraction) {
        return Err(format!(
            "trace-area precursor share {trace_fraction} is outside [0, 1]"
        ));
    }
    let maps = &lib.alignment.maps;
    let idx = EntryIndex::new(lib, stems.len());
    let bounds = Bounds {
        runs: stems.len(),
        precursors: lib.meta.len(),
    };
    let tmp_parent = match &limits.temp_dir {
        Some(d) => d.as_path(),
        None => params
            .mbr_path
            .parent()
            .filter(|d| !d.as_os_str().is_empty())
            .unwrap_or(Path::new(".")),
    };
    let tmp = tempfile::Builder::new()
        .prefix(".andes-finalize-")
        .tempdir_in(tmp_parent)
        .map_err(|e| format!("temporary directory in {}: {e}", tmp_parent.display()))?;
    let dir = tmp.path();
    let mut summary = FinalizeSummary::default();

    // Pass 1: donor summaries, in groups of precursors within the fit cap.
    let t0 = Instant::now();
    let upper: usize = lib.entries.iter().map(|e| e.donors.len()).sum();
    let groups = upper.div_ceil(limits.max_donor_fits.max(1)).max(1);
    let mut all_summaries = Summaries::new();
    for g in 0..groups {
        let items = Mutex::new(Vec::new());
        for_each_block(
            source,
            true,
            bounds,
            |_| Vec::new(),
            |local, _, _, rows| {
                for r in rows.iter().filter(|r| r.precursor as usize % groups == g) {
                    let w = consensus_scale_width(r, maps);
                    let f = donor_fit(r, maps);
                    let snr = (trace_fraction > 0.0).then(|| identified_snr(r)).flatten();
                    if w.is_some() || f.is_some() || snr.is_some() {
                        local.push((r.precursor, w, f, snr));
                    }
                }
                Ok(())
            },
            |_, local| {
                items.lock().expect("summary lock").push(local);
                Ok(())
            },
        )?;
        // One allocation of the exact size, the run vectors freed as they go.
        let mut runs: Vec<Vec<_>> = items.into_inner().expect("summary lock");
        let mut all = Vec::with_capacity(runs.iter().map(Vec::len).sum());
        while let Some(v) = runs.pop() {
            all.extend(v);
        }
        all_summaries.extend(summarise(all));
    }
    if trace_fraction > 0.0 {
        mark_trace_precursors(&mut all_summaries, trace_fraction);
    }
    eprintln!(
        "[mbr] finalize: {} extract tables, {} threads, temporary tables in {}; donor summaries of {} precursors in {groups} pass(es), {:.1} s",
        source.chunks(),
        rayon::current_num_threads(),
        dir.display(),
        all_summaries.len(),
        secs(t0)
    );
    let summaries = &all_summaries;

    // Pass 2: the training rows of both classifiers (one sample per member
    // of the transfer classifier), keyed by table and row so that they enter
    // training in table order.
    let t0 = Instant::now();
    type Key = (u32, u64);
    let members = limits.transfer_models.max(1);
    let t_samples: Vec<Mutex<TrainingSample<Key>>> = (0..members)
        .map(|m| {
            Mutex::new(TrainingSample::new(
                limits.max_train_rows,
                member_seed(seed, m),
            ))
        })
        .collect();
    let i_sample: Mutex<TrainingSample<Key>> =
        Mutex::new(TrainingSample::new(limits.max_train_rows, seed));
    let w_sample: Mutex<TrainingSample<Key>> =
        Mutex::new(TrainingSample::new(limits.max_train_rows, seed));
    let w_seen = AtomicUsize::new(0);
    let no_features = AtomicUsize::new(0);
    let identified_no_features = AtomicUsize::new(0);
    let rows_read = AtomicUsize::new(0);
    for_each_block(
        source,
        false,
        bounds,
        |_| (),
        |_, chunk, first, rows| {
            rows_read.fetch_add(rows.len(), Ordering::Relaxed);
            let ctx = BlockContext::new(rows, maps)?;
            let shapes = BlockShapes {
                summaries,
                own: &ctx.own,
            };
            let feats = table_features(rows, &idx, &shapes, maps, |t: Table| t != Table::Pseudo);
            let (mut t_rows, mut i_rows) = (Vec::new(), Vec::new());
            for (i, t, fam, x) in feats {
                let Some(x) = x else {
                    match t {
                        Table::Identified => &identified_no_features,
                        _ => &no_features,
                    }
                    .fetch_add(1, Ordering::Relaxed);
                    continue;
                };
                let r = &rows[i];
                let item = (
                    (chunk as u32, first + i as u64),
                    row_unit(r.precursor, r.run),
                    FdrRow {
                        precursor: r.precursor,
                        is_decoy: fam != 0,
                        family: fam,
                        x,
                    },
                );
                match t {
                    Table::Identified => i_rows.push(item),
                    _ => t_rows.push(item),
                }
            }
            for sample in &t_samples {
                let mut s = sample.lock().expect("training sample lock");
                for (k, u, row) in &t_rows {
                    if s.wants(*u) {
                        s.push(*k, *u, row.clone());
                    } else {
                        s.skip();
                    }
                }
            }
            let mut s = i_sample.lock().expect("training sample lock");
            for (k, u, row) in i_rows {
                s.push(k, u, row);
            }
            drop(s);
            // The wrong-peak model learns from the pseudo-transfers: each
            // judged against the identified fit of the same cell.
            let mut w_rows = Vec::new();
            for (i, r) in rows.iter().enumerate() {
                let Some((y, x)) = pseudo_example(r, &ctx, &idx) else {
                    continue;
                };
                w_rows.push((
                    (chunk as u32, first + i as u64),
                    row_unit(r.precursor, r.run),
                    FdrRow {
                        precursor: r.precursor,
                        is_decoy: y,
                        family: 0,
                        x,
                    },
                ));
            }
            w_seen.fetch_add(w_rows.len(), Ordering::Relaxed);
            let mut s = w_sample.lock().expect("training sample lock");
            for (k, u, row) in w_rows {
                s.push(k, u, row);
            }
            Ok(())
        },
        |_, ()| Ok(()),
    )?;
    summary.no_features = no_features.into_inner();
    summary.identified_no_features = identified_no_features.into_inner();
    let t_samples: Vec<TrainingSample<Key>> = t_samples
        .into_iter()
        .map(|s| s.into_inner().expect("training sample lock"))
        .collect();
    let i_sample = i_sample.into_inner().expect("training sample lock");
    for (name, s) in [
        ("transfer", &t_samples[0]),
        ("identified quality", &i_sample),
    ] {
        let (seen, kept, share) = s.stats();
        if kept < seen {
            eprintln!(
                "[mbr] {name} classifier: trained on {kept} of {seen} rows ((precursor, run) pairs, {share} of them; every row is scored)"
            );
        }
    }
    let t1 = Instant::now();
    // Both classifiers' folds train at once.
    let t_sets: Vec<Vec<FdrRow>> = t_samples
        .into_iter()
        .map(TrainingSample::into_items)
        .collect();
    let i_rows = i_sample.into_items();
    let w_rows = w_sample
        .into_inner()
        .expect("training sample lock")
        .into_items();
    let t_refs: Vec<&[FdrRow]> = t_sets.iter().map(Vec::as_slice).collect();
    let (t_models, i_models) = rayon::join(
        || train_bagged_folds(&t_refs, seed, &transfer_train_params()),
        || train_folds(&i_rows, seed),
    );
    drop(t_refs);
    drop((t_sets, i_rows));
    let w_models = train_folds(&w_rows, seed);
    eprintln!(
        "[mbr] wrong-peak model: {} pseudo-transfers ({:.1} % wrong against the identified fit), {} failed folds",
        w_rows.len(),
        100.0 * w_rows.iter().filter(|r| r.is_decoy).count() as f64 / w_rows.len().max(1) as f64,
        w_models.failed_folds.len()
    );
    drop(w_rows);
    eprintln!(
        "[mbr] finalize: training rows of {} extract rows {:.1} s, fold models {:.1} s ({} transfer classifier members per fold)",
        rows_read.into_inner(),
        (t1 - t0).as_secs_f64(),
        secs(t1),
        t_models.members()
    );

    // Pass 3: every row's held-out score, into the q-value step functions;
    // each run's transfer scores are also kept in a temporary table for
    // pass 5.
    let t0 = Instant::now();
    let w_stride = (w_seen.into_inner() / WRONG_CALIBRATION_ROWS + 1) as u64;
    let w_held: Mutex<Vec<(f32, bool)>> = Mutex::new(Vec::new());
    let strata = idx.strata();
    let t_store = ScoreStore::new("transfer", dir, limits.max_held_scores, seed, strata);
    let i_store = ScoreStore::new("identified", dir, limits.max_held_scores, seed, 1);
    for_each_block(
        source,
        false,
        bounds,
        |_| {
            (
                FoldScores::with_strata(seed, strata),
                FoldScores::new(seed),
                Vec::new(),
            )
        },
        |(ts, is, cache), _, _, rows| {
            let ctx = BlockContext::new(rows, maps)?;
            let shapes = BlockShapes {
                summaries,
                own: &ctx.own,
            };
            let want = |t: Table| t != Table::Pseudo;
            let mut block = vec![f32::NAN; rows.len()];
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
                        Table::Identified => {
                            is.add(p, fam, s);
                            block[i] = s.unwrap_or(NO_SCORE);
                        }
                        _ => {
                            ts.add(p, fam, s);
                            block[i] = s.unwrap_or(f32::NAN);
                        }
                    }
                },
            );
            cache.extend(
                rows.iter()
                    .zip(block)
                    .filter(|(r, _)| is_cached(&idx, r))
                    .map(|(_, s)| s),
            );
            // Held-out margins of the pseudo-transfers, for the calibration
            // of the wrong-peak model (whole units thinned above the cap).
            let mut labels = Vec::new();
            let mut xs = Vec::new();
            let mut precursors = Vec::new();
            for r in rows {
                if sample_hash(row_unit(r.precursor, r.run), seed) % w_stride != 0 {
                    continue;
                }
                if let Some((y, x)) = pseudo_example(r, &ctx, &idx) {
                    labels.push(y);
                    xs.push(x);
                    precursors.push(r.precursor);
                }
            }
            let x_refs: Vec<&[f32]> = xs.iter().map(Vec::as_slice).collect();
            let mut margins = Vec::new();
            w_models.score(&precursors, &x_refs, &mut margins);
            w_held.lock().expect("wrong-peak calibration lock").extend(
                margins
                    .iter()
                    .zip(&labels)
                    .filter_map(|(m, &y)| m.map(|m| (m, y))),
            );
            Ok(())
        },
        |chunk, (ts, is, cache)| {
            // One table's transfer scores in row order (NaN without a score).
            write_scores(&cached_scores_path(dir, chunk), &cache)?;
            t_store.add(ts)?;
            i_store.add(is)
        },
    )?;
    let wrong_cal =
        WrongCalibration::fit(w_held.into_inner().expect("wrong-peak calibration lock"));
    match &wrong_cal {
        Some(c) => eprintln!(
            "[mbr] wrong-peak calibration: {} blocks, probability {:.4} at the best margin, {:.4} at the worst",
            c.p.len(),
            c.p.last().copied().unwrap_or(f64::NAN),
            c.p.first().copied().unwrap_or(f64::NAN)
        ),
        None => eprintln!("[mbr] wrong-peak calibration: no pseudo-transfers; every transfer's quant_q is 1"),
    }
    // The validation thresholds, VALIDATION_Q then mbr_fdr: the transfer
    // estimates other than the pooled one are only counted at them, and
    // the identified rows use the pooled estimate only.
    let thresholds = [VALIDATION_Q, mbr_fdr];
    let t1 = Instant::now();
    let (t_done, i_done) = rayon::join(
        || t_store.finish(&t_models.failed_folds, &thresholds),
        || i_store.finish(&i_models.failed_folds, &[]),
    );
    let ((fdr, t_spills), (id_fdr, i_spills)) = (t_done?, i_done?);
    eprintln!(
        "[mbr] finalize: held-out scores {:.1} s, step functions {:.1} s ({} and {} steps, {} and {} pooled; {t_spills} and {i_spills} score spills)",
        (t1 - t0).as_secs_f64(),
        secs(t1),
        fdr.steps(),
        id_fdr.steps(),
        fdr.pooled_steps(),
        id_fdr.pooled_steps()
    );

    // Pass 4: the reference cells of the run normalisation, from the
    // identified rows' scores of pass 3.
    let t0 = Instant::now();
    let cells: Mutex<TrainingSample<Key, NormCell>> =
        Mutex::new(TrainingSample::new(limits.max_norm_cells, seed));
    for_each_block(
        source,
        false,
        bounds,
        |_| (None::<Vec<f32>>, 0usize),
        |(cache, at), chunk, first, rows| {
            if cache.is_none() {
                *cache = Some(read_cached_scores(&cached_scores_path(dir, chunk))?);
            }
            let cache = cache.as_deref().unwrap_or_default();
            let mut local = Vec::new();
            for (i, r) in rows.iter().enumerate() {
                let Some((t @ (Table::Transfer | Table::Identified), fam)) = idx.table(r) else {
                    continue;
                };
                let s = *cache.get(*at).ok_or_else(|| {
                    format!("{}: fewer cached scores than rows", source.name(chunk))
                })?;
                *at += 1;
                if t != Table::Identified || s.is_nan() {
                    continue;
                }
                let s = (s != NO_SCORE).then_some(s);
                if r.kind != Hypothesis::Target
                    || !id_fdr
                        .q(Estimate::Pooled, r.precursor, fam, s)
                        .is_some_and(|q| q <= VALIDATION_Q)
                {
                    continue;
                }
                // A fit without signal reports no quantity (see the outputs
                // below), so it is no reference cell either.
                let Some(e) = r.evidence.as_ref().filter(|e| e.amplitude > 0.0) else {
                    continue;
                };
                let s = summaries.get(&r.precursor);
                let (w_c, trace) = (s.and_then(|s| s.w_c), s.is_some_and(|s| s.trace));
                let q = reported_quantity(e, r.run, w_c, trace, maps).0;
                if q > 0.0 && q.is_finite() {
                    local.push((
                        (chunk as u32, first + i as u64),
                        u64::from(r.precursor),
                        NormCell {
                            precursor: r.precursor,
                            run: r.run,
                            log2_quantity: q.log2(),
                            rt: e.t0,
                        },
                    ));
                }
            }
            let mut c = cells.lock().expect("normalisation cells lock");
            for (k, u, cell) in local {
                c.push(k, u, cell);
            }
            Ok(())
        },
        |_, _| Ok(()),
    )?;
    let cells = cells.into_inner().expect("normalisation cells lock");
    let (seen, kept, share) = cells.stats();
    if kept < seen {
        eprintln!(
            "[mbr] run normalisation: fitted on {kept} of {seen} reference cells (whole precursors, {share} of them)"
        );
    }
    let t1 = Instant::now();
    let norm = fit_normalization(&cells.into_items(), stems.len(), &NormParams::default());
    eprintln!(
        "[mbr] finalize: normalisation cells {:.1} s, fit {:.1} s",
        (t1 - t0).as_secs_f64(),
        secs(t1)
    );

    // Pass 5: per run, the output rows into a temporary table, and the
    // validation counts. The transfer and identified rows' scores are those
    // of pass 3 (a transfer row without features or without a fold model,
    // and an identified row without features, has none); the pseudo-transfer
    // rows are scored here.
    let t0 = Instant::now();
    let shared_feature = shared::Model::shipped();
    // Decoy transfer rows with a q-value, and their acceptance, by family:
    // mass shift, RT shift, identity decoy.
    let count_decoy = |v: &mut MbrValidation, q: Option<f64>, f: usize| {
        let Some(q) = q else { return };
        v.decoy_transfers[f] += 1;
        for (a, &t) in [&mut v.at_q01, &mut v.at_threshold]
            .into_iter()
            .zip(&thresholds)
        {
            a.decoys[f] += usize::from(q <= t);
        }
    };
    let parts: Mutex<Vec<(usize, PathBuf)>> = Mutex::new(Vec::new());
    let totals: Mutex<MbrValidation> = Mutex::new(MbrValidation::default());
    for_each_block(
        source,
        false,
        bounds,
        |_| (RunOutput::default(), None, 0usize),
        |(out, cache, at), chunk, _, rows| {
            if cache.is_none() {
                *cache = Some(read_cached_scores(&cached_scores_path(dir, chunk))?);
            }
            let cache = cache.as_deref().unwrap_or_default();
            let ctx = BlockContext::new(rows, maps)?;
            let shapes = BlockShapes {
                summaries,
                own: &ctx.own,
            };
            let v = &mut out.validation;
            // Per row with a score: its family and held-out score.
            let mut scored: Vec<Option<(u8, Option<f32>)>> = vec![None; rows.len()];
            for (i, r) in rows.iter().enumerate() {
                let Some((t @ (Table::Transfer | Table::Identified), fam)) = idx.table(r) else {
                    continue;
                };
                let s = *cache.get(*at).ok_or_else(|| {
                    format!("{}: fewer cached scores than rows", source.name(chunk))
                })?;
                *at += 1;
                let fam = idx.packed(r, t, fam);
                if s.is_finite() {
                    scored[i] = Some((fam, Some(s)));
                } else if t == Table::Identified && s == NO_SCORE {
                    scored[i] = Some((fam, None));
                }
            }
            for_each_scored(
                rows,
                &idx,
                &shapes,
                maps,
                &t_models,
                &i_models,
                |t: Table| t == Table::Pseudo,
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
                let s = summaries.get(&r.precursor);
                let (w_c, trace) = (s.and_then(|s| s.w_c), s.is_some_and(|s| s.trace));
                r.evidence
                    .as_ref()
                    .filter(|e| e.amplitude > 0.0)
                    .map(|e| reported_quantity(e, r.run, w_c, trace, maps))
            };
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
                let x_refs: Vec<&[f32]> = ok
                    .iter()
                    .map(|&k| xs[k].as_deref().expect("filtered on features"))
                    .collect();
                let mut margins = Vec::new();
                w_models.score(&precursors, &x_refs, &mut margins);
                for (&k, m) in ok.iter().zip(margins) {
                    if let (Some(m), Some(c)) = (m, &wrong_cal) {
                        p_wrong[at[k]] = c.at(m);
                    }
                }
            }
            for (i, r) in rows.iter().enumerate() {
                if idx.is_identity(r.precursor) {
                    // Never reported: only the identity-decoy family's counts.
                    if (r.role, r.kind) == (Role::Transfer, Hypothesis::Target) {
                        count_decoy(v, q_of(&fdr, Estimate::Pooled, i), 2);
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
                        count_decoy(v, q_of(&fdr, Estimate::Pooled, i), f - 1);
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
                        let ok = |e: Estimate, t: f64| {
                            usize::from(
                                scored[i]
                                    .is_some_and(|(fam, s)| fdr.accepts(e, r.precursor, fam, s, t)),
                            )
                        };
                        for (a, &t) in [&mut v.at_q01, &mut v.at_threshold]
                            .into_iter()
                            .zip(&thresholds)
                        {
                            a.pooled += usize::from(q.is_some_and(|q| q <= t));
                            a.mass_only += ok(Estimate::MassOnly, t);
                            a.rt_only += ok(Estimate::RtOnly, t);
                            a.wrong_peak_only += ok(Estimate::WrongPeak, t);
                        }
                        let supported = r.ms2.is_some_and(|m| m.supports());
                        if let Some(m) = &r.ms2 {
                            v.ms2_rows[0] += 1;
                            v.ms2_supported[0] += usize::from(m.supports());
                        }
                        // A transfer without a quantity (zero-amplitude fit)
                        // has nothing to report and is not accepted. The
                        // transfer FDR pool and its decoys are unchanged,
                        // so the estimator keeps its symmetry; only the
                        // targets it would accept shrink.
                        let acc = transfer_accepted(q, mbr_fdr, quantity.is_some());
                        let downgrade = acc
                            && q.is_some_and(|q| {
                                shared_feature.is_shared(&shared::Row {
                                    claim: &r.claim,
                                    shared: r.shared.as_deref(),
                                    ms2_supported: supported,
                                    transfer_q: q,
                                    llr: e.map_or(0.0, |e| e.llr),
                                    width: e.map_or(0.0, |e| e.w),
                                })
                            });
                        let t = match (acc, supported, downgrade) {
                            (true, _, true) => {
                                v.shared_feature += 1;
                                EvidenceType::SharedFeature
                            }
                            (true, true, false) => {
                                v.transfer_ms2 += 1;
                                EvidenceType::TransferMs2
                            }
                            (true, false, false) => EvidenceType::TransferMs1,
                            (false, _, _) => EvidenceType::BelowThreshold,
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
                        let wrong = match (e, ctx.identified.get(&(r.precursor, r.run))) {
                            (Some(p), Some(id)) => Some((p.t0 - id.t0).abs() > id.w),
                            _ => None,
                        };
                        if let Some(w) = wrong {
                            v.wrong_peak_pairs += 1;
                            v.wrong_peak += usize::from(w);
                        }
                        for (a, &t) in [&mut v.at_q01, &mut v.at_threshold]
                            .into_iter()
                            .zip(&thresholds)
                        {
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
                if r.role == Role::Transfer {
                    out.wrong.push((out.rows.len(), p_wrong[i]));
                }
                out.rows.push(OutRow {
                    precursor: r.precursor,
                    run: r.run,
                    evidence_type,
                    quantity: quantity.map(|q| q.0),
                    quantity_se: quantity.map(|q| q.1),
                    quantity_normalized,
                    transfer_q,
                    quant_q,
                    apex_rt_s: e.map(|e| e.t0),
                    library_q: r.library_q,
                    mass_err_ppm: e.map(|e| e.mass_err_ppm),
                    claim: if r.role == Role::Transfer {
                        r.claim
                    } else {
                        Claim::default()
                    },
                });
            }
            Ok(())
        },
        |chunk, (mut out, cache, at)| {
            if cache.is_some_and(|c| c.len() != at) {
                return Err(format!(
                    "{}: more cached transfer scores than rows",
                    source.name(chunk)
                ));
            }
            totals
                .lock()
                .expect("validation lock")
                .add_counts(&out.validation);
            if out.rows.is_empty() {
                return Ok(());
            }
            // quant_q of an accepted transfer: the running mean, over the
            // run's accepted transfers in ascending probability of a wrong
            // measurement, of that probability (every row is kept; a
            // transfer below the threshold keeps its transfer_q).
            let mut pool: Vec<(f64, usize)> = out
                .wrong
                .iter()
                .filter(|&&(i, _)| {
                    matches!(
                        out.rows[i].evidence_type,
                        EvidenceType::TransferMs2
                            | EvidenceType::TransferMs1
                            | EvidenceType::SharedFeature
                    )
                })
                .map(|&(i, p)| (p, i))
                .collect();
            running_mean_q(&mut pool, |i, q| out.rows[i].quant_q = Some(q));
            out.rows.sort_unstable_by_key(|r| (r.precursor, r.run));
            if let Some(w) = out
                .rows
                .windows(2)
                .find(|w| (w[0].precursor, w[0].run) == (w[1].precursor, w[1].run))
            {
                return Err(format!(
                    "{}: two output rows of precursor {} in run {}",
                    source.name(chunk),
                    w[0].precursor,
                    w[0].run
                ));
            }
            let path = dir.join(format!("rows-{chunk:06}.parquet"));
            let mut w = OutWriter::create(&path)?;
            w.write(&out.rows)?;
            w.finish()?;
            parts.lock().expect("parts lock").push((chunk, path));
            Ok(())
        },
    )?;
    let identity_entries = lib.entries.iter().filter(|e| e.identity_decoy).count();
    let mut v = MbrValidation {
        q_threshold: mbr_fdr,
        quant_feature_threshold: quant_feature_fdr,
        target_entries: lib.entries.len() - identity_entries,
        identity_entries,
        ..MbrValidation::default()
    };
    v.add_counts(&totals.into_inner().expect("validation lock"));
    let mut parts = parts.into_inner().expect("parts lock");
    parts.sort_unstable();
    let paths: Vec<PathBuf> = parts.into_iter().map(|(_, p)| p).collect();
    let t1 = Instant::now();
    log_summary(
        stems,
        &v,
        &summary,
        &fdr,
        &id_fdr,
        &norm,
        mbr_fdr,
        quant_feature_fdr,
    );
    summary.fdr = fdr.counts.clone();
    summary.identified_fdr = id_fdr.counts.clone();
    // Only the output tables are left to read.
    drop((all_summaries, t_models, i_models, norm, fdr, id_fdr));

    // The per-run tables merged into the outputs: one thread merges, one
    // writes `<base>.mbr.parquet`, and this one hands the rows on.
    summary.rows_written = write_merged(&paths, dir, lib, stems, params.mbr_path, on_features)?;
    write_mbr_validation(params.validation_path, &v)?;
    drop(tmp);
    eprintln!(
        "[mbr] finalize: output rows {:.1} s, merged outputs {:.1} s",
        (t1 - t0).as_secs_f64(),
        secs(t1)
    );
    summary.validation = v;
    Ok(summary)
}

/// Merges the per-run output tables `paths` into `<base>.mbr.parquet` at
/// `mbr_path` and hands the rows to `on_features`, both in (precursor, run)
/// order; the merge, the table writer and `on_features` run in three
/// threads. Returns the rows written.
fn write_merged(
    paths: &[PathBuf],
    dir: &Path,
    lib: &Library,
    stems: &[String],
    mbr_path: &Path,
    on_features: &mut dyn FnMut(Vec<MbrFeature>) -> Result<(), String>,
) -> Result<usize, String> {
    use std::sync::mpsc::sync_channel;
    let stopped = || "the output merge stopped early".to_string();
    let (to_table, table_rows) = sync_channel::<Arc<Vec<OutRow>>>(4);
    let (to_features, feature_rows) = sync_channel::<Arc<Vec<OutRow>>>(4);
    std::thread::scope(|scope| {
        let merger = scope.spawn(move || {
            merge_out(paths, dir, 0, &mut |rows| {
                let b = Arc::new(rows.to_vec());
                to_table.send(b.clone()).map_err(|_| stopped())?;
                to_features.send(b).map_err(|_| stopped())
            })
        });
        let table = scope.spawn(move || -> Result<(MbrParquetWriter, usize), String> {
            let mut writer = MbrParquetWriter::create(mbr_path)?;
            let mut n = 0;
            for rows in table_rows {
                for r in rows.iter() {
                    writer.push(&MbrRow {
                        precursor_key: &lib.meta[r.precursor as usize].key,
                        run: &stems[r.run as usize],
                        evidence_type: r.evidence_type,
                        quantity: r.quantity,
                        quantity_se: r.quantity_se,
                        quantity_normalized: r.quantity_normalized,
                        transfer_q: r.transfer_q,
                        quant_q: r.quant_q,
                        apex_rt_s: r.apex_rt_s,
                        library_q: r.library_q,
                        claim: r.evidence_type.is_transfer().then(|| MbrClaim {
                            claimed_by: r.claim.claimant.map(|p| lib.meta[p as usize].key.as_str()),
                            dt_s: r.claim.dt_s.map(f64::from),
                            identified: r.claim.identified,
                            co_apex: r.claim.co_apex,
                            n_co_apex: u32::from(r.claim.n_co_apex),
                            n_identified_co_apex: u32::from(r.claim.n_identified_co_apex),
                        }),
                    })?;
                }
                n += rows.len();
            }
            Ok((writer, n))
        });
        let mut handed = Ok(());
        for rows in feature_rows.iter() {
            handed = on_features(rows.iter().map(OutRow::feature).collect());
            if handed.is_err() {
                break;
            }
        }
        // Unblocks the merge when the rows were not all taken.
        drop(feature_rows);
        let merged = merger.join().expect("output merge thread");
        let written = table.join().expect("output table thread");
        // The first error by cause: a failed writer or hand-off stops the
        // merge, which then reports the stop. The table is only finished
        // (renamed into place) when every row reached it.
        let (writer, n) = written?;
        handed?;
        merged?;
        writer.finish()?;
        Ok(n)
    })
}

/// Logs what finalize did.
#[allow(clippy::too_many_arguments)]
fn log_summary(
    stems: &[String],
    v: &MbrValidation,
    summary: &FinalizeSummary,
    fdr: &QValues,
    id_fdr: &QValues,
    norm: &quant::normalize::RunNormalization,
    mbr_fdr: f64,
    quant_feature_fdr: Option<f64>,
) {
    eprintln!(
        "[mbr] transfer FDR (internal quant-quality score) features: {}",
        TRANSFER_FEATURES.join(",")
    );
    eprintln!(
        "[mbr] library entries {} targets, {} identity decoys; {} transfer targets not transferred (best donor PSM q > {TRANSFER_DONOR_Q})",
        v.target_entries,
        v.identity_entries,
        v.donor_gated
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
    for f in &id_fdr.counts.failed_folds {
        eprintln!("[mbr] identified quality fold without a model: {f}");
    }
    eprintln!(
        "[mbr] identified rows {}: {} with a quality q-value (quant_q, internal), {} at q <= {VALIDATION_Q}; {} identified rows without features, {} unscored, {} scored without a q-value{}",
        v.identified,
        v.identified_scored,
        v.identified_q01,
        summary.identified_no_features,
        id_fdr.counts.unscored,
        id_fdr.counts.no_q,
        quant_feature_fdr.map_or(String::new(), |t| format!(
            "; {} identified_below_threshold at --quant-feature-fdr {t}",
            v.identified_below_threshold
        ))
    );
    eprintln!(
        "[mbr] pseudo-transfers {}: wrong peak {} of {}",
        v.pseudo_transfers, v.wrong_peak, v.wrong_peak_pairs
    );
    let rate = |a: usize, b: usize| {
        if b > 0 {
            format!("{:.4}", a as f64 / b as f64)
        } else {
            "NA".to_string()
        }
    };
    eprintln!(
        "[mbr] MS2 support (at least {} fragments at angle >= {}): transfer targets {} of {} ({}), mass-shift decoys {} of {} ({}), RT-shift decoys {} of {} ({}), pseudo-transfers {} of {} ({}); {} accepted transfers are transfer_ms2",
        quant::ms2::SUPPORT_MIN_MATCHED,
        quant::ms2::SUPPORT_MIN_ANGLE,
        v.ms2_supported[0],
        v.ms2_rows[0],
        rate(v.ms2_supported[0], v.ms2_rows[0]),
        v.ms2_supported[1],
        v.ms2_rows[1],
        rate(v.ms2_supported[1], v.ms2_rows[1]),
        v.ms2_supported[2],
        v.ms2_rows[2],
        rate(v.ms2_supported[2], v.ms2_rows[2]),
        v.ms2_pseudo_supported,
        v.ms2_pseudo_rows,
        rate(v.ms2_pseudo_supported, v.ms2_pseudo_rows),
        v.transfer_ms2
    );
    eprintln!(
        "[mbr] shared-feature classifier (probability above {}): {} accepted transfers reported as shared_feature",
        quant::shared::Model::shipped().threshold,
        v.shared_feature
    );
    eprintln!(
        "[mbr] run normalisation (quantity_normalized): {} reference precursors (identified, quant_q <= {VALIDATION_Q}, in >= 3 runs or all runs)",
        norm.precursors
    );
    for (run, stem) in stems.iter().enumerate() {
        match norm.summary(run as u16) {
            Some((level, sd, n)) => eprintln!(
                "[mbr]   {stem}: median offset {level:+.3} log2, sd along RT {sd:.3}, {n} reference cells"
            ),
            None => eprintln!("[mbr]   {stem}: not normalised (too few reference cells)"),
        }
    }
    let shown: &[(f64, &Accepted)] = if mbr_fdr == VALIDATION_Q {
        &[(VALIDATION_Q, &v.at_q01)]
    } else {
        &[(VALIDATION_Q, &v.at_q01), (mbr_fdr, &v.at_threshold)]
    };
    for &(t, a) in shown {
        eprintln!(
            "[mbr] at q <= {t}: decoy transfers accepted: mass-shift {} of {}, RT-shift {} of {}, identity decoys {} of {}; {} transfers would pass on the wrong-peak decoys alone",
            a.decoys[0],
            v.decoy_transfers[0],
            a.decoys[1],
            v.decoy_transfers[1],
            a.decoys[2],
            v.decoy_transfers[2],
            a.wrong_peak_only
        );
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
}

/// Whether a transfer with q-value `q` is accepted at `mbr_fdr`: it passes
/// the threshold and has a quantity (a fit with zero amplitude has none).
fn transfer_accepted(q: Option<f64>, mbr_fdr: f64, has_quantity: bool) -> bool {
    has_quantity && q.is_some_and(|q| q <= mbr_fdr)
}

#[cfg(test)]
mod tests {
    /// A transfer whose fit has zero amplitude reports no quantity and is
    /// not accepted whatever its q-value.
    #[test]
    fn zero_amplitude_transfer_is_not_accepted() {
        use super::transfer_accepted;
        assert!(transfer_accepted(Some(0.001), 0.01, true));
        assert!(!transfer_accepted(Some(0.001), 0.01, false));
        assert!(!transfer_accepted(Some(0.5), 0.01, true));
        assert!(!transfer_accepted(None, 0.01, true));
    }

    #[test]
    fn the_brightest_share_of_precursors_reports_trace_areas() {
        use super::{mark_trace_precursors, PrecursorSummary, Summaries};
        let summaries = |snrs: &[Option<f64>]| -> Summaries {
            snrs.iter()
                .enumerate()
                .map(|(p, &snr)| {
                    (
                        p as u32,
                        PrecursorSummary {
                            snr,
                            ..PrecursorSummary::default()
                        },
                    )
                })
                .collect()
        };
        let marked = |s: &Summaries| {
            let mut v: Vec<u32> = s.iter().filter(|(_, s)| s.trace).map(|(p, _)| *p).collect();
            v.sort_unstable();
            v
        };
        // Signal-to-noise 1..=10 for precursors 0..10, none for precursor 10.
        let mut snrs: Vec<Option<f64>> = (1..=10).map(|v| Some(v as f64)).collect();
        snrs.push(None);
        let mut s = summaries(&snrs);
        mark_trace_precursors(&mut s, 0.3);
        assert_eq!(marked(&s), vec![7, 8, 9]);
        // ceil(0.25 * 10) = 3 as well; a tie at the cut is kept.
        snrs[6] = Some(8.0);
        let mut s = summaries(&snrs);
        mark_trace_precursors(&mut s, 0.25);
        assert_eq!(marked(&s), vec![6, 7, 8, 9]);
        let mut s = summaries(&snrs);
        mark_trace_precursors(&mut s, 1.0);
        assert_eq!(marked(&s), (0..10).collect::<Vec<u32>>());
        // A small share still marks the brightest one.
        let mut s = summaries(&snrs);
        mark_trace_precursors(&mut s, 1e-9);
        assert_eq!(marked(&s), vec![9]);
        let mut s = summaries(&[None, None]);
        mark_trace_precursors(&mut s, 1.0);
        assert!(marked(&s).is_empty());
    }

    /// Two clusters of runs, entries identified in all runs of one cluster:
    /// the strata are cut from the prior of the pairs without an
    /// identification, with ascending edges.
    #[test]
    fn strata_edges_are_ascending_and_bounded() {
        use quant::library::{DonorPsm, LibraryEntry};
        let psm = DonorPsm {
            score: 1.0,
            q: 0.001,
            pep: 0.001,
        };
        let entries: Vec<LibraryEntry> = (0..400u32)
            .map(|p| LibraryEntry {
                precursor: p,
                identity_decoy: false,
                q: 0.001,
                donor_psm: psm,
                consensus_rt: 10.0,
                faims_cv: None,
                ion_mobility: None,
                // Cluster 0 (runs 0-5) or 1 (runs 6-11); a few misses.
                donors: (0..6u16)
                    .map(|r| r + 6 * (p % 2) as u16)
                    .filter(|r| (u32::from(*r) + p) % 7 != 0)
                    .collect(),
                best_rt: Vec::new(),
            })
            .collect();
        let entry_of: Vec<Option<&LibraryEntry>> = entries.iter().map(Some).collect();
        let neighbours: Vec<Vec<u16>> = (0..12u16)
            .map(|r| {
                let c = r / 6;
                let mut n: Vec<u16> = (6 * c..6 * c + 6).filter(|&x| x != r).collect();
                n.extend((0..12u16).filter(|x| x / 6 != c));
                n
            })
            .collect();
        let edges = strata_edges(&entry_of, &neighbours, 12, 5);
        assert!(!edges.is_empty());
        assert!(edges.iter().all(|&e| (0.0..1.0).contains(&e)), "{edges:?}");
        assert!(edges.windows(2).all(|w| w[0] < w[1]), "{edges:?}");
        assert!(edges.len() < MAX_STRATA);
    }

    /// The transfer gate reads the library's best-donor q-value, which is the
    /// pooled q-value of that PSM when the rescoring result has one: an entry
    /// whose donor PSM is good in its run but not pooled is gated, and the
    /// other way round.
    #[test]
    fn the_transfer_gate_reads_the_pooled_q_of_the_best_donor() {
        use quant::align::AlignParams;
        use quant::library::{build_library, PrecursorMeta, PsmIn};
        let meta: Vec<PrecursorMeta> = (0..200)
            .map(|i| PrecursorMeta {
                key: format!("P{i}/2"),
                mz: 500.0 + i as f64,
                charge: 2,
                envelope: [0.5, 0.3, 0.15, 0.05],
            })
            .collect();
        let psms = |q: f64, q_gate: f64| -> Vec<PsmIn> {
            (0..200u32)
                .flat_map(|p| {
                    (0..2u16).map(move |r| PsmIn {
                        run: r,
                        precursor: p,
                        rt: 10.0 + p as f64,
                        score: 100.0 + p as f64,
                        q,
                        q_gate,
                        pep: 0.001,
                        is_decoy: false,
                        faims_cv: None,
                        ion_mobility: None,
                    })
                })
                .collect()
        };
        let ap = AlignParams::default();
        // Per-run q passes the gate, pooled q does not.
        let lib = build_library(&psms(0.0005, 0.005), meta.clone(), 2, 1.0, 0.01, None, &ap);
        let idx = EntryIndex::new(&lib, 2);
        assert!(idx.gated(0));
        // Per-run q fails it, pooled q passes (the PSM is still a donor).
        let lib = build_library(&psms(0.005, 0.0005), meta.clone(), 2, 1.0, 0.01, None, &ap);
        let idx = EntryIndex::new(&lib, 2);
        assert!(!idx.gated(0));
        // Without a pooled q-value the gate reads the q-value.
        let lib = build_library(&psms(0.005, f64::NAN), meta, 2, 1.0, 0.01, None, &ap);
        assert!(EntryIndex::new(&lib, 2).gated(0));
    }

    /// An even histogram fills every one of the `MAX_STRATA` strata, and the
    /// last one survives packing into the family byte.
    #[test]
    fn an_even_prior_histogram_gives_max_strata() {
        let edges = edges_of_hist(&[1; 32]);
        assert_eq!(edges.len() + 1, MAX_STRATA, "{edges:?}");
        assert!(edges.windows(2).all(|w| w[0] < w[1]), "{edges:?}");
        // Packed and unpacked, the last stratum is still distinct.
        let last = stratum_of(1.0, &edges);
        assert_eq!(usize::from(last), MAX_STRATA - 1);
        let packed = crate::transfer_fdr::unpack_family(pack_family(TARGET, last));
        assert_eq!(packed, (TARGET, MAX_STRATA - 1));
        assert!(edges_of_hist(&[]).is_empty());
        assert!(edges_of_hist(&[0, 0, 0]).is_empty());
    }

    use super::*;

    struct Lcg(u64);
    impl Lcg {
        fn below(&mut self, n: u64) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) % n
        }
    }

    #[test]
    fn median_window_gives_every_median_and_leave_one_out_median() {
        let mut g = Lcg(3);
        for n in 0..40usize {
            for _ in 0..20 {
                // Coarse values with ties and both zeros.
                let v: Vec<f64> = (0..n)
                    .map(|_| match g.below(9) {
                        0 => 0.0,
                        1 => -0.0,
                        k => (k as f64 - 4.0) * 0.25,
                    })
                    .collect();
                let mut sorted = v.clone();
                sorted.sort_by(f64::total_cmp);
                let w = MedianWindow::new(&sorted);
                let full = (!v.is_empty()).then(|| median(v.clone()));
                assert_eq!(
                    w.median().map(f64::to_bits),
                    full.map(f64::to_bits),
                    "{v:?}"
                );
                for i in 0..n {
                    let mut rest = v.clone();
                    let x = rest.remove(i);
                    let want = (!rest.is_empty()).then(|| median(rest));
                    assert_eq!(
                        w.median_without(x).map(f64::to_bits),
                        want.map(f64::to_bits),
                        "{v:?} without {x}"
                    );
                }
            }
        }
    }

    fn out_row(precursor: u32, run: u16) -> OutRow {
        OutRow {
            precursor,
            run,
            evidence_type: EVIDENCE_CODES[(precursor as usize + run as usize) % 5],
            quantity: Some(precursor as f64 * 1.5),
            quantity_se: (run % 2 == 0).then_some(f64::NAN),
            quantity_normalized: None,
            transfer_q: Some(0.25),
            quant_q: Some(-0.0),
            apex_rt_s: Some(f64::INFINITY),
            library_q: 1e-3,
            mass_err_ppm: Some(-1.25),
            claim: Claim {
                claimant: (precursor % 2 == 0).then_some(precursor + 1),
                dt_s: (precursor % 2 == 0).then_some(-2.5),
                identified: precursor % 4 == 0,
                co_apex: precursor % 2 == 0,
                n: 2,
                n_co_apex: 1,
                n_identified_co_apex: u16::from(precursor % 4 == 0),
                ..Claim::default()
            },
        }
    }

    #[test]
    fn output_tables_merge_in_precursor_and_run_order() {
        let dir = tempfile::tempdir().unwrap();
        let runs = MAX_OPEN_PARTS as u16 + 9;
        let mut paths = Vec::new();
        let mut want = Vec::new();
        for run in 0..runs {
            let rows: Vec<OutRow> = (0..7u32)
                .filter(|p| (p + run as u32) % 3 != 0)
                .map(|p| out_row(p * 2 + u32::from(run % 2), run))
                .collect();
            want.extend(rows.iter().copied());
            let path = dir.path().join(format!("rows-{run:06}.parquet"));
            let mut w = OutWriter::create(&path).unwrap();
            w.write(&rows).unwrap();
            w.finish().unwrap();
            paths.push(path);
        }
        want.sort_by_key(|r| (r.precursor, r.run));
        let mut got = Vec::new();
        merge_out(&paths, dir.path(), 0, &mut |rows| {
            got.extend_from_slice(rows);
            Ok(())
        })
        .unwrap();
        assert_eq!(got.len(), want.len());
        for (a, b) in got.iter().zip(&want) {
            assert_eq!(
                (a.precursor, a.run, a.evidence_type),
                (b.precursor, b.run, b.evidence_type)
            );
            let bits = |r: &OutRow| {
                [
                    r.quantity,
                    r.quantity_se,
                    r.quantity_normalized,
                    r.transfer_q,
                    r.quant_q,
                    r.apex_rt_s,
                    Some(r.library_q),
                    r.mass_err_ppm,
                ]
                .map(|v| v.map(f64::to_bits))
            };
            assert_eq!(bits(a), bits(b));
            assert_eq!(a.claim, b.claim);
        }
    }

    #[test]
    fn spilled_scores_give_the_in_memory_q_values() {
        spilled_scores_match(1);
        spilled_scores_match(3);
    }

    /// Spilled and in-memory scores give one set of q-values, with `strata`
    /// prior strata.
    fn spilled_scores_match(strata: usize) {
        let dir = tempfile::tempdir().unwrap();
        let mut g = Lcg(17);
        let rows: Vec<(u32, u8, Option<f32>)> = (0..5000u32)
            .map(|i| {
                let fam = g.below(4) as u8 | ((g.below(strata as u64) as u8) << 4);
                let s = (g.below(41) as f32 - 20.0) / 4.0 + if fam & 0x0F == 0 { 1.0 } else { 0.0 };
                (
                    i % 700,
                    fam,
                    (i % 97 != 0).then_some(if s == 0.0 && i % 2 == 0 { -0.0 } else { s }),
                )
            })
            .collect();
        let mut one = FoldScores::with_strata(5, strata);
        for &(p, f, s) in &rows {
            one.add(p, f, s);
        }
        let store = ScoreStore::new("t", dir.path(), 300, 5, strata);
        for piece in rows.chunks(211) {
            let mut s = FoldScores::with_strata(5, strata);
            for &(p, f, sc) in piece {
                s.add(p, f, sc);
            }
            store.add(s).unwrap();
        }
        let a = one.finish(&[]);
        let (b, spills) = store.finish(&[], &[0.05]).unwrap();
        assert!(spills > 3, "{spills}");
        assert_eq!(a.counts, b.counts);
        assert_eq!(a.pooled_steps(), b.pooled_steps());
        for &(p, f, s) in &rows {
            for e in [
                Estimate::Pooled,
                Estimate::WrongPeak,
                Estimate::MassOnly,
                Estimate::RtOnly,
            ] {
                assert_eq!(a.accepts(e, p, f, s, 0.05), b.accepts(e, p, f, s, 0.05));
            }
            assert_eq!(
                a.q(Estimate::Pooled, p, f, s).map(f64::to_bits),
                b.q(Estimate::Pooled, p, f, s).map(f64::to_bits)
            );
            assert_eq!(
                a.pseudo_q(p, f, s).map(f64::to_bits),
                b.pseudo_q(p, f, s).map(f64::to_bits)
            );
        }
        let left = std::fs::read_dir(dir.path()).unwrap().count();
        assert_eq!(left, 0, "spilled score tables are removed");
    }

    #[test]
    fn wrong_calibration_is_a_non_increasing_step_function_of_the_margin() {
        let mut g = Lcg(5);
        let pairs: Vec<(f32, bool)> = (0..4000)
            .map(|_| {
                let s = g.below(2000) as f32 / 100.0 - 10.0;
                let p = 1.0 / (1.0 + f64::from(s).exp());
                (s, (g.below(1000) as f64) / 1000.0 < p)
            })
            .collect();
        let c = WrongCalibration::fit(pairs).unwrap();
        assert!(c.p.windows(2).all(|w| w[0] >= w[1]));
        assert!(c.lo.windows(2).all(|w| w[0] < w[1]));
        assert!(c.at(-100.0) >= c.at(0.0) && c.at(0.0) >= c.at(100.0));
        assert!(c.at(100.0) < 0.05 && c.at(-100.0) > 0.8);
        assert!(WrongCalibration::fit(Vec::new()).is_none());
        // Tied margins are one point: the share of wrong ones.
        let c = WrongCalibration::fit(vec![(1.0, true), (1.0, false), (1.0, false), (1.0, false)])
            .unwrap();
        assert_eq!(c.p, vec![0.25]);
    }

    fn quant_qs(rows: &[f64]) -> Vec<f64> {
        let mut pool: Vec<(f64, usize)> = rows.iter().copied().zip(0..).collect();
        let mut q = vec![f64::NAN; rows.len()];
        running_mean_q(&mut pool, |i, v| q[i] = v);
        q
    }

    #[test]
    fn tied_wrong_peak_probabilities_share_one_quant_q() {
        // Groups 0.1 x1, 0.2 x3, 0.5 x2: each group gets the mean over all
        // rows up to the end of the group.
        let rows = [0.2, 0.5, 0.1, 0.2, 0.5, 0.2];
        let q = quant_qs(&rows);
        let at = |p: f64| -> Vec<f64> {
            rows.iter()
                .zip(&q)
                .filter(|r| *r.0 == p)
                .map(|r| *r.1)
                .collect()
        };
        assert_eq!(at(0.1), vec![0.1]);
        let tie = (0.1 + 3.0 * 0.2) / 4.0;
        assert!(at(0.2).iter().all(|&v| (v - tie).abs() < 1e-12), "{q:?}");
        let all = (0.1 + 3.0 * 0.2 + 2.0 * 0.5) / 6.0;
        assert!(at(0.5).iter().all(|&v| (v - all).abs() < 1e-12), "{q:?}");
        // Any permutation of the rows gives every row the same quant_q.
        let mut g = Lcg(11);
        let probs: Vec<f64> = (0..500)
            .map(|_| f64::from(g.below(7) as u32) / 10.0)
            .collect();
        let q = quant_qs(&probs);
        for _ in 0..5 {
            let mut order: Vec<usize> = (0..probs.len()).collect();
            for i in (1..order.len()).rev() {
                order.swap(i, g.below(i as u64 + 1) as usize);
            }
            let permuted: Vec<f64> = order.iter().map(|&i| probs[i]).collect();
            let qp = quant_qs(&permuted);
            for (j, &i) in order.iter().enumerate() {
                assert_eq!(qp[j].to_bits(), q[i].to_bits());
            }
        }
    }
}
