//! `quant finalize`: score every transfer of the study, estimate its q-values, normalise the runs
//! and write the MBR table. The stream and its summary live here; each step has its own module.

mod extract;
mod index;
mod out;
mod scores;
mod shapes;
mod wrong;

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use output::mbr::{write_mbr_validation, Accepted, EvidenceType, MbrValidation, VALIDATION_Q};
use quant::align::RunMap;
use quant::claims::Claim;
use quant::decoys::Hypothesis;
use quant::library::Library;
use quant::normalize::{fit as fit_normalization, NormCell, NormParams, RunNormalization};
use quant::shared;

use crate::quant::mbr::{
    describe_width_weights, donor_fit, family, relative_fit_width, reported_quantity,
    row_width_weight, weigh_widths, whole_pattern, width_deviation, ExtractRow, Role,
};
use crate::transfers::{
    member_seed, row_unit, sample_hash, train_bagged_folds, train_folds, transfer_train_params,
    Estimate, FdrCounts, FdrRow, FoldModels, FoldScores, QValues, TrainingSample,
    TRANSFER_FEATURES,
};
#[cfg(test)]
pub(crate) use extract::ExtractRuns;
use extract::{for_each_block, Bounds};
pub(crate) use extract::{BlockContext, ExtractFiles, ExtractSource};
pub(crate) use index::{for_each_scored, table_features, EntryIndex, Table};
use out::{write_merged, OutRow, OutWriter};
use scores::{
    cached_scores_path, is_cached, read_cached_scores, write_scores, ScoreStore, NO_SCORE,
};
use shapes::{
    identified_snr, mark_trace_precursors, quantity_rule, summarise, BlockShapes, Summaries,
};
pub(crate) use wrong::{
    mismatch, pseudo_example, running_mean_q, transfer_quant_q, WrongCalibration,
};
#[cfg(test)]
pub(crate) use wrong::{wrong_features, wrong_label, Mismatch};
use wrong::{wrong_peak_probabilities, WRONG_CALIBRATION_ROWS};

#[derive(Debug, Clone, Default)]
pub(crate) struct FinalizeSummary {
    pub validation: MbrValidation,
    pub rows_written: usize,
    pub no_features: usize,
    pub fdr: FdrCounts,
    pub identified_no_features: usize,
    pub identified_fdr: FdrCounts,
    #[cfg(test)]
    pub features: Vec<MbrFeature>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MbrFeature {
    pub precursor: u32,
    pub run: u16,
    pub evidence_type: EvidenceType,
    pub quantity: Option<f64>,
    pub quantity_se: Option<f64>,
    pub quantity_normalized: Option<f64>,
    pub transfer_q: Option<f64>,
    pub quant_q: Option<f64>,
    pub apex_rt_s: Option<f64>,
    pub library_q: f64,
    pub mass_err_ppm: Option<f64>,
    pub claim: Claim,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FinalizeLimits {
    pub max_train_rows: usize,
    pub transfer_models: usize,
    pub max_norm_cells: usize,
    pub max_donor_fits: usize,
    pub max_held_scores: usize,
    pub temp_dir: Option<PathBuf>,
}

pub(crate) const DEFAULT_MAX_TRAIN_ROWS: usize = 3_000_000;
pub(crate) const DEFAULT_TRACE_FRACTION: f64 = 0.1;
pub(crate) const DEFAULT_TRANSFER_MODELS: usize = 8;
pub(crate) const DEFAULT_MAX_NORM_CELLS: usize = 20_000_000;
pub(crate) const DEFAULT_MAX_DONOR_FITS: usize = 12_000_000;
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

pub(crate) const SCORE_CHUNK: usize = 65_536;
const READ_BATCH_ROWS: usize = 65_536;
const SCORE_GROUP_ROWS: usize = 65_536;
const OUT_GROUP_ROWS: usize = 4_096;
const MAX_OPEN_PARTS: usize = 512;

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

pub(crate) struct FinalizeParams<'a> {
    pub mbr_fdr: f64,
    pub quant_feature_fdr: Option<f64>,
    pub seed: u64,
    pub limits: FinalizeLimits,
    pub mbr_path: &'a Path,
    pub validation_path: &'a Path,
    pub trace_fraction: f64,
}

#[derive(Default)]
struct RunOutput {
    rows: Vec<OutRow>,
    wrong: Vec<(usize, f64)>,
    validation: MbrValidation,
}

fn secs(t: Instant) -> f64 {
    t.elapsed().as_secs_f64()
}

type Key = (u32, u64);

fn check_thresholds(params: &FinalizeParams<'_>) -> Result<(), String> {
    let (mbr_fdr, quant_feature_fdr, trace_fraction) = (
        params.mbr_fdr,
        params.quant_feature_fdr,
        params.trace_fraction,
    );
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
    if !(0.0..=1.0).contains(&trace_fraction) {
        return Err(format!(
            "trace-area precursor share {trace_fraction} is outside [0, 1]"
        ));
    }
    Ok(())
}

fn donor_summaries(
    lib: &Library,
    source: &dyn ExtractSource,
    bounds: Bounds,
    trace_fraction: f64,
    limits: &FinalizeLimits,
    dir: &Path,
) -> Result<Summaries, String> {
    let maps = &lib.alignment.maps;
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
                    let k = relative_fit_width(r);
                    let f = donor_fit(r, maps);
                    let snr = (trace_fraction > 0.0).then(|| identified_snr(r)).flatten();
                    if k.is_some() || f.is_some() || snr.is_some() {
                        local.push((r.precursor, k, f, snr));
                    }
                }
                Ok(())
            },
            |_, local| {
                items.lock().expect("summary lock").push(local);
                Ok(())
            },
        )?;
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
    log::debug!(
        target: "finalize",
        "{} extract tables, {} threads, temporary tables in {}; donor summaries of {} precursors in {groups} pass(es), {:.1} s",
        source.chunks(),
        rayon::current_num_threads(),
        dir.display(),
        all_summaries.len(),
        secs(t0)
    );
    Ok(all_summaries)
}

#[allow(clippy::too_many_arguments)]
fn normalisation_cells(
    source: &dyn ExtractSource,
    bounds: Bounds,
    idx: &EntryIndex<'_>,
    dir: &Path,
    summaries: &Summaries,
    id_fdr: &QValues,
    limits: &FinalizeLimits,
    seed: u64,
) -> Result<TrainingSample<Key, (NormCell, Option<f64>)>, String> {
    let cells: Mutex<TrainingSample<Key, (NormCell, Option<f64>)>> =
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
                        .q(Estimate::Pooled, r.precursor, fam.into(), s)
                        .is_some_and(|q| q <= VALIDATION_Q)
                {
                    continue;
                }
                let Some(e) = r.evidence.as_ref().filter(|e| e.amplitude > 0.0) else {
                    continue;
                };
                let (k, trace) = quantity_rule(summaries, r.precursor);
                let q = reported_quantity(e, k, trace, 0.0).0;
                if q > 0.0 && q.is_finite() {
                    local.push((
                        (chunk as u32, first + i as u64),
                        u64::from(r.precursor),
                        (
                            NormCell {
                                precursor: r.precursor,
                                run: r.run,
                                log2_quantity: q.log2(),
                                rt: e.t0,
                            },
                            width_deviation(e, k, trace),
                        ),
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
        log::debug!(
            target: "mbr",
            "run normalisation: fitted on {kept} of {seen} reference cells (whole precursors, {share} of them)"
        );
    }
    Ok(cells)
}

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
    check_thresholds(params)?;
    let trace_fraction = params.trace_fraction;
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

    let all_summaries = donor_summaries(lib, source, bounds, trace_fraction, limits, dir)?;
    let summaries = &all_summaries;

    let inputs = StreamInputs {
        source,
        bounds,
        idx: &idx,
        maps,
        summaries,
        limits,
        seed,
    };
    let (t_models, i_models, w_models, w_seen) = train_transfer_models(&inputs, &mut summary)?;

    let t0 = Instant::now();
    let w_stride = (w_seen / WRONG_CALIBRATION_ROWS + 1) as u64;
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
            write_scores(&cached_scores_path(dir, chunk), &cache)?;
            t_store.add(ts)?;
            i_store.add(is)
        },
    )?;
    let wrong_cal =
        WrongCalibration::fit(w_held.into_inner().expect("wrong-peak calibration lock"));
    match &wrong_cal {
        Some(c) => log::debug!(
            target: "mbr",
            "wrong-peak calibration: {} blocks, probability {:.4} at the best margin, {:.4} at the worst",
            c.p.len(),
            c.p.last().copied().unwrap_or(f64::NAN),
            c.p.first().copied().unwrap_or(f64::NAN)
        ),
        None => {
            log::warn!(target: "mbr", "wrong-peak calibration: no pseudo-transfers; every transfer's quant_q is 1")
        }
    }
    let thresholds = [VALIDATION_Q, mbr_fdr];
    let t1 = Instant::now();
    let (t_done, i_done) = rayon::join(
        || t_store.finish(&t_models.failed_folds, &thresholds),
        || i_store.finish(&i_models.failed_folds, &[]),
    );
    let ((fdr, t_spills), (id_fdr, i_spills)) = (t_done?, i_done?);
    log::debug!(
        target: "finalize",
        "held-out scores {:.1} s, step functions {:.1} s ({} and {} steps, {} and {} pooled; {t_spills} and {i_spills} score spills)",
        (t1 - t0).as_secs_f64(),
        secs(t1),
        fdr.steps(),
        id_fdr.steps(),
        fdr.pooled_steps(),
        id_fdr.pooled_steps()
    );

    let t0 = Instant::now();
    let cells = normalisation_cells(source, bounds, &idx, dir, summaries, &id_fdr, limits, seed)?;
    let t1 = Instant::now();
    let (cells, width_weights) = weigh_widths(cells.into_items(), stems.len());
    log::debug!(
        target: "finalize",
        "run width weights (0 = relative width, 1 = own fitted width), runs per weight: {}",
        describe_width_weights(&width_weights)
    );
    let norm = fit_normalization(&cells, stems.len(), &NormParams::default());
    log::debug!(
        target: "finalize",
        "normalisation cells {:.1} s, fit {:.1} s",
        (t1 - t0).as_secs_f64(),
        secs(t1)
    );

    let t0 = Instant::now();
    let shared_feature = shared::Model::shipped();
    let scoring = OutputScoring {
        lib,
        dir,
        shared_feature,
        t_models: &t_models,
        i_models: &i_models,
        w_models: &w_models,
        wrong_cal: wrong_cal.as_ref(),
        fdr: &fdr,
        id_fdr: &id_fdr,
        norm: &norm,
        width_weights: &width_weights,
        thresholds,
        mbr_fdr,
        quant_feature_fdr,
    };
    let (totals, paths) = write_output_rows(&inputs, &scoring)?;
    let identity_entries = lib.entries.iter().filter(|e| e.identity_decoy).count();
    let mut v = MbrValidation {
        q_threshold: mbr_fdr,
        quant_feature_threshold: quant_feature_fdr,
        target_entries: lib.entries.len() - identity_entries,
        identity_entries,
        ..MbrValidation::default()
    };
    v.add_counts(&totals);
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
    drop((all_summaries, t_models, i_models, norm, fdr, id_fdr));

    summary.rows_written = write_merged(&paths, dir, lib, stems, params.mbr_path, on_features)?;
    write_mbr_validation(params.validation_path, &v)?;
    drop(tmp);
    log::debug!(
        target: "finalize",
        "output rows {:.1} s, merged outputs {:.1} s",
        (t1 - t0).as_secs_f64(),
        secs(t1)
    );
    summary.validation = v;
    Ok(summary)
}

/// The per-run output rows of one extract block in the output pass, with the partial checks
/// and the measurement-q pool.
type OutputState = (RunOutput, Option<Vec<f32>>, usize);

/// What the output pass reads besides the block inputs: the library and temp folder, the
/// shared-feature model, the fold models, the wrong-peak calibration, the two q-value step
/// functions, the run normalisation and width weights, and the thresholds.
#[derive(Clone, Copy)]
struct OutputScoring<'a> {
    lib: &'a Library,
    dir: &'a Path,
    shared_feature: &'a shared::Model,
    t_models: &'a FoldModels,
    i_models: &'a FoldModels,
    w_models: &'a FoldModels,
    wrong_cal: Option<&'a WrongCalibration>,
    fdr: &'a QValues,
    id_fdr: &'a QValues,
    norm: &'a RunNormalization,
    width_weights: &'a [f64],
    thresholds: [f64; 2],
    mbr_fdr: f64,
    quant_feature_fdr: Option<f64>,
}

/// The output pass: every block's rows are scored, typed and written to a sorted part file;
/// returns the validation counts summed over blocks and the part paths in block order.
fn write_output_rows(
    inputs: &StreamInputs<'_>,
    scoring: &OutputScoring<'_>,
) -> Result<(MbrValidation, Vec<PathBuf>), String> {
    let StreamInputs { source, bounds, .. } = *inputs;
    let dir = scoring.dir;
    let parts: Mutex<Vec<(usize, PathBuf)>> = Mutex::new(Vec::new());
    let totals: Mutex<MbrValidation> = Mutex::new(MbrValidation::default());
    for_each_block(
        source,
        false,
        bounds,
        |_| (RunOutput::default(), None, 0usize),
        |state, chunk, _, rows| output_block(inputs, scoring, state, chunk, rows),
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
            running_mean_q(&mut pool, |i, q| {
                let r = &mut out.rows[i];
                r.measurement_q = Some(q);
                r.quant_q = Some(transfer_quant_q(r.transfer_q, q));
            });
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
    let mut parts = parts.into_inner().expect("parts lock");
    parts.sort_unstable();
    let paths: Vec<PathBuf> = parts.into_iter().map(|(_, p)| p).collect();
    Ok((totals.into_inner().expect("validation lock"), paths))
}

/// One block of the output pass: reads the block's cached scores, scores the pseudo-transfers,
/// and appends a typed output row per target row to `state`.
fn output_block(
    inputs: &StreamInputs<'_>,
    scoring: &OutputScoring<'_>,
    state: &mut OutputState,
    chunk: usize,
    rows: &[ExtractRow],
) -> Result<(), String> {
    let StreamInputs {
        source,
        idx,
        maps,
        summaries,
        ..
    } = *inputs;
    let OutputScoring {
        lib,
        dir,
        shared_feature,
        t_models,
        i_models,
        w_models,
        wrong_cal,
        fdr,
        id_fdr,
        norm,
        width_weights,
        thresholds,
        mbr_fdr,
        quant_feature_fdr,
    } = *scoring;
    let (out, cache, at) = state;
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
    let mut scored: Vec<Option<(u32, Option<f32>)>> = vec![None; rows.len()];
    for (i, r) in rows.iter().enumerate() {
        let Some((t @ (Table::Transfer | Table::Identified), fam)) = idx.table(r) else {
            continue;
        };
        let s = *cache
            .get(*at)
            .ok_or_else(|| format!("{}: fewer cached scores than rows", source.name(chunk)))?;
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
        idx,
        &shapes,
        maps,
        t_models,
        i_models,
        |t: Table| t == Table::Pseudo,
        |i, _, fam, s| {
            scored[i] = Some((fam, s));
        },
    );
    let q_of = |q: &QValues, e: Estimate, i: usize| {
        scored[i].and_then(|(fam, s)| q.q(e, rows[i].precursor, fam, s))
    };
    let reported = |r: &ExtractRow| {
        let (k, trace) = quantity_rule(summaries, r.precursor);
        r.evidence.as_ref().filter(|e| e.amplitude > 0.0).map(|e| {
            whole_pattern(
                reported_quantity(e, k, trace, row_width_weight(r, width_weights)),
                lib.meta[r.precursor as usize].isotope_share,
            )
        })
    };
    let p_wrong = wrong_peak_probabilities(rows, idx, w_models, wrong_cal);
    for (i, r) in rows.iter().enumerate() {
        if idx.is_identity(r.precursor) {
            if (r.role, r.kind) == (Role::Transfer, Hypothesis::Target) {
                count_decoy(v, q_of(fdr, Estimate::Pooled, i), 2);
            }
            continue;
        }
        if r.kind != Hypothesis::Target {
            if r.role == Role::Transfer {
                let f = family(r.kind) as usize;
                count_decoy(v, q_of(fdr, Estimate::Pooled, i), f - 1);
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
                let q = q_of(id_fdr, Estimate::Pooled, i);
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
                let q = q_of(fdr, Estimate::Pooled, i);
                v.transfers += 1;
                v.transfers_scored += usize::from(q.is_some());
                let ok = |e: Estimate, t: f64| {
                    usize::from(
                        scored[i].is_some_and(|(fam, s)| fdr.accepts(e, r.precursor, fam, s, t)),
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
                let wrong = e
                    .zip(ctx.identified.get(&(r.precursor, r.run)))
                    .and_then(|(p, id)| mismatch(p, id));
                if let Some(m) = wrong {
                    v.wrong_peak_pairs += 1;
                    v.wrong_peak += usize::from(m.apex);
                    v.wrong_measurement += usize::from(m.wrong());
                }
                for (a, &t) in [&mut v.at_q01, &mut v.at_threshold]
                    .into_iter()
                    .zip(&thresholds)
                {
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
            measurement_q: None,
            library_q: r.library_q,
            mass_err_ppm: e.map(|e| e.mass_err_ppm),
            claim: if r.role == Role::Transfer {
                r.claim
            } else {
                Claim::default()
            },
            transfer_pep: if r.role == Role::Transfer {
                scored[i].and_then(|(fam, s)| fdr.pep(r.precursor, fam, s))
            } else {
                None
            },
        });
    }
    Ok(())
}

/// What the block passes of `finalize_stream` share: the extract source and its bounds, the
/// library index, the run maps, the donor summaries, the limits and the seed.
#[derive(Clone, Copy)]
struct StreamInputs<'a> {
    source: &'a dyn ExtractSource,
    bounds: Bounds,
    idx: &'a EntryIndex<'a>,
    maps: &'a [Option<RunMap>],
    summaries: &'a Summaries,
    limits: &'a FinalizeLimits,
    seed: u64,
}

/// The first pass of `finalize_stream`: samples the training rows of every block (transfers,
/// identified quality, wrong peaks) and trains the fold models on them. Records the rows without
/// features in `summary`; returns the transfer, identified and wrong-peak models and the number of
/// wrong-peak rows seen.
fn train_transfer_models(
    inputs: &StreamInputs<'_>,
    summary: &mut FinalizeSummary,
) -> Result<(FoldModels, FoldModels, FoldModels, usize), String> {
    let StreamInputs {
        source,
        bounds,
        idx,
        maps,
        summaries,
        limits,
        seed,
    } = *inputs;
    let t0 = Instant::now();
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
            let feats = table_features(rows, idx, &shapes, maps, |t: Table| t != Table::Pseudo);
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
            let mut w_rows = Vec::new();
            for (i, r) in rows.iter().enumerate() {
                let Some((y, x)) = pseudo_example(r, &ctx, idx) else {
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
            log::debug!(
                target: "mbr",
                "{name} classifier: trained on {kept} of {seen} rows ((precursor, run) pairs, {share} of them; every row is scored)"
            );
        }
    }
    let t1 = Instant::now();
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
    log::debug!(
        target: "mbr",
        "wrong-peak model: {} pseudo-transfers ({:.1} % wrong against the identified fit), {} failed folds",
        w_rows.len(),
        100.0 * w_rows.iter().filter(|r| r.is_decoy).count() as f64 / w_rows.len().max(1) as f64,
        w_models.failed_folds.len()
    );
    for f in &w_models.failed_folds {
        log::warn!(target: "mbr", "wrong-peak model fold without a model: {f}");
    }
    drop(w_rows);
    log::debug!(
        target: "finalize",
        "training rows of {} extract rows {:.1} s, fold models {:.1} s ({} transfer classifier members per fold)",
        rows_read.into_inner(),
        (t1 - t0).as_secs_f64(),
        secs(t1),
        t_models.members()
    );
    Ok((t_models, i_models, w_models, w_seen.into_inner()))
}

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
    log::debug!(
        target: "mbr",
        "transfer FDR (internal quant-quality score) features: {}",
        TRANSFER_FEATURES.join(",")
    );
    log::info!(
        target: "mbr",
        "library entries {} targets, {} identity decoys",
        v.target_entries,
        v.identity_entries
    );
    for f in &fdr.counts.failed_folds {
        log::warn!(target: "mbr", "transfer FDR fold without a model: {f}");
    }
    log::info!(
        target: "mbr",
        "transfers {}: {} with a q-value; {} transfer rows without features, {} unscored, {} scored without a q-value",
        v.transfers,
        v.transfers_scored,
        summary.no_features,
        fdr.counts.unscored,
        fdr.counts.no_q
    );
    for f in &id_fdr.counts.failed_folds {
        log::warn!(target: "mbr", "identified quality fold without a model: {f}");
    }
    log::info!(
        target: "mbr",
        "identified rows {}: {} with a quality q-value (quant_q, internal), {} at q <= {VALIDATION_Q}; {} identified rows without features, {} unscored, {} scored without a q-value{}",
        v.identified,
        v.identified_scored,
        v.identified_q01,
        summary.identified_no_features,
        id_fdr.counts.unscored,
        id_fdr.counts.no_q,
        quant_feature_fdr.map_or(String::new(), |t| format!(
            "; {} identified_below_threshold at quant_q {t}",
            v.identified_below_threshold
        ))
    );
    log::info!(
        target: "mbr",
        "pseudo-transfers {}: wrong peak {} of {}, wrong measurement (wrong peak, or quantity off by more than a factor of 2) {}",
        v.pseudo_transfers, v.wrong_peak, v.wrong_peak_pairs, v.wrong_measurement
    );
    let rate = |a: usize, b: usize| {
        if b > 0 {
            format!("{:.4}", a as f64 / b as f64)
        } else {
            "NA".to_string()
        }
    };
    log::info!(
        target: "mbr",
        "MS2 support (at least {} fragments at angle >= {}): transfer targets {} of {} ({}), mass-shift decoys {} of {} ({}), RT decoys {} of {} ({}), pseudo-transfers {} of {} ({}); {} accepted transfers are transfer_ms2",
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
    log::info!(
        target: "mbr",
        "shared-feature classifier (probability above {}): {} accepted transfers reported as shared_feature",
        quant::shared::Model::shipped().threshold,
        v.shared_feature
    );
    log::info!(
        target: "mbr",
        "run normalisation (quantity_normalized): {} reference precursors (identified, quant_q <= {VALIDATION_Q}, in >= 3 runs or all runs)",
        norm.precursors
    );
    for (run, stem) in stems.iter().enumerate() {
        match norm.summary(run as u16) {
            Some((level, sd, n)) => log::info!(
                target: "mbr",
                "  {stem}: median offset {level:+.3} log2, sd along RT {sd:.3}, {n} reference cells"
            ),
            None => log::info!(target: "mbr", "  {stem}: not normalised (too few reference cells)"),
        }
    }
    let shown: &[(f64, &Accepted)] = if mbr_fdr == VALIDATION_Q {
        &[(VALIDATION_Q, &v.at_q01)]
    } else {
        &[(VALIDATION_Q, &v.at_q01), (mbr_fdr, &v.at_threshold)]
    };
    for &(t, a) in shown {
        log::info!(
            target: "mbr",
            "at q <= {t}: decoy transfers accepted: mass-shift {} of {}, RT-shift {} of {}, identity decoys {} of {}; {} transfers would pass on the wrong-peak decoys alone",
            a.decoys[0],
            v.decoy_transfers[0],
            a.decoys[1],
            v.decoy_transfers[1],
            a.decoys[2],
            v.decoy_transfers[2],
            a.wrong_peak_only
        );
        log::info!(
            target: "mbr",
            "at q <= {t}: {} transfers (mass-shift decoys only {}, RT decoys only {}); {} pseudo-transfers (acceptance rate {}, optimistic for transfers: the alignment kept their runs' anchors), {} of {} accepted at the wrong peak, {} with a wrong measurement",
            a.pooled,
            a.mass_only,
            a.rt_only,
            a.pseudo,
            v.pseudo_acceptance(a).map_or("NA".into(), |s| format!("{s:.4}")),
            a.wrong_peak,
            a.pairs,
            a.wrong_measurement
        );
    }
}

fn transfer_accepted(q: Option<f64>, mbr_fdr: f64, has_quantity: bool) -> bool {
    has_quantity && q.is_some_and(|q| q <= mbr_fdr)
}

#[cfg(test)]
#[path = "finalize/tests/finalize.rs"]
mod tests;
