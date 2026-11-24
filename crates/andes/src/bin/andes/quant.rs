//! `andes quant`: the library, extract and finalize steps of a study, run as subcommands or in
//! one process after the search.

mod extract;
pub(crate) mod finalize;
mod groups;
mod library;
pub(crate) mod mbr;
pub(crate) mod run;
mod scores;
pub(crate) mod tmt;

use std::path::{Path, PathBuf};

use output::extract::read_extract_tags;
use output::library::{library_fingerprint, read_library, StoredLibrary, FINGERPRINT_KEY};
pub(crate) use output::qpx::parse_peptidoform;
use output::FeatureParquetWriter;
use quant::evidence::EvidenceParams;

use crate::cli::{QuantCommand, QuantExtractArgs, QuantFinalizeArgs, QuantLibraryArgs};
use crate::quant::finalize::{
    finalize_stream, ExtractFiles, ExtractSource, FinalizeLimits, FinalizeParams,
};
use crate::quant::mbr::{mbr_paths, TransferSelection};
use extract::extract_and_write;
pub(crate) use extract::{evidence_params, scan_of, split_key};
use groups::feature_records;
pub(crate) use groups::{GroupTally, ProteinGroups};
pub(crate) use library::build_and_write_library;
pub(crate) use scores::{rescore_scores, ScoreLoader, ScoreSet};
use scores::{score_loader, score_source};

pub(crate) const MBR_SEED: u64 = 42;

pub(crate) const TABLE_SUFFIX: &str = ".quant.parquet";
const MS1_SUFFIX: &str = ".ms1.parquet";
pub(crate) const MS2_SUFFIX: &str = ".ms2.parquet";
const EXTRACT_SUFFIX: &str = ".extract.parquet";
const TOL_KEY: &str = "andes.extract.tol_ppm";
const SELECTION_KEY: &str = "andes.extract.transfer_selection";
const RUN_SCALES_KEY: &str = "andes.extract.run_scales";
const ENVELOPE_LEN: usize = 4;

pub(crate) fn artifact(dir: &Path, stem: &str, suffix: &str) -> PathBuf {
    dir.join(format!("{stem}{suffix}"))
}

pub(crate) fn library_runs(dir: &Path, runs: &[String]) -> Result<Vec<String>, String> {
    let mut stems: Vec<String> = if runs.is_empty() {
        let entries = std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        let mut v = Vec::new();
        for e in entries {
            let name = e
                .map_err(|e| format!("{}: {e}", dir.display()))?
                .file_name()
                .to_string_lossy()
                .into_owned();
            if let Some(stem) = name.strip_suffix(TABLE_SUFFIX) {
                v.push(stem.to_string());
            }
        }
        v
    } else {
        runs.to_vec()
    };
    stems.sort();
    if stems.is_empty() {
        return Err(format!(
            "no runs: {} holds no *{TABLE_SUFFIX} and no --runs were given",
            dir.display()
        ));
    }
    if let Some(w) = stems.windows(2).find(|w| w[0] == w[1]) {
        return Err(format!(
            "run {:?} is given twice; run stems must be unique",
            w[0]
        ));
    }
    if stems.len() > u16::MAX as usize {
        return Err(format!("{} runs exceed the supported number", stems.len()));
    }
    let missing: Vec<String> = stems
        .iter()
        .flat_map(|s| {
            [
                artifact(dir, s, TABLE_SUFFIX),
                artifact(dir, s, MS1_SUFFIX),
                artifact(dir, s, MS2_SUFFIX),
            ]
        })
        .filter(|p| !p.is_file())
        .map(|p| p.display().to_string())
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "missing per-run artifacts (written by `andes search --quant-dir`): {}",
            missing.join(", ")
        ));
    }
    Ok(stems)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn finalize_and_write(
    stored: &StoredLibrary,
    source: &dyn ExtractSource,
    mbr_fdr: f64,
    seed: u64,
    limits: FinalizeLimits,
    trace_fraction: f64,
    base: &Path,
    parquet_dir: Option<&Path>,
    groups: Option<&ProteinGroups<'_>>,
) -> Result<(), String> {
    let (table, validation) = mbr_paths(base);
    let mut feature_out = match parquet_dir {
        Some(dir) => {
            std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
            let path = dir.join("quantms.feature.parquet");
            Some((FeatureParquetWriter::new_mbr(&path), path, 0i64))
        }
        None => None,
    };
    let params = FinalizeParams {
        mbr_fdr,
        quant_feature_fdr: None,
        seed,
        limits,
        mbr_path: &table,
        validation_path: &validation,
        trace_fraction,
    };
    let mut tally = GroupTally::default();
    let summary = finalize_stream(&stored.runs, &stored.library, source, &params, &mut |f| {
        if let Some((w, path, last_id)) = feature_out.as_mut() {
            let recs = feature_records(stored, &f, last_id, groups, &mut tally);
            for r in &recs {
                tally.quantify(r, &r.run_file_name)?;
            }
            w.push(recs)
                .map_err(|e| format!("write {}: {e}", path.display()))?;
        }
        Ok(())
    })?;
    log::info!(
        target: "mbr",
        "wrote {} ({} rows; transfer_q and quant_q are internal quant-quality q-values, not identification FDR) and {}",
        table.display(),
        summary.rows_written,
        validation.display()
    );
    if let Some((w, path, n)) = feature_out {
        w.finish()
            .map_err(|e| format!("write {}: {e}", path.display()))?;
        log::info!(
            target: "mbr",
            "wrote {} ({n} feature rows with evidence type; below-threshold transfers and identified rows are in the mbr.parquet only)",
            path.display(),
        );
    }
    if let (Some(dir), Some(pg)) = (parquet_dir, groups) {
        tally.write_pg(dir, pg, &stored.runs)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_in_process(
    dir: &Path,
    stems: &[String],
    scores: ScoreSet,
    ep: &EvidenceParams,
    max_q: f64,
    donor_q: f64,
    donor_pep: Option<f64>,
    mbr_fdr: f64,
    base: &Path,
    parquet_dir: Option<&Path>,
    groups: Option<&ProteinGroups<'_>>,
) -> Result<(), String> {
    let t0 = std::time::Instant::now();
    let stems = library_runs(dir, stems)?;
    let stored = build_and_write_library(
        dir,
        &stems,
        ScoreLoader::Merged(Box::new(move |_| Ok(scores))),
        max_q,
        donor_q,
        donor_pep,
    )?;
    let fingerprint = library_fingerprint(&stored);
    for r in 0..stems.len() {
        extract_and_write(
            dir,
            &stored,
            &fingerprint,
            r as u16,
            ep,
            &TransferSelection::default(),
        )?;
    }
    let files = ExtractFiles(
        stems
            .iter()
            .map(|s| artifact(dir, s, EXTRACT_SUFFIX))
            .collect(),
    );
    finalize_and_write(
        &stored,
        &files,
        mbr_fdr,
        MBR_SEED,
        FinalizeLimits::default(),
        crate::quant::finalize::DEFAULT_TRACE_FRACTION,
        base,
        parquet_dir,
        groups,
    )?;
    log::info!(target: "mbr", "done in {:.2} s", t0.elapsed().as_secs_f64());
    Ok(())
}

fn run_library(a: QuantLibraryArgs) -> Result<(), String> {
    let stems = library_runs(&a.quant_dir, &a.runs)?;
    let source = score_source(&a.psms, &a.decoy_psms, &stems)?;
    build_and_write_library(
        &a.quant_dir,
        &stems,
        score_loader(&source),
        a.mbr_donor_fdr(),
        a.psm_fdr,
        a.pep,
    )
    .map(|_| ())
}

fn run_extract(a: QuantExtractArgs) -> Result<(), String> {
    if crate::dia::study::is_study_run(&a.quant_dir, &a.run) {
        return crate::dia::transfer::extract_run(&a.quant_dir, &a.run);
    }
    let ep = evidence_params(a.lfq_tol)?;
    let stored = read_library(&a.quant_dir)?;
    let run = stored
        .runs
        .iter()
        .position(|s| *s == a.run)
        .ok_or_else(|| {
            format!(
                "run {:?} is not in the library of {} (runs: {})",
                a.run,
                a.quant_dir.display(),
                stored.runs.join(", ")
            )
        })?;
    let fingerprint = library_fingerprint(&stored);
    let selection = if a.mbr_preselect {
        TransferSelection {
            min_donor_fraction: a.mbr_min_donor_fraction,
        }
    } else {
        TransferSelection::OFF
    };
    extract_and_write(
        &a.quant_dir,
        &stored,
        &fingerprint,
        run as u16,
        &ep,
        &selection,
    )
    .map(|_| ())
}

fn run_finalize(a: QuantFinalizeArgs) -> Result<(), String> {
    let dia_runs = crate::dia::study::study_runs(&a.quant_dir)?;
    if !dia_runs.is_empty() {
        return crate::dia::study::finalize(&a, &dia_runs);
    }
    if !a.quant_dir.join(output::library::LIBRARY_FILE).is_file() {
        let tables = crate::quant::tmt::study_tables(&a.quant_dir, &a.output)?;
        if !tables.is_empty() {
            return crate::quant::tmt::finalize(&a, &tables);
        }
    }
    let stored = read_library(&a.quant_dir)?;
    let fingerprint = library_fingerprint(&stored);
    let mut paths = Vec::with_capacity(stored.runs.len());
    let mut tol: Option<(String, PathBuf)> = None;
    for stem in &stored.runs {
        let path = artifact(&a.quant_dir, stem, EXTRACT_SUFFIX);
        let tags = read_extract_tags(&path)?;
        let get = |k: &str| tags.iter().find(|t| t.0 == k).map(|t| t.1.clone());
        if get(FINGERPRINT_KEY).as_deref() != Some(fingerprint.as_str()) {
            return Err(format!(
                "{}: extracted with another library than {} (stale extract); run `andes quant extract` again",
                path.display(),
                a.quant_dir.join(output::library::LIBRARY_FILE).display()
            ));
        }
        let t = get(TOL_KEY).ok_or_else(|| format!("{}: no {TOL_KEY}", path.display()))?;
        match &tol {
            None => tol = Some((t, path.clone())),
            Some((t0, p0)) if *t0 != t => {
                return Err(format!(
                    "{}: extracted at {t} ppm, {} at {t0} ppm; every run needs the same --lfq-tol",
                    path.display(),
                    p0.display()
                ));
            }
            _ => {}
        }
        if let Some(scales) = get(RUN_SCALES_KEY) {
            log::debug!(target: "mbr", "run {stem}: {scales}");
        }
        paths.push(path);
    }
    let transfer_fdr = a.quant_transfer_fdr();
    let limits = FinalizeLimits {
        max_train_rows: a.max_train_rows,
        transfer_models: a.transfer_models,
        max_norm_cells: a.max_norm_cells,
        max_donor_fits: a.max_donor_fits,
        max_held_scores: a.max_held_scores,
        temp_dir: a.temp_dir,
    };
    let inf = match (&a.output_parquet, &a.database) {
        (Some(_), Some(database)) => Some(
            crate::inference::infer_from_files(database, &a.proteins, &a.psms, &a.decoy_psms)
                .map_err(|e| e.to_string())?
                .0,
        ),
        _ => None,
    };
    if inf.is_some() {
        crate::memory::release_free_memory();
    }
    let groups = inf.as_ref().map(|inf| ProteinGroups {
        inf,
        max_q: a.proteins.protein_fdr,
    });
    finalize_and_write(
        &stored,
        &ExtractFiles(paths),
        transfer_fdr,
        a.seed,
        limits,
        a.quant_trace_fraction,
        &a.output,
        a.output_parquet.as_deref(),
        groups.as_ref(),
    )
}

pub(crate) fn run_quant_command(cmd: QuantCommand) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(n) = cmd.threads() {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global();
    }
    match cmd {
        QuantCommand::Library(a) => run_library(a),
        QuantCommand::Extract(a) => run_extract(a),
        QuantCommand::Finalize(a) => run_finalize(*a),
    }
    .map_err(Into::into)
}

#[cfg(test)]
#[path = "quant/tests/quant.rs"]
mod tests;
