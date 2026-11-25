//! The run loop of a TMT study: the run tables are decoded a window at a time in parallel
//! while the previous window is joined, and each table is joined, written and quantified in
//! table order, so the study table, the QPX views and every count are those of a serial loop.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use output::features::FeatureRecord;
use output::tmt::{TmtStudyLayout, TmtStudyRow, TmtStudyWriter};
use output::FeatureParquetWriter;
use quant::isobaric::Plex;
use rayon::prelude::*;

use super::{
    feature_record, read_run_table, run_design, study_meta, study_rows, Pooled, RunRow, RunTable,
    StudyCounts,
};
use crate::quant::{GroupTally, ProteinGroups};
use crate::sdrf::IsobaricDesign;

/// Most run tables decoded at once ahead of the join (the window in memory is at most twice
/// this many tables).
const TABLE_WINDOW: usize = 8;

/// The study's state while its run tables are joined in table order.
pub(super) struct StudyJoin<'a> {
    pub a: &'a crate::cli::QuantFinalizeArgs,
    pub design: Option<&'a IsobaricDesign>,
    pub groups: &'a ProteinGroups<'a>,
    pub out_path: &'a Path,
    pub n_tables: usize,
    pub pooled: HashMap<u128, Pooled>,
    pub writer: Option<(TmtStudyWriter, TmtStudyLayout, String)>,
    pub features: Option<(FeatureParquetWriter, PathBuf)>,
    pub tally: GroupTally,
    pub counts: StudyCounts,
    pub runs: Vec<String>,
    pub no_mz: usize,
    pub feature_id: i64,
    /// Time spent joining, writing the study table, and building and writing features.
    pub spent: [Duration; 3],
}

impl StudyJoin<'_> {
    /// Joins every table of `tables` in order: up to [`TABLE_WINDOW`] tables are decoded in
    /// parallel while the tables decoded before them are joined. An error is the one the
    /// serial loop meets first.
    pub(super) fn run(&mut self, tables: &[PathBuf]) -> Result<(), String> {
        let window = rayon::current_num_threads().clamp(1, TABLE_WINDOW);
        let decode = |paths: &[PathBuf]| -> Vec<Result<RunTable, String>> {
            paths.par_iter().map(|p| read_run_table(p)).collect()
        };
        let mut chunks = tables.chunks(window);
        let mut next = chunks.next().map(|c| (c, decode(c)));
        while let Some((paths, decoded)) = next.take() {
            let ahead = chunks.next();
            let (done, read) = rayon::join(
                || -> Result<(), String> {
                    for (path, table) in paths.iter().zip(decoded) {
                        self.push_table(path, table?)?;
                    }
                    Ok(())
                },
                || ahead.map(|c| (c, decode(c))),
            );
            done?;
            next = read;
        }
        Ok(())
    }

    /// Joins one run table: opens the study table on the first one and checks the others
    /// match its layout, then joins, writes and quantifies each of its runs.
    fn push_table(&mut self, path: &Path, table: RunTable) -> Result<(), String> {
        match &self.writer {
            None => {
                let meta = study_meta(
                    self.a,
                    &table,
                    self.a.proteins.psm_fdr,
                    self.n_tables,
                    self.design,
                );
                let w = TmtStudyWriter::create(self.out_path, table.layout, &meta)?;
                self.writer = Some((w, table.layout, table.unmix_mode.clone()));
            }
            Some((_, layout, mode)) if *layout != table.layout || *mode != table.unmix_mode => {
                return Err(format!(
                    "{}: {} (correction {}, unmixing {}) differs from the study's first table \
                     ({}, correction {}, unmixing {}); search every run with the same --tmt, \
                     --tmt-correction and --tmt-unmix",
                    path.display(),
                    table.layout.plex.name(),
                    table.layout.corrected,
                    table.unmix_mode,
                    layout.plex.name(),
                    layout.corrected,
                    mode
                ));
            }
            Some(_) => {}
        }
        let mut by_run: Vec<(String, Vec<RunRow>)> = Vec::new();
        for r in table.rows {
            match by_run.iter_mut().find(|(n, _)| *n == r.filename) {
                Some((_, v)) => v.push(r),
                None => by_run.push((r.filename.clone(), vec![r])),
            }
        }
        for (run, rows) in by_run {
            if self.runs.contains(&run) {
                return Err(format!(
                    "{}: run {run} is also in another table of {}",
                    path.display(),
                    self.a.quant_dir.display()
                ));
            }
            self.push_run(run, rows)?;
        }
        Ok(())
    }

    /// Joins one run's rows, logs them, and writes them to the study table while their
    /// features are built and quantified.
    fn push_run(&mut self, run: String, rows: Vec<RunRow>) -> Result<(), String> {
        let psm_fdr = self.a.proteins.psm_fdr;
        let Some((w, layout, _)) = self.writer.as_mut() else {
            return Err(format!("run {run}: no study table open"));
        };
        let plex = layout.plex;
        let d = run_design(self.design, &run, plex)?;
        let before = self.counts.clone();
        let tic = Instant::now();
        let rows = study_rows(
            rows,
            &mut self.pooled,
            self.groups,
            &mut self.tally,
            &mut self.counts,
        )?;
        self.spent[0] += tic.elapsed();
        let counts = &self.counts;
        if rows.is_empty() && counts.run_rows > before.run_rows {
            log::warn!(
                "run {run}: none of its {} PSMs is in the pooled results at q <= {psm_fdr}; \
                 were its PIN rows rescored with the others?",
                counts.run_rows - before.run_rows
            );
        }
        log::info!(
            target: "tmt",
            "run {run}: {} PSMs with reporters, {} at pooled q <= {psm_fdr}, {} passing the \
             gate, {} in a protein group at q <= {}{}",
            counts.run_rows - before.run_rows,
            counts.kept - before.kept,
            counts.gate_passed - before.gate_passed,
            counts.grouped - before.grouped,
            self.a.proteins.protein_fdr,
            d.plex.map_or(String::new(), |p| format!(
                " (plex {p}, fraction {})",
                d.fraction.as_deref().unwrap_or("none")
            ))
        );
        let unit = d.plex.map_or_else(|| run.clone(), |p| format!("plex {p}"));
        let feats = Features {
            out: self.features.as_mut(),
            tally: &mut self.tally,
            feature_id: &mut self.feature_id,
            no_mz: &mut self.no_mz,
        };
        let ((wrote, t_table), (featured, t_feat)) = rayon::join(
            || {
                let tic = Instant::now();
                (w.push_run(&rows, &d), tic.elapsed())
            },
            || {
                let tic = Instant::now();
                (feats.push(&rows, plex, &unit), tic.elapsed())
            },
        );
        self.spent[1] += t_table;
        self.spent[2] += t_feat;
        wrote?;
        featured?;
        self.runs.push(run);
        Ok(())
    }
}

/// The feature view of the study and what it adds to.
struct Features<'a> {
    out: Option<&'a mut (FeatureParquetWriter, PathBuf)>,
    tally: &'a mut GroupTally,
    feature_id: &'a mut i64,
    no_mz: &'a mut usize,
}

impl Features<'_> {
    /// Builds the features of a run's grouped rows in parallel, then numbers them, adds them
    /// to the protein quantities and writes them in row order.
    fn push(self, rows: &[TmtStudyRow], plex: Plex, unit: &str) -> Result<(), String> {
        let Some((fw, fpath)) = self.out else {
            return Ok(());
        };
        let built: Vec<Option<FeatureRecord>> = rows
            .par_iter()
            .filter(|r| r.group.is_some())
            .map(|r| feature_record(r, plex, 0))
            .collect();
        let mut recs = Vec::with_capacity(built.len());
        for f in built {
            match f {
                Some(mut f) => {
                    *self.feature_id += 1;
                    f.feature_id = *self.feature_id;
                    self.tally.quantify(&f, unit)?;
                    recs.push(f);
                }
                None => *self.no_mz += 1,
            }
        }
        fw.push(recs)
            .map_err(|e| format!("write {}: {e}", fpath.display()))
    }
}
