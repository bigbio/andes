//! The run loop of a TMT study: the run tables are decoded a window at a time in parallel
//! while the previous window is joined, each table is joined and quantified in table order,
//! and two writer threads encode the study table and the feature view in that order, so the
//! files and every count are those of a serial loop.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{sync_channel, SyncSender};
use std::time::{Duration, Instant};

use output::features::FeatureRecord;
use output::tmt::{RunDesignColumns, TmtStudyLayout, TmtStudyRow, TmtStudyWriter};
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
/// Runs queued for each writer thread.
const WRITER_QUEUE: usize = 2;

/// The error a send gets when its writer thread has stopped; the writer's own error is
/// returned instead.
const WRITER_STOPPED: &str = "a TMT study writer stopped";

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
    /// Time spent joining, building and quantifying features, and waiting on the writers.
    pub spent: [Duration; 3],
}

/// The layout of the study table and the queues of the writer threads.
struct Sinks {
    layout: TmtStudyLayout,
    mode: String,
    study: SyncSender<(Vec<TmtStudyRow>, RunDesignColumns)>,
    features: Option<SyncSender<Vec<FeatureRecord>>>,
}

impl StudyJoin<'_> {
    /// Joins every table of `tables` in order: up to [`TABLE_WINDOW`] tables are decoded in
    /// parallel while the tables decoded before them are joined. An error is the one the
    /// serial loop meets first, or a writer's error.
    pub(super) fn run(&mut self, tables: &[PathBuf]) -> Result<(), String> {
        let window = rayon::current_num_threads().clamp(1, TABLE_WINDOW);
        let decode = |paths: &[PathBuf]| -> Vec<Result<RunTable, String>> {
            paths.par_iter().map(|p| read_run_table(p)).collect()
        };
        let mut chunks = tables.chunks(window);
        let Some(first) = chunks.next() else {
            return Ok(());
        };
        let decoded = decode(first);
        let (w, layout, mode) = match decoded.first() {
            Some(Ok(t)) => self.open(t)?,
            Some(Err(e)) => return Err(e.clone()),
            None => return Ok(()),
        };
        let features = self.features.take();
        let (study_tx, study_rx) =
            sync_channel::<(Vec<TmtStudyRow>, RunDesignColumns)>(WRITER_QUEUE);
        let (feature_tx, feature_rx) = sync_channel::<Vec<FeatureRecord>>(WRITER_QUEUE);
        std::thread::scope(|s| {
            let study = s.spawn(move || -> Result<TmtStudyWriter, String> {
                let mut w = w;
                for (rows, d) in study_rx {
                    w.push_run(&rows, &d)?;
                }
                Ok(w)
            });
            let has_features = features.is_some();
            let feats = features.map(|(mut fw, path)| {
                s.spawn(move || -> Result<(FeatureParquetWriter, PathBuf), String> {
                    for recs in feature_rx {
                        fw.push(recs)
                            .map_err(|e| format!("write {}: {e}", path.display()))?;
                    }
                    Ok((fw, path))
                })
            });
            let mut sinks = Sinks {
                layout,
                mode,
                study: study_tx,
                features: has_features.then_some(feature_tx),
            };
            let joined = self.join_tables(&mut sinks, decoded, first, &mut chunks, &decode);
            let Sinks {
                layout,
                mode,
                study: study_tx,
                features: feature_tx,
            } = sinks;
            drop((study_tx, feature_tx));
            let study = study
                .join()
                .unwrap_or_else(|e| std::panic::resume_unwind(e));
            let feats = feats.map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)));
            let written = study.and_then(|w| {
                let fw = feats.transpose()?;
                Ok((w, fw))
            });
            match (joined, written) {
                (_, Err(e)) | (Err(e), Ok(_)) => Err(e),
                (Ok(()), Ok((w, fw))) => {
                    self.writer = Some((w, layout, mode));
                    self.features = fw;
                    Ok(())
                }
            }
        })
    }

    /// Opens the study table on the first run table.
    fn open(&self, table: &RunTable) -> Result<(TmtStudyWriter, TmtStudyLayout, String), String> {
        let meta = study_meta(
            self.a,
            table,
            self.a.proteins.psm_fdr,
            self.n_tables,
            self.design,
        );
        let w = TmtStudyWriter::create(self.out_path, table.layout, &meta)?;
        Ok((w, table.layout, table.unmix_mode.clone()))
    }

    /// Joins the decoded tables of `first`, then those of the following chunks, each window
    /// decoded while the one before it is joined.
    fn join_tables<'p>(
        &mut self,
        sinks: &mut Sinks,
        decoded: Vec<Result<RunTable, String>>,
        first: &'p [PathBuf],
        chunks: &mut std::slice::Chunks<'p, PathBuf>,
        decode: &(impl Fn(&[PathBuf]) -> Vec<Result<RunTable, String>> + Sync),
    ) -> Result<(), String> {
        let mut next = Some((first, decoded));
        while let Some((paths, decoded)) = next.take() {
            let ahead = chunks.next();
            let (done, read) = rayon::join(
                || -> Result<(), String> {
                    for (path, table) in paths.iter().zip(decoded) {
                        self.push_table(sinks, path, table?)?;
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

    /// Joins one run table: checks it matches the layout of the study's first table, then
    /// joins, quantifies and queues each of its runs.
    fn push_table(&mut self, sinks: &Sinks, path: &Path, table: RunTable) -> Result<(), String> {
        if sinks.layout != table.layout || sinks.mode != table.unmix_mode {
            return Err(format!(
                "{}: {} (correction {}, unmixing {}) differs from the study's first table \
                 ({}, correction {}, unmixing {}); search every run with the same --tmt, \
                 --tmt-correction and --tmt-unmix",
                path.display(),
                table.layout.plex.name(),
                table.layout.corrected,
                table.unmix_mode,
                sinks.layout.plex.name(),
                sinks.layout.corrected,
                sinks.mode
            ));
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
            self.push_run(sinks, run, rows)?;
        }
        Ok(())
    }

    /// Joins one run's rows, logs them, builds and quantifies their features and queues the
    /// rows and the features for the writers.
    fn push_run(&mut self, sinks: &Sinks, run: String, rows: Vec<RunRow>) -> Result<(), String> {
        let psm_fdr = self.a.proteins.psm_fdr;
        let plex = sinks.layout.plex;
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
        let tic = Instant::now();
        let recs = match &sinks.features {
            Some(_) => {
                let unit = d.plex.map_or_else(|| run.clone(), |p| format!("plex {p}"));
                Some(self.features_of(&rows, plex, &unit)?)
            }
            None => None,
        };
        self.spent[1] += tic.elapsed();
        let tic = Instant::now();
        sinks
            .study
            .send((rows, d))
            .map_err(|_| WRITER_STOPPED.to_string())?;
        if let (Some(tx), Some(recs)) = (&sinks.features, recs) {
            tx.send(recs).map_err(|_| WRITER_STOPPED.to_string())?;
        }
        self.spent[2] += tic.elapsed();
        self.runs.push(run);
        Ok(())
    }

    /// The features of a run's grouped rows, built in parallel, then numbered and added to
    /// the protein quantities in row order.
    fn features_of(
        &mut self,
        rows: &[TmtStudyRow],
        plex: Plex,
        unit: &str,
    ) -> Result<Vec<FeatureRecord>, String> {
        let built: Vec<Option<FeatureRecord>> = rows
            .par_iter()
            .filter(|r| r.group.is_some())
            .map(|r| feature_record(r, plex, 0))
            .collect();
        let mut recs = Vec::with_capacity(built.len());
        for f in built {
            match f {
                Some(mut f) => {
                    self.feature_id += 1;
                    f.feature_id = self.feature_id;
                    self.tally.quantify(&f, unit)?;
                    recs.push(f);
                }
                None => self.no_mz += 1,
            }
        }
        Ok(recs)
    }
}
