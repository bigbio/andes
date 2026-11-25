//! The last steps of a search: quantification, the proteoform table, the MBR job, the search TSV
//! and the release of the search's memory before the in-process MBR run.

use std::collections::{HashMap, HashSet};

use std::path::PathBuf;

use model::Spectrum;
use rayon::prelude::*;
use search::TopNQueue;

use super::outputs::{release_at_exit, write_proteoform_table, write_search_tsv, PsmTables};
use super::*;
use crate::cli::Cli;
use crate::quant::run::QuantCollector;

/// What the last steps of a search own: everything they write from or drop.
pub(super) struct Finish {
    pub(super) cli: Cli,
    pub(super) bench_mode: bool,
    pub(super) quant: Option<QuantCollector>,
    pub(super) spectra: Vec<Spectrum>,
    pub(super) queues: Vec<TopNQueue>,
    pub(super) target_db: model::protein::ProteinDb,
    pub(super) mbr_decoys: HashMap<String, output::PercolatorPsm>,
    pub(super) report_base: PathBuf,
    pub(super) params: search::SearchParams,
    pub(super) is_mgf: bool,
    pub(super) rescore_map: Option<HashMap<String, output::PercolatorPsm>>,
    pub(super) mbr_source: String,
    pub(super) rescore_decoys: Option<HashSet<String>>,
    pub(super) pin_candidates: Vec<search::candidates::Candidate>,
    pub(super) pin_index: search::SearchIndex,
    pub(super) knockoff_rows: HashSet<String>,
    pub(super) protein_groups: Option<output::Inference>,
    pub(super) cutoffs: Option<output::Cutoffs>,
}

/// Quantifies, writes the proteoform table and the search TSV, then runs MBR in process once the
/// search's memory is released.
pub(super) fn finish_outputs(finish: Finish) -> Result<(), Box<dyn std::error::Error>> {
    let Finish {
        cli,
        bench_mode,
        quant,
        spectra,
        queues,
        target_db,
        mbr_decoys,
        report_base,
        params,
        is_mgf,
        mut rescore_map,
        mbr_source,
        rescore_decoys,
        pin_candidates,
        pin_index,
        knockoff_rows,
        protein_groups,
        cutoffs,
    } = finish;
    let spectrum_paths = &cli.spectrum;
    let omitted_psms = match &quant {
        Some(_) => output::omitted_psm_rows(
            &spectra,
            &queues,
            &pin_candidates,
            rescore_map.as_ref(),
            rescore_decoys.as_ref(),
            &knockoff_rows,
        )
        .map_err(|e| format!("quant: {e}"))?,
        None => std::collections::HashSet::new(),
    };
    if let Some(ref mut map) = rescore_map {
        map.retain(|id, _| !knockoff_rows.contains(id));
    }

    let qpx_groups = protein_groups
        .as_ref()
        .map(|inf| crate::quant::ProteinGroups {
            inf,
            max_q: cli.protein_fdr.unwrap_or(0.01),
        });
    let tmt_mbr_input = (cli.mbr && cli.tmt.is_some()).then(|| crate::quant::run::TmtMbrInput {
        decoys: &mbr_decoys,
        source: &mbr_source,
        library_fdr: cli.mbr_donor_fdr(),
        accept_q: cli.quant_transfer_fdr(),
        wrong_log2: cli.mbr_wrong_log2,
    });
    let quant_summary = match &quant {
        Some(q) => Some(crate::quant::run::run_quant(
            q,
            &spectra,
            &queues,
            &pin_candidates,
            &pin_index,
            rescore_map.as_ref(),
            &omitted_psms,
            cutoffs.as_ref(),
            &report_base,
            cli.output_parquet.as_deref(),
            qpx_groups.as_ref(),
            tmt_mbr_input.as_ref(),
            (cli.mbr && cli.label.is_some()).then_some(&mbr_decoys),
        )?),
        None => None,
    };

    if let (Some(map), Some(cutoffs), true) = (
        rescore_map.as_ref(),
        cutoffs.as_ref(),
        pin_index.has_alternative_entries,
    ) {
        write_proteoform_table(
            map,
            &pin_index,
            cutoffs,
            quant_summary.as_ref(),
            cli.protein_fdr.zip(protein_groups.as_ref()),
            &report_base,
        );
    }

    let mbr_job = if cli.mbr && cli.tmt.is_none() && cli.label.is_none() {
        let (Some(dir), Some(map)) = (cli.quant_dir.clone(), rescore_map.as_ref()) else {
            return Err("--mbr needs --quant-dir and rescoring".into());
        };
        let stems: Vec<String> = spectrum_paths
            .iter()
            .map(|p| crate::quant::run::run_stem(p))
            .collect();
        let scores = crate::quant::rescore_scores(map, &mbr_decoys, &mbr_source)?;
        let ep = crate::quant::evidence_params(cli.lfq_tol)?;
        let (donor_q, donor_pep) = quant
            .as_ref()
            .map(|q| (q.settings.psm_fdr, q.settings.pep_cap))
            .ok_or("--mbr needs --lfq")?;
        Some((dir, stems, scores, ep, donor_q, donor_pep))
    } else {
        None
    };
    drop(mbr_decoys);

    if bench_mode {
        log::info!(target: "bench", "skipping TSV write.");
    } else if let Some(ref tsv_path) = cli.output_tsv {
        write_search_tsv(
            tsv_path,
            &PsmTables {
                spectra: &spectra,
                queues: &queues,
                candidates: &pin_candidates,
                params: &params,
                index: &pin_index,
            },
            spectrum_paths,
            is_mgf,
        )?;
    }

    crate::memory::mark_stage("outputs");
    let Some((dir, stems, scores, ep, donor_q, donor_pep)) = mbr_job else {
        release_at_exit(pin_candidates);
        return Ok(());
    };
    drop(spectra);
    drop(queues);
    drop(pin_candidates);
    drop(pin_index);
    drop(quant);
    drop(rescore_map);
    drop(target_db);
    crate::quant::run_in_process(
        &dir,
        &stems,
        scores,
        &ep,
        cli.mbr_donor_fdr(),
        donor_q,
        donor_pep,
        cli.quant_transfer_fdr(),
        &report_base.with_extension(""),
        cli.output_parquet.as_deref(),
        qpx_groups.as_ref(),
    )?;
    Ok(())
}

/// What the steps after rescoring read and own: the search's results, its rescoring results and
/// what the report and outputs are written from.
pub(super) struct AfterRescoring<'a> {
    pub(super) cli: Cli,
    pub(super) spectrum_path: PathBuf,
    pub(super) bench_mode: bool,
    pub(super) quant: Option<QuantCollector>,
    pub(super) run_ends: Vec<usize>,
    pub(super) measurer: FragmentErrorMeasurer,
    pub(super) file_calibrations: Vec<search::FileCalibration>,
    pub(super) spectra: Vec<Spectrum>,
    pub(super) queues: Vec<TopNQueue>,
    pub(super) target_db: model::protein::ProteinDb,
    pub(super) ipf_scoring: Option<search::engine::IpfScoring<'a>>,
    pub(super) run_stats: output::RunStatistics,
    pub(super) fact_rows: Option<Vec<Box<[output::FactRow]>>>,
    pub(super) dia_pass2: Option<output::dia::pass2::Pass2Table>,
    pub(super) rerank_dropped: HashSet<String>,
    pub(super) fdr_threshold: f64,
    pub(super) proteoform_fdr: f64,
    pub(super) use_native: bool,
    pub(super) mbr_decoys: HashMap<String, output::PercolatorPsm>,
    pub(super) original_spectrum_paths: Vec<PathBuf>,
    pub(super) output_pin_path: PathBuf,
    pub(super) report_base: PathBuf,
    pub(super) params: search::SearchParams,
    pub(super) restricted_site_masses: Vec<f64>,
    pub(super) fixed_terminal_mods: Vec<(model::modification::ModLocation, f64)>,
    pub(super) scorer: &'a RankScorer,
    pub(super) is_mgf: bool,
    pub(super) rescore_map: Option<HashMap<String, output::PercolatorPsm>>,
    pub(super) mbr_source: String,
    pub(super) rescore_decoys: Option<HashSet<String>>,
    pub(super) rescore_decoy_psms: HashMap<String, output::PercolatorPsm>,
    pub(super) row_scores: Option<output::RowScores>,
    pub(super) dia_replaced: HashSet<String>,
    pub(super) dia_precursor_ids: HashSet<String>,
    pub(super) pin_candidates: Vec<search::candidates::Candidate>,
    pub(super) pin_index: search::SearchIndex,
}

/// The steps after rescoring: the DIA study tables, the refined rows and IPF, protein groups and
/// cutoffs, the identification summary and statistics log, the open-search offsets, the knockoff
/// rows, the QPX bundle, then the last steps in `finish_outputs`.
pub(super) fn after_rescoring(state: AfterRescoring<'_>) -> Result<(), Box<dyn std::error::Error>> {
    let AfterRescoring {
        cli,
        spectrum_path,
        bench_mode,
        quant,
        run_ends,
        measurer,
        file_calibrations,
        spectra,
        queues,
        target_db,
        ipf_scoring,
        run_stats,
        fact_rows,
        mut dia_pass2,
        rerank_dropped,
        fdr_threshold,
        proteoform_fdr,
        use_native,
        mut mbr_decoys,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        params,
        restricted_site_masses,
        fixed_terminal_mods,
        scorer,
        is_mgf,
        mut rescore_map,
        mbr_source,
        rescore_decoys,
        mut rescore_decoy_psms,
        row_scores,
        dia_replaced,
        dia_precursor_ids,
        pin_candidates,
        pin_index,
    } = state;
    let spectrum_paths = &cli.spectrum;

    if let (Some(t), Some(dir)) = (dia_pass2.take(), cli.quant_dir.as_deref()) {
        crate::dia::study::store_run(
            dir,
            &spectrum_path,
            &output_pin_path,
            &t,
            None,
            &Default::default(),
        )?;
    }
    if cli.dia {
        crate::memory::release_free_memory();
    }

    let mut knockoff_rows: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut ipf_qpx: Option<std::collections::HashMap<String, output::ipf::IpfQpx>> = None;
    if let Some(ref mut map) = rescore_map {
        (knockoff_rows, ipf_qpx) = refine_rescored(
            &RescoredInputs {
                cli: &cli,
                pin_index: &pin_index,
                pin_candidates: &pin_candidates,
                rescore_decoy_psms: &rescore_decoy_psms,
                row_scores: &row_scores,
                queues: &queues,
                spectra: &spectra,
                restricted_site_masses: &restricted_site_masses,
                fixed_terminal_mods: &fixed_terminal_mods,
                report_base: &report_base,
                ipf_scoring: &ipf_scoring,
                run_ends: &run_ends,
                measurer: &measurer,
                fdr_threshold,
                proteoform_fdr,
            },
            map,
        )?;
    }
    let protein_groups = rescore_map.as_ref().map(|map| {
        infer_and_write_groups(
            map,
            &rescore_decoy_psms,
            &pin_index,
            &knockoff_rows,
            &report_base,
            (fdr_threshold, proteoform_fdr),
        )
    });
    if cli.mbr {
        mbr_decoys = std::mem::take(&mut rescore_decoy_psms);
    }
    let cutoffs = rescore_map
        .as_ref()
        .map(|m| output::Cutoffs::new(m, &pin_index, fdr_threshold, proteoform_fdr));

    let mut identified: Option<(Vec<&output::PercolatorPsm>, String)> = None;
    let mut open_search_summary: Option<String> = None;
    if let (Some(map), Some(cutoffs)) = (rescore_map.as_ref(), cutoffs.as_ref()) {
        let fdr = fdr_threshold;
        let (stem, dir) = report_stem(&report_base);
        let q_tag = format!("{fdr}").replace('.', "p");
        let filtered_path = dir.join(format!("{stem}.q{q_tag}.tsv"));
        let accepted = accepted_psms(&cli, map, cutoffs, &knockoff_rows, protein_groups.as_ref());
        let pep_note = cut_note(&cli, &pin_index, fdr, proteoform_fdr);
        match write_filtered_tsv(&filtered_path, &accepted, protein_groups.as_ref()) {
            Ok(()) => log::info!(
                target: "output",
                "Wrote filtered TSV: {} ({} PSMs at q<={}{})",
                filtered_path.display(),
                accepted.len(),
                fdr,
                pep_note
            ),
            Err(e) => {
                log::warn!(target: "output", "could not write {}: {e}", filtered_path.display())
            }
        }
        if cli.open_search {
            open_search_summary =
                write_open_offsets(&output_pin_path, dir, &stem, &accepted, &rescore_decoy_psms);
        }
        identified = Some((accepted, pep_note));
    }

    let mut id_summary = identified.as_ref().map(|(accepted, cut_note)| {
        let wanted: std::collections::HashSet<&str> =
            accepted.iter().map(|p| p.psm_id.as_str()).collect();
        let facts = match &fact_rows {
            Some(rows) => output::psm_facts_of_rows(
                &spectra,
                rows,
                &pin_candidates,
                &params,
                &run_ends,
                &file_calibrations,
                &wanted,
            ),
            None => output::psm_facts(
                &spectra,
                &queues,
                &pin_candidates,
                &params,
                &run_ends,
                &file_calibrations,
                &wanted,
            ),
        };
        let linked: Vec<output::PsmFacts> = accepted
            .iter()
            .filter_map(|p| facts.get(p.psm_id.as_str()).cloned())
            .collect();
        let decoys = if cli.mbr {
            &mbr_decoys
        } else {
            &rescore_decoy_psms
        };
        let decoys_at_cut = output::Cutoffs::new(decoys, &pin_index, fdr_threshold, proteoform_fdr)
            .count_passing(decoys, cli.pep);
        let run_names: Vec<String> = spectrum_paths
            .iter()
            .map(|p| crate::quant::run::run_stem(p))
            .collect();
        let fixed = output::FixedModSpec::declared(&params);
        output::IdentificationSummary::compute(
            &output::SpectraTally::of(&spectra, &run_ends),
            &run_names,
            &linked,
            accepted.len() - linked.len(),
            &output::SummaryContext {
                rescorer: if use_native {
                    "Percolator (native engine)"
                } else {
                    "Percolator"
                },
                cut: &format!("q ≤ {fdr_threshold}{cut_note}"),
                decoys_at_cut,
                protein_groups: protein_groups.as_ref().map(|inf| {
                    (
                        fdr_threshold,
                        inf.groups
                            .iter()
                            .filter(|g| g.q_value <= fdr_threshold)
                            .count(),
                    )
                }),
                fixed: &fixed,
                calibration_shift_ppm: match file_calibrations.len() {
                    0 | 1 => params.precursor_mass_shift_ppm,
                    n => file_calibrations.iter().map(|f| f.shift_ppm).sum::<f64>() / n as f64,
                },
            },
        )
    });
    if let Some(summary) = id_summary.as_mut() {
        summary.open_search = open_search_summary.take();
    }
    log::info!(target: "summary", "{}", run_stats.render_with(id_summary.as_ref()));
    write_statistics(&run_stats, id_summary.as_ref(), &report_base);

    knockoff_rows.extend(rerank_dropped);
    if !dia_precursor_ids.is_empty() {
        knockoff_rows.extend(dia_replaced);
        knockoff_rows.retain(|id| !dia_precursor_ids.contains(id));
    }
    if let Some(ref parquet_dir) = cli.output_parquet {
        write_qpx_bundle(
            parquet_dir,
            &PsmTables {
                spectra: &spectra,
                queues: &queues,
                candidates: &pin_candidates,
                params: &params,
                index: &pin_index,
            },
            &scorer.feature_match_tolerance(),
            (&spectrum_path, &original_spectrum_paths, &run_ends),
            (
                rescore_map.as_ref(),
                rescore_decoys.as_ref(),
                &knockoff_rows,
                ipf_qpx.as_ref(),
            ),
        )?;
    }

    finish_outputs(Finish {
        cli,
        bench_mode,
        quant,
        spectra,
        queues,
        target_db,
        mbr_decoys,
        report_base,
        params,
        is_mgf,
        rescore_map,
        mbr_source,
        rescore_decoys,
        pin_candidates,
        pin_index,
        knockoff_rows,
        protein_groups,
        cutoffs,
    })
}

/// The report's file stem (`andes` without one) and its directory.
fn report_stem(report_base: &Path) -> (String, &Path) {
    let stem = report_base
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
    (stem, dir)
}

/// The rescored PSMs that pass their q-value cutoff, `--pep` and `--protein-fdr` and are not
/// knockoff rows, in report order.
fn accepted_psms<'m>(
    cli: &Cli,
    map: &'m HashMap<String, output::PercolatorPsm>,
    cutoffs: &output::Cutoffs,
    knockoff_rows: &HashSet<String>,
    protein_groups: Option<&output::Inference>,
) -> Vec<&'m output::PercolatorPsm> {
    let pep_cap = cli.pep;
    let mut accepted: Vec<&'m output::PercolatorPsm> = map
        .iter()
        .filter(|(id, p)| {
            p.q_value <= cutoffs.of(id)
                && pep_cap.is_none_or(|t| p.pep <= t)
                && !knockoff_rows.contains(*id)
                && cli.protein_fdr.is_none_or(|t| {
                    protein_groups
                        .and_then(|inf| inf.group_of(p))
                        .is_some_and(|g| g.q_value <= t)
                })
        })
        .map(|(_, p)| p)
        .collect();
    accepted.sort_by(|a, b| a.report_cmp(b));
    accepted
}

/// The cut beyond the PSM q-value, as appended to the summary: PEP, proteoform and protein
/// group cutoffs when set.
fn cut_note(cli: &Cli, pin_index: &search::SearchIndex, fdr: f64, proteoform_fdr: f64) -> String {
    let pep_cap = cli.pep;
    let mut pep_note = pep_cap
        .map(|t| format!(" and PEP<={t}"))
        .unwrap_or_default();
    if pin_index.has_alternative_entries && proteoform_fdr != fdr {
        pep_note.push_str(&format!(", proteoform claims at q<={proteoform_fdr}"));
    }
    if let Some(t) = cli.protein_fdr {
        pep_note.push_str(&format!(", protein groups at q<={t}"));
    }
    pep_note
}

/// Writes `<stem>.statistics.log` beside the report.
fn write_statistics(
    run_stats: &output::RunStatistics,
    id_summary: Option<&output::IdentificationSummary>,
    report_base: &Path,
) {
    let stats_path = {
        let (stem, dir) = report_stem(report_base);
        dir.join(format!("{stem}.statistics.log"))
    };
    match output::write_statistics_log(run_stats, id_summary, &stats_path) {
        Ok(()) => log::info!(target: "output", "Wrote statistics: {}", stats_path.display()),
        Err(e) => log::warn!(target: "output", "could not write {}: {e}", stats_path.display()),
    }
}

/// What writing the PIN and rescoring it read and own; the scorer is borrowed, as the IPF scorer
/// borrows it.
pub(super) struct PinAndRescore<'a> {
    pub(super) cli: Cli,
    pub(super) t_total: std::time::Instant,
    pub(super) spectrum_path: PathBuf,
    pub(super) bench_mode: bool,
    pub(super) quant: Option<QuantCollector>,
    pub(super) predictors: Option<crate::dia::rt::Predictors>,
    pub(super) run_ends: Vec<usize>,
    pub(super) dia_ms1: Vec<Vec<model::scan::Ms1Scan>>,
    pub(super) measurer: FragmentErrorMeasurer,
    pub(super) file_calibrations: Vec<search::FileCalibration>,
    pub(super) spectra: Vec<Spectrum>,
    pub(super) queues: Vec<TopNQueue>,
    pub(super) target_db: model::protein::ProteinDb,
    pub(super) scratch_dir: PathBuf,
    pub(super) original_spectrum_paths: Vec<PathBuf>,
    pub(super) output_pin_path: PathBuf,
    pub(super) report_base: PathBuf,
    pub(super) params: search::SearchParams,
    pub(super) restricted_site_masses: Vec<f64>,
    pub(super) fixed_terminal_mods: Vec<(model::modification::ModLocation, f64)>,
    pub(super) scorer: &'a RankScorer,
    pub(super) is_mgf: bool,
    pub(super) t_phase: std::time::Instant,
    pub(super) refine_merged: bool,
    pub(super) ipf_scoring: Option<search::engine::IpfScoring<'a>>,
    pub(super) pin_candidates: Vec<search::candidates::Candidate>,
    pub(super) pin_index: search::SearchIndex,
}

/// Writes the PIN (with the RT and provenance features), the run statistics, the DIA predictor pass
/// and pooled re-rank, rescores it, then runs the steps after rescoring.
pub(super) fn write_pin_and_rescore(
    state: PinAndRescore<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    let PinAndRescore {
        cli,
        t_total,
        spectrum_path,
        bench_mode,
        quant,
        predictors,
        run_ends,
        dia_ms1,
        measurer,
        file_calibrations,
        mut spectra,
        mut queues,
        target_db,
        scratch_dir,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        params,
        restricted_site_masses,
        fixed_terminal_mods,
        scorer,
        is_mgf,
        t_phase,
        refine_merged,
        ipf_scoring,
        pin_candidates,
        pin_index,
    } = state;
    output::populate_rt_features_per_run(&spectra, &mut queues, &pin_candidates, &run_ends);

    output::populate_provenance_features(&mut queues, &pin_candidates, &pin_index);

    output::write_pin(
        &output_pin_path,
        &spectra,
        &queues,
        &pin_candidates,
        &params,
        &pin_index,
    )?;
    log::info!(
        target: "output",
        "Wrote PIN: {} in {:.2} s",
        output_pin_path.display(),
        t_phase.elapsed().as_secs_f64()
    );
    crate::memory::mark_stage("PIN written");
    let run_stats = output::RunStatistics::compute(
        &queues,
        &pin_candidates,
        &params,
        &scorer.feature_match_tolerance(),
        &scorer.pin_match_tolerance(),
        &scorer.node_score_window(),
    );
    let has_site_claims = pin_index.iter_target_proteins().any(|p| {
        p.annotation
            .as_deref()
            .is_some_and(|a| !a.site_claims.is_empty())
    });
    let fact_rows = (cli.dia
        && !has_site_claims
        && ipf_scoring.is_none()
        && quant.is_none()
        && cli.output_parquet.is_none()
        && cli.output_tsv.is_none())
    .then(|| output::fact_rows(&spectra, &queues, &pin_candidates));
    if fact_rows.is_some() {
        std::mem::take(&mut queues).into_par_iter().for_each(drop);
        crate::memory::release_free_memory();
    }
    log::info!(
        target: "total",
        "done in {:.2} s",
        t_total.elapsed().as_secs_f64()
    );
    let pooled = cli.dia && params.top_n_psms_per_spectrum > cli.top_n;
    let n_dia_ms1: usize = dia_ms1.iter().map(Vec::len).sum();
    let dia_ms1 = if cli.dia {
        crate::dia::ms1_indexes(dia_ms1)
    } else {
        Vec::new()
    };
    let mut dia_pass2: Option<output::dia::pass2::Pass2Table> = None;
    if let Some(p) = &predictors {
        let kept = p.add_to_pin(&output_pin_path, &spectra, &run_ends, pooled)?;
        crate::memory::mark_stage("predictions");
        if let Some((predicted, table)) = kept.filter(|_| pooled) {
            let pool = crate::dia::llr::add_to_pin(
                &output_pin_path,
                &spectra,
                &run_ends,
                &dia_ms1,
                table,
                &predicted,
                p.tolerance(),
            )?;
            crate::memory::mark_stage("pool scores");
            let nominees_out = cli
                .quant_dir
                .as_deref()
                .filter(|_| run_ends.len() == 1)
                .map(|d| {
                    crate::dia::transfer::nominee_path(
                        d,
                        &crate::quant::run::run_stem(&spectrum_path),
                    )
                });
            dia_pass2 = Some(crate::dia::pass2::run(
                &output_pin_path,
                &spectra,
                &run_ends,
                &dia_ms1,
                pool,
                predicted,
                nominees_out
                    .as_deref()
                    .map(|path| crate::dia::transfer::NomineesOut {
                        path,
                        spectrum: spectrum_path.as_path(),
                        mobility_slices: cli.dia_mobility_slices,
                    }),
            )?);
            crate::memory::mark_stage("pass 2");
        }
    }
    let early_learned = rescoring::start_learned(&cli, &mut dia_pass2);
    if cli.dia {
        crate::memory::release_free_memory();
    }
    crate::memory::mark_stage("pool released");
    let mut rerank_dropped = std::collections::HashSet::new();
    if pooled {
        rerank_dropped = crate::dia::rerank_pin(
            &output_pin_path,
            cli.top_n as usize,
            cli.output_parquet.is_some() || quant.is_some(),
        )?;
        crate::memory::release_free_memory();
        crate::memory::mark_stage("per-scan cut");
    }
    if cli.dia {
        crate::dia::add_to_pin(&output_pin_path, &spectra, &run_ends, &dia_ms1, n_dia_ms1)?;
        crate::memory::mark_stage("co-elution");
    }
    drop(dia_ms1);
    if cli.dia {
        crate::memory::release_free_memory();
    }
    if cli.refine {
        if refine_merged {
            log::info!("Refinement PSMs merged into unified PIN (Pass-1 ⊕ Pass-2).");
        } else {
            log::info!("Refinement produced no Pass-2 PSMs; unified PIN contains Pass-1 only.");
        }
    }

    let fdr_threshold = cli.psm_fdr();
    let proteoform_fdr = cli.proteoform_fdr();
    let use_native = !cli.external_percolator();
    let use_percolator = cli.rescore;
    if !use_percolator && (cli.psm_fdr.is_some() || cli.pep.is_some()) {
        let which = match (cli.psm_fdr.is_some(), cli.pep.is_some()) {
            (true, true) => "--psm-fdr and --pep",
            (true, false) => "--psm-fdr",
            _ => "--pep",
        };
        log::warn!(
            "{which} set without --rescore, so no rescoring runs and no filtered .tsv \
             is written. andes does not compute FDR on its own; the PIN is written for \
             Percolator. Add --rescore to rescore with the native Percolator engine."
        );
    }
    let mbr_decoys: std::collections::HashMap<String, output::PercolatorPsm> =
        std::collections::HashMap::new();
    let rescoring = rescoring::rescore_pin(
        &rescoring::RescoreInputs {
            cli: &cli,
            output_pin_path: &output_pin_path,
            spectrum_path: &spectrum_path,
            scratch_dir: &scratch_dir,
            use_native,
            use_percolator,
            free_peaks: fact_rows.is_some(),
            has_site_claims,
        },
        &mut spectra,
        &mut dia_pass2,
        early_learned,
    )?;
    let rescoring::Rescoring {
        rescore_map,
        mbr_source,
        rescore_decoys,
        rescore_decoy_psms,
        row_scores,
        dia_replaced,
        dia_precursor_ids,
    } = rescoring;
    after_rescoring(AfterRescoring {
        cli,
        spectrum_path,
        bench_mode,
        quant,
        run_ends,
        measurer,
        file_calibrations,
        spectra,
        queues,
        target_db,
        ipf_scoring,
        run_stats,
        fact_rows,
        dia_pass2,
        rerank_dropped,
        fdr_threshold,
        proteoform_fdr,
        use_native,
        mbr_decoys,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        params,
        restricted_site_masses,
        fixed_terminal_mods,
        scorer,
        is_mgf,
        rescore_map,
        mbr_source,
        rescore_decoys,
        rescore_decoy_psms,
        row_scores,
        dia_replaced,
        dia_precursor_ids,
        pin_candidates,
        pin_index,
    })
}

/// The spectra, their PSM queues and the refine output, as the search leaves them.
pub(super) type Searched = (
    Vec<Spectrum>,
    Vec<TopNQueue>,
    Option<search::refinement::RefinementOutput>,
);

/// What the steps after the search read.
#[derive(Clone, Copy)]
pub(super) struct PostSearch<'a> {
    pub(super) t_phase: std::time::Instant,
    pub(super) xl_snapped: usize,
    pub(super) cli: &'a Cli,
    pub(super) acquisition: &'a crate::acquisition::Acquisition,
    pub(super) decoy_strategy: search::decoy::DecoyStrategy,
    pub(super) fragment_tol_da: f64,
    pub(super) mono_active: bool,
    pub(super) mono_params: &'a search::precursor::monoisotopic::MonoParams,
    pub(super) routed_activation: ActivationMethod,
    pub(super) spectrum_paths: &'a Vec<PathBuf>,
    pub(super) standalone: bool,
    pub(super) t_search_start: std::time::Instant,
    pub(super) xl_snap_files: usize,
    pub(super) bench_mode: bool,
    pub(super) file_calibrations: &'a Vec<search::FileCalibration>,
    pub(super) idx: &'a search::SearchIndex,
    pub(super) output_pin_path: &'a PathBuf,
    pub(super) params: &'a search::SearchParams,
    pub(super) quant: &'a Option<QuantCollector>,
    pub(super) report_base: &'a PathBuf,
    pub(super) run_ends: &'a Vec<usize>,
    pub(super) scorer: &'a RankScorer,
    pub(super) spectrum_path: &'a PathBuf,
    pub(super) t_total: std::time::Instant,
}

/// What the steps after the search consume.
pub(super) struct PostSearchOwned {
    pub(super) all_spectra: Vec<Spectrum>,
    pub(super) all_queues: Vec<TopNQueue>,
    pub(super) mono_table: Vec<Option<MonoCorrection>>,
    pub(super) mono_dump: Option<MonoDump>,
    pub(super) mono_stats: MonoStats,
    pub(super) xl_ms3: Vec<Vec<model::scan::ProductScan>>,
}

/// The steps between the search and the PIN: the search and routing log, the monoisotopic
/// correction, the standalone and DIA logs, open-search and cross-link summaries, the glyco pass and
/// the refine pass. Returns the spectra, their PSM queues and the refine output, or none when
/// the run ends there (a bench run or a glyco or cross-link search).
pub(super) fn post_search(
    inputs: &PostSearch<'_>,
    owned: PostSearchOwned,
    prepared: &mut PreparedSearch<'_>,
) -> Result<Option<Searched>, Box<dyn std::error::Error>> {
    let PostSearch {
        t_phase,
        xl_snapped,
        cli,
        acquisition,
        decoy_strategy,
        fragment_tol_da,
        mono_active,
        mono_params,
        routed_activation,
        spectrum_paths,
        standalone,
        t_search_start,
        xl_snap_files,
        bench_mode,
        file_calibrations,
        idx,
        output_pin_path,
        params,
        quant,
        report_base,
        run_ends,
        scorer,
        spectrum_path,
        t_total,
    } = *inputs;
    let PostSearchOwned {
        all_spectra,
        all_queues,
        mut mono_table,
        mut mono_dump,
        mono_stats,
        xl_ms3,
    } = owned;
    let search_elapsed = t_search_start.elapsed();
    log::info!(
        target: "stream_search",
        "Loaded and scored {} spectra from {} in chunks of {}, done in {:.2} s",
        all_spectra.len(),
        spectrum_path.display(),
        CHUNK_SIZE,
        t_phase.elapsed().as_secs_f64()
    );
    log::info!(
        target: "chunk_candidate_prep",
        "done in {:.2} s (part of stream_search: fragment-index build or record expansion, summed over chunks)",
        prepared.chunk_prep_seconds()
    );
    for (activation, own_model, scans) in prepared.per_scan_routes() {
        log::info!(
            target: "routing",
            "{scans} {} scan scorings with {} (run routed {})",
            activation.name(),
            if own_model {
                "the model of their own activation"
            } else {
                "the run's model"
            },
            routed_activation.name()
        );
    }
    if bench_mode {
        log::info!(target: "bench", "capped at {} spectra", cli.max_spectra);
    }

    prepared.sync_materialized_candidates();

    let mono_columns = mono_active;
    let (mono_active, glyco_isotope_override) = finish_mono_correction(
        cli,
        params,
        mono_params,
        mono_active,
        &mono_stats,
        &mut mono_dump,
        &mut mono_table,
    )?;
    if mono_active {
        debug_assert_eq!(mono_table.len(), all_spectra.len());
    }

    let spectra = all_spectra;
    let queues = all_queues;

    if standalone {
        log::info!(
            "Standard peptide search skipped in --{} mode (PSMs unused; the driver \
             derives its own candidates); {} spectra loaded ({:.2}s)",
            if cli.glyco { "glyco" } else { "xl" },
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    } else {
        let non_empty = queues.iter().filter(|q| !q.is_empty()).count();
        log::info!(
            "Search complete: {non_empty} / {} spectra have PSMs (search wall: {:.2}s)",
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    }
    if cli.dia {
        crate::memory::release_free_memory();
    }
    crate::memory::mark_stage("search");
    if cli.open_search {
        let (attempts, wins) = prepared.open_search_stats();
        log::info!(
            target: "open-search",
            "{attempts} modification-aware re-scores ({wins} improved the node \
             score)"
        );
    }

    if cli.xl.is_some() {
        if xl_snap_files > 0 {
            log::info!(
                target: "xl",
                "{xl_snapped} precursor m/z values moved to the dominant MS1 peak within {} ppm",
                search::crosslinks::PRECURSOR_SNAP_PPM
            );
        }
        crate::xl::run_xl(
            cli,
            &spectra,
            &xl_ms3,
            prepared,
            report_base,
            spectrum_paths,
            run_ends,
            t_total,
        )?;
        return Ok(None);
    }

    if cli.glyco {
        let unfitted: Vec<Option<MonoCorrection>>;
        let mono_for_pin = if mono_active {
            Some(mono_table.as_slice())
        } else if mono_columns {
            unfitted = vec![None; spectra.len()];
            Some(unfitted.as_slice())
        } else {
            None
        };
        run_glyco(
            cli,
            &spectra,
            prepared,
            params,
            idx,
            output_pin_path,
            spectrum_paths,
            acquisition.high_res,
            mono_for_pin,
            glyco_isotope_override,
            quant.as_ref(),
            t_total,
        )?;
        return Ok(None);
    }

    let refine_output = if cli.refine {
        let base_cfg = match &cli.refine_config {
            Some(p) => search::RefineConfig::from_yaml_str(&std::fs::read_to_string(p)?)
                .map_err(|e| format!("parsing --refine-config {}: {e}", p.display()))?,
            None => search::RefineConfig::default_tier(),
        };
        let cfg = base_cfg;

        let high_res = acquisition.high_res;

        let high_res_only = cfg.high_res_only;
        if high_res_only && !high_res {
            log::warn!("refine is high-res-only and the data is low-res; skipping refinement.");
            None
        } else {
            let target_db = idx.target_db();
            let t_refine = std::time::Instant::now();
            let out = search::refinement::run_refinement(
                &queues,
                &spectra,
                &prepared.candidates,
                &target_db,
                params,
                scorer,
                &cfg,
                0.01,
                high_res,
                fragment_tol_da,
                &cli.db.decoy_prefix,
                decoy_strategy,
                cli.decoy_seed,
                file_calibrations,
            );
            log::info!(
                target: "refinement",
                "done in {:.2} s",
                t_refine.elapsed().as_secs_f64()
            );
            out
        }
    } else {
        None
    };

    Ok(Some((spectra, queues, refine_output)))
}
