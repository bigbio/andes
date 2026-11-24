//! The last steps of a search: quantification, the proteoform table, the MBR job, the search TSV
//! and the release of the search's memory before the in-process MBR run.

use std::collections::{HashMap, HashSet};

use std::path::PathBuf;

use model::Spectrum;
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
        crate::dia::study::store_run(dir, &spectrum_path, &output_pin_path, &t, None)?;
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
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        let q_tag = format!("{fdr}").replace('.', "p");
        let filtered_path = dir.join(format!("{stem}.q{q_tag}.tsv"));
        let pep_cap = cli.pep;
        let mut accepted: Vec<&output::PercolatorPsm> = map
            .iter()
            .filter(|(id, p)| {
                p.q_value <= cutoffs.of(id)
                    && pep_cap.is_none_or(|t| p.pep <= t)
                    && !knockoff_rows.contains(*id)
                    && cli.protein_fdr.is_none_or(|t| {
                        protein_groups
                            .as_ref()
                            .and_then(|inf| inf.group_of(p))
                            .is_some_and(|g| g.q_value <= t)
                    })
            })
            .map(|(_, p)| p)
            .collect();
        accepted.sort_by(|a, b| a.report_cmp(b));
        let mut pep_note = pep_cap
            .map(|t| format!(" and PEP<={t}"))
            .unwrap_or_default();
        if pin_index.has_alternative_entries && proteoform_fdr != fdr {
            pep_note.push_str(&format!(", proteoform claims at q<={proteoform_fdr}"));
        }
        if let Some(t) = cli.protein_fdr {
            pep_note.push_str(&format!(", protein groups at q<={t}"));
        }
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
    let stats_path = {
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        dir.join(format!("{stem}.statistics.log"))
    };
    match output::write_statistics_log(&run_stats, id_summary.as_ref(), &stats_path) {
        Ok(()) => log::info!(target: "output", "Wrote statistics: {}", stats_path.display()),
        Err(e) => log::warn!(target: "output", "could not write {}: {e}", stats_path.display()),
    }

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
