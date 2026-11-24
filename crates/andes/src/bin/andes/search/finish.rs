//! The last steps of a search: quantification, the proteoform table, the MBR job, the search TSV
//! and the release of the search's memory before the in-process MBR run.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use model::Spectrum;
use search::TopNQueue;

use super::outputs::{release_at_exit, write_proteoform_table, write_search_tsv, PsmTables};
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
