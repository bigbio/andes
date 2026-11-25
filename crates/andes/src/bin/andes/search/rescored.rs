//! Post-processing of the rescored PSMs: FDR within each variant class, knockoff q-values, site
//! claims and IPF localization.

use std::path::Path;

use model::Spectrum;
use search::TopNQueue;

use super::calibrate::FragmentErrorMeasurer;
use super::RescoredOutputs;
use crate::cli::Cli;

/// What the post-processing of the rescored PSMs reads: the flags, the PIN candidates and their
/// index, the decoy PSMs and row scores, the searched spectra and queues, the site-claim inputs, the
/// IPF scorer and the per-file fragment errors.
#[derive(Clone, Copy)]
pub(super) struct RescoredInputs<'a, 'b> {
    pub(super) cli: &'a Cli,
    pub(super) pin_index: &'a search::SearchIndex,
    pub(super) pin_candidates: &'a [search::candidates::Candidate],
    pub(super) rescore_decoy_psms: &'a std::collections::HashMap<String, output::PercolatorPsm>,
    pub(super) row_scores: &'a Option<output::RowScores>,
    pub(super) queues: &'a [TopNQueue],
    pub(super) spectra: &'a [Spectrum],
    pub(super) restricted_site_masses: &'a [f64],
    pub(super) fixed_terminal_mods: &'a [(model::ModLocation, f64)],
    pub(super) report_base: &'a Path,
    pub(super) ipf_scoring: &'a Option<search::engine::IpfScoring<'b>>,
    pub(super) run_ends: &'a [usize],
    pub(super) measurer: &'a FragmentErrorMeasurer,
    pub(super) fdr_threshold: f64,
    pub(super) proteoform_fdr: f64,
}

/// FDR within each variant class, knockoff q-values, site claims (written to `.sites.tsv`) and, with
/// `--localize ipf`, site localization (`.ipf.tsv`) of the rescored PSMs in `map`. Returns the
/// knockoff rows and the IPF table for the QPX output.
pub(super) fn refine_rescored(
    inputs: &RescoredInputs<'_, '_>,
    map: &mut std::collections::HashMap<String, output::PercolatorPsm>,
) -> Result<RescoredOutputs, Box<dyn std::error::Error>> {
    let RescoredInputs {
        cli,
        pin_index,
        pin_candidates,
        rescore_decoy_psms,
        row_scores,
        queues,
        spectra,
        restricted_site_masses,
        fixed_terminal_mods,
        report_base,
        ipf_scoring,
        run_ends,
        measurer,
        fdr_threshold,
        proteoform_fdr,
    } = *inputs;
    let mut ipf_qpx: Option<std::collections::HashMap<String, output::ipf::IpfQpx>> = None;
    if pin_index.has_alternative_entries {
        let shared =
            output::share_isobaric_references(map, pin_index, fdr_threshold.max(proteoform_fdr));
        if shared > 0 {
            log::info!(
                target: "rescore",
                "{shared} variant-class PSMs share their peptide with a reference \
                 protein once I and L are taken as equal; they now name it"
            );
        }
    }
    let counts = output::apply_class_fdr(map, rescore_decoy_psms, pin_index);
    if counts.len() > 1 {
        let line: Vec<String> = counts
            .iter()
            .map(|c| {
                format!(
                    "{} {}/{} {} pi0 {:.2}",
                    c.class,
                    c.targets,
                    c.decoys,
                    c.method.name(),
                    c.pi0
                )
            })
            .collect();
        log::info!(
            target: "rescore",
            "FDR within each class (targets/decoys, method, pi0): {}",
            line.join(", ")
        );
    }
    let outcome = output::apply_knockoff_fdr(map, pin_index, proteoform_fdr, cli.db.knockoffs > 0);
    let mut knockoff_rows = outcome.knockoff_rows;
    let untested: usize = outcome
        .counts
        .iter()
        .map(|c| c.real_before - c.tested_before)
        .sum();
    if untested > 0 {
        log::warn!(
            "{untested} variant PSMs at q <= {proteoform_fdr} have no knockoff \
             (the enzyme or the variant kind leaves no position for one) and are not \
             accepted; --knockoffs 0 leaves variants to the class FDR alone"
        );
    }
    if !outcome.counts.is_empty() {
        let line: Vec<String> = outcome
            .counts
            .iter()
            .map(|c| {
                format!(
                    "{} {} -> {} ({} tested; knockoffs {:.1}, false share {:.2})",
                    c.stratum,
                    c.real_before,
                    c.real_after,
                    c.tested_before,
                    c.knockoff_weighted,
                    c.false_share()
                )
            })
            .collect();
        log::info!(
            target: "rescore",
            "variant PSMs at q<={proteoform_fdr} per substitution class, before -> \
             after knockoff q-values: {}",
            line.join(", ")
        );
    }
    let posteriors = row_scores.as_ref().map(|(target_scores, decoy_scores)| {
        let rows: std::collections::HashMap<String, output::fdr::identity::RowEvidence> = queues
            .iter()
            .enumerate()
            .filter(|(_, q)| q.iter_psms().any(|p| p.site_evidence.is_some()))
            .flat_map(|(i, q)| {
                output::emitted_row_evidence(&spectra[i], q, pin_candidates)
                    .into_iter()
                    .map(move |(id, isotope_offset, evidence)| {
                        (
                            id,
                            output::fdr::identity::RowEvidence {
                                spectrum: i as u32,
                                isotope_offset,
                                evidence,
                            },
                        )
                    })
            })
            .collect();
        output::fdr::identity::row_posteriors(target_scores, decoy_scores, &rows)
    });
    let sites = output::site_claims(
        map,
        pin_index,
        fdr_threshold,
        proteoform_fdr,
        restricted_site_masses,
        posteriors.as_ref(),
        fixed_terminal_mods,
    );
    knockoff_rows.extend(sites.knockoff_rows.iter().cloned());
    if !sites.counts.is_empty() {
        let line: Vec<String> = sites
            .counts
            .iter()
            .map(|c| {
                format!(
                    "{} {} -> {} ({} tested; knockoffs {:.1})",
                    c.stratum, c.claims, c.accepted, c.tested, c.knockoff_weighted
                )
            })
            .collect();
        log::info!(
            target: "rescore",
            "site-claim PSMs at q<={fdr_threshold} per modification and residue, \
             before -> after site q<={proteoform_fdr}: {}",
            line.join(", ")
        );
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        let path = dir.join(format!("{stem}.sites.tsv"));
        match output::write_sites_tsv(&path, &sites.rows) {
            Ok(()) => log::info!(
                target: "output",
                "Wrote site claims: {} ({} sites)",
                path.display(),
                sites.rows.len()
            ),
            Err(e) => log::warn!(target: "output", "could not write {}: {e}", path.display()),
        }
    }
    if let Some(scoring) = ipf_scoring.as_ref() {
        use rayon::prelude::*;
        let t_ipf = std::time::Instant::now();
        let map_ref: &std::collections::HashMap<String, output::percolator::PercolatorPsm> = map;
        let cutoffs = output::Cutoffs::new(map_ref, pin_index, fdr_threshold, proteoform_fdr);
        let searched_with = scoring.fragment_error();
        let mut evidence: std::collections::HashMap<String, search::ipf::IpfEvidence> =
            std::collections::HashMap::new();
        let mut lo = 0usize;
        for (f, &hi) in run_ends.iter().enumerate() {
            scoring.set_fragment_error(measurer.per_file.get(f).copied().flatten());
            let part: Vec<(String, search::ipf::IpfEvidence)> = queues[lo..hi]
                .par_iter()
                .enumerate()
                .map(|(j, q)| (lo + j, q))
                .flat_map_iter(|(i, q)| {
                    let rows: Vec<(String, u32, u8)> =
                        output::emitted_row_candidates(&spectra[i], q, pin_candidates)
                            .into_iter()
                            .filter(|(id, c, _)| {
                                !pin_candidates[*c as usize].is_decoy
                                    && !knockoff_rows.contains(id)
                                    && map_ref.get(id).is_some_and(|p| p.q_value <= cutoffs.of(id))
                            })
                            .collect();
                    let peptides: Vec<(&model::peptide::Peptide, u8)> = rows
                        .iter()
                        .map(|(_, c, z)| (&pin_candidates[*c as usize].peptide, *z))
                        .collect();
                    let found = if peptides.is_empty() {
                        Vec::new()
                    } else {
                        scoring.evidence(&spectra[i], &peptides)
                    };
                    rows.into_iter()
                        .zip(found)
                        .filter_map(|((id, _, _), e)| e.map(|e| (id, e)))
                        .collect::<Vec<_>>()
                })
                .collect();
            evidence.extend(part);
            lo = hi;
        }
        scoring.set_fragment_error(searched_with);
        let (rows, model) = output::ipf::localize_psms(map_ref, &evidence, cli.ipf_debug);
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        let path = dir.join(format!("{stem}.ipf.tsv"));
        output::ipf::write_ipf_tsv(&path, &rows)
            .map_err(|e| format!("could not write {}: {e}", path.display()))?;
        if rows.is_empty() {
            log::warn!(
                "--localize ipf: no identified target PSM has a variable \
                 modification with more than one possible site; {} is empty",
                path.display()
            );
        }
        let truncated = rows.iter().filter(|r| r.truncated).count();
        log::info!(
            target: "output",
            "Wrote IPF localization: {} ({} PSMs, {} with placements limited to the \
             {} nearest the search's; ion PEP by intensity quartile {:?}, unmatched \
             {:.3}) [{:.2}s]",
            path.display(),
            rows.len(),
            truncated,
            search::ipf::MAX_PLACEMENTS,
            model.matched.map(|p| (p * 1000.0).round() / 1000.0),
            model.unmatched,
            t_ipf.elapsed().as_secs_f64()
        );
        if cli.output_parquet.is_some() {
            ipf_qpx = Some(output::ipf::qpx_map(&rows));
        }
    }
    Ok((knockoff_rows, ipf_qpx))
}
