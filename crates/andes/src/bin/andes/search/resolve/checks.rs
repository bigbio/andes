//! Checks of the search flags before anything runs: incompatible combinations, MBR arguments,
//! and warnings for deuterated tags, restricted artefact masses and other risky settings.

use super::set_flags;
use crate::cli::Cli;

/// Rejects flag combinations that cannot work, using the flags alone.
pub(crate) fn check_flag_combinations(cli: &Cli) -> Result<(), String> {
    if cli.top_n == 0 {
        return Err("--top-n must be at least 1".into());
    }
    if cli.dia_mobility_slices == 0 {
        return Err("--dia-mobility-slices must be at least 1".into());
    }
    if !cli.dia_ms1_tol_ppm.is_finite() || cli.dia_ms1_tol_ppm <= 0.0 {
        return Err("--dia-ms1-tol-ppm must be a positive number".into());
    }
    if cli.localize.is_some() && !cli.rescore {
        return Err("--localize ipf needs --rescore".into());
    }
    if cli.xl.is_none() && !cli.external_percolator() {
        if !cli.percolator_args.trim().is_empty() {
            return Err(
                "native rescoring takes no percolator_args; add --percolator-bin or --percolator-docker to run an external Percolator"
                    .into(),
            );
        }
        if cli.percolator_image != output::DEFAULT_PERCOLATOR_IMAGE {
            return Err(
                "--percolator-image has no effect without --percolator-docker; add --percolator-docker to run Percolator from that image"
                    .into(),
            );
        }
    }
    if cli.dia {
        if cli.chimeric {
            return Err("--dia and --chimeric cannot be combined: --dia searches every scan across its whole window".into());
        }
        if cli.open_search {
            return Err("--dia and --open-search cannot be combined: a DIA scan has no selected precursor to offset from".into());
        }
        if cli.refine {
            return Err(
                "--dia and --refine cannot be combined: the refinement pass has no MS1 gate".into(),
            );
        }
        if cli.glyco {
            return Err("--dia and --glyco cannot be combined: the glycopeptide search has no isolation-window or MS1 gate".into());
        }
        if cli.xl.is_some() {
            return Err("--dia and --xl cannot be combined: the cross-link search has no isolation-window or MS1 gate".into());
        }
    }
    if cli.add_entrapment {
        let unsupported = set_flags(&[
            (cli.glyco, "--glyco"),
            (cli.xl.is_some(), "--xl"),
            (cli.open_search, "--open-search"),
        ]);
        if !unsupported.is_empty() {
            return Err(format!(
                "--add-entrapment is not combined with {}: these modes write their own results \
                 and FDR, which the entrapment estimate does not read",
                unsupported.join(", ")
            ));
        }
        if !cli.rescore {
            log::warn!(target: "entrapment",
                "--add-entrapment without --rescore: the entrapment proteins are searched and \
                 written to the PIN, but no FDP is estimated (it needs the rescored q-values)"
            );
        }
    }
    if cli.glyco_index_sequon_only && !cli.glyco {
        return Err(
            "--glyco-index-sequon-only requires --glyco (the standard peptide \
                    search scores every candidate)"
                .into(),
        );
    }
    if cli.open_search {
        if cli.glyco {
            return Err("--open-search is not supported with --glyco".into());
        }
        if cli.chimeric {
            return Err("--open-search is not combined with --chimeric".into());
        }
        if cli.label.is_some() {
            return Err("--open-search is not combined with --label".into());
        }
        if cli.tmt.is_some() {
            return Err("--open-search is not combined with --tmt".into());
        }
        if cli.lfq || cli.mbr {
            return Err("--open-search is not combined with --lfq or --mbr".into());
        }
        if cli.refine {
            return Err("--open-search is not combined with --refine".into());
        }
        if cli.xl.is_some() {
            return Err("--open-search is not combined with --xl".into());
        }
        if cli.isotope_error.is_some() {
            return Err(
                "--open-search chooses the isotope offset from the fragments (range 0..=2); \
                 remove --isotope-error"
                    .into(),
            );
        }
    }
    if cli.xl.is_some() {
        let unsupported = set_flags(&[
            (cli.output_tsv.is_some(), "--output-tsv"),
            (cli.output_parquet.is_some(), "--output-parquet"),
            (cli.refine, "--refine"),
            (cli.psm_fdr.is_some(), "--psm-fdr"),
            (cli.pep.is_some(), "--pep"),
            (cli.lfq, "--lfq"),
            (cli.mbr, "--mbr"),
            (cli.quant_dir.is_some(), "--quant-dir"),
            (cli.chimeric, "--chimeric"),
            (cli.localize.is_some(), "--localize"),
            (cli.tmt.is_some(), "--tmt"),
            (cli.label.is_some(), "--label"),
        ]);
        if !unsupported.is_empty() {
            return Err(format!(
                "--xl does not support: {} (crosslink mode writes its own .xl.pin and \
                 tables, scored by the pair score)",
                unsupported.join(", ")
            ));
        }
    }
    if cli.glyco {
        let unsupported = set_flags(&[
            (cli.output_tsv.is_some(), "--output-tsv"),
            (cli.output_parquet.is_some(), "--output-parquet"),
            (cli.rescore, "--rescore"),
            (cli.refine, "--refine"),
            (cli.psm_fdr.is_some(), "--psm-fdr"),
            (cli.pep.is_some(), "--pep"),
        ]);
        if !unsupported.is_empty() {
            return Err(format!(
                "--glyco does not support: {} (glyco mode writes a standalone \
                 .glyco.pin and skips the standard PIN/rescore/TSV/Parquet/refine \
                 pipeline; run Percolator on the .glyco.pin separately)",
                unsupported.join(", ")
            ));
        }
        if cli.mbr && !cli.lfq {
            return Err(
                "--glyco --mbr needs --lfq: match-between-runs transfers the identified intact \
                 glycopeptides into the other runs"
                    .into(),
            );
        }
        if cli.optimize {
            return Err("--optimize is not supported with --glyco".into());
        }
    }
    Ok(())
}

/// Logs a pre-formatted warning under the `search` target without its `WARN: ` prefix.
pub(super) fn warn_search(warning: &str) {
    log::warn!(target: "search", "{}", warning.strip_prefix("WARN: ").unwrap_or(warning));
}

const ARTEFACT_MASSES: [f64; 4] = [0.984016, 15.994915, 43.005814, 14.015650];

pub(super) fn warn_restricted_artefacts(aa: &model::AminoAcidSet, variants_on: bool) {
    if !variants_on {
        return;
    }
    for m in aa.distinct_mods() {
        if m.annotated_site_mass.is_some()
            && ARTEFACT_MASSES
                .iter()
                .any(|a| (m.mass_delta - a).abs() <= model::protein::SITE_MASS_TOLERANCE)
        {
            log::warn!(target: "search",
                "{} ({:+.4} Da) is restricted to annotated sites while \
                 --uniprot-variants is on: variants of the same mass lose the \
                 competing explanation at unannotated residues.",
                m.name,
                m.mass_delta
            );
        }
    }
}

pub(super) fn check_mbr_args(cli: &Cli) -> Result<(), String> {
    if cli.spectrum.len() < 2 {
        return Err(format!(
            "--mbr transfers identifications between runs and needs two or more --spectrum inputs (got {})",
            cli.spectrum.len()
        ));
    }
    if !cli.lfq && cli.tmt.is_none() && cli.label.is_none() {
        return Err(
            "--mbr needs --lfq (label-free), --tmt (reporter ions) or --label (labelled pairs)"
                .into(),
        );
    }
    if !cli.rescore {
        return Err(
            "--mbr needs --rescore: the library ranks precursors by the rescored PSM scores".into(),
        );
    }
    if cli.tmt.is_some() {
        if cli.lfq {
            return Err(
                "--mbr transfers either label-free MS1 quantities (--lfq) or reporter ions \
                 (--tmt), not both; run them separately"
                    .into(),
            );
        }
        if cli.tmt_level.is_some_and(|level| level >= 3) {
            return Err(
                "--mbr with --tmt needs the reporter ions in the MS2 scans and does not apply \
                 to SPS-MS3 reporters (--tmt-level 3)"
                    .into(),
            );
        }
        return Ok(());
    }
    crate::quant::evidence_params(cli.lfq_tol)?;
    Ok(())
}

pub(super) fn deuterated_tag_warning(
    plex: Option<quant::Plex>,
    precursor_tol: model::Tolerance,
) -> Option<String> {
    use model::Tolerance;
    let plex = plex.filter(|p| p.has_deuterated_channels())?;
    let shift = 3.0 * quant::isobaric::DEUTERATED_TAG_SHIFT_DA;
    if precursor_tol.as_da(1000.0) >= shift {
        return None;
    }
    let tol = match precursor_tol {
        Tolerance::Ppm(v) => format!("{v}ppm"),
        Tolerance::Da(v) => format!("{v}da"),
    };
    Some(format!(
        "WARN: --tmt {}: each deuterated TMTpro tag is {:.1} mDa heavier than the TMTpro mass \
         (304.2071 Da), so a 1000 Da peptide with three of them sits {:.1} ppm above its \
         candidate, outside --precursor-tol {tol}. Widen it to 10ppm or more, or the \
         deuterated channels' precursors are missed.",
        plex.name(),
        quant::isobaric::DEUTERATED_TAG_SHIFT_DA * 1e3,
        shift / 1000.0 * 1e6,
    ))
}
