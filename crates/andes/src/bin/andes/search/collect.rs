//! The quantification collector of a search and the checks that its inputs can feed it.

use std::path::PathBuf;

use scoring_crate::RankScorer;

use crate::cli::Cli;
use crate::quant::run::{QuantCollector, QuantSettings};

/// The quantification collector of the run (`--lfq`, `--tmt`, `--label`, `--quant-dir`), after the
/// checks that the inputs and the other flags can feed it; `None` when no quantification is asked.
pub(super) fn quant_collector(
    cli: &Cli,
    high_res: bool,
    tmt_level: u8,
    scorer: &RankScorer,
    is_mgf: bool,
    is_d: bool,
    spectrum_paths: &[PathBuf],
) -> Result<Option<QuantCollector>, Box<dyn std::error::Error>> {
    if cli.glyco && cli.mbr && !cli.lfq {
        return Err(
            "--glyco --mbr needs --lfq: match-between-runs transfers the identified intact \
             glycopeptides into the other runs"
                .into(),
        );
    }
    let quant: Option<QuantCollector> =
        QuantSettings::from_cli(cli, high_res, tmt_level)?.map(|mut settings| {
            settings.fragment_tol = Some(scorer.feature_match_tolerance());
            QuantCollector::new(settings)
        });
    if let Some(q) = &quant {
        let s = &q.settings;
        if is_mgf && (s.needs_ms1_series() || s.needs_ms3()) {
            return Err(
                "--lfq, --quant-dir and --tmt-level 3 need MS1/MS3 scans; MGF input carries MS2 only".into(),
            );
        }
        if is_d && s.needs_ms3() {
            return Err(
                "--tmt-level 3 reads reporter ions from MS3 scans, and a Bruker .d has none \
                 (a TDF records MS1 and PASEF MS2 frames only); quantify the reporters from \
                 the MS2 (--tmt-level 2)"
                    .into(),
            );
        }
        if cli.chimeric && s.needs_ms3() {
            return Err(
                "--tmt-level 3 cannot be combined with --chimeric: the cascade reads MS2 only"
                    .into(),
            );
        }
        if cli.glyco {
            let mut unsupported: Vec<&str> = Vec::new();
            if s.plex.is_some() {
                unsupported.push("--tmt");
            }
            if cli.quant_dir.is_some() {
                unsupported.push("--quant-dir");
            }
            if !unsupported.is_empty() {
                return Err(format!(
                    "--glyco does not support: {} (label-free quantification of intact \
                     glycopeptides is `--glyco --lfq`; the isobaric and match-between-runs \
                     stages are not implemented for glyco mode)",
                    unsupported.join(", ")
                )
                .into());
            }
        }
        if !cli.glyco && !cli.rescore && (s.plex.is_some() || s.lfq) {
            log::warn!(
                "--tmt/--lfq without --rescore: every rank-1 target \
                 PSM is quantified; filter the tables by the q-values you compute downstream, \
                 or add --rescore so only PSMs at q<=--psm-fdr are quantified."
            );
        }
        if !cli.glyco {
            log::info!(target: "quant", "{}", s.describe());
        }
    }
    if cli.dia && cli.quant_dir.is_some() && spectrum_paths.len() > 1 {
        return Err(
            "--dia --quant-dir takes one input per process (each run's pass-2 table and PIN go \
             to the quant directory); run andes once per run, then `andes quant finalize`"
                .into(),
        );
    }
    Ok(quant)
}
