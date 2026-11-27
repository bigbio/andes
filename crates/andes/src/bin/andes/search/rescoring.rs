//! Rescoring a search's PIN with Percolator: the native engine or an external binary, the DIA
//! RT refit and precursor-level q, the DIA precursor PIN, and the per-row scores for site
//! posteriors.

use std::path::{Path, PathBuf};

use model::Spectrum;

use crate::cli::Cli;

/// What rescoring reads from the search.
#[derive(Clone, Copy)]
pub(super) struct RescoreInputs<'a> {
    pub(super) cli: &'a Cli,
    pub(super) output_pin_path: &'a PathBuf,
    pub(super) spectrum_path: &'a PathBuf,
    pub(super) scratch_dir: &'a PathBuf,
    pub(super) use_native: bool,
    pub(super) use_percolator: bool,
    /// The search no longer needs the spectra's peaks after precursor-level q (DIA fact rows).
    pub(super) free_peaks: bool,
    pub(super) has_site_claims: bool,
}

/// What rescoring yields: the target PSMs by id (none without `--rescore`), the decoy ids and PSMs,
/// the per-row scores, the DIA precursor ids and the rows the precursor PIN replaced, and where the
/// MBR scores come from.
pub(super) struct Rescoring {
    pub(super) rescore_map: Option<std::collections::HashMap<String, output::PercolatorPsm>>,
    pub(super) mbr_source: String,
    pub(super) rescore_decoys: Option<std::collections::HashSet<String>>,
    pub(super) rescore_decoy_psms: std::collections::HashMap<String, output::PercolatorPsm>,
    pub(super) row_scores: Option<output::RowScores>,
    pub(super) dia_replaced: std::collections::HashSet<String>,
    pub(super) dia_precursor_ids: std::collections::HashSet<String>,
    pub(super) dia_precursor_charges: std::collections::HashMap<String, u8>,
}

/// The DIA precursor learned score trained beside the steps before rescoring, and the pass-2
/// table it was trained on.
pub(super) type EarlyLearned =
    std::thread::JoinHandle<(output::dia::pass2::Pass2Table, crate::dia::score::Learned)>;

/// The Percolator arguments of the run: `--percolator-args`, on DIA with the default train FDR.
fn user_args(cli: &Cli) -> Vec<String> {
    let mut args: Vec<String> = cli
        .percolator_args
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();
    if cli.dia {
        crate::dia::default_train_fdr(&mut args);
    }
    args
}

/// Starts the DIA precursor learned score on its own thread as soon as pass 2 has its table,
/// when rescoring will use it (`--rescore`, DIA, no `--quant-dir`): the table moves into the
/// thread and comes back with the score when [`rescore_pin`] joins it. The score depends only
/// on the table and the seed, so it is the one [`rescore_pin`] would train itself.
pub(super) fn start_learned(
    cli: &Cli,
    dia_pass2: &mut Option<output::dia::pass2::Pass2Table>,
) -> Option<EarlyLearned> {
    if !(cli.rescore && cli.dia && cli.quant_dir.is_none()) {
        return None;
    }
    let t = dia_pass2.take()?;
    let seed = crate::dia::score::percolator_seed(&user_args(cli));
    Some(std::thread::spawn(move || {
        let learned = crate::dia::score::learn(&t, &[], seed);
        (t, learned)
    }))
}

/// Rescores the search's PIN when `--rescore` is set; takes `dia_pass2` when it is used.
/// `early` is the learned score [`start_learned`] started, joined after pass 1.
pub(super) fn rescore_pin(
    inputs: &RescoreInputs<'_>,
    spectra: &mut [Spectrum],
    dia_pass2: &mut Option<output::dia::pass2::Pass2Table>,
    early: Option<EarlyLearned>,
) -> Result<Rescoring, Box<dyn std::error::Error>> {
    let RescoreInputs {
        cli,
        output_pin_path,
        spectrum_path,
        scratch_dir,
        use_native,
        use_percolator,
        free_peaks,
        has_site_claims,
    } = *inputs;
    let mut mbr_source = String::new();
    let mut rescore_decoys: Option<std::collections::HashSet<String>> = None;
    let mut rescore_decoy_psms: std::collections::HashMap<String, output::PercolatorPsm> =
        std::collections::HashMap::new();
    let mut row_scores: Option<output::RowScores> = None;
    let mut dia_replaced: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut dia_precursor_ids: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut dia_precursor_charges: std::collections::HashMap<String, u8> =
        std::collections::HashMap::new();
    let rescore_map = if use_percolator {
        let backend = if use_native {
            crate::cli::ensure_scratch_dir(scratch_dir)?;
            output::PercolatorBackend::Native {
                memory: crate::memory::rescore_memory(cli.rescore_memory_gb, scratch_dir.clone()),
            }
        } else {
            output::resolve_backend(
                cli.percolator_bin.as_deref(),
                cli.percolator_docker,
                &cli.percolator_image,
            )
            .map_err(|e| format!("{e}"))?
        };
        log::info!(target: "rescore", "Percolator backend = {}", backend.describe());
        let user_args = user_args(cli);
        let extra = output::percolator_args_with_threads(&user_args, cli.threads);
        let t_resc = std::time::Instant::now();
        if cli.dia && free_peaks {
            for s in spectra.iter_mut() {
                s.peaks = Vec::new();
            }
            crate::memory::release_free_memory();
        }
        let seed = crate::dia::score::percolator_seed(&user_args);
        let pass1 = rescore_pass1(cli, output_pin_path, spectra, &backend, &extra);
        let learned = match early {
            Some(h) => {
                let (t, learned) = h
                    .join()
                    .map_err(|_| "the precursor learned score panicked".to_string())?;
                *dia_pass2 = Some(t);
                Some(learned)
            }
            None => None,
        };
        let pass1 = pass1?;
        let Pass1 {
            targets: mut map,
            decoy_ids: mut decoys,
            decoys: pass1_decoys,
            q_count: pass1_q,
            best: pass1_best,
            removed: pass1_removed,
        } = pass1;
        rescore_decoy_psms = pass1_decoys;
        let mut results_pin = output_pin_path.clone();
        if let (Some(t), Some(dir)) = (dia_pass2.as_ref(), cli.quant_dir.as_deref()) {
            crate::dia::study::store_run(
                dir,
                spectrum_path,
                output_pin_path,
                t,
                Some((&backend, &user_args)),
                &pass1_removed,
            )?;
            crate::memory::mark_stage("study tables stored");
            *dia_pass2 = None;
        }
        if let Some(t) = dia_pass2.take() {
            let rescore = |pin: &Path| -> Result<crate::dia::score::Rescored, String> {
                let (targets, decoy_ids) =
                    output::run_percolator_with_decoys(&backend, pin, &extra)
                        .map_err(|e| format!("{e}"))?;
                let path = output::decoy_results_path(pin);
                let text = std::fs::read_to_string(&path)
                    .map_err(|e| format!("{}: {e}", path.display()))?;
                let decoys = output::parse_psm_results(&text)
                    .map_err(|e| format!("{}: {e}", path.display()))?;
                Ok(crate::dia::score::Rescored {
                    targets,
                    decoy_ids,
                    decoys,
                })
            };
            match crate::dia::score::rescore_precursors(
                output_pin_path,
                &t,
                seed,
                &pass1_best,
                &pass1_removed,
                rescore,
                learned,
            ) {
                Ok(r) => {
                    dia_precursor_charges = t
                        .spec_id
                        .iter()
                        .cloned()
                        .zip(t.charge.iter().copied())
                        .collect();
                    log::info!(
                        target: "dia",
                        "precursor PIN used: {} target precursors at q <= 0.01 (pass 1: \
                         {pass1_q}), {} degenerate learned-score folds",
                        r.n_q,
                        r.ids.degenerate_folds
                    );
                    map = r.rescored.targets;
                    decoys = r.rescored.decoy_ids;
                    rescore_decoy_psms = r.rescored.decoys;
                    dia_replaced = r.replaced;
                    dia_precursor_ids = r.ids.targets.into_iter().chain(r.ids.decoys).collect();
                    results_pin = r.pin;
                }
                Err(why) => log::info!(
                    target: "dia",
                    "pass-1 results kept, the precursor PIN is not used: {why}"
                ),
            }
            crate::memory::mark_stage("precursor rescoring");
        }
        rescore_decoys = Some(decoys);
        log::info!(
            target: "rescore",
            "Percolator produced {} target PSMs [{:.2}s]",
            map.len(),
            t_resc.elapsed().as_secs_f64()
        );
        if has_site_claims {
            row_scores = Some(
                output::score_all_rows(&backend, &results_pin, &extra)
                    .map_err(|e| format!("scoring every PIN row for the site posteriors: {e}"))?,
            );
        }
        if cli.mbr {
            mbr_source = format!("the Percolator results of {}", results_pin.display());
        }
        Some(map)
    } else {
        None
    };
    Ok(Rescoring {
        rescore_map,
        mbr_source,
        rescore_decoys,
        rescore_decoy_psms,
        row_scores,
        dia_replaced,
        dia_precursor_ids,
        dia_precursor_charges,
    })
}

/// Pass 1's rescored PSMs: target rows, decoy ids and decoy rows by SpecId, and on DIA the
/// precursor-level target count at q <= 0.01 and the best `(score, is decoy)` of every precursor.
struct Pass1 {
    targets: std::collections::HashMap<String, output::PercolatorPsm>,
    decoy_ids: std::collections::HashSet<String>,
    decoys: std::collections::HashMap<String, output::PercolatorPsm>,
    q_count: usize,
    best: Vec<(f64, bool)>,
    removed: std::collections::HashSet<(String, String, u8)>,
}

/// Rescores the pass-1 PIN: Percolator, on DIA the RT refit and a second Percolator run when
/// the refit changed the PIN, then on DIA the precursor-level q-values.
fn rescore_pass1(
    cli: &Cli,
    pin: &Path,
    spectra: &[Spectrum],
    backend: &output::PercolatorBackend,
    extra: &[String],
) -> Result<Pass1, String> {
    let (mut map, mut decoys) =
        output::run_percolator_with_decoys(backend, pin, extra).map_err(|e| format!("{e}"))?;
    crate::memory::mark_stage("rescoring");
    if cli.dia
        && crate::dia::refine_rt(pin, spectra, &|id: &str| {
            map.get(id).is_some_and(|p| p.q_value <= 0.01)
        })?
    {
        (map, decoys) =
            output::run_percolator_with_decoys(backend, pin, extra).map_err(|e| format!("{e}"))?;
        crate::memory::mark_stage("rescoring after the RT refit");
    }
    let path = output::decoy_results_path(pin);
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let decoy_psms =
        output::parse_psm_results(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut q_count = 0;
    let mut best = Vec::new();
    let mut removed = std::collections::HashSet::new();
    if cli.dia {
        let ((t, d, n), b, r) =
            crate::dia::precursor_q(pin, spectra, &mut map, &decoy_psms, cli.dia_interference)?;
        q_count = n;
        best = b;
        removed = r;
        crate::memory::mark_stage("precursor q");
        log::info!(
            target: "rescore",
            "Precursor-level q (DIA): {t} target and {d} decoy precursors, {n} at q <= 0.01"
        );
    }
    Ok(Pass1 {
        targets: map,
        decoy_ids: decoys,
        decoys: decoy_psms,
        q_count,
        best,
        removed,
    })
}
