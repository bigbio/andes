pub(crate) mod matching;
pub(crate) mod quantities;

use std::path::{Path, PathBuf};

use crate::cli::Cli;
use model::Spectrum;
use search::{PreparedSearch, SearchIndex, SearchParams};

const DEFAULT_GLYCO_TOL_PPM: f64 = 20.0;
const LOW_RES_GLYCO_TOL_PPM: f64 = 500.0;

fn glyco_fragment_tol_ppm(configured: f64, explicit: bool, acquisition_high_res: bool) -> f64 {
    if explicit || acquisition_high_res {
        configured
    } else {
        LOW_RES_GLYCO_TOL_PPM
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_glyco(
    cli: &Cli,
    spectra: &[Spectrum],
    prepared: &PreparedSearch,
    params: &SearchParams,
    idx: &SearchIndex,
    output_pin_path: &Path,
    spectrum_paths: &[PathBuf],
    acquisition_high_res: bool,
    mono: Option<&[Option<search::precursor::monoisotopic::MonoCorrection>]>,
    isotope_error_override: Option<(std::ops::RangeInclusive<i8>, Vec<bool>)>,
    quant: Option<&crate::quant::run::QuantCollector>,
    t_total: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let t_glyco = std::time::Instant::now();
    scoring_crate::scoring::init_cz_settings(scoring_crate::scoring::CzSettings {
        zmax_override: cli.glyco_cz_max_charge,
    });
    glyco::backbone::init_y_max_charge(cli.glyco_y_max_charge);

    let (glycan_list, source_desc) = match (
        cli.glyco_glycan_gdb.as_ref(),
        cli.glyco_glycan_wurcs.as_ref(),
        cli.glyco_species.as_ref(),
    ) {
        (Some(path), _, _) => {
            let content = std::fs::read_to_string(path)
                .map_err(|e| format!("--glyco-glycan-gdb {}: {e}", path.display()))?;
            let list = glyco::glycan::db::load_glycan_gdb(&content)?;
            (list, format!("{} (file)", path.display()))
        }
        (None, Some(path), _) => {
            let content = std::fs::read_to_string(path)
                .map_err(|e| format!("--glyco-glycan-wurcs {}: {e}", path.display()))?;
            let list = glyco::glycan::db::load_glycan_wurcs_csv(&content)?;
            (list, format!("{} (WURCS)", path.display()))
        }
        (None, None, Some(species)) => {
            let list = glyco::glycan::db::load_species_glycan_db(species.key())?;
            (list, format!("--glyco-species {}", species.key()))
        }
        (None, None, None) => {
            return Err(
                "--glyco requires --glyco-glycan-gdb, --glyco-glycan-wurcs or --glyco-species: \
                 the built-in composition enumerator was removed; supply a `.gdb` glycan \
                 database, a WURCS 2.0 CSV, or a bundled species"
                    .into(),
            );
        }
    };
    if glycan_list.is_empty() {
        return Err(format!(
            "glycan database {source_desc} contains no glycans: `--glyco` has no search space"
        )
        .into());
    }
    log::info!(
        target: "glyco",
        "glycan list: {} compositions loaded from {} (structure preserved; NeuGc content used as-is)",
        glycan_list.len(),
        source_desc
    );
    let glyco_tol_ppm = glyco_fragment_tol_ppm(
        cli.glyco_tol_ppm,
        crate::arg_present("--glyco-tol-ppm") || cli.glyco_tol_ppm != DEFAULT_GLYCO_TOL_PPM,
        acquisition_high_res,
    );
    if glyco_tol_ppm != cli.glyco_tol_ppm {
        log::info!(
            target: "glyco",
            "fragment tolerance: {glyco_tol_ppm} ppm for low-resolution MS2 \
             (--glyco-tol-ppm sets it)"
        );
    }
    if glyco_tol_ppm < DEFAULT_GLYCO_TOL_PPM {
        log::warn!(
            target: "glyco",
            "--glyco-tol-ppm {glyco_tol_ppm} is tighter than the 20 ppm the \
             glyco defaults were validated at; oxonium, core-Y and c/z matching may \
             under-fire"
        );
    }
    let spectra_for_glyco: &[_] = if cli.max_spectra > 0 {
        &spectra[..spectra.len().min(cli.max_spectra)]
    } else {
        spectra
    };
    let retrieval_ppm: Option<f64> = cli
        .glyco_retrieval_tol_ppm
        .or_else(|| acquisition_high_res.then_some(glyco_tol_ppm));
    match retrieval_ppm {
        Some(ppm) => log::info!(
            target: "glyco",
            "retrieval window: {ppm} ppm ({})",
            if cli.glyco_retrieval_tol_ppm.is_some() {
                "--glyco-retrieval-tol-ppm"
            } else {
                "high-resolution acquisition"
            }
        ),
        None => log::info!(
            target: "glyco",
            "retrieval window: model fragment tolerance in Da (low-resolution acquisition)"
        ),
    }
    let glyco_cfg = search::glycopeptides::GlycoConfig {
        core: cli.glyco_core.into(),
        min_core_y: cli.glyco_min_core_y,
        max_gen_peaks: cli.glyco_max_peaks,
        scan_filter_path: cli.glyco_scans.clone(),
        isotope_error_override_mask: isotope_error_override
            .as_ref()
            .map(|(_, mask)| mask.clone()),
        isotope_error_override: isotope_error_override.map(|(window, _)| window),
        retrieval_tol_ppm: retrieval_ppm,
        full_glycan_db: !cli.glyco_peptide_first,
        debug: cli.debug_glyco,
        hcd_pair: cli.glyco_hcd_pair && spectrum_paths.len() == 1,
    };
    let t_table = std::time::Instant::now();
    let table =
        search::glycopeptides::GlycoCandidates::from_index(prepared, cli.glyco_core.into());
    log::info!(
        target: "glyco",
        "candidate table: {} {}-bearing candidates from the out-of-core index \
         ({} base peptides) [{:.2}s]",
        table.candidates().len(),
        match cli.glyco_core {
            crate::cli::GlycoCoreFlag::N => "N-X-S/T",
            crate::cli::GlycoCoreFlag::O => "S/T",
        },
        prepared.mmap_index.len(),
        t_table.elapsed().as_secs_f64()
    );
    let pass1 = search::glycopeptides::glyco_search_run(
        spectra_for_glyco,
        prepared,
        &table,
        &glycan_list,
        glyco_tol_ppm,
        glyco_cfg,
    );
    let total_pass1_rows: usize = pass1.iter().map(|r| r.hits.len()).sum();
    log::info!(
        target: "glyco",
        "scored {} spectra → {} glyco-PSM rows [{:.2}s]",
        pass1.len(),
        total_pass1_rows,
        t_glyco.elapsed().as_secs_f64()
    );

    let glyco_pin_path = {
        let mut s = if output_pin_path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("pin"))
        {
            output_pin_path.with_extension("").into_os_string()
        } else {
            output_pin_path.to_path_buf().into_os_string()
        };
        s.push(".glyco.pin");
        std::path::PathBuf::from(s)
    };
    log::info!(target: "glyco", "PIN will be written to: {}", glyco_pin_path.display());
    let mut glyco_results = pass1;

    let total_glyco_rows: usize = glyco_results.iter().map(|r| r.hits.len()).sum();
    output::populate_glyco_rt_features(
        spectra,
        &mut glyco_results,
        table.candidates(),
        &glycan_list,
    );

    output::write_glyco_pin(
        &glyco_pin_path,
        spectra,
        &glyco_results,
        table.candidates(),
        params,
        idx,
        cli.debug_glyco,
        mono,
        cli.glyco_core.into(),
    )?;
    log::info!(
        target: "glyco",
        "wrote glyco PIN {} ({} PSM rows); total {:.2} s",
        glyco_pin_path.display(),
        total_glyco_rows,
        t_total.elapsed().as_secs_f64()
    );

    if let Some(collector) = quant {
        let runs = collector.run_ms1();
        if runs.is_empty() {
            return Err(
                "--glyco --lfq: no run's MS1 scans were captured, so no intact glycopeptide \
                 can be quantified. Label-free quantification needs input that carries MS1 \
                 scans (mzML, Thermo .raw or Bruker .d; not MGF)"
                    .into(),
            );
        }
        crate::glyco::quantities::run(
            &glyco_pin_path,
            spectra,
            &glyco_results,
            table.candidates(),
            idx,
            &runs,
            cli,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod glyco_tol_tests {
    use super::*;

    #[test]
    fn low_res_acquisition_widens_the_default_glyco_tolerance_only() {
        assert_eq!(glyco_fragment_tol_ppm(20.0, false, true), 20.0);
        assert_eq!(
            glyco_fragment_tol_ppm(20.0, false, false),
            LOW_RES_GLYCO_TOL_PPM
        );
        assert_eq!(glyco_fragment_tol_ppm(20.0, true, false), 20.0);
        assert_eq!(glyco_fragment_tol_ppm(30.0, true, true), 30.0);
    }
}
