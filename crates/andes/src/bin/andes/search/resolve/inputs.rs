//! Run inputs and paths: scratch directories, input detection, `.d` archives and the output
//! and rescoring destinations.

use std::path::{Path, PathBuf};

use super::super::scratch_tempdir;
use super::checks::{check_mbr_args, deuterated_tag_warning, warn_search};
use super::{isobaric_label_of_plex, RunInputs, RunPaths, TempDirs};
use crate::acquisition::Acquisition;
use crate::cli::{Cli, EthcdActivationFlag, Protocol};
use crate::models::{label_from_mods, routing_protocol};
use crate::quant::run::QuantSettings;
use crate::spectra::{
    compare_with_first_input, detect_acquisition_metadata, detect_isobaric_sampled,
    input_format_flags, sample_ms2_meta, sample_spectra,
};

/// Creates the scratch directory, unpacks archived inputs, validates the quantification
/// arguments and decides where the PIN goes.
pub(crate) fn prepare_run(cli: &mut Cli) -> Result<RunPaths, Box<dyn std::error::Error>> {
    let scratch_dir = crate::cli::dir_or_output_dir(
        cli.temp_dir.as_deref(),
        &output_base(
            cli.output_pin.as_deref(),
            cli.output_tsv.as_deref(),
            cli.output_parquet.as_deref(),
        ),
    );
    cli.temp_dir = Some(scratch_dir.clone());
    let original_spectrum_paths = cli.spectrum.clone();
    let _d_archives = unpack_d_archives(&mut cli.spectrum, &scratch_dir)?;
    let mut _mbr_tmp: Option<tempfile::TempDir> = None;
    if cli.mbr && !cli.glyco {
        check_mbr_args(cli)?;
        crate::quant::run::check_unique_stems(&cli.spectrum)?;
        if cli.quant_dir.is_none() && cli.label.is_none() {
            let dir = scratch_tempdir("andes-quant-", &scratch_dir)?;
            cli.quant_dir = Some(dir.path().to_path_buf());
            _mbr_tmp = Some(dir);
        }
    }
    if cli.quant_dir.is_some() {
        crate::quant::run::check_unique_stems(&cli.spectrum)?;
    }
    QuantSettings::check_args(cli)?;
    if let Some(warning) = deuterated_tag_warning(cli.tmt, cli.precursor_tol) {
        warn_search(&warning);
    }
    let mut _rescore_tmp: Option<tempfile::TempDir> = None;
    let (output_pin_path, report_base) = rescore_pin_destination(
        cli.output_pin.as_deref(),
        cli.output_tsv.as_deref(),
        cli.output_parquet.as_deref(),
        cli.keep_pin,
        &scratch_dir,
        &mut _rescore_tmp,
    )?;
    if cli.spectrum.len() > 1 {
        log::info!(target: "search",
            "Multi-spectrum search: {} inputs → one PIN",
            cli.spectrum.len()
        );
    }
    Ok(RunPaths {
        scratch_dir,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        temp_dirs: TempDirs {
            _rescore: _rescore_tmp,
            _mbr: _mbr_tmp,
            _d_archives,
        },
    })
}

/// Detects the input format, the acquisition metadata and any isobaric label.
pub(crate) fn detect_inputs(
    cli: &Cli,
    aa: &model::AminoAcidSet,
) -> Result<RunInputs, Box<dyn std::error::Error>> {
    let spectrum_paths = &cli.spectrum;
    let spectrum_path: PathBuf = spectrum_paths[0].clone();
    let (is_mzml, is_raw, is_d, is_mgf) = input_format_flags(&spectrum_path);

    input::mzml::init_ethcd_as_etd(cli.ethcd_activation == EthcdActivationFlag::Etd);

    if is_raw {
        #[cfg(feature = "thermo")]
        input::ThermoRawReader::open(&spectrum_path).map_err(|e| e.to_string())?;
        #[cfg(not(feature = "thermo"))]
        return Err("this andes build has no Thermo .raw support; rebuild with \
                    `--features thermo` (and run with the .NET 8 runtime installed). \
                    mzML/MGF inputs work without it."
            .into());
    }

    let (detected_activation, detected_instrument) = detect_acquisition_metadata(&spectrum_path);
    let instrument_was_detected = detected_instrument.is_some();
    for other_path in &spectrum_paths[1..] {
        let other = detect_acquisition_metadata(other_path);
        if let Some(warning) = compare_with_first_input(
            &spectrum_path,
            (detected_activation, detected_instrument),
            other_path,
            other,
        )? {
            warn_search(&warning);
        }
    }

    let peak_evidence = if !instrument_was_detected
        && cli.fragment_tol_ppm.is_none()
        && cli.fragment_tol_da.is_none()
    {
        let ev = sample_spectra(&spectrum_path, is_mzml, is_mgf, is_raw, cli.ms_level as u32)
            .map(|sample| input::infer_fragment_resolution(&sample));
        if let Some(e) = &ev {
            log::info!(
                target: "acquisition",
                "no analyzer metadata; fragment resolution from the peaks: {} -> {}",
                e.describe(),
                match e.high_res {
                    Some(true) => "high resolution (HCD high-res model)",
                    Some(false) => "low resolution",
                    None => "undecided",
                }
            );
        }
        ev
    } else {
        None
    };
    let acquisition = Acquisition::resolve(
        detected_activation,
        detected_instrument,
        cli.fragment_tol_ppm,
        cli.fragment_tol_da,
        peak_evidence,
        sample_ms2_meta(&spectrum_path, is_mzml, is_raw, is_d),
    );
    if let (Some(said), Some(read)) = (cli.sdrf_dissociation, acquisition.activation) {
        if said != read {
            log::warn!(target: "search",
                "the SDRF dissociation method is {} but {} records {}; andes follows the file",
                said.name(),
                spectrum_path.display(),
                read.name()
            );
        }
    }

    let detect = || {
        detect_isobaric_sampled(
            &spectrum_path,
            is_mzml,
            is_mgf,
            is_raw,
            cli.ms_level as u32,
            acquisition.high_res,
        )
    };
    let label: Option<input::IsobaricLabel> = match (cli.tmt, cli.protocol) {
        (Some(plex), _) => Some(isobaric_label_of_plex(plex)),
        (None, Protocol::Auto | Protocol::Tmt | Protocol::Itraq | Protocol::ItraqPhospho) => {
            detect()
        }
        (None, _) => None,
    };
    let label = label.or_else(|| {
        if cli.protocol == Protocol::Auto {
            label_from_mods(aa)
        } else {
            None
        }
    });
    let route_protocol = routing_protocol(cli.protocol, label);
    Ok(RunInputs {
        is_mzml,
        is_raw,
        is_d,
        is_mgf,
        detected_activation,
        detected_instrument,
        instrument_was_detected,
        acquisition,
        label,
        route_protocol,
    })
}

fn output_base(
    output_pin: Option<&Path>,
    output_tsv: Option<&Path>,
    output_parquet: Option<&Path>,
) -> PathBuf {
    output_pin
        .or(output_tsv)
        .or(output_parquet)
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("andes.pin"))
}

pub(super) fn rescore_pin_destination(
    output_pin: Option<&Path>,
    output_tsv: Option<&Path>,
    output_parquet: Option<&Path>,
    keep_pin: bool,
    scratch_dir: &Path,
    tmp: &mut Option<tempfile::TempDir>,
) -> Result<(PathBuf, PathBuf), Box<dyn std::error::Error>> {
    let base = output_base(output_pin, output_tsv, output_parquet);
    if output_pin.is_some() {
        return Ok((base.clone(), base));
    }
    if keep_pin {
        return Ok((base.with_extension("pin"), base));
    }
    let dir = scratch_tempdir("andes-pin-", scratch_dir)?;
    let path = dir.path().join("andes.pin");
    *tmp = Some(dir);
    Ok((path, base))
}

#[cfg(feature = "timstof")]
pub(super) fn check_distinct_d_archive_runs(inputs: &[PathBuf]) -> Result<(), String> {
    let mut seen: std::collections::HashMap<String, &PathBuf> = std::collections::HashMap::new();
    for p in inputs {
        let Some((_, run)) = input::d_archive_kind(p) else {
            continue;
        };
        if let Some(first) = seen.insert(run.clone(), p) {
            return Err(format!(
                "the archives {} and {} both hold run '{run}'; give each run one archive",
                first.display(),
                p.display()
            ));
        }
    }
    Ok(())
}

fn unpack_d_archives(
    inputs: &mut [PathBuf],
    scratch_dir: &Path,
) -> Result<Option<tempfile::TempDir>, Box<dyn std::error::Error>> {
    if !inputs.iter().any(|p| crate::spectra::is_d_archive(p)) {
        return Ok(None);
    }
    #[cfg(feature = "timstof")]
    {
        check_distinct_d_archive_runs(inputs)?;
        let dir = scratch_tempdir("andes-d-", scratch_dir)?;
        for p in inputs
            .iter_mut()
            .filter(|p| crate::spectra::is_d_archive(p))
        {
            let t = std::time::Instant::now();
            let run = input::unpack_d_archive(p, dir.path())?;
            log::info!(target: "search",
                "Bruker .d archive {}: unpacked to {} ({:.1}s)",
                p.display(),
                run.display(),
                t.elapsed().as_secs_f64()
            );
            *p = run;
        }
        Ok(Some(dir))
    }
    #[cfg(not(feature = "timstof"))]
    {
        let _ = scratch_dir;
        Err(
            "this andes build has no Bruker .d (timsTOF) support; rebuild with \
             `--features timstof` to read archived .d inputs"
                .into(),
        )
    }
}
