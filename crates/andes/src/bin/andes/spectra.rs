use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::mpsc::SyncSender;

use input::{detect_instrument_type, MgfReader, Ms1Link, MzMLReader};
use model::scan::Ms1Scan;
use model::{activation::ActivationMethod, InstrumentType, Spectrum, Tolerance};
use search::precursor_cal::{constants as cal_constants, sample_every_nth};
use search::{
    learn_calibration_stats, CalibrationStats, PrecursorCalMode, PreparedSearch, SearchParams,
    SpecKey,
};

#[derive(Debug, Default)]
pub(crate) struct ParseStats {
    pub(crate) error_count: usize,
    pub(crate) first_errors: Vec<String>,
}

pub(crate) fn spectrum_ext_lower(path: &std::path::Path) -> Option<String> {
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    let effective: std::path::PathBuf = if is_gz {
        path.file_stem()
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|| path.to_path_buf())
    } else {
        path.to_path_buf()
    };
    effective
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
}

pub(crate) fn is_d_archive(path: &Path) -> bool {
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    let lower = name.to_ascii_lowercase();
    [".d.zip", ".d.tar", ".d.tar.gz", ".d.tgz", ".d.gz"]
        .iter()
        .any(|s| lower.len() > s.len() && lower.ends_with(s))
}

pub(crate) fn input_format_flags(path: &Path) -> (bool, bool, bool, bool) {
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    let effective: std::path::PathBuf = if is_gz {
        path.file_stem()
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|| path.to_path_buf())
    } else {
        path.to_path_buf()
    };
    let ext = effective
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase());
    let is_mzml = matches!(ext.as_deref(), Some("mzml"));
    let is_raw = matches!(ext.as_deref(), Some("raw"));
    let is_d = matches!(ext.as_deref(), Some("d")) || is_d_archive(path);
    let is_mgf = !is_mzml && !is_raw && !is_d;
    (is_mzml, is_raw, is_d, is_mgf)
}

pub(crate) fn title_prefix_for(file_stem: &str) -> Option<String> {
    Some(format!("{file_stem}/"))
}

pub(crate) fn prefix_spectrum_titles(chunk: &mut [Spectrum], prefix: &str) {
    const SEPARATORS: [char; 3] = ['\t', '\n', '\r'];
    for spec in chunk.iter_mut() {
        let title = if spec.title.is_empty() {
            format!("{prefix}scan={}", spec.scan.unwrap_or(0))
        } else {
            format!("{prefix}{}", spec.title)
        };
        spec.title = if title.contains(SEPARATORS) {
            title.replace(SEPARATORS, " ")
        } else {
            title
        };
    }
}

pub(crate) fn add_mzml_skipped<R: std::io::BufRead>(
    stats: &mut ParseStats,
    reader: &MzMLReader<R>,
) {
    stats.error_count += reader.skipped_count();
    for e in reader.skipped_errors() {
        if stats.first_errors.len() < 3 {
            stats.first_errors.push(e.clone());
        }
    }
}

pub(crate) fn merge_parse_stats(acc: &mut ParseStats, part: ParseStats) {
    acc.error_count += part.error_count;
    for e in part.first_errors {
        if acc.first_errors.len() < 10 {
            acc.first_errors.push(e);
        }
    }
}

pub(crate) fn with_first_read_error(msg: String, first_errors: &[String]) -> String {
    match first_errors.first() {
        Some(e) => format!("{msg} First read error: {e}"),
        None => msg,
    }
}

pub(crate) fn require_spectra(path: &Path, ms_level: u8, found: usize) -> Result<(), String> {
    if found > 0 {
        return Ok(());
    }
    let (is_mzml, is_raw, is_d, _) = input_format_flags(path);
    let format = if is_mzml {
        "mzML"
    } else if is_raw {
        "Thermo .raw"
    } else if is_d {
        "Bruker .d"
    } else {
        "MGF"
    };
    Err(format!(
        "no MS{ms_level} spectra could be read from {} ({format}): the file is empty or \
         corrupt, or has no MS{ms_level} scans with a precursor m/z",
        path.display()
    ))
}

pub(crate) fn send_chunks<R, E>(
    reader: R,
    chunk_size: usize,
    bench_cap: usize,
    tx: SyncSender<Vec<Spectrum>>,
) -> ParseStats
where
    R: Iterator<Item = Result<Spectrum, E>>,
    E: std::fmt::Display,
{
    let mut stats = ParseStats::default();
    let mut chunk: Vec<Spectrum> = Vec::with_capacity(chunk_size);
    let mut total = 0usize;
    for result in reader {
        if total >= bench_cap {
            break;
        }
        match result {
            Ok(s) => {
                chunk.push(s);
                total += 1;
                if chunk.len() >= chunk_size {
                    let payload = std::mem::replace(&mut chunk, Vec::with_capacity(chunk_size));
                    if tx.send(payload).is_err() {
                        return stats;
                    }
                }
            }
            Err(e) => {
                stats.error_count += 1;
                if stats.first_errors.len() < 3 {
                    stats.first_errors.push(format!("{e}"));
                }
            }
        }
    }
    if !chunk.is_empty() {
        let _ = tx.send(chunk);
    }
    stats
}

pub(crate) fn send_linked_chunks(
    mut spectra: Vec<Spectrum>,
    ms2_to_ms1: &[Option<usize>],
    ms1: &[Ms1Scan],
    chunk_size: usize,
    tx: &SyncSender<(Vec<Spectrum>, Ms1Link)>,
) {
    for links in ms2_to_ms1.chunks(chunk_size) {
        let chunk: Vec<Spectrum> = spectra.drain(..links.len()).collect();
        let mut link = Ms1Link::default();
        let mut last: Option<usize> = None;
        for &l in links {
            link.ms2_to_ms1.push(l.map(|m| {
                if last != Some(m) {
                    link.ms1_peaks.push(ms1[m].peaks.clone());
                    last = Some(m);
                }
                link.ms1_peaks.len() - 1
            }));
        }
        if tx.send((chunk, link)).is_err() {
            return;
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct SpectrumMeta {
    pub(crate) precursor_charge: Option<i32>,
    pub(crate) num_peaks: usize,
}

pub(crate) fn scan_spectrum_metadata(
    path: &Path,
    is_mzml: bool,
    ms_level: u32,
    bench_cap: usize,
) -> Result<Vec<SpectrumMeta>, Box<dyn std::error::Error>> {
    let mut out = Vec::new();
    if input_format_flags(path).2 {
        return Ok(d_metadata(open_d_ms2(path)?, bench_cap));
    } else if is_mzml {
        let reader = MzMLReader::new(input::open_buf_maybe_gz(path)?)
            .with_ms_level_range(ms_level, ms_level);
        for result in reader {
            if out.len() >= bench_cap {
                break;
            }
            let spec = result.map_err(|e| format!("mzML parse: {e}"))?;
            out.push(SpectrumMeta {
                precursor_charge: spec.precursor_charge,
                num_peaks: spec.peaks.len(),
            });
        }
    } else {
        let reader = MgfReader::new(input::open_buf_maybe_gz(path)?);
        for result in reader {
            if out.len() >= bench_cap {
                break;
            }
            let spec = result.map_err(|e| format!("MGF parse: {e}"))?;
            out.push(SpectrumMeta {
                precursor_charge: spec.precursor_charge,
                num_peaks: spec.peaks.len(),
            });
        }
    }
    Ok(out)
}

fn readable_spectra<E>(
    items: impl Iterator<Item = Result<Spectrum, E>>,
    bench_cap: usize,
) -> impl Iterator<Item = (usize, Spectrum)> {
    items.filter_map(Result::ok).take(bench_cap).enumerate()
}

fn d_metadata<E>(
    items: impl Iterator<Item = Result<Spectrum, E>>,
    bench_cap: usize,
) -> Vec<SpectrumMeta> {
    readable_spectra(items, bench_cap)
        .map(|(_, spec)| SpectrumMeta {
            precursor_charge: spec.precursor_charge,
            num_peaks: spec.peaks.len(),
        })
        .collect()
}

fn d_load_by_index<E>(
    items: impl Iterator<Item = Result<Spectrum, E>>,
    indices: &HashSet<usize>,
    bench_cap: usize,
) -> HashMap<usize, Spectrum> {
    let mut loaded = HashMap::new();
    for (idx, spec) in readable_spectra(items, bench_cap) {
        if indices.contains(&idx) {
            loaded.insert(idx, spec);
            if loaded.len() == indices.len() {
                break;
            }
        }
    }
    loaded
}

type SpectrumIter = Box<dyn Iterator<Item = Result<Spectrum, Box<dyn std::error::Error>>>>;

fn open_d_ms2(path: &Path) -> Result<SpectrumIter, Box<dyn std::error::Error>> {
    #[cfg(feature = "timstof")]
    {
        let reader = input::TimsTofReader::open(path)?;
        Ok(Box::new(reader.map(|r| r.map_err(|e| e.into()))))
    }
    #[cfg(not(feature = "timstof"))]
    {
        Err(format!(
            "this andes build has no Bruker .d (timsTOF) support; rebuild with \
             `--features timstof`: {}",
            path.display()
        )
        .into())
    }
}

pub(crate) fn build_spec_keys_from_metadata(
    meta: &[SpectrumMeta],
    charge_range: std::ops::RangeInclusive<u8>,
    min_peaks: u32,
) -> Vec<SpecKey> {
    search::mass_calibrator::build_spec_keys_from_counts(
        meta.iter().map(|m| (m.precursor_charge, m.num_peaks)),
        &charge_range,
        min_peaks,
    )
}

pub(crate) fn load_spectra_by_index(
    path: &Path,
    is_mzml: bool,
    ms_level: u32,
    indices: &HashSet<usize>,
    bench_cap: usize,
) -> Result<HashMap<usize, Spectrum>, Box<dyn std::error::Error>> {
    let mut loaded = HashMap::new();
    if indices.is_empty() {
        return Ok(loaded);
    }
    if input_format_flags(path).2 {
        return Ok(d_load_by_index(open_d_ms2(path)?, indices, bench_cap));
    } else if is_mzml {
        let reader = MzMLReader::new(input::open_buf_maybe_gz(path)?)
            .with_ms_level_range(ms_level, ms_level);
        for (idx, result) in reader.enumerate() {
            if idx >= bench_cap {
                break;
            }
            if !indices.contains(&idx) {
                continue;
            }
            let spec = result.map_err(|e| format!("mzML parse: {e}"))?;
            loaded.insert(idx, spec);
            if loaded.len() == indices.len() {
                break;
            }
        }
    } else {
        let reader = MgfReader::new(input::open_buf_maybe_gz(path)?);
        for (idx, result) in reader.enumerate() {
            if idx >= bench_cap {
                break;
            }
            if !indices.contains(&idx) {
                continue;
            }
            let spec = result.map_err(|e| format!("MGF parse: {e}"))?;
            loaded.insert(idx, spec);
            if loaded.len() == indices.len() {
                break;
            }
        }
    }
    Ok(loaded)
}

pub(crate) fn detect_isobaric_sampled(
    path: &Path,
    is_mzml: bool,
    is_mgf: bool,
    is_raw: bool,
    ms_level: u32,
    high_res: bool,
) -> Option<input::IsobaricLabel> {
    let found = sample_spectra(path, is_mzml, is_mgf, is_raw, ms_level)
        .and_then(|sample| input::detect_isobaric(&sample, high_res));
    if found.is_some() || ms_level != 2 || !(is_mzml || is_raw) {
        return found;
    }
    let ms3 = sample_spectra(path, is_mzml, is_mgf, is_raw, 3)?;
    let label = input::detect_isobaric(&ms3, true);
    if label.is_some() {
        log::info!(target: "protocol", "isobaric reporters found in the MS3 scans (SPS-MS3)");
    }
    label
}

pub(crate) fn ms3_carries_reporters(path: &Path, is_mzml: bool, is_raw: bool) -> bool {
    if !(is_mzml || is_raw) {
        return false;
    }
    let Some(ms3) = sample_spectra(path, is_mzml, false, is_raw, 3) else {
        return false;
    };
    let found = input::detect_isobaric(&ms3, true).is_some();
    if !found {
        log::info!(
            target: "protocol",
            "no isobaric reporters in a sample of {} MS3 scans",
            ms3.len()
        );
    }
    found
}

pub(crate) fn sample_ms2_meta(
    path: &Path,
    is_mzml: bool,
    is_raw: bool,
    is_d: bool,
) -> input::Ms2MetaSample {
    const SAMPLE_N: usize = 1000;
    if is_mzml {
        return input::open_buf_maybe_gz(path)
            .map(|r| input::sample_mzml_ms2_meta(r, SAMPLE_N))
            .unwrap_or_default();
    }
    #[cfg(feature = "thermo")]
    if is_raw {
        return input::thermo::sample_ms2_meta(path, SAMPLE_N).unwrap_or_default();
    }
    #[cfg(feature = "timstof")]
    if is_d {
        return input::TimsTofReader::open(path)
            .map(|r| r.sample_ms2_meta(SAMPLE_N))
            .unwrap_or_default();
    }
    let _ = (is_raw, is_d);
    input::Ms2MetaSample::default()
}

pub(crate) fn sample_spectra(
    path: &Path,
    is_mzml: bool,
    is_mgf: bool,
    is_raw: bool,
    ms_level: u32,
) -> Option<Vec<Spectrum>> {
    const SAMPLE_N: usize = 1000;
    if is_raw {
        #[cfg(feature = "thermo")]
        {
            let level = u8::try_from(ms_level).ok()?;
            let mut reader = input::ThermoRawReader::open(path).ok()?;
            if level >= 3 {
                return Some(reader.sample_level(level, SAMPLE_N));
            }
            let reader = reader.with_ms_level(Some(level));
            return Some(reader.filter_map(Result::ok).take(SAMPLE_N).collect());
        }
        #[cfg(not(feature = "thermo"))]
        return None;
    }
    if input_format_flags(path).2 {
        #[cfg(feature = "timstof")]
        {
            if ms_level != 2 {
                return None;
            }
            let reader = input::TimsTofReader::open(path).ok()?;
            return Some(reader.filter_map(Result::ok).take(SAMPLE_N).collect());
        }
        #[cfg(not(feature = "timstof"))]
        return None;
    }
    if !(is_mzml || is_mgf) {
        return None;
    }
    let indices: HashSet<usize> = (0..SAMPLE_N).collect();
    let loaded = load_spectra_by_index(path, is_mzml, ms_level, &indices, usize::MAX).ok()?;
    Some(loaded.into_values().collect())
}

pub(crate) fn sample_spectra_strided(
    path: &Path,
    is_mzml: bool,
    ms_level: u32,
    want: usize,
) -> Option<Vec<Spectrum>> {
    if !is_mzml || want == 0 {
        return None;
    }
    let meta = scan_spectrum_metadata(path, true, ms_level, usize::MAX).ok()?;
    let n = meta.len();
    if n == 0 {
        return None;
    }
    let want = want.min(n);
    let stride = (n / want).max(1);
    let indices: HashSet<usize> = (0..n).step_by(stride).take(want).collect();
    let loaded = load_spectra_by_index(path, true, ms_level, &indices, usize::MAX).ok()?;
    let mut ordered: Vec<(usize, Spectrum)> = loaded.into_iter().collect();
    ordered.sort_by_key(|(i, _)| *i);
    Some(ordered.into_iter().map(|(_, s)| s).collect())
}

pub(crate) fn tolerance_ppm_display(t: Tolerance) -> Option<f64> {
    match t {
        Tolerance::Ppm(v) => Some(v),
        Tolerance::Da(_) => None,
    }
}

pub(crate) fn run_precursor_calibration(
    spectrum_path: &Path,
    is_mzml: bool,
    ms_level: u32,
    bench_cap: usize,
    params: &SearchParams,
    prepared: &mut PreparedSearch<'_>,
) -> Result<CalibrationStats, Box<dyn std::error::Error>> {
    if params.precursor_cal_mode == PrecursorCalMode::Off {
        return Ok(CalibrationStats::default());
    }

    let t_cal = std::time::Instant::now();
    let meta = scan_spectrum_metadata(spectrum_path, is_mzml, ms_level, bench_cap)?;
    let spec_keys =
        build_spec_keys_from_metadata(&meta, params.charge_range.clone(), params.min_peaks);

    if spec_keys.len() < params.cal_min_spec_keys {
        log::info!(
            target: "calibration",
            "calibration skipped ({} SpecKeys < {} threshold; elapsed: {:.2}s). \
             The sample is too small for a reliable calibration pre-pass.",
            spec_keys.len(),
            params.cal_min_spec_keys,
            t_cal.elapsed().as_secs_f64()
        );
        return Ok(CalibrationStats::default());
    }

    let sampled = sample_every_nth(
        &spec_keys,
        cal_constants::SAMPLING_STRIDE,
        cal_constants::MAX_SAMPLED,
    );
    let needed: HashSet<usize> = sampled.iter().map(|k| k.spectrum_idx).collect();
    let originals = load_spectra_by_index(spectrum_path, is_mzml, ms_level, &needed, bench_cap)?;

    let stats = learn_calibration_stats(&spec_keys, &originals, prepared, params);

    if stats.has_reliable_stats() {
        log::info!(
            target: "calibration",
            "shift learned: {:.3} ppm from {} confident PSMs (robust sigma {:.3} ppm; elapsed: {:.2}s)",
            stats.shift_ppm,
            stats.confident_psm_count,
            stats.robust_sigma_ppm,
            t_cal.elapsed().as_secs_f64()
        );
    } else {
        log::info!(
            target: "calibration",
            "calibration skipped (insufficient confident PSMs: {} with PSMs, {} below RawScore floor, {} failed |residual|>50ppm; elapsed: {:.2}s)",
            stats.queues_with_psm,
            stats.rejected_low_score,
            stats.rejected_residual,
            t_cal.elapsed().as_secs_f64()
        );
    }
    Ok(stats)
}

pub(crate) fn detect_dominant_activation(
    spectrum_path: &std::path::Path,
) -> Option<ActivationMethod> {
    let ext_lower = spectrum_ext_lower(spectrum_path);
    if ext_lower.as_deref() != Some("mzml") {
        return None;
    }

    const MAX_PEEK: usize = 64;

    let reader = MzMLReader::new(input::open_buf_maybe_gz(spectrum_path).ok()?);

    let mut counts: std::collections::HashMap<ActivationMethod, usize> =
        std::collections::HashMap::new();
    for (seen, item) in reader.enumerate() {
        if seen >= MAX_PEEK {
            break;
        }
        if let Ok(spec) = item {
            if let Some(m) = spec.activation_method {
                *counts.entry(m).or_insert(0) += 1;
            }
        }
    }

    if counts.is_empty() {
        return None;
    }

    let dominant = counts
        .iter()
        .max_by_key(|(&m, &n)| (n, std::cmp::Reverse(m as u8)))
        .map(|(&m, _)| m)?;

    if counts.len() > 1 {
        let mut other_pairs: Vec<(ActivationMethod, usize)> = counts
            .iter()
            .filter(|(&m, _)| m != dominant)
            .map(|(&m, &n)| (m, n))
            .collect();
        other_pairs.sort_by(|a, b| b.1.cmp(&a.1));
        let total: usize = counts.values().sum();
        let dominant_count = counts[&dominant];
        log::warn!(
            target: "model",
            "mixed activation methods in input ({} different methods \
             across {} peeked MS2 spectra). The run is routed by the dominant = {} ({}/{}); \
             each scan is scored with the model of its own activation (not with --model or \
             --fragmentation). Other methods present: {}",
            counts.len(),
            total,
            dominant.name(),
            dominant_count,
            total,
            other_pairs
                .iter()
                .map(|(m, n)| format!("{}={}", m.name(), n))
                .collect::<Vec<_>>()
                .join(", "),
        );
    }

    Some(dominant)
}

pub(crate) fn detect_instrument_type_for_path(
    spectrum_path: &std::path::Path,
) -> Option<InstrumentType> {
    let ext_lower = spectrum_ext_lower(spectrum_path);
    if ext_lower.as_deref() != Some("mzml") {
        return None;
    }

    detect_instrument_type(input::open_buf_maybe_gz(spectrum_path).ok()?)
}

pub(crate) fn detect_acquisition_metadata(
    spectrum_path: &Path,
) -> (Option<ActivationMethod>, Option<InstrumentType>) {
    let (is_mzml, is_raw, is_d, _is_mgf) = input_format_flags(spectrum_path);
    if is_mzml {
        (
            detect_dominant_activation(spectrum_path),
            detect_instrument_type_for_path(spectrum_path),
        )
    } else if is_raw {
        #[cfg(feature = "thermo")]
        {
            input::thermo::detect_activation_instrument(spectrum_path, 64)
        }
        #[cfg(not(feature = "thermo"))]
        {
            (None, None)
        }
    } else if is_d {
        (Some(ActivationMethod::CID), Some(InstrumentType::TimsTOF))
    } else {
        (None, None)
    }
}

pub(crate) fn compare_with_first_input(
    first_path: &Path,
    first: (Option<ActivationMethod>, Option<InstrumentType>),
    other_path: &Path,
    other: (Option<ActivationMethod>, Option<InstrumentType>),
) -> Result<Option<String>, String> {
    let class = |i: InstrumentType| {
        if i.is_high_resolution() {
            "high"
        } else {
            "low"
        }
    };
    let analyzer_name = |i: Option<InstrumentType>| i.map_or("not recorded", |i| i.name());
    let activation_name = |a: Option<ActivationMethod>| a.map_or("not recorded", |a| a.name());
    if let (Some(a), Some(b)) = (first.1, other.1) {
        if a.is_high_resolution() != b.is_high_resolution() {
            return Err(format!(
                "--spectrum inputs differ in fragment resolution: {} is {}-resolution MS2 ({}), \
                 {} is {}-resolution MS2 ({}). The model, the fragment tolerance and the score \
                 features are chosen once per search, from the first input; search these files \
                 separately.",
                first_path.display(),
                class(a),
                a.name(),
                other_path.display(),
                class(b),
                b.name(),
            ));
        }
    }
    if first == other {
        return Ok(None);
    }
    Ok(Some(format!(
        "WARN: {} differs from the first input {}: analyzer {} vs {}, dominant activation {} vs \
         {}. The run-level model and tolerances come from the first input.",
        other_path.display(),
        first_path.display(),
        analyzer_name(other.1),
        analyzer_name(first.1),
        activation_name(other.0),
        activation_name(first.0),
    )))
}

#[cfg(test)]
mod acquisition_metadata_tests {
    use super::{compare_with_first_input, detect_acquisition_metadata};
    use model::{activation::ActivationMethod, InstrumentType};
    use std::path::Path;

    const HCD: Option<ActivationMethod> = Some(ActivationMethod::HCD);
    const CID: Option<ActivationMethod> = Some(ActivationMethod::CID);

    fn compare(
        first: (Option<ActivationMethod>, Option<InstrumentType>),
        other: (Option<ActivationMethod>, Option<InstrumentType>),
    ) -> Result<Option<String>, String> {
        compare_with_first_input(Path::new("a.mzML"), first, Path::new("b.mzML"), other)
    }

    #[test]
    fn identical_inputs_are_silent() {
        let q = (HCD, Some(InstrumentType::QExactive));
        assert_eq!(compare(q, q), Ok(None));
        assert_eq!(compare((None, None), (None, None)), Ok(None));
    }

    #[test]
    fn different_resolution_classes_are_an_error() {
        let err = compare(
            (HCD, Some(InstrumentType::QExactive)),
            (CID, Some(InstrumentType::LowRes)),
        )
        .unwrap_err();
        assert!(err.contains("a.mzML is high-resolution MS2"), "{err}");
        assert!(err.contains("b.mzML is low-resolution MS2"), "{err}");
        assert!(compare(
            (None, Some(InstrumentType::LowRes)),
            (None, Some(InstrumentType::OrbitrapAstral)),
        )
        .is_err());
    }

    #[test]
    fn same_class_differences_are_a_warning() {
        let w = compare(
            (HCD, Some(InstrumentType::QExactive)),
            (HCD, Some(InstrumentType::OrbitrapAstral)),
        )
        .unwrap()
        .unwrap();
        assert!(
            w.contains("b.mzML differs from the first input a.mzML"),
            "{w}"
        );
        let w = compare(
            (HCD, Some(InstrumentType::QExactive)),
            (CID, Some(InstrumentType::QExactive)),
        )
        .unwrap()
        .unwrap();
        assert!(w.contains("dominant activation CID vs HCD"), "{w}");
        let w = compare((HCD, Some(InstrumentType::QExactive)), (None, None))
            .unwrap()
            .unwrap();
        assert!(w.contains("analyzer not recorded vs"), "{w}");
    }

    #[test]
    fn mgf_has_no_acquisition_metadata() {
        assert_eq!(
            detect_acquisition_metadata(Path::new("no-such-file.mgf")),
            (None, None)
        );
    }
}

#[cfg(test)]
mod require_spectra_tests {
    use super::{require_spectra, with_first_read_error};

    #[test]
    fn the_first_read_error_is_appended() {
        assert_eq!(with_first_read_error("none.".into(), &[]), "none.");
        let m = with_first_read_error("none.".into(), &["bad table".into(), "x".into()]);
        assert_eq!(m, "none. First read error: bad table");
    }
    use std::path::Path;

    #[test]
    fn a_file_with_spectra_passes() {
        assert!(require_spectra(Path::new("run.raw"), 2, 1).is_ok());
    }

    #[test]
    fn an_empty_file_is_an_error_naming_the_file_and_format() {
        for (file, format) in [
            ("dir/run.raw", "Thermo .raw"),
            ("dir/run.mzML", "mzML"),
            ("dir/run.mzML.gz", "mzML"),
            ("dir/run.mgf", "MGF"),
            ("dir/run.d", "Bruker .d"),
        ] {
            let err = require_spectra(Path::new(file), 2, 0).unwrap_err();
            assert!(err.contains(file) && err.contains(format), "{err}");
            assert!(err.contains("MS2"), "{err}");
        }
        let err = require_spectra(Path::new("run.mzML"), 3, 0).unwrap_err();
        assert!(err.contains("MS3"), "{err}");
    }
}

#[cfg(test)]
mod d_index_tests {
    use super::{d_load_by_index, d_metadata, Spectrum};
    use std::collections::HashSet;

    fn items() -> impl Iterator<Item = Result<Spectrum, String>> {
        (1..=5).map(|z| {
            if z == 2 {
                return Err("frame could not be decoded".to_string());
            }
            Ok(Spectrum {
                precursor_charge: Some(z),
                peaks: vec![(100.0, 1.0); z as usize],
                ..Default::default()
            })
        })
    }

    #[test]
    fn both_pre_passes_number_the_readable_spectra_alike() {
        let meta = d_metadata(items(), usize::MAX);
        let charges: Vec<Option<i32>> = meta.iter().map(|m| m.precursor_charge).collect();
        assert_eq!(charges, vec![Some(1), Some(3), Some(4), Some(5)]);
        let all: HashSet<usize> = (0..meta.len()).collect();
        let loaded = d_load_by_index(items(), &all, usize::MAX);
        assert_eq!(loaded.len(), meta.len());
        for (idx, m) in meta.iter().enumerate() {
            assert_eq!(loaded[&idx].precursor_charge, m.precursor_charge);
            assert_eq!(loaded[&idx].peaks.len(), m.num_peaks);
        }
    }

    #[test]
    fn the_bench_cap_counts_readable_spectra() {
        let meta = d_metadata(items(), 2);
        assert_eq!(meta.len(), 2);
        let loaded = d_load_by_index(items(), &HashSet::from([1]), 2);
        assert_eq!(loaded[&1].precursor_charge, Some(3));
        assert!(d_load_by_index(items(), &HashSet::from([2]), 2).is_empty());
    }
}

#[cfg(test)]
mod format_routing_tests {
    use super::input_format_flags;
    use std::path::Path;

    #[test]
    fn gz_is_routed_by_the_underlying_extension() {
        assert_eq!(
            input_format_flags(Path::new("x/foo.mzML.gz")),
            (true, false, false, false)
        );
        assert_eq!(
            input_format_flags(Path::new("foo.MGF.GZ")),
            (false, false, false, true)
        );
        assert_eq!(
            input_format_flags(Path::new("foo.mzML")),
            (true, false, false, false)
        );
        assert_eq!(
            input_format_flags(Path::new("foo.raw")),
            (false, true, false, false)
        );
        assert_eq!(
            input_format_flags(Path::new("run.d")),
            (false, false, true, false)
        );
        assert_eq!(
            input_format_flags(Path::new("foo.mgf")),
            (false, false, false, true)
        );
    }
}

#[cfg(test)]
mod title_tests {
    use super::prefix_spectrum_titles;
    use model::Spectrum;

    fn titled(title: &str, scan: Option<i32>) -> Spectrum {
        Spectrum {
            title: title.to_string(),
            scan,
            ..Default::default()
        }
    }

    #[test]
    fn titles_get_the_prefix_and_untitled_spectra_their_scan() {
        let mut spectra = vec![titled("Scan 42", Some(42)), titled("", Some(7))];
        prefix_spectrum_titles(&mut spectra, "run/");
        assert_eq!(spectra[0].title, "run/Scan 42");
        assert_eq!(spectra[1].title, "run/scan=7");
    }

    #[test]
    fn column_and_row_separators_in_a_title_become_spaces() {
        let mut spectra = vec![titled("a\tb\r\nc", Some(1))];
        prefix_spectrum_titles(&mut spectra, "ru\tn/");
        assert_eq!(spectra[0].title, "ru n/a b  c");
    }
}
