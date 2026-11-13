//! Spectrum-file loading, format routing, metadata scanning and precursor calibration.

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

/// Statistics returned by the parser-thread helper.
#[derive(Debug, Default)]
pub(crate) struct ParseStats {
    pub(crate) error_count: usize,
    pub(crate) first_errors: Vec<String>,
}

/// Lowercased spectrum-file extension with a trailing `.gz` stripped, so
/// `run.mzML.gz` reports `mzml` rather than `gz`.
///
/// `Path::extension` returns `gz` for a double extension, which silently
/// defeated the `== "mzml"` guards on the metadata-detection helpers below:
/// a gzipped mzML skipped instrument and activation detection entirely and
/// fell back to the low-res default model. The readers those guards protect
/// all use `open_buf_maybe_gz`, so the guard -- not the reader -- was the
/// limitation. Mirrors the `.gz` handling in `input_format_flags`.
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

/// Whether `path` names an archived Bruker `.d` folder (`X.d.zip`,
/// `X.d.tar`, `X.d.tar.gz`, `X.d.tgz`, `X.d.gz`; case-insensitive).
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
    // Strip a trailing `.gz` so the format is detected from the underlying
    // extension (`spectra.mzML.gz` → mzML, `spectra.mgf.gz` → MGF). `.raw`/`.d`
    // are binary/directory and never gzipped, so this only affects mzML/MGF.
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

/// The per-file SpecId/title prefix, `<stem>/`. Every search applies it, also
/// with one input: native spectrum titles repeat across files, and PINs from
/// separate per-file searches are pooled into one Percolator run, so the
/// SpecIds must be unique across runs, not only within one search.
pub(crate) fn title_prefix_for(file_stem: &str) -> Option<String> {
    Some(format!("{file_stem}/"))
}

/// Put `prefix` before every title; a spectrum without a title is named by
/// its scan. The title is the `SpecId` of the tab-separated outputs (PIN,
/// TSVs), so a tab or a line break in it becomes a space.
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

/// Fold the spectra a tolerant mzML reader skipped (malformed or in an
/// unsupported binary encoding) into `stats`, so the run reports them.
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

/// `msg` followed by the first read error the parser recorded, when there is one.
pub(crate) fn with_first_read_error(msg: String, first_errors: &[String]) -> String {
    match first_errors.first() {
        Some(e) => format!("{msg} First read error: {e}"),
        None => msg,
    }
}

/// An input file that yielded no spectra at the searched MS level is an error,
/// not an empty search: a search of zero spectra would only run on the default
/// acquisition parameters and write an empty PIN. `found` counts the spectra
/// the file contributed.
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

/// Producer helper: drains `reader` into fixed-size chunks of `Spectrum`
/// and sends them through `tx`. Stops at `bench_cap` total spectra (or
/// `usize::MAX` for unbounded). Parse errors are counted and the first few
/// captured for downstream reporting; the channel is closed when the
/// reader is exhausted or the consumer hangs up.
///
/// Generic over the reader's error type so the same helper serves both
/// MGF and mzML.
///
/// Runs on a dedicated thread so chunk N+1 is PARSED while chunk N is SCORED.
/// Channel capacity is 2 (one in-flight + one queued) so the producer stays at
/// most one chunk ahead.
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
                    // If the consumer hung up, stop. Sender is moved into the
                    // function, so dropping returns `Err(SendError(chunk))`.
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

/// Sends `spectra` (file order; `ms2_to_ms1` indexes `ms1`) through `tx` in
/// chunks of `chunk_size`, each with the [`Ms1Link`]
/// [`MzMLReader::read_with_ms1_chunked`] gives it: the peaks of every MS1 scan
/// the chunk links to, once each, in file order. Stops when the consumer hangs
/// up.
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

/// Lightweight metadata collected in one linear file scan for precursorCal.
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

/// The readable spectra of `items`, numbered as the search reader counts
/// them: an unreadable spectrum is skipped and takes no index, and
/// `bench_cap` counts readable spectra only.
fn readable_spectra<E>(
    items: impl Iterator<Item = Result<Spectrum, E>>,
    bench_cap: usize,
) -> impl Iterator<Item = (usize, Spectrum)> {
    items.filter_map(Result::ok).take(bench_cap).enumerate()
}

/// Metadata of the readable `.d` spectra, indexed as [`readable_spectra`].
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

/// The readable `.d` spectra at `indices`, indexed as [`readable_spectra`]
/// (the same numbering as [`d_metadata`]).
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

/// The MS2 spectra of a Bruker `.d`, in the order the search reads them.
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
    // Metadata only: no placeholder peak lists, which on a large file were a
    // transient allocation of one tuple per peak for every spectrum.
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

/// Auto-detect an isobaric label (TMT/TMTpro/iTRAQ 4/8) by sampling the first
/// `SAMPLE_N` spectra at `ms_level` and inspecting their reporter-ion region.
///
/// An SPS-MS3 run carries its reporters in the MS3 scans, not the (ion-trap)
/// MS2 it searches, so when an MS2 sample shows no label the MS3 scans of an
/// mzML or `.raw` are sampled too (high-res: SPS-MS3 reporters are read in the
/// Orbitrap). Returns `None` for label-free data.
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
        eprintln!("Protocol resolver: isobaric reporters found in the MS3 scans (SPS-MS3)");
    }
    label
}

/// Whether the MS3 scans of an mzML or Thermo `.raw` carry isobaric reporter
/// ions (an SPS-MS3 run).
pub(crate) fn ms3_carries_reporters(path: &Path, is_mzml: bool, is_raw: bool) -> bool {
    if !(is_mzml || is_raw) {
        return false;
    }
    let Some(ms3) = sample_spectra(path, is_mzml, false, is_raw, 3) else {
        return false;
    };
    let found = input::detect_isobaric(&ms3, true).is_some();
    if !found {
        eprintln!(
            "Protocol resolver: no isobaric reporters in a sample of {} MS3 scans",
            ms3.len()
        );
    }
    found
}

/// Collision energy, isolation width and FAIMS CV of the first MS2 scans of an
/// mzML, Thermo `.raw` or Bruker `.d`. Empty for MGF and for a file that cannot
/// be read.
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

/// The first spectra at `ms_level` of an mzML, MGF, Thermo `.raw` or Bruker
/// `.d` (MS2 only; a `.d` has no other level with a precursor).
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
            // MS3 scans are sampled for their reporter ions only, so they are
            // kept without a precursor m/z; MS2 goes through the search path.
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
        eprintln!(
            "Precursor mass calibration skipped ({} SpecKeys < {} threshold; elapsed: {:.2}s). \
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
        eprintln!(
            "Precursor mass shift learned: {:.3} ppm from {} confident PSMs (robust sigma {:.3} ppm; elapsed: {:.2}s)",
            stats.shift_ppm,
            stats.confident_psm_count,
            stats.robust_sigma_ppm,
            t_cal.elapsed().as_secs_f64()
        );
    } else {
        eprintln!(
            "Precursor mass calibration skipped (insufficient confident PSMs: {} with PSMs, {} below RawScore floor, {} failed |residual|>50ppm; elapsed: {:.2}s)",
            stats.queues_with_psm,
            stats.rejected_low_score,
            stats.rejected_residual,
            t_cal.elapsed().as_secs_f64()
        );
    }
    Ok(stats)
}

/// Peek the spectrum file and return the dominant
/// `ActivationMethod` across the first several MS2 spectra.
///
/// Reads up to `MAX_PEEK` spectra (early-exit) and tallies a histogram of
/// activation methods. Returns the most-common method, or `None` when no
/// spectra carry an activation cvParam (older mzMLs, MGF, etc.).
///
/// Currently only mzML files (`.mzml` / `.mzML` extension) carry an
/// `<activation>` block. For anything else (MGF, unknown extension) we
/// return `None` and the caller falls back to the historical default.
///
/// When multiple activation methods are present, prints a single
/// `eprintln!` warning naming the runner-up and its count.
pub(crate) fn detect_dominant_activation(
    spectrum_path: &std::path::Path,
) -> Option<ActivationMethod> {
    // Only mzML carries `<activation>`. Other formats: caller falls back.
    let ext_lower = spectrum_ext_lower(spectrum_path);
    if ext_lower.as_deref() != Some("mzml") {
        return None;
    }

    const MAX_PEEK: usize = 64;

    let reader = MzMLReader::new(input::open_buf_maybe_gz(spectrum_path).ok()?);

    // Tally counts keyed by ActivationMethod variant.
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

    // Find the dominant method. Ties are broken by ActivationMethod's
    // declaration order via match below, which is stable.
    let dominant = counts
        .iter()
        // Deterministic on ties: HashMap iteration order is randomised per
        // process, so a bare `max_by_key(count)` picks an ARBITRARY maximum. A
        // 1:1 interleaved HCD/ETD acquisition ties exactly at 32/32 over the
        // peeked window and flipped the selected model run-to-run on identical
        // input. Tie-break on the discriminant so the choice is reproducible.
        .max_by_key(|(&m, &n)| (n, std::cmp::Reverse(m as u8)))
        .map(|(&m, _)| m)?;

    // Warn on mixed activation. The dominant method still wins; this is
    // purely informational so the user can spot heterogeneous mzMLs.
    if counts.len() > 1 {
        let mut other_pairs: Vec<(ActivationMethod, usize)> = counts
            .iter()
            .filter(|(&m, _)| m != dominant)
            .map(|(&m, &n)| (m, n))
            .collect();
        other_pairs.sort_by(|a, b| b.1.cmp(&a.1));
        let total: usize = counts.values().sum();
        let dominant_count = counts[&dominant];
        eprintln!(
            "Param resolver: mixed activation methods in input ({} different methods \
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

/// Helper to call `input::detect_instrument_type` on an mzML path.
///
/// Mirrors the structure of `detect_dominant_activation` so the two
/// detection passes look symmetric at the call site. Returns `None` for
/// non-mzML inputs or when the mzML has no recoverable instrument metadata.
pub(crate) fn detect_instrument_type_for_path(
    spectrum_path: &std::path::Path,
) -> Option<InstrumentType> {
    let ext_lower = spectrum_ext_lower(spectrum_path);
    if ext_lower.as_deref() != Some("mzml") {
        return None;
    }

    detect_instrument_type(input::open_buf_maybe_gz(spectrum_path).ok()?)
}

/// Dominant activation and analyzer of one input, each read from the metadata
/// on its own: an mzML or a `.raw` whose scans name no activation still names
/// its analyzer. MGF and unreadable inputs give `(None, None)`.
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
        // timsTOF DDA-PASEF: CID fragmentation on a TOF analyzer.
        (Some(ActivationMethod::CID), Some(InstrumentType::TimsTOF))
    } else {
        (None, None)
    }
}

/// Compare a further `--spectrum` input with the first one, whose metadata
/// selects the model, the fragment tolerance and the score features of the
/// whole search.
///
/// `Err` when the two analyzers are of different resolution classes: the
/// run-level tolerance would be wrong for one of the files. `Ok(Some(_))` is a
/// warning for a difference the search tolerates (another analyzer of the same
/// class, another dominant activation, metadata missing on one side).
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
        // Also without an activation on either side.
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
        // Metadata on one side only cannot show a class conflict.
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

    /// Five spectra (charge 1..=5, one peak per unit of charge) with the
    /// second one unreadable.
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

    // A gzipped spectrum is read transparently (input::open_maybe_gz) and must be
    // routed by its UNDERLYING extension, not the bare `.gz` (finding 2.5).
    #[test]
    fn gz_is_routed_by_the_underlying_extension() {
        // (is_mzml, is_raw, is_d, is_mgf)
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
