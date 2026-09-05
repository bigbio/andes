//! Spectrum-file loading, format routing, metadata scanning and precursor calibration.

use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::mpsc::SyncSender;

use input::{detect_instrument_type, MgfReader, MzMLReader};
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
    let is_d = matches!(ext.as_deref(), Some("d"));
    let is_mgf = !is_mzml && !is_raw && !is_d;
    (is_mzml, is_raw, is_d, is_mgf)
}

/// Prefix spectrum titles so pooled multi-file PIN SpecIds stay unique.
/// Decide the per-file SpecId/title prefix. Returns `None` for a single-file
/// search so its PIN output stays byte-identical to the pre-multi-file path;
/// returns `Some("<stem>/")` only when disambiguating across multiple inputs.
pub(crate) fn title_prefix_for(num_files: usize, file_stem: &str) -> Option<String> {
    (num_files > 1).then(|| format!("{file_stem}/"))
}

pub(crate) fn prefix_spectrum_titles(chunk: &mut [Spectrum], prefix: &str) {
    for spec in chunk.iter_mut() {
        if spec.title.is_empty() {
            spec.title = format!("{prefix}scan={}", spec.scan.unwrap_or(0));
        } else {
            spec.title = format!("{prefix}{}", spec.title);
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
    if is_mzml {
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
    if is_mzml {
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

/// Auto-detect an isobaric label (TMT/iTRAQ) by sampling the first `SAMPLE_N`
/// MS2 spectra and inspecting their reporter-ion region. Used only when
/// `--protocol auto` is left at its default, to engage the isobaric windowed
/// peak filter with zero config.
///
/// Returns `None` for `.raw`/`.d` (the sampling reader here is mzML/MGF only —
/// the protocol then stays as-is, byte-identical) and for label-free data, so
/// non-isobaric runs are unchanged. The mzML benchmark datasets (Astral, UPS1,
/// TMT) all flow through the mzML branch.
pub(crate) fn detect_isobaric_sampled(
    path: &Path,
    is_mzml: bool,
    is_mgf: bool,
    ms_level: u32,
    high_res: bool,
) -> Option<input::IsobaricLabel> {
    const SAMPLE_N: usize = 1000;
    if !(is_mzml || is_mgf) {
        return None;
    }
    let indices: HashSet<usize> = (0..SAMPLE_N).collect();
    let loaded = load_spectra_by_index(path, is_mzml, ms_level, &indices, usize::MAX).ok()?;
    let sample: Vec<Spectrum> = loaded.into_values().collect();
    input::detect_isobaric(&sample, high_res)
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
    prepared: &PreparedSearch<'_>,
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

/// Warn BEFORE scoring when the in-RAM candidate index is unlikely to fit.
///
/// Without this, a large database (a whole human proteome at three missed cleavages)
/// runs for half an hour and is then killed by the OOM killer with no message from
/// andes at all — the user sees only a dead process and no output. Measured on a
/// 20k-protein human FASTA: ~0.65 KB resident per candidate for a plain search and
/// ~0.92 KB under `--glyco`, which holds per-spectrum glyco state on top.
///
/// This warns rather than aborts: the estimate is a linear fit, machines differ, and
/// refusing to start a run that would have succeeded is worse than a noisy warning.
/// Only Linux exposes MemAvailable cheaply; elsewhere the check is skipped.
pub(crate) fn warn_if_index_will_not_fit(n_candidates: usize, glyco: bool) {
    const BYTES_PER_CANDIDATE_PLAIN: f64 = 665.0;
    const BYTES_PER_CANDIDATE_GLYCO: f64 = 940.0;

    let available = match std::fs::read_to_string("/proc/meminfo") {
        Ok(text) => text
            .lines()
            .find_map(|l| l.strip_prefix("MemAvailable:"))
            .and_then(|v| v.split_whitespace().next()?.parse::<u64>().ok())
            .map(|kb| kb * 1024),
        Err(_) => None,
    };
    let Some(available) = available else { return };

    let per = if glyco {
        BYTES_PER_CANDIDATE_GLYCO
    } else {
        BYTES_PER_CANDIDATE_PLAIN
    };
    let estimate = (n_candidates as f64 * per) as u64;
    if estimate <= available {
        return;
    }
    let gb = |b: u64| b as f64 / (1024.0 * 1024.0 * 1024.0);
    eprintln!(
        "WARNING: this search needs roughly {:.1} GB for the in-RAM candidate index \
         ({} candidates) but only {:.1} GB is available. The process is likely to be \
         killed by the operating system partway through, with no result written.",
        gb(estimate),
        n_candidates,
        gb(available)
    );
    if glyco {
        eprintln!(
            "  --glyco cannot use the out-of-core index yet (--candidate-index mmap is \
             rejected in glyco mode), so reduce the search space instead: pass \
             --max-missed-cleavages 1 or 2 (glyco defaults to 3), restrict the FASTA to \
             the proteins of interest, or split the database and merge the .glyco.pin \
             files afterwards."
        );
    } else {
        eprintln!("  Re-run with --candidate-index mmap to page the index from disk instead.");
    }
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
             across {} peeked MS2 spectra). Using dominant = {} ({}/{}); other methods \
             present: {}",
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
