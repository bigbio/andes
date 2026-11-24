//! Scan lookups: run stems, native ids, the survey scan of an MS2 and its isolation window.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use model::scan::ScanContext;
use model::spectrum::Spectrum;
use model::stats::{median_f32_mut, Even};
use model::Tolerance;
use quant::isobaric::Plex;
use quant::ms1::Ms1RunIndex;
use quant::unmix::{window_record, Survey, WindowRecord};

use super::collect::RunCapture;
use super::REPORTER_REGION_MARGIN;

pub(crate) fn run_stem(path: &Path) -> String {
    let name = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string());
    let name = name.strip_suffix(".gz").map(str::to_string).unwrap_or(name);
    match name.rfind('.') {
        Some(i) if i > 0 => name[..i].to_string(),
        _ => name,
    }
}

pub(crate) fn check_unique_stems(paths: &[PathBuf]) -> Result<(), String> {
    let mut seen: HashMap<String, &Path> = HashMap::new();
    for p in paths {
        if let Some(first) = seen.insert(run_stem(p), p) {
            return Err(format!(
                "--quant-dir: inputs {} and {} share the run name {:?}; the per-run artifacts \
                 are named by it, so rename one of the files",
                first.display(),
                p.display(),
                run_stem(p)
            ));
        }
    }
    Ok(())
}

pub(super) fn native_id<'a>(run: &RunCapture, title: &'a str) -> &'a str {
    match &run.title_prefix {
        Some(p) => title.strip_prefix(p.as_str()).unwrap_or(title),
        None => title,
    }
}

pub(super) fn survey_scan(
    ms1: &Ms1RunIndex,
    context: Option<&ScanContext>,
    rt: f64,
) -> Option<usize> {
    const MAX_BACKTRACK: usize = 64;
    if let Some(i) = context
        .and_then(|c| c.parent_id.as_deref())
        .and_then(|id| ms1.position_of(id))
    {
        return Some(i);
    }
    let last = ms1.scan_at_or_before(rt)?;
    let Some(cv) = context.and_then(|c| c.faims_cv) else {
        return Some(last);
    };
    (last.saturating_sub(MAX_BACKTRACK)..=last)
        .rev()
        .find(|&i| ms1.scan(i).faims_cv == Some(cv))
}

pub(super) fn covers_reporters(
    scan_window: Option<(f64, f64)>,
    plex: Plex,
    tol: Tolerance,
) -> bool {
    let Some((lo, hi)) = scan_window else {
        return true;
    };
    let (first, last) = plex
        .channels()
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), c| {
            (a.min(c.mz), b.max(c.mz))
        });
    lo <= first - tol.as_da(first) && hi >= last + tol.as_da(last)
}

pub(super) fn next_survey(ms1: &Ms1RunIndex, after: usize, faims_cv: Option<f32>) -> Option<usize> {
    const MAX_FORWARD: usize = 64;
    (after + 1..ms1.len().min(after + 1 + MAX_FORWARD))
        .find(|&i| faims_cv.is_none() || ms1.scan(i).faims_cv == faims_cv)
}

pub(super) fn isolation_window_record(
    ms1: &Ms1RunIndex,
    context: Option<&ScanContext>,
    spec: &Spectrum,
    mono_mz: f64,
    charge: u8,
) -> Option<WindowRecord> {
    let rt = spec.rt_seconds?;
    let centre = spec.isolation_target_mz.unwrap_or(spec.precursor_mz);
    let lo = centre - spec.isolation_lower_offset?;
    let hi = centre + spec.isolation_upper_offset?;
    let survey = |i: usize| Survey {
        rt: ms1.rt(i),
        peaks: ms1.peaks(i),
        detection_limit: ms1.detection_limit(i),
    };
    let prev = survey_scan(ms1, context, rt);
    let cv = context.and_then(|c| c.faims_cv);
    let next = match prev {
        Some(p) => next_survey(ms1, p, cv),
        None => (0..ms1.len().min(64)).find(|&i| cv.is_none() || ms1.scan(i).faims_cv == cv),
    };
    window_record(
        prev.map(survey),
        next.map(survey),
        rt,
        lo,
        hi,
        mono_mz,
        charge,
    )
}

pub(super) fn peptide_ion_current(peaks: &[(f64, f32)], plex: Plex) -> f32 {
    let top = plex
        .channels()
        .iter()
        .fold(f64::NEG_INFINITY, |a, c| a.max(c.mz))
        + REPORTER_REGION_MARGIN;
    peaks
        .iter()
        .filter(|p| p.0 > top && p.1.is_finite())
        .map(|p| p.1.max(0.0))
        .sum()
}

pub(super) fn median_f32(v: &[f32]) -> f32 {
    median_f32_mut(&mut v.to_vec(), Even::Upper).unwrap_or(f32::NAN)
}
