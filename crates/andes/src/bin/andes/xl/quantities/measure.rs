//! MS1 measurement of cross-link precursors: around the retention time of an identification,
//! or around the predicted retention time of a transfer, each next to a decoy feature that
//! the q-values are computed against. The width of the elution peaks of a run is estimated
//! from the run itself and sets the smoothing, the prior on the apex and the extraction window.

use std::collections::BTreeMap;

use model::hash::splitmix64;
use model::stats::{median_mut, Even};
use quant::decoys::{probes, Hypothesis, Probe};
use quant::envelope::{elution_sd_s, pattern_trace, quantify_envelope};
use quant::evidence::Prior;
use quant::lfq::{FeatureQuant, LfqParams, LfqTarget};
use quant::ms1::Ms1RunIndex;
use rayon::prelude::*;

use super::align::{LinearAnchor, LinearKey};
use super::library::Library;
use super::{Evidence, Row, Settings};
use crate::quant::run::RunMs1;

/// A precursor to look for in a run that did not identify it.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Transfer {
    pub(super) precursor: usize,
    pub(super) run: usize,
    /// The predicted retention time and its standard deviation, in seconds.
    pub(super) rt: f64,
    pub(super) sd_s: f64,
}

/// The elution peaks of one run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct PeakScale {
    /// The standard deviation of a typical elution peak, in seconds.
    pub(super) sd_s: f64,
    /// The median time between MS1 scans, in seconds.
    pub(super) step_s: f64,
}

/// The full width at half maximum of a Gaussian peak in standard deviations: 2 sqrt(2 ln 2).
const FWHM_PER_SD: f64 = 2.354_820_045_030_949;
/// The pattern signal is smoothed with a Gaussian of this many elution-peak standard
/// deviations: narrow enough to keep neighbouring peaks apart, wide enough to average the
/// scans of one peak.
const SMOOTHING_PEAK_SDS: f64 = 0.5;
/// `gaussian_smooth` takes the width that holds this many standard deviations of its kernel.
const SMOOTHING_WIDTH_SDS: f64 = 8.0;
/// The MS2 spectrum of an identification was taken while its precursor eluted, so the prior
/// on the apex around it has the full width at half maximum of a peak as standard deviation.
const IDENTIFIED_PRIOR_PEAK_SDS: f64 = FWHM_PER_SD;
/// Retention-time half-window, in seconds, in which the elution width is estimated around an
/// identification.
const WIDTH_WINDOW_S: f64 = 60.0;

/// What the measurement of the identifications yields.
pub(super) struct Identified {
    pub(super) rows: Vec<Row>,
    /// Per run, the apex retention time of each linear anchor that has a feature.
    pub(super) linear: Vec<BTreeMap<LinearKey, f64>>,
    /// Per run, the width of its elution peaks; `None` for a run without any identification.
    pub(super) scales: Vec<Option<PeakScale>>,
}

fn target_of(mono_mz: f64, charge: u8, envelope: &[f64], share: f64, rt: f64) -> LfqTarget {
    LfqTarget {
        mono_mz,
        charge,
        anchor_rt: rt,
        envelope: envelope.to_vec(),
        isotope_share: share,
    }
}

/// The extraction settings of a probe in a run of elution width `scale`.
fn params_at(settings: &Settings, scale: PeakScale, window_s: f64) -> LfqParams {
    let climb = (scale.sd_s / scale.step_s.max(f64::MIN_POSITIVE)).ceil();
    LfqParams {
        rt_window_s: window_s,
        gaussian_width_s: settings
            .params
            .gaussian_width_s
            .or(Some(SMOOTHING_WIDTH_SDS * SMOOTHING_PEAK_SDS * scale.sd_s)),
        max_apex_climb: if climb.is_finite() { climb as usize } else { 0 },
        ..settings.params.clone()
    }
}

/// The standard deviation of the prior on the apex of an identification and the extraction
/// half-window, both in seconds.
fn identified_prior(settings: &Settings, scale: PeakScale) -> (f64, f64) {
    let sd = IDENTIFIED_PRIOR_PEAK_SDS * scale.sd_s;
    let window = settings
        .rt_window_s
        .unwrap_or_else(|| half_window(sd, scale));
    (sd, window)
}

fn half_window(prior_sd_s: f64, scale: PeakScale) -> f64 {
    Prior {
        t_hat: 0.0,
        s: prior_sd_s,
        w_prior: scale.sd_s,
    }
    .half_window()
}

/// The width of the elution peaks of a run: the median over the pattern traces around its
/// identifications.
fn peak_scale(ms1: &Ms1RunIndex, targets: &[LfqTarget], settings: &Settings) -> Option<PeakScale> {
    let step_s = ms1.median_scan_spacing()?;
    let params = LfqParams {
        rt_window_s: WIDTH_WINDOW_S,
        ..settings.params.clone()
    };
    let mut widths: Vec<f64> = targets
        .par_iter()
        .filter_map(|t| {
            let (rt, signal) = pattern_trace(ms1, t, &params)?;
            elution_sd_s(&rt, &signal, step_s)
        })
        .collect();
    let sd_s = median_mut(&mut widths, Even::Upper)?;
    (sd_s.is_finite() && sd_s > 0.0).then_some(PeakScale { sd_s, step_s })
}

/// The feature at `probe`: the precursor's own m/z or a decoy's, around the probe's time.
fn feature_at(
    ms1: &Ms1RunIndex,
    target: &LfqTarget,
    probe: &Probe,
    settings: &Settings,
    scale: PeakScale,
    window_s: f64,
) -> Option<FeatureQuant> {
    let at = LfqTarget {
        anchor_rt: probe.t_hat,
        ..target.clone()
    };
    let params = params_at(settings, scale, window_s);
    quantify_envelope(ms1, &at, probe.mz, &params, probe.s)
}

/// The bit of a row's hash that picks the kind of decoy of a transfer; `quant::decoys::probes`
/// reads the three bits below it.
const DECOY_KIND_BIT: u32 = 3;

/// The decoy a row competes with: for an identification the same envelope at a shifted mass;
/// for a transfer that or, for every other one, the same mass at a shifted retention time.
fn decoy_of<'a>(row: &Row, all: &'a [Probe]) -> Option<&'a Probe> {
    let mass_shift = all
        .iter()
        .find(|p| matches!(p.kind, Hypothesis::MassShiftDecoy(_)));
    let hash = splitmix64(((row.precursor as u64) << 16) | row.run as u64);
    if row.evidence == Evidence::Transferred && (hash >> DECOY_KIND_BIT) & 1 == 1 {
        all.iter()
            .find(|p| p.kind == Hypothesis::RtShiftDecoy)
            .or(mass_shift)
    } else {
        mass_shift
    }
}

/// The target feature of `row` and the feature of its decoy.
fn measure_row(
    ms1: &Ms1RunIndex,
    lib: &Library,
    row: &Row,
    settings: &Settings,
    scale: PeakScale,
) -> (Option<FeatureQuant>, Option<FeatureQuant>) {
    let p = &lib.precursors[row.precursor];
    let target = target_of(
        p.mono_mz,
        p.charge,
        &p.envelope,
        p.isotope_share,
        row.anchor_rt,
    );
    let Some(run_rt) = ms1.rt_range() else {
        return (None, None);
    };
    let all = probes(
        row.precursor as u32,
        row.run as u16,
        p.mono_mz,
        p.charge.max(1),
        row.anchor_rt,
        row.anchor_sd_s,
        scale.sd_s,
        run_rt,
    );
    let at = |probe: Option<&Probe>| {
        probe.and_then(|probe| feature_at(ms1, &target, probe, settings, scale, row.window_s))
    };
    (
        at(all.iter().find(|p| p.kind == Hypothesis::Target)),
        at(decoy_of(row, &all)),
    )
}

/// The apex retention time of a linear anchor, when it has a feature whose isotope pattern
/// fits.
fn linear_apex(
    ms1: &Ms1RunIndex,
    anchor: &LinearAnchor,
    settings: &Settings,
    scale: PeakScale,
) -> Option<f64> {
    let (sd, window) = identified_prior(settings, scale);
    let target = target_of(
        anchor.mono_mz,
        anchor.charge,
        &anchor.envelope,
        1.0,
        anchor.rt,
    );
    let params = params_at(settings, scale, window);
    let f = quantify_envelope(ms1, &target, anchor.mono_mz, &params, sd)?;
    (f.cosine.is_finite() && f.cosine >= settings.params.min_cosine).then_some(f.apex_rt)
}

fn measure_rows(
    ms1: &Ms1RunIndex,
    lib: &Library,
    rows: &mut [Row],
    ids: &[usize],
    settings: &Settings,
    scale: PeakScale,
) {
    let measured: Vec<_> = ids
        .par_iter()
        .map(|&i| (i, measure_row(ms1, lib, &rows[i], settings, scale)))
        .collect();
    for (i, (feature, decoy)) in measured {
        rows[i].feature = feature;
        rows[i].decoy = decoy;
    }
}

fn rows_by_run(rows: &[Row], n_runs: usize) -> Vec<Vec<usize>> {
    let mut by_run: Vec<Vec<usize>> = vec![Vec::new(); n_runs];
    for (i, row) in rows.iter().enumerate() {
        by_run[row.run].push(i);
    }
    by_run
}

/// Measures every identification in its run, and in the same pass the linear anchors of the
/// run, whose apex retention times the alignment is fitted on.
pub(super) fn identified(
    lib: &Library,
    anchors: &[Vec<LinearAnchor>],
    runs: &[RunMs1<'_>],
    settings: &Settings,
) -> Result<Identified, String> {
    let mut rows: Vec<Row> = lib
        .identifications
        .iter()
        .map(|id| Row {
            precursor: id.precursor,
            run: id.run,
            evidence: Evidence::Identified,
            anchor_rt: id.rt,
            anchor_sd_s: 0.0,
            window_s: 0.0,
            scan: id.scan,
            csm_score: id.score,
            n_csms: id.n_csms,
            feature: None,
            decoy: None,
            score: None,
            q: None,
        })
        .collect();
    let by_run = rows_by_run(&rows, runs.len());
    let mut linear = vec![BTreeMap::new(); runs.len()];
    let mut scales = vec![None; runs.len()];
    for (r, run) in runs.iter().enumerate() {
        let ids = &by_run[r];
        let run_anchors = anchors.get(r).map_or(&[][..], Vec::as_slice);
        if ids.is_empty() && run_anchors.is_empty() {
            continue;
        }
        let t0 = std::time::Instant::now();
        let ms1 = run.load().map_err(|e| format!("--xl --lfq: {e}"))?;
        let mut targets: Vec<LfqTarget> = ids
            .iter()
            .map(|&i| {
                let p = &lib.precursors[rows[i].precursor];
                target_of(
                    p.mono_mz,
                    p.charge,
                    &p.envelope,
                    p.isotope_share,
                    rows[i].anchor_rt,
                )
            })
            .collect();
        targets.extend(
            run_anchors
                .iter()
                .map(|a| target_of(a.mono_mz, a.charge, &a.envelope, 1.0, a.rt)),
        );
        let Some(scale) = peak_scale(&ms1, &targets, settings) else {
            log::warn!(
                target: "xl-quant",
                "{}: no elution peak width could be estimated around its {} identifications, \
                 so none of them is quantified",
                run.stem,
                targets.len()
            );
            continue;
        };
        scales[r] = Some(scale);
        let (sd, window) = identified_prior(settings, scale);
        for &i in ids {
            rows[i].anchor_sd_s = sd;
            rows[i].window_s = window;
        }
        measure_rows(&ms1, lib, &mut rows, ids, settings, scale);
        linear[r] = run_anchors
            .par_iter()
            .filter_map(|a| Some((a.key, linear_apex(&ms1, a, settings, scale)?)))
            .collect();
        log::debug!(
            target: "xl-quant",
            "{}: elution peak sd {:.1} s at {:.2} s between MS1 scans; {} features in {} \
             identifications, {} of {} linear anchors with a feature [{:.2}s]",
            run.stem,
            scale.sd_s,
            scale.step_s,
            ids.iter().filter(|&&i| rows[i].feature.is_some()).count(),
            ids.len(),
            linear[r].len(),
            run_anchors.len(),
            t0.elapsed().as_secs_f64()
        );
    }
    Ok(Identified {
        rows,
        linear,
        scales,
    })
}

/// The row of a transfer in a run of elution width `scale`: the prior on the apex combines
/// the uncertainty of the predicted retention time with the width of an elution peak.
fn transfer_row(t: &Transfer, scale: PeakScale) -> Row {
    let sd = t.sd_s.hypot(scale.sd_s);
    Row {
        precursor: t.precursor,
        run: t.run,
        evidence: Evidence::Transferred,
        anchor_rt: t.rt,
        anchor_sd_s: sd,
        window_s: half_window(sd, scale),
        scan: 0,
        csm_score: f64::NAN,
        n_csms: 0,
        feature: None,
        decoy: None,
        score: None,
        q: None,
    }
}

/// Measures every transfer around its predicted retention time, in the runs whose elution
/// width is known.
pub(super) fn transferred(
    lib: &Library,
    transfers: Vec<Transfer>,
    scales: &[Option<PeakScale>],
    runs: &[RunMs1<'_>],
    settings: &Settings,
) -> Result<Vec<Row>, String> {
    let mut rows: Vec<Row> = transfers
        .iter()
        .filter_map(|t| Some(transfer_row(t, scales.get(t.run).copied().flatten()?)))
        .collect();
    let by_run = rows_by_run(&rows, runs.len());
    for (r, run) in runs.iter().enumerate() {
        let (ids, Some(scale)) = (&by_run[r], scales.get(r).copied().flatten()) else {
            continue;
        };
        if ids.is_empty() {
            continue;
        }
        let ms1 = run.load().map_err(|e| format!("--xl --lfq --mbr: {e}"))?;
        measure_rows(&ms1, lib, &mut rows, ids, settings, scale);
    }
    Ok(rows)
}

#[cfg(test)]
#[path = "tests/measure.rs"]
mod tests;
