use std::io::Write;
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};

use model::spectrum::Spectrum;
use model::tolerance::{PrecursorTolerance, Tolerance};
use scoring_crate::RankScorer;
use search::{PreparedSearch, SearchIndex, SearchParams};

const TOLERANCES_PPM: [f64; 4] = [5.0, 10.0, 20.0, 50.0];
const TOP_NS: [u32; 4] = [1, 2, 3, 5];
const OPTIMIZE_Q: f64 = 0.01;

pub(crate) struct SampleRescorer {
    pub(crate) backend: output::PercolatorBackend,
    pub(crate) args: Vec<String>,
    pub(crate) dir: PathBuf,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn optimize_parameters(
    idx: &SearchIndex,
    scorer: &RankScorer,
    fragment_tol_da: f64,
    decoy_prefix: &str,
    index_path: &Path,
    params: &mut SearchParams,
    sample: &[Spectrum],
    rescorer: &SampleRescorer,
) {
    if sample.is_empty() {
        log::warn!("the sample is empty; keeping the resolved parameters");
        return;
    }
    let base = params.clone();
    let resolved_top_n = base.top_n_psms_per_spectrum;
    let resolved_range = base.charge_range.clone();
    let mut prepared = match PreparedSearch::prepare_mmap(
        idx,
        &base,
        scorer,
        fragment_tol_da,
        decoy_prefix,
        index_path,
    ) {
        Ok(p) => p,
        Err(e) => {
            log::warn!(
                "cannot prepare the sample search: {e}; keeping the resolved parameters"
            );
            return;
        }
    };
    log::info!("{} spectra sampled from the run", sample.len());

    struct Row {
        stage: &'static str,
        tol: Option<f64>,
        top_n: u32,
        range: RangeInclusive<u8>,
        targets: usize,
    }
    let mut rows: Vec<Row> = Vec::new();

    struct BestTol {
        key: (usize, bool),
        params: SearchParams,
        tol: f64,
    }
    struct Best {
        key: (usize, bool, bool),
        params: SearchParams,
        targets: usize,
        top_n: u32,
        range: RangeInclusive<u8>,
    }
    let mut best_tol: Option<BestTol> = None;
    if base.open_search.is_none()
        && matches!(base.precursor_tolerance.left, Tolerance::Ppm(_))
        && matches!(base.precursor_tolerance.right, Tolerance::Ppm(_))
    {
        for &tol in &TOLERANCES_PPM {
            let mut p = base.clone();
            p.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(tol));
            let Some(targets) = evaluate(&mut prepared, sample, idx, &p, rescorer) else {
                continue;
            };
            rows.push(Row {
                stage: "tol",
                tol: Some(tol),
                top_n: resolved_top_n,
                range: resolved_range.clone(),
                targets,
            });
            let key = (targets, p.precursor_tolerance == base.precursor_tolerance);
            if best_tol.as_ref().is_none_or(|b| key > b.key) {
                best_tol = Some(BestTol {
                    key,
                    params: p,
                    tol,
                });
            }
        }
    }
    let (tol_label, stage1_params) = match best_tol {
        Some(b) => (b.tol, b.params),
        None => (f64::NAN, base.clone()),
    };

    let mut top_ns: Vec<u32> = TOP_NS.to_vec();
    if !top_ns.contains(&resolved_top_n) {
        top_ns.push(resolved_top_n);
    }
    let mut ranges: Vec<RangeInclusive<u8>> = vec![resolved_range.clone()];
    if *resolved_range.start() != 2 || *resolved_range.end() != 3 {
        ranges.push(2..=3);
    }
    let mut best: Option<Best> = None;
    for &top_n in &top_ns {
        for range in &ranges {
            let mut p = stage1_params.clone();
            p.top_n_psms_per_spectrum = top_n;
            p.charge_range = range.clone();
            let Some(targets) = evaluate(&mut prepared, sample, idx, &p, rescorer) else {
                continue;
            };
            rows.push(Row {
                stage: "topN",
                tol: None,
                top_n,
                range: range.clone(),
                targets,
            });
            let key = (
                targets,
                top_n == resolved_top_n,
                *range.start() == *resolved_range.start() && *range.end() == *resolved_range.end(),
            );
            if best.as_ref().is_none_or(|b| key > b.key) {
                best = Some(Best {
                    key,
                    params: p,
                    targets,
                    top_n,
                    range: range.clone(),
                });
            }
        }
    }

    log::info!("  stage  tol(ppm) topN charge   targets");
    for row in &rows {
        log::info!(
            "  {:>5} {:>8} {:>4} {:>2}..={:<2} {:>8}",
            row.stage,
            row.tol.map_or("-".to_string(), |t| format!("{t:.0}")),
            row.top_n,
            row.range.start(),
            row.range.end(),
            row.targets
        );
    }
    let Some(Best {
        params: chosen,
        targets,
        top_n,
        range,
        ..
    }) = best
    else {
        log::warn!("no combination could be scored; keeping the resolved parameters");
        return;
    };
    log::info!(
        "chosen --precursor-tol {}{} --top-n {top_n} --charge {}..={} \
         ({targets} targets at q <= {OPTIMIZE_Q})",
        if tol_label.is_nan() {
            "resolved".to_string()
        } else {
            format!("{tol_label:.0} ppm")
        },
        if chosen.precursor_tolerance == base.precursor_tolerance {
            " (unchanged)"
        } else {
            ""
        },
        range.start(),
        range.end()
    );
    *params = chosen;
}

fn evaluate(
    prepared: &mut PreparedSearch,
    sample: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    rescorer: &SampleRescorer,
) -> Option<usize> {
    let queues = prepared.run_chunk_with_params(sample, 0, params);
    prepared.sync_materialized_candidates();
    let pin = rescorer.dir.join("optimize.pin");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&pin).ok()?);
    output::write_pin_to(&mut file, sample, &queues, &prepared.candidates, params, idx).ok()?;
    file.flush().ok()?;
    drop(file);
    let targets = output::run_percolator_pooled(&rescorer.backend, &pin, &rescorer.args).ok()?;
    Some(
        targets
            .values()
            .filter(|p| p.q_value <= OPTIMIZE_Q)
            .count(),
    )
}
