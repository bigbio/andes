//! Parameter grid search (`--optimize`).
//!
//! Scores a representative sample of the run under a grid of search parameters
//! (precursor tolerance x top-N x charge range), rescores each grid's PIN with
//! the native (Percolator-free) rescorer, keeps the combination with the most
//! target PSMs at `q <= 0.01`, and applies it to the whole run. This is the
//! parameter-optimization step of MSFragger's open search pipeline
//! (Yu et al., *Nat. Commun.* 2020), scoped to the parameters andes does not
//! derive from the data:
//!
//! * the **fragment** tolerance is not searched: the per-run fragment-error
//!   calibration estimates it from the confident matches, and a count-ranked
//!   grid would only push it wider;
//! * the precursor ppm calibration (`--precursor-cal`) is a separate, optional
//!   step and is not part of the grid.
//!
//! The grid is **staged**: the tolerance first (at the resolved top-N and
//! charge range), then top-N and charge range at the chosen tolerance. That
//! keeps the sample searches to ~14 instead of the full cross product, which
//! matters because a sample large enough to rank tolerances reliably (thousands
//! of spectra across the file) makes each search expensive.
//!
//! The ranking uses the native rescorer, not raw target-decoy counts: it sees
//! the same PIN features and competition groups as Percolator. On a tie the
//! resolved value wins, so the grid never changes a setting it cannot separate.

use std::ops::RangeInclusive;
use std::path::Path;

use model::spectrum::Spectrum;
use model::tolerance::{PrecursorTolerance, Tolerance};
use scoring_crate::RankScorer;
use search::{PreparedSearch, SearchIndex, SearchParams};

/// Precursor tolerances scored, in ppm (symmetric). Only used when the
/// resolved tolerance is in ppm; a Da tolerance is kept as is.
const TOLERANCES_PPM: [f64; 4] = [5.0, 10.0, 20.0, 50.0];
/// Top-N PSMs per spectrum scored. The resolved value joins the list so a tie
/// never changes it.
const TOP_NS: [u32; 4] = [1, 2, 3, 5];
/// Target q the grid counts identifications at.
const OPTIMIZE_Q: f64 = 0.01;
/// The native rescorer's seed, fixed so the grid is reproducible.
const OPTIMIZE_SEED: u64 = 42;

/// Run the grid on `sample` and overwrite the grid dimensions of `params` with
/// the best combination. Prints the whole grid and the choice.
pub(crate) fn optimize_parameters(
    idx: &SearchIndex,
    scorer: &RankScorer,
    fragment_tol_da: f64,
    decoy_prefix: &str,
    index_path: &Path,
    params: &mut SearchParams,
    sample: &[Spectrum],
) {
    if sample.is_empty() {
        eprintln!("optimize: the sample is empty; keeping the resolved parameters");
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
            eprintln!(
                "optimize: cannot prepare the sample search: {e}; keeping the resolved parameters"
            );
            return;
        }
    };
    eprintln!("optimize: {} spectra sampled from the run", sample.len());

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
    // ── Stage 1: precursor tolerance at the resolved top-N / charge range ──
    // Skipped under --open-search, where the precursor tolerance IS the mass
    // offset window the user asked for; and for a Da tolerance, which a ppm
    // grid cannot express.
    let mut best_tol: Option<BestTol> = None;
    if base.open_search.is_none()
        && matches!(base.precursor_tolerance.left, Tolerance::Ppm(_))
        && matches!(base.precursor_tolerance.right, Tolerance::Ppm(_))
    {
        for &tol in &TOLERANCES_PPM {
            let mut p = base.clone();
            p.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(tol));
            let Some(targets) = evaluate(&mut prepared, sample, idx, &p) else {
                continue;
            };
            rows.push(Row { stage: "tol", tol: Some(tol), top_n: resolved_top_n, range: resolved_range.clone(), targets });
            let key = (targets, p.precursor_tolerance == base.precursor_tolerance);
            if best_tol.as_ref().is_none_or(|b| key > b.key) {
                best_tol = Some(BestTol { key, params: p, tol });
            }
        }
    }
    let (tol_label, stage1_params) = match best_tol {
        Some(b) => (b.tol, b.params),
        None => (f64::NAN, base.clone()),
    };

    // ── Stage 2: top-N x charge range at the chosen tolerance ──
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
            let Some(targets) = evaluate(&mut prepared, sample, idx, &p) else {
                continue;
            };
            rows.push(Row { stage: "topN", tol: None, top_n, range: range.clone(), targets });
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

    eprintln!("  stage  tol(ppm) topN charge   targets");
    for row in &rows {
        eprintln!(
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
        eprintln!("optimize: no combination could be scored; keeping the resolved parameters");
        return;
    };
    eprintln!(
        "optimize: chosen --precursor-tol {}{} --top-n {top_n} --charge {}..={} \
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

/// Target PSMs at `q <= OPTIMIZE_Q` for one grid point: write the sample's PIN
/// and run the native rescorer over it.
fn evaluate(
    prepared: &mut PreparedSearch,
    sample: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
) -> Option<usize> {
    let queues = prepared.run_chunk_with_params(sample, 0, params);
    // Resolve this run's candidate indices before the next one replaces the
    // materialized list.
    prepared.sync_materialized_candidates();
    let mut buf: Vec<u8> = Vec::new();
    output::write_pin_to(
        &mut buf,
        sample,
        &queues,
        &prepared.candidates,
        params,
        idx,
    )
    .ok()?;
    let text = String::from_utf8_lossy(&buf);
    let rescored = crate::rescore::native_rescore_pin_all(&text, OPTIMIZE_SEED).ok()?;
    Some(
        rescored
            .targets
            .values()
            .filter(|p| p.q_value <= OPTIMIZE_Q)
            .count(),
    )
}
