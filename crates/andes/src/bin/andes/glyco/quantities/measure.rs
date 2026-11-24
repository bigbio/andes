//! MS1 measurement of glycopeptide targets: the extracted isotope envelope, identified and
//! transferred measurements, and the monoisotopic area.

use std::collections::{BTreeSet, HashMap};

use quant::lfq::{quantify, quantify_decoy, FeatureQuant, LfqParams, LfqTarget};
use rayon::prelude::*;

use super::{FeatureRow, Settings, Target, ISOTOPE_COVERAGE, MAX_ISOTOPES, MIN_ISOTOPES};
use crate::glyco::matching;
use crate::quant::run::RunMs1;

pub(super) fn extracted_envelope(full: &[f64], offset: usize) -> Vec<f64> {
    let tail = full.get(offset.min(full.len())..).unwrap_or(&[]);
    if tail.is_empty() {
        return vec![1.0];
    }
    let total: f64 = tail.iter().sum();
    let mut k = tail.len();
    let mut cum = 0.0;
    for (i, &v) in tail.iter().enumerate() {
        cum += v;
        if i + 1 >= MIN_ISOTOPES && cum >= ISOTOPE_COVERAGE * total {
            k = i + 1;
            break;
        }
    }
    let k = k.clamp(MIN_ISOTOPES, MAX_ISOTOPES).min(tail.len());
    let mut env: Vec<f64> = tail[..k].to_vec();
    let sum: f64 = env.iter().sum();
    if sum > 0.0 {
        for v in &mut env {
            *v /= sum;
        }
    }
    env
}

#[derive(Debug, Clone, Default)]
pub(super) struct Measurement {
    pub(super) feature: Option<FeatureQuant>,
    pub(super) decoy: Option<FeatureQuant>,
    pub(super) match_score: Option<f32>,
    pub(super) decoy_match: Option<f32>,
}

pub(super) fn measure(
    targets: &[Target],
    runs: &[RunMs1<'_>],
    settings: &Settings,
    phase: &str,
) -> Result<(usize, Vec<Measurement>), String> {
    measure_with(targets, runs, settings, phase, measure_identified)
}

pub(super) fn measure_identified(
    ms1: &quant::ms1::Ms1RunIndex,
    t: &Target,
    params: &LfqParams,
) -> Measurement {
    let lfq = t.lfq_target();
    let params = LfqParams {
        n_isotopes: t.envelope.len(),
        rt_window_s: t.window_s,
        ..params.clone()
    };
    let mut feature = quantify(ms1, &lfq, &params);
    if let Some(f) = feature.as_mut() {
        f.area = mono_area(f, &lfq.envelope);
    }
    Measurement {
        feature,
        decoy: quantify_decoy(ms1, &lfq, &params),
        match_score: None,
        decoy_match: None,
    }
}

fn measure_with(
    targets: &[Target],
    runs: &[RunMs1<'_>],
    settings: &Settings,
    phase: &str,
    f: impl Fn(&quant::ms1::Ms1RunIndex, &Target, &LfqParams) -> Measurement + Sync,
) -> Result<(usize, Vec<Measurement>), String> {
    let mut by_run: Vec<Vec<usize>> = vec![Vec::new(); runs.len()];
    for (i, t) in targets.iter().enumerate() {
        by_run[t.run_idx].push(i);
    }
    let mut slot: Vec<Option<Measurement>> = vec![None; targets.len()];
    let mut found = 0usize;
    for (ri, run) in runs.iter().enumerate() {
        if by_run[ri].is_empty() {
            continue;
        }
        let t0 = std::time::Instant::now();
        let ms1 = run.load().map_err(|e| format!("--glyco --lfq: {e}"))?;
        let part: Vec<(usize, Measurement)> = by_run[ri]
            .par_iter()
            .map(|&i| (i, f(&ms1, &targets[i], &settings.params)))
            .collect();
        let run_found = part.iter().filter(|(_, m)| m.feature.is_some()).count();
        found += run_found;
        log::info!(
            target: "glyco-quant",
            "{phase} {}: {} features in {} targets, {:.2}s",
            run.stem,
            run_found,
            part.len(),
            t0.elapsed().as_secs_f64()
        );
        for (i, m) in part {
            slot[i] = Some(m);
        }
    }
    Ok((
        found,
        slot.into_iter().map(Option::unwrap_or_default).collect(),
    ))
}

pub(super) struct MatchReferences {
    pub(super) target: matching::Reference,
    pub(super) decoy: Option<matching::Reference>,
}

pub(super) fn build_match_references(
    rows: &[FeatureRow],
    transfers: &[Target],
    runs: &[RunMs1<'_>],
    settings: &Settings,
) -> Result<HashMap<usize, MatchReferences>, String> {
    let mut by_run: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); runs.len()];
    for t in transfers {
        if let Some(d) = t.donor {
            by_run[rows[d].target.run_idx].insert(d);
        }
    }
    let mut refs: HashMap<usize, MatchReferences> = HashMap::new();
    for (ri, ids) in by_run.iter().enumerate() {
        if ids.is_empty() {
            continue;
        }
        let ms1 = runs[ri]
            .load()
            .map_err(|e| format!("--glyco --mbr match reference: {e}"))?;
        for &id in ids {
            let row = &rows[id];
            let Some(f) = &row.feature else {
                continue;
            };
            let t = &row.target;
            let build = |mono_mz: f64| {
                matching::Reference::build(
                    &ms1,
                    mono_mz,
                    t.charge,
                    t.envelope.len(),
                    t.anchor_rt,
                    t.window_s,
                    f.rt_start,
                    f.rt_stop,
                    settings.params.tol,
                )
            };
            let Some(target) = build(t.mono_mz) else {
                continue;
            };
            let decoy = build(t.mono_mz + settings.params.decoy_mz_shift);
            refs.insert(id, MatchReferences { target, decoy });
        }
    }
    Ok(refs)
}

pub(super) fn measure_transfers(
    targets: &[Target],
    refs: &HashMap<usize, MatchReferences>,
    runs: &[RunMs1<'_>],
    settings: &Settings,
) -> Result<(usize, Vec<Measurement>), String> {
    measure_with(targets, runs, settings, "transferred", |ms1, t, params| {
        let reference = t.donor.and_then(|d| refs.get(&d));
        let mut anchor = t.anchor_rt;
        let mut window = t.window_s;
        let mut decoy_anchor = anchor;
        let mut decoy_window = window;
        let mut match_score = None;
        let mut decoy_match = None;
        if let Some(r) = reference {
            let predicted = r.target.window_start(t.anchor_rt);
            if let Some(m) = matching::best_match(ms1, &r.target, predicted, t.search_s, params.tol)
            {
                anchor = 0.5 * (m.start_rt + m.stop_rt);
                window = m.window_s;
                match_score = Some(m.score);
            }
            if let Some(d) = &r.decoy {
                let predicted = d.window_start(t.anchor_rt);
                if let Some(m) = matching::best_match(ms1, d, predicted, t.search_s, params.tol) {
                    decoy_anchor = 0.5 * (m.start_rt + m.stop_rt);
                    decoy_window = m.window_s;
                    decoy_match = Some(m.score);
                }
            }
        }
        if match_score.is_none() || decoy_match.is_none() {
            match_score = None;
            decoy_match = None;
            anchor = t.anchor_rt;
            window = t.window_s;
        }
        let lfq = LfqTarget {
            anchor_rt: anchor,
            ..t.lfq_target()
        };
        let p = LfqParams {
            n_isotopes: t.envelope.len(),
            rt_window_s: window,
            ..params.clone()
        };
        let mut feature = quantify(ms1, &lfq, &p);
        if let Some(f) = feature.as_mut() {
            f.area = mono_area(f, &t.envelope);
        }
        let decoy = if match_score.is_some() {
            let twin = LfqTarget {
                mono_mz: t.mono_mz + params.decoy_mz_shift,
                anchor_rt: decoy_anchor,
                ..t.lfq_target()
            };
            quantify(
                ms1,
                &twin,
                &LfqParams {
                    n_isotopes: t.envelope.len(),
                    rt_window_s: decoy_window,
                    ..params.clone()
                },
            )
        } else {
            quantify_decoy(ms1, &lfq, &p)
        };
        Measurement {
            feature,
            decoy,
            match_score,
            decoy_match,
        }
    })
}

fn mono_area(f: &FeatureQuant, envelope: &[f64]) -> f64 {
    match (f.isotope_areas.first(), envelope.first()) {
        (Some(&mono), Some(&share)) if share > 0.0 => mono / share,
        _ => 0.0,
    }
}
