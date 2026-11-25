//! Match between runs for glycopeptides: piecewise RT maps between runs, alignment, the
//! transfers they allow and the optional match-training dump.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use model::stats::{median_mad, Even, MAD_TO_SD};

use super::{accepted, Evidence, FeatureRow, Settings, Target};
use crate::glyco::matching;
use crate::quant::run::RunMs1;

pub(super) const MBR_MIN_ANCHORS: usize = 50;
const MBR_MAX_RESIDUAL_S: f64 = 45.0;
const MBR_ANCHORS_PER_KNOT: usize = 300;
const MBR_MIN_KNOTS: usize = 3;
const MBR_MAX_KNOTS: usize = 10;
pub(super) const MBR_MIN_WINDOW_S: f64 = 10.0;
const MBR_MAX_WINDOW_S: f64 = 60.0;
const MBR_MATCH_MIN_SEARCH_S: f64 = 60.0;
const MBR_MATCH_MAX_SEARCH_S: f64 = 120.0;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct TransferCounts {
    pub(super) maps: usize,
    pub(super) anchors: usize,
    pub(super) candidates: usize,
    pub(super) measured: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct RtMap {
    xs: Vec<f64>,
    ys: Vec<f64>,
    pub(super) sigma_s: f64,
}

impl RtMap {
    fn identity() -> RtMap {
        RtMap {
            xs: vec![0.0, 1.0],
            ys: vec![0.0, 1.0],
            sigma_s: 0.0,
        }
    }

    pub(super) fn forward(&self, rt: f64) -> f64 {
        interpolate(&self.xs, &self.ys, rt)
    }

    pub(super) fn inverse(&self, rt_ref: f64) -> Option<f64> {
        let ys_are_increasing = self.ys.windows(2).all(|w| w[1] > w[0]);
        ys_are_increasing.then(|| interpolate(&self.ys, &self.xs, rt_ref))
    }

    fn inverse_scale(&self, rt_ref: f64) -> f64 {
        let i = self
            .ys
            .partition_point(|&k| k <= rt_ref)
            .saturating_sub(1)
            .min(self.ys.len() - 2);
        (self.xs[i + 1] - self.xs[i]) / (self.ys[i + 1] - self.ys[i])
    }
}

fn interpolate(a: &[f64], b: &[f64], v: f64) -> f64 {
    debug_assert!(a.len() >= 2 && a.len() == b.len());
    let i = a
        .partition_point(|&k| k <= v)
        .saturating_sub(1)
        .min(a.len() - 2);
    let slope = (b[i + 1] - b[i]) / (a[i + 1] - a[i]);
    b[i] + slope * (v - a[i])
}

fn fit_line(pairs: &[(f64, f64)]) -> Option<(f64, f64)> {
    if pairs.len() < 3 {
        return None;
    }
    let n = pairs.len() as f64;
    let mx = pairs.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pairs.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = pairs.iter().map(|p| (p.0 - mx).powi(2)).sum();
    if sxx.is_nan() || sxx <= 0.0 {
        return None;
    }
    let sxy: f64 = pairs.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let slope = sxy / sxx;
    if !slope.is_finite() || slope <= 0.0 {
        return None;
    }
    Some((my - slope * mx, slope))
}

pub(super) fn fit_rt_map(pairs: &[(f64, f64)]) -> Option<RtMap> {
    let (intercept, slope) = fit_line(pairs)?;
    let residuals: Vec<f64> = pairs
        .iter()
        .map(|&(x, y)| y - (intercept + slope * x))
        .collect();
    let (center, mad) = median_mad(&residuals, Even::Mean).unwrap_or((f64::NAN, f64::NAN));
    let mad = MAD_TO_SD * mad;
    let mut kept: Vec<(f64, f64)> = pairs
        .iter()
        .zip(&residuals)
        .filter(|(_, r)| mad <= 0.0 || ((*r - center).abs()) <= 3.0 * mad)
        .map(|(p, _)| *p)
        .collect();
    if kept.len() < 3 {
        return None;
    }
    kept.sort_by(|a, b| a.0.total_cmp(&b.0));

    let bins = (kept.len() / MBR_ANCHORS_PER_KNOT).clamp(MBR_MIN_KNOTS, MBR_MAX_KNOTS);
    let (mut xs, mut ys) = (Vec::with_capacity(bins), Vec::with_capacity(bins));
    for k in 0..bins {
        let bin = &kept[k * kept.len() / bins..(k + 1) * kept.len() / bins];
        let x = crate::quant::mbr::median(bin.iter().map(|p| p.0).collect());
        let y = crate::quant::mbr::median(bin.iter().map(|p| p.1).collect());
        let advances = xs.last().is_none_or(|&px| x > px) && ys.last().is_none_or(|&py| y > py);
        if advances {
            xs.push(x);
            ys.push(y);
        }
    }
    if xs.len() < 2 {
        return None;
    }
    let map = RtMap {
        xs,
        ys,
        sigma_s: 0.0,
    };
    let dev: Vec<f64> = kept.iter().map(|&(x, y)| y - map.forward(x)).collect();
    let sigma = MAD_TO_SD * median_mad(&dev, Even::Mean).map_or(f64::NAN, |(_, mad)| mad);
    let map = RtMap {
        sigma_s: sigma,
        ..map
    };
    map.sigma_s.is_finite().then_some(map)
}

pub(super) type GlycoKey = (String, String, u8);

fn glyco_key(t: &Target) -> GlycoKey {
    (t.peptide.clone(), t.glycan.clone(), t.charge)
}

pub(super) struct Alignment {
    pub(super) identified: Vec<std::collections::HashSet<GlycoKey>>,
    pub(super) donor: Vec<HashMap<GlycoKey, usize>>,
    pub(super) maps: Vec<Option<RtMap>>,
    counts: TransferCounts,
}

impl Alignment {
    pub(super) fn build(rows: &[FeatureRow], run_names: &[&str], settings: &Settings) -> Self {
        let mut counts = TransferCounts {
            maps: 1,
            ..TransferCounts::default()
        };
        let mut identified: Vec<std::collections::HashSet<GlycoKey>> =
            vec![std::collections::HashSet::new(); run_names.len()];
        let mut donor: Vec<HashMap<GlycoKey, usize>> = vec![HashMap::new(); run_names.len()];
        for (i, row) in rows.iter().enumerate() {
            if row.target.evidence != Evidence::Identified {
                continue;
            }
            let key = glyco_key(&row.target);
            identified[row.target.run_idx].insert(key.clone());
            if !accepted(row, settings) {
                continue;
            }
            let score = row.score.unwrap_or(f32::NEG_INFINITY);
            match donor[row.target.run_idx].get(&key) {
                Some(&j) if rows[j].score.unwrap_or(f32::NEG_INFINITY) >= score => {}
                _ => {
                    donor[row.target.run_idx].insert(key, i);
                }
            }
        }

        let mut maps: Vec<Option<RtMap>> = vec![None; run_names.len()];
        maps[0] = Some(RtMap::identity());
        for (r, map) in maps.iter_mut().enumerate().skip(1) {
            let pairs: Vec<(f64, f64)> = donor[r]
                .iter()
                .filter_map(|(key, &j)| {
                    let other = *donor[0].get(key)?;
                    let rt = rows[j].feature.as_ref()?.apex_rt;
                    let rt_ref = rows[other].feature.as_ref()?.apex_rt;
                    Some((rt, rt_ref))
                })
                .collect();
            if pairs.len() < MBR_MIN_ANCHORS {
                log::warn!(
                    target: "glyco-mbr",
                    "{}: {} shared glycopeptides, fewer than {MBR_MIN_ANCHORS} \
                     needed for an RT map; no transfers into or out of it",
                    run_names[r],
                    pairs.len()
                );
                continue;
            }
            let Some(fitted) = fit_rt_map(&pairs) else {
                continue;
            };
            if fitted.sigma_s > MBR_MAX_RESIDUAL_S {
                log::warn!(
                    target: "glyco-mbr",
                    "{}: RT map residual {:.1}s is above {MBR_MAX_RESIDUAL_S}s; \
                     no transfers into or out of it",
                    run_names[r], fitted.sigma_s
                );
                continue;
            }
            counts.maps += 1;
            counts.anchors += pairs.len();
            *map = Some(fitted);
        }
        Self {
            identified,
            donor,
            maps,
            counts,
        }
    }

    fn predict(&self, from_idx: usize, to_idx: usize, rt: f64) -> Option<(f64, f64)> {
        let from = self.maps[from_idx].as_ref()?;
        let to = self.maps[to_idx].as_ref()?;
        let t_ref = from.forward(rt);
        let predicted = to.inverse(t_ref)?;
        let sigma = (from.sigma_s.powi(2) + (to.sigma_s * to.inverse_scale(t_ref)).powi(2)).sqrt();
        Some((predicted, sigma))
    }
}

pub(super) fn build_transfers(
    rows: &[FeatureRow],
    run_names: &[&str],
    settings: &Settings,
    alignment: &Alignment,
) -> (Vec<Target>, TransferCounts) {
    let Alignment {
        identified,
        donor,
        maps,
        counts,
    } = alignment;
    let mut counts = *counts;
    let mut best: HashMap<(GlycoKey, usize), (f32, Target)> = HashMap::new();
    for (i, row) in rows.iter().enumerate() {
        if row.target.evidence != Evidence::Identified
            || !accepted(row, settings)
            || donor[row.target.run_idx].get(&glyco_key(&row.target)) != Some(&i)
        {
            continue;
        }
        let Some(feature) = &row.feature else {
            continue;
        };
        if maps[row.target.run_idx].is_none() {
            continue;
        }
        let key = glyco_key(&row.target);
        for (r, to) in maps.iter().enumerate() {
            if r == row.target.run_idx {
                continue;
            }
            if to.is_none() {
                continue;
            }
            if identified[r].contains(&key) {
                continue;
            }
            let Some((rt, sigma)) = alignment.predict(row.target.run_idx, r, feature.apex_rt)
            else {
                continue;
            };
            let window = (3.0 * sigma).clamp(MBR_MIN_WINDOW_S, MBR_MAX_WINDOW_S);
            let search = (3.0 * sigma).clamp(MBR_MATCH_MIN_SEARCH_S, MBR_MATCH_MAX_SEARCH_S);
            let score = row.score.unwrap_or(f32::NEG_INFINITY);
            let slot = (key.clone(), r);
            if best
                .get(&slot)
                .is_some_and(|(existing, _)| *existing >= score)
            {
                continue;
            }
            let mut target = row.target.clone();
            target.run_idx = r;
            target.run = run_names[r].to_string();
            target.spec_id = format!("{}_mbr_{}", target.spec_id, run_names[r]);
            target.scan = 0;
            target.evidence = Evidence::Transferred;
            target.anchor_rt = rt;
            target.window_s = window;
            target.search_s = search;
            target.donor = Some(i);
            best.insert(slot, (score, target));
        }
    }
    let mut transfers: Vec<Target> = best.into_values().map(|(_, t)| t).collect();
    transfers.sort_by(|a, b| {
        (
            a.run_idx,
            &a.peptide,
            &a.glycan,
            &a.topology,
            a.charge,
            &a.spec_id,
        )
            .cmp(&(
                b.run_idx,
                &b.peptide,
                &b.glycan,
                &b.topology,
                b.charge,
                &b.spec_id,
            ))
    });
    counts.candidates = transfers.len();
    (transfers, counts)
}

pub(super) fn dump_match_training(
    rows: &[FeatureRow],
    runs: &[RunMs1<'_>],
    settings: &Settings,
    alignment: &Alignment,
    path: &Path,
) -> Result<(), String> {
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(path).map_err(|e| format!("{}: {e}", path.display()))?,
    );
    writeln!(
        out,
        "# glyco-mbr match training table: label, match-model features \
         f0..f9, peptide, donor/target run"
    )
    .map_err(|e| format!("{}: {e}", path.display()))?;
    let mut pairs = 0usize;
    let mut positives = 0usize;
    let mut negatives = 0usize;
    for d in 0..runs.len() {
        if alignment.maps[d].is_none() || alignment.donor[d].is_empty() {
            continue;
        }
        let ms1_d = runs[d]
            .load()
            .map_err(|e| format!("--glyco-mbr-match-train: {e}"))?;
        let mut refs: Vec<(
            GlycoKey,
            usize,
            matching::Reference,
            Option<matching::Reference>,
        )> = Vec::with_capacity(alignment.donor[d].len());
        for (key, &j) in &alignment.donor[d] {
            let row = &rows[j];
            let Some(f) = &row.feature else { continue };
            let t = &row.target;
            let build = |mono_mz: f64| {
                matching::Reference::build(
                    &ms1_d,
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
            refs.push((key.clone(), j, target, decoy));
        }
        drop(ms1_d);
        for (t_run, run) in runs.iter().enumerate() {
            if t_run == d || alignment.maps[t_run].is_none() {
                continue;
            }
            let mut truth: HashMap<&GlycoKey, usize> = HashMap::new();
            for (key, &j) in &alignment.donor[t_run] {
                truth.insert(key, j);
            }
            if truth.is_empty() {
                continue;
            }
            let ms1 = run
                .load()
                .map_err(|e| format!("--glyco-mbr-match-train: {e}"))?;
            let (mut pair_pos, mut pair_neg) = (0usize, 0usize);
            for (key, donor_row, reference, decoy) in &refs {
                let Some(&truth_row) = truth.get(key) else {
                    continue;
                };
                let donor = &rows[*donor_row];
                let true_row = &rows[truth_row];
                let (Some(df), Some(tf)) = (&donor.feature, &true_row.feature) else {
                    continue;
                };
                let Some((predicted_rt, sigma)) = alignment.predict(d, t_run, df.apex_rt) else {
                    continue;
                };
                let search = (3.0 * sigma).clamp(MBR_MATCH_MIN_SEARCH_S, MBR_MATCH_MAX_SEARCH_S);
                let Some(candidates) = matching::score_candidates(
                    &ms1,
                    reference,
                    reference.window_start(predicted_rt),
                    search,
                    settings.params.tol,
                ) else {
                    continue;
                };
                let scored = matching::scored_top(&candidates);
                if scored.is_empty() {
                    continue;
                }
                let truth_width = (tf.rt_stop - tf.rt_start).max(0.0);
                let min_overlap = 0.5 * truth_width;
                let overlaps: Vec<f64> = scored
                    .iter()
                    .map(|m| (m.stop_rt.min(tf.rt_stop) - m.start_rt.max(tf.rt_start)).max(0.0))
                    .collect();
                if !overlaps.iter().any(|&o| o >= min_overlap && o > 0.0) {
                    continue;
                }
                let group = &key.0;
                for (m, &overlap) in scored.iter().zip(&overlaps) {
                    let label = u8::from(overlap >= min_overlap && overlap > 0.0);
                    write_training_row(&mut out, label, m, group, d, t_run)
                        .map_err(|e| format!("{}: {e}", path.display()))?;
                    if label == 1 {
                        pair_pos += 1;
                    } else {
                        pair_neg += 1;
                    }
                }
                if let Some(decoy) = decoy {
                    if let Some(candidates) = matching::score_candidates(
                        &ms1,
                        decoy,
                        decoy.window_start(predicted_rt),
                        search,
                        settings.params.tol,
                    ) {
                        for m in matching::scored_top(&candidates) {
                            write_training_row(&mut out, 0, &m, group, d, t_run)
                                .map_err(|e| format!("{}: {e}", path.display()))?;
                            pair_neg += 1;
                        }
                    }
                }
                pairs += 1;
            }
            positives += pair_pos;
            negatives += pair_neg;
            log::debug!(
                target: "glyco-mbr-match-train",
                "runs {d}->{t_run}: {pair_pos} positive, \
                 {pair_neg} negative candidate windows"
            );
        }
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", path.display()))?;
    log::info!(
        target: "glyco-mbr-match-train",
        "{} glycopeptide pairs, {positives} positive and \
         {negatives} negative rows written to {}",
        pairs,
        path.display()
    );
    Ok(())
}

fn write_training_row(
    out: &mut impl Write,
    label: u8,
    m: &matching::Match,
    group: &str,
    donor_run: usize,
    target_run: usize,
) -> std::io::Result<()> {
    write!(out, "{label}")?;
    for f in m.features {
        write!(out, "\t{f:.6}")?;
    }
    writeln!(out, "\t{group}\t{donor_run}>{target_run}")
}
