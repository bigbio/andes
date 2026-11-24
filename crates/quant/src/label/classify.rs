//! Channel classification: decoy shifts, channel scores and q-values, identified or
//! transferred status, log ratios and retention-time shifts between channels.

use model::fdr::{ranked_q, Counts, Eval};
use model::mass::ISOTOPE;
use model::stats::{median_mut, nan_low, Even};

use super::fit::ChannelFit;

const DECOY_UNITS: [f64; 12] = [
    7.0, -7.0, 11.0, -11.0, 13.0, -13.0, 17.0, -17.0, 19.0, -19.0, 23.0, -23.0,
];

pub fn decoy_shifts(
    moved: &[usize],
    deltas: &[f64],
    min_sep: f64,
    seed: u64,
    n: usize,
) -> Vec<f64> {
    let ok: Vec<f64> = DECOY_UNITS
        .iter()
        .map(|u| u * crate::decoys::PEPTIDE_MASS_SLOPE)
        .filter(|&s| {
            moved.iter().all(|&m| {
                deltas.iter().all(|&d| {
                    (-3..=6).all(|k| (deltas[m] + s - d - k as f64 * ISOTOPE).abs() >= min_sep)
                })
            })
        })
        .collect();
    let k = ok.len().max(1);
    (0..n.min(ok.len()))
        .map(|i| ok[(seed as usize).wrapping_add(i) % k])
        .collect()
}

pub fn decoy_min_sep(mass: f64, tol_ppm: f64) -> f64 {
    3.0 * tol_ppm * 1e-6 * mass
}

pub fn mass_factor(err_ppm: Option<f64>, expected_ppm: f64, sigma_ppm: f64) -> f64 {
    let Some(err_ppm) = err_ppm else {
        return 0.0;
    };
    let d = (err_ppm - expected_ppm) / sigma_ppm.max(1e-3);
    (-0.5 * d * d).exp()
}

pub fn coelution_factor(offset: f64, sigma_s: f64) -> f64 {
    let d = offset / sigma_s.max(1e-3);
    (-0.5 * d * d).exp()
}

pub fn channel_score(fit: &ChannelFit) -> f64 {
    (2.0 * fit.llr.max(0.0)).sqrt() * fit.envelope_cos.max(0.0) * fit.coherence.max(0.0)
}

pub fn channel_qvalues(rows: &[(f64, bool)], decoys_per_target: usize) -> Vec<f64> {
    let k = decoys_per_target.max(1) as f64;
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        nan_low(rows[b].0)
            .total_cmp(&nan_low(rows[a].0))
            .then(rows[b].1.cmp(&rows[a].1))
    });
    let ranked = ranked_q::<Counts>(
        rows.len(),
        |c, rank| c.add(rows[order[rank]].1),
        |_, _| false,
        Eval::PerRow,
        |c| c.ratio_plus_one_over(k),
        1.0,
    );
    let mut q = vec![1.0; rows.len()];
    for (&i, qi) in order.iter().zip(ranked) {
        q[i] = qi;
    }
    q
}

pub fn weighted_qvalues(rows: &[(f64, f64)]) -> Vec<f64> {
    let decoys: Vec<f64> = rows.iter().filter(|r| r.1 > 0.0).map(|r| r.1).collect();
    let w_bar = if decoys.is_empty() {
        1.0
    } else {
        decoys.iter().sum::<f64>() / decoys.len() as f64
    };
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        nan_low(rows[b].0)
            .total_cmp(&nan_low(rows[a].0))
            .then((rows[b].1 > 0.0).cmp(&(rows[a].1 > 0.0)))
    });
    let ranked = ranked_q::<(usize, f64)>(
        rows.len(),
        |(t, d), rank| {
            let w = rows[order[rank]].1;
            if w > 0.0 {
                *d += w;
            } else {
                *t += 1;
            }
        },
        |_, _| false,
        Eval::PerRow,
        |&(t, d)| (d + w_bar) / t.max(1) as f64,
        1.0,
    );
    let mut q = vec![1.0; rows.len()];
    for (&i, qi) in order.iter().zip(ranked) {
        q[i] = qi;
    }
    q
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChannelStatus {
    Identified,
    Matched,
    Transferred,
    BelowDetection,
    NotQuantifiable,
}

impl ChannelStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            ChannelStatus::Identified => "identified",
            ChannelStatus::Matched => "match_channel",
            ChannelStatus::Transferred => "transfer",
            ChannelStatus::BelowDetection => "below_detection",
            ChannelStatus::NotQuantifiable => "not_quantifiable",
        }
    }

    pub fn is_quantified(self) -> bool {
        matches!(
            self,
            ChannelStatus::Identified | ChannelStatus::Matched | ChannelStatus::Transferred
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ClassifyParams {
    pub fdr: f64,
    pub min_signal_llr: f64,
    pub min_cosine: f64,
    pub weak_fraction: f64,
    pub min_weak_coherence: f64,
    pub max_m1_ratio: f64,
}

impl Default for ClassifyParams {
    fn default() -> Self {
        Self {
            fdr: 0.01,
            min_signal_llr: 4.0,
            min_cosine: 0.7,
            weak_fraction: 0.2,
            min_weak_coherence: 0.5,
            max_m1_ratio: 0.3,
        }
    }
}

pub fn classify_channel(
    fit: &ChannelFit,
    q: f64,
    anchor_amplitude: f64,
    claimed: bool,
    p: &ClassifyParams,
) -> ChannelStatus {
    if claimed {
        return ChannelStatus::NotQuantifiable;
    }
    if fit.amplitude <= 0.0 || fit.llr < p.min_signal_llr {
        return ChannelStatus::BelowDetection;
    }
    let weak = fit.amplitude < p.weak_fraction * anchor_amplitude;
    if q > p.fdr
        || fit.envelope_cos < p.min_cosine
        || fit.m1_ratio > p.max_m1_ratio
        || (weak && fit.coherence < p.min_weak_coherence)
    {
        return ChannelStatus::NotQuantifiable;
    }
    ChannelStatus::Matched
}

pub fn classify_identified(fit: &ChannelFit, p: &ClassifyParams) -> ChannelStatus {
    if fit.amplitude > 0.0 && (fit.envelope_cos < p.min_cosine || fit.m1_ratio > p.max_m1_ratio) {
        ChannelStatus::NotQuantifiable
    } else {
        ChannelStatus::Identified
    }
}

pub fn classify_transferred(fit: &ChannelFit, strongest: f64, p: &ClassifyParams) -> ChannelStatus {
    match classify_channel(fit, 0.0, strongest, false, p) {
        ChannelStatus::Matched => ChannelStatus::Transferred,
        other => other,
    }
}

pub fn log2_ratio(c: &ChannelFit, reference: &ChannelFit) -> Option<(f64, f64)> {
    if c.amplitude <= 0.0 || reference.amplitude <= 0.0 {
        return None;
    }
    let r = (c.amplitude / reference.amplitude).log2();
    let rel = (c.amplitude_se / c.amplitude).powi(2)
        + (reference.amplitude_se / reference.amplitude).powi(2);
    Some((r, rel.sqrt() / std::f64::consts::LN_2))
}

pub fn learn_rt_shifts(pairs: &[(usize, f64)], n_channels: usize, min_pairs: usize) -> Vec<f64> {
    (0..n_channels)
        .map(|c| {
            let mut v: Vec<f64> = pairs
                .iter()
                .filter(|p| p.0 == c && p.1.is_finite())
                .map(|p| p.1)
                .collect();
            if v.len() < min_pairs.max(1) {
                return 0.0;
            }
            median_mut(&mut v, Even::Mean).unwrap_or(0.0)
        })
        .collect()
}

pub fn median(mut v: Vec<f64>) -> f64 {
    v.retain(|x| x.is_finite());
    median_mut(&mut v, Even::Mean).unwrap_or(f64::NAN)
}
