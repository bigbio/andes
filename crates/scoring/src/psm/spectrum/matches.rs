//! Walks over the theoretical ions of a peptide that find, or miss, a peak: node and
//! neutral-loss ions, their rank and mass error, and the nearest-peak lookups behind them.

use model::peptide::Peptide;

use super::SegmentPartitionSlice;
use crate::params::{IonType, Partition};
use crate::psm::rank::RankScorer;

#[allow(
    clippy::too_many_arguments,
    reason = "mirrors the scoring loop's argument bundle"
)]
#[allow(
    clippy::too_many_arguments,
    reason = "private inner matcher for the node-scoring hot path; every argument is an orthogonal, already-resolved input and bundling them would add a struct build per call"
)]
pub(super) fn visit_directional_node_ion_matches<F>(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    segment_partition_cache: SegmentPartitionSlice<'_>,
    scorer: &RankScorer,
    nominal_mass: f64,
    exact_mass: Option<f64>,
    is_prefix: bool,
    charge: u8,
    parent_mass: f64,
    mut visit: F,
) where
    F: FnMut(Partition, IonType, Option<(u32, f64)>, &[f32], f64, f64),
{
    use crate::params::IonType;

    let param = scorer.param();
    let node_window = scorer.node_score_window();
    let mme = &node_window;
    let num_segs = param.num_segments as usize;
    let use_cache = !segment_partition_cache.is_empty();
    #[allow(clippy::needless_range_loop)]
    for seg in 0..num_segs {
        let (partition, ion_logs_slice): (Partition, &[(IonType, Vec<f32>)]) = if use_cache {
            (
                segment_partition_cache[seg].0,
                segment_partition_cache[seg].1.as_slice(),
            )
        } else {
            let p = param.partition_for(charge, parent_mass, seg);
            let logs = scorer.partition_ion_logs(&p);
            (p, logs)
        };
        if !ion_logs_slice.is_empty() {
            for (ion, logs) in ion_logs_slice {
                let theo_mz = match (is_prefix, *ion) {
                    (true, IonType::Prefix { .. }) => match exact_mass {
                        Some(m) => ion.mz_exact(m),
                        None => ion.mz(nominal_mass),
                    },
                    (false, IonType::Suffix { .. }) => match exact_mass {
                        Some(m) => ion.mz_exact(m),
                        None => ion.mz(nominal_mass),
                    },
                    _ => continue,
                };
                if param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let matched = nearest_peak_rank_and_mz_in(peaks, ranks, theo_mz, tol_da);
                visit(partition, *ion, matched, logs, theo_mz, tol_da);
            }
        } else {
            for &ion in param.ion_types_for_partition_slice(charge, parent_mass, seg) {
                let theo_mz = match (is_prefix, ion) {
                    (true, IonType::Prefix { .. }) => match exact_mass {
                        Some(m) => ion.mz_exact(m),
                        None => ion.mz(nominal_mass),
                    },
                    (false, IonType::Suffix { .. }) => match exact_mass {
                        Some(m) => ion.mz_exact(m),
                        None => ion.mz(nominal_mass),
                    },
                    _ => continue,
                };
                if param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let matched = nearest_peak_rank_and_mz_in(peaks, ranks, theo_mz, tol_da);
                visit(partition, ion, matched, &[], theo_mz, tol_da);
            }
        }
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "mirrors visit_directional_node_ion_matches plus the active-loss slice"
)]
pub(super) fn visit_directional_loss_ion_matches<F>(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    segment_partition_cache: SegmentPartitionSlice<'_>,
    scorer: &RankScorer,
    nominal_mass: f64,
    is_prefix: bool,
    charge: u8,
    parent_mass: f64,
    active_losses: &[(f64, u8)],
    mut visit: F,
) where
    F: FnMut(Partition, IonType, Option<u32>, &[f32], f64, f64),
{
    use crate::params::IonType;
    if active_losses.is_empty() {
        return;
    }
    let param = scorer.param();
    let node_window = scorer.node_score_window();
    let mme = &node_window;
    let num_segs = param.num_segments as usize;
    let use_cache = !segment_partition_cache.is_empty();
    #[allow(clippy::needless_range_loop)]
    for seg in 0..num_segs {
        let partition = if use_cache {
            segment_partition_cache[seg].0
        } else {
            param.partition_for(charge, parent_mass, seg)
        };
        let loss_logs = scorer.partition_loss_ion_logs(&partition);
        if loss_logs.is_empty() {
            continue;
        }
        for (ion, logs) in loss_logs {
            match (is_prefix, ion) {
                (true, IonType::Prefix { .. }) | (false, IonType::Suffix { .. }) => {}
                _ => continue,
            }
            let cls = ion.loss_class();
            let ion_charge = ion.charge().unwrap_or(1).max(1) as f64;
            let base_mz = ion.mz(nominal_mass);
            for &(loss, lcls) in active_losses {
                if lcls != cls {
                    continue;
                }
                let theo_mz = base_mz - loss / ion_charge;
                if theo_mz <= 0.0 || param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let rank = nearest_peak_rank_in(peaks, ranks, theo_mz, tol_da);
                visit(partition, *ion, rank, logs, theo_mz, tol_da);
            }
        }
    }
}

pub(crate) fn span_losses(peptide: &Peptide, range: std::ops::Range<usize>) -> Vec<(f64, u8)> {
    let mut v: Vec<(f64, u8)> = Vec::new();
    for i in range {
        if let Some(m) = peptide.residues[i].mod_.as_ref() {
            if m.loss_class != 0 {
                for &l in &m.neutral_losses {
                    if !v
                        .iter()
                        .any(|&(el, ec)| ec == m.loss_class && (el - l).abs() < 1e-9)
                    {
                        v.push((l, m.loss_class));
                    }
                }
            }
        }
    }
    v
}

pub(super) fn ion_fact_rank_and_error(
    rank: Option<u32>,
    peaks: &[(f64, f32)],
    ranks: &[u32],
    theo_mz: f64,
    tol_da: f64,
    max_rank: u32,
    esf: i32,
) -> (Option<u32>, Option<u32>) {
    let matched = rank.map(|r| {
        let mz = matched_peak_mz(peaks, ranks, theo_mz, tol_da).unwrap_or(theo_mz);
        (r, mz)
    });
    ion_fact_rank_and_error_from_match(matched, theo_mz, max_rank, esf)
}

pub(super) fn ion_fact_rank_and_error_from_match(
    matched: Option<(u32, f64)>,
    theo_mz: f64,
    max_rank: u32,
    esf: i32,
) -> (Option<u32>, Option<u32>) {
    match matched {
        Some((r, peak_mz)) => {
            let clamped = r.min(max_rank).max(1);
            let ebin = if esf > 0 {
                let error_da = (peak_mz - theo_mz) as f32;
                let mut idx = (error_da * esf as f32).round() as i32;
                if idx > esf {
                    idx = esf;
                } else if idx < -esf {
                    idx = -esf;
                }
                idx += esf;
                Some(idx as u32)
            } else {
                None
            };
            (Some(clamped), ebin)
        }
        None => (None, None),
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "mirrors visit_directional_node_ion_matches plus the active-loss slice"
)]
pub(super) fn visit_directional_loss_ion_facts<F>(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    segment_partition_cache: SegmentPartitionSlice<'_>,
    scorer: &RankScorer,
    nominal_mass: f64,
    exact_mass: Option<f64>,
    is_prefix: bool,
    charge: u8,
    parent_mass: f64,
    active_losses: &[(f64, u8)],
    mut visit: F,
) where
    F: FnMut(Partition, IonType, Option<u32>, f64, f64),
{
    use crate::params::IonType;
    if active_losses.is_empty() {
        return;
    }
    let param = scorer.param();
    let node_window = scorer.node_score_window();
    let mme = &node_window;
    let num_segs = param.num_segments as usize;
    let use_cache = !segment_partition_cache.is_empty();
    #[allow(clippy::needless_range_loop)]
    for seg in 0..num_segs {
        let (partition, intact_logs): (Partition, &[(IonType, Vec<f32>)]) = if use_cache {
            (
                segment_partition_cache[seg].0,
                segment_partition_cache[seg].1.as_slice(),
            )
        } else {
            let p = param.partition_for(charge, parent_mass, seg);
            (p, scorer.partition_ion_logs(&p))
        };
        for (ion, _logs) in intact_logs {
            let (icharge, ioff) = match (is_prefix, ion) {
                (
                    true,
                    IonType::Prefix {
                        charge,
                        offset_bits,
                        ..
                    },
                ) => (*charge, *offset_bits),
                (
                    false,
                    IonType::Suffix {
                        charge,
                        offset_bits,
                        ..
                    },
                ) => (*charge, *offset_bits),
                _ => continue,
            };
            let ion_charge = icharge.max(1) as f64;
            let base_mz = match exact_mass {
                Some(m) => ion.mz_exact(m),
                None => ion.mz(nominal_mass),
            };
            for &(loss, cls) in active_losses {
                let theo_mz = base_mz - loss / ion_charge;
                if theo_mz <= 0.0 || param.segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = mme.as_da(theo_mz);
                let rank = nearest_peak_rank_in(peaks, ranks, theo_mz, tol_da);
                let loss_ion = if is_prefix {
                    IonType::Prefix {
                        charge: icharge,
                        offset_bits: ioff,
                        loss_class: cls,
                    }
                } else {
                    IonType::Suffix {
                        charge: icharge,
                        offset_bits: ioff,
                        loss_class: cls,
                    }
                };
                visit(partition, loss_ion, rank, theo_mz, tol_da);
            }
        }
    }
}

pub(super) fn matched_peak_mz(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    target_mz: f64,
    tolerance_da: f64,
) -> Option<f64> {
    if peaks.is_empty() {
        return None;
    }
    let lo_mz = target_mz - tolerance_da;
    let hi_mz = target_mz + tolerance_da;
    let start = peaks.partition_point(|&(mz, _)| mz < lo_mz);
    let mut best: Option<(f64, f32)> = None;
    for i in start..peaks.len() {
        let (mz, intensity) = peaks[i];
        if mz > hi_mz {
            break;
        }
        if ranks[i] == u32::MAX {
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_int)| intensity > *best_int)
        {
            best = Some((mz, intensity));
        }
    }
    best.map(|(mz, _)| mz)
}

pub(super) fn nearest_peak_rank_in(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    target_mz: f64,
    tolerance_da: f64,
) -> Option<u32> {
    nearest_peak_rank_and_mz_in(peaks, ranks, target_mz, tolerance_da).map(|(rank, _mz)| rank)
}

pub(super) fn nearest_peak_rank_and_mz_in(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    target_mz: f64,
    tolerance_da: f64,
) -> Option<(u32, f64)> {
    if peaks.is_empty() {
        return None;
    }
    let lo_mz = target_mz - tolerance_da;
    let hi_mz = target_mz + tolerance_da;
    let start = peaks.partition_point(|&(mz, _)| mz < lo_mz);
    let mut best: Option<(usize, f32)> = None;
    for i in start..peaks.len() {
        let (mz, intensity) = peaks[i];
        if mz > hi_mz {
            break;
        }
        if ranks[i] == u32::MAX {
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_int)| intensity > *best_int)
        {
            best = Some((i, intensity));
        }
    }
    best.map(|(i, _)| (ranks[i], peaks[i].0))
}

pub(super) fn max_active_intensity_of(peaks: &[(f64, f32)], ranks: &[u32]) -> f64 {
    peaks
        .iter()
        .zip(ranks.iter())
        .filter(|(_, &rank)| rank != u32::MAX)
        .map(|(&(_, intensity), _)| intensity as f64)
        .fold(0.0_f64, f64::max)
}
