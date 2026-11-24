//! Charge deconvolution of a spectrum before scoring, and the main ion type of a scoring
//! parameter set.

use super::PROTON;
use crate::params::{IonType, Param};

pub(super) type DeconvResult = (Option<Vec<(f64, f32)>>, Option<Vec<u32>>);

pub(super) fn deconvolute_spectrum(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    precursor_charge: u8,
    tol: f64,
) -> (Vec<(f64, f32)>, Vec<u32>) {
    const ISOTOPE: f64 = model::mass::ISOTOPE;

    let n = peaks.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut ignore = vec![false; n];
    let mut out: Vec<(f64, f32, u32)> = Vec::with_capacity(n);
    let charge_i32 = precursor_charge as i32;

    for i in 0..n {
        if ignore[i] {
            continue;
        }
        let (mut p_mz, p_int) = peaks[i];
        let p_rank = ranks[i];

        for ion_charge_i in 2..charge_i32.min(4) {
            let ion_charge = ion_charge_i as f64;
            let expected_diff = ISOTOPE / ion_charge;
            let mut is_deconvoluted = false;
            for j in (i + 1)..n {
                let (p2_mz, p2_int) = peaks[j];
                let diff = p2_mz - p_mz - expected_diff;
                if diff > -tol && diff < tol {
                    ignore[j] = true;
                    let p_new_mz = ion_charge * p_mz - (ion_charge - 1.0) * PROTON;
                    let p2_new_mz = ion_charge * p2_mz - (ion_charge - 1.0) * PROTON;
                    p_mz = p_new_mz;
                    is_deconvoluted = true;

                    let neutral_mass = ion_charge * (p2_mz - PROTON);
                    let p3_diff_expected =
                        model::isotope::averagine_m1_m2_spacing(neutral_mass) / ion_charge;
                    for k in (j + 1)..n {
                        let (p3_mz, p3_int) = peaks[k];
                        let diff2 = p3_mz - p2_mz - p3_diff_expected;
                        if diff2 > -tol && diff2 < tol {
                            ignore[k] = true;
                            let p3_new_mz = ion_charge * p3_mz - (ion_charge - 1.0) * PROTON;
                            out.push((p3_new_mz, p3_int, ranks[k]));
                            break;
                        } else if diff2 > tol {
                            break;
                        }
                    }
                    out.push((p2_new_mz, p2_int, ranks[j]));
                    break;
                } else if diff > tol {
                    break;
                }
            }
            if is_deconvoluted {
                break;
            }
        }
        out.push((p_mz, p_int, p_rank));
    }

    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut out_peaks: Vec<(f64, f32)> = Vec::with_capacity(out.len());
    let mut out_ranks: Vec<u32> = Vec::with_capacity(out.len());
    for (mz, intensity, rank) in out {
        out_peaks.push((mz, intensity));
        out_ranks.push(rank);
    }
    (out_peaks, out_ranks)
}

pub(super) fn main_ion_from_param(param: &Param, partition: crate::params::Partition) -> IonType {
    let fallback = IonType::Prefix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 0,
    };
    let num_segments = param.num_segments.max(1) as usize;
    let mut ion_freq: std::collections::HashMap<IonType, f32> = std::collections::HashMap::new();
    for seg in 0..num_segments {
        let part = crate::params::Partition {
            charge: partition.charge,
            parent_mass: partition.parent_mass,
            seg_num: seg as i32,
        };
        if let Some(frag_list) = param.frag_off_table.get(&part) {
            for f in frag_list {
                if matches!(f.ion_type, IonType::Noise) {
                    continue;
                }
                *ion_freq.entry(f.ion_type).or_insert(0.0) += f.frequency;
            }
        }
    }
    let tiebreak = |ion: IonType| -> (u8, i32, u32, u8) {
        match ion {
            IonType::Suffix {
                charge,
                offset_bits,
                loss_class,
            } => (0, charge, offset_bits, loss_class),
            IonType::Prefix {
                charge,
                offset_bits,
                loss_class,
            } => (1, charge, offset_bits, loss_class),
            IonType::Noise => (2, 0, 0, 0),
        }
    };
    let mut best_ion: Option<IonType> = None;
    let mut best_freq = f32::NEG_INFINITY;
    let mut best_key = (u8::MAX, i32::MAX, u32::MAX, u8::MAX);
    for (&ion, &freq) in &ion_freq {
        let key = tiebreak(ion);
        if freq > best_freq || (freq == best_freq && key < best_key) {
            best_freq = freq;
            best_key = key;
            best_ion = Some(ion);
        }
    }
    best_ion.unwrap_or(fallback)
}
