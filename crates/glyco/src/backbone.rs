use crate::glycan_first::GlycanCore;
use crate::glycan_mass::{CORE_Y_STEPS, O_CORE_Y_STEPS, PROTON};

#[inline]
fn core_y_steps(core: GlycanCore) -> &'static [f64] {
    match core {
        GlycanCore::NGlycan => &CORE_Y_STEPS,
        GlycanCore::OGlycan => &O_CORE_Y_STEPS,
    }
}

pub(crate) const H2O: f64 = 18.010565;

const BY_COMPLEMENT_OFFSET: f64 = 2.0 * PROTON;

#[derive(Debug, Clone)]
pub struct BackboneCandidate {
    pub backbone_mass: f64,
    pub core_y_hits: u8,
    pub votes: u32,
    pub intensity_score: f64,
    pub complement_score: f64,
}

fn complement_score(peaks: &[(f64, f32)], norm_sqrt: &[f64], bb: f64, tol_ppm: f64) -> f64 {
    let target = bb + BY_COMPLEMENT_OFFSET;
    let tol = (target * (tol_ppm / 1e6)).max(0.02);

    let n = peaks.len();
    if n < 2 {
        return 0.0;
    }

    let b_lo = 50.0_f64;
    let b_hi = bb / 2.0;

    let mut score = 0.0_f64;

    let mut lo_cur = n;
    let mut hi_cur = n;
    for (i, &(b_mz, _)) in peaks.iter().enumerate() {
        if b_mz < b_lo {
            continue;
        }
        if b_mz > b_hi {
            break;
        }
        let y_mz = target - b_mz;
        if y_mz <= b_mz || y_mz <= 0.0 {
            continue;
        }
        let y_lo = y_mz - tol;
        let y_hi = y_mz + tol;
        while hi_cur > 0 && peaks[hi_cur - 1].0 > y_hi {
            hi_cur -= 1;
        }
        if lo_cur > hi_cur {
            lo_cur = hi_cur;
        }
        while lo_cur > 0 && peaks[lo_cur - 1].0 >= y_lo {
            lo_cur -= 1;
        }
        for j in lo_cur..hi_cur {
            if j == i {
                continue;
            }
            let w = norm_sqrt[i].min(norm_sqrt[j]);
            score += w;
        }
    }
    score
}

pub fn solve_backbone(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    tol_ppm: f64,
    top_k: usize,
) -> Vec<BackboneCandidate> {
    solve_backbone_min(peaks, precursor_neutral, precursor_z, tol_ppm, top_k, 2)
}

pub fn solve_backbone_min(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    tol_ppm: f64,
    top_k: usize,
    min_core_y: u8,
) -> Vec<BackboneCandidate> {
    use std::collections::HashMap;

    let min_core_y = min_core_y.max(1);

    if peaks.is_empty() {
        return Vec::new();
    }

    let base_peak_intensity = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0_f32, f32::max)
        .max(1.0_f32);

    let norm_sqrt: Vec<f64> = peaks
        .iter()
        .map(|&(_, i)| (i / base_peak_intensity).max(0.0).sqrt() as f64)
        .collect();

    const MIN_BB: f64 = 500.0;
    const MIN_GLYCAN: f64 = 406.0;

    let prec_z_max = precursor_z.max(1);
    let neutral = |mz: f64, z: f64| (mz - PROTON) * z;

    let mut rung_best: HashMap<(i64, u8), (f64, f64, u32)> = HashMap::new();

    const RUNG_WEIGHT: [f64; 6] = [2.0, 2.0, 1.5, 1.0, 1.0, 1.0];

    for (peak_idx, &(pmz, _pi)) in peaks.iter().enumerate() {
        let w = norm_sqrt[peak_idx];
        for z in 1..=prec_z_max {
            let pn = neutral(pmz, z as f64);
            if pn <= 0.0 || pn > precursor_neutral {
                continue;
            }
            for (ri, step) in std::iter::once(0.0)
                .chain(CORE_Y_STEPS.iter().copied())
                .enumerate()
            {
                let bb = pn - step;
                let bb_res = bb - H2O;
                if bb < MIN_BB || bb_res >= precursor_neutral {
                    continue;
                }
                if precursor_neutral - bb_res < MIN_GLYCAN {
                    continue;
                }
                let key = (bb * 100.0).round() as i64;
                let tol_key = ((bb * tol_ppm / 1e6).max(0.01) * 100.0).round() as i64;
                let ri_u8 = ri.min(5) as u8;
                let rw = RUNG_WEIGHT[ri_u8 as usize];
                let weighted_w = w * rw;
                for k in (key - tol_key)..=(key + tol_key) {
                    let e = rung_best.entry((k, ri_u8)).or_insert((0.0, 0.0, 0));
                    if weighted_w > e.0 {
                        e.0 = weighted_w;
                    }
                    e.1 += bb;
                    e.2 += 1;
                }
            }
        }
    }

    let mut bins: HashMap<i64, (u32, [bool; 6], f64, f64)> = HashMap::new();

    for ((k, ri), (best_w, mass_sum, vote_count)) in rung_best {
        let e = bins.entry(k).or_insert((0, [false; 6], 0.0, 0.0));
        e.0 += vote_count;
        if (ri as usize) < 6 {
            e.1[ri as usize] = true;
        }
        e.2 += mass_sum;
        e.3 += best_w;
    }

    let sorted_peaks: Vec<(f64, f32)>;
    let sorted_norm: Vec<f64>;
    let (sp, sn): (&[(f64, f32)], &[f64]) = {
        let already_sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
        if already_sorted {
            (peaks, &norm_sqrt)
        } else {
            let mut idx: Vec<usize> = (0..peaks.len()).collect();
            idx.sort_by(|&a, &b| {
                peaks[a]
                    .0
                    .partial_cmp(&peaks[b].0)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            sorted_peaks = idx.iter().map(|&i| peaks[i]).collect();
            sorted_norm = idx.iter().map(|&i| norm_sqrt[i]).collect();
            (&sorted_peaks, &sorted_norm)
        }
    };

    let mut cands: Vec<BackboneCandidate> = bins
        .into_iter()
        .filter_map(|(_k, (v, hits, mass_sum, int_score))| {
            let core_y_hits = hits.iter().filter(|&&h| h).count() as u8;
            if core_y_hits < min_core_y {
                return None;
            }
            let bb = mass_sum / v as f64;
            let cscore = complement_score(sp, sn, bb, tol_ppm);
            if core_y_hits >= 2 || cscore > 0.0 {
                Some(BackboneCandidate {
                    backbone_mass: bb,
                    core_y_hits,
                    votes: v,
                    intensity_score: int_score,
                    complement_score: cscore,
                })
            } else {
                None
            }
        })
        .collect();

    const COMPLEMENT_WEIGHT: f64 = 0.3;
    cands.sort_by(|a, b| {
        b.core_y_hits
            .cmp(&a.core_y_hits)
            .then({
                let sa = a.intensity_score + a.complement_score * COMPLEMENT_WEIGHT;
                let sb = b.intensity_score + b.complement_score * COMPLEMENT_WEIGHT;
                sb.partial_cmp(&sa).unwrap_or(std::cmp::Ordering::Equal)
            })
            .then(
                a.backbone_mass
                    .partial_cmp(&b.backbone_mass)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
    });

    cands.dedup_by(|cur, kept| {
        if (kept.backbone_mass - cur.backbone_mass).abs() < 0.05 {
            let total = kept.votes as f64 + cur.votes as f64;
            kept.backbone_mass = (kept.backbone_mass * kept.votes as f64
                + cur.backbone_mass * cur.votes as f64)
                / total;
            kept.votes += cur.votes;
            kept.intensity_score += cur.intensity_score;
            kept.complement_score = kept.complement_score.max(cur.complement_score);
            kept.core_y_hits = kept.core_y_hits.max(cur.core_y_hits);
            true
        } else {
            false
        }
    });
    cands.truncate(top_k);
    cands
}

#[derive(Clone, Copy)]
pub struct SpectrumStats {
    pub base: f64,
    pub sorted: bool,
}

impl SpectrumStats {
    #[inline]
    pub fn new(peaks: &[(f64, f32)]) -> Self {
        let base = peaks
            .iter()
            .map(|&(_, i)| i)
            .fold(0.0f32, f32::max)
            .max(1e-9) as f64;
        let sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
        Self { base, sorted }
    }
}

static Y_MAX_CHARGE: std::sync::OnceLock<u8> = std::sync::OnceLock::new();

pub fn init_y_max_charge(zmax: u8) {
    let _ = Y_MAX_CHARGE.set(zmax.max(1));
}

fn y_hicharge_enabled() -> bool {
    *Y_MAX_CHARGE.get().unwrap_or(&3) > 3
}

#[inline]
fn max_intensity_in_window(peaks: &[(f64, f32)], sorted: bool, lo: f64, hi: f64) -> f32 {
    if sorted {
        let start = peaks.partition_point(|&(m, _)| m < lo);
        peaks[start..]
            .iter()
            .take_while(|&&(m, _)| m <= hi)
            .map(|&(_, i)| i)
            .fold(0.0f32, f32::max)
    } else {
        peaks
            .iter()
            .filter(|&&(m, _)| m >= lo && m <= hi)
            .map(|&(_, i)| i)
            .fold(0.0f32, f32::max)
    }
}

#[inline]
pub fn best_frag_intensity(
    peaks: &[(f64, f32)],
    sorted: bool,
    neutral: f64,
    tol_ppm: f64,
    max_charge: u8,
) -> f32 {
    const ISOTOPE: f64 = 1.003_354_83;
    let hi_charge = y_hicharge_enabled();
    let zmax = if hi_charge {
        max_charge.clamp(1, 5)
    } else {
        max_charge.clamp(1, 3)
    };
    let mut best = 0.0f32;
    for z in 1..=zmax {
        let ion_mz = (neutral + z as f64 * PROTON) / z as f64;
        let tol = (ion_mz * tol_ppm / 1e6).max(0.01);
        let found = max_intensity_in_window(peaks, sorted, ion_mz - tol, ion_mz + tol);
        if hi_charge && z > 3 && found > 0.0 {
            let iso_mz = ion_mz + ISOTOPE / z as f64;
            let iso_tol = (iso_mz * tol_ppm / 1e6).max(0.01);
            if max_intensity_in_window(peaks, sorted, iso_mz - iso_tol, iso_mz + iso_tol) <= 0.0 {
                continue;
            }
        }
        best = best.max(found);
    }
    best
}

pub fn core_y_intensity(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb: f64,
    tol_ppm: f64,
    max_charge: u8,
    core: GlycanCore,
) -> f64 {
    let steps = core_y_steps(core);
    let mut score =
        best_frag_intensity(peaks, stats.sorted, bb, tol_ppm, max_charge) as f64 / stats.base;
    for &step in steps {
        score += best_frag_intensity(peaks, stats.sorted, bb + step, tol_ppm, max_charge) as f64
            / stats.base;
    }
    score
}

pub fn partial_glycan_by_intensity(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    residues: &[f64],
    tol_ppm: f64,
    max_charge: u8,
) -> f64 {
    use crate::glycan_mass::{HEX, HEXNAC};
    let n = residues.len();
    if n < 2 {
        return 0.0;
    }
    let total: f64 = residues.iter().sum();
    let core_adds = [
        HEXNAC,
        2.0 * HEXNAC,
        2.0 * HEXNAC + HEX,
        2.0 * HEXNAC + 2.0 * HEX,
        2.0 * HEXNAC + 3.0 * HEX,
    ];
    let mut score = 0.0f64;
    let mut prefix = 0.0f64;
    for i in 1..n {
        prefix += residues[i - 1];
        let b = prefix;
        let y = total - prefix + H2O;
        for &add in &core_adds {
            score += best_frag_intensity(peaks, stats.sorted, b + add, tol_ppm, max_charge) as f64
                / stats.base;
            score += best_frag_intensity(peaks, stats.sorted, y + add, tol_ppm, max_charge) as f64
                / stats.base;
        }
    }
    score
}

pub fn glycan_y_intensity(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb_neutral: f64,
    comp: &crate::glycan_db::GlycanComp,
    tol_ppm: f64,
    max_charge: u8,
) -> f64 {
    let retained = comp.y_ion_retained_masses();

    let match_int = |neutral: f64| -> f64 {
        best_frag_intensity(peaks, stats.sorted, neutral, tol_ppm, max_charge) as f64 / stats.base
    };

    let mut score = match_int(bb_neutral);
    let mut n_rungs = 1usize;
    for r in &retained {
        score += match_int(bb_neutral + r);
        n_rungs += 1;
    }
    ladder_norm_scale(score, n_rungs)
}

const CORE_RUNGS: usize = 6;

fn ladder_norm_enabled() -> bool {
    use std::sync::OnceLock;
    static CELL: OnceLock<bool> = OnceLock::new();
    *CELL.get_or_init(|| true)
}

fn ladder_norm_scale(score: f64, n_rungs: usize) -> f64 {
    if ladder_norm_enabled() && n_rungs > 0 {
        score * (CORE_RUNGS as f64) / (n_rungs as f64)
    } else {
        score
    }
}

#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[allow(clippy::too_many_arguments)]
pub fn glycan_y_intensity_decoy(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb_neutral: f64,
    comp: &crate::glycan_db::GlycanComp,
    tol_ppm: f64,
    max_charge: u8,
    seed: u64,
) -> f64 {
    let retained = comp.y_ion_retained_masses();

    let match_int = |neutral: f64| -> f64 {
        best_frag_intensity(peaks, stats.sorted, neutral, tol_ppm, max_charge) as f64 / stats.base
    };

    let mut score = match_int(bb_neutral);
    let last = retained.len().saturating_sub(1);
    for (ai, &r) in retained.iter().enumerate() {
        let shift = if ai == 0 || ai == last {
            0.0
        } else {
            let h = splitmix64(seed ^ (ai as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            1.0 + (h % 29_000) as f64 / 1000.0
        };
        score += match_int(bb_neutral + r + shift);
    }
    ladder_norm_scale(score, 1 + retained.len())
}

#[allow(clippy::too_many_arguments)]
pub fn glycan_y_hit_frac(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb_neutral: f64,
    comp: &crate::glycan_db::GlycanComp,
    tol_ppm: f64,
    max_charge: u8,
    decoy_seed: Option<u64>,
) -> f64 {
    let retained = comp.y_ion_retained_masses();
    let hit = |neutral: f64| -> bool {
        best_frag_intensity(peaks, stats.sorted, neutral, tol_ppm, max_charge) > 0.0
    };
    let total = retained.len() + 1;
    let mut hits = usize::from(hit(bb_neutral));
    let last = retained.len().saturating_sub(1);
    for (ai, &r) in retained.iter().enumerate() {
        let shift = match decoy_seed {
            Some(seed) if ai != 0 && ai != last => {
                let h = splitmix64(seed ^ (ai as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
                1.0 + (h % 29_000) as f64 / 1000.0
            }
            _ => 0.0,
        };
        if hit(bb_neutral + r + shift) {
            hits += 1;
        }
    }
    hits as f64 / total as f64
}

pub fn y0y1_anchor_intensity(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb_neutral: f64,
    prec_z: u8,
    tol_ppm: f64,
) -> f64 {
    use crate::glycan_mass::HEXNAC;
    let match_int = |ion_mz: f64| -> f64 {
        let tol = (ion_mz * tol_ppm / 1e6).max(0.01);
        let (lo, hi) = (ion_mz - tol, ion_mz + tol);
        let best = if stats.sorted {
            let start = peaks.partition_point(|&(m, _)| m < lo);
            peaks[start..]
                .iter()
                .take_while(|&&(m, _)| m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        } else {
            peaks
                .iter()
                .filter(|&&(m, _)| m >= lo && m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        };
        (best as f64) / stats.base
    };

    let zmax = prec_z.clamp(1, 2);
    let y0_neutral = bb_neutral;
    let y1_neutral = bb_neutral + HEXNAC;
    let mut score = 0.0;
    for z in 1..=zmax {
        let zf = z as f64;
        score += match_int((y0_neutral + zf * PROTON) / zf);
        score += match_int((y1_neutral + zf * PROTON) / zf);
    }
    score
}

pub fn count_core_y_hits(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb: f64,
    tol_ppm: f64,
    max_charge: u8,
    core: GlycanCore,
) -> u8 {
    let steps = core_y_steps(core);
    let mut hits: u8 = 0;
    if best_frag_intensity(peaks, stats.sorted, bb, tol_ppm, max_charge) > 0.0 {
        hits += 1;
    }
    for &step in steps {
        if best_frag_intensity(peaks, stats.sorted, bb + step, tol_ppm, max_charge) > 0.0 {
            hits += 1;
        }
    }
    hits
}

pub fn acceptor_core_y_gate(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb: f64,
    tol_ppm: f64,
    max_charge: u8,
    min_hits: u8,
    core: GlycanCore,
) -> bool {
    if count_core_y_hits(peaks, stats, bb, tol_ppm, max_charge, core) < min_hits {
        return false;
    }
    best_frag_intensity(
        peaks,
        stats.sorted,
        bb + core_y_steps(core)[0],
        tol_ppm,
        max_charge,
    ) > 0.0
}

#[cfg(test)]
mod tests {
    #[test]
    fn complement_score_window_follows_the_tolerance() {
        let bb = 1500.0;
        let target = bb + BY_COMPLEMENT_OFFSET;
        let peaks = vec![(300.0, 100.0), (target - 300.0 + 0.1, 100.0)];
        let norm: Vec<f64> = peaks.iter().map(|&(_, it)| (it as f64).sqrt()).collect();
        assert_eq!(complement_score(&peaks, &norm, bb, 20.0), 0.0);
        assert!(complement_score(&peaks, &norm, bb, 500.0) > 0.0);
    }

    #[test]
    fn complement_score_finds_the_by_pairs_of_a_known_peptide() {
        let proton = crate::glycan_mass::PROTON;
        let residues = [
            129.042593, 113.084064, 99.068414, 113.084064, 87.032028, 113.084064, 113.084064,
            99.068414, 129.042593, 87.032028, 128.094963,
        ];
        let residue_sum: f64 = residues.iter().sum();
        let neutral = residue_sum + H2O;

        let mut peaks: Vec<(f64, f32)> = Vec::new();
        for i in 2..=6 {
            let prefix: f64 = residues[..i].iter().sum();
            let suffix: f64 = residues[i..].iter().sum();
            peaks.push((prefix + proton, 500.0));
            peaks.push((suffix + H2O + proton, 500.0));
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let norm: Vec<f64> = peaks.iter().map(|&(_, it)| (it as f64).sqrt()).collect();

        let score = complement_score(&peaks, &norm, neutral, 20.0);
        assert!(
            score > 0.0,
            "the b/y pairs of the true peptide must be found, got {score}"
        );
        assert_eq!(
            complement_score(&peaks, &norm, residue_sum, 20.0),
            0.0,
            "the water-free residue mass is not the peptide neutral mass"
        );
    }

    #[test]
    fn complement_score_two_pointer_matches_bruteforce() {
        use super::{complement_score, BY_COMPLEMENT_OFFSET};
        fn reference(peaks: &[(f64, f32)], norm_sqrt: &[f64], bb: f64) -> f64 {
            let target = bb + BY_COMPLEMENT_OFFSET;
            let tol = (target * 20e-6).max(0.02);
            let n = peaks.len();
            if n < 2 {
                return 0.0;
            }
            let (b_lo, b_hi) = (50.0_f64, bb / 2.0);
            let mut score = 0.0_f64;
            for (i, &(b_mz, _)) in peaks.iter().enumerate() {
                if b_mz < b_lo {
                    continue;
                }
                if b_mz > b_hi {
                    break;
                }
                let y_mz = target - b_mz;
                if y_mz <= b_mz || y_mz <= 0.0 {
                    continue;
                }
                let start = peaks.partition_point(|&(m, _)| m < y_mz - tol);
                let end = peaks.partition_point(|&(m, _)| m <= y_mz + tol);
                for j in start..end {
                    if j == i {
                        continue;
                    }
                    score += norm_sqrt[i].min(norm_sqrt[j]);
                }
            }
            score
        }
        let mut seed = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for &npk in &[2usize, 5, 50, 400] {
            for &bb in &[800.0_f64, 1500.0, 2600.0, 4200.0] {
                let mut peaks: Vec<(f64, f32)> = (0..npk)
                    .map(|_| {
                        let m = 40.0 + (next() % 5_000_000) as f64 / 1000.0;
                        (m, 1.0 + (next() % 1000) as f32)
                    })
                    .collect();
                let target = bb + BY_COMPLEMENT_OFFSET;
                for k in 1..=3 {
                    let b = 100.0 * k as f64;
                    peaks.push((b, 500.0));
                    peaks.push((target - b, 500.0));
                }
                peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
                let norm: Vec<f64> = peaks.iter().map(|&(_, it)| (it as f64).sqrt()).collect();
                let got = complement_score(&peaks, &norm, bb, 20.0);
                let want = reference(&peaks, &norm, bb);
                assert_eq!(
                    got.to_bits(),
                    want.to_bits(),
                    "npk={npk} bb={bb}: two-pointer {got} != bruteforce {want}"
                );
            }
        }
    }

    use super::*;

    #[test]
    fn core_y_matches_multiply_charged_ions() {
        let bb = 2000.0_f64;
        let mz2 = |neutral: f64| (neutral + 2.0 * PROTON) / 2.0;
        let mut peaks = vec![
            (mz2(bb), 100.0f32),
            (mz2(bb + CORE_Y_STEPS[0]), 100.0f32),
            (mz2(bb + CORE_Y_STEPS[1]), 100.0f32),
        ];
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let stats = SpectrumStats::new(&peaks);
        assert_eq!(
            count_core_y_hits(&peaks, &stats, bb, 20.0, 1, GlycanCore::NGlycan),
            0,
            "z=1-only must miss 2+ Y-ions"
        );
        assert_eq!(
            count_core_y_hits(&peaks, &stats, bb, 20.0, 2, GlycanCore::NGlycan),
            3,
            "z<=2 must match the 2+ Y-ions"
        );
        assert!(core_y_intensity(&peaks, &stats, bb, 20.0, 2, GlycanCore::NGlycan) > 0.0);
    }

    #[test]
    fn solve_finds_backbone_with_minimal_2hexnac_core_glycan() {
        use crate::glycan_mass::HEXNAC;
        let bb_neutral = 1500.0_f64;
        let bb_res = bb_neutral - H2O;
        let glycan = 2.0 * HEXNAC;
        let precursor_neutral = bb_res + glycan;
        let mut peaks = vec![
            (bb_neutral + PROTON, 1000.0f32),
            (bb_neutral + HEXNAC + PROTON, 900.0f32),
            (bb_neutral + 2.0 * HEXNAC + PROTON, 500.0f32),
        ];
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let cands = solve_backbone_min(&peaks, precursor_neutral, 2, 20.0, 50, 2);
        assert!(
            cands
                .iter()
                .any(|c| (c.backbone_mass - bb_neutral).abs() < 0.05),
            "minimal 2xHexNAc-core backbone dropped by the off-by-H2O gate; got {:?}",
            cands.iter().map(|c| c.backbone_mass).collect::<Vec<_>>()
        );
    }

    #[test]
    fn solve_backbone_recovers_known_mass_from_core_ladder() {
        let bb = 1500.0;
        let mut peaks: Vec<(f64, f32)> = crate::glycan_mass::CORE_Y_STEPS
            .iter()
            .map(|&s| (bb + s + crate::glycan_mass::PROTON, 50.0))
            .collect();
        peaks.push((bb + crate::glycan_mass::PROTON, 40.0));
        peaks.push((999.9, 5.0));
        let precursor = bb + 1444.53;
        let out = solve_backbone(&peaks, precursor, 2, 20.0, 5);
        assert!(!out.is_empty());
        assert!(
            (out[0].backbone_mass - bb).abs() < 0.05,
            "got {}",
            out[0].backbone_mass
        );
        assert!(out[0].core_y_hits >= 2);
    }

    #[test]
    fn solve_backbone_empty_without_core_quorum() {
        let peaks = vec![(700.0, 50.0_f32), (1234.5, 50.0)];
        assert!(solve_backbone(&peaks, 2500.0, 2, 20.0, 5).is_empty());
    }

    #[test]
    fn glycan_y_intensity_prefers_true_composition_over_same_mass_decoy() {
        use crate::glycan_db::GlycanComp;
        use crate::glycan_mass::{HEX, HEXNAC};
        let proton = PROTON;
        let bb = 1500.0_f64;

        let truth = GlycanComp::from_tree_str("(N(N(H(H)(H))))").unwrap();
        let mut cum = 0.0;
        let mut peaks: Vec<(f64, f32)> = vec![(bb + proton, 1000.0)];
        for m in [HEXNAC, HEXNAC, HEX, HEX, HEX] {
            cum += m;
            peaks.push((bb + cum + proton, 500.0));
        }
        peaks.push((277.0, 20.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let stats = SpectrumStats::new(&peaks);

        let decoy = GlycanComp::from_tree_str("(N(H(H)(H)(H)))").unwrap();

        let s_true = glycan_y_intensity(&peaks, &stats, bb, &truth, 20.0, 3);
        let s_decoy = glycan_y_intensity(&peaks, &stats, bb, &decoy, 20.0, 3);
        assert!(
            s_true > s_decoy,
            "true composition Y-ladder ({s_true}) must beat the decoy ({s_decoy})"
        );
        assert!(
            s_true > 0.0,
            "true composition must have positive Y-ladder intensity"
        );
    }

    #[test]
    fn y0y1_anchor_scores_and_discriminates() {
        use crate::glycan_mass::HEXNAC;
        let bb_a = 1500.0_f64;
        let bb_b = 1650.0_f64;
        let peaks: Vec<(f64, f32)> = vec![
            (bb_a + PROTON, 1000.0),
            (bb_a + HEXNAC + PROTON, 800.0),
            ((bb_a + 2.0 * PROTON) / 2.0, 500.0),
            (204.087, 300.0),
        ];
        let stats = SpectrumStats::new(&peaks);
        let a = y0y1_anchor_intensity(&peaks, &stats, bb_a, 2, 20.0);
        let b = y0y1_anchor_intensity(&peaks, &stats, bb_b, 2, 20.0);
        assert!(a > 0.0, "anchor must fire for the true backbone, got {a}");
        assert!(a > b, "anchor must discriminate: A ({a}) > B ({b})");

        let empty = vec![(204.087, 100.0), (366.14, 80.0)];
        let stats_empty = SpectrumStats::new(&empty);
        assert_eq!(
            y0y1_anchor_intensity(&empty, &stats_empty, bb_a, 2, 20.0),
            0.0
        );
    }

    #[test]
    fn glycan_y_intensity_decoy_scores_below_target_on_true_spectrum() {
        use crate::glycan_db::GlycanComp;
        use crate::glycan_mass::{HEX, HEXNAC};
        let proton = PROTON;
        let bb = 1500.0_f64;
        let truth = GlycanComp::from_tree_str("(N(N(H(H)(H))))").unwrap();
        let mut cum = 0.0;
        let mut peaks: Vec<(f64, f32)> = vec![(bb + proton, 1000.0)];
        for m in [HEXNAC, HEXNAC, HEX, HEX, HEX] {
            cum += m;
            peaks.push((bb + cum + proton, 500.0));
        }
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let stats = SpectrumStats::new(&peaks);

        let s_true = glycan_y_intensity(&peaks, &stats, bb, &truth, 20.0, 3);
        let s_decoy = glycan_y_intensity_decoy(&peaks, &stats, bb, &truth, 20.0, 3, 12345);
        assert!(
            s_decoy < s_true,
            "decoy ladder ({s_decoy}) must score below the target ({s_true})"
        );
        assert!(s_decoy > 0.0, "decoy keeps Y0/Y1 anchor, should be > 0");
        let s_decoy2 = glycan_y_intensity_decoy(&peaks, &stats, bb, &truth, 20.0, 3, 12345);
        assert!(
            (s_decoy - s_decoy2).abs() < 1e-12,
            "decoy must be deterministic"
        );
    }

    #[test]
    fn glycan_y_intensity_decoy_matches_multicharge_anchors() {
        use crate::glycan_db::GlycanComp;
        use crate::glycan_mass::{HEX, HEXNAC};
        let bb = 3200.0_f64;
        let truth = GlycanComp::from_tree_str("(N(N(H(H)(H))))").unwrap();
        let mz2 = |neutral: f64| (neutral + 2.0 * PROTON) / 2.0;
        let mut cum = 0.0;
        let mut peaks: Vec<(f64, f32)> = vec![(mz2(bb), 1000.0)];
        for m in [HEXNAC, HEXNAC, HEX, HEX, HEX] {
            cum += m;
            peaks.push((mz2(bb + cum), 500.0));
        }
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let stats = SpectrumStats::new(&peaks);

        let s_true = glycan_y_intensity(&peaks, &stats, bb, &truth, 20.0, 3);
        assert!(s_true > 0.0, "target must match the z=2 ladder");
        let s_decoy = glycan_y_intensity_decoy(&peaks, &stats, bb, &truth, 20.0, 3, 777);
        assert!(
            s_decoy > 0.0,
            "decoy must match multi-charge Y0/Y1 anchors (exchangeable null); got {s_decoy}"
        );
        assert!(
            s_decoy < s_true,
            "decoy ({s_decoy}) still below target ({s_true}): interior rungs shifted away"
        );
        let s_decoy_z1 = glycan_y_intensity_decoy(&peaks, &stats, bb, &truth, 20.0, 1, 777);
        assert!(
            s_decoy_z1 == 0.0,
            "z=1 matcher cannot see z=2 anchors; got {s_decoy_z1}"
        );
    }

    #[test]
    fn glycan_y_intensity_decoy_keeps_top_rung() {
        use crate::glycan_db::GlycanComp;
        use crate::glycan_mass::{HEX, HEXNAC};
        let bb = 1500.0_f64;
        let truth = GlycanComp::from_tree_str("(N(N(H(H)(H))))").unwrap();
        let full = 2.0 * HEXNAC + 3.0 * HEX;
        let hexnac = HEXNAC;
        let mut peaks: Vec<(f64, f32)> = vec![
            (bb + PROTON, 1000.0),
            (bb + hexnac + PROTON, 500.0),
            (bb + full + PROTON, 400.0),
        ];
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let stats = SpectrumStats::new(&peaks);
        let s_true = glycan_y_intensity(&peaks, &stats, bb, &truth, 20.0, 3);
        let s_decoy = glycan_y_intensity_decoy(&peaks, &stats, bb, &truth, 20.0, 3, 42);
        assert!(
            (s_true - s_decoy).abs() < 1e-12,
            "target ({s_true}) and decoy ({s_decoy}) must match when only kept rungs are present"
        );
    }

    #[test]
    fn core_y_intensity_and_count_use_neutral_convention() {
        let proton = PROTON;
        let steps = CORE_Y_STEPS;
        let neutral = 1500.0_f64;
        let mut peaks: Vec<(f64, f32)> = vec![(neutral + proton, 1000.0)];
        for &s in steps.iter() {
            peaks.push((neutral + s + proton, 500.0));
        }
        peaks.push((321.0, 10.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let stats = SpectrumStats::new(&peaks);

        assert_eq!(
            count_core_y_hits(&peaks, &stats, neutral, 20.0, 3, GlycanCore::NGlycan),
            6,
            "neutral must match all 6 rungs"
        );
        let inten = core_y_intensity(&peaks, &stats, neutral, 20.0, 3, GlycanCore::NGlycan);
        assert!(
            inten > 0.0,
            "core_y_intensity must be positive for a real ladder, got {inten}"
        );
        assert!(
            (inten - 3.5).abs() < 1e-6,
            "expected 3.5 normalised intensity, got {inten}"
        );

        let residue = neutral - H2O;
        assert!(
            count_core_y_hits(&peaks, &stats, residue, 20.0, 3, GlycanCore::NGlycan) < 6,
            "residue mass (H2O too low) must NOT match the neutral-anchored ladder"
        );
    }

    #[test]
    fn glycan_y_hit_frac_full_vs_partial_vs_decoy() {
        use crate::glycan_db::GlycanComp;
        use crate::glycan_mass::{HEX, HEXNAC};
        let bb = 1500.0_f64;
        let g = GlycanComp::from_tree_str("(N(N(H(H)(H))))").unwrap();
        let mut full: Vec<(f64, f32)> = vec![(bb + PROTON, 1000.0)];
        let mut cum = 0.0;
        for m in [HEXNAC, HEXNAC, HEX, HEX, HEX] {
            cum += m;
            full.push((bb + cum + PROTON, 500.0));
        }
        full.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let sf = SpectrumStats::new(&full);
        assert!(
            (glycan_y_hit_frac(&full, &sf, bb, &g, 20.0, 3, None) - 1.0).abs() < 1e-9,
            "full ladder → 1.0"
        );
        let p0 = vec![(bb + PROTON, 1000.0)];
        let s0 = SpectrumStats::new(&p0);
        let f0 = glycan_y_hit_frac(&p0, &s0, bb, &g, 20.0, 3, None);
        assert!((f0 - 1.0 / 5.0).abs() < 1e-9, "Y0 only → 1/5, got {f0}");
        let dec = glycan_y_hit_frac(&full, &sf, bb, &g, 20.0, 3, Some(123));
        assert!(
            dec < 1.0,
            "decoy ladder must miss shifted interior rungs, got {dec}"
        );
    }

    #[test]
    fn acceptor_core_y_gate_requires_min_hits_and_y1() {
        let bb = 1500.0_f64;
        let mk = |rungs: &[usize]| {
            let mut peaks: Vec<(f64, f32)> = rungs
                .iter()
                .map(|&k| {
                    let neutral = if k == 0 { bb } else { bb + CORE_Y_STEPS[k - 1] };
                    (neutral + PROTON, 500.0)
                })
                .collect();
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            peaks
        };
        let p = mk(&[0, 1, 2]);
        let s = SpectrumStats::new(&p);
        assert!(
            acceptor_core_y_gate(&p, &s, bb, 20.0, 3, 3, GlycanCore::NGlycan),
            "3 rungs incl Y1 must pass"
        );
        let p2 = mk(&[0, 1]);
        let s2 = SpectrumStats::new(&p2);
        assert!(
            !acceptor_core_y_gate(&p2, &s2, bb, 20.0, 3, 3, GlycanCore::NGlycan),
            "<3 rungs must fail"
        );
        let p3 = mk(&[0, 2, 3]);
        let s3 = SpectrumStats::new(&p3);
        assert!(
            !acceptor_core_y_gate(&p3, &s3, bb, 20.0, 3, 3, GlycanCore::NGlycan),
            "missing Y1 must fail even with 3 rungs"
        );
    }

    #[test]
    fn quorum1_rejects_single_rung_without_complement() {
        let peaks = vec![
            (700.0, 50.0_f32),
            (911.3, 40.0),
            (1301.7, 45.0),
            (1502.9, 30.0),
        ];
        let out = solve_backbone_min(&peaks, 3000.0, 2, 20.0, 50, 1);
        assert!(
            out.is_empty(),
            "quorum-1 must reject single-rung candidates lacking b/y complement, got {out:?}"
        );
    }

    #[test]
    fn solve_backbone_reports_cluster_centroid_not_low_edge() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;
        let bb = 1500.005_f64;
        let mut peaks: Vec<(f64, f32)> = vec![(bb + proton, 50.0)];
        for &s in steps.iter() {
            peaks.push((bb + s + proton, 50.0));
        }
        let precursor = bb + 1444.53;
        let out = solve_backbone(&peaks, precursor, 2, 20.0, 5);
        assert!(!out.is_empty());
        assert!(
            (out[0].backbone_mass - bb).abs() < 0.002,
            "centroid off by {:.4} (got {})",
            (out[0].backbone_mass - bb).abs(),
            out[0].backbone_mass
        );
    }

    #[test]
    fn solve_backbone_is_deterministic_under_ties() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        let bb1 = 1200.0_f64;
        let bb2 = 1800.0_f64;

        let mut peaks: Vec<(f64, f32)> = Vec::new();
        for &bb in &[bb1, bb2] {
            peaks.push((bb + proton, 50.0));
            for &s in steps.iter() {
                peaks.push((bb + s + proton, 50.0));
            }
        }

        let precursor_neutral = bb2 + steps.last().copied().unwrap_or(0.0) + 500.0;

        let run1 = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);
        let run2 = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);

        assert!(!run1.is_empty(), "expected candidates in run1");
        assert_eq!(run1.len(), run2.len(), "result length differs between runs");
        for (c1, c2) in run1.iter().zip(run2.iter()) {
            assert_eq!(
                c1.core_y_hits, c2.core_y_hits,
                "core_y_hits differ between runs"
            );
            assert_eq!(c1.votes, c2.votes, "votes differ between runs");
            assert!(
                (c1.backbone_mass - c2.backbone_mass).abs() < 1e-9,
                "backbone_mass differs between runs: {} vs {}",
                c1.backbone_mass,
                c2.backbone_mass
            );
        }
    }

    #[test]
    fn solve_backbone_complement_pairs_promote_true_backbone() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        let true_bb = 1500.0_f64;
        let spurious_bb = 1300.0_f64;

        let intensity = 500.0_f32;
        let mut peaks: Vec<(f64, f32)> = Vec::new();

        peaks.push((true_bb + proton, intensity));
        for &s in steps.iter() {
            peaks.push((true_bb + s + proton, intensity));
        }

        peaks.push((spurious_bb + proton, intensity));
        for &s in steps.iter() {
            peaks.push((spurious_bb + s + proton, intensity));
        }

        let comp_target = true_bb + BY_COMPLEMENT_OFFSET;
        for &b_mz in &[300.0_f64, 450.0, 600.0] {
            let y_mz = comp_target - b_mz;
            if y_mz > 0.0 {
                peaks.push((b_mz, intensity));
                peaks.push((y_mz, intensity));
            }
        }

        let precursor_neutral = true_bb + 1444.53;

        let out = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);

        assert!(!out.is_empty(), "expected candidates");
        assert!(
            (out[0].backbone_mass - true_bb).abs() < 0.05,
            "expected true_bb ({}) first, got backbone_mass={:.4} complement_score={:.4}",
            true_bb,
            out[0].backbone_mass,
            out[0].complement_score,
        );
        assert!(
            out[0].complement_score > 0.0,
            "true_bb must have positive complement_score, got {}",
            out[0].complement_score
        );
    }

    #[test]
    fn solve_backbone_intensity_weighted_prefers_bright_ladder() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        let true_bb = 1500.0_f64;
        let mut peaks: Vec<(f64, f32)> = vec![(true_bb + proton, 10000.0)];
        for &s in steps.iter() {
            peaks.push((true_bb + s + proton, 8000.0));
        }

        let spurious_bb = 800.0_f64;
        peaks.push((spurious_bb + proton, 1.0));
        for &s in steps.iter() {
            peaks.push((spurious_bb + s + proton, 1.0));
        }

        let precursor_neutral = true_bb + 1444.53;
        let out = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);
        assert!(!out.is_empty());
        assert!(
            (out[0].backbone_mass - true_bb).abs() < 0.05,
            "expected true backbone first, got {}",
            out[0].backbone_mass
        );
    }

    #[test]
    fn o_core_y_intensity_matches_single_galnac_rung() {
        use crate::glycan_mass::{HEXNAC, PROTON};
        let bb = 1500.0_f64;
        let mut peaks = vec![(bb + PROTON, 1000.0f32), (bb + HEXNAC + PROTON, 900.0f32)];
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let stats = SpectrumStats::new(&peaks);
        assert_eq!(
            count_core_y_hits(&peaks, &stats, bb, 20.0, 1, GlycanCore::OGlycan),
            2,
            "O-glycan core-Y must count Y0 + the single GalNAc rung"
        );
        assert!(
            core_y_intensity(&peaks, &stats, bb, 20.0, 1, GlycanCore::OGlycan) > 0.0,
            "O-glycan core-Y intensity must be positive for the GalNAc rung"
        );
    }
}
