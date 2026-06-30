use crate::glycan_mass::{CORE_Y_STEPS, PROTON};

/// Water molecule mass (monoisotopic).
const H2O: f64 = 18.010565;

/// Constant for b/y complement-pair sum: for a backbone of neutral mass `bb`,
/// a singly-charged b-ion and its complementary singly-charged y-ion satisfy:
///   b_mz + y_mz = bb + H2O + 2*PROTON
/// so we check |p1_mz + p2_mz - (bb + BY_COMPLEMENT_OFFSET)| <= tol.
const BY_COMPLEMENT_OFFSET: f64 = H2O + 2.0 * PROTON; // 18.010565 + 2*1.0072765 = 20.025118

#[derive(Debug, Clone)]
pub struct BackboneCandidate {
    pub backbone_mass: f64,
    pub core_y_hits: u8,
    pub votes: u32,
    /// Sum of sqrt-compressed, base-peak-normalised intensities of all peaks
    /// that voted for this backbone cluster. Primary ranking after core_y_hits.
    pub intensity_score: f64,
    /// Complement-pair (b/y) confirmation score: sum of min(i1, i2) for peak
    /// pairs (p1_mz, p2_mz) where |p1_mz + p2_mz - (bb + BY_COMPLEMENT_OFFSET)|
    /// ≤ max(|bb + BY_COMPLEMENT_OFFSET| * 20 ppm, 0.02). Spurious clusters have
    /// few/no complement pairs; real backbones with b/y ladders score higher.
    pub complement_score: f64,
}

/// Count complement (b/y) pairs for a candidate backbone neutral mass `bb`.
///
/// Returns the sum of min(sqrt_norm_i1, sqrt_norm_i2) over unique singly-charged
/// (b, y) peak pairs where:
///   - b_mz ∈ [50, bb/2]  (singly-charged b-ion, below midpoint — avoids glycan-heavy region)
///   - y_mz = target - b_mz  (complementary y-ion)
///   - |b_mz + y_mz - target| ≤ tol
///
/// Restricting b_mz < bb/2 ensures each pair is counted only once and avoids the
/// large-mz glycan-ion region (glycan Y-ions dominate the upper half of the spectrum
/// and create incidental pair noise for every candidate backbone).
///
/// `peaks` must be sorted by m/z ascending.
/// `norm_sqrt` is the sqrt-compressed base-peak-normalised intensity array
/// (parallel to `peaks`).
fn complement_score(peaks: &[(f64, f32)], norm_sqrt: &[f64], bb: f64) -> f64 {
    let target = bb + BY_COMPLEMENT_OFFSET;
    // tolerance: ±20 ppm of the target sum, floor 0.02 Da
    let tol = (target * 20e-6).max(0.02);

    let n = peaks.len();
    if n < 2 {
        return 0.0;
    }

    // b-ion window: singly-charged peptide backbone b-ions are in [50, bb/2].
    // This cuts the glycan-heavy region (y-ion series > bb/2) which creates
    // noise pairs for every candidate regardless of whether it is the true backbone.
    let b_lo = 50.0_f64;
    let b_hi = bb / 2.0; // midpoint — each pair (b, y) has b < y when b < bb/2

    let mut score = 0.0_f64;

    // For each candidate b-ion (mz ∈ [b_lo, b_hi]), look up its complement y at
    // y_mz = target - b_mz using a binary search into the sorted peak list.
    for (i, &(b_mz, _)) in peaks.iter().enumerate() {
        if b_mz < b_lo {
            continue;
        }
        if b_mz > b_hi {
            break; // peaks are sorted; nothing further qualifies as b-ion
        }
        let y_mz = target - b_mz;
        // y_mz must be > b_mz (ensured by b_mz < bb/2 since target ≈ bb+20), > 0
        if y_mz <= b_mz || y_mz <= 0.0 {
            continue;
        }
        // Binary search for y_mz ± tol in sorted peaks (skip self)
        let y_lo = y_mz - tol;
        let y_hi = y_mz + tol;
        let start = peaks.partition_point(|&(m, _)| m < y_lo);
        let end = peaks.partition_point(|&(m, _)| m <= y_hi);
        for j in start..end {
            if j == i {
                continue; // self-pair guard (shouldn't happen given b<y constraint)
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
    use std::collections::HashMap;

    if peaks.is_empty() {
        return Vec::new();
    }

    // --- intensity normalisation: sqrt-compress then normalise to base peak ---
    // sqrt compression prevents a single huge peak from dominating (e.g. oxonium
    // at m/z 204 that happens to vote for a spurious low-mass backbone), while
    // still weighting true Y-ions (which tend to be higher-intensity) over noise.
    let base_peak_intensity = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0_f32, f32::max)
        .max(1.0_f32); // avoid divide-by-zero on all-zero spectra

    // Precompute normalised sqrt intensities once.
    let norm_sqrt: Vec<f64> = peaks
        .iter()
        .map(|&(_, i)| (i / base_peak_intensity).max(0.0).sqrt() as f64)
        .collect();

    // Minimum plausible backbone mass: real N-glycopeptides are >=500 Da
    // (smallest tryptic glycopeptide with NXS/T sequon). This cuts the very common
    // spurious cluster at ~203 Da (HexNAc oxonium m/z 204 at z=1 votes as Y0).
    const MIN_BB: f64 = 500.0;
    // Minimum implied glycan mass for N-glycopeptides: the core requires at least
    // 2×HexNAc (~406 Da). This cuts spurious high-mass backbone candidates that
    // leave no room for a real glycan.
    const MIN_GLYCAN: f64 = 406.0; // 2×HexNAc = minimum N-glycan core

    let prec_z_max = precursor_z.max(1);
    let neutral = |mz: f64, z: f64| (mz - PROTON) * z;

    // Per bin key, accumulate:
    //   .0  raw vote count (u32)             — for centroid/dedup
    //   .1  per-rung best-intensity seen      — [f64; 6], one slot per Y-ion rung
    //   .2  vote-weighted mass sum (f64)      — for centroid
    //   .3  sum of normalised sqrt intensities — intensity ranking score
    //
    // KEY CHANGE: each rung slot records the MAXIMUM intensity that voted for it
    // (de-duplicated across charges). This means a single peak seen at z=2 and z=3
    // does NOT double-count; only its best weight contributes to the intensity score.
    // This collapses the spurious inflation of low-mass clusters.
    //
    // Per-key per-rung accumulation: we first collect per (key, rung) the best
    // weight from any charge, then aggregate into the bucket.
    //
    // Data layout: outer HashMap key = backbone bin key; inner = rung index (0–5),
    // value = (best_w_for_rung, mass_accum, vote_count). We compute the final
    // intensity_score as the sum of per-rung best weights.

    // Two-pass: collect best-weight per (backbone_bin, rung) then aggregate.
    // We use a flat HashMap keyed by (backbone_bin, rung_idx) for efficiency.
    let mut rung_best: HashMap<(i64, u8), (f64 /*best_w*/, f64 /*mass for centroid*/, u32 /*vote count*/)> =
        HashMap::new();

    // Per-rung weights: Y0 (bare backbone) and Y1 (first HexNAc) are the most
    // diagnostic evidence for a specific backbone mass; they are unique to the
    // backbone and not part of a repeated glycan motif. Y3-Y5 (Hex additions)
    // are less specific because hexose residues are common in many glycan
    // antennae structures and can arise coincidentally.
    //   ri: 0=Y0, 1=Y1, 2=Y2, 3=Y3, 4=Y4, 5=Y5
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
                // Plausibility gates:
                //  1. backbone must be within [MIN_BB, precursor_neutral)
                //  2. implied glycan (precursor − backbone) must be ≥ MIN_GLYCAN
                //     to avoid spurious candidates where the backbone consumes nearly
                //     the entire precursor mass, leaving no room for a real glycan.
                if bb < MIN_BB || bb >= precursor_neutral {
                    continue;
                }
                if precursor_neutral - bb < MIN_GLYCAN {
                    continue;
                }
                let key = (bb * 100.0).round() as i64;
                let tol_key = ((bb * tol_ppm / 1e6).max(0.01) * 100.0).round() as i64;
                let ri_u8 = ri.min(5) as u8;
                // Rung-weighted intensity: Y0 and Y1 are more diagnostic than Y3-Y5.
                let rw = RUNG_WEIGHT[ri_u8 as usize];
                let weighted_w = w * rw;
                for k in (key - tol_key)..=(key + tol_key) {
                    let e = rung_best
                        .entry((k, ri_u8))
                        .or_insert((0.0, 0.0, 0));
                    // keep the best rung-weighted intensity seen (charge-dedup)
                    if weighted_w > e.0 {
                        e.0 = weighted_w;
                    }
                    e.1 += bb; // for vote-weighted centroid
                    e.2 += 1;
                }
            }
        }
    }

    // Aggregate per (bin, rung) into per-bin stats.
    // Per-bin stats: (raw_votes: u32, rung_hit_mask: [bool;6], mass_sum: f64, intensity_score: f64)
    let mut bins: HashMap<i64, (u32, [bool; 6], f64, f64)> = HashMap::new();

    for ((k, ri), (best_w, mass_sum, vote_count)) in rung_best {
        let e = bins.entry(k).or_insert((0, [false; 6], 0.0, 0.0));
        e.0 += vote_count;
        if (ri as usize) < 6 {
            e.1[ri as usize] = true;
        }
        e.2 += mass_sum;
        // intensity score = sum of per-rung best weights; each rung contributes
        // at most once (the best weight seen across all charges + all peaks
        // that could have generated this rung).
        e.3 += best_w;
    }

    // Sort peaks by m/z for the two-pointer complement-score sweep.
    // Note: peaks passed in may not be sorted; we need a sorted copy.
    // Avoid cloning the whole peaks vec if already sorted (most mzML parsers sort).
    let sorted_peaks: Vec<(f64, f32)>;
    let sorted_norm: Vec<f64>;
    let (sp, sn): (&[(f64, f32)], &[f64]) = {
        let already_sorted = peaks
            .windows(2)
            .all(|w| w[0].0 <= w[1].0);
        if already_sorted {
            (peaks, &norm_sqrt)
        } else {
            let mut idx: Vec<usize> = (0..peaks.len()).collect();
            idx.sort_by(|&a, &b| peaks[a].0.partial_cmp(&peaks[b].0).unwrap_or(std::cmp::Ordering::Equal));
            sorted_peaks = idx.iter().map(|&i| peaks[i]).collect();
            sorted_norm = idx.iter().map(|&i| norm_sqrt[i]).collect();
            (&sorted_peaks, &sorted_norm)
        }
    };

    // Build candidates from bins.
    let mut cands: Vec<BackboneCandidate> = bins
        .into_iter()
        .map(|(_k, (v, hits, mass_sum, int_score))| {
            let bb = mass_sum / v as f64;
            let cscore = complement_score(sp, sn, bb);
            BackboneCandidate {
                backbone_mass: bb,
                core_y_hits: hits.iter().filter(|&&h| h).count() as u8,
                votes: v,
                intensity_score: int_score,
                complement_score: cscore,
            }
        })
        .filter(|c| c.core_y_hits >= 2)
        .collect(); // core-Y quorum

    // Sort: PRIMARY = core_y_hits (more distinct rungs = stronger evidence),
    // SECONDARY = combined score: intensity_score + complement_score * COMPLEMENT_WEIGHT
    //   (complement pairs from real b/y ladders break ties; spurious clusters have none),
    // TERTIARY = backbone_mass ascending (deterministic tiebreak).
    //
    // COMPLEMENT_WEIGHT tuned so complement evidence has meaningful influence but does
    // not override a large core_y_hits gap. A weight of 0.5 means a candidate needs
    // ~2 extra complement pairs per unit of intensity score deficit to win the tie.
    const COMPLEMENT_WEIGHT: f64 = 0.3;
    cands.sort_by(|a, b| {
        b.core_y_hits
            .cmp(&a.core_y_hits)
            .then({
                let sa = a.intensity_score + a.complement_score * COMPLEMENT_WEIGHT;
                let sb = b.intensity_score + b.complement_score * COMPLEMENT_WEIGHT;
                sb.partial_cmp(&sa)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then(
                a.backbone_mass
                    .partial_cmp(&b.backbone_mass)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
    });

    // Merge near-mass clusters into a single vote-weighted centroid. We keep the
    // higher-ranked representative (kept) and fold the removed neighbour (cur) into it.
    // dedup_by calls the closure as (cur, kept): `cur` is the element under
    // inspection and `kept` is the surviving representative. Returning true
    // removes `cur`; we first fold its votes/mass into `kept` so the survivor
    // reports the cluster CENTROID rather than its lowest-mass edge.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solve_backbone_recovers_known_mass_from_core_ladder() {
        let bb = 1500.0;
        let mut peaks: Vec<(f64, f32)> = crate::glycan_mass::CORE_Y_STEPS
            .iter()
            .map(|&s| (bb + s + crate::glycan_mass::PROTON, 50.0))
            .collect();
        peaks.push((bb + crate::glycan_mass::PROTON, 40.0)); // Y0 = bare backbone+H
        peaks.push((999.9, 5.0)); // noise
        let precursor = bb + 1444.53; // backbone + a HexNAc2Hex5 glycan (~1444.53)
        let out = solve_backbone(&peaks, precursor, 2, 20.0, 5);
        assert!(!out.is_empty());
        // dedup_by merges candidates within 0.05 Da; the surviving rep may be up to 0.05 Da off
        assert!((out[0].backbone_mass - bb).abs() < 0.05, "got {}", out[0].backbone_mass);
        assert!(out[0].core_y_hits >= 2);
    }

    #[test]
    fn solve_backbone_empty_without_core_quorum() {
        let peaks = vec![(700.0, 50.0_f32), (1234.5, 50.0)]; // no core-Y ladder
        assert!(solve_backbone(&peaks, 2500.0, 2, 20.0, 5).is_empty());
    }

    /// The reported backbone mass must be the vote-weighted CENTROID of a merged
    /// near-mass cluster, not its lowest-mass edge. We build a Y-ladder anchored
    /// at a backbone whose mass is offset by a non-integer amount from the 0.01-Da
    /// bin grid; the tolerance window spreads votes symmetrically across several
    /// bins. With a centroid, the recovered mass must land within ~1 mDa of the
    /// true center; a lowest-edge tiebreak would bias it low by the tol window.
    #[test]
    fn solve_backbone_reports_cluster_centroid_not_low_edge() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;
        // True backbone chosen to sit between 0.01-Da bins (…x.xx5 region).
        let bb = 1500.005_f64;
        let mut peaks: Vec<(f64, f32)> = vec![(bb + proton, 50.0)]; // Y0
        for &s in steps.iter() {
            peaks.push((bb + s + proton, 50.0));
        }
        let precursor = bb + 1444.53;
        let out = solve_backbone(&peaks, precursor, 2, 20.0, 5);
        assert!(!out.is_empty());
        // Centroid must be within 2 mDa of the true backbone (unbiased).
        // The old lowest-edge behaviour skewed low by up to the 0.03-Da window.
        assert!(
            (out[0].backbone_mass - bb).abs() < 0.002,
            "centroid off by {:.4} (got {})",
            (out[0].backbone_mass - bb).abs(),
            out[0].backbone_mass
        );
    }

    /// Regression test: solve_backbone must return the same ordered result on repeated calls.
    /// We build two synthetic peptide+core-Y ladders with backbone masses 1200.0 and 1800.0 Da,
    /// each yielding the same number of core-Y hits and votes, producing a tie on
    /// (core_y_hits, intensity_score). Without the backbone_mass tiebreaker the HashMap
    /// iteration order could permute the two candidates across runs.
    #[test]
    fn solve_backbone_is_deterministic_under_ties() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // Build identical core-Y ladders anchored at two different backbone masses.
        // Both backbones get exactly the same number of Y-ladder peaks → tied hits & scores.
        let bb1 = 1200.0_f64;
        let bb2 = 1800.0_f64;

        let mut peaks: Vec<(f64, f32)> = Vec::new();
        for &bb in &[bb1, bb2] {
            // Y0 (bare backbone)
            peaks.push((bb + proton, 50.0));
            // Y_r for each core step
            for &s in steps.iter() {
                peaks.push((bb + s + proton, 50.0));
            }
        }

        // Precursor must be larger than any peak neutral mass.
        let precursor_neutral = bb2 + steps.last().copied().unwrap_or(0.0) + 500.0;

        // Call solve_backbone twice with identical args.
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

    /// Complement-pair (b/y) confirmation: a candidate WITH complement-pair support
    /// must outrank an otherwise-equal candidate WITHOUT. We build two synthetic
    /// backbones at the same precursor:
    ///  - `true_bb`: full core-Y ladder (same rung count and intensity as spurious)
    ///    PLUS synthetic b/y fragment ion pairs that satisfy b+y ≈ bb+20.025118.
    ///  - `spurious_bb`: same core-Y ladder (identical rung count and intensity)
    ///    but NO complement b/y pairs.
    /// After solve_backbone, `true_bb` must be first (complement_score tips the tie).
    #[test]
    fn solve_backbone_complement_pairs_promote_true_backbone() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // Two backbones at different masses; same number of core-Y rungs (6 each)
        // and same intensity — so core_y_hits and intensity_score are tied.
        // Use intensity 500.0 for both so sqrt-norm is equal.
        let true_bb = 1500.0_f64;
        let spurious_bb = 1300.0_f64;

        let intensity = 500.0_f32;
        let mut peaks: Vec<(f64, f32)> = Vec::new();

        // Core-Y ladder for true_bb
        peaks.push((true_bb + proton, intensity)); // Y0
        for &s in steps.iter() {
            peaks.push((true_bb + s + proton, intensity));
        }

        // Core-Y ladder for spurious_bb (identical count and intensity)
        peaks.push((spurious_bb + proton, intensity)); // Y0
        for &s in steps.iter() {
            peaks.push((spurious_bb + s + proton, intensity));
        }

        // Add complement pairs only for true_bb.
        // A valid b-ion at b_mz and its complement y at y_mz = (true_bb + BY_COMPLEMENT_OFFSET) - b_mz.
        // Choose b_mz values far from the core-Y peaks to avoid confusion.
        let comp_target = true_bb + BY_COMPLEMENT_OFFSET; // 1500 + 20.025118 = 1520.025118
        // Three b/y pairs: b at 300, 450, 600 Da → y = comp_target - b
        for &b_mz in &[300.0_f64, 450.0, 600.0] {
            let y_mz = comp_target - b_mz;
            // Only add if y_mz > 0 and distinct from existing peaks
            if y_mz > 0.0 {
                peaks.push((b_mz, intensity));
                peaks.push((y_mz, intensity));
            }
        }

        // Precursor must be large enough: true_bb + glycan (1444.53 ≥ MIN_GLYCAN=406)
        let precursor_neutral = true_bb + 1444.53;

        // Both spurious_bb (1300) and true_bb (1500) are within [MIN_BB, precursor_neutral)
        // and leave room for MIN_GLYCAN:
        //   precursor_neutral - spurious_bb = 1500+1444.53 - 1300 = 1644.53 ≥ 406 ✓
        let out = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);

        assert!(!out.is_empty(), "expected candidates");
        // true_bb must be ranked first because its complement_score > 0 while spurious_bb = 0.
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

    /// High-intensity true Y-ions should rank above equal-rung-count low-intensity noise.
    /// Build a known backbone ladder with high intensity, and a spurious same-rung-count
    /// cluster built from very low-intensity peaks. The true backbone must rank first.
    #[test]
    fn solve_backbone_intensity_weighted_prefers_bright_ladder() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // True backbone at 1500 Da — bright peaks
        let true_bb = 1500.0_f64;
        let mut peaks: Vec<(f64, f32)> = vec![(true_bb + proton, 10000.0)]; // Y0 bright
        for &s in steps.iter() {
            peaks.push((true_bb + s + proton, 8000.0));
        }

        // Spurious backbone at 800 Da — dim peaks that by flat-vote coincidence
        // would accumulate the same rung count but much lower intensity.
        let spurious_bb = 800.0_f64;
        peaks.push((spurious_bb + proton, 1.0)); // Y0 dim
        for &s in steps.iter() {
            peaks.push((spurious_bb + s + proton, 1.0));
        }

        let precursor_neutral = true_bb + 1444.53;
        let out = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);
        assert!(!out.is_empty());
        // The bright true-backbone cluster must be ranked first.
        assert!(
            (out[0].backbone_mass - true_bb).abs() < 0.05,
            "expected true backbone first, got {}",
            out[0].backbone_mass
        );
    }
}
