use crate::glycan_mass::{CORE_Y_STEPS, PROTON};

#[derive(Debug, Clone)]
pub struct BackboneCandidate {
    pub backbone_mass: f64,
    pub core_y_hits: u8,
    pub votes: u32,
    /// Sum of sqrt-compressed, base-peak-normalised intensities of all peaks
    /// that voted for this backbone cluster. Primary ranking after core_y_hits.
    pub intensity_score: f64,
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

    // Build candidates from bins.
    let mut cands: Vec<BackboneCandidate> = bins
        .into_iter()
        .map(|(_k, (v, hits, mass_sum, int_score))| BackboneCandidate {
            backbone_mass: mass_sum / v as f64,
            core_y_hits: hits.iter().filter(|&&h| h).count() as u8,
            votes: v,
            intensity_score: int_score,
        })
        .filter(|c| c.core_y_hits >= 2)
        .collect(); // core-Y quorum

    // Sort: PRIMARY = core_y_hits (more distinct rungs = stronger evidence),
    // SECONDARY = intensity_score (sum of rung-weighted sqrt-normalised intensities,
    // charge-deduplicated — true Y-ions tend to be brighter than spurious matches),
    // TERTIARY = backbone_mass ascending (deterministic tiebreak).
    cands.sort_by(|a, b| {
        b.core_y_hits
            .cmp(&a.core_y_hits)
            .then(
                b.intensity_score
                    .partial_cmp(&a.intensity_score)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
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
