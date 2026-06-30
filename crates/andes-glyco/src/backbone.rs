use crate::glycan_mass::{CORE_Y_STEPS, PROTON};

#[derive(Debug, Clone)]
pub struct BackboneCandidate {
    pub backbone_mass: f64,
    pub core_y_hits: u8,
    pub votes: u32,
}

pub fn solve_backbone(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    tol_ppm: f64,
    top_k: usize,
) -> Vec<BackboneCandidate> {
    use std::collections::HashMap;
    // Bucket votes for a backbone mass at 0.01-Da resolution.
    // Per key we track: vote count, which core-Y steps contributed, and the
    // vote-weighted sum of the *true* (unrounded) backbone-mass estimates that
    // landed in this bin. The weighted sum lets us report the CENTROID mass of a
    // merged near-mass cluster rather than the lowest-mass edge — without the
    // centroid the reported mass is biased low by up to the match tolerance,
    // which would invalidate the ±20 ppm searchable-backbone gate.
    let mut votes: HashMap<i64, (u32, [bool; 6], f64)> = HashMap::new();
    let neutral = |mz: f64, z: f64| (mz - PROTON) * z; // peak neutral mass at charge z
    for &(pmz, _pi) in peaks {
        for z in 1..=precursor_z.max(1) {
            let pn = neutral(pmz, z as f64);
            if pn <= 0.0 || pn > precursor_neutral {
                continue;
            }
            // this peak could be Y0 (bare backbone) or Y_r (backbone + core step r)
            // candidate backbone = pn - core_step (Y0: step 0)
            for (ri, step) in std::iter::once(0.0).chain(CORE_Y_STEPS.iter().copied()).enumerate() {
                let bb = pn - step;
                if bb <= 0.0 {
                    continue;
                }
                let key = (bb * 100.0).round() as i64;
                let tol_key = ((bb * tol_ppm / 1e6).max(0.01) * 100.0).round() as i64;
                // accumulate into nearby keys within tolerance
                for k in (key - tol_key)..=(key + tol_key) {
                    let e = votes.entry(k).or_insert((0, [false; 6], 0.0));
                    e.0 += 1;
                    if ri < 6 {
                        e.1[ri] = true;
                    }
                    // accumulate the unrounded backbone estimate, weighted by 1
                    // vote, so backbone_mass = mass_sum / votes is the centroid.
                    e.2 += bb;
                }
            }
        }
    }
    // Per-key candidates: backbone_mass is the vote-weighted centroid (mass_sum /
    // votes), not the bin's rounded center, so it is unbiased within the bin.
    let mut cands: Vec<BackboneCandidate> = votes
        .into_iter()
        .map(|(_k, (v, hits, mass_sum))| BackboneCandidate {
            backbone_mass: mass_sum / v as f64,
            core_y_hits: hits.iter().filter(|&&h| h).count() as u8,
            votes: v,
        })
        .filter(|c| c.core_y_hits >= 2)
        .collect(); // core-Y quorum
    cands.sort_by(|a, b| {
        b.core_y_hits.cmp(&a.core_y_hits)
            .then(b.votes.cmp(&a.votes))
            .then(a.backbone_mass.partial_cmp(&b.backbone_mass).unwrap_or(std::cmp::Ordering::Equal))
    });
    // Merge near-mass clusters into a single vote-weighted centroid. We keep the
    // higher-ranked representative (a) and fold the merged neighbour's (b) votes
    // and mass into it so the survivor reports the cluster CENTROID rather than
    // its lowest-mass edge. dedup_by removes b and keeps the mutated a.
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
    /// (core_y_hits, votes). Without the backbone_mass tiebreaker the HashMap iteration
    /// order could permute the two candidates across runs.
    #[test]
    fn solve_backbone_is_deterministic_under_ties() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // Build identical core-Y ladders anchored at two different backbone masses.
        // Both backbones get exactly the same number of Y-ladder peaks → tied hits & votes.
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
}
