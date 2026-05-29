//! Chimeric Phase 3 — greedy shared-fragment competition.
//!
//! When `--chimeric` emits multiple distinct peptides per scan, PSM-level FDR
//! over the inflated multi-PSM set is structurally untrustworthy (chimeric
//! Phase 2 bench). The 2026-05-29 Astral overlap diagnostic confirmed the
//! "fragment-theft" premise: a substantial fraction of co-emitted runner-up
//! peptides match the SAME MS2 peaks as the more-confident top peptide.
//!
//! This module implements the discriminator: peptides compete for fragment
//! peaks **most-confident-first** (greedy, à la Yu et al. 2023). Each peptide
//! claims its matched peaks; a later, less-confident peptide is credited only
//! for the peaks not already claimed ("unique" evidence). The unique-evidence
//! metrics here (a) feed the residual SpecEValue re-score in `match_engine`
//! (the in-engine discriminator) and (b) are emitted as additive PIN columns
//! for Percolator.
//!
//! This file holds the **pure** core (set arithmetic over peak keys +
//! intensities), independent of the scorer; the extraction of a peptide's
//! matched (peak-key, intensity) list from a `ScoredSpectrum` and the residual
//! re-score live in `match_engine` where the scoring context is available.

use rustc_hash::FxHashSet;

/// A peptide's matched charge-1 b/y peaks: `(quantized m/z key, intensity)`.
/// The key is `round(peak_mz * 1000)` — identical to `matched_peak_keys` and
/// the validated overlap diagnostic, so a peak claimed by one peptide is
/// recognised as the same peak by the next.
pub(crate) type MatchedPeak = (i64, f32);

/// Unique-fragment evidence for one peptide, computed against the set of peaks
/// already claimed by more-confident peptides on the same scan.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct UniqueEvidence {
    /// Count of matched peaks NOT already claimed.
    pub unique_matched_ions: u32,
    /// Σ intensity(unique peaks) / Σ intensity(all matched peaks). 0.0 if the
    /// peptide matched no peaks.
    pub unique_explained_fraction: f32,
    /// |matched ∩ claimed| / |matched|. 0.0 if the peptide matched no peaks.
    pub shared_frac_claimed: f32,
}

impl UniqueEvidence {
    /// Evidence for a peptide that matched no peaks (degenerate).
    pub(crate) const EMPTY: UniqueEvidence = UniqueEvidence {
        unique_matched_ions: 0,
        unique_explained_fraction: 0.0,
        shared_frac_claimed: 0.0,
    };
}

/// Compute one peptide's unique-evidence metrics given the peaks already
/// claimed by more-confident peptides. Does **not** mutate `claimed` — the
/// caller claims afterwards via [`claim`] (so a peptide never competes against
/// its own peaks).
pub(crate) fn unique_evidence(matched: &[MatchedPeak], claimed: &FxHashSet<i64>) -> UniqueEvidence {
    if matched.is_empty() {
        return UniqueEvidence::EMPTY;
    }
    let mut unique_count: u32 = 0;
    let mut shared_count: u32 = 0;
    let mut unique_intensity: f64 = 0.0;
    let mut total_intensity: f64 = 0.0;
    for &(key, intensity) in matched {
        let i = intensity as f64;
        total_intensity += i;
        if claimed.contains(&key) {
            shared_count += 1;
        } else {
            unique_count += 1;
            unique_intensity += i;
        }
    }
    let n = matched.len() as f32;
    UniqueEvidence {
        unique_matched_ions: unique_count,
        unique_explained_fraction: if total_intensity > 0.0 {
            (unique_intensity / total_intensity) as f32
        } else {
            0.0
        },
        shared_frac_claimed: shared_count as f32 / n,
    }
}

/// Insert all of a peptide's matched peaks into the claimed set.
pub(crate) fn claim(matched: &[MatchedPeak], claimed: &mut FxHashSet<i64>) {
    claimed.extend(matched.iter().map(|&(key, _)| key));
}

/// Run the full greedy competition over peptides already ordered
/// most-confident-first. Returns one [`UniqueEvidence`] per input peptide,
/// aligned by index. The first peptide (rank-1) always sees an empty claimed
/// set, so its evidence reflects the full spectrum (unique == all matched,
/// `shared_frac_claimed == 0`). Pure; used directly by the unit tests and
/// mirrored by the `match_engine` hook (which additionally re-scores).
#[cfg(test)]
pub(crate) fn compete(ordered: &[Vec<MatchedPeak>]) -> Vec<UniqueEvidence> {
    let mut claimed: FxHashSet<i64> = FxHashSet::default();
    let mut out = Vec::with_capacity(ordered.len());
    for matched in ordered {
        out.push(unique_evidence(matched, &claimed));
        claim(matched, &mut claimed);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(keys: &[(i64, f32)]) -> Vec<MatchedPeak> {
        keys.to_vec()
    }

    #[test]
    fn rank1_sees_full_spectrum() {
        // The first (most-confident) peptide competes against nothing.
        let a = p(&[(1, 10.0), (2, 20.0), (3, 30.0)]);
        let ev = compete(&[a]);
        assert_eq!(ev[0].unique_matched_ions, 3);
        assert_eq!(ev[0].shared_frac_claimed, 0.0);
        assert_eq!(ev[0].unique_explained_fraction, 1.0);
    }

    #[test]
    fn subset_theft_leaves_runner_up_with_no_unique() {
        // B's peaks ⊂ A's peaks: after A claims, B has nothing of its own.
        let a = p(&[(1, 10.0), (2, 20.0), (3, 30.0)]);
        let b = p(&[(1, 10.0), (2, 20.0)]);
        let ev = compete(&[a, b]);
        assert_eq!(ev[1].unique_matched_ions, 0);
        assert_eq!(ev[1].shared_frac_claimed, 1.0);
        assert_eq!(ev[1].unique_explained_fraction, 0.0);
    }

    #[test]
    fn disjoint_runner_up_keeps_all_its_evidence() {
        let a = p(&[(1, 10.0), (2, 20.0), (3, 30.0)]);
        let c = p(&[(4, 40.0), (5, 50.0)]);
        let ev = compete(&[a, c]);
        assert_eq!(ev[1].unique_matched_ions, 2);
        assert_eq!(ev[1].shared_frac_claimed, 0.0);
        assert_eq!(ev[1].unique_explained_fraction, 1.0);
    }

    #[test]
    fn partial_overlap_splits_unique_and_shared() {
        // D shares peak 3 with A, uniquely owns peak 4.
        let a = p(&[(1, 10.0), (2, 20.0), (3, 30.0)]);
        let d = p(&[(3, 30.0), (4, 40.0)]);
        let ev = compete(&[a, d]);
        assert_eq!(ev[1].unique_matched_ions, 1);
        assert_eq!(ev[1].shared_frac_claimed, 0.5);
        // unique intensity 40 / total matched 70.
        assert!((ev[1].unique_explained_fraction - (40.0 / 70.0)).abs() < 1e-6);
    }

    #[test]
    fn greedy_is_transitive_across_three_peptides() {
        // A claims {1,2}; B claims {2,3} but 2 already taken → unique {3};
        // C matches {1,3} both already claimed → 0 unique.
        let a = p(&[(1, 10.0), (2, 20.0)]);
        let b = p(&[(2, 20.0), (3, 30.0)]);
        let c = p(&[(1, 10.0), (3, 30.0)]);
        let ev = compete(&[a, b, c]);
        assert_eq!(ev[0].unique_matched_ions, 2);
        assert_eq!(ev[1].unique_matched_ions, 1); // only peak 3
        assert_eq!(ev[2].unique_matched_ions, 0); // 1 and 3 both claimed
        assert_eq!(ev[2].shared_frac_claimed, 1.0);
    }

    #[test]
    fn no_matched_peaks_is_empty_evidence() {
        let ev = unique_evidence(&[], &FxHashSet::default());
        assert_eq!(ev, UniqueEvidence::EMPTY);
    }

    #[test]
    fn identical_evidence_is_order_symmetric_for_target_and_decoy() {
        // A target runner-up and a decoy runner-up with identical matched peaks
        // get identical evidence when competing after the same rank-1 — the
        // decoy-symmetry property the FDR relies on.
        let a = p(&[(1, 10.0), (2, 20.0)]);
        let target_b = p(&[(2, 20.0), (3, 30.0)]);
        let decoy_b = p(&[(2, 20.0), (3, 30.0)]);
        let ev_t = compete(&[a.clone(), target_b]);
        let ev_d = compete(&[a, decoy_b]);
        assert_eq!(ev_t[1], ev_d[1]);
    }
}
