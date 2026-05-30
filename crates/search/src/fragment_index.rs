//! Chimeric Phase-4 fragment-evidence prefilter: an inverted
//! `fragment-bin -> [candidate_id]` index used as a candidate generator under
//! `--chimeric`, so only candidates with real fragment evidence are scored.
//!
//! Task 3 wires the index into the chimeric hot path via `match_engine`.
//!
//! SUPERSEDED by `sage_index.rs` (Approach B): this Approach-A FragmentIndex /
//! FragmentVoter is the failed algorithm and is no longer wired into the search
//! (the `PreparedSearch.fragment_index` build is forced to `None`). Kept in the
//! tree as a record per the plan; `dead_code` is allowed until it is removed.
#![allow(dead_code)]

use crate::candidate_gen::Candidate;
use scoring_crate::scoring::fragment_ions::predict_by_ions;
use std::cmp::Ordering;

/// Inverted index: fragment-m/z bin -> candidate ids that have a charge-1 b/y
/// fragment in that bin. CSR layout (offsets + concatenated ids) keeps it
/// compact (`u32` ids) — memory was the failure mode of the abandoned Java
/// attempt, so this stays as tight as possible.
pub(crate) struct FragmentIndex {
    bin_width: f64,
    min_mz: f64,
    n_bins: usize,
    /// `bucket_offsets[b]..bucket_offsets[b+1]` indexes into `bucket_candidates`.
    bucket_offsets: Vec<u32>,
    bucket_candidates: Vec<u32>,
}

impl FragmentIndex {
    /// Build over the full candidate set (target+decoy, mod-expanded).
    /// `bin_width` is the fragment-m/z bin in Da (caller picks ~tolerance:
    /// 0.02 for high-res, 0.5 for low-res).
    pub(crate) fn build(candidates: &[Candidate], bin_width: f64) -> Self {
        // Pass A: bounds.
        let mut min_mz = f64::INFINITY;
        let mut max_mz = f64::NEG_INFINITY;
        for c in candidates {
            for ion in predict_by_ions(&c.peptide, 1..=1) {
                if ion.mz < min_mz {
                    min_mz = ion.mz;
                }
                if ion.mz > max_mz {
                    max_mz = ion.mz;
                }
            }
        }
        if !min_mz.is_finite() {
            return FragmentIndex {
                bin_width,
                min_mz: 0.0,
                n_bins: 0,
                bucket_offsets: vec![0],
                bucket_candidates: Vec::new(),
            };
        }
        let n_bins = (((max_mz - min_mz) / bin_width).floor() as usize) + 1;

        // Pass B: per-bin counts.
        let mut counts = vec![0u32; n_bins];
        let bin_of = |mz: f64| -> Option<usize> {
            if mz < min_mz {
                return None;
            }
            let b = ((mz - min_mz) / bin_width).floor() as usize;
            if b < n_bins {
                Some(b)
            } else {
                None
            }
        };
        for c in candidates {
            for ion in predict_by_ions(&c.peptide, 1..=1) {
                if let Some(b) = bin_of(ion.mz) {
                    counts[b] += 1;
                }
            }
        }

        // Prefix sum -> offsets.
        let mut bucket_offsets = vec![0u32; n_bins + 1];
        let mut acc = 0u32;
        for b in 0..n_bins {
            bucket_offsets[b] = acc;
            acc += counts[b];
        }
        bucket_offsets[n_bins] = acc;

        // Pass C: fill via a moving cursor copy of offsets.
        let mut cursor: Vec<u32> = bucket_offsets[..n_bins].to_vec();
        let mut bucket_candidates = vec![0u32; acc as usize];
        for (cid, c) in candidates.iter().enumerate() {
            for ion in predict_by_ions(&c.peptide, 1..=1) {
                if let Some(b) = bin_of(ion.mz) {
                    let pos = cursor[b] as usize;
                    bucket_candidates[pos] = cid as u32;
                    cursor[b] += 1;
                }
            }
        }

        FragmentIndex {
            bin_width,
            min_mz,
            n_bins,
            bucket_offsets,
            bucket_candidates,
        }
    }

    #[inline]
    fn bin_index(&self, mz: f64) -> Option<usize> {
        if self.n_bins == 0 || mz < self.min_mz {
            return None;
        }
        let b = ((mz - self.min_mz) / self.bin_width).floor() as usize;
        if b < self.n_bins {
            Some(b)
        } else {
            None
        }
    }

    /// Candidate ids whose charge-1 b/y fragment falls in the bin containing `mz`.
    /// (Callers also probe `mz ± bin_width` to cover tolerance at bin edges.)
    pub(crate) fn candidates_in_bin(&self, mz: f64) -> &[u32] {
        match self.bin_index(mz) {
            Some(b) => {
                let lo = self.bucket_offsets[b] as usize;
                let hi = self.bucket_offsets[b + 1] as usize;
                &self.bucket_candidates[lo..hi]
            }
            None => &[],
        }
    }

    /// Total indexed (fragment, candidate) entries — for memory accounting.
    pub(crate) fn n_entries(&self) -> usize {
        self.bucket_candidates.len()
    }
}

/// Per-thread reusable scratch for the per-spectrum vote/top-K step. `votes` is
/// sized to the candidate count; `touched` records which entries were written so
/// reset is O(touched), never O(n_candidates) (the Sage pattern — avoids the
/// global per-spectrum allocation that OOM'd the Java attempt).
pub(crate) struct FragmentVoter {
    votes: Vec<f32>,
    touched: Vec<u32>,
}

impl FragmentVoter {
    pub(crate) fn new(n_candidates: usize) -> Self {
        FragmentVoter {
            votes: vec![0.0; n_candidates],
            touched: Vec::with_capacity(4096),
        }
    }

    /// Accumulate one vote per matched fragment bin and return up to `k`
    /// in-window candidate ids ranked by vote (descending; ties broken by
    /// ascending id for determinism). `peaks` are `(rank, mz)` for the active
    /// observed peaks; `in_window(cid)` gates by precursor-mass eligibility.
    /// Probes the peak's bin and both neighbours to cover ±bin_width tolerance.
    pub(crate) fn top_k<F: Fn(u32) -> bool>(
        &mut self,
        idx: &FragmentIndex,
        peaks: &[(u32, f64)],
        in_window: F,
        k: usize,
    ) -> Vec<u32> {
        // Reset prior votes.
        for &c in &self.touched {
            self.votes[c as usize] = 0.0;
        }
        self.touched.clear();

        for &(_rank, mz) in peaks {
            // weight = 1.0 (matched-fragment count). Rank-weighting is a P4
            // tuning knob; count is a strong, deterministic baseline.
            for probe in [mz - idx.bin_width, mz, mz + idx.bin_width] {
                for &cid in idx.candidates_in_bin(probe) {
                    let v = &mut self.votes[cid as usize];
                    if *v == 0.0 {
                        self.touched.push(cid);
                    }
                    *v += 1.0;
                }
            }
        }

        // Collect in-window touched candidates with their votes, partial-sort top-k.
        let mut scored: Vec<(f32, u32)> = self
            .touched
            .iter()
            .copied()
            .filter(|&c| in_window(c))
            .map(|c| (self.votes[c as usize], c))
            .collect();
        scored.sort_by(|a, b| {
            b.0.partial_cmp(&a.0)
                .unwrap_or(Ordering::Equal)
                .then(a.1.cmp(&b.1))
        });
        scored.truncate(k);
        scored.into_iter().map(|(_, c)| c).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidate_gen::Candidate;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;

    fn cand(seq: &str) -> Candidate {
        let residues = seq
            .bytes()
            .map(|r| AminoAcid::standard(r).unwrap())
            .collect();
        Candidate {
            peptide: Peptide::new(residues, b'-', b'-'),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    #[test]
    fn build_indexes_every_candidate_fragment() {
        let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
        let idx = FragmentIndex::build(&cands, 0.02);
        // Every charge-1 b/y fragment of every candidate must be retrievable:
        // querying the bin of a known fragment m/z returns that candidate id.
        let frags = scoring_crate::scoring::fragment_ions::predict_by_ions(&cands[0].peptide, 1..=1);
        let probe = frags[0].mz;
        let hits = idx.candidates_in_bin(probe);
        assert!(
            hits.contains(&0u32),
            "candidate 0 must be indexed at its own fragment m/z"
        );
    }

    #[test]
    fn csr_fill_packs_colliding_candidates_and_counts_entries() {
        let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
        // Deliberately wide bin so fragments from BOTH candidates collide into
        // the same few bins, exercising the cursor-fill packing across cands.
        let idx = FragmentIndex::build(&cands, 1000.0);

        // A fragment m/z of candidate 0; with bin_width 1000 candidate 1's
        // nearby fragments share the bin, so the bin must contain BOTH ids.
        let frags0 = predict_by_ions(&cands[0].peptide, 1..=1);
        let probe = frags0[0].mz;
        let hits = idx.candidates_in_bin(probe);
        let both_here = hits.contains(&0u32) && hits.contains(&1u32);
        // Fallback: at minimum some bin across the whole index packs both ids.
        let both_somewhere = (0..idx.n_bins).any(|b| {
            let lo = idx.bucket_offsets[b] as usize;
            let hi = idx.bucket_offsets[b + 1] as usize;
            let bin = &idx.bucket_candidates[lo..hi];
            bin.contains(&0u32) && bin.contains(&1u32)
        });
        assert!(
            both_here || both_somewhere,
            "wide bin_width must pack both candidate ids into a shared bin"
        );

        // Cursor fill must pack every charge-1 b/y fragment exactly once.
        let expected: usize = cands
            .iter()
            .map(|c| predict_by_ions(&c.peptide, 1..=1).len())
            .sum();
        assert_eq!(
            idx.n_entries(),
            expected,
            "n_entries must equal total charge-1 b/y fragments across both candidates"
        );
    }

    #[test]
    fn unknown_mz_returns_empty() {
        let cands = vec![cand("PEPTIDEK")];
        let idx = FragmentIndex::build(&cands, 0.02);
        assert!(idx.candidates_in_bin(5.0).is_empty());
        assert!(idx.candidates_in_bin(99999.0).is_empty());
    }

    #[test]
    fn voter_ranks_candidate_with_most_matched_fragments_first() {
        // B shares 3 fragments with the observed peaks; A shares 1. B must rank above A.
        let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
        let idx = FragmentIndex::build(&cands, 0.02);
        let b_frags = predict_by_ions(&cands[1].peptide, 1..=1);
        let a_frags = predict_by_ions(&cands[0].peptide, 1..=1);
        // observed peaks (rank, mz): 3 of B's fragments + 1 of A's.
        let peaks = vec![
            (1u32, b_frags[0].mz),
            (2, b_frags[1].mz),
            (3, b_frags[2].mz),
            (4, a_frags[0].mz),
        ];
        let mut voter = FragmentVoter::new(cands.len());
        // in_window = both candidates eligible.
        let topk = voter.top_k(&idx, &peaks, |_cid| true, 2);
        assert_eq!(topk[0], 1u32, "candidate B (3 matches) ranks first");
        assert!(topk.contains(&0u32));
    }

    #[test]
    fn voter_excludes_out_of_window_candidates() {
        let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
        let idx = FragmentIndex::build(&cands, 0.02);
        let b_frags = predict_by_ions(&cands[1].peptide, 1..=1);
        let peaks = vec![(1u32, b_frags[0].mz), (2, b_frags[1].mz)];
        let mut voter = FragmentVoter::new(cands.len());
        // window excludes candidate 1 -> it must not appear even though it has the votes.
        let topk = voter.top_k(&idx, &peaks, |cid| cid == 0, 2);
        assert!(!topk.contains(&1u32));
    }

    #[test]
    fn voter_resets_between_calls() {
        let cands = vec![cand("PEPTIDEK"), cand("ACDEFGHIK")];
        let idx = FragmentIndex::build(&cands, 0.02);
        let b_frags = predict_by_ions(&cands[1].peptide, 1..=1);
        let mut voter = FragmentVoter::new(cands.len());
        let _ = voter.top_k(&idx, &[(1, b_frags[0].mz)], |_| true, 2);
        // second call with NO peaks must yield no votes (scratch cleared).
        let topk = voter.top_k(&idx, &[], |_| true, 2);
        assert!(topk.is_empty());
    }
}
