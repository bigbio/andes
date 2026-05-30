//! Chimeric Sage-style fragment index (Approach B). Peptides sorted by precursor
//! mass (so a mass window is a contiguous index range); fragments sorted by m/z in
//! fixed buckets, each bucket re-sorted by peptide index. The per-spectrum query
//! (Task 2) does a dual binary search bounded to the precursor window — the bound
//! Approach A's vote-all-touched prefilter lacked.
//!
//! Fields below (`sorted_cand`, `sorted_mass`, `fragments`, `bucket_min_mz`, and
//! `Frag.mz`) are populated by `build` here but only read by the per-spectrum
//! `query` added in Task 2; allow dead_code until that lands.
#![allow(dead_code)]

use crate::candidate_gen::Candidate;
use scoring_crate::scoring::fragment_ions::predict_by_ions;

/// Fixed fragment bucket size (Sage uses 8192; power-of-two).
const BUCKET: usize = 8192;

#[derive(Clone, Copy)]
struct Frag {
    mz: f32,
    /// index into `sorted_cand` / `sorted_mass` (Sage's PeptideIx).
    pidx: u32,
}

pub(crate) struct SageIndex {
    /// candidate ids ordered by ascending peptide neutral mass.
    sorted_cand: Vec<u32>,
    /// parallel to `sorted_cand`: ascending neutral masses (for precursor binary search).
    sorted_mass: Vec<f64>,
    /// all candidates' charge-1 b/y fragments; globally m/z-sorted, then each
    /// `BUCKET`-sized chunk re-sorted by `pidx`.
    fragments: Vec<Frag>,
    /// min fragment m/z per bucket (len = ceil(fragments.len()/BUCKET)).
    bucket_min_mz: Vec<f32>,
}

impl SageIndex {
    /// Build over the full candidate set (target+decoy, mod-expanded).
    pub(crate) fn build(candidates: &[Candidate]) -> Self {
        // 1. mass-sorted candidate order.
        let mut order: Vec<u32> = (0..candidates.len() as u32).collect();
        order.sort_by(|&a, &b| {
            candidates[a as usize]
                .peptide
                .mass()
                .partial_cmp(&candidates[b as usize].peptide.mass())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let sorted_mass: Vec<f64> = order
            .iter()
            .map(|&c| candidates[c as usize].peptide.mass())
            .collect();

        // 2. fragments at each candidate's pidx (compute predict_by_ions ONCE/candidate).
        let mut fragments: Vec<Frag> = Vec::new();
        for (pidx, &cid) in order.iter().enumerate() {
            for ion in predict_by_ions(&candidates[cid as usize].peptide, 1..=1) {
                fragments.push(Frag {
                    mz: ion.mz as f32,
                    pidx: pidx as u32,
                });
            }
        }

        // 3. global m/z sort, then per-bucket re-sort by pidx.
        fragments.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));
        let mut bucket_min_mz = Vec::with_capacity(fragments.len() / BUCKET + 1);
        let mut start = 0;
        while start < fragments.len() {
            let end = (start + BUCKET).min(fragments.len());
            bucket_min_mz.push(fragments[start].mz);
            fragments[start..end].sort_by_key(|f| f.pidx);
            start = end;
        }

        SageIndex {
            sorted_cand: order,
            sorted_mass,
            fragments,
            bucket_min_mz,
        }
    }

    /// Total indexed fragment entries (memory accounting: 8 B each).
    pub(crate) fn n_fragments(&self) -> usize {
        self.fragments.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;

    fn cand(seq: &str) -> Candidate {
        let residues = seq.bytes().map(|b| AminoAcid::standard(b).unwrap()).collect();
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
    fn build_sorts_by_mass_and_indexes_fragments() {
        // ACDEFGHIK is heavier than AAK -> mass order must be [AAK, ACDEFGHIK].
        let cands = vec![cand("ACDEFGHIK"), cand("AAK")];
        let idx = SageIndex::build(&cands);
        // sorted_mass ascending:
        assert!(idx.sorted_mass[0] <= idx.sorted_mass[1]);
        // lighter peptide (AAK = candidate id 1) is pidx 0:
        assert_eq!(idx.sorted_cand[0], 1u32);
        // total fragments = sum of charge-1 b/y across both candidates:
        let expect: usize = cands
            .iter()
            .map(|c| predict_by_ions(&c.peptide, 1..=1).len())
            .sum();
        assert_eq!(idx.n_fragments(), expect);
        // bucket_min_mz is non-decreasing.
        for w in idx.bucket_min_mz.windows(2) {
            assert!(w[0] <= w[1]);
        }
    }
}
