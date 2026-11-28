//! Building a `PreparedSearch`: over a memory-mapped or in-memory index, with the open-search
//! index when it is on.

use std::sync::atomic::{AtomicU64, AtomicUsize};
use std::sync::{Mutex, RwLock};

use model::amino_acid::set::AminoAcidSet;
use model::enzyme::Enzyme;
use scoring_crate::psm::RankScorer;

use super::dedup::*;
use super::{MmapAccumulator, PreparedSearch};
use crate::candidates::records::MmapCandidateIndex;
use crate::fragments::index::ChunkFragmentIndex;
use crate::index::SearchIndex;
use crate::params::SearchParams;

impl<'a> PreparedSearch<'a> {
    pub fn prepare_mmap(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        decoy_prefix: &str,
        index_path: &std::path::Path,
    ) -> std::io::Result<Self> {
        let (mmap_index, was_built) =
            MmapCandidateIndex::open_or_build(index_path, idx, params, decoy_prefix)?;
        if mmap_index.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "the digest built no candidate peptides (enzyme {}, extra enzymes {}, \
                     length {}..{}, max missed cleavages {}, termini {}); check --enzyme, \
                     --min-length, --max-length, --max-missed-cleavages and the database",
                    params.enzyme.name(),
                    params.extra_enzymes.len(),
                    params.min_length,
                    params.max_length,
                    params.max_missed_cleavages,
                    params.num_tolerable_termini,
                ),
            ));
        }
        if was_built {
            log::info!(
                target: "candidate-index",
                "built and cached {} records → {}",
                mmap_index.len(),
                index_path.display()
            );
        } else {
            log::info!(
                target: "candidate-index",
                "reusing cached index ({} records) from {}",
                mmap_index.len(),
                index_path.display()
            );
        }
        let mut prepared = Self::over_index(idx, params, scorer, fragment_tolerance_da, mmap_index);
        prepared.open_index = Self::build_open_index(idx, params, scorer, &prepared.mmap_index);
        Ok(prepared)
    }

    pub fn prepare_in_memory(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        decoy_prefix: &str,
    ) -> std::io::Result<Self> {
        let mmap_index = MmapCandidateIndex::build_in_memory(idx, params, decoy_prefix)?;
        let mut prepared = Self::over_index(idx, params, scorer, fragment_tolerance_da, mmap_index);
        prepared.open_index = Self::build_open_index(idx, params, scorer, &prepared.mmap_index);
        Ok(prepared)
    }

    fn build_open_index(
        idx: &SearchIndex,
        params: &SearchParams,
        scorer: &RankScorer,
        mmap_index: &MmapCandidateIndex,
    ) -> Option<Box<ChunkFragmentIndex>> {
        params.open_search?;
        let t = std::time::Instant::now();
        let records = mmap_index.records();
        let (lo, hi) = match (records.first(), records.last()) {
            (Some(a), Some(b)) => (a.mass_milli as f64 / 1000.0, b.mass_milli as f64 / 1000.0),
            _ => return None,
        };
        let max_mod = params
            .aa_set
            .distinct_mods()
            .iter()
            .filter(|m| !m.fixed)
            .map(|m| m.mass_delta.max(0.0))
            .fold(0.0_f64, f64::max);
        let max_sites = params.max_variable_mods_per_peptide + params.max_label_sites;
        let hi = hi + max_mod * f64::from(max_sites) + 1.0;
        let index = ChunkFragmentIndex::build(
            records,
            idx,
            params,
            scorer.feature_match_tolerance(),
            lo,
            hi,
        );
        log::info!(
            target: "open-search-index",
            "{:.1}-{:.1} Da, {} forms, built in {:.1} s",
            lo,
            hi,
            index.n_forms(),
            t.elapsed().as_secs_f64()
        );
        Some(Box::new(index))
    }

    fn over_index(
        idx: &'a SearchIndex,
        params: &'a SearchParams,
        scorer: &'a RankScorer,
        fragment_tolerance_da: f64,
        mmap_index: MmapCandidateIndex,
    ) -> Self {
        let mut aa_set_for_scoring: AminoAcidSet = params.aa_set.clone();
        if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
            aa_set_for_scoring.register_enzyme(params.enzyme, 0.95, 0.95);
        }
        let target_bare_seqs = build_target_bare_seqs(&mmap_index, idx);
        PreparedSearch {
            idx,
            params,
            scorer,
            fragment_tolerance_da,
            candidates: Vec::new(),
            mmap_index,
            mmap_accum: Mutex::new(MmapAccumulator::default()),
            chunk_prep_nanos: AtomicU64::new(0),
            open_rescore_attempts: AtomicUsize::new(0),
            open_rescore_wins: AtomicUsize::new(0),
            aa_set_for_scoring,
            ms1_link: None,
            dia_ms1: Vec::new(),
            dia_ms1_sorted: Vec::new(),
            dia_scan_mobility: input::DiaScanMobility::default(),
            dia_mobility_margin: std::sync::RwLock::new(None),
            target_bare_seqs,
            candidate_filter: None,
            open_index: None,
            per_scan: None,
            form_rt_model: None,
            dia_rt_prior: RwLock::new(None),
            form_rt_cache: crate::dia::FormRtCache::default(),
            frag_pred_cache: crate::fragments::predictions::FragPredCache::default(),
        }
    }
}
