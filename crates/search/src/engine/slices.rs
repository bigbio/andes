//! Search of a whole run over the memory-mapped candidate index, one precursor-mass slice
//! at a time. Part of the match engine; `PreparedSearch` runs each slice.

use super::*;

pub const INDEX_SLICE_DA: f64 = 80.0;

pub fn match_spectra(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    scorer: &RankScorer,
    fragment_tolerance_da: f64,
    decoy_prefix: &str,
    cache_dir: &std::path::Path,
) -> std::io::Result<(Vec<TopNQueue>, Vec<Candidate>)> {
    let mut params = params.clone();
    params.fragment_index_top_k = crate::params::FRAGMENT_INDEX_TOP_K;
    let cache = crate::candidates::records::index_cache_path(idx, &params, cache_dir);
    let mut prepared = PreparedSearch::prepare_mmap(
        idx,
        &params,
        scorer,
        fragment_tolerance_da,
        decoy_prefix,
        &cache,
    )?;
    let mut order: Vec<usize> = (0..spectra.len()).collect();
    order.sort_by(|&a, &b| {
        crate::fragments::index::index_order_mass(&spectra[a])
            .total_cmp(&crate::fragments::index::index_order_mass(&spectra[b]))
            .then(a.cmp(&b))
    });
    let mut queues: Vec<Option<TopNQueue>> = (0..spectra.len()).map(|_| None).collect();
    let mut rest: &[usize] = &order;
    while !rest.is_empty() {
        let take = crate::fragments::index::index_chunk_len(
            rest,
            |&i| crate::fragments::index::index_order_mass(&spectra[i]),
            INDEX_SLICE_DA,
        );
        let (head, tail) = rest.split_at(take);
        let chunk: Vec<Spectrum> = head.iter().map(|&i| spectra[i].clone()).collect();
        for (mut q, &orig) in prepared.run_chunk(&chunk, 0).into_iter().zip(head) {
            q.update_in_place(|psm| psm.spectrum_idx = orig);
            queues[orig] = Some(q);
        }
        rest = tail;
    }
    prepared.sync_materialized_candidates();
    let queues = queues
        .into_iter()
        .map(|q| q.expect("every spectrum is in exactly one chunk"))
        .collect();
    Ok((queues, prepared.candidates))
}
