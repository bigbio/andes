//! Out-of-core search over the memory-mapped candidate index, and the second (co-isolation)
//! pass that rescores a spectrum's precursor window for co-isolated peptides. Part of the
//! match engine; `PreparedSearch` drives both.

use super::*;

pub const OUT_OF_CORE_SLICE_DA: f64 = 80.0;

pub fn match_spectra_out_of_core(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    scorer: &RankScorer,
    fragment_tolerance_da: f64,
    decoy_prefix: &str,
    cache_dir: &std::path::Path,
) -> std::io::Result<(Vec<TopNQueue>, Vec<Candidate>)> {
    let mut params = params.clone();
    params.fragment_index_top_k = crate::search_params::FRAGMENT_INDEX_TOP_K;
    let cache = crate::candidate_index::index_cache_path(idx, &params, cache_dir);
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
        crate::fragment_index::index_order_mass(&spectra[a])
            .total_cmp(&crate::fragment_index::index_order_mass(&spectra[b]))
            .then(a.cmp(&b))
    });
    let mut queues: Vec<Option<TopNQueue>> = (0..spectra.len()).map(|_| None).collect();
    let mut rest: &[usize] = &order;
    while !rest.is_empty() {
        let take = crate::fragment_index::index_chunk_len(
            rest,
            |&i| crate::fragment_index::index_order_mass(&spectra[i]),
            OUT_OF_CORE_SLICE_DA,
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

/// The isolation window around its target m/z (the precursor m/z when the scan records no
/// target), and the precursor tolerance in Da.
pub(super) fn coisolation_window(spec: &Spectrum, params: &SearchParams) -> (f64, f64, f64) {
    let (lo, hi) = crate::precursor_matching::isolation_window_mz(
        spec,
        params.chimeric_isolation_halfwidth_da,
    );
    let tol = params
        .precursor_tolerance
        .left
        .as_da(spec.precursor_mz)
        .max(0.01);
    (lo, hi, tol)
}

pub fn coisolation_ms1_bounds(spec: &Spectrum, params: &SearchParams) -> (f64, f64) {
    const MARGIN_DA: f64 = 0.01;
    let (lo, hi, tol) = coisolation_window(spec, params);
    (
        lo - tol - MARGIN_DA,
        hi + 3.0 * model::mass::ISOTOPE + tol + MARGIN_DA,
    )
}

pub(super) struct Pass2Scan {
    primary: Peptide,
    coisolated: Vec<crate::coisolation::CoIsolated>,
}

pub fn run_pass2_coisolation(
    prepared: &PreparedSearch,
    spectra: &[Spectrum],
    queues: &mut [TopNQueue],
    params: &SearchParams,
    link: &Ms1Link,
    offset: usize,
) {
    if !params.chimeric {
        return;
    }
    let shifted;
    let spectra = match recentred_spectra(spectra, prepared.scorer.fragment_recentre_ppm()) {
        Some(v) => {
            shifted = v;
            &shifted[..]
        }
        None => spectra,
    };
    let scorer_set = prepared.scorer_set(spectra);

    let scans: Vec<Option<Pass2Scan>> = queues
        .par_iter()
        .enumerate()
        .map(|(spec_idx, q)| {
            if q.is_empty() {
                return None;
            }
            let spec = spectra.get(spec_idx)?;
            let ms1_idx = (*link.ms2_to_ms1.get(spec_idx)?)?;
            let ms1 = link.ms1_peaks.get(ms1_idx)?;
            let (lo, hi, tol) = coisolation_window(spec, params);
            let coisolated = crate::coisolation::detect_coisolated(
                ms1,
                lo,
                hi,
                spec.precursor_mz,
                *params.charge_range.start()..=*params.charge_range.end(),
                tol,
                params.chimeric_max_kl,
                params.chimeric_max_coisolated,
            );
            if coisolated.is_empty() {
                return None;
            }
            let primary = prepared.candidate_peptide(q.peek_top()?.primary_candidate_idx());
            Some(Pass2Scan {
                primary,
                coisolated,
            })
        })
        .collect();

    let found = mmap_secondaries(prepared, &scorer_set, spectra, &scans, params);

    let mut accum = prepared.mmap_accum.lock().unwrap();
    for (spec_idx, (q, secondaries)) in queues.iter_mut().zip(found).enumerate() {
        if secondaries.is_empty() {
            continue;
        }
        let mut emitted: std::collections::BTreeSet<(PepDedupKey, u8, bool)> = q
            .iter_psms()
            .filter_map(|p| {
                let i = p.primary_candidate_idx() as usize;
                let c = accum
                    .candidates
                    .get(i)
                    .or_else(|| prepared.candidates.get(i))?;
                Some((
                    PepDedupKey::from_peptide(&c.peptide),
                    p.charge_used,
                    c.is_decoy,
                ))
            })
            .collect();
        for (mut psm, cand) in secondaries {
            let key = (
                PepDedupKey::from_peptide(&cand.peptide),
                psm.charge_used,
                cand.is_decoy,
            );
            if !emitted.insert(key) {
                continue;
            }
            psm.candidate_idxs = smallvec::smallvec![accum.intern(&cand)];
            psm.spectrum_idx = offset + spec_idx;
            q.force_push(psm);
        }
    }
}

pub(super) const PASS2_POOL_SPAN_NOMINAL: i32 = 32;

pub(super) fn mmap_secondaries(
    prepared: &PreparedSearch,
    scorer_set: &ScorerSet<'_>,
    spectra: &[Spectrum],
    scans: &[Option<Pass2Scan>],
    params: &SearchParams,
) -> Vec<Vec<(PsmMatch, Candidate)>> {
    use crate::coisolation::{search_secondary, secondary_nominal_window, CoIsolated};
    type Hit = (PsmMatch, Candidate, std::collections::HashSet<i64>);
    type Item = (usize, CoIsolated, (i32, i32));

    let mut claimed: Vec<std::collections::HashSet<i64>> =
        vec![std::collections::HashSet::new(); scans.len()];
    let mut found: Vec<Vec<(PsmMatch, Candidate)>> = vec![Vec::new(); scans.len()];
    let rounds = scans
        .iter()
        .flatten()
        .map(|s| s.coisolated.len())
        .max()
        .unwrap_or(0);
    for k in 0..rounds {
        let mut items: Vec<Item> = scans
            .iter()
            .enumerate()
            .filter_map(|(i, s)| {
                let co = *s.as_ref()?.coisolated.get(k)?;
                Some((i, co, secondary_nominal_window(co, params)))
            })
            .collect();
        items.sort_by_key(|&(i, _, w)| (w, i));
        let mut batches: Vec<&[Item]> = Vec::new();
        let mut rest = items.as_slice();
        while let Some(&(_, _, (lo, _))) = rest.first() {
            let n = rest
                .iter()
                .take_while(|&&(_, _, (_, hi))| hi - lo <= PASS2_POOL_SPAN_NOMINAL)
                .count()
                .max(1);
            let (batch, tail) = rest.split_at(n);
            batches.push(batch);
            rest = tail;
        }
        let claimed_ref = &claimed;
        let hits: Vec<(usize, Option<Hit>)> = batches
            .par_iter()
            .flat_map_iter(|batch| {
                let mut windows: Vec<(i32, i32)> = Vec::new();
                for &(_, _, (lo, hi)) in batch.iter() {
                    match windows.last_mut() {
                        Some(last) if lo <= last.1 + 1 => last.1 = last.1.max(hi),
                        _ => windows.push((lo, hi)),
                    }
                }
                let pool = prepared.expand_mmap_window_candidates(prepared.params, &windows, None);
                let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
                for (i, c) in pool.iter().enumerate() {
                    bucket_index
                        .entry(c.peptide.nominal_residue_mass())
                        .or_default()
                        .push(i);
                }
                batch
                    .iter()
                    .map(|&(spec_idx, co, _)| {
                        let scan = scans[spec_idx]
                            .as_ref()
                            .expect("only scans with co-isolated precursors are batched");
                        let hit = search_secondary(
                            &spectra[spec_idx],
                            &scan.primary,
                            &claimed_ref[spec_idx],
                            co,
                            &pool,
                            &bucket_index,
                            scorer_set.for_spec(&spectra[spec_idx]),
                            &prepared.aa_set_for_scoring,
                            params,
                            prepared.fragment_tolerance_da,
                        )
                        .map(|(psm, winner_claimed)| {
                            let cand = pool[psm.primary_candidate_idx() as usize].clone();
                            (psm, cand, winner_claimed)
                        });
                        (spec_idx, hit)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        for (spec_idx, hit) in hits {
            if let Some((psm, cand, winner_claimed)) = hit {
                claimed[spec_idx].extend(winner_claimed);
                found[spec_idx].push((psm, cand));
            }
        }
    }
    found
}
