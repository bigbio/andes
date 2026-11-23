use std::borrow::Cow;

use model::spectrum::Spectrum;
use search::candidates::Candidate;
use search::isomer::collapse_isomers;
use search::psm::PsmMatch;
use search::index::SearchIndex;

pub(crate) struct RowContext {
    pub scan: i32,
    pub spec_id: String,
    pub accession: String,
}

impl RowContext {
    pub(crate) fn new(spec: &Spectrum, cand: &Candidate, search_index: &SearchIndex) -> Self {
        let scan = spec.scan.unwrap_or(0);
        let spec_id = if spec.title.is_empty() {
            format!("scan={scan}")
        } else {
            spec.title.clone()
        };
        let accession = resolve_accession(cand, search_index);
        Self {
            scan,
            spec_id,
            accession,
        }
    }
}

pub(crate) fn resolve_accession(cand: &Candidate, search_index: &SearchIndex) -> String {
    let idx = cand.protein_index;
    match search_index.protein_at(idx) {
        Some(prot) => prot.accession.clone(),
        None => format!("PROT_{idx}"),
    }
}

pub(crate) fn iter_ranked_by_rank_score<'a, I>(
    queue_sorted: I,
) -> impl Iterator<Item = (u32, &'a PsmMatch)>
where
    I: IntoIterator<Item = &'a PsmMatch>,
{
    let mut rank = 0u32;
    let mut prev_rs: Option<f32> = None;
    queue_sorted.into_iter().map(move |psm| {
        let ties_prev = prev_rs.is_some_and(|prev| {
            psm.rank_score == prev || (psm.rank_score.is_nan() && prev.is_nan())
        });
        if !ties_prev {
            rank += 1;
            prev_rs = Some(psm.rank_score);
        }
        (rank, psm)
    })
}

pub fn emitted_rows<'a>(
    spec: &Spectrum,
    psms_sorted: &'a [PsmMatch],
    candidates: &[Candidate],
) -> Vec<(u32, Cow<'a, PsmMatch>)> {
    let kept: Vec<Cow<'a, PsmMatch>> = if psms_sorted.len() <= 1 {
        psms_sorted.iter().map(Cow::Borrowed).collect()
    } else {
        collapse_isomers(spec.precursor_mz, psms_sorted, candidates)
            .into_iter()
            .map(|(i, collapsed)| {
                with_collapsed_proteins(&psms_sorted[i], &collapsed, psms_sorted, candidates)
            })
            .collect()
    };
    let ranks: Vec<u32> = iter_ranked_by_rank_score(kept.iter().map(|p| &**p))
        .map(|(rank, _)| rank)
        .collect();
    ranks.into_iter().zip(kept).collect()
}

fn with_collapsed_proteins<'a>(
    kept: &'a PsmMatch,
    collapsed: &[usize],
    psms_sorted: &[PsmMatch],
    candidates: &[Candidate],
) -> Cow<'a, PsmMatch> {
    let mut proteins: Vec<usize> = kept
        .candidate_idxs
        .iter()
        .map(|&c| candidates[c as usize].protein_index)
        .collect();
    let mut extra: Vec<u32> = Vec::new();
    for &j in collapsed {
        for &c in &psms_sorted[j].candidate_idxs {
            let protein = candidates[c as usize].protein_index;
            if !proteins.contains(&protein) {
                proteins.push(protein);
                extra.push(c);
            }
        }
    }
    if extra.is_empty() {
        Cow::Borrowed(kept)
    } else {
        let mut merged = kept.clone();
        merged.candidate_idxs.extend(extra);
        Cow::Owned(merged)
    }
}
