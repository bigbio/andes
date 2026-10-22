//! Shared per-PSM row context used by both PIN and TSV writers.
//!
//! Computes spectrum- and PSM-level fields that both formats need (rank,
//! accession string, scan number, spec_id) once per PSM, so each format
//! only has to format columns from a stable struct.

use std::borrow::Cow;

use model::spectrum::Spectrum;
use search::candidate_gen::Candidate;
use search::isomer::collapse_isomers;
use search::psm::PsmMatch;
use search::search_index::SearchIndex;

/// Fields derived once per PSM that are used by both PIN and TSV writers.
///
/// Format-specific fields (e.g. PIN's `exp_mass`/`dm`, TSV's `frag_method`)
/// are computed in the per-writer code; only the intersection lives here.
pub(crate) struct RowContext {
    /// Raw scan number (`spec.scan.unwrap_or(0)`).
    pub scan: i32,
    /// Spectrum identifier string: `spec.title` if non-empty, else `"scan=N"`.
    pub spec_id: String,
    /// Resolved protein accession (decoy accessions already carry their prefix).
    pub accession: String,
}

impl RowContext {
    /// Build a `RowContext` for one PSM. Caller passes the resolved
    /// `Candidate` (looked up via `psm.primary_candidate_idx()`) so this layer doesn't
    /// need its own `candidates` slice reference.
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

/// Resolve a protein accession from the `SearchIndex` for a given `Candidate`.
///
/// The combined target+decoy `ProteinDb` inside `search_index` already carries
/// decoy prefixes on decoy accessions (set by `search::decoy::decoy_db`), so no
/// prefix arithmetic is needed here. Falls back to `"PROT_{idx}"` if the
/// index is out of range.
pub(crate) fn resolve_accession(cand: &Candidate, search_index: &SearchIndex) -> String {
    let idx = cand.protein_index;
    match search_index.protein_at(idx) {
        Some(prot) => prot.accession.clone(),
        None => format!("PROT_{idx}"),
    }
}

/// Iterate PSMs (pre-sorted best-first by `rank_score` descending) yielding
/// `(rank, psm)`. Rank is 1-based and increments when `rank_score` changes —
/// ties share the same rank. `rank_score` (RawScore) is the sole ranking
/// signal now that the generating function is removed.
pub(crate) fn iter_ranked_by_rank_score<'a, I>(
    queue_sorted: I,
) -> impl Iterator<Item = (u32, &'a PsmMatch)>
where
    I: IntoIterator<Item = &'a PsmMatch>,
{
    let mut rank = 0u32;
    let mut prev_rs: Option<f32> = None;
    queue_sorted.into_iter().map(move |psm| {
        // Two NaN scores tie (NaNs collapse into one worst-score bucket in the
        // training/gate paths), so a run of NaN-scored PSMs must share one rank
        // rather than each starting a new one as `!= f32::NAN` would do.
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

/// The rows one spectrum emits, as `(rank, psm)` in emission order.
///
/// `psms_sorted` is the spectrum's queue sorted best-first
/// (`TopNQueue::into_rank_sorted_vec`). Positional isomers are collapsed
/// ([`search::isomer::collapse_isomers`]): of the PSMs sharing an
/// [`search::isomer::IsomerKey`] (same precursor, charge, target/decoy side,
/// residues and modification-mass multiset), only the first in rank order is
/// emitted. Site alternatives of one peptide score almost identically and would
/// otherwise crowd the per-spectrum competition downstream (Percolator keeps the
/// best PSM of each spectrum), so one row carries the peptide. Chimeric
/// secondaries are matched at their own precursor and never collapse into the
/// primary. Ranks are assigned over the emitted rows (dense, 1-based, ties
/// shared), so they stay contiguous after the collapse.
///
/// A collapsed PSM can hold the same residues in a protein the kept one does not
/// list (the same peptide with different flanking residues scores differently,
/// so the per-charge dedup keeps the copies apart). The kept row takes the union
/// of their proteins, one candidate per protein, so its `Proteins` column still
/// names every protein the retained copies came from. A row with no isomer is
/// borrowed unchanged.
///
/// Every writer that emits per-PSM rows (PIN, TSV, QPX) goes through this, so
/// the same rows — and the same `SpecId`s — appear in each.
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

/// `kept` with the candidates of the PSMs collapsed into it appended, one per
/// protein not already listed; borrowed when nothing is added.
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
