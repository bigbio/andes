//! Shared per-PSM row context used by both PIN and TSV writers.
//!
//! Computes spectrum- and PSM-level fields that both formats need (rank,
//! accession string, scan number, spec_id) once per PSM, so each format
//! only has to format columns from a stable struct.

use search::candidate_gen::Candidate;
use search::psm::PsmMatch;
use search::search_index::SearchIndex;
use model::spectrum::Spectrum;

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
        Self { scan, spec_id, accession }
    }
}

/// Resolve a protein accession from the `SearchIndex` for a given `Candidate`.
///
/// The combined target+decoy `ProteinDb` inside `search_index` already carries
/// decoy prefixes on decoy accessions (set by `target_plus_decoy`), so no
/// prefix arithmetic is needed here. Falls back to `"PROT_{idx}"` if the
/// index is out of range.
pub(crate) fn resolve_accession(cand: &Candidate, search_index: &SearchIndex) -> String {
    let idx = cand.protein_index;
    match search_index.protein_at(idx) {
        Some(prot) => prot.accession.clone(),
        None => format!("PROT_{idx}"),
    }
}

/// Iterate a slice of PSMs (pre-sorted best-first by `rank_score` descending)
/// yielding `(rank, psm)`. Rank is 1-based and increments when `rank_score`
/// changes — ties share the same rank. `rank_score` (RawScore) is the sole
/// ranking signal now that the generating function is removed.
pub(crate) fn iter_ranked_by_rank_score(
    queue_sorted: &[PsmMatch],
) -> impl Iterator<Item = (u32, &PsmMatch)> {
    let mut rank = 0u32;
    let mut prev_rs: Option<f32> = None;
    queue_sorted.iter().map(move |psm| {
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
