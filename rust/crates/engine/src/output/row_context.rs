//! Shared per-PSM row context used by both PIN and TSV writers.
//!
//! Computes spectrum- and PSM-level fields that both formats need (rank,
//! accession string, scan number, spec_id) once per PSM, so each format
//! only has to format columns from a stable struct.

use crate::psm::PsmMatch;
use crate::search_index::SearchIndex;
use crate::spectrum::Spectrum;

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
    /// Build a `RowContext` for one PSM.
    pub(crate) fn new(spec: &Spectrum, psm: &PsmMatch, search_index: &SearchIndex) -> Self {
        let scan = spec.scan.unwrap_or(0);
        let spec_id = if spec.title.is_empty() {
            format!("scan={scan}")
        } else {
            spec.title.clone()
        };
        let accession = resolve_accession(psm, search_index);
        Self { scan, spec_id, accession }
    }
}

/// Resolve a protein accession from the `SearchIndex` for a given PSM.
///
/// The combined target+decoy `ProteinDb` inside `search_index` already carries
/// decoy prefixes on decoy accessions (set by `target_plus_decoy`), so no
/// prefix arithmetic is needed here. Falls back to `"PROT_{idx}"` if the
/// index is out of range.
pub(crate) fn resolve_accession(psm: &PsmMatch, search_index: &SearchIndex) -> String {
    let idx = psm.candidate.protein_index;
    match search_index.protein_at(idx) {
        Some(prot) => prot.accession.clone(),
        None => format!("PROT_{idx}"),
    }
}

/// Iterate a slice of PSMs (pre-sorted best-first) yielding `(rank, psm)`.
///
/// Rank is 1-based and increments only when `spec_e_value` changes — ties
/// share the same rank.
pub(crate) fn iter_ranked(queue_sorted: &[PsmMatch]) -> impl Iterator<Item = (u32, &PsmMatch)> {
    let mut rank = 0u32;
    let mut prev_sev = f64::NAN;
    queue_sorted.iter().map(move |psm| {
        if psm.spec_e_value != prev_sev {
            rank += 1;
            prev_sev = psm.spec_e_value;
        }
        (rank, psm)
    })
}
