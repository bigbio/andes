//! Conservative target-substring rescue of observed entrapment hits, independent of digestion.

use std::collections::HashSet;

use output::PercolatorPsm;
use search::decoy::entrapment::{is_entrapment_accession, target_substrings};
use search::SearchIndex;

pub(super) fn peptides<'a>(
    index: &SearchIndex,
    rows: impl Iterator<Item = &'a PercolatorPsm>,
) -> Result<HashSet<String>, String> {
    let queries: Vec<String> = rows
        .filter(|p| {
            let mut proteins = p.proteins.split(';').filter(|a| !a.is_empty()).peekable();
            proteins.peek().is_some() && proteins.all(is_entrapment_accession)
        })
        .map(|p| output::offsets::stripped_sequence(&p.peptide).replace('I', "L"))
        .filter(|p| !p.is_empty())
        .collect();
    let targets = index
        .iter_target_proteins()
        .filter(|p| !is_entrapment_accession(&p.accession))
        .map(|p| p.sequence.as_slice());
    target_substrings(targets, queries).map_err(|e| {
        format!("--add-entrapment: cannot index observed peptides for target rescue: {e}")
    })
}
