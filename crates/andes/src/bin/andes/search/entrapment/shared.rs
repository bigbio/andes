//! Conservative target-substring rescue of observed entrapment hits, independent of digestion.

use std::collections::{BTreeSet, HashSet};

use output::PercolatorPsm;
use search::decoy::entrapment::is_entrapment_accession;
use search::SearchIndex;

pub(super) fn peptides<'a>(
    index: &SearchIndex,
    rows: impl Iterator<Item = &'a PercolatorPsm>,
) -> Result<HashSet<String>, String> {
    let queries: BTreeSet<String> = rows
        .filter(|p| {
            let mut proteins = p.proteins.split(';').filter(|a| !a.is_empty()).peekable();
            proteins.peek().is_some() && proteins.all(is_entrapment_accession)
        })
        .map(|p| output::offsets::stripped_sequence(&p.peptide).replace('I', "L"))
        .filter(|p| !p.is_empty())
        .collect();
    if queries.is_empty() {
        return Ok(HashSet::new());
    }
    let queries: Vec<String> = queries.into_iter().collect();
    let matcher = aho_corasick::AhoCorasick::new(&queries).map_err(|e| {
        format!("--add-entrapment: cannot index observed peptides for target rescue: {e}")
    })?;
    let mut found = vec![false; queries.len()];
    let mut sequence = Vec::new();
    for protein in index
        .iter_target_proteins()
        .filter(|p| !is_entrapment_accession(&p.accession))
    {
        sequence.clear();
        sequence.extend(
            protein
                .sequence
                .iter()
                .map(|&r| if r == b'I' { b'L' } else { r }),
        );
        for hit in matcher.find_overlapping_iter(&sequence) {
            found[hit.pattern().as_usize()] = true;
        }
    }
    Ok(queries
        .into_iter()
        .zip(found)
        .filter_map(|(p, found)| found.then_some(p))
        .collect())
}
