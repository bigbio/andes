//! Candidate peptide enumeration via protein-walk. Mirrors Java
//! `CandidatePeptideGrid` but operates per-protein rather than over
//! the suffix array. The candidate set produced is equivalent (only
//! iteration order differs).
//!
//! Phase 4d Task 2 ships the basic enumerator: enzyme-cleaved spans
//! within length range. NO variable-mod expansion (Task 4).
//! Missed-cleavage handling already works correctly via the
//! `missed_count` check; Task 3 adds tests for it.

use crate::amino_acid::AminoAcid;
use crate::enzyme::Enzyme;
use crate::peptide::Peptide;
use crate::protein::Protein;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;

#[derive(Debug, Clone)]
pub struct Candidate {
    pub peptide: Peptide,
    pub protein_index: usize,
    pub start_offset_in_protein: usize,
    pub is_decoy: bool,
}

/// Enumerate every candidate peptide from `idx` matching `params`.
/// Order: by `(protein_index, start_offset, mod_combination_index)`.
pub fn enumerate_candidates<'a>(
    idx: &'a SearchIndex,
    params: &'a SearchParams,
    decoy_prefix: &'a str,
) -> impl Iterator<Item = Candidate> + 'a {
    let normalized_prefix = format!("{}_", decoy_prefix.trim_end_matches('_'));
    idx.db.proteins.iter().enumerate().flat_map(move |(p_idx, protein)| {
        let is_decoy = protein.accession.starts_with(&normalized_prefix);
        enumerate_protein(protein, p_idx, is_decoy, params).into_iter()
    })
}

fn enumerate_protein(
    protein: &Protein,
    protein_index: usize,
    is_decoy: bool,
    params: &SearchParams,
) -> Vec<Candidate> {
    let seq = &protein.sequence;
    let n = seq.len() as u32;
    if n < params.min_length {
        return Vec::new();
    }

    let cleavage_positions = compute_cleavage_positions(seq, params.enzyme);

    let mut out = Vec::new();
    for (i, &start) in cleavage_positions.iter().enumerate() {
        for (offset, &end) in cleavage_positions[i + 1..].iter().enumerate() {
            let len = end - start;
            if len > params.max_length {
                break;  // future values give larger lengths only
            }
            if len < params.min_length {
                continue;
            }
            // Missed cleavages = number of cleavage positions strictly between start and end.
            // offset 0 means adjacent positions → 0 missed; offset k → k missed.
            let missed = offset as u32;
            if missed > params.max_missed_cleavages {
                continue;
            }

            let span = &seq[start as usize..end as usize];
            // Skip spans containing non-standard residues.
            if span.iter().any(|&r| AminoAcid::standard(r).is_none()) {
                continue;
            }

            let pre = if start == 0 { b'_' } else { seq[start as usize - 1] };
            let post = if end as usize == seq.len() { b'-' } else { seq[end as usize] };

            let mod_combinations = expand_mod_combinations(span, params);
            for residues in mod_combinations {
                let peptide = Peptide::new(residues, pre, post);
                out.push(Candidate {
                    peptide,
                    protein_index,
                    start_offset_in_protein: start as usize,
                    is_decoy,
                });
            }
        }
    }

    out
}

/// Generate every combination of variable-mod applications for `span`,
/// up to `params.max_variable_mods_per_peptide` mods total.
fn expand_mod_combinations(span: &[u8], params: &SearchParams) -> Vec<Vec<AminoAcid>> {
    use crate::modification::ModLocation;

    // For each position, the list of variants at that residue.
    let position_variants: Vec<Vec<AminoAcid>> = span.iter().map(|&r| {
        params.aa_set.variants_for(r, ModLocation::Anywhere).to_vec()
    }).collect();

    let mut out = Vec::new();
    let mut current = Vec::with_capacity(span.len());
    expand_recursive(
        &position_variants, 0, &mut current, 0,
        params.max_variable_mods_per_peptide, &mut out,
    );
    out
}

fn expand_recursive(
    position_variants: &[Vec<AminoAcid>],
    pos: usize,
    current: &mut Vec<AminoAcid>,
    mods_used: u32,
    max_mods: u32,
    out: &mut Vec<Vec<AminoAcid>>,
) {
    if pos == position_variants.len() {
        out.push(current.clone());
        return;
    }
    for variant in &position_variants[pos] {
        let new_mods = mods_used + if variant.is_modified() { 1 } else { 0 };
        if new_mods > max_mods {
            continue;
        }
        current.push(variant.clone());
        expand_recursive(
            position_variants, pos + 1, current, new_mods, max_mods, out,
        );
        current.pop();
    }
}

/// Cleavage positions: 0 (start of protein), n (end of protein), and
/// every i in 1..n where `enzyme.is_cleavable_after(seq[i-1])` (for
/// C-term cutters like Trypsin) OR `enzyme.is_cleavable_before(seq[i])`
/// (for N-term cutters like AspN/LysN).
fn compute_cleavage_positions(seq: &[u8], enzyme: Enzyme) -> Vec<u32> {
    let n = seq.len() as u32;

    if matches!(enzyme, Enzyme::NoCleavage) {
        return vec![0, n];
    }

    if matches!(enzyme, Enzyme::NonSpecific) {
        return (0..=n).collect();
    }

    let mut positions = vec![0u32];
    for i in 1..n {
        let prev = seq[(i - 1) as usize];
        let here = seq[i as usize];
        if enzyme.is_cleavable_after(prev) || enzyme.is_cleavable_before(here) {
            positions.push(i);
        }
    }
    if *positions.last().unwrap() != n {
        positions.push(n);
    }
    positions
}
