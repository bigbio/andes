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

            let is_protein_n_term = start == 0;
            let is_protein_c_term = end as usize == seq.len();
            let mod_combinations =
                expand_mod_combinations(span, params, is_protein_n_term, is_protein_c_term);
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
///
/// `is_protein_n_term`: the span begins at position 0 of the protein sequence.
/// `is_protein_c_term`: the span ends at the last residue of the protein sequence.
///
/// These flags control which terminal-location mod variants are consulted:
/// - Position 0: Protein_N_Term (if is_protein_n_term) or N_Term variants are
///   merged in addition to Anywhere variants.
/// - Position n-1: Protein_C_Term (if is_protein_c_term) or C_Term variants are
///   merged in addition to Anywhere variants.
/// - All other positions: Anywhere only (unchanged).
///
/// Mirrors Java `CandidatePeptideGrid.java:43` which maintains separate cached
/// AA-set arrays per terminal context (aaSetN, aaSetC, aaSetProtN, aaSetProtC).
fn expand_mod_combinations(
    span: &[u8],
    params: &SearchParams,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
) -> Vec<Vec<AminoAcid>> {
    use crate::modification::ModLocation;

    let n = span.len();
    // For each position, the list of variants at that residue.
    let position_variants: Vec<Vec<AminoAcid>> = span.iter().enumerate().map(|(i, &r)| {
        let mut variants = params.aa_set.variants_for(r, ModLocation::Anywhere).to_vec();

        // Position 0: add N-Term or Protein_N_Term variants.
        if i == 0 {
            let term_loc = if is_protein_n_term {
                ModLocation::ProtNTerm
            } else {
                ModLocation::NTerm
            };
            for v in params.aa_set.variants_for(r, term_loc) {
                if !variants.contains(v) {
                    variants.push(v.clone());
                }
            }
        }

        // Position n-1: add C-Term or Protein_C_Term variants.
        // (A single-residue span has i==0 and i==n-1 simultaneously — both sets apply.)
        if i == n - 1 {
            let term_loc = if is_protein_c_term {
                ModLocation::ProtCTerm
            } else {
                ModLocation::CTerm
            };
            for v in params.aa_set.variants_for(r, term_loc) {
                if !variants.contains(v) {
                    variants.push(v.clone());
                }
            }
        }

        variants
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
