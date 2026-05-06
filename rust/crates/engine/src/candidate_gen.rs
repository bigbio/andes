//! Candidate peptide enumeration via protein-walk. Mirrors Java
//! `CandidatePeptideGrid` but operates per-protein rather than over
//! the suffix array. The candidate set produced is equivalent (only
//! iteration order differs).
//!
//! Phase 4d Task 2 ships the basic enumerator: enzyme-cleaved spans
//! within length range. NO variable-mod expansion (Task 4).
//! Missed-cleavage handling already works correctly via the
//! `missed_count` check; Task 3 adds tests for it.
//!
//! Track B5: N-terminal Met cleavage support.
//! Mirrors Java `CandidatePeptideGridConsideringMetCleavage.java:16`.
//! When a protein starts with M, a parallel enumeration treats
//! `sequence[1..]` as the effective protein sequence (initial Met loss).
//! Both enumerations run concurrently; Met-cleaved candidates differ by
//! `is_protein_n_term=true` at offset 1 of the original sequence and are
//! NOT deduplicated — they have a distinct search space (protein-N-term
//! mod variants apply) matching Java's behaviour.

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

    // Standard enumeration: full sequence from offset 0.
    let mut out = enumerate_protein_from_offset(seq, 0, protein_index, is_decoy, params);

    // N-terminal Met cleavage: mirrors Java CandidatePeptideGridConsideringMetCleavage.java:16.
    // When the protein starts with M (and has >1 residue), also enumerate candidates
    // treating sequence[1..] as the effective start. The Met-cleaved peptides still carry
    // is_protein_n_term=true (the post-Met residue is the new biological N-terminus) and
    // are NOT deduplicated — they differ by terminal-mod search space.
    if seq.first() == Some(&b'M') && seq.len() > 1 {
        out.extend(enumerate_protein_from_offset(seq, 1, protein_index, is_decoy, params));
    }

    out
}

/// Enumerate candidates starting from `seq_offset` into `seq`.
///
/// `seq_offset = 0` → normal full-protein walk.
/// `seq_offset = 1` → Met-cleaved walk: `seq[1..]` is the effective protein
///   sequence. Cleavage positions, lengths, and missed-cleavage counts are
///   computed over the sub-sequence. The `start_offset_in_protein` stored on
///   each `Candidate` is adjusted back to the original protein coordinates
///   (i.e. `sub_start + seq_offset`). When `sub_start == 0`, `is_protein_n_term`
///   is set to `true` — the post-Met residue is the effective protein N-terminus.
///   The `pre` context residue for sub_start == 0 is `b'M'` (the cleaved Met).
fn enumerate_protein_from_offset(
    seq: &[u8],
    seq_offset: usize,
    protein_index: usize,
    is_decoy: bool,
    params: &SearchParams,
) -> Vec<Candidate> {
    let sub_seq = &seq[seq_offset..];
    let n = sub_seq.len() as u32;
    if n < params.min_length {
        return Vec::new();
    }

    let cleavage_positions = compute_cleavage_positions(sub_seq, params.enzyme);

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

            let span = &sub_seq[start as usize..end as usize];
            // Skip spans containing non-standard residues.
            if span.iter().any(|&r| AminoAcid::standard(r).is_none()) {
                continue;
            }

            // Context residues in original sequence coordinates.
            // For the Met-cleaved pass (seq_offset == 1) and sub_start == 0, the
            // context residue to the left is the cleaved M (seq[0]).
            let abs_start = start as usize + seq_offset;
            let abs_end = end as usize + seq_offset;
            let pre = if abs_start == 0 {
                b'_'
            } else {
                seq[abs_start - 1]
            };
            let post = if abs_end == seq.len() { b'-' } else { seq[abs_end] };

            // is_protein_n_term: true when the span begins at the effective protein start.
            // For seq_offset == 0: start == 0 (same as before).
            // For seq_offset == 1 (Met-cleaved): start == 0 in sub_seq means the
            //   post-Met residue is at the biological N-terminus — still protein-N-term.
            let is_protein_n_term = start == 0;
            let is_protein_c_term = abs_end == seq.len();
            let mod_combinations =
                expand_mod_combinations(span, params, is_protein_n_term, is_protein_c_term);
            for residues in mod_combinations {
                let peptide = Peptide::new(residues, pre, post);
                out.push(Candidate {
                    peptide,
                    protein_index,
                    start_offset_in_protein: abs_start,
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
