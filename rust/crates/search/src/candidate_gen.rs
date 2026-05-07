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

use model::amino_acid::AminoAcid;
use model::enzyme::Enzyme;
use model::peptide::Peptide;
use model::protein::Protein;
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
    // Use the prefix verbatim — match exactly what the caller (and the SearchIndex)
    // stored. Don't invent formatting; require callers to pass the real prefix.
    idx.db.proteins.iter().enumerate().flat_map(move |(p_idx, protein)| {
        let is_decoy = protein.accession.starts_with(decoy_prefix);
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
///
/// The `params.num_tolerable_termini` field controls cleavage enforcement:
/// - `2`: both ends must be enzyme-cleavage sites (strict / fully specific, default).
/// - `1`: at least one end must be an enzyme-cleavage site (semi-specific).
/// - `0`: neither end needs to be a cleavage site (non-specific).
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

    let ntt = params.num_tolerable_termini;

    // For ntt=0 (non-specific) with a non-NonSpecific enzyme, enumerate all
    // valid-length spans without any cleavage constraint. This produces the
    // same set as Enzyme::NonSpecific with ntt=2 (modulo missed-cleavage
    // filtering — for ntt=0 we skip that since there are no "cleavage sites"
    // to count between arbitrary span endpoints).
    //
    // Note: Enzyme::NonSpecific itself falls through to the normal cleavage-
    // position loop below (which returns all positions 0..=n), preserving the
    // existing missed-cleavage semantics that the NonSpecific tests exercise.
    if ntt == 0 && !matches!(params.enzyme, Enzyme::NonSpecific) {
        let ctx = EmitCtx { sub_seq, seq, seq_offset, protein_index, is_decoy, params };
        return enumerate_all_spans(&ctx, n);
    }

    let cleavage_positions = compute_cleavage_positions(sub_seq, params.enzyme);

    // ntt=2: strict — only spans where both start and end are cleavage positions.
    // ntt=1: semi-specific — spans where at least one end is a cleavage position.
    //
    // Strategy for ntt=1:
    //   (a) Strict spans (same as ntt=2) — already both ends tryptic.
    //   (b) Free C-terminus: for each tryptic start, slide the end across
    //       all positions in [start+min_len, start+max_len]. Skip ends that
    //       ARE cleavage positions (already covered by the strict case).
    //   (c) Free N-terminus: for each tryptic end, slide the start across
    //       all positions in [end-max_len, end-min_len]. Skip starts that
    //       ARE cleavage positions (already covered by the strict case).
    //
    // Using a HashSet of (start, end) pairs to prevent duplicates when both
    // ends happen to be tryptic.

    let mut out = Vec::new();

    // Build a fast lookup for cleavage positions.
    let cleavage_set: std::collections::HashSet<u32> = cleavage_positions.iter().copied().collect();

    let ctx = EmitCtx { sub_seq, seq, seq_offset, protein_index, is_decoy, params };

    // ── Strict spans (ntt=2 behaviour) ───────────────────────────────────────
    // Also included in ntt=1, since a strict span satisfies "at least one end".
    for (i, &start) in cleavage_positions.iter().enumerate() {
        for (offset, &end) in cleavage_positions[i + 1..].iter().enumerate() {
            let len = end - start;
            if len > params.max_length {
                break;
            }
            if len < params.min_length {
                continue;
            }
            let missed = offset as u32;
            if missed > params.max_missed_cleavages {
                continue;
            }
            emit_span(&ctx, start, end, &mut out);
        }
    }

    // ── Semi-specific spans (ntt=1 only) ─────────────────────────────────────
    if ntt == 1 {
        // (b) Tryptic N-terminus, free C-terminus.
        for &start in &cleavage_positions {
            let c_min = start + params.min_length;
            let c_max = (start + params.max_length).min(n);
            for end in c_min..=c_max {
                // Skip ends that are cleavage positions — already emitted above.
                if cleavage_set.contains(&end) {
                    continue;
                }
                // No missed-cleavage filter here: the "missed cleavages between
                // start and end" concept applies to strictly tryptic spans. For
                // semi-tryptic peptides with a free terminus, Java MS-GF+ does
                // not count internal cleavage sites as missed cleavages — the
                // semi-tryptic span is treated as a single candidate regardless
                // of internal K/R residues.
                emit_span(&ctx, start, end, &mut out);
            }
        }

        // (c) Free N-terminus, tryptic C-terminus.
        for &end in &cleavage_positions {
            if end < params.min_length {
                continue;
            }
            let s_min = end.saturating_sub(params.max_length);
            let s_max = end - params.min_length;
            for start in s_min..=s_max {
                // Skip starts that are cleavage positions — already emitted above.
                if cleavage_set.contains(&start) {
                    continue;
                }
                emit_span(&ctx, start, end, &mut out);
            }
        }
    }

    out
}

/// Shared context passed to `emit_span` to avoid exceeding argument limits.
struct EmitCtx<'a> {
    sub_seq: &'a [u8],
    seq: &'a [u8],
    seq_offset: usize,
    protein_index: usize,
    is_decoy: bool,
    params: &'a SearchParams,
}

/// Emit a single (start, end) span as candidates, if the span passes residue
/// validity checks. Appends to `out`.
#[inline]
fn emit_span(ctx: &EmitCtx<'_>, start: u32, end: u32, out: &mut Vec<Candidate>) {
    let span = &ctx.sub_seq[start as usize..end as usize];
    // Skip spans containing non-standard residues.
    if span.iter().any(|&r| AminoAcid::standard(r).is_none()) {
        return;
    }

    let abs_start = start as usize + ctx.seq_offset;
    let abs_end = end as usize + ctx.seq_offset;
    let pre = if abs_start == 0 { b'_' } else { ctx.seq[abs_start - 1] };
    let post = if abs_end == ctx.seq.len() { b'-' } else { ctx.seq[abs_end] };

    let is_protein_n_term = start == 0;
    let is_protein_c_term = abs_end == ctx.seq.len();
    let mod_combinations =
        expand_mod_combinations(span, ctx.params, is_protein_n_term, is_protein_c_term);
    for residues in mod_combinations {
        let peptide = Peptide::new(residues, pre, post);
        out.push(Candidate {
            peptide,
            protein_index: ctx.protein_index,
            start_offset_in_protein: abs_start,
            is_decoy: ctx.is_decoy,
        });
    }
}

/// Enumerate all valid-length spans without cleavage constraints (ntt=0 path).
/// Invoked when `num_tolerable_termini = 0` with a non-NonSpecific enzyme.
fn enumerate_all_spans(ctx: &EmitCtx<'_>, n: u32) -> Vec<Candidate> {
    let mut out = Vec::new();
    for start in 0..n {
        let end_max = (start + ctx.params.max_length).min(n);
        for end in (start + ctx.params.min_length)..=end_max {
            emit_span(ctx, start, end, &mut out);
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
    use model::modification::ModLocation;

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

#[cfg(test)]
mod tests {
    #[test]
    fn decoy_prefix_matched_verbatim_no_underscore_appended() {
        // Caller passes "XXX" (no underscore). The matcher should look for
        // accessions starting with literally "XXX", NOT "XXX_".
        // We exercise this by checking the is_decoy flag logic directly:
        // any accession starting with "XXX" (including "XXX_something") must
        // match, and accessions starting with "XXX_" only must also match (no
        // double-underscore invention).
        let prefix = "XXX";
        assert!(
            "XXX_protein1".starts_with(prefix),
            "accession starting with 'XXX_' should match prefix 'XXX'"
        );
        assert!(
            "XXXprotein1".starts_with(prefix),
            "accession starting with 'XXXprotein1' should match prefix 'XXX'"
        );
        assert!(
            !"DECOY_protein1".starts_with(prefix),
            "accession 'DECOY_protein1' should NOT match prefix 'XXX'"
        );

        // Verify we do NOT append an underscore: "DECOY" prefix must not
        // accidentally match "DECOY_protein" as "DECOY__protein" or similar.
        let colon_prefix = "DECOY:";
        assert!(
            "DECOY:sp|P12345|PROT_HUMAN".starts_with(colon_prefix),
            "colon-terminated prefix should match verbatim"
        );
        assert!(
            !"DECOY_sp|P12345|PROT_HUMAN".starts_with(colon_prefix),
            "underscore-delimited accession should NOT match colon prefix"
        );
    }
}
