//! Candidate peptide enumeration via per-protein walk.
//!
//! Enumerates enzyme-cleaved spans within the configured length range,
//! including missed-cleavage spans (governed by the `missed_count` check).
//!
//! ## N-terminal Met cleavage
//!
//! When a protein starts with M, a parallel enumeration treats
//! `sequence[1..]` as the effective protein sequence (initial Met loss).
//! Both enumerations run concurrently; Met-cleaved candidates differ by
//! `is_protein_n_term=true` at offset 1 of the original sequence and are
//! NOT deduplicated — they have a distinct search space (protein-N-term
//! mod variants apply).

use model::amino_acid::AminoAcid;
use model::enzyme::Enzyme;
use model::modification::ModLocation;
use model::peptide::Peptide;
use model::protein::Protein;
use crate::decoy::is_decoy_accession_affix;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;

#[derive(Debug, Clone)]
pub struct Candidate {
    pub peptide: Peptide,
    pub protein_index: usize,
    pub start_offset_in_protein: usize,
    pub is_decoy: bool,
    /// True when this peptide spans the protein's biological N-terminus.
    /// For Met-cleaved peptides, this is true even though `start_offset_in_protein > 0`.
    pub is_protein_n_term: bool,
    /// True when this peptide spans the protein's C-terminus
    /// (`abs_end == sequence_length`).
    pub is_protein_c_term: bool,
}

/// Enumerate every candidate peptide from `idx` matching `params`.
/// Order: by `(protein_index, start_offset, mod_combination_index)`.
pub fn enumerate_candidates<'a>(
    idx: &'a SearchIndex,
    params: &'a SearchParams,
    decoy_prefix: &'a str,
) -> impl Iterator<Item = Candidate> + 'a {
    // Decoy membership uses the SHARED affix test (prefix needle "<prefix>_" OR
    // the index's optional decoy suffix) — the SAME test SearchIndex uses — so
    // candidate generation and the index can never diverge. Native andes decoys
    // are built as "<prefix>_<orig>" (the "_" delimiter is required, so a target
    // accession merely starting with the bare prefix is NOT a decoy); the suffix
    // form recognizes externally-built decoys (e.g. quantms "<orig>_rev").
    let suffix = idx.decoy_suffix.clone();
    idx.db.proteins.iter().enumerate().flat_map(move |(p_idx, protein)| {
        let is_decoy = is_decoy_accession_affix(&protein.accession, decoy_prefix, suffix.as_deref());
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

    // N-terminal Met cleavage: when the protein starts with M (and has >1
    // residue), also enumerate candidates treating sequence[1..] as the
    // effective start. The Met-cleaved peptides still carry
    // is_protein_n_term=true (the post-Met residue is the new biological
    // N-terminus) and are NOT deduplicated — they differ by terminal-mod
    // search space.
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

    let cleavage_positions =
        compute_cleavage_positions(sub_seq, params.enzyme, &params.extra_enzymes);

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
                // start and end" concept applies to strictly tryptic spans.
                // For semi-tryptic peptides with a free terminus, the
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
            is_protein_n_term,
            is_protein_c_term,
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
/// - All other positions: Anywhere only — borrowed directly from AminoAcidSet,
///   no clone.
pub(crate) fn expand_mod_combinations(
    span: &[u8],
    params: &SearchParams,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
) -> Vec<Vec<AminoAcid>> {
    let n = span.len();

    // Build owned merged-variant vecs only for the (up to two) terminal
    // positions. Interior positions will borrow the Anywhere slice directly,
    // eliminating the per-position `to_vec` clone that showed up in perf
    // traces (~87% of positions on real tryptic peptides).
    let pos0_owned: Option<Vec<AminoAcid>> = (n > 0).then(|| {
        build_terminal_variants(params, span[0], 0, n, is_protein_n_term, is_protein_c_term)
    });
    let pos_last_owned: Option<Vec<AminoAcid>> = (n > 1).then(|| {
        build_terminal_variants(params, span[n - 1], n - 1, n, is_protein_n_term, is_protein_c_term)
    });

    // Collect per-position variant slices. Terminal positions reference the
    // owned vecs above; interior positions borrow directly from AminoAcidSet.
    // All borrows are valid for the duration of this function.
    let position_variants_refs: Vec<&[AminoAcid]> = span.iter().enumerate().map(|(i, &r)| {
        if i == 0 {
            pos0_owned.as_ref().unwrap().as_slice()
        } else if i == n - 1 {
            // n > 1 guaranteed here because n == 1 means i == 0 == n-1,
            // which is already handled by the first branch.
            pos_last_owned.as_ref().unwrap().as_slice()
        } else {
            // Interior position: borrow Anywhere variants — no clone.
            params.aa_set.variants_for(r, ModLocation::Anywhere)
        }
    }).collect();

    let mut out = Vec::new();
    let mut current = Vec::with_capacity(n);
    expand_recursive(
        &position_variants_refs, 0, &mut current, 0,
        params.max_variable_mods_per_peptide, &mut out,
    );
    out
}

/// Build the merged variant list for a terminal position (pos 0 or pos n-1).
///
/// Only called for the 1-2 terminal positions per span.
fn build_terminal_variants(
    params: &SearchParams,
    residue: u8,
    pos: usize,
    span_len: usize,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
) -> Vec<AminoAcid> {
    let anywhere_variants = params.aa_set.variants_for(residue, ModLocation::Anywhere);

    // Helper: returns true if `term_variants` contains a FIXED mod variant
    // for this residue. When a fixed terminal mod applies, the residue
    // MUST carry it — the unmodified Anywhere variant is not a valid
    // candidate. Fixed mods are mandatory (Kim et al., Nat Commun 5:5277, 2014).
    let has_fixed_in = |term_variants: &[AminoAcid]| -> bool {
        term_variants.iter().any(|aa| aa.mod_.as_ref().map(|m| m.fixed).unwrap_or(false))
    };

    let n_term_variants: &[AminoAcid] = if pos == 0 {
        let loc = if is_protein_n_term { ModLocation::ProtNTerm } else { ModLocation::NTerm };
        params.aa_set.variants_for(residue, loc)
    } else {
        &[]
    };
    let c_term_variants: &[AminoAcid] = if pos == span_len - 1 {
        let loc = if is_protein_c_term { ModLocation::ProtCTerm } else { ModLocation::CTerm };
        params.aa_set.variants_for(residue, loc)
    } else {
        &[]
    };

    let has_fixed_n = has_fixed_in(n_term_variants);
    let has_fixed_c = has_fixed_in(c_term_variants);

    // If a fixed terminal mod is mandatory at this position, the
    // unmodified Anywhere variant is not a legal candidate. Drop the
    // Anywhere variants in that case; otherwise include them. This
    // prevents the candidate explosion that wildcard fixed N-term TMT
    // would otherwise cause (every peptide would be enumerated twice
    // at position 0: once unmodded, once TMT-modded).
    //
    // Note: Anywhere variants always include the residue's own fixed
    // mods folded in (e.g. K-anywhere already carries K-TMT), so this
    // rule applies only to terminal mods.
    let mut variants: Vec<AminoAcid> = if has_fixed_n || has_fixed_c {
        Vec::new()
    } else {
        anywhere_variants.to_vec()
    };

    // Append all terminal variants (fixed + variable). When a fixed
    // mod is present, the modded variant is the only legal one for
    // that mod's residue/location slot; variable mods stack on top
    // by adding additional explored variants.
    for v in n_term_variants {
        if !variants.contains(v) {
            variants.push(v.clone());
        }
    }
    for v in c_term_variants {
        if !variants.contains(v) {
            variants.push(v.clone());
        }
    }

    variants
}

fn expand_recursive(
    position_variants: &[&[AminoAcid]],
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
    for variant in position_variants[pos] {
        // Only VARIABLE mods consume slots against the per-peptide cap.
        // Fixed mods are unconditionally applied by the AminoAcidSet (e.g.
        // CAM-on-C, TMT-on-K, TMT-on-N-term-wildcard) and must not count
        // against max_variable_mods_per_peptide — otherwise a peptide with
        // two fixed mods (e.g. TQAHTQQNMVEK + N-term-TMT + K-TMT) is pruned
        // when NumMods=1, which is exactly the bug that caused 86% of TMT
        // top-1 PSMs to diverge from the reference implementation.
        //
        // Kim et al. (Nat Commun 5:5277, 2014): only variable mods consume
        // slots against the per-peptide modification cap.
        let consumes_slot = variant
            .mod_
            .as_ref()
            .map(|m| !m.fixed)
            .unwrap_or(false);
        let new_mods = mods_used + if consumes_slot { 1 } else { 0 };
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
/// every i in 1..n where any configured protease cuts — `is_cleavable_after(seq[i-1])`
/// (for C-term cutters like Trypsin) OR `is_cleavable_before(seq[i])` (for
/// N-term cutters like AspN/LysN).
///
/// `primary` is the configured `enzyme`; `extras` are additional proteases for
/// a multi-protease digest. A position is a cleavage site if the union of all
/// proteases cuts there. With `extras` empty this is the single-enzyme path
/// (bit-identical to before multi-enzyme support).
///
/// NonSpecific/NoCleavage are whole-protease modes, so they only short-circuit
/// when they apply to *every* configured protease — a NonSpecific extra makes
/// the whole digest non-specific; a NoCleavage extra contributes no cut sites
/// (it never short-circuits unless it's the only protease).
fn compute_cleavage_positions(seq: &[u8], primary: Enzyme, extras: &[Enzyme]) -> Vec<u32> {
    let n = seq.len() as u32;

    let all = || std::iter::once(primary).chain(extras.iter().copied());

    // Any non-specific protease makes the whole digest non-specific.
    if all().any(|e| matches!(e, Enzyme::NonSpecific)) {
        return (0..=n).collect();
    }
    // Only when *every* protease is NoCleavage are there no internal cut sites.
    if all().all(|e| matches!(e, Enzyme::NoCleavage)) {
        return vec![0, n];
    }

    let mut positions = vec![0u32];
    for i in 1..n {
        let prev = seq[(i - 1) as usize];
        let here = seq[i as usize];
        if all().any(|e| e.is_cleavable_after(prev) || e.is_cleavable_before(here)) {
            positions.push(i);
        }
    }
    if *positions.last().unwrap() != n {
        positions.push(n);
    }
    positions
}

/// Distinct variable-mod mass deltas configured in `params.aa_set`.
///
/// Walks every variant in the set, keeps the non-fixed (variable) mods, and
/// de-duplicates by `mass_delta.to_bits()`. Fixed mods are excluded — they are
/// already folded into each residue's mass (and into the base-peptide index),
/// so they never contribute to the lazy mod-mass offset `Δ`.
fn variable_mod_deltas(params: &SearchParams) -> Vec<f64> {
    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut out: Vec<f64> = Vec::new();
    for aa in params.aa_set.iter_variants() {
        if let Some(m) = aa.mod_.as_ref() {
            if !m.fixed && seen.insert(m.mass_delta.to_bits()) {
                out.push(m.mass_delta);
            }
        }
    }
    out
}

/// Enumerate the distinct total mod-mass offsets `Δ` reachable by applying up
/// to `max_variable_mods_per_peptide` variable mods (with repetition — the same
/// mod can apply at multiple residues).
///
/// Returns a list of `(delta_milli, delta)` pairs, including `Δ = 0` (the
/// unmodified backbone), de-duplicated by rounded milli-Dalton so two mods that
/// sum to the same shift collapse to one window probe.
///
/// This is intentionally a **superset** of the offsets any single base peptide
/// can actually achieve: a peptide may lack the target residue for a given mod,
/// in which case that `Δ` simply yields no surviving peptidoform after the
/// per-record mass filter. Over-generating `Δ` only costs extra (filtered-out)
/// window probes; it can never drop a candidate. The decisive correctness
/// guarantee comes from re-running [`expand_mod_combinations`] per base record
/// and filtering on the resulting peptidoform mass — not from `Δ` precision.
fn distinct_mod_mass_offsets(params: &SearchParams) -> Vec<(u64, f64)> {
    let deltas = variable_mod_deltas(params);
    let max_mods = params.max_variable_mods_per_peptide as usize;

    // BFS over multisets of mod deltas: level k = applying exactly k variable
    // mods. We track each combination's summed delta and de-dup by milli-Da.
    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut out: Vec<(u64, f64)> = Vec::new();

    // Δ = 0 (no variable mods) is always present.
    let zero_milli = (0.0_f64 * 1000.0).round() as i64 as u64;
    seen.insert(zero_milli);
    out.push((zero_milli, 0.0));

    if deltas.is_empty() || max_mods == 0 {
        return out;
    }

    // `frontier` holds (summed_delta, start_index) where start_index enforces a
    // non-decreasing choice of mod indices so each multiset is generated once.
    let mut frontier: Vec<(f64, usize)> = vec![(0.0, 0)];
    for _ in 0..max_mods {
        let mut next: Vec<(f64, usize)> = Vec::new();
        for &(sum, start) in &frontier {
            for (i, &d) in deltas.iter().enumerate().skip(start) {
                let new_sum = sum + d;
                next.push((new_sum, i));
                let key = milli_of(new_sum);
                if seen.insert(key) {
                    out.push((key, new_sum));
                }
            }
        }
        frontier = next;
        if frontier.is_empty() {
            break;
        }
    }
    out
}

/// Round a neutral mass (Da) to milli-Daltons, consistent with Task 1's
/// `(mass * 1000.0).round()` index key. Uses `i64` intermediately so negative
/// offsets (e.g. pyro-glu) round correctly before the `u64` window arithmetic.
#[inline]
fn milli_of(mass: f64) -> u64 {
    (mass * 1000.0).round() as i64 as u64
}

/// Lazily enumerate the modified candidate set for a precursor mass window,
/// reading base peptides from the out-of-core index `mi` and placing variable
/// mods on the fly.
///
/// **Decisive correctness property** (guarded by
/// `lazy_enum_matches_inram_candidate_set_for_window`): the SET of `Candidate`s
/// returned here equals
/// `{ c ∈ enumerate_candidates(db, params, decoy_prefix) : |c.peptide.mass() − precursor_mass| ≤ tol_da }`,
/// compared at the peptidoform level (residues + mods, mass, protein_index,
/// start_offset, decoy).
///
/// ## Algorithm
/// 1. Enumerate the distinct mod-mass offsets `Δ` (incl. 0) reachable with up
///    to `params.max_variable_mods_per_peptide` variable mods
///    ([`distinct_mod_mass_offsets`]).
/// 2. For each `Δ`, binary-search the index for base peptides whose **base**
///    mass lands in `[precursor − Δ − tol, precursor − Δ + tol]`. Collect the
///    union of base records across all `Δ`, de-duplicated by
///    `(protein_index, start_offset, length, flags)` so a record reached via
///    several `Δ` windows is expanded only once.
/// 3. For each unique base record, reconstruct the residue span from
///    `db.proteins[protein_index].sequence[start_offset..start_offset+length]`
///    and the terminal context (`pre`/`post` from absolute coordinates;
///    `is_protein_n_term`/`is_protein_c_term` from the stored flags — these are
///    NOT derivable from coordinates alone for Met-cleaved peptides). Then call
///    the SAME [`expand_mod_combinations`] enumerate_candidates uses, yielding
///    every legal peptidoform, and keep those whose actual
///    `peptide.mass()` is within `tol_da` of `precursor_mass`.
///
/// Step 3 is why the result is bit-identical to the in-RAM set: mod placement
/// and mass computation flow through the exact same code path. The per-`Δ`
/// windows in step 2 are purely an out-of-core retrieval optimisation — they
/// only decide which base records to *look at*, never which peptidoforms survive.
pub fn lazy_candidates_for_precursor(
    mi: &crate::candidate_index::MmapCandidateIndex,
    db: &SearchIndex,
    params: &SearchParams,
    precursor_mass: f64,
    tol_left_da: f64,
    tol_right_da: f64,
) -> Vec<Candidate> {
    use crate::candidate_index::{flags, IndexRecord};

    // (1) Distinct mod-mass offsets, including Δ = 0.
    let offsets = distinct_mod_mass_offsets(params);

    // (2) Union of base records across all Δ windows, de-duplicated.
    //
    // DIRECTIONAL tolerance (mirrors `matches_precursor`): a candidate of mass
    // `m` has signed error `m − precursor`; when LIGHTER (error < 0) the bound is
    // `tol_left`, when HEAVIER (error ≥ 0) the bound is `tol_right`. So the base
    // window is `[center − tol_left, center + tol_right]`, NOT a symmetric `±max`.
    // Using `max(left,right)` would pre-fetch extra candidates on the wider side
    // which — even though re-gated by the final mass filter below — perturb the
    // per-spectrum scored-candidate set (tailor histogram / null-pool) and diverge
    // from the in-RAM path's scored set.
    let tol_left_milli = (tol_left_da * 1000.0).round() as u64;
    let tol_right_milli = (tol_right_da * 1000.0).round() as u64;
    let mut seen_records: std::collections::HashSet<(u32, u32, u16, u16)> =
        std::collections::HashSet::new();
    let mut base_records: Vec<IndexRecord> = Vec::new();
    for (_, delta) in &offsets {
        // Base mass = precursor − Δ. Window asymmetric: −tol_left .. +tol_right.
        let center = precursor_mass - delta;
        let lo = (center * 1000.0).round() as i64 - tol_left_milli as i64;
        let hi = (center * 1000.0).round() as i64 + tol_right_milli as i64;
        // A base mass is always positive; clamp the lower bound at 0.
        let lo_milli = (lo - 1).max(0) as u64;
        if hi < 0 {
            continue;
        }
        let hi_milli = (hi + 1) as u64;
        for rec in mi.mass_window(lo_milli, hi_milli) {
            let key = (rec.protein_index, rec.start_offset, rec.length, rec.flags);
            if seen_records.insert(key) {
                base_records.push(rec);
            }
        }
    }

    // (3) Reconstruct each base peptide span + terminal context, place mods via
    //     the SAME expand_mod_combinations path, keep peptidoforms in window.
    let mut out: Vec<Candidate> = Vec::new();
    for rec in &base_records {
        let protein = &db.db.proteins[rec.protein_index as usize];
        let seq = &protein.sequence;
        let abs_start = rec.start_offset as usize;
        let abs_end = abs_start + rec.length as usize;
        // Defensive: a corrupt record could point past the sequence.
        if abs_end > seq.len() {
            continue;
        }
        let span = &seq[abs_start..abs_end];

        let is_protein_n_term = rec.flags & flags::IS_PROTEIN_N_TERM != 0;
        let is_protein_c_term = rec.flags & flags::IS_PROTEIN_C_TERM != 0;
        let is_decoy = rec.flags & flags::IS_DECOY != 0;

        // pre/post are determined by absolute coordinates exactly as emit_span
        // computes them (the N-term-Met case: abs_start==1 ⇒ pre = seq[0] = 'M').
        let pre = if abs_start == 0 { b'_' } else { seq[abs_start - 1] };
        let post = if abs_end == seq.len() { b'-' } else { seq[abs_end] };

        let mod_combinations =
            expand_mod_combinations(span, params, is_protein_n_term, is_protein_c_term);
        for residues in mod_combinations {
            let peptide = Peptide::new(residues, pre, post);
            // DIRECTIONAL final gate, mirroring `matches_precursor`: lighter →
            // tol_left, heavier → tol_right. Keeps the lazily-materialized set
            // EXACTLY the set the in-RAM `matches_precursor` admits per spectrum.
            let err = peptide.mass() - precursor_mass;
            let allowed = if err < 0.0 { tol_left_da } else { tol_right_da };
            if err.abs() <= allowed {
                out.push(Candidate {
                    peptide,
                    protein_index: rec.protein_index as usize,
                    start_offset_in_protein: abs_start,
                    is_decoy,
                    is_protein_n_term,
                    is_protein_c_term,
                });
            }
        }
    }
    out
}

/// Lazily enumerate the modified candidate set whose **nominal residue mass**
/// falls in the inclusive integer-bucket window `[min_nominal, max_nominal]`,
/// reproducing the in-RAM `bucket_index.range(candidate_nominal_bounds(..))`
/// membership predicate BIT-FOR-BIT.
///
/// The in-RAM (`Ram`) backing scores exactly the candidates whose
/// `peptide.nominal_residue_mass()` (= `nominal_from(peptide.mass() − H2O)`)
/// lands in the coarse integer bucket window. That window is LOSSY at bucket
/// boundaries: the f64-milli-Da directional window used by
/// [`lazy_candidates_for_precursor`] is a DIFFERENT set at the edges, so the two
/// backings diverge on boundary candidates (proven on scan 33203). To keep the
/// `Mmap` backing result-identical to `Ram`, the per-spectrum scored set must be
/// gated on the SAME nominal-bucket predicate, NOT the f64 window.
///
/// ## Algorithm
/// Identical to [`lazy_candidates_for_precursor`] (per-Δ base-record fetch →
/// `expand_mod_combinations`) EXCEPT:
/// - The base-record fetch window is derived from the nominal bounds, converted
///   back to a conservative milli-Da mass window that is WIDE ENOUGH to cover
///   every base mass that could expand to an in-window peptidoform (over-fetch
///   is fine — surplus records are filtered out by the acceptance predicate).
/// - Each expanded peptidoform is accepted iff
///   `nominal_from(peptide.mass() − H2O) ∈ [min_nominal, max_nominal]` — the
///   EXACT RAM membership test — replacing the directional f64 final gate.
///
/// The lazy per-Δ enumeration is retained for bounded-RSS (O(window), not full
/// materialization); only the per-candidate ACCEPTANCE predicate changes.
pub fn lazy_candidates_for_nominal_window(
    mi: &crate::candidate_index::MmapCandidateIndex,
    db: &SearchIndex,
    params: &SearchParams,
    min_nominal: i32,
    max_nominal: i32,
) -> Vec<Candidate> {
    use crate::candidate_index::{flags, IndexRecord};
    use model::mass::{H2O, INTEGER_MASS_SCALER};

    if min_nominal > max_nominal {
        return Vec::new();
    }

    // (1) Distinct mod-mass offsets, including Δ = 0 (same superset as the
    //     precursor-window path).
    let offsets = distinct_mod_mass_offsets(params);

    // (2) Conservative full-peptide (incl. H2O) mass window covering the nominal
    //     bounds. A peptidoform with bucket key `k = round(SCALER * residue_mass)`
    //     has `residue_mass ∈ ((k − 0.5)/SCALER, (k + 0.5)/SCALER)`. Pad by one
    //     extra bucket-width on each side to be safe against float rounding, then
    //     add H2O to move from residue-mass space to full-peptide mass space.
    let scaler = INTEGER_MASS_SCALER as f64;
    // One nominal unit spans ~1/SCALER Da of real mass; pad generously.
    let bucket_da = 1.0 / scaler;
    let residue_lo = (min_nominal as f64 - 0.5) / scaler - bucket_da;
    let residue_hi = (max_nominal as f64 + 0.5) / scaler + bucket_da;
    let full_lo = residue_lo + H2O;
    let full_hi = residue_hi + H2O;

    // (3) Union of base records across all Δ windows, de-duplicated. Base mass =
    //     full_peptide_mass − Δ, so the base-mass fetch window is
    //     `[full_lo − Δ, full_hi − Δ]`.
    let mut seen_records: std::collections::HashSet<(u32, u32, u16, u16)> =
        std::collections::HashSet::new();
    let mut base_records: Vec<IndexRecord> = Vec::new();
    for (_, delta) in &offsets {
        let lo = ((full_lo - delta) * 1000.0).floor() as i64;
        let hi = ((full_hi - delta) * 1000.0).ceil() as i64;
        if hi < 0 {
            continue;
        }
        let lo_milli = lo.max(0) as u64;
        let hi_milli = hi as u64;
        for rec in mi.mass_window(lo_milli, hi_milli) {
            let key = (rec.protein_index, rec.start_offset, rec.length, rec.flags);
            if seen_records.insert(key) {
                base_records.push(rec);
            }
        }
    }

    // Sort base_records to match enumerate_candidates order (protein_index ASC,
    // start_offset ASC, flags ASC) before expanding mod combinations. This is
    // required for the MMAP path to visit candidates in the same order as the
    // RAM path (which iterates candidates in enumerate_candidates order via
    // ascending global indices). A HashSet was used for dedup in step (3), so
    // the vec order is non-deterministic without this sort.
    base_records.sort_unstable_by_key(|r| (r.protein_index, r.start_offset, r.flags));

    // (4) Reconstruct + place mods exactly as the precursor-window path, but gate
    //     acceptance on the RAM nominal-bucket predicate.
    let mut out: Vec<Candidate> = Vec::new();
    for rec in &base_records {
        let protein = &db.db.proteins[rec.protein_index as usize];
        let seq = &protein.sequence;
        let abs_start = rec.start_offset as usize;
        let abs_end = abs_start + rec.length as usize;
        if abs_end > seq.len() {
            continue;
        }
        let span = &seq[abs_start..abs_end];

        let is_protein_n_term = rec.flags & flags::IS_PROTEIN_N_TERM != 0;
        let is_protein_c_term = rec.flags & flags::IS_PROTEIN_C_TERM != 0;
        let is_decoy = rec.flags & flags::IS_DECOY != 0;

        let pre = if abs_start == 0 { b'_' } else { seq[abs_start - 1] };
        let post = if abs_end == seq.len() { b'-' } else { seq[abs_end] };

        let mod_combinations =
            expand_mod_combinations(span, params, is_protein_n_term, is_protein_c_term);
        for residues in mod_combinations {
            let peptide = Peptide::new(residues, pre, post);
            // RAM membership: nominal_from(peptide.mass() − H2O) in bounds.
            // `Peptide::nominal_residue_mass()` caches exactly this value.
            let key = peptide.nominal_residue_mass();
            if key >= min_nominal && key <= max_nominal {
                out.push(Candidate {
                    peptide,
                    protein_index: rec.protein_index as usize,
                    start_offset_in_protein: abs_start,
                    is_decoy,
                    is_protein_n_term,
                    is_protein_c_term,
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    #[test]
    fn decoy_membership_requires_underscore_delimiter() {
        use crate::decoy::is_decoy_accession;
        // Decoys are built as "<prefix>_<orig>" (reverse_db/shuffle_db), so
        // membership requires the "_" delimiter, and candidate generation MUST
        // use the SAME needle as SearchIndex (they share decoy_accession_needle)
        // so they cannot diverge. A target accession that merely starts with the
        // bare prefix (e.g. a real protein "XXXfoo") is NOT a decoy.
        assert!(is_decoy_accession("XXX_protein1", "XXX"));
        assert!(
            !is_decoy_accession("XXXprotein1", "XXX"),
            "no '_' delimiter ⇒ a real 'XXX...' target must NOT be read as decoy"
        );
        assert!(!is_decoy_accession("DECOY_protein1", "XXX"));
        // Short custom prefix: a genuine "rev..." target is not a decoy.
        assert!(is_decoy_accession("rev_sp|P1", "rev"));
        assert!(!is_decoy_accession("reverse_kinase", "rev"));
    }

    #[test]
    fn multi_enzyme_cleavage_is_union_of_rules() {
        use super::compute_cleavage_positions;
        use model::enzyme::Enzyme;
        // M K A E S R P E  (n=8). Trypsin cuts after K,R; GluC cuts after E.
        let seq = b"MKAESRPE";

        let trypsin = compute_cleavage_positions(seq, Enzyme::Trypsin, &[]);
        assert_eq!(trypsin, vec![0, 2, 6, 8], "Trypsin alone: after K(2), R(6)");

        let gluc = compute_cleavage_positions(seq, Enzyme::GluC, &[]);
        assert_eq!(gluc, vec![0, 4, 8], "GluC alone: after first E(4)");

        // Union: every Trypsin site AND every GluC site is a cut point.
        let both = compute_cleavage_positions(seq, Enzyme::GluC, &[Enzyme::Trypsin]);
        assert_eq!(both, vec![0, 2, 4, 6, 8], "GluC+Trypsin = union of cut sites");
        // Order of primary vs extra must not matter for the cut-site set.
        let both_rev = compute_cleavage_positions(seq, Enzyme::Trypsin, &[Enzyme::GluC]);
        assert_eq!(both, both_rev, "union is order-independent");
    }

    #[test]
    fn single_enzyme_path_unchanged_with_empty_extras() {
        use super::compute_cleavage_positions;
        use model::enzyme::Enzyme;
        // Empty extras must be bit-identical to the legacy single-enzyme result
        // for every protease mode, including the NonSpecific/NoCleavage cases.
        let seq = b"MKAESRPEDK";
        for e in [Enzyme::Trypsin, Enzyme::LysC, Enzyme::GluC, Enzyme::AspN,
                  Enzyme::NonSpecific, Enzyme::NoCleavage] {
            let single = compute_cleavage_positions(seq, e, &[]);
            // NoCleavage extras contribute no sites; result must equal the primary alone.
            let with_nocleave = compute_cleavage_positions(seq, e, &[Enzyme::NoCleavage]);
            assert_eq!(single, with_nocleave, "{e:?}: NoCleavage extra adds no sites");
        }
        // A NonSpecific protease anywhere in the set makes the whole digest non-specific.
        let ns = compute_cleavage_positions(seq, Enzyme::Trypsin, &[Enzyme::NonSpecific]);
        assert_eq!(ns, (0..=seq.len() as u32).collect::<Vec<_>>());
    }

    // ───────────────────────── Task 3: lazy modified enumeration ─────────────
    //
    // Decisive correctness property: for a precursor window, the SET of
    // Candidates from `lazy_candidates_for_precursor` (which reads BASE peptides
    // from the out-of-core index and places mods lazily) must EXACTLY equal the
    // mass-filtered subset of the in-RAM `enumerate_candidates`.

    use super::{
        enumerate_candidates, lazy_candidates_for_nominal_window, lazy_candidates_for_precursor,
        Candidate,
    };
    use crate::candidate_index::{build_base_peptide_index, MmapCandidateIndex};
    use crate::search_index::SearchIndex;
    use crate::search_params::SearchParams;
    use model::aa_set::AminoAcidSetBuilder;
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use model::protein::{Protein, ProteinDb};

    /// Canonical identity of a candidate peptidoform — independent of any
    /// enumeration order. Captures decoy, protein coordinates, the residue+mod
    /// signature (via the Display form `pre.SEQ±delta.post`), and the rounded
    /// mass. Two candidates with this key equal are the same peptidoform.
    fn canonical_key(c: &Candidate) -> (bool, usize, usize, String, u64) {
        (
            c.is_decoy,
            c.protein_index,
            c.start_offset_in_protein,
            c.peptide.to_string(),
            (c.peptide.mass() * 1000.0).round() as u64,
        )
    }

    fn base_mass_of(db: &SearchIndex, params: &SearchParams, seq: &[u8]) -> f64 {
        // Find the unmodified candidate whose residues spell `seq` and return
        // its base mass (Δ = 0). Uses enumerate_candidates with NumMods forced
        // to 0 so we get exactly the unmodified backbone.
        let base_params = SearchParams {
            max_variable_mods_per_peptide: 0,
            ..params.clone()
        };
        let cands: Vec<_> = enumerate_candidates(db, &base_params, "XXX").collect();
        cands
            .into_iter()
            .filter(|c| !c.is_decoy)
            .find(|c| {
                c.peptide.residues.iter().map(|aa| aa.residue).eq(seq.iter().copied())
            })
            .unwrap_or_else(|| panic!("no unmodified candidate spelling {:?}", std::str::from_utf8(seq)))
            .peptide
            .mass()
    }

    fn build_and_open(
        db: &SearchIndex,
        params: &SearchParams,
    ) -> (tempfile::NamedTempFile, MmapCandidateIndex) {
        use std::io::BufWriter;
        let tmp = tempfile::NamedTempFile::new().expect("tempfile");
        {
            let mut bw = BufWriter::new(tmp.as_file());
            build_base_peptide_index(db, params, "XXX", &mut bw)
                .expect("build_base_peptide_index failed");
        }
        let mi = MmapCandidateIndex::open(tmp.path()).expect("open index");
        (tmp, mi)
    }

    fn oxidation_m() -> Modification {
        Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.994915,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        }
    }

    fn deamidation_n() -> Modification {
        Modification {
            name: "Deamidation".to_string(),
            mass_delta: 0.984016,
            residue: ResidueSpec::Specific(b'N'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        }
    }

    /// SearchIndex + params with Oxidation-M as the only variable mod, NumMods=1.
    /// Protein contains M (PEPTMK has one M) plus other tryptic peptides.
    fn oxm_fixture() -> (SearchIndex, SearchParams) {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "ox-m fixture".into(),
                // PEPTMK | DALMR | QPK  → tryptic peptides with M present.
                sequence: b"PEPTMKDALMRQPK".to_vec(),
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(oxidation_m())
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 1;
        (db, params)
    }

    /// Set-equality assertion shared by both tests: the lazy candidate set for
    /// `(precursor, tol)` must equal the mass-filtered in-RAM set.
    fn assert_lazy_matches_inram(
        db: &SearchIndex,
        params: &SearchParams,
        mi: &MmapCandidateIndex,
        precursor: f64,
        tol: f64,
    ) {
        use std::collections::HashSet;
        let lazy: HashSet<_> = lazy_candidates_for_precursor(mi, db, params, precursor, tol, tol)
            .iter()
            .map(canonical_key)
            .collect();
        let inram: HashSet<_> = enumerate_candidates(db, params, "XXX")
            .filter(|c| (c.peptide.mass() - precursor).abs() <= tol)
            .map(|c| canonical_key(&c))
            .collect();
        assert_eq!(
            lazy, inram,
            "lazy enumeration must match the in-RAM candidate set\n\
             only-in-lazy: {:?}\n only-in-inram: {:?}",
            lazy.difference(&inram).collect::<Vec<_>>(),
            inram.difference(&lazy).collect::<Vec<_>>(),
        );
    }

    // ───────── Task 4 fix: mmap nominal window == RAM bucket window ──────────
    //
    // RESULT-IDENTITY regression guard. The in-RAM (`Ram`) backing scores
    // exactly `{ c : nominal_from(c.mass − H2O) ∈ candidate_nominal_bounds(spec,z) }`
    // (its `bucket_index.range(min..=max)` set). The `Mmap` backing previously
    // gated on an f64 directional milli-Da window, which is a DIFFERENT set at
    // integer-bucket boundaries — the proven scan-33203 divergence. This test
    // asserts the new `lazy_candidates_for_nominal_window` reproduces RAM's
    // bucket-window membership BIT-FOR-BIT, and (fail-before proof) that the OLD
    // f64-window path diverged from RAM on a boundary fixture.

    /// RAM membership set for a spectrum at charge `z`: every `enumerate_candidates`
    /// candidate whose nominal residue mass is in `candidate_nominal_bounds`.
    /// This is precisely what `bucket_index.range(min..=max)` selects.
    fn ram_window_keys(
        db: &SearchIndex,
        params: &SearchParams,
        spec: &model::Spectrum,
        z: u8,
    ) -> std::collections::HashSet<(bool, usize, usize, String, u64)> {
        let (min_nominal, max_nominal) =
            crate::match_engine::candidate_nominal_bounds(spec, z, params, 0.0);
        enumerate_candidates(db, params, "XXX")
            .filter(|c| {
                let k = c.peptide.nominal_residue_mass();
                k >= min_nominal && k <= max_nominal
            })
            .map(|c| canonical_key(&c))
            .collect()
    }

    /// A spectrum whose precursor m/z corresponds (at the given charge, iso 0) to
    /// the supplied neutral peptide mass.
    fn spectrum_at_mass(neutral_mass: f64, z: u8) -> model::Spectrum {
        let mz = (neutral_mass + z as f64 * model::mass::PROTON) / z as f64;
        model::Spectrum {
            title: "boundary".into(),
            precursor_mz: mz,
            precursor_intensity: None,
            precursor_charge: Some(z as i32),
            rt_seconds: None,
            scan: None,
            peaks: Vec::new(),
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        }
    }

    #[test]
    fn mmap_nominal_window_equals_ram_bucket_window_on_boundary() {
        use model::tolerance::{PrecursorTolerance, Tolerance};
        use std::collections::HashSet;

        let (db, mut params) = oxm_fixture();
        // No isotope offsets: isolate the window-breadth divergence at a single
        // bucket edge (offsets only shift the bounds; the boundary effect is the
        // integer-vs-f64 widening at the edge).
        params.isotope_error_range = 0..=0;
        // A fractional Da tolerance whose `round(tol − 0.4999)` widening differs
        // from the exact f64 window at bucket boundaries. 0.5 Da → widen by 0,
        // so the nominal window is a single bucket while the f64 window spans
        // ±0.5 Da and can dip into the neighboring bucket.
        params.precursor_tolerance =
            PrecursorTolerance::symmetric(Tolerance::Da(0.5));
        let (_tmp, mi) = build_and_open(&db, &params);

        let z = 1u8;

        // Collect every candidate's true neutral mass so we can place precursors
        // on bucket boundaries (where exact f64 mass crosses an integer edge).
        let masses: Vec<f64> = enumerate_candidates(&db, &params, "XXX")
            .map(|c| c.peptide.mass())
            .collect();
        assert!(!masses.is_empty(), "fixture produced no candidates");

        // Scan a fine grid of precursor masses straddling each candidate's mass.
        // We assert (1) the NEW nominal-window path always equals RAM's bucket
        // window, and (2) detect at least one precursor where the OLD f64 window
        // DIVERGED from RAM — proving the boundary case is exercised (fail-before).
        let mut boundary_seen = false;
        for &m in &masses {
            for step in -8i32..=8 {
                let precursor = m + (step as f64) * 0.125; // 0.125 Da grid
                let spec = spectrum_at_mass(precursor, z);

                let ram = ram_window_keys(&db, &params, &spec, z);

                // NEW path: nominal-bucket gated mmap window.
                let (min_nominal, max_nominal) =
                    crate::match_engine::candidate_nominal_bounds(&spec, z, &params, 0.0);
                let mmap_nominal: HashSet<_> =
                    lazy_candidates_for_nominal_window(&mi, &db, &params, min_nominal, max_nominal)
                        .iter()
                        .map(canonical_key)
                        .collect();
                assert_eq!(
                    mmap_nominal, ram,
                    "NEW mmap nominal window must equal RAM bucket window \
                     (precursor={precursor})\n only-in-mmap: {:?}\n only-in-ram: {:?}",
                    mmap_nominal.difference(&ram).collect::<Vec<_>>(),
                    ram.difference(&mmap_nominal).collect::<Vec<_>>(),
                );

                // OLD path: f64 directional milli-Da window (what shipped before
                // this fix). Reproduce exactly how the mmap loop called it.
                let neutral_obs = precursor; // shift_ppm = 0, iso 0
                let tol_left = params.precursor_tolerance.left.as_da(neutral_obs);
                let tol_right = params.precursor_tolerance.right.as_da(neutral_obs);
                let old_f64: HashSet<_> = lazy_candidates_for_precursor(
                    &mi, &db, &params, neutral_obs, tol_left, tol_right,
                )
                .iter()
                .map(canonical_key)
                .collect();
                if old_f64 != ram {
                    boundary_seen = true;
                }
            }
        }
        assert!(
            boundary_seen,
            "fixture never exercised a bucket boundary where the OLD f64 window \
             diverged from RAM — the regression guard would not have caught the bug"
        );
    }

    #[test]
    fn lazy_enum_matches_inram_candidate_set_for_window() {
        let (db, params) = oxm_fixture();
        let (_tmp, mi) = build_and_open(&db, &params);

        // (a) Ox-M precursor: base mass of PEPTMK + one oxidation. The Ox-M
        //     peptidoform must appear (and the unmodified form must NOT, since
        //     its base mass is 15.99 Da lighter and outside a 0.01 window).
        let base = base_mass_of(&db, &params, b"PEPTMK");
        let p_oxm = base + 15.994915;
        assert_lazy_matches_inram(&db, &params, &mi, p_oxm, 0.01);

        // Spot-check: the Ox-M peptidoform is actually present.
        let lazy = lazy_candidates_for_precursor(&mi, &db, &params, p_oxm, 0.01, 0.01);
        assert!(
            lazy.iter().any(|c| {
                !c.is_decoy
                    && c.peptide.residues.iter().map(|aa| aa.residue).eq(b"PEPTMK".iter().copied())
                    && c.peptide.residues.iter().any(|aa| aa.is_modified())
            }),
            "expected the Ox-M PEPTMK peptidoform"
        );

        // (b) Base (unmodified) precursor: the unmodified PEPTMK must appear.
        assert_lazy_matches_inram(&db, &params, &mi, base, 0.01);
        let lazy_base = lazy_candidates_for_precursor(&mi, &db, &params, base, 0.01, 0.01);
        assert!(
            lazy_base.iter().any(|c| {
                !c.is_decoy
                    && c.peptide.residues.iter().map(|aa| aa.residue).eq(b"PEPTMK".iter().copied())
                    && c.peptide.residues.iter().all(|aa| !aa.is_modified())
            }),
            "expected the unmodified PEPTMK peptidoform"
        );

        // (c) A wide window covering many masses — stresses the union/dedup of
        //     base records across Δ windows.
        assert_lazy_matches_inram(&db, &params, &mi, base, 200.0);
    }

    #[test]
    fn lazy_enum_matches_inram_two_variable_mods() {
        // Ox-M + Deamidation-N, NumMods=2 — stresses multi-Δ combination logic
        // (Δ ∈ {0, Ox, Deam, 2·Ox, Ox+Deam, 2·Deam}).
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P2".into(),
                description: "two-mod fixture".into(),
                // Tryptic peptides containing M and N (sometimes both).
                sequence: b"MNPEMTKNMDALERQNMPK".to_vec(),
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(oxidation_m())
            .add_variable_mod(deamidation_n())
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 2;

        let (_tmp, mi) = build_and_open(&db, &params);

        // Sweep a series of precursor masses across the full enumerated range so
        // every Δ combination (incl. 2 mods on different/same residues, decoys)
        // is exercised. For each enumerated candidate mass, assert set-equality
        // of the lazy vs in-RAM window centred there.
        let all: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
        assert!(!all.is_empty(), "fixture must enumerate candidates");
        // Distinct rounded masses to probe.
        let mut masses: Vec<f64> = all.iter().map(|c| c.peptide.mass()).collect();
        masses.sort_by(|a, b| a.partial_cmp(b).unwrap());
        masses.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        for &m in &masses {
            assert_lazy_matches_inram(&db, &params, &mi, m, 0.01);
        }

        // And one wide window covering everything (max stress on dedup).
        let max_mass = masses.last().copied().unwrap();
        let min_mass = masses.first().copied().unwrap();
        let mid = (min_mass + max_mass) / 2.0;
        let half = (max_mass - min_mass) / 2.0 + 1.0;
        assert_lazy_matches_inram(&db, &params, &mi, mid, half);
    }
    /// Regression test for the tolerance-boundary set-equality gap:
    ///
    /// `lazy_candidates_for_precursor` silently drops genuine candidates when a
    /// base peptide's milli key falls just outside the integer-rounded fetch
    /// window even though the peptidoform mass satisfies the float filter.
    ///
    /// ## Root cause
    ///
    /// The FETCH window bounds are computed as:
    ///
    /// ```text
    /// lo = round(center × 1000) − round(tol × 1000)
    /// hi = round(center × 1000) + round(tol × 1000)
    /// ```
    ///
    /// A base record is stored at key K = round(base_mass × 1000). When the
    /// fractional parts of `center × 1000` and `tol × 1000` combine adversely,
    /// K can land 1 milli-Da outside [lo, hi] even though |base_mass − center|
    /// ≤ tol.
    ///
    /// ## Concrete construction (hi-boundary failure)
    ///
    /// PEPTMK has base mass B where frac(B × 1000) ≈ 0.812.  With Ox-M delta
    /// (Δ = 15.994915 Da) and tol constructed so that:
    ///
    /// ```text
    /// frac(tol × 1000) ∈ (frac(B × 1000) − 0.5, 0.5)   →   tol rounds DOWN
    /// center = B − tol * 0.999   (B is near the high edge of the Δ=Ox window)
    /// frac(center × 1000)  =  frac(B × 1000) − frac(tol × 1000) * 0.999 < 0.5
    ///                                                              → rounds DOWN
    /// ```
    ///
    /// Under these conditions:
    ///
    /// ```text
    /// hi      = round(center × 1000) + round(tol × 1000)
    ///         = (floor(B × 1000) − 1) + round(tol × 1000)    (both round DOWN)
    ///         = K − 1                                          ← K > hi
    /// K       = round(B × 1000)                               ← MISSED by integer fetch
    /// ```
    ///
    /// The float filter still passes: |B − center| = tol × 0.999 < tol.
    ///
    /// ## Fix
    ///
    /// Widen the fetch window by ±1 milli-Da:
    ///
    /// ```text
    /// lo_milli = (lo − 1).max(0) as u64
    /// hi_milli = (hi + 1)        as u64
    /// ```
    ///
    /// The extra milli is benign: the float filter in step 3 discards any false
    /// positives introduced by the wider window.
    #[test]
    fn lazy_enum_matches_inram_at_tolerance_edge() {
        // Use the oxm_fixture (PEPTMKDALMRQPK, Oxidation-M variable, NumMods=1).
        let (db, params) = oxm_fixture();
        let (_tmp, mi) = build_and_open(&db, &params);

        // Actual monoisotopic base mass of PEPTMK (no variable mods).
        let base = base_mass_of(&db, &params, b"PEPTMK");

        // frac(base × 1000) for PEPTMK is ≈ 0.812.  The hi-boundary bug fires
        // when frac(tol × 1000) ∈ (frac(base × 1000) − 0.5, 0.5), i.e., when
        // frac(base × 1000) − frac(tol × 1000) < 0.5 so that center (= base − tol)
        // rounds DOWN and hi = round(center×1000) + round(tol×1000) < round(base×1000).
        //
        // Pick frac_tol = 0.40, which satisfies (0.812 − 0.5, 0.5) = (0.312, 0.5).
        let frac_tol = 0.40_f64; // ∈ (frac_base − 0.5, 0.5); verified analytically.
        let tol = (10.0 + frac_tol) / 1000.0; // 0.0104 — frac_tol < 0.5 → rounds DOWN

        // Place center so PEPTMK sits just inside the high edge of the Δ=Ox window.
        // Using 0.999 × tol keeps the float distance < tol (avoids float-equality
        // hazards at the exact boundary with PEPTMK+Ox peptidoform mass).
        let ox_delta = 15.994915_f64;
        let center = base - tol * 0.999;
        let precursor = center + ox_delta;

        // Verify hi-boundary mismatch holds (documents that the bug is real here).
        {
            let tol_milli = (tol * 1000.0).round() as i64;
            let center_milli = (center * 1000.0).round() as i64;
            let base_key = (base * 1000.0).round() as i64;
            let hi = center_milli + tol_milli;
            assert!(
                base_key > hi,
                "expected hi-boundary mismatch (base_key={base_key} > hi={hi})                  — rounding conditions may not hold for this mass; verify frac_tol"
            );
        }

        // Sanity: the in-RAM candidate set must be non-empty (PEPTMK+Ox is within tol).
        let inram_count: usize = enumerate_candidates(&db, &params, "XXX")
            .filter(|c| (c.peptide.mass() - precursor).abs() <= tol)
            .count();
        assert!(
            inram_count > 0,
            "fixture must yield at least one in-RAM candidate at the engineered precursor              (base={base:.6}, tol={tol:.6}, center={center:.6}, precursor={precursor:.6})"
        );

        // Decisive assertion — must match in-RAM exactly.
        // FAILS before the ±1 milli fetch-window fix; PASSES after.
        assert_lazy_matches_inram(&db, &params, &mi, precursor, tol);
    }
}
