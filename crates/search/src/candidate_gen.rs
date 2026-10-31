//! Candidate peptide enumeration via per-protein walk.
//!
//! Enumerates enzyme-cleaved spans within the configured length range,
//! including missed-cleavage spans (governed by the `missed_count` check).
//!
//! ## N-terminal Met cleavage
//!
//! When a protein starts with M, a second walk treats `sequence[1..]` as the
//! effective protein sequence (initial Met loss) and emits only the peptides
//! starting at its first residue, with `is_protein_n_term=true` (so the
//! protein-N-term mod forms apply). The full walk leaves those spans to it,
//! so every span of the protein is enumerated once.
//! A mature-chain entry (UniProt XML) is never Met-cleaved: its first residue
//! is not an initiator Met.
//!
//! ## Database annotation (UniProt XML)
//!
//! A sequence-variant protein ([`Protein::variant_parts`]) only yields spans
//! that carry its variant; the rest of it reads like the canonical protein,
//! which is searched already. A mature-chain entry ([`Protein::chain_ends`])
//! only yields the spans at its new termini. A variable mod declared
//! `sites=annotated` ([`model::Modification::annotated_site_mass`]) is only
//! placed on residues the span's protein annotates with a modification of
//! that mass ([`Protein::mod_sites`]); every expansion path below applies the
//! same test through [`SpanSites`], so the in-RAM, lazy and fragment-index
//! candidate sets stay identical.

use crate::decoy::is_decoy_accession_affix;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use model::aa_set::{form_index, AminoAcidSet, Terminus, FORM_TABLE_LEN};
use model::amino_acid::AminoAcid;
use model::enzyme::Enzyme;
use model::peptide::Peptide;
use model::protein::{ChainEnds, ModSite, Protein, Substitution, VariantRegion};

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
    let require_sequon = params.require_nxst_sequon;
    idx.db
        .proteins
        .iter()
        .enumerate()
        .flat_map(move |(p_idx, protein)| {
            let is_decoy =
                is_decoy_accession_affix(&protein.accession, decoy_prefix, suffix.as_deref());
            enumerate_protein(protein, p_idx, is_decoy, params).into_iter()
        })
        .filter(move |c| !require_sequon || candidate_has_nxst_sequon(c))
}

/// The glyco path's sequon membership test (`glyco_search::GlycoCandidates::from_out_of_core`):
/// an internal N-X-S/T, or one completed by the residue after the peptide.
/// Kept in one place so `SearchParams::require_nxst_sequon` can never admit a
/// different candidate set than the scorer would use.
pub fn candidate_has_nxst_sequon(c: &Candidate) -> bool {
    let res: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
    span_has_nxst_sequon(&res, c.peptide.post)
}

/// [`candidate_has_nxst_sequon`] on a bare residue span and the residue after
/// it (`b'-'` at the protein C-terminus). Variable mods never change a residue
/// letter, so every peptidoform of a span shares the span's answer.
pub fn span_has_nxst_sequon(residues: &[u8], post: u8) -> bool {
    use andes_glyco::sequon::{boundary_nxst_site, has_nxst_sequon};
    has_nxst_sequon(residues) || boundary_nxst_site(residues, post).is_some()
}

/// The O-glycan path's candidate membership test: any Ser/Thr in the span. Unlike
/// N-glycan there is no sequon and no boundary case — a bare S/T is the motif.
pub fn candidate_has_ser_or_thr(c: &Candidate) -> bool {
    use andes_glyco::sequon::has_ser_or_thr;
    let res: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
    has_ser_or_thr(&res)
}

/// [`candidate_has_ser_or_thr`] on a bare residue span. The `post` residue is
/// unused (no boundary case for O-glycans) but the signature matches
/// [`span_has_nxst_sequon`] so both slot into `out_of_core_candidates_where`.
pub fn span_has_ser_or_thr(residues: &[u8], _post: u8) -> bool {
    use andes_glyco::sequon::has_ser_or_thr;
    has_ser_or_thr(residues)
}

fn enumerate_protein(
    protein: &Protein,
    protein_index: usize,
    is_decoy: bool,
    params: &SearchParams,
) -> Vec<Candidate> {
    // Standard enumeration: full sequence from offset 0.
    let mut out = enumerate_protein_from_offset(protein, 0, protein_index, is_decoy, params);

    // N-terminal Met cleavage: when the protein starts with M (and has >1
    // residue), also enumerate the peptides that start at the post-Met
    // residue, the new biological N-terminus (is_protein_n_term=true). Every
    // other span of sequence[1..] is emitted by the walk above, which leaves
    // the post-Met spans to this one. A mature chain was cut out of its
    // precursor; its first residue is not an initiator Met.
    if has_met_cleaved_walk(protein) {
        out.extend(enumerate_protein_from_offset(
            protein,
            1,
            protein_index,
            is_decoy,
            params,
        ));
    }

    out
}

/// Whether `protein` gets the Met-cleaved walk: it starts with M, has more
/// than one residue, and is not a mature chain.
#[inline]
pub(crate) fn has_met_cleaved_walk(protein: &Protein) -> bool {
    let seq = &protein.sequence;
    seq.len() > 1 && seq[0] == b'M' && protein.chain_ends().is_none()
}

/// Enumerate candidates starting from `seq_offset` into `seq`.
///
/// `seq_offset = 0` → normal full-protein walk.
/// `seq_offset = 1` → Met-cleaved walk: `seq[1..]` is the effective protein
///   sequence, and only the spans starting at its first residue are emitted
///   (with `is_protein_n_term = true`; the `pre` residue is the cleaved `M`).
///   Cleavage positions, lengths, and missed-cleavage counts are computed over
///   the sub-sequence. The `start_offset_in_protein` stored on each
///   `Candidate` is in the original protein coordinates (`sub_start + 1`).
///
/// The `params.num_tolerable_termini` field controls cleavage enforcement:
/// - `2`: both ends must be enzyme-cleavage sites (strict / fully specific, default).
/// - `1`: at least one end must be an enzyme-cleavage site (semi-specific).
/// - `0`: neither end needs to be a cleavage site (non-specific).
fn enumerate_protein_from_offset(
    protein: &Protein,
    seq_offset: usize,
    protein_index: usize,
    is_decoy: bool,
    params: &SearchParams,
) -> Vec<Candidate> {
    let seq = &protein.sequence;
    let sub_seq = &seq[seq_offset..];
    let n = sub_seq.len() as u32;
    if n < params.min_length {
        return Vec::new();
    }

    // A variant protein only yields spans that cover its variant (`emit_span`),
    // and those lie within `max_length` of it: walk just that window of
    // `sub_seq` instead of the whole protein.
    let window = match protein.variant_region() {
        Some(v) => {
            let reach = params.max_length as usize;
            // `start > end` for an indel in a repeat (see `VariantRegion`).
            let lo = (v.start.min(v.end) as usize).saturating_sub(reach);
            let hi = v.start.max(v.end) as usize + reach;
            (
                lo.saturating_sub(seq_offset) as u32,
                hi.saturating_sub(seq_offset).min(n as usize) as u32,
            )
        }
        None => (0, n),
    };

    let ntt = params.num_tolerable_termini;
    // The Met-cleaved walk emits only spans that start at the new N-terminus.
    let n_term_only = seq_offset > 0;

    // A non-specific digest (ntt=0, or `Enzyme::NonSpecific` as the primary or an
    // extra protease) emits every valid-length span. The missed-cleavage cap does
    // not apply: every position is a cleavage site, so a span of length L has
    // L - 1 internal sites and the cap would only act as a second length limit
    // (with the default cap of 1 and `min_length` 8 it removed every candidate).
    if ntt == 0
        || params.enzyme == Enzyme::NonSpecific
        || params.extra_enzymes.contains(&Enzyme::NonSpecific)
    {
        let ctx = EmitCtx {
            sub_seq,
            seq,
            seq_offset,
            protein_index,
            is_decoy,
            params,
            variant: protein.variant_parts(),
            chain: protein.chain_ends(),
            protein,
            met_cleaved_walk: has_met_cleaved_walk(protein),
        };
        let starts_end = if n_term_only {
            window.1.min(1)
        } else {
            window.1
        };
        return enumerate_all_spans(&ctx, window.0..starts_end, n);
    }

    let cleavage_positions =
        compute_cleavage_positions(sub_seq, params.enzyme, &params.extra_enzymes);
    // The cleavage positions inside the window, anchoring every span walked
    // below. A span inside the window has all of its internal cleavage
    // positions inside it too, so missed-cleavage counts are unchanged.
    let anchors = &cleavage_positions[cleavage_positions.partition_point(|&c| c < window.0)
        ..cleavage_positions.partition_point(|&c| c <= window.1)];

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

    let ctx = EmitCtx {
        sub_seq,
        seq,
        seq_offset,
        protein_index,
        is_decoy,
        params,
        variant: protein.variant_parts(),
        chain: protein.chain_ends(),
        protein,
        met_cleaved_walk: has_met_cleaved_walk(protein),
    };

    // ── Strict spans (ntt=2 behaviour) ───────────────────────────────────────
    // Also included in ntt=1, since a strict span satisfies "at least one end".
    for (i, &start) in anchors.iter().enumerate() {
        if n_term_only && start != 0 {
            break;
        }
        for (offset, &end) in anchors[i + 1..].iter().enumerate() {
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
        for &start in anchors {
            if n_term_only && start != 0 {
                break;
            }
            let c_min = start + params.min_length;
            let c_max = (start + params.max_length).min(n);
            for end in c_min..=c_max {
                // Skip ends that are cleavage positions — already emitted above.
                if cleavage_set.contains(&end) {
                    continue;
                }
                // Bound internal missed cleavages the same way strict spans are
                // bounded: a semi-tryptic peptide with more than
                // `max_missed_cleavages` internal K/R sites is implausible and
                // bloats the candidate space.
                if missed_cleavages_between(&cleavage_positions, start, end)
                    > params.max_missed_cleavages
                {
                    continue;
                }
                emit_span(&ctx, start, end, &mut out);
            }
        }

        // (c) Free N-terminus, tryptic C-terminus. Its starts are never
        // cleavage positions, so never the N-terminus of a Met-cleaved walk.
        for &end in anchors.iter().filter(|_| !n_term_only) {
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
                // Same internal missed-cleavage bound as the strict path.
                if missed_cleavages_between(&cleavage_positions, start, end)
                    > params.max_missed_cleavages
                {
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
    /// The parts of the protein's variant, when it is a sequence-variant
    /// entry (empty otherwise).
    variant: &'a [VariantRegion],
    /// The protein's new termini, when it is a mature-chain entry.
    chain: Option<ChainEnds>,
    /// The protein (annotated sites and variant substitutions).
    protein: &'a Protein,
    /// Whether the protein gets the Met-cleaved walk ([`has_met_cleaved_walk`]).
    met_cleaved_walk: bool,
}

/// Emit a single (start, end) span as candidates, if the span passes residue
/// validity checks. Appends to `out`.
#[inline]
fn emit_span(ctx: &EmitCtx<'_>, start: u32, end: u32, out: &mut Vec<Candidate>) {
    let abs_start = start as usize + ctx.seq_offset;
    let abs_end = end as usize + ctx.seq_offset;
    // A variant protein contributes only the spans that carry its variant, a
    // mature chain only those at its new termini.
    if (!ctx.variant.is_empty()
        && !ctx
            .variant
            .iter()
            .any(|v| v.covered_by(abs_start as u32, abs_end as u32)))
        || ctx
            .chain
            .is_some_and(|c| !c.covered_by(abs_start as u32, abs_end as u32, ctx.seq.len() as u32))
    {
        return;
    }
    let span = &ctx.sub_seq[start as usize..end as usize];
    // Skip spans containing non-standard residues.
    if span.iter().any(|&r| AminoAcid::standard(r).is_none()) {
        return;
    }
    // A span starting at the post-Met residue of an M-initial protein belongs
    // to the Met-cleaved walk (as a protein N-terminal peptide), not this one.
    if ctx.seq_offset == 0 && abs_start == 1 && ctx.met_cleaved_walk {
        return;
    }
    let pre = if abs_start == 0 {
        b'_'
    } else {
        ctx.seq[abs_start - 1]
    };
    let post = if abs_end == ctx.seq.len() {
        b'-'
    } else {
        ctx.seq[abs_end]
    };

    let is_protein_n_term = start == 0;
    let is_protein_c_term = abs_end == ctx.seq.len();
    let sites = SpanSites::new(ctx.protein, abs_start);
    let mod_combinations = expand_mod_combinations(
        span,
        ctx.params,
        is_protein_n_term,
        is_protein_c_term,
        sites,
    );
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

/// Enumerate all valid-length spans starting in `starts` without cleavage
/// constraints or a missed-cleavage cap (non-specific digest).
fn enumerate_all_spans(ctx: &EmitCtx<'_>, starts: std::ops::Range<u32>, n: u32) -> Vec<Candidate> {
    let mut out = Vec::new();
    for start in starts {
        let end_max = (start + ctx.params.max_length).min(n);
        for end in (start + ctx.params.min_length)..=end_max {
            emit_span(ctx, start, end, &mut out);
        }
    }
    out
}

/// The (N-terminal, C-terminal) context of position `i` of a span of length
/// `n`, for [`model::AminoAcidSet::forms`].
#[inline]
fn position_termini(
    i: usize,
    n: usize,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
) -> (Terminus, Terminus) {
    (
        Terminus::at(i == 0, is_protein_n_term),
        Terminus::at(i + 1 == n, is_protein_c_term),
    )
}

/// What a span's protein annotation says about its residues: the annotated
/// modification sites and, for a variant entry, its single-residue
/// substitutions, with the span's start in protein coordinates. Decides where
/// a mod restricted to annotated sites
/// ([`model::Modification::annotated_site_mass`]) may go, and leaves out a
/// substituted residue whose form weighs the same as a searchable form of the
/// residue it replaced (D, or N deamidated, where the canonical has N): that
/// form is the canonical peptide's, not evidence for the variant.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct SpanSites<'a> {
    sites: &'a [ModSite],
    substitutions: &'a [Substitution],
    start: u32,
    /// The search's residue forms and the span's length and protein-terminal
    /// flags; set by [`SpanSites::with_forms`] before a walk. The forms are
    /// needed for a span with substitutions, the rest also for mods that
    /// block cleavage.
    forms: Option<&'a AminoAcidSet>,
    len: usize,
    protein_n_term: bool,
    protein_c_term: bool,
    /// The search's enzyme, for mods that block cleavage.
    enzyme: Option<Enzyme>,
}

impl<'a> SpanSites<'a> {
    pub(crate) fn new(protein: &'a Protein, span_start: usize) -> Self {
        Self {
            sites: protein.mod_sites(),
            substitutions: protein.substitutions(),
            start: span_start as u32,
            forms: None,
            len: 0,
            protein_n_term: false,
            protein_c_term: false,
            enzyme: None,
        }
    }

    /// These sites, judging substitutions against the residue forms of `aa`
    /// for a span of `len` residues with the given protein-terminal flags: a
    /// substituted position is compared with the forms the replaced residue
    /// takes in the same terminal context. `enzyme` decides where a mod that
    /// blocks cleavage may not sit.
    pub(crate) fn with_forms<'b>(
        self,
        aa: &'b AminoAcidSet,
        len: usize,
        protein_n_term: bool,
        protein_c_term: bool,
        enzyme: Enzyme,
    ) -> SpanSites<'b>
    where
        'a: 'b,
    {
        SpanSites {
            forms: Some(aa),
            len,
            protein_n_term,
            protein_c_term,
            enzyme: Some(enzyme),
            ..self
        }
    }

    /// Whether the span was cleaved at its residue `i` (`residue`): its last
    /// residue for an enzyme that cuts after it, its first for one that cuts
    /// before, unless that end is the protein's own terminus.
    #[inline]
    fn cut_at(&self, i: usize, residue: u8) -> bool {
        let Some(enzyme) = self.enzyme else {
            return false;
        };
        if enzyme.is_n_term() {
            i == 0 && !self.protein_n_term && enzyme.is_cleavable_before(residue)
        } else {
            i + 1 == self.len && !self.protein_c_term && enzyme.is_cleavable_after(residue)
        }
    }

    /// Whether `variant` may sit at position `i` of the span.
    #[inline]
    fn allows(&self, i: usize, variant: &AminoAcid) -> bool {
        let position = self.start + i as u32;
        if !variant
            .mod_
            .is_none_or(|m| m.allowed_at(self.sites, position))
        {
            return false;
        }
        if variant.mod_.is_some_and(|m| m.blocks_cleavage) && self.cut_at(i, variant.residue) {
            return false;
        }
        let Some(aa) = self.forms else {
            return true;
        };
        if self.substitutions.is_empty() {
            return true;
        }
        let Ok(k) = self
            .substitutions
            .binary_search_by_key(&position, |s| s.position)
        else {
            return true;
        };
        let mass = |a: &AminoAcid| a.mass + a.mod_.map_or(0.0, |m| m.mass_delta);
        let here = mass(variant);
        let (nt, ct) = position_termini(i, self.len, self.protein_n_term, self.protein_c_term);
        !aa.forms(self.substitutions[k].original, nt, ct)
            .iter()
            .any(|f| (mass(f) - here).abs() <= model::protein::SITE_MASS_TOLERANCE)
    }
}

/// Generate every combination of variable-mod applications for `span`,
/// up to `params.max_variable_mods_per_peptide` mods total. A mod restricted
/// to annotated sites is only placed where `sites` allows it.
///
/// `is_protein_n_term`: the span begins at the protein's (biological)
/// N-terminus. `is_protein_c_term`: the span ends at the protein's C-terminus.
///
/// Each position takes the forms [`model::AminoAcidSet::forms`] lists for its
/// terminal context: the first and last positions carry the peptide- (or
/// protein-) terminal mods stacked with the residue's side-chain mods;
/// interior positions carry side-chain mods only.
pub(crate) fn expand_mod_combinations(
    span: &[u8],
    params: &SearchParams,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
    sites: SpanSites<'_>,
) -> Vec<Vec<AminoAcid>> {
    let n = span.len();
    let mut position_variants: Vec<&[AminoAcid]> = Vec::with_capacity(n);
    let mut position_slots: Vec<&[u8]> = Vec::with_capacity(n);
    for (i, &r) in span.iter().enumerate() {
        let (nt, ct) = position_termini(i, n, is_protein_n_term, is_protein_c_term);
        position_variants.push(params.aa_set.forms(r, nt, ct));
        position_slots.push(params.aa_set.form_slots(r, nt, ct));
    }

    let mut out = Vec::new();
    let mut current = Vec::with_capacity(n);
    expand_recursive(
        &position_variants,
        &position_slots,
        0,
        &mut current,
        0,
        params.max_variable_mods_per_peptide,
        sites.with_forms(
            &params.aa_set,
            n,
            is_protein_n_term,
            is_protein_c_term,
            params.enzyme,
        ),
        &mut out,
    );
    out
}

/// Per-search mass bounds for [`walk_forms_bounded`]: the (lightest,
/// heaviest) form mass of every `(residue, N context, C context)` entry of
/// [`model::AminoAcidSet::forms`], indexed by [`form_index`];
/// `(f64::MAX, f64::MIN)` for an entry with no form.
pub struct FormWalkTables {
    bounds: Vec<(f64, f64)>,
}

/// Mass of one variant: residue plus its modification delta.
#[inline]
fn variant_mass(aa: &AminoAcid) -> f64 {
    aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta)
}

/// (lightest, heaviest) variant mass of `variants`, or `(f64::MAX, f64::MIN)`
/// when there are none.
fn variant_bounds(variants: &[AminoAcid]) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for v in variants {
        let m = variant_mass(v);
        lo = lo.min(m);
        hi = hi.max(m);
    }
    (lo, hi)
}

impl FormWalkTables {
    pub fn new(params: &SearchParams) -> Self {
        let mut bounds = vec![(f64::MAX, f64::MIN); FORM_TABLE_LEN];
        for r in 0..=255u8 {
            for nt in Terminus::ALL {
                for ct in Terminus::ALL {
                    bounds[form_index(r, nt, ct)] = variant_bounds(params.aa_set.forms(r, nt, ct));
                }
            }
        }
        Self { bounds }
    }

    /// (lightest, heaviest) variant mass at position `i` of a span of length `n`.
    #[inline]
    fn position_bounds(&self, span: &[u8], i: usize, prot_n: bool, prot_c: bool) -> (f64, f64) {
        let (nt, ct) = position_termini(i, span.len(), prot_n, prot_c);
        self.bounds[form_index(span[i], nt, ct)]
    }
}

/// Peptide lengths up to this many residues keep the walk's per-position
/// buffers on the stack.
const WALK_INLINE_LEN: usize = 64;

/// Per-position variant lists of one walk.
type WalkVariants<'t> = smallvec::SmallVec<[&'t [AminoAcid]; WALK_INLINE_LEN]>;
/// Per-position variable-mod counts of each variant of one walk.
type WalkSlots<'t> = smallvec::SmallVec<[&'t [u8]; WALK_INLINE_LEN]>;
/// Suffix sums of the per-position minimum or maximum variant mass.
type WalkRest = smallvec::SmallVec<[f64; WALK_INLINE_LEN]>;

/// Fill `position_variants` / `position_slots` with the forms (and their
/// variable-mod counts) of every position of `span`, and `min_rest` /
/// `max_rest` with the suffix sums of the per-position minimum and maximum
/// variant mass; `false` when the span is empty or a position has no variant.
/// The buffers are the caller's, so nothing is returned by value.
#[allow(clippy::too_many_arguments)]
fn walk_bounds<'t>(
    span: &[u8],
    params: &'t SearchParams,
    tables: &FormWalkTables,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
    position_variants: &mut WalkVariants<'t>,
    position_slots: &mut WalkSlots<'t>,
    min_rest: &mut WalkRest,
    max_rest: &mut WalkRest,
) -> bool {
    let n = span.len();
    if n == 0 {
        return false;
    }
    position_variants.clear();
    position_slots.clear();
    for (i, &r) in span.iter().enumerate() {
        let (nt, ct) = position_termini(i, n, is_protein_n_term, is_protein_c_term);
        position_variants.push(params.aa_set.forms(r, nt, ct));
        position_slots.push(params.aa_set.form_slots(r, nt, ct));
    }
    min_rest.clear();
    min_rest.resize(n + 1, 0.0);
    max_rest.clear();
    max_rest.resize(n + 1, 0.0);
    for i in (0..n).rev() {
        if position_variants[i].is_empty() {
            return false;
        }
        let (lo, hi) = tables.position_bounds(span, i, is_protein_n_term, is_protein_c_term);
        min_rest[i] = min_rest[i + 1] + lo;
        max_rest[i] = max_rest[i + 1] + hi;
    }
    true
}

/// The lightest and heaviest neutral mass (without `H2O`) any form of `span`
/// can take: `min_rest[0]` and `max_rest[0]` of [`walk_bounds`], summed in the
/// same order so the values are identical; `None` exactly when it fails.
fn span_mass_bounds(
    span: &[u8],
    tables: &FormWalkTables,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
) -> Option<(f64, f64)> {
    let n = span.len();
    if n == 0 {
        return None;
    }
    let (mut min_rest, mut max_rest) = (0.0f64, 0.0f64);
    for i in (0..n).rev() {
        let (lo, hi) = tables.position_bounds(span, i, is_protein_n_term, is_protein_c_term);
        if lo > hi {
            return None;
        }
        min_rest += lo;
        max_rest += hi;
    }
    Some((min_rest, max_rest))
}

/// Whether the bounded walk over `rec` with window `[mass_lo, mass_hi]` can
/// visit any form: `false` exactly when [`for_each_record_form_masses_bounded`]
/// would return without visiting one, decided from the same bounds.
pub fn record_may_have_forms_in(
    db: &SearchIndex,
    tables: &FormWalkTables,
    rec: &crate::candidate_index::IndexRecord,
    mass_lo: f64,
    mass_hi: f64,
) -> bool {
    use model::mass::H2O;
    let Some(RecordSpan {
        span,
        n_term,
        c_term,
        ..
    }) = record_span(db, rec)
    else {
        return false;
    };
    match span_mass_bounds(span, tables, n_term, c_term) {
        Some((min_rest, max_rest)) => !(min_rest + H2O > mass_hi || max_rest + H2O < mass_lo),
        None => false,
    }
}

/// Keep only the records of `recs` that [`record_may_have_forms_in`] accepts,
/// in order. The records arrive in mass order, so each one's protein and span
/// are a cache miss; they are touched for a batch of records in tight loops
/// first, letting the misses overlap, and then tested.
pub fn retain_records_with_forms_in(
    db: &SearchIndex,
    tables: &FormWalkTables,
    recs: &mut Vec<crate::candidate_index::IndexRecord>,
    mass_lo: f64,
    mass_hi: f64,
) {
    const BATCH: usize = 32;
    let proteins = &db.db.proteins;
    let mut kept = 0usize;
    let mut start = 0usize;
    while start < recs.len() {
        let end = (start + BATCH).min(recs.len());
        let mut spans: [&[u8]; BATCH] = [&[]; BATCH];
        for (slot, rec) in spans.iter_mut().zip(&recs[start..end]) {
            if let Some(p) = proteins.get(rec.protein_index as usize) {
                *slot = &p.sequence;
            }
        }
        let mut touched = 0u8;
        for (span, rec) in spans.iter().zip(&recs[start..end]) {
            touched ^= span.get(rec.start_offset as usize).copied().unwrap_or(0);
        }
        std::hint::black_box(touched);
        for i in start..end {
            if record_may_have_forms_in(db, tables, &recs[i], mass_lo, mass_hi) {
                recs.swap(kept, i);
                kept += 1;
            }
        }
        start = end;
    }
    recs.truncate(kept);
}

/// Shared bounded walk over the peptidoforms of `span`: the same recursion
/// and iteration order as [`expand_mod_combinations`], but subtrees whose
/// neutral mass cannot land in `[mass_lo, mass_hi]` are pruned. Variable
/// modifications only add mass here, so a branch's running sum plus the
/// minimum (or maximum) of the remaining positions bounds its final mass; the
/// maximum ignores the per-peptide mod cap, which can only keep a branch, never
/// drop a feasible one. `visit` sees the pruned index `k` (dense over the
/// forms that survive), the chosen variants, and their masses. Both the
/// fragment-ion index build and its materialisation use this walk, so `k`
/// means the same thing on both sides. Mods restricted to annotated sites are
/// placed where `sites` allows them, as in [`expand_mod_combinations`]; the
/// bounds still count them everywhere, which only loosens the pruning.
#[allow(clippy::too_many_arguments)]
fn walk_forms_bounded(
    span: &[u8],
    params: &SearchParams,
    tables: &FormWalkTables,
    is_protein_n_term: bool,
    is_protein_c_term: bool,
    sites: SpanSites<'_>,
    mass_lo: f64,
    mass_hi: f64,
    mut visit: impl FnMut(usize, &[&AminoAcid], &[f64]),
) {
    use model::mass::H2O;
    use smallvec::SmallVec;
    let mut position_variants = WalkVariants::new();
    let mut position_slots = WalkSlots::new();
    let mut min_rest = WalkRest::new();
    let mut max_rest = WalkRest::new();
    if !walk_bounds(
        span,
        params,
        tables,
        is_protein_n_term,
        is_protein_c_term,
        &mut position_variants,
        &mut position_slots,
        &mut min_rest,
        &mut max_rest,
    ) {
        return;
    }
    if min_rest[0] + H2O > mass_hi || max_rest[0] + H2O < mass_lo {
        return;
    }
    let n = span.len();
    let mut masses: SmallVec<[f64; WALK_INLINE_LEN]> = SmallVec::with_capacity(n);
    let mut chosen: SmallVec<[&AminoAcid; WALK_INLINE_LEN]> = SmallVec::with_capacity(n);
    let mut k = 0usize;
    #[allow(clippy::too_many_arguments)]
    fn rec<'a>(
        pv: &[&'a [AminoAcid]],
        ps: &[&[u8]],
        pos: usize,
        sum: f64,
        masses: &mut SmallVec<[f64; WALK_INLINE_LEN]>,
        chosen: &mut SmallVec<[&'a AminoAcid; WALK_INLINE_LEN]>,
        mods_used: u32,
        max_mods: u32,
        sites: SpanSites<'_>,
        min_rest: &[f64],
        max_rest: &[f64],
        lo: f64,
        hi: f64,
        k: &mut usize,
        visit: &mut impl FnMut(usize, &[&AminoAcid], &[f64]),
    ) {
        if pos == pv.len() {
            visit(*k, chosen, masses);
            *k += 1;
            return;
        }
        for (variant, &slots) in pv[pos].iter().zip(ps[pos]) {
            let new_mods = mods_used + u32::from(slots);
            if new_mods > max_mods || !sites.allows(pos, variant) {
                continue;
            }
            let m = variant.mass + variant.mod_.as_ref().map_or(0.0, |md| md.mass_delta);
            let next = sum + m;
            // Bounds for everything under this choice.
            if next + min_rest[pos + 1] + H2O > hi || next + max_rest[pos + 1] + H2O < lo {
                continue;
            }
            masses.push(m);
            chosen.push(variant);
            rec(
                pv,
                ps,
                pos + 1,
                next,
                masses,
                chosen,
                new_mods,
                max_mods,
                sites,
                min_rest,
                max_rest,
                lo,
                hi,
                k,
                visit,
            );
            chosen.pop();
            masses.pop();
        }
    }
    rec(
        &position_variants,
        &position_slots,
        0,
        0.0,
        &mut masses,
        &mut chosen,
        0,
        params.max_variable_mods_per_peptide,
        sites.with_forms(
            &params.aa_set,
            n,
            is_protein_n_term,
            is_protein_c_term,
            params.enzyme,
        ),
        &min_rest,
        &max_rest,
        mass_lo,
        mass_hi,
        &mut k,
        &mut visit,
    );
}

/// An index record resolved against its protein.
struct RecordSpan<'a> {
    span: &'a [u8],
    n_term: bool,
    c_term: bool,
    pre: u8,
    post: u8,
    sites: SpanSites<'a>,
}

/// Span, terminal flags, flanks and annotated sites of an index record;
/// `None` when the record does not fit its protein.
fn record_span<'a>(
    db: &'a SearchIndex,
    rec: &crate::candidate_index::IndexRecord,
) -> Option<RecordSpan<'a>> {
    use crate::candidate_index::flags;
    let protein = db.db.proteins.get(rec.protein_index as usize)?;
    let seq = &protein.sequence;
    let abs_start = rec.start_offset as usize;
    let abs_end = abs_start + rec.length as usize;
    if abs_end > seq.len() {
        return None;
    }
    let pre = if abs_start == 0 {
        b'_'
    } else {
        seq[abs_start - 1]
    };
    let post = if abs_end == seq.len() {
        b'-'
    } else {
        seq[abs_end]
    };
    Some(RecordSpan {
        span: &seq[abs_start..abs_end],
        n_term: rec.flags & flags::IS_PROTEIN_N_TERM != 0,
        c_term: rec.flags & flags::IS_PROTEIN_C_TERM != 0,
        pre,
        post,
        sites: SpanSites::new(protein, abs_start),
    })
}

/// Visit the peptidoforms of one index record whose neutral mass lies in
/// `[mass_lo, mass_hi]`, as per-residue masses, with the pruned index `k`.
/// Nothing per form is allocated.
pub fn for_each_record_form_masses_bounded(
    db: &SearchIndex,
    params: &SearchParams,
    tables: &FormWalkTables,
    rec: &crate::candidate_index::IndexRecord,
    mass_lo: f64,
    mass_hi: f64,
    mut f: impl FnMut(usize, &[f64]),
) {
    let Some(RecordSpan {
        span,
        n_term,
        c_term,
        sites,
        ..
    }) = record_span(db, rec)
    else {
        return;
    };
    walk_forms_bounded(
        span,
        params,
        tables,
        n_term,
        c_term,
        sites,
        mass_lo,
        mass_hi,
        |k, _, masses| f(k, masses),
    );
}

/// The peptidoforms of one index record whose neutral mass lies in
/// `[mass_lo, mass_hi]` and whose pruned index (the `k`
/// [`for_each_record_form_masses_bounded`] reports) is in `ks`, materialised
/// as candidates in ascending `k` order. `ks` must be sorted ascending;
/// indices past the record's last in-window form are ignored.
#[allow(clippy::too_many_arguments)]
pub fn expand_base_record_bounded(
    db: &SearchIndex,
    params: &SearchParams,
    tables: &FormWalkTables,
    rec: &crate::candidate_index::IndexRecord,
    mass_lo: f64,
    mass_hi: f64,
    ks: &[u32],
) -> Vec<Candidate> {
    use crate::candidate_index::flags;
    let Some(RecordSpan {
        span,
        n_term,
        c_term,
        pre,
        post,
        sites,
    }) = record_span(db, rec)
    else {
        return Vec::new();
    };
    let is_decoy = rec.flags & flags::IS_DECOY != 0;
    let mut out = Vec::with_capacity(ks.len());
    let mut next = 0usize;
    walk_forms_bounded(
        span,
        params,
        tables,
        n_term,
        c_term,
        sites,
        mass_lo,
        mass_hi,
        |k, chosen, _| {
            if next >= ks.len() || ks[next] as usize != k {
                return;
            }
            next += 1;
            let residues: Vec<AminoAcid> = chosen.iter().map(|aa| (*aa).clone()).collect();
            out.push(Candidate {
                peptide: Peptide::new(residues, pre, post),
                protein_index: rec.protein_index as usize,
                start_offset_in_protein: rec.start_offset as usize,
                is_decoy,
                is_protein_n_term: n_term,
                is_protein_c_term: c_term,
            });
        },
    );
    out
}

/// Every combination of one form per position of `position_variants`
/// carrying at most `max_mods` variable modifications in total, counted by
/// `position_slots` (fixed modifications never count).
#[allow(clippy::too_many_arguments)]
fn expand_recursive(
    position_variants: &[&[AminoAcid]],
    position_slots: &[&[u8]],
    pos: usize,
    current: &mut Vec<AminoAcid>,
    mods_used: u32,
    max_mods: u32,
    sites: SpanSites<'_>,
    out: &mut Vec<Vec<AminoAcid>>,
) {
    if pos == position_variants.len() {
        out.push(current.clone());
        return;
    }
    for (variant, &slots) in position_variants[pos].iter().zip(position_slots[pos]) {
        let new_mods = mods_used + u32::from(slots);
        if new_mods > max_mods || !sites.allows(pos, variant) {
            continue;
        }
        current.push(variant.clone());
        expand_recursive(
            position_variants,
            position_slots,
            pos + 1,
            current,
            new_mods,
            max_mods,
            sites,
            out,
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
pub(crate) fn compute_cleavage_positions(
    seq: &[u8],
    primary: Enzyme,
    extras: &[Enzyme],
) -> Vec<u32> {
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

/// Distinct mass offsets a single position's form can add to the base
/// peptide: for every form list of `params.aa_set`, each form's mass minus the
/// mass of the list's form without variable modifications (the form the
/// base-peptide index folds in), de-duplicated by `to_bits()`. A form that
/// stacks two variable mods, or replaces a fixed terminal mod with a variable
/// one, contributes its whole offset; every such form carries at least one
/// variable mod, so a peptidoform's offset from its base is a sum of at most
/// `max_variable_mods_per_peptide` of these.
fn variable_mod_deltas(params: &SearchParams) -> Vec<f64> {
    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut out: Vec<f64> = Vec::new();
    for (variants, slots) in params.aa_set.form_lists() {
        let Some(base) = variants
            .iter()
            .zip(slots)
            .find(|(_, &s)| s == 0)
            .map(|(aa, _)| variant_mass(aa))
        else {
            continue;
        };
        for (aa, _) in variants.iter().zip(slots).filter(|(_, &s)| s > 0) {
            let d = variant_mass(aa) - base;
            if seen.insert(d.to_bits()) {
                out.push(d);
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
/// Step 3 is why the result is bit-identical to the mass-filtered
/// [`enumerate_candidates`] set: mod placement
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
    // per-spectrum scored-candidate set (tailor histogram / null-pool).
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
        let pre = if abs_start == 0 {
            b'_'
        } else {
            seq[abs_start - 1]
        };
        let post = if abs_end == seq.len() {
            b'-'
        } else {
            seq[abs_end]
        };

        let sites = SpanSites::new(protein, abs_start);
        let mod_combinations =
            expand_mod_combinations(span, params, is_protein_n_term, is_protein_c_term, sites);
        for residues in mod_combinations {
            let peptide = Peptide::new(residues, pre, post);
            // DIRECTIONAL final gate, mirroring `matches_precursor`: lighter →
            // tol_left, heavier → tol_right. Keeps the lazily-materialized set
            // EXACTLY the set `matches_precursor` admits per spectrum.
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
/// gated on the nominal-bucket membership predicate of
/// `candidate_nominal_bounds(..)`.
///
/// A spectrum scores exactly the candidates whose
/// `peptide.nominal_residue_mass()` (= `nominal_from(peptide.mass() − H2O)`)
/// lands in the coarse integer bucket window. That window is LOSSY at bucket
/// boundaries: the f64-milli-Da directional window used by
/// [`lazy_candidates_for_precursor`] is a DIFFERENT set at the edges, so the
/// per-spectrum scored set must be gated on the nominal-bucket predicate, NOT the
/// f64 window.
///
/// ## Algorithm
/// Identical to [`lazy_candidates_for_precursor`] (per-Δ base-record fetch →
/// `expand_mod_combinations`) EXCEPT:
/// - The base-record fetch window is derived from the nominal bounds, converted
///   back to a conservative milli-Da mass window that is WIDE ENOUGH to cover
///   every base mass that could expand to an in-window peptidoform (over-fetch
///   is fine — surplus records are filtered out by the acceptance predicate).
/// - Each expanded peptidoform is accepted iff
///   `nominal_from(peptide.mass() − H2O) ∈ [min_nominal, max_nominal]`,
///   replacing the directional f64 final gate.
///
/// The lazy per-Δ enumeration is retained for bounded-RSS (O(window), not full
/// materialization); only the per-candidate ACCEPTANCE predicate changes.
/// Identity of one base record in the on-disk index:
/// `(protein_index, start_offset, length, flags)`. Identical copies of a span
/// (the N-terminal-Met re-enumeration) share a key; their count is the
/// enumeration multiplicity.
pub type BaseRecordKey = (u32, u32, u16, u16);

/// The [`BaseRecordKey`] of an index record.
pub fn base_record_key(rec: &crate::candidate_index::IndexRecord) -> BaseRecordKey {
    (rec.protein_index, rec.start_offset, rec.length, rec.flags)
}

/// Base records whose peptidoforms MAY land in the nominal-bucket window
/// `[min_nominal, max_nominal]`, over-fetched across every distinct variable-mod
/// mass offset, each with the number of identical copies the index holds for it.
/// Sorted by key, which is `enumerate_candidates` order (protein ASC, start
/// ASC, length ASC, flags ASC). Cheap: index range scans only, no expansion.
pub fn base_records_for_nominal_window(
    mi: &crate::candidate_index::MmapCandidateIndex,
    params: &SearchParams,
    min_nominal: i32,
    max_nominal: i32,
) -> Vec<(crate::candidate_index::IndexRecord, u32)> {
    use crate::candidate_index::IndexRecord;
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
    let bucket_da = 1.0 / scaler;
    let residue_lo = (min_nominal as f64 - 0.5) / scaler - bucket_da;
    let residue_hi = (max_nominal as f64 + 0.5) / scaler + bucket_da;
    let full_lo = residue_lo + H2O;
    let full_hi = residue_hi + H2O;

    // (3) Union of base records across all Δ windows. Base mass =
    //     full_peptide_mass − Δ, so the base-mass fetch window is
    //     `[full_lo − Δ, full_hi − Δ]`. Copies of one span sit at ONE mass, so
    //     every Δ window that reaches the span returns all of its copies; the
    //     copy count is therefore the same in every window that has the key.
    let mut records: rustc_hash::FxHashMap<BaseRecordKey, (IndexRecord, u32)> =
        rustc_hash::FxHashMap::default();
    let mut counts_this_window: rustc_hash::FxHashMap<BaseRecordKey, u32> =
        rustc_hash::FxHashMap::default();
    for (_, delta) in &offsets {
        let lo = ((full_lo - delta) * 1000.0).floor() as i64;
        let hi = ((full_hi - delta) * 1000.0).ceil() as i64;
        if hi < 0 {
            continue;
        }
        let lo_milli = lo.max(0) as u64;
        let hi_milli = hi as u64;
        counts_this_window.clear();
        for rec in mi.mass_window(lo_milli, hi_milli) {
            let key = base_record_key(&rec);
            *counts_this_window.entry(key).or_insert(0) += 1;
            records.entry(key).or_insert((rec, 0));
        }
        for (key, n) in counts_this_window.drain() {
            let slot = records.get_mut(&key).expect("inserted above");
            slot.1 = slot.1.max(n);
        }
    }

    let mut out: Vec<(IndexRecord, u32)> = records.into_values().collect();
    out.sort_unstable_by_key(|(r, _)| base_record_key(r));
    out
}

/// The distinct base records [`base_records_for_nominal_window`] returns for
/// `[min_nominal, max_nominal]` (one per [`BaseRecordKey`], sorted by key),
/// without the copy counts, restricted to those `keep` retains. The union of
/// the per-Δ mass windows is scanned in parallel slices, so one wide window
/// does not run on a single thread; `keep` filters each slice's records in
/// place inside the scan, so only the kept records are sorted.
pub fn distinct_base_records_for_nominal_window(
    mi: &crate::candidate_index::MmapCandidateIndex,
    params: &SearchParams,
    min_nominal: i32,
    max_nominal: i32,
    keep: impl Fn(&mut Vec<crate::candidate_index::IndexRecord>) + Sync,
) -> Vec<crate::candidate_index::IndexRecord> {
    use model::mass::{H2O, INTEGER_MASS_SCALER};
    use rayon::prelude::*;

    if min_nominal > max_nominal {
        return Vec::new();
    }
    // The same per-Δ milli-Da windows as `base_records_for_nominal_window`.
    let scaler = INTEGER_MASS_SCALER as f64;
    let bucket_da = 1.0 / scaler;
    let full_lo = (min_nominal as f64 - 0.5) / scaler - bucket_da + H2O;
    let full_hi = (max_nominal as f64 + 0.5) / scaler + bucket_da + H2O;
    let windows: Vec<(u64, u64)> = distinct_mod_mass_offsets(params)
        .iter()
        .filter_map(|(_, delta)| {
            let lo = ((full_lo - delta) * 1000.0).floor() as i64;
            let hi = ((full_hi - delta) * 1000.0).ceil() as i64;
            (hi >= 0).then_some((lo.max(0) as u64, hi as u64))
        })
        .collect();
    let (Some(span_lo), Some(span_hi)) = (
        windows.iter().map(|w| w.0).min(),
        windows.iter().map(|w| w.1).max(),
    ) else {
        return Vec::new();
    };
    // Disjoint inclusive slices covering [span_lo, span_hi].
    let n_slices = (rayon::current_num_threads() * 4) as u64;
    let step = ((span_hi - span_lo) / n_slices).max(1);
    let slices: Vec<(u64, u64)> = (0..)
        .map(|i| span_lo + i * step)
        .take_while(|&a| a <= span_hi)
        .map(|a| {
            (
                a,
                if a + step > span_hi {
                    span_hi
                } else {
                    a + step - 1
                },
            )
        })
        .collect();
    let mut out: Vec<crate::candidate_index::IndexRecord> = slices
        .par_iter()
        .flat_map_iter(|&(a, b)| {
            let mut recs = mi.mass_window(a, b);
            recs.retain(|rec| {
                windows
                    .iter()
                    .any(|&(lo, hi)| rec.mass_milli >= lo && rec.mass_milli <= hi)
            });
            keep(&mut recs);
            recs
        })
        .collect();
    out.par_sort_unstable_by_key(base_record_key);
    out.dedup_by_key(|r| base_record_key(r));
    out
}

/// Every peptidoform of one base record — the span reconstructed from the
/// protein, variable mods placed exactly as the precursor-window path places
/// them — with NO window filter. Empty when the record does not fit its
/// protein (a corrupt index).
pub fn expand_base_record(
    db: &SearchIndex,
    params: &SearchParams,
    rec: &crate::candidate_index::IndexRecord,
) -> Vec<Candidate> {
    use crate::candidate_index::flags;

    let Some(protein) = db.db.proteins.get(rec.protein_index as usize) else {
        return Vec::new();
    };
    let seq = &protein.sequence;
    let abs_start = rec.start_offset as usize;
    let abs_end = abs_start + rec.length as usize;
    if abs_end > seq.len() {
        return Vec::new();
    }
    let span = &seq[abs_start..abs_end];

    let is_protein_n_term = rec.flags & flags::IS_PROTEIN_N_TERM != 0;
    let is_protein_c_term = rec.flags & flags::IS_PROTEIN_C_TERM != 0;
    let is_decoy = rec.flags & flags::IS_DECOY != 0;

    let pre = if abs_start == 0 {
        b'_'
    } else {
        seq[abs_start - 1]
    };
    let post = if abs_end == seq.len() {
        b'-'
    } else {
        seq[abs_end]
    };

    let sites = SpanSites::new(protein, abs_start);
    expand_mod_combinations(span, params, is_protein_n_term, is_protein_c_term, sites)
        .into_iter()
        .map(|residues| Candidate {
            peptide: Peptide::new(residues, pre, post),
            protein_index: rec.protein_index as usize,
            start_offset_in_protein: abs_start,
            is_decoy,
            is_protein_n_term,
            is_protein_c_term,
        })
        .collect()
}

/// Peptidoforms in the nominal-bucket window `[min_nominal, max_nominal]`, in
/// the order `enumerate_candidates` visits them: base records in
/// `enumerate_candidates` order, each expanded in `expand_mod_combinations`
/// order, and gated on the membership predicate
/// `nominal_from(peptide.mass() − H2O) ∈ [min, max]` (cached as
/// `Peptide::nominal_residue_mass()`).
pub fn lazy_candidates_for_nominal_window(
    mi: &crate::candidate_index::MmapCandidateIndex,
    db: &SearchIndex,
    params: &SearchParams,
    min_nominal: i32,
    max_nominal: i32,
) -> Vec<Candidate> {
    base_records_for_nominal_window(mi, params, min_nominal, max_nominal)
        .into_iter()
        .flat_map(|(rec, _)| expand_base_record(db, params, &rec))
        .filter(|c| {
            let key = c.peptide.nominal_residue_mass();
            key >= min_nominal && key <= max_nominal
        })
        .collect()
}

/// Count enzymatic cleavage sites strictly inside `(start, end)` — the missed
/// cleavages internal to a candidate peptide spanning `[start, end]`.
/// `cleavage_positions` is sorted ascending. This matches the strict path's
/// `offset` count (positions strictly between the two tryptic ends), so the
/// semi-specific path can apply the same `max_missed_cleavages` bound.
fn missed_cleavages_between(cleavage_positions: &[u32], start: u32, end: u32) -> u32 {
    let lo = cleavage_positions.partition_point(|&c| c <= start);
    let hi = cleavage_positions.partition_point(|&c| c < end);
    (hi - lo) as u32
}

#[cfg(test)]
mod missed_cleavage_tests {
    use super::missed_cleavages_between;

    #[test]
    fn counts_only_strictly_internal_sites() {
        let cp = [0u32, 5, 10, 15, 20];
        // endpoints excluded; only sites strictly between are missed cleavages
        assert_eq!(missed_cleavages_between(&cp, 0, 20), 3); // 5,10,15
        assert_eq!(missed_cleavages_between(&cp, 5, 10), 0); // adjacent
        assert_eq!(missed_cleavages_between(&cp, 0, 10), 1); // 5
        assert_eq!(missed_cleavages_between(&cp, 3, 17), 3); // 5,10,15
        assert_eq!(missed_cleavages_between(&cp, 10, 11), 0);
    }

    #[test]
    fn matches_strict_offset_semantics() {
        // For strict spans cp[i]..cp[j], missed == j-i-1 (the strict path's offset).
        let cp = [2u32, 9, 14, 30];
        assert_eq!(missed_cleavages_between(&cp, cp[0], cp[1]), 0); // offset 0
        assert_eq!(missed_cleavages_between(&cp, cp[0], cp[2]), 1); // offset 1
        assert_eq!(missed_cleavages_between(&cp, cp[0], cp[3]), 2); // offset 2
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn decoy_membership_requires_underscore_delimiter() {
        use crate::decoy::is_decoy_accession;
        // Decoys are built as "<prefix>_<orig>" (decoy_db), so
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
        assert_eq!(
            both,
            vec![0, 2, 4, 6, 8],
            "GluC+Trypsin = union of cut sites"
        );
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
        for e in [
            Enzyme::Trypsin,
            Enzyme::LysC,
            Enzyme::GluC,
            Enzyme::AspN,
            Enzyme::NonSpecific,
            Enzyme::NoCleavage,
        ] {
            let single = compute_cleavage_positions(seq, e, &[]);
            // NoCleavage extras contribute no sites; result must equal the primary alone.
            let with_nocleave = compute_cleavage_positions(seq, e, &[Enzyme::NoCleavage]);
            assert_eq!(
                single, with_nocleave,
                "{e:?}: NoCleavage extra adds no sites"
            );
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
    // mass-filtered subset of the full `enumerate_candidates`.

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
                c.peptide
                    .residues
                    .iter()
                    .map(|aa| aa.residue)
                    .eq(seq.iter().copied())
            })
            .unwrap_or_else(|| {
                panic!(
                    "no unmodified candidate spelling {:?}",
                    std::str::from_utf8(seq)
                )
            })
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
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
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
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
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
                annotation: None,
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

    // A FIXED peptide-N-term mod (e.g. the TMT/iTRAQ tag) must be applied at
    // EVERY pos-0 residue, INCLUDING protein-N-term peptides; otherwise a
    // protein-N-term peptide carries the wrong precursor/fragment mass.
    #[test]
    fn fixed_n_term_tag_applies_to_protein_n_term_peptides() {
        use model::modification::ResidueSpec;
        let tmt = Modification {
            name: "TMT6plex".into(),
            mass_delta: 229.162932,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::NTerm,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
        };
        let acetyl = Modification {
            name: "Acetyl".into(),
            mass_delta: 42.010565,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::ProtNTerm,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
        };
        // First tryptic peptide AAAPEPTIDEK sits at the protein N-terminus (offset 0).
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "tmt protein-nterm".into(),
                sequence: b"AAAPEPTIDEKQQQR".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(tmt)
            .add_variable_mod(acetyl)
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 1;

        let pnt: Vec<_> = enumerate_candidates(&db, &params, "XXX")
            .filter(|c| !c.is_decoy && c.start_offset_in_protein == 0)
            .collect();
        assert!(
            !pnt.is_empty(),
            "expected a protein-N-term peptide candidate"
        );
        let has_tmt = pnt.iter().any(|c| {
            c.peptide
                .residues
                .first()
                .and_then(|aa| aa.mod_.as_ref())
                .map(|m| (m.mass_delta - 229.162932).abs() < 1e-3)
                .unwrap_or(false)
        });
        assert!(
            has_tmt,
            "protein-N-term peptide must carry the fixed N-term TMT tag (finding 2.3)"
        );
    }

    /// A substituted residue at a peptide N-terminus is judged against the
    /// forms the replaced residue takes there: with a fixed N-term tag and
    /// deamidation searched, D+tag at the start of a peptide weighs the same as
    /// the canonical N deamidated + tag, so the variant form is left out.
    /// Acetyl-K that blocks cleavage: never on the K a peptide was cleaved
    /// after, but on an internal K (a missed cleavage) and on the protein's
    /// own C-terminal K.
    #[test]
    fn a_mod_that_blocks_cleavage_never_sits_on_the_cleaved_residue() {
        let seq = b"AAAAKGGGGKCCCCK".to_vec();
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: String::new(),
                sequence: seq.clone(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let run = |blocks: bool| {
            let acetyl = Modification {
                name: "Acetyl".into(),
                mass_delta: 42.010565,
                residue: ResidueSpec::Specific(b'K'),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
                annotated_site_mass: None,
                blocks_cleavage: blocks,
                diagnostic_ions: Vec::new(),
            };
            let aa = AminoAcidSetBuilder::new_standard()
                .add_variable_mod(acetyl)
                .build()
                .unwrap();
            let mut params = SearchParams::default_tryptic(aa);
            params.min_length = 3;
            params.max_variable_mods_per_peptide = 2;
            // (start, length, modified positions) of every target candidate.
            enumerate_candidates(&db, &params, "XXX")
                .filter(|c| !c.is_decoy)
                .map(|c| {
                    let mods: Vec<usize> = c
                        .peptide
                        .residues
                        .iter()
                        .enumerate()
                        .filter(|(_, a)| a.mod_.is_some())
                        .map(|(i, _)| i)
                        .collect();
                    (c.start_offset_in_protein, c.peptide.length(), mods)
                })
                .collect::<Vec<_>>()
        };
        let blocked = run(true);
        for (start, len, mods) in &blocked {
            let protein_end = start + len == seq.len();
            assert!(
                protein_end || !mods.contains(&(len - 1)),
                "acetyl on the cleaved K of {start}+{len}"
            );
        }
        assert!(
            blocked
                .iter()
                .any(|(s, l, m)| *s == 0 && *l == 10 && m == &vec![4]),
            "acetyl on the internal K of a missed cleavage"
        );
        assert!(
            blocked
                .iter()
                .any(|(s, l, m)| *s == 10 && *l == 5 && m == &vec![4]),
            "acetyl on the protein's own C-terminal K"
        );
        // Without the attribute the cleaved K may carry it.
        assert!(run(false)
            .iter()
            .any(|(s, l, m)| *s == 0 && *l == 5 && m == &vec![4]));
    }

    #[test]
    fn substitution_at_a_peptide_terminus_uses_the_terminal_forms() {
        use model::protein::{ProteinAnnotation, Substitution, VariantRegion};
        let annotation = ProteinAnnotation {
            variant: vec![VariantRegion { start: 4, end: 5 }],
            ..Default::default()
        }
        .with_substitutions(vec![Substitution {
            position: 4,
            original: b'N',
        }]);
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "variant fixture".into(),
                sequence: b"PEPKDAGTRK".to_vec(),
                annotation: Some(std::sync::Arc::new(annotation)),
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let has_dagtr = |deamidation: bool| {
            let mut b = AminoAcidSetBuilder::new_standard().add_fixed_mod(Modification {
                name: "TMT6plex".into(),
                mass_delta: 229.162932,
                residue: ResidueSpec::Wildcard,
                location: ModLocation::NTerm,
                fixed: true,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
                annotated_site_mass: None,
                blocks_cleavage: false,
                diagnostic_ions: Vec::new(),
            });
            if deamidation {
                b = b.add_variable_mod(deamidation_n());
            }
            let mut params = SearchParams::default_tryptic(b.build().unwrap());
            params.min_length = 3;
            params.max_variable_mods_per_peptide = 1;
            let found = enumerate_candidates(&db, &params, "XXX")
                .any(|c| !c.is_decoy && c.start_offset_in_protein == 4 && c.peptide.length() == 5);
            found
        };
        assert!(
            has_dagtr(false),
            "without deamidation D is evidence for the variant"
        );
        assert!(
            !has_dagtr(true),
            "D+TMT is the canonical N deamidated + TMT"
        );
    }

    /// A mature chain starting with M gets no Met-cleaved walk, so the full
    /// walk keeps its spans that start at the second residue: a C-terminal
    /// chain's span from there to its end is enumerated once, not as a
    /// protein N-terminal peptide.
    #[test]
    fn mature_chain_starting_with_met_keeps_its_second_residue_spans() {
        use model::protein::{ChainEnds, ProteinAnnotation};
        let annotation = ProteinAnnotation {
            chain: Some(ChainEnds {
                n_term: false,
                c_term: true,
            }),
            ..Default::default()
        };
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "chain fixture".into(),
                sequence: b"MAKAGERLLK".to_vec(),
                annotation: Some(std::sync::Arc::new(annotation)),
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 2;
        params.num_tolerable_termini = 0;
        let starting_at_1: Vec<_> = enumerate_candidates(&db, &params, "XXX")
            .filter(|c| !c.is_decoy && c.start_offset_in_protein == 1)
            .collect();
        assert_eq!(starting_at_1.len(), 1, "AKAGERLLK once");
        assert_eq!(starting_at_1[0].peptide.residues.len(), 9);
        assert!(!starting_at_1[0].is_protein_n_term);
    }

    /// The Met-cleaved walk adds only the peptides starting at the post-Met
    /// residue: every other span of an M-initial protein is enumerated once.
    #[test]
    fn met_cleaved_walk_adds_only_post_met_n_term_peptides() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "met fixture".into(),
                sequence: b"MAKAGERLLKPEPTIDER".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 2;
        params.max_missed_cleavages = 1;
        for ntt in [2u8, 1, 0] {
            params.num_tolerable_termini = ntt;
            let cands: Vec<_> = enumerate_candidates(&db, &params, "XXX")
                .filter(|c| !c.is_decoy)
                .collect();
            let mut keys: Vec<(usize, usize)> = cands
                .iter()
                .map(|c| (c.start_offset_in_protein, c.peptide.residues.len()))
                .collect();
            let n = keys.len();
            keys.sort_unstable();
            keys.dedup();
            assert_eq!(keys.len(), n, "ntt={ntt}: a span is enumerated twice");
            assert!(
                cands
                    .iter()
                    .any(|c| c.start_offset_in_protein == 1 && c.is_protein_n_term),
                "ntt={ntt}: missing the Met-cleaved N-terminal peptide"
            );
            assert!(
                cands
                    .iter()
                    .all(|c| c.is_protein_n_term == (c.start_offset_in_protein <= 1)),
                "ntt={ntt}: only spans starting at 0 or 1 are protein N-terminal"
            );
        }
    }

    /// TMT fixed on K and the N-terminus, CAM fixed on C, Ox-M and protein
    /// N-term Acetyl variable: the N-terminal residue carries its side-chain
    /// mod stacked with the terminal one, and the lazy out-of-core path (whose
    /// Δ offsets include Acetyl replacing the N-terminal tag) matches the
    /// in-RAM candidate set.
    #[test]
    fn stacked_terminal_forms_enumerate_and_match_lazy() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "229.162932,K,fix,any,TMT6plex\n\
             229.162932,*,fix,N-term,TMT6plex\n\
             57.021464,C,fix,any,Carbamidomethyl\n\
             15.994915,M,opt,any,Oxidation\n\
             42.010565,*,opt,Prot-N-term,Acetyl\n",
        )
        .unwrap();
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_mods_from_file(tmp.path())
            .unwrap()
            .build()
            .unwrap();
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "tmt stacking fixture".into(),
                sequence: b"MCKAGERKLLRMAEPK".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 2;
        params.max_variable_mods_per_peptide = 2;
        let all: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
        let first_delta = |c: &Candidate| {
            c.peptide.residues[0]
                .mod_
                .as_ref()
                .map_or(0.0, |m| m.mass_delta)
        };
        let has = |start: usize, len: usize, delta: f64| {
            all.iter().any(|c| {
                !c.is_decoy
                    && c.start_offset_in_protein == start
                    && c.peptide.residues.len() == len
                    && (first_delta(c) - delta).abs() < 1e-3
            })
        };
        // MCK at the protein N-terminus: M with Ox + Ac (two variable mods).
        assert!(has(0, 3, 15.994915 + 42.010565));
        // CK after Met cleavage: C keeps CAM under the TMT tag or under Acetyl.
        assert!(has(1, 2, 57.021464 + 229.162932));
        assert!(has(1, 2, 57.021464 + 42.010565));
        assert!(!has(1, 2, 42.010565));
        // KLLR: an N-terminal K carries both TMT tags.
        assert!(has(7, 4, 2.0 * 229.162932));
        // MAEPK inside the protein: N-terminal M is tagged, optionally oxidised.
        assert!(has(11, 5, 229.162932 + 15.994915));

        let (_tmp, mi) = build_and_open(&db, &params);
        let mut masses: Vec<f64> = all.iter().map(|c| c.peptide.mass()).collect();
        masses.sort_by(|a, b| a.partial_cmp(b).unwrap());
        masses.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        for &m in &masses {
            assert_lazy_matches_inram(&db, &params, &mi, m, 0.01);
        }
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
            lazy,
            inram,
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
            isolation_target_mz: None,
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
        params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Da(0.5));
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
                    mmap_nominal,
                    ram,
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
                    &mi,
                    &db,
                    &params,
                    neutral_obs,
                    tol_left,
                    tol_right,
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
                    && c.peptide
                        .residues
                        .iter()
                        .map(|aa| aa.residue)
                        .eq(b"PEPTMK".iter().copied())
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
                    && c.peptide
                        .residues
                        .iter()
                        .map(|aa| aa.residue)
                        .eq(b"PEPTMK".iter().copied())
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
                annotation: None,
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
