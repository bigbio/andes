pub mod records;

use crate::decoy::is_decoy_accession_affix;
use crate::index::SearchIndex;
use crate::params::SearchParams;
use model::aa_set::label_code::LabelState;
use model::aa_set::{form_index, label_code, AminoAcidSet, Terminus, FORM_TABLE_LEN};
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
    pub is_protein_n_term: bool,
    pub is_protein_c_term: bool,
}

pub fn enumerate_candidates<'a>(
    idx: &'a SearchIndex,
    params: &'a SearchParams,
    decoy_prefix: &'a str,
) -> impl Iterator<Item = Candidate> + 'a {
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

pub fn candidate_has_nxst_sequon(c: &Candidate) -> bool {
    let res: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
    span_has_nxst_sequon(&res, c.peptide.post)
}

pub fn span_has_nxst_sequon(residues: &[u8], post: u8) -> bool {
    use glyco::sequon::{boundary_nxst_site, has_nxst_sequon};
    has_nxst_sequon(residues) || boundary_nxst_site(residues, post).is_some()
}

pub fn candidate_has_ser_or_thr(c: &Candidate) -> bool {
    use glyco::sequon::has_ser_or_thr;
    let res: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
    has_ser_or_thr(&res)
}

pub fn span_has_ser_or_thr(residues: &[u8], _post: u8) -> bool {
    use glyco::sequon::has_ser_or_thr;
    has_ser_or_thr(residues)
}

fn enumerate_protein(
    protein: &Protein,
    protein_index: usize,
    is_decoy: bool,
    params: &SearchParams,
) -> Vec<Candidate> {
    let mut out = enumerate_protein_from_offset(protein, 0, protein_index, is_decoy, params);

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

#[inline]
pub(crate) fn has_met_cleaved_walk(protein: &Protein) -> bool {
    let seq = &protein.sequence;
    seq.len() > 1 && seq[0] == b'M' && protein.chain_ends().is_none()
}

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

    let window = match protein.variant_region() {
        Some(v) => {
            let reach = params.max_length as usize;
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
    let n_term_only = seq_offset > 0;

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
    let anchors = &cleavage_positions[cleavage_positions.partition_point(|&c| c < window.0)
        ..cleavage_positions.partition_point(|&c| c <= window.1)];

    let mut out = Vec::new();

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

    if ntt == 1 {
        for &start in anchors {
            if n_term_only && start != 0 {
                break;
            }
            let c_min = start + params.min_length;
            let c_max = (start + params.max_length).min(n);
            for end in c_min..=c_max {
                if cleavage_set.contains(&end) {
                    continue;
                }
                if missed_cleavages_between(&cleavage_positions, start, end)
                    > params.max_missed_cleavages
                {
                    continue;
                }
                emit_span(&ctx, start, end, &mut out);
            }
        }

        for &end in anchors.iter().filter(|_| !n_term_only) {
            if end < params.min_length {
                continue;
            }
            let s_min = end.saturating_sub(params.max_length);
            let s_max = end - params.min_length;
            for start in s_min..=s_max {
                if cleavage_set.contains(&start) {
                    continue;
                }
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

struct EmitCtx<'a> {
    sub_seq: &'a [u8],
    seq: &'a [u8],
    seq_offset: usize,
    protein_index: usize,
    is_decoy: bool,
    params: &'a SearchParams,
    variant: &'a [VariantRegion],
    chain: Option<ChainEnds>,
    protein: &'a Protein,
    met_cleaved_walk: bool,
}

#[inline]
fn emit_span(ctx: &EmitCtx<'_>, start: u32, end: u32, out: &mut Vec<Candidate>) {
    let abs_start = start as usize + ctx.seq_offset;
    let abs_end = end as usize + ctx.seq_offset;
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
    if span.iter().any(|&r| AminoAcid::standard(r).is_none()) {
        return;
    }
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

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct SpanSites<'a> {
    sites: &'a [ModSite],
    substitutions: &'a [Substitution],
    start: u32,
    forms: Option<&'a AminoAcidSet>,
    len: usize,
    protein_n_term: bool,
    protein_c_term: bool,
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
    let mut position_labels: Vec<Option<&[u8]>> = Vec::new();
    let mut position_label_slots: Vec<&[u8]> = Vec::new();
    let labelled = params.aa_set.has_label_channels();
    for (i, &r) in span.iter().enumerate() {
        let (nt, ct) = position_termini(i, n, is_protein_n_term, is_protein_c_term);
        position_variants.push(params.aa_set.forms(r, nt, ct));
        position_slots.push(params.aa_set.form_slots(r, nt, ct));
        if labelled {
            position_labels.push(params.aa_set.form_labels(r, nt, ct));
            position_label_slots.push(params.aa_set.form_label_slots(r, nt, ct));
        }
    }

    let mut out = Vec::new();
    let mut current = Vec::with_capacity(n);
    expand_recursive(
        &position_variants,
        &position_slots,
        &position_labels,
        &position_label_slots,
        LabelState::default(),
        0,
        &mut current,
        ModBudget::default(),
        ModBudget::max_of(params),
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

pub struct FormWalkTables {
    bounds: Vec<(f64, f64)>,
}

#[inline]
fn variant_mass(aa: &AminoAcid) -> f64 {
    aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta)
}

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

    #[inline]
    fn position_bounds(&self, span: &[u8], i: usize, prot_n: bool, prot_c: bool) -> (f64, f64) {
        let (nt, ct) = position_termini(i, span.len(), prot_n, prot_c);
        self.bounds[form_index(span[i], nt, ct)]
    }
}

const WALK_INLINE_LEN: usize = 64;

type WalkVariants<'t> = smallvec::SmallVec<[&'t [AminoAcid]; WALK_INLINE_LEN]>;
type WalkSlots<'t> = smallvec::SmallVec<[&'t [u8]; WALK_INLINE_LEN]>;
type WalkRest = smallvec::SmallVec<[f64; WALK_INLINE_LEN]>;

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

pub fn record_may_have_forms_in(
    db: &SearchIndex,
    tables: &FormWalkTables,
    rec: &crate::candidates::records::IndexRecord,
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

pub fn retain_records_with_forms_in(
    db: &SearchIndex,
    tables: &FormWalkTables,
    recs: &mut Vec<crate::candidates::records::IndexRecord>,
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
        pl: &[Option<&[u8]>],
        pls: &[&[u8]],
        label: LabelState,
        pos: usize,
        sum: f64,
        masses: &mut SmallVec<[f64; WALK_INLINE_LEN]>,
        chosen: &mut SmallVec<[&'a AminoAcid; WALK_INLINE_LEN]>,
        mods_used: ModBudget,
        max_mods: ModBudget,
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
        for (j, (variant, &slots)) in pv[pos].iter().zip(ps[pos]).enumerate() {
            let Some(new_mods) = mods_used.take(slots, pls, pos, j, max_mods) else {
                continue;
            };
            if !sites.allows(pos, variant) {
                continue;
            }
            let Some(next_label) = step_label(pl, pos, j, label) else {
                continue;
            };
            let m = variant.mass + variant.mod_.as_ref().map_or(0.0, |md| md.mass_delta);
            let next = sum + m;
            if next + min_rest[pos + 1] + H2O > hi || next + max_rest[pos + 1] + H2O < lo {
                continue;
            }
            masses.push(m);
            chosen.push(variant);
            rec(
                pv,
                ps,
                pl,
                pls,
                next_label,
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
    let mut position_labels: smallvec::SmallVec<[Option<&[u8]>; WALK_INLINE_LEN]> =
        smallvec::SmallVec::new();
    let mut position_label_slots: smallvec::SmallVec<[&[u8]; WALK_INLINE_LEN]> =
        smallvec::SmallVec::new();
    if params.aa_set.has_label_channels() {
        for (i, &r) in span.iter().enumerate() {
            let (nt, ct) = position_termini(i, n, is_protein_n_term, is_protein_c_term);
            position_labels.push(params.aa_set.form_labels(r, nt, ct));
            position_label_slots.push(params.aa_set.form_label_slots(r, nt, ct));
        }
    }
    rec(
        &position_variants,
        &position_slots,
        &position_labels,
        &position_label_slots,
        LabelState::default(),
        0,
        0.0,
        &mut masses,
        &mut chosen,
        ModBudget::default(),
        ModBudget::max_of(params),
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

struct RecordSpan<'a> {
    span: &'a [u8],
    n_term: bool,
    c_term: bool,
    pre: u8,
    post: u8,
    sites: SpanSites<'a>,
}

fn record_span<'a>(
    db: &'a SearchIndex,
    rec: &crate::candidates::records::IndexRecord,
) -> Option<RecordSpan<'a>> {
    use crate::candidates::records::flags;
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

pub fn for_each_record_form_masses_bounded(
    db: &SearchIndex,
    params: &SearchParams,
    tables: &FormWalkTables,
    rec: &crate::candidates::records::IndexRecord,
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

#[allow(clippy::too_many_arguments)]
pub fn expand_base_record_bounded(
    db: &SearchIndex,
    params: &SearchParams,
    tables: &FormWalkTables,
    rec: &crate::candidates::records::IndexRecord,
    mass_lo: f64,
    mass_hi: f64,
    ks: &[u32],
) -> Vec<Candidate> {
    use crate::candidates::records::flags;
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

#[allow(clippy::too_many_arguments)]
fn expand_recursive(
    position_variants: &[&[AminoAcid]],
    position_slots: &[&[u8]],
    position_labels: &[Option<&[u8]>],
    position_label_slots: &[&[u8]],
    label: LabelState,
    pos: usize,
    current: &mut Vec<AminoAcid>,
    mods_used: ModBudget,
    max_mods: ModBudget,
    sites: SpanSites<'_>,
    out: &mut Vec<Vec<AminoAcid>>,
) {
    if pos == position_variants.len() {
        out.push(current.clone());
        return;
    }
    for (j, (variant, &slots)) in position_variants[pos]
        .iter()
        .zip(position_slots[pos])
        .enumerate()
    {
        let Some(new_mods) = mods_used.take(slots, position_label_slots, pos, j, max_mods) else {
            continue;
        };
        if !sites.allows(pos, variant) {
            continue;
        }
        let Some(next_label) = step_label(position_labels, pos, j, label) else {
            continue;
        };
        current.push(variant.clone());
        expand_recursive(
            position_variants,
            position_slots,
            position_labels,
            position_label_slots,
            next_label,
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

#[derive(Debug, Clone, Copy, Default)]
struct ModBudget {
    mods: u32,
    labels: u32,
}

impl ModBudget {
    fn max_of(params: &SearchParams) -> Self {
        Self {
            mods: params.max_variable_mods_per_peptide,
            labels: params.max_label_sites,
        }
    }

    #[inline]
    fn take(
        self,
        slots: u8,
        label_slots: &[&[u8]],
        pos: usize,
        j: usize,
        max: Self,
    ) -> Option<Self> {
        let labels = label_slots
            .get(pos)
            .and_then(|l| l.get(j))
            .map_or(0, |&l| u32::from(l));
        let next = Self {
            mods: self.mods + u32::from(slots) - labels,
            labels: self.labels + labels,
        };
        (next.mods <= max.mods && next.labels <= max.labels).then_some(next)
    }
}

#[inline]
fn step_label(
    labels: &[Option<&[u8]>],
    pos: usize,
    j: usize,
    state: LabelState,
) -> Option<LabelState> {
    match labels.get(pos) {
        Some(Some(codes)) => state.step(codes[j]),
        _ => Some(state),
    }
}

pub(crate) fn compute_cleavage_positions(
    seq: &[u8],
    primary: Enzyme,
    extras: &[Enzyme],
) -> Vec<u32> {
    let n = seq.len() as u32;

    let all = || std::iter::once(primary).chain(extras.iter().copied());

    if all().any(|e| matches!(e, Enzyme::NonSpecific)) {
        return (0..=n).collect();
    }
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

#[derive(Debug, Clone, Copy, PartialEq)]
struct FormDelta {
    delta: f64,
    mods: u32,
    labels: u32,
    channel: u8,
}

fn variable_mod_deltas(params: &SearchParams) -> Vec<FormDelta> {
    let mut seen: std::collections::HashSet<(u64, u32, u32, u8)> = std::collections::HashSet::new();
    let mut out: Vec<FormDelta> = Vec::new();
    for (variants, slots, label_slots, codes) in params.aa_set.form_lists_with_labels() {
        let Some(base) = variants
            .iter()
            .zip(slots)
            .find(|(_, &s)| s == 0)
            .map(|(aa, _)| variant_mass(aa))
        else {
            continue;
        };
        for (j, (aa, &s)) in variants.iter().zip(slots).enumerate() {
            if s == 0 {
                continue;
            }
            let labels = u32::from(label_slots.get(j).copied().unwrap_or(0));
            let channel = if labels == 0 {
                0
            } else {
                match codes.get(j).copied() {
                    Some(c) if (1..=label_code::MAX_CHANNEL).contains(&c) => c,
                    _ => continue,
                }
            };
            let f = FormDelta {
                delta: variant_mass(aa) - base,
                mods: u32::from(s) - labels,
                labels,
                channel,
            };
            if seen.insert((f.delta.to_bits(), f.mods, f.labels, f.channel)) {
                out.push(f);
            }
        }
    }
    out
}

fn distinct_mod_mass_offsets(params: &SearchParams) -> std::sync::Arc<Vec<(u64, f64)>> {
    use std::collections::HashMap;
    use std::sync::{Arc, Mutex, OnceLock};
    type Key = (Vec<(u64, u32, u32, u8)>, u32, u32);
    type Offsets = Arc<Vec<(u64, f64)>>;
    static CACHE: OnceLock<Mutex<HashMap<Key, Offsets>>> = OnceLock::new();
    let deltas = variable_mod_deltas(params);
    let key: Key = (
        deltas
            .iter()
            .map(|f| (f.delta.to_bits(), f.mods, f.labels, f.channel))
            .collect(),
        params.max_variable_mods_per_peptide,
        params.max_label_sites,
    );
    let cache = CACHE.get_or_init(Default::default);
    if let Some(hit) = cache.lock().ok().and_then(|c| c.get(&key).cloned()) {
        return hit;
    }
    let offsets = Arc::new(enumerate_mod_mass_offsets(&deltas, key.1, key.2));
    if let Ok(mut c) = cache.lock() {
        if c.len() >= 16 {
            c.clear();
        }
        c.insert(key, Arc::clone(&offsets));
    }
    offsets
}

fn enumerate_mod_mass_offsets(
    deltas: &[FormDelta],
    max_mods: u32,
    max_labels: u32,
) -> Vec<(u64, f64)> {
    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut out: Vec<(u64, f64)> = Vec::new();
    seen.insert(milli_of(0.0));
    out.push((milli_of(0.0), 0.0));

    type State = (usize, u64, u32, u32);
    fn walk(
        items: &[FormDelta],
        start: usize,
        sum: f64,
        budget: (u32, u32),
        visited: &mut std::collections::HashSet<State>,
        seen: &mut std::collections::HashSet<u64>,
        out: &mut Vec<(u64, f64)>,
    ) {
        let state = (start, sum.to_bits(), budget.0, budget.1);
        if !visited.insert(state) {
            return;
        }
        for (i, f) in items.iter().enumerate().skip(start) {
            if f.mods > budget.0 || f.labels > budget.1 {
                continue;
            }
            let s = sum + f.delta;
            let key = milli_of(s);
            if seen.insert(key) {
                out.push((key, s));
            }
            walk(
                items,
                i,
                s,
                (budget.0 - f.mods, budget.1 - f.labels),
                visited,
                seen,
                out,
            );
        }
    }

    let budget = (max_mods, max_labels);
    let mut channels: Vec<u8> = deltas.iter().map(|f| f.channel).collect();
    channels.push(0);
    channels.sort_unstable();
    channels.dedup();
    for c in channels {
        let items: Vec<FormDelta> = deltas
            .iter()
            .filter(|f| f.channel == 0 || f.channel == c)
            .copied()
            .collect();
        let mut visited = std::collections::HashSet::new();
        walk(&items, 0, 0.0, budget, &mut visited, &mut seen, &mut out);
    }
    out
}

#[inline]
fn milli_of(mass: f64) -> u64 {
    (mass * 1000.0).round() as i64 as u64
}

pub fn lazy_candidates_for_precursor(
    mi: &crate::candidates::records::MmapCandidateIndex,
    db: &SearchIndex,
    params: &SearchParams,
    precursor_mass: f64,
    tol_left_da: f64,
    tol_right_da: f64,
) -> Vec<Candidate> {
    use crate::candidates::records::{flags, IndexRecord};

    let offsets = distinct_mod_mass_offsets(params);

    let tol_left_milli = (tol_left_da * 1000.0).round() as u64;
    let tol_right_milli = (tol_right_da * 1000.0).round() as u64;
    let mut seen_records: std::collections::HashSet<(u32, u32, u16, u16)> =
        std::collections::HashSet::new();
    let mut base_records: Vec<IndexRecord> = Vec::new();
    for (_, delta) in offsets.iter() {
        let center = precursor_mass - delta;
        let lo = (center * 1000.0).round() as i64 - tol_left_milli as i64;
        let hi = (center * 1000.0).round() as i64 + tol_right_milli as i64;
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

pub type BaseRecordKey = (u32, u32, u16, u16);

pub fn base_record_key(rec: &crate::candidates::records::IndexRecord) -> BaseRecordKey {
    (rec.protein_index, rec.start_offset, rec.length, rec.flags)
}

pub fn base_records_for_nominal_window(
    mi: &crate::candidates::records::MmapCandidateIndex,
    params: &SearchParams,
    min_nominal: i32,
    max_nominal: i32,
) -> Vec<(crate::candidates::records::IndexRecord, u32)> {
    use crate::candidates::records::IndexRecord;
    use model::mass::{H2O, INTEGER_MASS_SCALER};

    if min_nominal > max_nominal {
        return Vec::new();
    }

    let offsets = distinct_mod_mass_offsets(params);

    let scaler = INTEGER_MASS_SCALER as f64;
    let bucket_da = 1.0 / scaler;
    let residue_lo = (min_nominal as f64 - 0.5) / scaler - bucket_da;
    let residue_hi = (max_nominal as f64 + 0.5) / scaler + bucket_da;
    let full_lo = residue_lo + H2O;
    let full_hi = residue_hi + H2O;

    let mut records: rustc_hash::FxHashMap<BaseRecordKey, (IndexRecord, u32)> =
        rustc_hash::FxHashMap::default();
    let mut counts_this_window: rustc_hash::FxHashMap<BaseRecordKey, u32> =
        rustc_hash::FxHashMap::default();
    for (_, delta) in offsets.iter() {
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

pub fn distinct_base_records_for_nominal_window(
    mi: &crate::candidates::records::MmapCandidateIndex,
    params: &SearchParams,
    min_nominal: i32,
    max_nominal: i32,
    keep: impl Fn(&mut Vec<crate::candidates::records::IndexRecord>) + Sync,
) -> Vec<crate::candidates::records::IndexRecord> {
    use model::mass::{H2O, INTEGER_MASS_SCALER};
    use rayon::prelude::*;

    if min_nominal > max_nominal {
        return Vec::new();
    }
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
    let mut out: Vec<crate::candidates::records::IndexRecord> = slices
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

pub fn expand_base_record(
    db: &SearchIndex,
    params: &SearchParams,
    rec: &crate::candidates::records::IndexRecord,
) -> Vec<Candidate> {
    use crate::candidates::records::flags;

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

pub fn lazy_candidates_for_nominal_window(
    mi: &crate::candidates::records::MmapCandidateIndex,
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

pub(crate) fn missed_cleavages_between(cleavage_positions: &[u32], start: u32, end: u32) -> u32 {
    let lo = cleavage_positions.partition_point(|&c| c <= start);
    let hi = cleavage_positions.partition_point(|&c| c < end);
    (hi - lo) as u32
}

#[cfg(test)]
#[path = "candidates/tests/missed_cleavage.rs"]
mod missed_cleavage_tests;

#[cfg(test)]
#[path = "candidates/tests/candidates.rs"]
mod tests;
