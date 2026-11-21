use crate::decoy::is_decoy_accession_affix;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
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
    mi: &crate::candidate_index::MmapCandidateIndex,
    db: &SearchIndex,
    params: &SearchParams,
    precursor_mass: f64,
    tol_left_da: f64,
    tol_right_da: f64,
) -> Vec<Candidate> {
    use crate::candidate_index::{flags, IndexRecord};

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

pub fn base_record_key(rec: &crate::candidate_index::IndexRecord) -> BaseRecordKey {
    (rec.protein_index, rec.start_offset, rec.length, rec.flags)
}

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

pub(crate) fn missed_cleavages_between(cleavage_positions: &[u32], start: u32, end: u32) -> u32 {
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
        assert_eq!(missed_cleavages_between(&cp, 0, 20), 3);
        assert_eq!(missed_cleavages_between(&cp, 5, 10), 0);
        assert_eq!(missed_cleavages_between(&cp, 0, 10), 1);
        assert_eq!(missed_cleavages_between(&cp, 3, 17), 3);
        assert_eq!(missed_cleavages_between(&cp, 10, 11), 0);
    }

    #[test]
    fn matches_strict_offset_semantics() {
        let cp = [2u32, 9, 14, 30];
        assert_eq!(missed_cleavages_between(&cp, cp[0], cp[1]), 0);
        assert_eq!(missed_cleavages_between(&cp, cp[0], cp[2]), 1);
        assert_eq!(missed_cleavages_between(&cp, cp[0], cp[3]), 2);
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn decoy_membership_requires_underscore_delimiter() {
        use crate::decoy::is_decoy_accession;
        assert!(is_decoy_accession("XXX_protein1", "XXX"));
        assert!(
            !is_decoy_accession("XXXprotein1", "XXX"),
            "no '_' delimiter ⇒ a real 'XXX...' target must NOT be read as decoy"
        );
        assert!(!is_decoy_accession("DECOY_protein1", "XXX"));
        assert!(is_decoy_accession("rev_sp|P1", "rev"));
        assert!(!is_decoy_accession("reverse_kinase", "rev"));
    }

    #[test]
    fn multi_enzyme_cleavage_is_union_of_rules() {
        use super::compute_cleavage_positions;
        use model::enzyme::Enzyme;
        let seq = b"MKAESRPE";

        let trypsin = compute_cleavage_positions(seq, Enzyme::Trypsin, &[]);
        assert_eq!(trypsin, vec![0, 2, 6, 8], "Trypsin alone: after K(2), R(6)");

        let gluc = compute_cleavage_positions(seq, Enzyme::GluC, &[]);
        assert_eq!(gluc, vec![0, 4, 8], "GluC alone: after first E(4)");

        let both = compute_cleavage_positions(seq, Enzyme::GluC, &[Enzyme::Trypsin]);
        assert_eq!(
            both,
            vec![0, 2, 4, 6, 8],
            "GluC+Trypsin = union of cut sites"
        );
        let both_rev = compute_cleavage_positions(seq, Enzyme::Trypsin, &[Enzyme::GluC]);
        assert_eq!(both, both_rev, "union is order-independent");
    }

    #[test]
    fn single_enzyme_path_unchanged_with_empty_extras() {
        use super::compute_cleavage_positions;
        use model::enzyme::Enzyme;
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
            let with_nocleave = compute_cleavage_positions(seq, e, &[Enzyme::NoCleavage]);
            assert_eq!(
                single, with_nocleave,
                "{e:?}: NoCleavage extra adds no sites"
            );
        }
        let ns = compute_cleavage_positions(seq, Enzyme::Trypsin, &[Enzyme::NonSpecific]);
        assert_eq!(ns, (0..=seq.len() as u32).collect::<Vec<_>>());
    }

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
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        }
    }

    fn label_mod(name: &str, delta: f64, residue: u8, location: ModLocation) -> Modification {
        use model::modification::ResidueSpec;
        Modification {
            name: name.into(),
            mass_delta: delta,
            residue: if residue == b'*' {
                ResidueSpec::Wildcard
            } else {
                ResidueSpec::Specific(residue)
            },
            location,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
            parts: Vec::new(),
        }
    }

    fn label_names(c: &Candidate) -> Vec<String> {
        c.peptide
            .residues
            .iter()
            .filter_map(|aa| aa.mod_.map(|m| m.name.clone()))
            .collect()
    }

    #[test]
    fn non_specific_digest_labels_every_site() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "silac".into(),
                sequence: b"AKKKKA".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let k8 = label_mod("Label:13C(6)15N(2)", 8.014199, b'K', ModLocation::Anywhere);
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(k8)
            .with_label_channels(vec![("Label:13C(6)15N(2)".into(), 1)])
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 6;
        params.num_tolerable_termini = 0;
        params.max_missed_cleavages = 0;
        params.max_variable_mods_per_peptide = 0;
        params.max_label_sites = 8;
        let full: Vec<_> = enumerate_candidates(&db, &params, "XXX")
            .filter(|c| c.peptide.residues.len() == 6)
            .map(|c| label_names(&c).len())
            .collect();
        assert!(full.contains(&0), "{full:?}");
        assert!(full.contains(&4), "{full:?}");
        assert!(full.iter().all(|&n| n == 0 || n == 4), "{full:?}");
    }

    #[test]
    fn fixed_mod_on_a_labelable_site_is_neutral() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "fixed".into(),
                sequence: b"MAGLKPEPTIDEK".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let mut b = AminoAcidSetBuilder::new_standard();
        let mut channels = Vec::new();
        for (c, (name, d)) in [("Dimethyl", 28.0313), ("Dimethyl:2H(4)", 32.056407)]
            .into_iter()
            .enumerate()
        {
            for (res, loc, suffix) in [
                (b'K', ModLocation::Anywhere, ""),
                (b'*', ModLocation::NTerm, "-N"),
            ] {
                let nm = format!("{name}{suffix}");
                b = b.add_variable_mod(label_mod(&nm, d, res, loc));
                channels.push((nm, c as u8 + 1));
            }
        }
        let mut acetyl = label_mod("Acetyl", 42.010565, b'*', ModLocation::ProtNTerm);
        acetyl.fixed = true;
        let aa_set = b
            .add_fixed_mod(acetyl)
            .with_label_channels(channels)
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params.max_label_sites = 3;
        let kept: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
        let maglk: Vec<Vec<String>> = kept
            .iter()
            .filter(|c| c.peptide.residues.iter().map(|a| a.residue).eq(*b"MAGLK"))
            .map(label_names)
            .collect();
        for ch in ["Dimethyl", "Dimethyl:2H(4)"] {
            assert!(
                maglk.iter().any(|n| n == &["Acetyl", ch]),
                "{ch}: {maglk:?}"
            );
        }
    }

    #[test]
    fn mass_offsets_of_a_large_label_allowance() {
        let mut b = AminoAcidSetBuilder::new_standard();
        let mut channels = Vec::new();
        for (r, d) in [
            (b'K', 8.0),
            (b'R', 9.0),
            (b'L', 10.0),
            (b'I', 11.0),
            (b'V', 12.0),
        ] {
            let name = format!("Label:{}", r as char);
            b = b.add_variable_mod(label_mod(&name, d, r, ModLocation::Anywhere));
            channels.push((name, 1));
        }
        let mut params =
            SearchParams::default_tryptic(b.with_label_channels(channels).build().unwrap());
        params.max_variable_mods_per_peptide = 0;
        params.max_label_sites = 100;
        let mut got: Vec<i64> = crate::candidate_gen::distinct_mod_mass_offsets(&params)
            .iter()
            .map(|&(_, d)| d.round() as i64)
            .collect();
        got.sort_unstable();
        let mut want: Vec<i64> = (0..=100i64).flat_map(|k| 8 * k..=12 * k).collect();
        want.sort_unstable();
        want.dedup();
        assert_eq!(got, want);
        assert!(std::sync::Arc::ptr_eq(
            &crate::candidate_gen::distinct_mod_mass_offsets(&params),
            &crate::candidate_gen::distinct_mod_mass_offsets(&params),
        ));
    }

    #[test]
    fn memoized_offsets_equal_the_plain_walk() {
        fn plain(params: &SearchParams) -> Vec<u64> {
            let deltas = crate::candidate_gen::variable_mod_deltas(params);
            fn rec(
                items: &[crate::candidate_gen::FormDelta],
                start: usize,
                sum: f64,
                budget: (u32, u32),
                out: &mut std::collections::BTreeSet<u64>,
            ) {
                for (i, f) in items.iter().enumerate().skip(start) {
                    if f.mods > budget.0 || f.labels > budget.1 {
                        continue;
                    }
                    let s = sum + f.delta;
                    out.insert(crate::candidate_gen::milli_of(s));
                    rec(items, i, s, (budget.0 - f.mods, budget.1 - f.labels), out);
                }
            }
            let mut out = std::collections::BTreeSet::from([crate::candidate_gen::milli_of(0.0)]);
            let mut channels: Vec<u8> = deltas.iter().map(|f| f.channel).collect();
            channels.push(0);
            channels.sort_unstable();
            channels.dedup();
            for c in channels {
                let items: Vec<_> = deltas
                    .iter()
                    .filter(|f| f.channel == 0 || f.channel == c)
                    .copied()
                    .collect();
                rec(
                    &items,
                    0,
                    0.0,
                    (params.max_variable_mods_per_peptide, params.max_label_sites),
                    &mut out,
                );
            }
            out.into_iter().collect()
        }
        let memo = |params: &SearchParams| -> Vec<u64> {
            let mut v: Vec<u64> = crate::candidate_gen::distinct_mod_mass_offsets(params)
                .iter()
                .map(|&(k, _)| k)
                .collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        let build = |labels: &[(&str, f64, u8, u8)], ox: bool, mods: u32, sites: u32| {
            let mut b = AminoAcidSetBuilder::new_standard();
            for &(name, d, r, _) in labels {
                b = b.add_variable_mod(label_mod(name, d, r, ModLocation::Anywhere));
            }
            if ox {
                b = b.add_variable_mod(oxidation_m());
            }
            let channels = labels
                .iter()
                .map(|&(n, _, _, c)| (n.to_string(), c))
                .collect();
            let mut p =
                SearchParams::default_tryptic(b.with_label_channels(channels).build().unwrap());
            p.max_variable_mods_per_peptide = mods;
            p.max_label_sites = sites;
            p
        };
        let coderabbit = build(
            &[
                ("A", 1.0, b'K', 1),
                ("B", 1.000_000_24, b'L', 1),
                ("C", 10.000_249_9, b'R', 1),
            ],
            false,
            0,
            4,
        );
        let got = memo(&coderabbit);
        assert_eq!(got, plain(&coderabbit));
        assert!(got.contains(&22001), "{got:?}");
        let silac = build(
            &[
                ("Label:2H(4)", 4.025107, b'K', 1),
                ("Label:13C(6)", 6.020129, b'R', 1),
                ("Label:13C(6)15N(2)", 8.014199, b'K', 2),
                ("Label:13C(6)15N(4)", 10.008269, b'R', 2),
            ],
            true,
            1,
            6,
        );
        assert_eq!(memo(&silac), plain(&silac));
    }

    #[test]
    fn mass_offsets_keep_channels_apart() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "offsets".into(),
                sequence: b"PEMKTKIDERAAKGGLK".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let aa = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(label_mod(
                "Label:2H(4)",
                4.025107,
                b'K',
                ModLocation::Anywhere,
            ))
            .add_variable_mod(label_mod(
                "Label:13C(6)15N(2)",
                8.014199,
                b'K',
                ModLocation::Anywhere,
            ))
            .add_variable_mod(oxidation_m())
            .with_label_channels(vec![
                ("Label:2H(4)".into(), 1),
                ("Label:13C(6)15N(2)".into(), 2),
            ])
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 1;
        params.max_label_sites = 3;
        let offsets: Vec<f64> = crate::candidate_gen::distinct_mod_mass_offsets(&params)
            .iter()
            .map(|&(_, d)| d)
            .collect();
        let has = |d: f64| offsets.iter().any(|o| (o - d).abs() < 1e-6);
        assert!(!has(4.025107 + 8.014199), "{offsets:?}");
        assert!(has(3.0 * 8.014199 + 15.994915));
        assert!(!has(2.0 * 15.994915));
        for c in enumerate_candidates(&db, &params, "XXX") {
            let d: f64 = c
                .peptide
                .residues
                .iter()
                .filter_map(|a| a.mod_.map(|m| m.mass_delta))
                .sum();
            assert!(has(d), "{d} missing from {offsets:?}");
        }
    }

    #[test]
    fn silac_labels_are_channel_consistent() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "silac".into(),
                sequence: b"PEPKTIDERAAGGLK".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let k8 = label_mod("Label:13C(6)15N(2)", 8.014199, b'K', ModLocation::Anywhere);
        let r10 = label_mod("Label:13C(6)15N(4)", 10.008269, b'R', ModLocation::Anywhere);
        let build = |labels: bool| {
            let mut b = AminoAcidSetBuilder::new_standard()
                .add_variable_mod(k8.clone())
                .add_variable_mod(r10.clone())
                .add_variable_mod(oxidation_m());
            if labels {
                b = b.with_label_channels(vec![
                    ("Label:13C(6)15N(2)".into(), 1),
                    ("Label:13C(6)15N(4)".into(), 1),
                ]);
            }
            let mut params = SearchParams::default_tryptic(b.build().unwrap());
            params.min_length = 3;
            params.max_variable_mods_per_peptide = 3;
            params.max_label_sites = 3;
            params
        };
        let free = build(false);
        let labelled = build(true);
        let all: Vec<_> = enumerate_candidates(&db, &free, "XXX").collect();
        let kept: Vec<_> = enumerate_candidates(&db, &labelled, "XXX").collect();
        let span = |c: &Candidate| {
            c.peptide
                .residues
                .iter()
                .map(|a| a.residue as char)
                .collect::<String>()
        };
        let n = |v: &[Candidate], s: &str| v.iter().filter(|c| span(c) == s).count();
        assert_eq!(n(&all, "PEPKTIDER"), 4);
        assert_eq!(n(&kept, "PEPKTIDER"), 2);
        for c in kept.iter().filter(|c| span(c) == "PEPKTIDER") {
            let names = label_names(c);
            assert!(names.is_empty() || names.len() == 2, "{names:?}");
        }
        assert_eq!(n(&kept, "TIDER"), 2);
        let mut no_mods = build(true);
        no_mods.max_variable_mods_per_peptide = 0;
        let strict: Vec<_> = enumerate_candidates(&db, &no_mods, "XXX").collect();
        assert_eq!(n(&strict, "PEPKTIDER"), 2);
        assert!(strict.iter().all(|c| c
            .peptide
            .residues
            .iter()
            .all(|a| a.mod_.is_none_or(|m| !m.name.contains("Oxidation")))));
        let tables = crate::candidate_gen::FormWalkTables::new(&labelled);
        let mut walked = 0;
        crate::candidate_gen::walk_forms_bounded(
            b"PEPKTIDER",
            &labelled,
            &tables,
            true,
            false,
            crate::candidate_gen::SpanSites::new(&target.proteins[0], 0),
            0.0,
            f64::MAX,
            |_, _, _| walked += 1,
        );
        assert_eq!(walked, 2);
    }

    #[test]
    fn max_mods_bounds_only_non_label_mods() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "silac".into(),
                sequence: b"PEPMKTIDERAAGGLK".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let k8 = label_mod("Label:13C(6)15N(2)", 8.014199, b'K', ModLocation::Anywhere);
        let r10 = label_mod("Label:13C(6)15N(4)", 10.008269, b'R', ModLocation::Anywhere);
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(k8)
            .add_variable_mod(r10)
            .add_variable_mod(oxidation_m())
            .with_label_channels(vec![
                ("Label:13C(6)15N(2)".into(), 1),
                ("Label:13C(6)15N(4)".into(), 1),
            ])
            .build()
            .unwrap();
        let oxidised = |c: &[&model::amino_acid::AminoAcid]| {
            c.iter()
                .any(|a| a.mod_.is_some_and(|m| m.name.contains("Oxidation")))
        };
        for (max_mods, expect, with_ox) in [(1u32, 4usize, 2usize), (0, 2, 0)] {
            let mut params = SearchParams::default_tryptic(aa_set.clone());
            params.min_length = 3;
            params.max_variable_mods_per_peptide = max_mods;
            params.max_label_sites = 3;
            let forms: Vec<Vec<model::amino_acid::AminoAcid>> =
                enumerate_candidates(&db, &params, "XXX")
                    .filter(|c| {
                        c.peptide
                            .residues
                            .iter()
                            .map(|a| a.residue)
                            .collect::<Vec<_>>()
                            == b"PEPMKTIDER"
                    })
                    .map(|c| c.peptide.residues)
                    .collect();
            assert_eq!(forms.len(), expect, "max_mods {max_mods}");
            let refs: Vec<Vec<&model::amino_acid::AminoAcid>> =
                forms.iter().map(|f| f.iter().collect()).collect();
            assert_eq!(refs.iter().filter(|f| oxidised(f)).count(), with_ox);
            let tables = crate::candidate_gen::FormWalkTables::new(&params);
            let (mut walked, mut walked_ox) = (0, 0);
            crate::candidate_gen::walk_forms_bounded(
                b"PEPMKTIDER",
                &params,
                &tables,
                true,
                false,
                crate::candidate_gen::SpanSites::new(&target.proteins[0], 0),
                0.0,
                f64::MAX,
                |_, chosen, _| {
                    walked += 1;
                    walked_ox += usize::from(oxidised(chosen));
                },
            );
            assert_eq!(
                (walked, walked_ox),
                (expect, with_ox),
                "max_mods {max_mods}"
            );
        }
    }

    #[test]
    fn dimethyl_labels_cover_n_terminus_and_lysine() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "dimethyl".into(),
                sequence: b"MAGLKPEPTIDEK".to_vec(),
                annotation: None,
            }],
        };
        let db = SearchIndex::from_target_db(&target, "XXX");
        let mut b = AminoAcidSetBuilder::new_standard();
        let mut channels = Vec::new();
        for (c, (name, d)) in [
            ("Dimethyl", 28.0313),
            ("Dimethyl:2H(4)", 32.056407),
            ("Dimethyl:2H(6)13C(2)", 36.07567),
        ]
        .into_iter()
        .enumerate()
        {
            for (res, loc, suffix) in [
                (b'K', ModLocation::Anywhere, ""),
                (b'*', ModLocation::NTerm, "-N"),
            ] {
                let nm = format!("{name}{suffix}");
                b = b.add_variable_mod(label_mod(&nm, d, res, loc));
                channels.push((nm, c as u8 + 1));
            }
        }
        let acetyl = label_mod("Acetyl", 42.010565, b'*', ModLocation::ProtNTerm);
        let aa_set = b
            .add_variable_mod(acetyl)
            .with_label_channels(channels)
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 1;
        params.max_label_sites = 3;
        let kept: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
        let of = |s: &str| -> Vec<Vec<String>> {
            kept.iter()
                .filter(|c| {
                    c.peptide
                        .residues
                        .iter()
                        .map(|a| a.residue as char)
                        .collect::<String>()
                        == s
                })
                .map(label_names)
                .collect()
        };
        let pk = of("PEPTIDEK");
        assert_eq!(pk.len(), 4, "{pk:?}");
        for names in &pk {
            assert!(names.is_empty() || names.len() == 2, "{names:?}");
        }
        let mk = of("MAGLK");
        assert_eq!(mk.len(), 8, "{mk:?}");
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
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        }
    }

    fn oxm_fixture() -> (SearchIndex, SearchParams) {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "ox-m fixture".into(),
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
            parts: Vec::new(),
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
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
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
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
            };
            let aa = AminoAcidSetBuilder::new_standard()
                .add_variable_mod(acetyl)
                .build()
                .unwrap();
            let mut params = SearchParams::default_tryptic(aa);
            params.min_length = 3;
            params.max_variable_mods_per_peptide = 2;
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
                parts: Vec::new(),
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

    #[test]
    fn stacked_terminal_forms_enumerate_and_match_lazy() {
        use model::mod_spec::parse_mod_list;
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_mods(
                parse_mod_list(
                    "TMT6plex (K), TMT6plex (N-term), Carbamidomethyl (C)",
                    true,
                )
                .unwrap(),
            )
            .add_mods(parse_mod_list("Oxidation (M), Acetyl (Protein N-term)", false).unwrap())
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
        assert!(has(0, 3, 15.994915 + 42.010565));
        assert!(has(1, 2, 57.021464 + 229.162932));
        assert!(has(1, 2, 57.021464 + 42.010565));
        assert!(!has(1, 2, 42.010565));
        assert!(has(7, 4, 2.0 * 229.162932));
        assert!(has(11, 5, 229.162932 + 15.994915));

        let (_tmp, mi) = build_and_open(&db, &params);
        let mut masses: Vec<f64> = all.iter().map(|c| c.peptide.mass()).collect();
        masses.sort_by(|a, b| a.partial_cmp(b).unwrap());
        masses.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        for &m in &masses {
            assert_lazy_matches_inram(&db, &params, &mi, m, 0.01);
        }
    }

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
            mobility_window: None,
        }
    }

    #[test]
    fn mmap_nominal_window_equals_ram_bucket_window_on_boundary() {
        use model::tolerance::{PrecursorTolerance, Tolerance};
        use std::collections::HashSet;

        let (db, mut params) = oxm_fixture();
        params.isotope_error_range = 0..=0;
        params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Da(0.5));
        let (_tmp, mi) = build_and_open(&db, &params);

        let z = 1u8;

        let masses: Vec<f64> = enumerate_candidates(&db, &params, "XXX")
            .map(|c| c.peptide.mass())
            .collect();
        assert!(!masses.is_empty(), "fixture produced no candidates");

        let mut boundary_seen = false;
        for &m in &masses {
            for step in -8i32..=8 {
                let precursor = m + (step as f64) * 0.125;
                let spec = spectrum_at_mass(precursor, z);

                let ram = ram_window_keys(&db, &params, &spec, z);

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

                let neutral_obs = precursor;
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

        let base = base_mass_of(&db, &params, b"PEPTMK");
        let p_oxm = base + 15.994915;
        assert_lazy_matches_inram(&db, &params, &mi, p_oxm, 0.01);

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

        assert_lazy_matches_inram(&db, &params, &mi, base, 200.0);
    }

    #[test]
    fn lazy_enum_matches_inram_two_variable_mods() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P2".into(),
                description: "two-mod fixture".into(),
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

        let all: Vec<_> = enumerate_candidates(&db, &params, "XXX").collect();
        assert!(!all.is_empty(), "fixture must enumerate candidates");
        let mut masses: Vec<f64> = all.iter().map(|c| c.peptide.mass()).collect();
        masses.sort_by(|a, b| a.partial_cmp(b).unwrap());
        masses.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        for &m in &masses {
            assert_lazy_matches_inram(&db, &params, &mi, m, 0.01);
        }

        let max_mass = masses.last().copied().unwrap();
        let min_mass = masses.first().copied().unwrap();
        let mid = (min_mass + max_mass) / 2.0;
        let half = (max_mass - min_mass) / 2.0 + 1.0;
        assert_lazy_matches_inram(&db, &params, &mi, mid, half);
    }
    #[test]
    fn lazy_enum_matches_inram_at_tolerance_edge() {
        let (db, params) = oxm_fixture();
        let (_tmp, mi) = build_and_open(&db, &params);

        let base = base_mass_of(&db, &params, b"PEPTMK");

        let frac_tol = 0.40_f64;
        let tol = (10.0 + frac_tol) / 1000.0;

        let ox_delta = 15.994915_f64;
        let center = base - tol * 0.999;
        let precursor = center + ox_delta;

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

        let inram_count: usize = enumerate_candidates(&db, &params, "XXX")
            .filter(|c| (c.peptide.mass() - precursor).abs() <= tol)
            .count();
        assert!(
            inram_count > 0,
            "fixture must yield at least one in-RAM candidate at the engineered precursor              (base={base:.6}, tol={tol:.6}, center={center:.6}, precursor={precursor:.6})"
        );

        assert_lazy_matches_inram(&db, &params, &mi, precursor, tol);
    }
}
