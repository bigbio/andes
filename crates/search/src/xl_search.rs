//! Crosslink search driver (`--xl`): doublet-gated (MS-cleavable linkers) or
//! open-mass fragment-index (non-cleavable linkers) pair lookup and alpha/beta
//! co-scoring with the andes core scorer.
//!
//! The design is described in the `andes_xl` crate docs. In short: every
//! peptide is scored alone with `score_psm + psm_edge_score` on the
//! spectrum's [`ScoredSpectrum`], carrying its partner (or the cleaved linker
//! stub) as a mass on the link residue; a CSM is the pair with the best
//! `min(alpha, beta)` among the peptide masses the reporter doublets (or the
//! open-mass retrieval) and the precursor allow.

use std::collections::HashMap;
use std::ops::RangeInclusive;
use std::sync::Mutex;

use andes_xl::crosslinker::{Crosslinker, LinkForm, ISOTOPE_SPACING, N_FORMS, PROTON};
use andes_xl::decoy::PairClass;
use andes_xl::doublet::{
    find_doublets, group_by_mass, stub_ions_present, DoubletMass, DoubletParams,
};
use andes_xl::pair::{pairs_by_sum, MassIndex};
use andes_xl::sites::{link_sites, LinkSite};
use model::amino_acid::AminoAcid;
use model::modification::{leak_mod, ModLocation, Modification, ResidueSpec};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use rayon::prelude::*;
use scoring_crate::scoring::{psm_edge_score, score_psm, RankScorer, ScoredSpectrum};

use crate::candidate_gen::Candidate;
use crate::fragment_index::ChunkFragmentIndex;
use crate::match_engine::PreparedSearch;

/// Settings of the crosslink search.
#[derive(Debug, Clone)]
pub struct XlConfig {
    /// The crosslinking reagent.
    pub linker: Crosslinker,
    /// Reporter-doublet detection.
    pub doublet: DoubletParams,
    /// Tolerance of a doublet-derived peptide mass, in ppm: grouping of
    /// doublets across charges and the peptide lookup.
    pub doublet_mass_tol_ppm: f64,
    /// Precursor isotope offsets tried (0 = monoisotopic; negative = the
    /// precursor picked below the monoisotope).
    pub isotope_offsets: RangeInclusive<i8>,
    /// Lowest precursor charge searched.
    pub min_precursor_charge: u8,
    /// Lightest peptide mass looked up.
    pub min_peptide_mass: f64,
    /// Peptides kept per mass query.
    pub top_per_mass: usize,
    /// Spectra (in input order, among those searched) also run through the
    /// exhaustive pair oracle; 0 disables it.
    pub oracle_spectra: usize,
    /// Open-mass retrieval (non-cleavable linkers): peptides kept per spectrum
    /// from the fragment-index votes.
    pub open_top_k: usize,
    /// Open-mass retrieval: fragment ions a peptide needs to be kept.
    pub open_min_votes: u16,
}

impl XlConfig {
    /// The defaults for `linker`.
    pub fn new(linker: Crosslinker) -> Self {
        XlConfig {
            linker,
            doublet: DoubletParams::default(),
            doublet_mass_tol_ppm: 10.0,
            isotope_offsets: 0..=1,
            min_precursor_charge: 3,
            min_peptide_mass: 400.0,
            top_per_mass: 4,
            oracle_spectra: 0,
            open_top_k: 50,
            open_min_votes: crate::search_params::FRAGMENT_INDEX_MIN_MATCHED,
        }
    }
}

/// A linkable peptidoform, merged over the proteins that contain it.
#[derive(Debug, Clone)]
pub struct XlPeptide {
    /// The peptidoform (flanks of its first occurrence).
    pub peptide: Peptide,
    /// Protein indices of every occurrence, ascending.
    pub proteins: Vec<u32>,
    /// Every occurrence as (protein index, 0-based start in the protein),
    /// ascending.
    pub occurrences: Vec<(u32, u32)>,
    /// The occurrences where the peptide's N-terminus is a link site: protein
    /// N-terminal ones, or all when the linker reacts with every peptide
    /// N-terminus. Ascending.
    pub n_term_occurrences: Vec<(u32, u32)>,
    /// The protein C-terminal occurrences, ascending.
    pub c_term_occurrences: Vec<(u32, u32)>,
    /// True when every occurrence is a decoy.
    pub is_decoy: bool,
    /// Link sites of the linker's first end (of both ends when the linker is
    /// homobifunctional).
    pub sites: Vec<LinkSite>,
    /// Link sites of the second end of a heterobifunctional linker; `None`
    /// for a homobifunctional one.
    pub sites_end2: Option<Vec<LinkSite>>,
}

impl XlPeptide {
    /// Residues as one-letter codes.
    pub fn bare(&self) -> Vec<u8> {
        self.peptide.residues.iter().map(|aa| aa.residue).collect()
    }

    /// The occurrences where `site` is a link site: a terminal site holds only
    /// where the peptide is at that terminus, a residue site everywhere.
    pub fn occurrences_at(&self, site: &LinkSite) -> &[(u32, u32)] {
        if site.n_term {
            &self.n_term_occurrences
        } else if site.c_term {
            &self.c_term_occurrences
        } else {
            &self.occurrences
        }
    }

    /// The sites of linker end `end` (0 or 1).
    pub fn sites_of_end(&self, end: usize) -> &[LinkSite] {
        match (end, &self.sites_end2) {
            (1, Some(v)) => v,
            _ => &self.sites,
        }
    }
}

/// The ends a peptide pair is scored at: `(alpha end, beta end)` for every
/// assignment of the linker's ends to the two peptides.
fn end_assignments(linker: &Crosslinker) -> &'static [(usize, usize)] {
    if linker.is_homobifunctional() {
        &[(0, 0)]
    } else {
        &[(0, 1), (1, 0)]
    }
}

/// The open-mass retrieval index of a non-cleavable search: the fragment-ion
/// index over the linkable base records, and the table row of each indexed
/// peptidoform.
pub struct OpenIndex {
    index: ChunkFragmentIndex,
    row_of: HashMap<String, u32>,
}

/// The linkable peptide table of a crosslink search.
pub struct XlCandidates {
    peptides: Vec<XlPeptide>,
    index: MassIndex,
    open: Option<OpenIndex>,
}

/// Peptidoform identity without flanks: residues and modification masses.
fn peptidoform_key(p: &Peptide) -> String {
    use std::fmt::Write;
    let mut s = String::with_capacity(p.residues.len() * 2);
    for aa in &p.residues {
        s.push(aa.residue as char);
        if let Some(m) = aa.mod_ {
            let _ = write!(s, "[{:+.4}]", m.mass_delta);
        }
    }
    s
}

/// Whether residue 0 carries an N-terminal modification (the amine is taken).
fn n_term_blocked(p: &Peptide) -> bool {
    p.residues.first().is_some_and(|aa| {
        aa.mod_
            .is_some_and(|m| matches!(m.location, ModLocation::NTerm | ModLocation::ProtNTerm))
    })
}

fn sites_of(c: &Candidate, linker: &Crosslinker, end: usize) -> Vec<LinkSite> {
    let end = &linker.ends[end];
    let residues: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
    let mut modified: Vec<bool> = c
        .peptide
        .residues
        .iter()
        .map(|aa| aa.mod_.is_some())
        .collect();
    // Residue 0 is open to an N-terminal link unless an N-terminal mod sits on
    // it; a side-chain mod on residue 0 (e.g. carbamidomethyl Cys) leaves the
    // amine free. Lys side chains with any mod are taken.
    if let Some(first) = modified.first_mut() {
        let side_chain_reactive = end.reacts_with(residues[0]);
        if !side_chain_reactive {
            *first = n_term_blocked(&c.peptide);
        }
    }
    link_sites(
        &residues,
        &modified,
        c.is_protein_n_term,
        c.is_protein_c_term,
        end,
    )
}

impl XlCandidates {
    /// The linkable candidates of the out-of-core index of `prepared`; for a
    /// non-cleavable linker also the open-mass retrieval index over them.
    pub fn from_out_of_core(prepared: &PreparedSearch<'_>, linker: &Crosslinker) -> Self {
        // A span is kept when it has a reactive residue or could be
        // N-terminal; the exact site rule needs the candidate (protein
        // termini, mods) and runs in `from_candidates`.
        let keep = |lk: Crosslinker| {
            move |residues: &[u8], _post: u8| {
                lk.n_term()
                    || lk.ends.iter().any(|e| e.protein_c_term)
                    || residues.iter().any(|&r| lk.reacts_with(r))
            }
        };
        let cands = prepared.out_of_core_candidates_where(keep(linker.clone()));
        let mut table = Self::from_candidates(cands, linker);
        if linker.cleavage.is_none() && !table.peptides.is_empty() {
            let records: Vec<_> = prepared
                .out_of_core_records_where(keep(linker.clone()))
                .into_iter()
                .map(|(rec, _)| rec)
                .collect();
            let masses = &table.index;
            let lo = (0..masses.len() as u32)
                .map(|i| masses.mass(i))
                .fold(f64::MAX, f64::min);
            let hi = (0..masses.len() as u32)
                .map(|i| masses.mass(i))
                .fold(f64::MIN, f64::max);
            let index = ChunkFragmentIndex::build(
                records,
                prepared.idx,
                prepared.params,
                prepared.scorer.feature_match_tolerance(),
                lo - 1.0,
                hi + 1.0,
            );
            let row_of = table
                .peptides
                .iter()
                .enumerate()
                .map(|(i, p)| (peptidoform_key(&p.peptide), i as u32))
                .collect();
            table.open = Some(OpenIndex { index, row_of });
        }
        table
    }

    /// The table over `candidates`, keeping those with at least one link site
    /// on either end of the linker.
    pub fn from_candidates(candidates: Vec<Candidate>, linker: &Crosslinker) -> Self {
        let hetero = !linker.is_homobifunctional();
        let mut by_key: HashMap<String, usize> = HashMap::new();
        let mut peptides: Vec<XlPeptide> = Vec::new();
        let merge = |into: &mut Vec<LinkSite>, from: Vec<LinkSite>| {
            for s in from {
                if !into.contains(&s) {
                    into.push(s);
                }
            }
        };
        let any_n_term = linker.ends.iter().any(|e| e.peptide_n_term);
        for c in candidates {
            let occ = (c.protein_index as u32, c.start_offset_in_protein as u32);
            let n_occ = (any_n_term || c.is_protein_n_term).then_some(occ);
            let c_occ = c.is_protein_c_term.then_some(occ);
            let sites = sites_of(&c, linker, 0);
            let sites_end2 = hetero.then(|| sites_of(&c, linker, 1));
            if sites.is_empty() && sites_end2.as_ref().is_none_or(|v| v.is_empty()) {
                continue;
            }
            let key = peptidoform_key(&c.peptide);
            match by_key.get(&key) {
                Some(&i) => {
                    let p = &mut peptides[i];
                    p.proteins.push(c.protein_index as u32);
                    p.occurrences.push(occ);
                    p.n_term_occurrences.extend(n_occ);
                    p.c_term_occurrences.extend(c_occ);
                    p.is_decoy &= c.is_decoy;
                    merge(&mut p.sites, sites);
                    if let (Some(into), Some(from)) = (p.sites_end2.as_mut(), sites_end2) {
                        merge(into, from);
                    }
                }
                None => {
                    by_key.insert(key, peptides.len());
                    peptides.push(XlPeptide {
                        peptide: c.peptide,
                        proteins: vec![c.protein_index as u32],
                        occurrences: vec![occ],
                        n_term_occurrences: n_occ.into_iter().collect(),
                        c_term_occurrences: c_occ.into_iter().collect(),
                        is_decoy: c.is_decoy,
                        sites,
                        sites_end2,
                    });
                }
            }
        }
        for p in &mut peptides {
            p.proteins.sort_unstable();
            p.proteins.dedup();
            for v in [
                &mut p.occurrences,
                &mut p.n_term_occurrences,
                &mut p.c_term_occurrences,
            ] {
                v.sort_unstable();
                v.dedup();
            }
            p.sites.sort_by_key(|s| (s.pos, s.n_term, s.c_term));
            if let Some(v) = p.sites_end2.as_mut() {
                v.sort_by_key(|s| (s.pos, s.n_term, s.c_term));
            }
        }
        let index = MassIndex::new(peptides.iter().map(|p| p.peptide.mass()).collect());
        XlCandidates {
            peptides,
            index,
            open: None,
        }
    }

    /// The linkable peptides.
    pub fn peptides(&self) -> &[XlPeptide] {
        &self.peptides
    }

    /// The mass index over [`Self::peptides`].
    pub fn index(&self) -> &MassIndex {
        &self.index
    }

    /// Whether the open-mass retrieval index was built.
    pub fn has_open_index(&self) -> bool {
        self.open.is_some()
    }

    /// Table rows of the peptides the fragment-index votes retrieve for
    /// `spec` with neutral mass in `[lo, hi]`, best votes first, at most
    /// `top_k` forms; empty without an open index.
    pub fn open_retrieve(
        &self,
        spec: &Spectrum,
        prepared: &PreparedSearch<'_>,
        lo: f64,
        hi: f64,
        top_k: usize,
        min_votes: u16,
    ) -> Vec<u32> {
        let Some(open) = &self.open else {
            return Vec::new();
        };
        let tol = prepared.scorer.feature_match_tolerance();
        let sel = open.index.query_mass_range(
            spec,
            lo,
            hi,
            tol,
            top_k,
            min_votes,
            crate::fragment_index::vote_peaks_for(tol),
        );
        let cands = open.index.materialise(
            &sel,
            prepared.idx,
            prepared.params,
            None,
            &|_: &mut [Candidate]| {},
        );
        let mut rows: Vec<u32> = Vec::new();
        for c in &cands {
            if let Some(&r) = open.row_of.get(&peptidoform_key(&c.peptide)) {
                if !rows.contains(&r) {
                    rows.push(r);
                }
            }
        }
        rows
    }
}

/// How a candidate mass pair was proposed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairKind {
    /// Both peptide masses carry a reporter doublet.
    Full,
    /// One mass carries a doublet; the other is the precursor complement with
    /// at least one stub ion observed.
    Half,
    /// Open-mass retrieval: one peptide from the fragment-index votes, the
    /// other the precursor complement (non-cleavable linkers).
    Open,
}

impl PairKind {
    pub fn name(self) -> &'static str {
        match self {
            PairKind::Full => "full",
            PairKind::Half => "half",
            PairKind::Open => "open",
        }
    }
}

/// A pair of peptide masses to look up.
#[derive(Debug, Clone, Copy)]
struct MassPair {
    a_mass: f64,
    a_tol: f64,
    b_mass: f64,
    b_tol: f64,
    kind: PairKind,
    iso: i8,
    /// Stub ions of the complement seen (half pairs; the linker's number of
    /// stubs for full pairs).
    b_stub_ions: u8,
}

/// The best placement of the linker on one peptide against one spectrum.
#[derive(Debug, Clone, Copy)]
pub struct PeptideHit {
    /// Index into [`XlCandidates::peptides`].
    pub pep: u32,
    /// The best site.
    pub site: LinkSite,
    /// The best form at that site.
    pub form: LinkForm,
    /// `score_psm + psm_edge_score` of the best (site, form).
    pub score: f32,
    /// Best score of each [`LinkForm`] over the sites (indexed by
    /// [`LinkForm::index`]; `-inf` for a form the linker does not have).
    pub form_scores: [f32; N_FORMS],
    /// Number of sites scored.
    pub n_sites: u8,
    /// Score of the best peptide of another sequence for the same mass query
    /// minus this score (`None` when there is none).
    pub gap_next: Option<f32>,
    /// Peptide mass error against the queried mass, in ppm.
    pub mass_err_ppm: f64,
}

/// Reporter-doublet summary of one spectrum.
#[derive(Debug, Clone, Default)]
pub struct DoubletSummary {
    pub n_peaks: usize,
    pub base_peak: f64,
    pub tic: f64,
    pub n_doublets: usize,
    pub n_doublets_iso: usize,
    /// Distinct doublet peptide masses, most intense first.
    pub masses: Vec<DoubletMass>,
    /// Highest `min(light, heavy)` member intensity of any doublet, relative
    /// to the base peak.
    pub top_rel_intensity: f64,
    pub n_full_pairs: usize,
    pub n_half_pairs: usize,
    /// Doublet-mass counts for the doublet probe, `[floor][iso][spacing]`: floors
    /// [`PROBE_FLOORS`] of the base peak; any doublet, or isotope-confirmed
    /// only; the linker's spacing, or the null spacings [`PROBE_NULL_SHIFTS`]
    /// (mean count over the two).
    pub probe: [[[f32; 2]; 2]; 3],
}

/// Intensity floors (fraction of the base peak) of the doublet probe.
pub const PROBE_FLOORS: [f64; 3] = [0.0, 0.012, 0.05];
/// Shifts of the doublet spacing (Da) that give the null doublet counts: whole
/// Da keep the mass defect of real fragment pairs, and five Da avoid isotope
/// copies of real doublets.
pub const PROBE_NULL_SHIFTS: [f64; 2] = [-5.0, 5.0];

/// The crosslink-spectrum match of one spectrum.
#[derive(Debug, Clone)]
pub struct Csm {
    /// The higher-scoring peptide.
    pub alpha: PeptideHit,
    /// The lower-scoring peptide; its score is the pair score.
    pub beta: PeptideHit,
    pub pair_score: f32,
    /// Best pair score among pairs that differ from this one in at least one
    /// peptide sequence.
    pub second_pair_score: Option<f32>,
    pub kind: PairKind,
    pub iso: i8,
    pub class: PairClass,
    pub intra: bool,
    pub precursor_err_ppm: f64,
    pub b_stub_ions: u8,
    pub n_mass_pairs: usize,
    /// The reporter doublet at alpha's mass, if any.
    pub alpha_doublet: Option<DoubletMass>,
    /// The reporter doublet at beta's mass, if any.
    pub beta_doublet: Option<DoubletMass>,
    /// Re-scores of the chosen (site, form) of each peptide on reduced
    /// spectra, see [`Rescores`].
    pub rescores: Rescores,
}

/// Alpha and beta re-scored at their chosen site and form on spectra with
/// peaks removed, as rescoring features.
#[derive(Debug, Clone, Copy, Default)]
pub struct Rescores {
    /// Without the stub ions of both peptides (the doublet members and the
    /// sulfenic form, at every charge below the precursor's).
    pub alpha_no_stubs: f32,
    pub beta_no_stubs: f32,
    /// Decoupled: alpha without the peaks that beta's b/y ions explain, and
    /// beta without those of alpha.
    pub alpha_decoupled: f32,
    pub beta_decoupled: f32,
}

/// Exhaustive-oracle result of one spectrum.
#[derive(Debug, Clone)]
pub struct OracleResult {
    pub n_pairs: usize,
    /// For a spectrum with a doublet-path CSM: whether the oracle enumerated
    /// its peptide pair, and the oracle's score of that pair.
    pub csm_pair: Option<(bool, f32)>,
    pub best: Option<(u32, u32, f32)>,
    pub class: Option<PairClass>,
    /// Whether each peptide of the best pair has a doublet mass within the
    /// doublet tolerance.
    pub doublet_support: (bool, bool),
}

/// The crosslink-search result of one spectrum.
#[derive(Debug, Clone)]
pub struct XlSpectrumResult {
    pub spec_idx: usize,
    pub charge: u8,
    pub precursor_mass: f64,
    pub doublets: DoubletSummary,
    pub csm: Option<Csm>,
    pub oracle: Option<OracleResult>,
}

/// Leaked stub modifications, one per (form, existing residue mod delta).
/// Bounded: the existing delta comes from the search's modification set.
struct StubMods {
    /// Modification name of each form, by [`LinkForm::index`].
    names: Vec<String>,
    cache: Mutex<HashMap<(LinkForm, u64), &'static Modification>>,
}

impl StubMods {
    fn new(linker: &Crosslinker) -> Self {
        let mut names = vec![String::new(); N_FORMS];
        for form in linker.forms() {
            names[form.index()] = format!("XL:{}", linker.form_name(form));
        }
        StubMods {
            names,
            cache: Mutex::new(HashMap::new()),
        }
    }

    fn get(&self, form: LinkForm, stub: f64, existing: f64) -> &'static Modification {
        let key = (form, existing.to_bits());
        let mut c = self.cache.lock().expect("stub mod cache");
        c.entry(key).or_insert_with(|| {
            let mut m = Modification::new(
                self.names[form.index()].clone(),
                existing + stub,
                ResidueSpec::Wildcard,
                ModLocation::Anywhere,
                false,
            );
            m.blocks_cleavage = true;
            leak_mod(m)
        })
    }
}

/// `peptide` with `delta` on residue `pos`. A stub form is a modification on
/// the residue; the intact form (a per-spectrum partner mass) is folded into
/// the residue mass, which the core scorer reads the same way.
fn with_link_mass(
    peptide: &Peptide,
    pos: usize,
    form: LinkForm,
    delta: f64,
    mods: &StubMods,
) -> Peptide {
    let mut residues: Vec<AminoAcid> = peptide.residues.clone();
    let aa = &mut residues[pos];
    let existing = aa.mod_.map_or(0.0, |m| m.mass_delta);
    match form {
        LinkForm::Intact => aa.mass += delta,
        _ => aa.mod_ = Some(mods.get(form, delta, existing)),
    }
    Peptide::new(residues, peptide.pre, peptide.post)
}

struct SpecCtx<'a> {
    ss: ScoredSpectrum<'a>,
    scorer: &'a RankScorer,
    z: u8,
    frag_tol: f64,
}

/// The best placement of the linker on peptide `pep_idx` at the sites of
/// linker end `end`; `None` when the peptide has no site there.
#[allow(clippy::too_many_arguments)]
fn score_peptide(
    ctx: &SpecCtx<'_>,
    table: &XlCandidates,
    pep_idx: u32,
    partner_mass: f64,
    linker: &Crosslinker,
    mods: &StubMods,
    end: usize,
) -> Option<PeptideHit> {
    let p = &table.peptides[pep_idx as usize];
    let sites = p.sites_of_end(end);
    let mut best: Option<(f32, LinkSite, LinkForm)> = None;
    let mut form_scores = [f32::NEG_INFINITY; N_FORMS];
    for &site in sites {
        for form in linker.forms() {
            let Some(delta) = form.delta(linker, partner_mass) else {
                continue;
            };
            let scored = with_link_mass(&p.peptide, site.pos, form, delta, mods);
            let s = score_psm(&ctx.ss, &scored, ctx.scorer, ctx.z, ctx.frag_tol)
                + psm_edge_score(&ctx.ss, &scored, ctx.scorer, ctx.z) as f32;
            let fi = form.index();
            if s > form_scores[fi] {
                form_scores[fi] = s;
            }
            if best.is_none_or(|(b, _, _)| s > b) {
                best = Some((s, site, form));
            }
        }
    }
    let (score, site, form) = best?;
    Some(PeptideHit {
        pep: pep_idx,
        site,
        form,
        score,
        form_scores,
        n_sites: sites.len().min(255) as u8,
        gap_next: None,
        mass_err_ppm: 0.0,
    })
}

/// Score every peptide within `tol` of `mass`; best first, ties by index.
/// Each hit's `gap_next` is filled against the best hit of another sequence.
#[allow(clippy::too_many_arguments)]
fn query(
    ctx: &SpecCtx<'_>,
    table: &XlCandidates,
    mass: f64,
    tol: f64,
    partner_mass: f64,
    linker: &Crosslinker,
    mods: &StubMods,
    keep: usize,
    end: usize,
) -> Vec<PeptideHit> {
    let hits: Vec<PeptideHit> = table
        .index
        .within(mass, tol)
        .iter()
        .filter_map(|&i| {
            let mut h = score_peptide(ctx, table, i, partner_mass, linker, mods, end)?;
            h.mass_err_ppm = (table.index.mass(i) - mass) / mass * 1e6;
            Some(h)
        })
        .collect();
    rank_hits(hits, table, keep)
}

/// `hits` best first (ties by index), each with its `gap_next` against the
/// best hit of another sequence, cut to `keep`.
fn rank_hits(mut hits: Vec<PeptideHit>, table: &XlCandidates, keep: usize) -> Vec<PeptideHit> {
    hits.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.pep.cmp(&b.pep)));
    let bare: Vec<Vec<u8>> = hits
        .iter()
        .map(|h| table.peptides[h.pep as usize].bare())
        .collect();
    for k in 0..hits.len() {
        hits[k].gap_next = (0..hits.len())
            .find(|&j| j != k && bare[j] != bare[k])
            .map(|j| hits[k].score - hits[j].score);
    }
    hits.truncate(keep);
    hits
}

/// Accession of protein `pi` with the decoy affix removed: a target protein
/// and its decoy share it.
fn protein_family<'p>(prepared: &'p PreparedSearch<'_>, pi: u32) -> &'p str {
    let idx = prepared.idx;
    let Some(prot) = idx.db.proteins.get(pi as usize) else {
        return "";
    };
    let acc = prot.accession.as_str();
    let needle = crate::decoy::decoy_accession_needle(&idx.decoy_prefix);
    if let Some(rest) = acc.strip_prefix(needle.as_str()) {
        return rest;
    }
    match idx.decoy_suffix.as_deref() {
        Some(sfx) if !sfx.is_empty() => acc.strip_suffix(sfx).unwrap_or(acc),
        _ => acc,
    }
}

/// The intra/inter class of a pair, by one rule for targets and decoys:
/// intra when the peptides' proteins share a family (a protein and its decoy
/// are one family, so a target peptide linked to a decoy peptide of its own
/// protein is an intra TD).
fn is_intra(prepared: &PreparedSearch<'_>, a: &XlPeptide, b: &XlPeptide) -> bool {
    let fam_b: Vec<&str> = b
        .proteins
        .iter()
        .map(|&p| protein_family(prepared, p))
        .collect();
    a.proteins
        .iter()
        .any(|&p| fam_b.contains(&protein_family(prepared, p)))
}

fn precursor_window_da(params: &crate::search_params::SearchParams, mass: f64) -> f64 {
    let t = &params.precursor_tolerance;
    t.left.as_da(mass).max(t.right.as_da(mass))
}

/// m/z-sorted peaks of `spec` (a copy only when the input is unsorted).
fn sorted_peaks(spec: &Spectrum) -> std::borrow::Cow<'_, [(f64, f32)]> {
    if spec.peaks.windows(2).all(|w| w[0].0 <= w[1].0) {
        std::borrow::Cow::Borrowed(&spec.peaks)
    } else {
        let mut p = spec.peaks.clone();
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        std::borrow::Cow::Owned(p)
    }
}

fn mass_pairs(
    peaks: &[(f64, f32)],
    masses: &[DoubletMass],
    precursor_mass: f64,
    prec_tol: f64,
    z: u8,
    cfg: &XlConfig,
) -> Vec<MassPair> {
    let Some(cleavage) = cfg.linker.cleavage.as_ref() else {
        return Vec::new();
    };
    let n_stub_forms = cleavage.stubs.len() as u8;
    let dtol = |m: f64| m * cfg.doublet_mass_tol_ppm * 1e-6;
    let mut out: Vec<MassPair> = Vec::new();
    for iso in cfg.isotope_offsets.clone() {
        let target = precursor_mass - iso as f64 * ISOTOPE_SPACING - cfg.linker.intact;
        for (i, a) in masses.iter().enumerate() {
            for b in &masses[i..] {
                let tol = prec_tol + dtol(a.mass) + dtol(b.mass);
                if (a.mass + b.mass - target).abs() <= tol {
                    let (hi, lo) = if a.mass >= b.mass { (a, b) } else { (b, a) };
                    out.push(MassPair {
                        a_mass: hi.mass,
                        a_tol: dtol(hi.mass),
                        b_mass: lo.mass,
                        b_tol: dtol(lo.mass),
                        kind: PairKind::Full,
                        iso,
                        b_stub_ions: n_stub_forms,
                    });
                }
            }
        }
    }
    for iso in cfg.isotope_offsets.clone() {
        let target = precursor_mass - iso as f64 * ISOTOPE_SPACING - cfg.linker.intact;
        for a in masses {
            let b_mass = target - a.mass;
            if b_mass < cfg.min_peptide_mass {
                continue;
            }
            let b_tol = prec_tol + dtol(a.mass);
            // A complement that is itself a doublet mass is a full pair.
            if out.iter().any(|p| {
                p.kind == PairKind::Full && p.iso == iso && {
                    let close = |x: f64, y: f64| (x - y).abs() <= b_tol;
                    (close(p.a_mass, a.mass) && close(p.b_mass, b_mass))
                        || (close(p.b_mass, a.mass) && close(p.a_mass, b_mass))
                }
            }) {
                continue;
            }
            let stubs = stub_ions_present(
                peaks,
                b_mass,
                cleavage,
                z.saturating_sub(1).max(1),
                cfg.doublet.tol_ppm,
            );
            let n_stubs = stubs.iter().filter(|&&s| s).count() as u8;
            if n_stubs == 0 {
                continue;
            }
            out.push(MassPair {
                a_mass: a.mass,
                a_tol: dtol(a.mass),
                b_mass,
                b_tol,
                kind: PairKind::Half,
                iso,
                b_stub_ions: n_stubs,
            });
        }
    }
    out.retain(|p| p.a_mass >= cfg.min_peptide_mass && p.b_mass >= cfg.min_peptide_mass);
    out
}

/// The scoring peptide of a hit: its peptide with the hit's form at its site.
fn hit_peptide(
    table: &XlCandidates,
    h: &PeptideHit,
    partner_mass: f64,
    linker: &Crosslinker,
    mods: &StubMods,
) -> Option<Peptide> {
    let p = &table.peptides[h.pep as usize];
    let delta = h.form.delta(linker, partner_mass)?;
    Some(with_link_mass(&p.peptide, h.site.pos, h.form, delta, mods))
}

/// m/z of the b and y ions of `p` at charges `1..=max_z`.
fn by_ions(p: &Peptide, max_z: u8) -> Vec<f64> {
    let total: f64 = p
        .residues
        .iter()
        .map(|aa| aa.mass + aa.mod_.map_or(0.0, |m| m.mass_delta))
        .sum();
    let mut out = Vec::new();
    let mut prefix = 0.0;
    for aa in &p.residues[..p.residues.len().saturating_sub(1)] {
        prefix += aa.mass + aa.mod_.map_or(0.0, |m| m.mass_delta);
        let b = prefix;
        let y = total - prefix + model::mass::H2O;
        for z in 1..=max_z.max(1) {
            let zf = z as f64;
            out.push(b / zf + PROTON);
            out.push(y / zf + PROTON);
        }
    }
    out
}

/// `spec` without the peaks within `tol` Da of any of `mzs`.
fn without_peaks(spec: &Spectrum, mzs: &[f64], tol: f64) -> Spectrum {
    let mut sorted = mzs.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mut out = spec.clone();
    out.peaks.retain(|&(mz, _)| {
        let i = sorted.partition_point(|&m| m < mz - tol);
        !sorted.get(i).is_some_and(|&m| m <= mz + tol)
    });
    out
}

fn rescore(
    spec: &Spectrum,
    prepared: &PreparedSearch<'_>,
    z: u8,
    peptides: &[&Peptide],
) -> Vec<f32> {
    let ss = ScoredSpectrum::new(spec, prepared.scorer, z);
    peptides
        .iter()
        .map(|p| {
            score_psm(&ss, p, prepared.scorer, z, prepared.fragment_tolerance_da)
                + psm_edge_score(&ss, p, prepared.scorer, z) as f32
        })
        .collect()
}

/// [`Rescores`] of the chosen pair.
#[allow(clippy::too_many_arguments)]
fn rescores(
    spec: &Spectrum,
    prepared: &PreparedSearch<'_>,
    z: u8,
    table: &XlCandidates,
    alpha: &PeptideHit,
    beta: &PeptideHit,
    linker: &Crosslinker,
    mods: &StubMods,
) -> Rescores {
    let ma = table.index.mass(alpha.pep);
    let mb = table.index.mass(beta.pep);
    let (Some(pa), Some(pb)) = (
        hit_peptide(table, alpha, mb, linker, mods),
        hit_peptide(table, beta, ma, linker, mods),
    ) else {
        return Rescores::default();
    };
    let max_z = z.saturating_sub(1).max(1);
    let tol = prepared.fragment_tolerance_da;
    let mut stubs = Vec::new();
    if let Some(c) = linker.cleavage.as_ref() {
        for m in [ma, mb] {
            for stub in c.stub_masses() {
                for q in 1..=max_z {
                    stubs.push((m + stub) / q as f64 + PROTON);
                }
            }
        }
    }
    let ns = rescore(&without_peaks(spec, &stubs, tol), prepared, z, &[&pa, &pb]);
    let frag_z = max_z.min(2);
    let da = rescore(
        &without_peaks(spec, &by_ions(&pb, frag_z), tol),
        prepared,
        z,
        &[&pa],
    );
    let db = rescore(
        &without_peaks(spec, &by_ions(&pa, frag_z), tol),
        prepared,
        z,
        &[&pb],
    );
    Rescores {
        alpha_no_stubs: ns[0],
        beta_no_stubs: ns[1],
        alpha_decoupled: da[0],
        beta_decoupled: db[0],
    }
}

/// The doublet mass within the doublet tolerance of `mass`, most intense
/// first.
fn doublet_at(masses: &[DoubletMass], mass: f64, cfg: &XlConfig) -> Option<DoubletMass> {
    let tol = mass * cfg.doublet_mass_tol_ppm * 1e-6;
    masses
        .iter()
        .find(|d| (d.mass - mass).abs() <= tol)
        .cloned()
}

fn summarize_doublets(peaks: &[(f64, f32)], z: u8, cfg: &XlConfig) -> DoubletSummary {
    let mut s = DoubletSummary {
        n_peaks: peaks.len(),
        base_peak: peaks.iter().map(|p| p.1 as f64).fold(0.0, f64::max),
        tic: peaks.iter().map(|p| p.1 as f64).sum(),
        ..Default::default()
    };
    let Some(cleavage) = cfg.linker.cleavage.as_ref() else {
        return s;
    };
    let ds = find_doublets(peaks, z, cleavage, &cfg.doublet);
    s.n_doublets = ds.len();
    s.n_doublets_iso = ds.iter().filter(|d| d.isotope_confirmed).count();
    if s.base_peak > 0.0 {
        s.top_rel_intensity = ds
            .iter()
            .map(|d| peaks[d.light_idx].1.min(peaks[d.heavy_idx].1) as f64 / s.base_peak)
            .fold(0.0, f64::max);
    }
    s.masses = group_by_mass(&ds, cfg.doublet_mass_tol_ppm);
    for (fi, &floor) in PROBE_FLOORS.iter().enumerate() {
        let p = DoubletParams {
            min_rel_intensity: floor,
            ..cfg.doublet
        };
        let count = |c: &andes_xl::crosslinker::Cleavage| {
            let ds = find_doublets(peaks, z, c, &p);
            let iso: Vec<_> = ds.iter().copied().filter(|d| d.isotope_confirmed).collect();
            [
                group_by_mass(&ds, cfg.doublet_mass_tol_ppm).len() as f32,
                group_by_mass(&iso, cfg.doublet_mass_tol_ppm).len() as f32,
            ]
        };
        let real = count(cleavage);
        let mut null = [0.0_f32; 2];
        for shift in PROBE_NULL_SHIFTS {
            let n = count(&cleavage.with_heavy_shift(shift));
            null[0] += n[0] / PROBE_NULL_SHIFTS.len() as f32;
            null[1] += n[1] / PROBE_NULL_SHIFTS.len() as f32;
        }
        s.probe[fi] = [[real[0], null[0]], [real[1], null[1]]];
    }
    s
}

/// The crosslink search of one spectrum at its precursor charge. A spectrum
/// without one is searched at every charge of the search's charge range from
/// [`XlConfig::min_precursor_charge`] up, and the charge whose CSM has the best
/// pair score is kept (the lowest charge when none has a CSM).
#[allow(clippy::too_many_arguments)]
fn search_spectrum(
    spec_idx: usize,
    spec: &Spectrum,
    prepared: &PreparedSearch<'_>,
    table: &XlCandidates,
    cfg: &XlConfig,
    mods: &StubMods,
    run_oracle: bool,
) -> Option<XlSpectrumResult> {
    match spec.precursor_charge {
        Some(z) => {
            let z = u8::try_from(z).ok()?;
            if z < cfg.min_precursor_charge {
                return None;
            }
            search_spectrum_at(spec_idx, spec, z, prepared, table, cfg, mods, run_oracle)
        }
        None => {
            let range = &prepared.params.charge_range;
            let lo = (*range.start()).max(cfg.min_precursor_charge);
            let mut best: Option<XlSpectrumResult> = None;
            for z in lo..=*range.end() {
                let Some(r) =
                    search_spectrum_at(spec_idx, spec, z, prepared, table, cfg, mods, false)
                else {
                    continue;
                };
                let score = |r: &XlSpectrumResult| r.csm.as_ref().map(|c| c.pair_score);
                let better = match &best {
                    None => true,
                    Some(b) => match (score(&r), score(b)) {
                        (Some(x), Some(y)) => x > y,
                        (Some(_), None) => true,
                        _ => false,
                    },
                };
                if better {
                    best = Some(r);
                }
            }
            best
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn search_spectrum_at(
    spec_idx: usize,
    spec: &Spectrum,
    z: u8,
    prepared: &PreparedSearch<'_>,
    table: &XlCandidates,
    cfg: &XlConfig,
    mods: &StubMods,
    run_oracle: bool,
) -> Option<XlSpectrumResult> {
    let precursor_mass = crate::precursor_cal::adjusted_observed_neutral_mass(
        (spec.precursor_mz - PROTON) * z as f64,
        prepared.params.precursor_mass_shift_ppm,
    );
    let peaks = sorted_peaks(spec);
    let mut doublets = summarize_doublets(&peaks, z, cfg);
    let prec_tol = precursor_window_da(prepared.params, precursor_mass);
    let pairs = mass_pairs(&peaks, &doublets.masses, precursor_mass, prec_tol, z, cfg);
    doublets.n_full_pairs = pairs.iter().filter(|p| p.kind == PairKind::Full).count();
    doublets.n_half_pairs = pairs.len() - doublets.n_full_pairs;

    let linker = &cfg.linker;
    // Open-mass retrieval (non-cleavable linkers): the peptides the fragment
    // ions vote for, at any mass that leaves room for a partner.
    let open_rows = if table.has_open_index() {
        table.open_retrieve(
            spec,
            prepared,
            cfg.min_peptide_mass - prec_tol,
            precursor_mass - linker.intact - cfg.min_peptide_mass + prec_tol,
            cfg.open_top_k,
            cfg.open_min_votes,
        )
    } else {
        Vec::new()
    };

    let mut result = XlSpectrumResult {
        spec_idx,
        charge: z,
        precursor_mass,
        doublets,
        csm: None,
        oracle: None,
    };
    if pairs.is_empty() && open_rows.is_empty() && !run_oracle {
        return Some(result);
    }
    let ctx = SpecCtx {
        ss: ScoredSpectrum::new(spec, prepared.scorer, z),
        scorer: prepared.scorer,
        z,
        frag_tol: prepared.fragment_tolerance_da,
    };

    // Every (alpha, beta) combination of the kept hits of every mass pair.
    struct Cand {
        a: PeptideHit,
        b: PeptideHit,
        score: f32,
        pair: MassPair,
    }
    let mut combos: Vec<Cand> = Vec::new();
    for (mp, &(ea, eb)) in pairs
        .iter()
        .flat_map(|mp| end_assignments(linker).iter().map(move |e| (mp, e)))
    {
        let a_hits = query(
            &ctx,
            table,
            mp.a_mass,
            mp.a_tol,
            mp.b_mass,
            linker,
            mods,
            cfg.top_per_mass,
            ea,
        );
        if a_hits.is_empty() {
            continue;
        }
        let b_hits = query(
            &ctx,
            table,
            mp.b_mass,
            mp.b_tol,
            mp.a_mass,
            linker,
            mods,
            cfg.top_per_mass,
            eb,
        );
        let iso_mass = mp.iso as f64 * ISOTOPE_SPACING;
        for a in &a_hits {
            for b in &b_hits {
                // The pair must fit the precursor within the search's precursor
                // tolerance, as a linear peptide does; the per-mass windows
                // above are only the lookup.
                let theo =
                    table.index.mass(a.pep) + table.index.mass(b.pep) + linker.intact + iso_mass;
                if (precursor_mass - theo).abs() > prec_tol {
                    continue;
                }
                combos.push(Cand {
                    a: *a,
                    b: *b,
                    score: a.score.min(b.score),
                    pair: *mp,
                });
            }
        }
    }
    // Open-mass pairs: each retrieved peptide is scored with the precursor
    // complement as its partner, and the partner is looked up by that mass.
    for iso in cfg.isotope_offsets.clone() {
        if open_rows.is_empty() {
            break;
        }
        let target = precursor_mass - iso as f64 * ISOTOPE_SPACING - linker.intact;
        for &(ea, eb) in end_assignments(linker) {
            let a_hits: Vec<PeptideHit> = open_rows
                .iter()
                .filter_map(|&r| {
                    let b_mass = target - table.index.mass(r);
                    if b_mass < cfg.min_peptide_mass - prec_tol {
                        return None;
                    }
                    score_peptide(&ctx, table, r, b_mass, linker, mods, ea)
                })
                .collect();
            for a in rank_hits(a_hits, table, usize::MAX) {
                let a_mass = table.index.mass(a.pep);
                let b_mass = target - a_mass;
                let b_hits = query(
                    &ctx,
                    table,
                    b_mass,
                    prec_tol,
                    a_mass,
                    linker,
                    mods,
                    cfg.top_per_mass,
                    eb,
                );
                for b in b_hits {
                    let theo = a_mass
                        + table.index.mass(b.pep)
                        + linker.intact
                        + iso as f64 * ISOTOPE_SPACING;
                    if (precursor_mass - theo).abs() > prec_tol {
                        continue;
                    }
                    combos.push(Cand {
                        a,
                        b,
                        score: a.score.min(b.score),
                        pair: MassPair {
                            a_mass,
                            a_tol: 0.0,
                            b_mass,
                            b_tol: prec_tol,
                            kind: PairKind::Open,
                            iso,
                            b_stub_ions: 0,
                        },
                    });
                }
            }
        }
    }
    let seq_pair = |c: &Cand| {
        let mut v = [
            table.peptides[c.a.pep as usize].bare(),
            table.peptides[c.b.pep as usize].bare(),
        ];
        v.sort();
        v
    };
    // Best by pair score; ties by full pairs first, then candidate indices.
    let best = combos.iter().enumerate().max_by(|(_, x), (_, y)| {
        x.score
            .total_cmp(&y.score)
            .then((x.pair.kind == PairKind::Full).cmp(&(y.pair.kind == PairKind::Full)))
            .then(y.a.pep.min(y.b.pep).cmp(&x.a.pep.min(x.b.pep)))
            .then(y.a.pep.max(y.b.pep).cmp(&x.a.pep.max(x.b.pep)))
    });
    if let Some((bi, c)) = best {
        let best_seqs = seq_pair(c);
        let second = combos
            .iter()
            .enumerate()
            .filter(|(i, o)| *i != bi && seq_pair(o) != best_seqs)
            .map(|(_, o)| o.score)
            .max_by(f32::total_cmp);
        let (alpha, beta) = if c.a.score >= c.b.score {
            (c.a, c.b)
        } else {
            (c.b, c.a)
        };
        let pa = &table.peptides[alpha.pep as usize];
        let pb = &table.peptides[beta.pep as usize];
        let theo = pa.peptide.mass()
            + pb.peptide.mass()
            + linker.intact
            + c.pair.iso as f64 * ISOTOPE_SPACING;
        result.csm = Some(Csm {
            alpha,
            beta,
            pair_score: c.score,
            second_pair_score: second,
            kind: c.pair.kind,
            iso: c.pair.iso,
            class: PairClass::of(pa.is_decoy, pb.is_decoy),
            intra: is_intra(prepared, pa, pb),
            precursor_err_ppm: (precursor_mass - theo) / theo * 1e6,
            b_stub_ions: c.pair.b_stub_ions,
            n_mass_pairs: pairs.len() + open_rows.len(),
            alpha_doublet: doublet_at(&result.doublets.masses, pa.peptide.mass(), cfg),
            beta_doublet: doublet_at(&result.doublets.masses, pb.peptide.mass(), cfg),
            rescores: rescores(spec, prepared, z, table, &alpha, &beta, linker, mods),
        });
    }

    if run_oracle {
        result.oracle = Some(oracle(
            &ctx,
            table,
            precursor_mass,
            prec_tol,
            &result,
            cfg,
            mods,
        ));
    }
    Some(result)
}

/// Exhaustive search of one spectrum: every pair of linkable peptides whose
/// summed mass fits the precursor (any isotope offset of the config), no
/// doublet gate. Each peptide is scored with its partner taken as the exact
/// precursor complement.
fn oracle(
    ctx: &SpecCtx<'_>,
    table: &XlCandidates,
    precursor_mass: f64,
    prec_tol: f64,
    result: &XlSpectrumResult,
    cfg: &XlConfig,
    mods: &StubMods,
) -> OracleResult {
    let linker = &cfg.linker;
    let mut n_pairs = 0;
    let mut best: Option<(u32, u32, f32, i8)> = None;
    let mut scored: HashMap<(u32, i8, usize), f32> = HashMap::new();
    let csm_key = result.csm.as_ref().map(|c| {
        let (i, j) = (c.alpha.pep, c.beta.pep);
        ((i.min(j), i.max(j)), c.iso)
    });
    let mut csm_pair: Option<(bool, f32)> = csm_key.map(|_| (false, f32::NEG_INFINITY));
    for iso in cfg.isotope_offsets.clone() {
        let target = precursor_mass - iso as f64 * ISOTOPE_SPACING - linker.intact;
        let pairs = pairs_by_sum(&table.index, target, prec_tol);
        n_pairs += pairs.len();
        for (i, j) in pairs {
            let mut s = |k: u32, end: usize| -> f32 {
                *scored.entry((k, iso, end)).or_insert_with(|| {
                    let partner = target - table.index.mass(k);
                    score_peptide(ctx, table, k, partner, linker, mods, end)
                        .map_or(f32::NEG_INFINITY, |h| h.score)
                })
            };
            let pair = end_assignments(linker)
                .iter()
                .map(|&(ea, eb)| s(i, ea).min(s(j, eb)))
                .fold(f32::NEG_INFINITY, f32::max);
            if csm_key == Some(((i, j), iso)) {
                csm_pair = Some((true, pair));
            }
            if best.is_none_or(|(bi, bj, bs, _)| pair > bs || (pair == bs && (i, j) < (bi, bj))) {
                best = Some((i, j, pair, iso));
            }
        }
    }
    let support = |k: u32| {
        let m = table.index.mass(k);
        result
            .doublets
            .masses
            .iter()
            .any(|d| (d.mass - m).abs() <= m * cfg.doublet_mass_tol_ppm * 1e-6)
    };
    OracleResult {
        n_pairs,
        csm_pair,
        best: best.map(|(i, j, s, _)| (i, j, s)),
        class: best.map(|(i, j, _, _)| {
            PairClass::of(
                table.peptides[i as usize].is_decoy,
                table.peptides[j as usize].is_decoy,
            )
        }),
        doublet_support: best.map_or((false, false), |(i, j, _, _)| (support(i), support(j))),
    }
}

/// Search `spectra` for crosslinks of `cfg.linker` against `table`.
pub fn xl_search_run(
    spectra: &[Spectrum],
    prepared: &PreparedSearch<'_>,
    table: &XlCandidates,
    cfg: &XlConfig,
) -> Vec<XlSpectrumResult> {
    let mods = StubMods::new(&cfg.linker);
    // The oracle runs on the first `oracle_spectra` spectra that pass the
    // charge gate, in input order.
    let eligible: Vec<usize> = spectra
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            s.precursor_charge
                .is_some_and(|z| z >= i32::from(cfg.min_precursor_charge))
        })
        .map(|(i, _)| i)
        .take(cfg.oracle_spectra)
        .collect();
    spectra
        .par_iter()
        .enumerate()
        .filter_map(|(i, spec)| {
            let run_oracle = eligible.binary_search(&i).is_ok();
            search_spectrum(i, spec, prepared, table, cfg, &mods, run_oracle)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'K', b'A')
    }

    #[test]
    fn by_ions_and_peak_removal() {
        let p = pep(b"GAK");
        let ions = by_ions(&p, 1);
        // b1 = G + proton, y1 = K + water + proton.
        let g = AminoAcid::standard(b'G').unwrap().mass;
        let k = AminoAcid::standard(b'K').unwrap().mass;
        assert!(ions.iter().any(|&m| (m - (g + PROTON)).abs() < 1e-9));
        assert!(ions
            .iter()
            .any(|&m| (m - (k + model::mass::H2O + PROTON)).abs() < 1e-9));
        assert_eq!(ions.len(), 4);
        let spec = Spectrum {
            title: "t".into(),
            precursor_mz: 0.0,
            precursor_intensity: None,
            precursor_charge: Some(3),
            rt_seconds: None,
            scan: None,
            peaks: vec![(g + PROTON + 0.001, 1.0), (500.0, 2.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        };
        let kept = without_peaks(&spec, &ions, 0.01);
        assert_eq!(kept.peaks, vec![(500.0, 2.0)]);
    }

    #[test]
    fn link_mass_is_on_the_site_and_in_the_peptide_mass() {
        let x = Crosslinker::dsso();
        let mods = StubMods::new(&x);
        let p = pep(b"PEPKTIDE");
        let stub = LinkForm::Stub(0).delta(&x, 0.0).unwrap();
        let q = with_link_mass(&p, 3, LinkForm::Stub(0), stub, &mods);
        assert_eq!(q.residues[3].mod_.unwrap().name, "XL:alkene");
        assert!((q.mass() - p.mass() - stub).abs() < 1e-9);
        assert!((q.residues[3].mod_.unwrap().mass_delta - stub).abs() < 1e-12);
        // The same form reuses one modification record.
        let r = with_link_mass(&p, 3, LinkForm::Stub(0), stub, &mods);
        assert!(std::ptr::eq(
            q.residues[3].mod_.unwrap(),
            r.residues[3].mod_.unwrap()
        ));
        // The intact form folds the partner into the residue mass.
        let intact = LinkForm::Intact.delta(&x, 1000.0).unwrap();
        let s = with_link_mass(&p, 3, LinkForm::Intact, intact, &mods);
        assert!((s.mass() - p.mass() - intact).abs() < 1e-9);
        assert!(s.residues[3].mod_.is_none());
        assert_eq!(
            s.nominal_residue_mass(),
            model::mass::nominal_from(s.mass() - model::mass::H2O)
        );
    }

    #[test]
    fn a_terminal_site_is_reported_only_where_the_peptide_is_terminal() {
        let cand = |protein: usize, start: usize, n_term: bool| Candidate {
            peptide: pep(b"AKR"),
            protein_index: protein,
            start_offset_in_protein: start,
            is_decoy: false,
            is_protein_n_term: n_term,
            is_protein_c_term: false,
        };
        // AKR is protein N-terminal in protein 0 and internal in protein 1.
        let t = XlCandidates::from_candidates(
            vec![cand(0, 0, true), cand(1, 7, false)],
            &Crosslinker::dsso(),
        );
        assert_eq!(t.peptides().len(), 1);
        let p = &t.peptides()[0];
        let n = p.sites.iter().find(|s| s.n_term).unwrap();
        let k = p.sites.iter().find(|s| !s.n_term).unwrap();
        assert_eq!(p.occurrences_at(n), &[(0, 0)]);
        assert_eq!(p.occurrences_at(k), &[(0, 0), (1, 7)]);
    }

    #[test]
    fn mass_pairs_find_full_and_half_pairs() {
        let cfg = XlConfig::new(Crosslinker::dsso());
        let thiol = cfg.linker.cleavage.as_ref().unwrap().stubs[1].mass;
        let (ma, mb) = (1500.0, 900.0);
        let prec = ma + mb + cfg.linker.intact;
        let dm = |m: f64| DoubletMass {
            mass: m,
            intensity: 1.0,
            n_charges: 1,
            isotope_confirmed: false,
            min_charge: 1,
        };
        // Both doublets: one full pair (the half pairs it implies are dropped).
        let full = mass_pairs(&[], &[dm(ma), dm(mb)], prec, 0.01, 4, &cfg);
        assert_eq!(full.len(), 1);
        assert_eq!(full[0].kind, PairKind::Full);
        assert_eq!(full[0].b_stub_ions, 3);
        assert!((full[0].a_mass - ma).abs() < 1e-9 && (full[0].b_mass - mb).abs() < 1e-9);
        // One doublet: a half pair only with a stub ion of the complement.
        assert!(mass_pairs(&[], &[dm(ma)], prec, 0.01, 4, &cfg).is_empty());
        let stub_peak = [((mb + thiol) / 2.0 + PROTON, 10.0_f32)];
        let half = mass_pairs(&stub_peak, &[dm(ma)], prec, 0.01, 4, &cfg);
        assert_eq!(half.len(), 1);
        assert_eq!(half[0].kind, PairKind::Half);
        assert_eq!(half[0].iso, 0);
        assert_eq!(half[0].b_stub_ions, 1);
    }
}
