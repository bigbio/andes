use std::collections::HashMap;
use std::ops::RangeInclusive;
use std::sync::Mutex;

use model::amino_acid::AminoAcid;
use model::modification::{leak_mod, ModLocation, Modification, ResidueSpec};
use model::peptide::Peptide;
use model::scan::ProductScan;
use model::spectrum::Spectrum;
use rayon::prelude::*;
use scoring_crate::scoring::{psm_edge_score, score_psm, RankScorer, ScoredSpectrum};
use xl::crosslinker::{Crosslinker, LinkForm, ISOTOPE_SPACING, N_FORMS, PROTON};
use xl::decoy::PairClass;
use xl::doublet::{find_doublets, group_by_mass, stub_ions_present, DoubletMass, DoubletParams};
use xl::pair::{pairs_by_sum, MassIndex};
use xl::sites::{link_sites, LinkSite};

use crate::candidate_gen::Candidate;
use crate::fragment_index::ChunkFragmentIndex;
use crate::match_engine::PreparedSearch;

#[derive(Debug, Clone)]
pub struct XlConfig {
    pub linker: Crosslinker,
    pub doublet: DoubletParams,
    pub doublet_mass_tol_ppm: f64,
    pub isotope_offsets: RangeInclusive<i8>,
    pub min_precursor_charge: u8,
    pub min_peptide_mass: f64,
    pub top_per_mass: usize,
    pub oracle_spectra: usize,
    pub open_top_k: usize,
    pub open_min_votes: u16,
    pub stub_weight: f32,
    pub secondary_site_penalty: f32,
    /// A spectrum without a precursor charge keeps the charge of its best pair score instead
    /// of its best combined score: set when several linkers are searched, whose results are
    /// merged by the pair score.
    pub charge_by_pair_score: bool,
}

pub const STUB_ION_WEIGHT: f32 = 3.0;

pub const SECONDARY_SITE_PENALTY: f32 = 12.0;

pub const PAIR_DELTA_WEIGHT: f32 = 0.3;
/// A crosslink precursor m/z moves to the most intense MS1 peak within this many ppm when that
/// peak is at least `PRECURSOR_SNAP_RATIO` times the reported one: an Orbitrap at times reports
/// a weak side peak 15-25 ppm off an intense ion as the precursor.
pub const PRECURSOR_SNAP_PPM: f64 = 25.0;
pub const PRECURSOR_SNAP_RATIO: f32 = 3.0;

pub const LOW_CHARGE_HALF_PAIR_PENALTY: f32 = 15.0;
/// An MS3 scan matches a stub peptide when its precursor m/z is within this many Th of the
/// stub peptide's m/z: a Thermo MS3 precursor is the MSn filter's target, given to 0.01.
pub const MS3_PRECURSOR_MZ_TOL: f64 = 0.01;

impl XlConfig {
    pub fn new(linker: Crosslinker) -> Self {
        XlConfig {
            linker,
            doublet: DoubletParams::default(),
            doublet_mass_tol_ppm: 10.0,
            isotope_offsets: 0..=2,
            min_precursor_charge: 3,
            min_peptide_mass: 400.0,
            top_per_mass: 4,
            oracle_spectra: 0,
            open_top_k: 50,
            open_min_votes: crate::search_params::FRAGMENT_INDEX_MIN_MATCHED,
            stub_weight: STUB_ION_WEIGHT,
            secondary_site_penalty: SECONDARY_SITE_PENALTY,
            charge_by_pair_score: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct XlPeptide {
    pub peptide: Peptide,
    pub proteins: Vec<u32>,
    pub occurrences: Vec<(u32, u32)>,
    pub n_term_occurrences: Vec<(u32, u32)>,
    pub c_term_occurrences: Vec<(u32, u32)>,
    pub is_decoy: bool,
    pub sites: Vec<LinkSite>,
    pub sites_end2: Option<Vec<LinkSite>>,
}

impl XlPeptide {
    pub fn bare(&self) -> Vec<u8> {
        self.peptide.residues.iter().map(|aa| aa.residue).collect()
    }

    pub fn occurrences_at(&self, site: &LinkSite) -> &[(u32, u32)] {
        if site.n_term {
            &self.n_term_occurrences
        } else if site.c_term {
            &self.c_term_occurrences
        } else {
            &self.occurrences
        }
    }

    pub fn sites_of_end(&self, end: usize) -> &[LinkSite] {
        match (end, &self.sites_end2) {
            (1, Some(v)) => v,
            _ => &self.sites,
        }
    }
}

fn end_assignments(linker: &Crosslinker) -> &'static [(usize, usize)] {
    if linker.is_homobifunctional() {
        &[(0, 0)]
    } else {
        &[(0, 1), (1, 0)]
    }
}

pub struct OpenIndex {
    index: ChunkFragmentIndex,
    row_of: HashMap<String, u32>,
}

pub struct XlCandidates {
    peptides: Vec<XlPeptide>,
    index: MassIndex,
    open: Option<OpenIndex>,
}

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
    pub fn from_out_of_core(prepared: &PreparedSearch<'_>, linker: &Crosslinker) -> Self {
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

    pub fn peptides(&self) -> &[XlPeptide] {
        &self.peptides
    }

    pub fn index(&self) -> &MassIndex {
        &self.index
    }

    pub fn has_open_index(&self) -> bool {
        self.open.is_some()
    }

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairKind {
    Full,
    Half,
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

#[derive(Debug, Clone, Copy)]
struct MassPair {
    a_mass: f64,
    a_tol: f64,
    b_mass: f64,
    b_tol: f64,
    kind: PairKind,
    iso: i8,
    b_stub_ions: u8,
}

#[derive(Debug, Clone, Copy)]
pub struct PeptideHit {
    pub pep: u32,
    pub site: LinkSite,
    pub form: LinkForm,
    pub score: f32,
    pub form_scores: [f32; N_FORMS],
    pub n_sites: u8,
    pub gap_next: Option<f32>,
    pub mass_err_ppm: f64,
    /// The best score of an MS3 scan of the spectrum on this peptide with its cleaved stub,
    /// already part of `score`; `None` when no MS3 precursor matches the stub peptide.
    pub ms3: Option<f32>,
}

#[derive(Debug, Clone, Default)]
pub struct DoubletSummary {
    pub n_peaks: usize,
    pub base_peak: f64,
    pub tic: f64,
    pub n_doublets: usize,
    pub n_doublets_iso: usize,
    pub masses: Vec<DoubletMass>,
    pub top_rel_intensity: f64,
    pub n_full_pairs: usize,
    pub n_half_pairs: usize,
    pub probe: [[[f32; 2]; 2]; 3],
    /// How many of `masses` (the most intense) were paired; see `paired_masses`.
    pub n_paired_masses: usize,
}

/// The doublet masses of a spectrum worth pairing: the most intense ones, as many as the
/// doublets in excess of chance. The chance count is the number of doublet masses found with
/// the heavy stub shifted off its true mass (`probe[0][0][1]`, the mean over
/// `PROBE_NULL_SHIFTS`); a spectrum of a few real doublets has almost none, while one with
/// thousands of noise peaks finds hundreds of doublets at any spacing. At least
/// `MIN_PAIRED_DOUBLET_MASSES` are kept. `masses` is sorted by intensity, highest first.
pub fn paired_masses(d: &DoubletSummary) -> &[DoubletMass] {
    let [real, chance] = d.probe[0][0];
    let excess = (real - chance).max(0.0).ceil() as usize;
    let keep = excess.max(MIN_PAIRED_DOUBLET_MASSES).min(d.masses.len());
    &d.masses[..keep]
}

pub const PROBE_FLOORS: [f64; 3] = [0.0, 0.012, 0.05];
/// The fewest doublet masses of a spectrum that are paired, however many of its doublets the
/// chance estimate explains.
pub const MIN_PAIRED_DOUBLET_MASSES: usize = 16;
pub const PROBE_NULL_SHIFTS: [f64; 2] = [-5.0, 5.0];

#[derive(Debug, Clone)]
pub struct Csm {
    pub alpha: PeptideHit,
    pub beta: PeptideHit,
    /// The peptide pair as its two bare sequences, sorted.
    pub pair_key: PairKey,
    pub pair_score: f32,
    pub second_pair_score: Option<f32>,
    pub kind: PairKind,
    pub iso: i8,
    pub class: PairClass,
    pub intra: bool,
    pub homotypic: bool,
    pub precursor_err_ppm: f64,
    pub b_stub_ions: u8,
    pub n_mass_pairs: usize,
    pub alpha_doublet: Option<DoubletMass>,
    pub beta_doublet: Option<DoubletMass>,
    pub rescores: Rescores,
    pub charge: u8,
}

/// A peptide pair as its two bare sequences, sorted: pairs that differ only in
/// modifications or link sites are the same pair.
pub type PairKey = [Vec<u8>; 2];

/// The best pair score among `pairs` of a peptide pair other than `winner`.
pub fn runner_up(pairs: &[(PairKey, f32)], winner: &PairKey) -> Option<f32> {
    pairs
        .iter()
        .filter(|(k, _)| k != winner)
        .map(|(_, s)| *s)
        .max_by(f32::total_cmp)
}

pub fn combined_pair_score(pair: f32, second: Option<f32>) -> f32 {
    let delta = second.map_or(0.0, |s| (pair - s).max(0.0));
    pair + PAIR_DELTA_WEIGHT * delta
}

pub fn half_pair_penalty(kind: PairKind, charge: u8) -> f32 {
    if kind == PairKind::Half && charge <= 3 {
        LOW_CHARGE_HALF_PAIR_PENALTY
    } else {
        0.0
    }
}

impl Csm {
    pub fn score(&self) -> f32 {
        combined_pair_score(self.pair_score, self.second_pair_score)
            - half_pair_penalty(self.kind, self.charge)
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Rescores {
    pub alpha_no_stubs: f32,
    pub beta_no_stubs: f32,
    pub alpha_decoupled: f32,
    pub beta_decoupled: f32,
}

#[derive(Debug, Clone)]
pub struct OracleResult {
    pub n_pairs: usize,
    pub csm_pair: Option<(bool, f32)>,
    pub best: Option<(u32, u32, f32)>,
    pub class: Option<PairClass>,
    pub doublet_support: (bool, bool),
}

#[derive(Debug, Clone)]
pub struct XlSpectrumResult {
    pub spec_idx: usize,
    pub charge: u8,
    pub precursor_mass: f64,
    pub doublets: DoubletSummary,
    pub csm: Option<Csm>,
    pub oracle: Option<OracleResult>,
    pub linker: usize,
    /// The best and the next-best peptide pair at every charge searched, with their pair
    /// scores: the candidates of the CSM's runner-up across charges and linkers.
    pub pairs: Vec<(PairKey, f32)>,
}

struct StubMods {
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
    peaks: &'a [(f64, f32)],
    stub_weight: f32,
    stub_tol_ppm: f64,
    secondary_site_penalty: f32,
    ms3: &'a [Ms3Ctx<'a>],
}

/// One MS3 scan of the spectrum, scored at every charge a stub peptide of the precursor can
/// carry (1 to the precursor charge - 1).
struct Ms3Ctx<'a> {
    mz: f64,
    by_charge: Vec<(u8, ScoredSpectrum<'a>)>,
}

fn ms3_spectrum(p: &ProductScan) -> Spectrum {
    let mut peaks = p.peaks.clone();
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    Spectrum {
        title: p.id.clone(),
        precursor_mz: p.precursor_mz,
        scan: p.scan,
        rt_seconds: p.rt_seconds,
        peaks,
        ..Default::default()
    }
}

fn ms3_contexts<'a>(specs: &'a [Spectrum], scorer: &'a RankScorer, z: u8) -> Vec<Ms3Ctx<'a>> {
    specs
        .iter()
        .filter(|s| s.precursor_mz > 0.0 && !s.peaks.is_empty())
        .map(|s| Ms3Ctx {
            mz: s.precursor_mz,
            by_charge: (1..z.max(2))
                .map(|q| (q, ScoredSpectrum::new(s, scorer, q)))
                .collect(),
        })
        .collect()
}

/// Whether an MS3 precursor at `mz` is a stub peptide of neutral mass `mass` at `charge`,
/// on any isotope peak the instrument may have isolated (`max_isotope`, from
/// `model::isotope::max_isolated_isotope`).
fn ms3_precursor_matches(mz: f64, charge: u8, mass: f64, max_isotope: usize) -> bool {
    model::isotope::isolated_isotope_peak(mz, mass, charge, max_isotope, MS3_PRECURSOR_MZ_TOL)
        .is_some()
}

/// The best score of an MS3 scan whose precursor is `stub_peptide` at some charge.
fn ms3_score(ctx: &SpecCtx<'_>, stub_peptide: &Peptide) -> Option<f32> {
    let m = stub_peptide.mass();
    let max_isotope = model::isotope::max_isolated_isotope(m);
    let mut best: Option<f32> = None;
    for c in ctx.ms3 {
        for (q, ss) in &c.by_charge {
            if !ms3_precursor_matches(c.mz, *q, m, max_isotope) {
                continue;
            }
            let s = score_psm(ss, stub_peptide, ctx.scorer, *q, ctx.frag_tol)
                + psm_edge_score(ss, stub_peptide, ctx.scorer, *q) as f32;
            best = Some(best.map_or(s, |b| b.max(s)));
        }
    }
    best
}

fn site_penalty(
    penalty: f32,
    p: &Peptide,
    site: &LinkSite,
    linker: &Crosslinker,
    end: usize,
) -> f32 {
    let e = &linker.ends[end];
    let selective_amine = e.reacts_with(b'K') && !e.reacts_with(b'G');
    if penalty == 0.0 || site.n_term || site.c_term || !selective_amine {
        return 0.0;
    }
    match p.residues.get(site.pos).map(|aa| aa.residue) {
        Some(b'S' | b'T' | b'Y') => penalty,
        _ => 0.0,
    }
}

fn stub_bonus(ctx: &SpecCtx<'_>, mass: f64, linker: &Crosslinker) -> f32 {
    if ctx.stub_weight == 0.0 {
        return 0.0;
    }
    let Some(c) = linker.cleavage.as_ref() else {
        return 0.0;
    };
    let n = stub_ions_present(
        ctx.peaks,
        mass,
        c,
        ctx.z.saturating_sub(1).max(1),
        ctx.stub_tol_ppm,
    )
    .iter()
    .filter(|&&v| v)
    .count();
    ctx.stub_weight * n as f32
}

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
    let mut best: Option<(f32, LinkSite, LinkForm, Option<f32>)> = None;
    let mut form_scores = [f32::NEG_INFINITY; N_FORMS];
    for &site in sites {
        for form in linker.forms() {
            let Some(delta) = form.delta(linker, partner_mass) else {
                continue;
            };
            let scored = with_link_mass(&p.peptide, site.pos, form, delta, mods);
            let ms3 = match form {
                LinkForm::Stub(_) if !ctx.ms3.is_empty() => ms3_score(ctx, &scored),
                _ => None,
            };
            let s = score_psm(&ctx.ss, &scored, ctx.scorer, ctx.z, ctx.frag_tol)
                + psm_edge_score(&ctx.ss, &scored, ctx.scorer, ctx.z) as f32
                - site_penalty(ctx.secondary_site_penalty, &p.peptide, &site, linker, end)
                + ms3.map_or(0.0, |v| v.max(0.0));
            let fi = form.index();
            if s > form_scores[fi] {
                form_scores[fi] = s;
            }
            if best.is_none_or(|(b, _, _, _)| s > b) {
                best = Some((s, site, form, ms3));
            }
        }
    }
    let (score, site, form, ms3) = best?;
    let score = score + stub_bonus(ctx, p.peptide.mass(), linker);
    Some(PeptideHit {
        pep: pep_idx,
        site,
        form,
        score,
        form_scores,
        n_sites: sites.len().min(255) as u8,
        gap_next: None,
        mass_err_ppm: 0.0,
        ms3,
    })
}

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
        let count = |c: &xl::crosslinker::Cleavage| {
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

#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_arguments)]
fn search_spectrum(
    spec_idx: usize,
    spec: &Spectrum,
    ms3: &[ProductScan],
    prepared: &PreparedSearch<'_>,
    table: &XlCandidates,
    cfg: &XlConfig,
    mods: &StubMods,
    run_oracle: bool,
) -> Option<XlSpectrumResult> {
    let ms3: Vec<Spectrum> = if cfg.linker.cleavage.is_some() {
        ms3.iter().map(ms3_spectrum).collect()
    } else {
        Vec::new()
    };
    let ms3 = ms3.as_slice();
    match spec.precursor_charge {
        Some(z) => {
            let z = u8::try_from(z).ok()?;
            if z < cfg.min_precursor_charge {
                return None;
            }
            search_spectrum_at(
                spec_idx, spec, ms3, z, prepared, table, cfg, mods, run_oracle,
            )
        }
        None => {
            let range = &prepared.params.charge_range;
            let lo = (*range.start()).max(cfg.min_precursor_charge);
            choose_charge(
                (lo..=*range.end()).filter_map(|z| {
                    search_spectrum_at(spec_idx, spec, ms3, z, prepared, table, cfg, mods, false)
                }),
                cfg.charge_by_pair_score,
            )
        }
    }
}

/// The result of a spectrum without a precursor charge: the charge with the best combined
/// pair score (or with `by_pair_score`, the best pair score), which keeps its own runner-up.
/// The peptide pairs of every charge stay with it; a search over several linkers takes the
/// runner-up from them once the linkers are merged.
fn choose_charge(
    trials: impl Iterator<Item = XlSpectrumResult>,
    by_pair_score: bool,
) -> Option<XlSpectrumResult> {
    let score = |r: &XlSpectrumResult| {
        r.csm.as_ref().map(|c| {
            if by_pair_score {
                c.pair_score
            } else {
                combined_pair_score(c.pair_score, c.second_pair_score)
            }
        })
    };
    let mut best: Option<XlSpectrumResult> = None;
    let mut pairs: Vec<(PairKey, f32)> = Vec::new();
    for r in trials {
        let better = match &best {
            None => true,
            Some(b) => match (score(&r), score(b)) {
                (Some(x), Some(y)) => x > y,
                (Some(_), None) => true,
                _ => false,
            },
        };
        pairs.extend(r.pairs.iter().cloned());
        if better {
            best = Some(r);
        }
    }
    if let Some(b) = best.as_mut() {
        b.pairs = pairs;
    }
    best
}

#[allow(clippy::too_many_arguments)]
fn search_spectrum_at(
    spec_idx: usize,
    spec: &Spectrum,
    ms3: &[Spectrum],
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
    let paired = paired_masses(&doublets);
    let n_paired = paired.len();
    let pairs = mass_pairs(&peaks, paired, precursor_mass, prec_tol, z, cfg);
    doublets.n_paired_masses = n_paired;
    doublets.n_full_pairs = pairs.iter().filter(|p| p.kind == PairKind::Full).count();
    doublets.n_half_pairs = pairs.len() - doublets.n_full_pairs;

    let linker = &cfg.linker;
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
        linker: 0,
        pairs: Vec::new(),
    };
    if pairs.is_empty() && open_rows.is_empty() && !run_oracle {
        return Some(result);
    }
    let ms3_ctx = ms3_contexts(ms3, prepared.scorer, z);
    let ctx = SpecCtx {
        ss: ScoredSpectrum::new(spec, prepared.scorer, z),
        scorer: prepared.scorer,
        z,
        frag_tol: prepared.fragment_tolerance_da,
        peaks: &peaks,
        stub_weight: cfg.stub_weight,
        stub_tol_ppm: cfg.doublet.tol_ppm,
        secondary_site_penalty: cfg.secondary_site_penalty,
        ms3: &ms3_ctx,
    };

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
            .map(|(_, o)| o)
            .max_by(|x, y| x.score.total_cmp(&y.score));
        result.pairs.push((best_seqs.clone(), c.score));
        if let Some(o) = second {
            result.pairs.push((seq_pair(o), o.score));
        }
        let second = second.map(|o| o.score);
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
            pair_key: best_seqs,
            pair_score: c.score,
            second_pair_score: second,
            kind: c.pair.kind,
            iso: c.pair.iso,
            class: PairClass::of(pa.is_decoy, pb.is_decoy),
            intra: is_intra(prepared, pa, pb),
            homotypic: pa.bare() == pb.bare(),
            precursor_err_ppm: (precursor_mass - theo) / theo * 1e6,
            b_stub_ions: c.pair.b_stub_ions,
            n_mass_pairs: pairs.len() + open_rows.len(),
            alpha_doublet: doublet_at(&result.doublets.masses, pa.peptide.mass(), cfg),
            beta_doublet: doublet_at(&result.doublets.masses, pb.peptide.mass(), cfg),
            rescores: rescores(spec, prepared, z, table, &alpha, &beta, linker, mods),
            charge: z,
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

/// Search every spectrum; `ms3[i]` holds the MS3 scans of spectrum `i` (missing or empty
/// when it has none).
pub fn xl_search_run(
    spectra: &[Spectrum],
    ms3: &[Vec<ProductScan>],
    prepared: &PreparedSearch<'_>,
    table: &XlCandidates,
    cfg: &XlConfig,
) -> Vec<XlSpectrumResult> {
    let mods = StubMods::new(&cfg.linker);
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
            let children = ms3.get(i).map_or(&[][..], Vec::as_slice);
            search_spectrum(i, spec, children, prepared, table, cfg, &mods, run_oracle)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::isotope::max_isolated_isotope;

    pub(crate) fn linker(spec: &str) -> Crosslinker {
        xl::spec::XlDecl::parse(spec)
            .and_then(|d| d.resolve(model::unimod::bundled()))
            .unwrap()
    }

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'K', b'A')
    }

    #[test]
    fn a_peptide_of_two_proteins_lists_both_whatever_its_flanks() {
        use model::aa_set::AminoAcidSetBuilder;
        use model::protein::{Protein, ProteinDb};
        let prot = |acc: &str, seq: &[u8]| Protein {
            accession: acc.into(),
            description: String::new(),
            sequence: seq.to_vec(),
            annotation: None,
        };
        let db = ProteinDb {
            proteins: vec![
                prot("CETN2", b"MSEKAFKLFDDDETGKISFKNLK"),
                prot("CETN1", b"MSEKAFRLFDDDETGKISFKNLK"),
            ],
        };
        let idx = crate::search_index::SearchIndex::from_target_db(&db, "XXX_");
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        let params = crate::search_params::SearchParams::default_tryptic(aa);
        let cands: Vec<Candidate> =
            crate::candidate_gen::enumerate_candidates(&idx, &params, "XXX_").collect();
        let table = XlCandidates::from_candidates(cands, &linker("DSSO"));
        let shared = |t: &XlCandidates| {
            t.peptides()
                .iter()
                .find(|p| p.bare() == b"LFDDDETGKISFK" && !p.is_decoy)
                .map(|p| p.proteins.clone())
                .unwrap()
        };
        assert_eq!(shared(&table).len(), 2);
        let scorer = crate::coisolation::tests::tiny_scorer();
        let prepared =
            PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX_").unwrap();
        let table = XlCandidates::from_out_of_core(&prepared, &linker("DSSO"));
        assert_eq!(shared(&table).len(), 2, "{:?}", shared(&table));
    }

    #[test]
    fn only_the_doublets_in_excess_of_chance_are_paired() {
        let mass = |m: f64, i: f64| DoubletMass {
            mass: m,
            intensity: i,
            n_charges: 1,
            isotope_confirmed: false,
            min_charge: 2,
        };
        let summary = |n: usize, chance: f32| {
            let mut d = DoubletSummary {
                masses: (0..n)
                    .map(|k| mass(500.0 + k as f64, (n - k) as f64))
                    .collect(),
                ..Default::default()
            };
            d.probe[0][0] = [n as f32, chance];
            d
        };
        let clean = summary(30, 0.0);
        assert_eq!(paired_masses(&clean).len(), 30);
        let noisy = summary(540, 500.5);
        let kept = paired_masses(&noisy);
        assert_eq!(kept.len(), 40);
        assert!(kept
            .iter()
            .all(|m| m.intensity >= noisy.masses[39].intensity));
        assert_eq!(
            paired_masses(&summary(300, 299.0)).len(),
            MIN_PAIRED_DOUBLET_MASSES
        );
        assert_eq!(paired_masses(&summary(5, 9.0)).len(), 5);
    }

    #[test]
    fn an_ms3_precursor_matches_its_stub_peptide_at_its_charge() {
        let dsso = linker("DSSO");
        let alkene = dsso.cleavage.as_ref().unwrap().stubs[0].mass;
        let stub_peptide = 1000.5 + alkene;
        let mz2 = stub_peptide / 2.0 + PROTON;
        assert!(ms3_precursor_matches(
            mz2,
            2,
            stub_peptide,
            max_isolated_isotope(stub_peptide)
        ));
        assert!(ms3_precursor_matches(
            mz2 + 0.009,
            2,
            stub_peptide,
            max_isolated_isotope(stub_peptide)
        ));
        assert!(!ms3_precursor_matches(
            mz2 + 0.02,
            2,
            stub_peptide,
            max_isolated_isotope(stub_peptide)
        ));
        assert!(!ms3_precursor_matches(
            mz2,
            3,
            stub_peptide,
            max_isolated_isotope(stub_peptide)
        ));
        assert!(!ms3_precursor_matches(
            mz2,
            0,
            stub_peptide,
            max_isolated_isotope(stub_peptide)
        ));
        let heavy = 2496.2496 + alkene;
        let m1 = (heavy + ISOTOPE_SPACING) / 2.0 + PROTON;
        assert!(
            ms3_precursor_matches(m1, 2, heavy, max_isolated_isotope(heavy)),
            "M+1 of a 2.5 kDa stub peptide"
        );
        assert!(!ms3_precursor_matches(
            m1 + 0.02,
            2,
            heavy,
            max_isolated_isotope(heavy)
        ));
        let light = 600.0 + alkene;
        let far = (light + 3.0 * ISOTOPE_SPACING) / 2.0 + PROTON;
        assert!(
            !ms3_precursor_matches(far, 2, light, max_isolated_isotope(light)),
            "M+3 of a small stub peptide"
        );
        let s = ms3_spectrum(&ProductScan {
            id: "scan=9".into(),
            precursor_mz: mz2,
            peaks: vec![(300.0, 1.0), (200.0, 2.0)],
            ..Default::default()
        });
        assert_eq!(s.peaks, vec![(200.0, 2.0), (300.0, 1.0)]);
        assert_eq!(s.precursor_charge, None);
    }

    #[test]
    fn by_ions_and_peak_removal() {
        let p = pep(b"GAK");
        let ions = by_ions(&p, 1);
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
            mobility_window: None,
        };
        let kept = without_peaks(&spec, &ions, 0.01);
        assert_eq!(kept.peaks, vec![(500.0, 2.0)]);
    }

    #[test]
    fn the_csm_score_adds_a_share_of_the_gap_to_the_next_pair() {
        assert_eq!(combined_pair_score(40.0, None), 40.0);
        assert!(
            (combined_pair_score(40.0, Some(30.0)) - (40.0 + PAIR_DELTA_WEIGHT * 10.0)).abs()
                < 1e-6
        );
        assert_eq!(combined_pair_score(40.0, Some(45.0)), 40.0);
        assert!(combined_pair_score(40.0, Some(10.0)) > combined_pair_score(40.0, Some(39.0)));
    }

    fn trial(z: u8, pair: &str, score: f32, second: Option<(&str, f32)>) -> XlSpectrumResult {
        let key = |p: &str| [p.as_bytes().to_vec(), b"PEPTIDEK".to_vec()];
        let hit = PeptideHit {
            pep: 0,
            site: xl::sites::LinkSite {
                pos: 0,
                n_term: false,
                c_term: false,
            },
            form: xl::crosslinker::LinkForm::Intact,
            score,
            form_scores: Default::default(),
            n_sites: 1,
            gap_next: None,
            mass_err_ppm: 0.0,
            ms3: None,
        };
        let mut pairs = vec![(key(pair), score)];
        pairs.extend(second.map(|(p, s)| (key(p), s)));
        XlSpectrumResult {
            spec_idx: 0,
            charge: z,
            precursor_mass: 0.0,
            doublets: Default::default(),
            csm: Some(Csm {
                alpha: hit,
                beta: hit,
                pair_key: key(pair),
                pair_score: score,
                second_pair_score: second.map(|(_, s)| s),
                kind: PairKind::Full,
                iso: 0,
                class: PairClass::TT,
                intra: true,
                homotypic: false,
                precursor_err_ppm: 0.0,
                b_stub_ions: 0,
                n_mass_pairs: 1,
                alpha_doublet: None,
                beta_doublet: None,
                rescores: Default::default(),
                charge: z,
            }),
            oracle: None,
            linker: 0,
            pairs,
        }
    }

    #[test]
    fn an_uncharged_spectrum_keeps_the_runner_up_of_its_chosen_charge() {
        // Charge 3 (100 / 50, combined 115) beats charge 4 (110 / 109, combined 110.3).
        let trials = || {
            [
                trial(3, "AK", 100.0, Some(("BK", 50.0))),
                trial(4, "CK", 110.0, Some(("DK", 109.0))),
            ]
            .into_iter()
        };
        // Several linkers: the charge of the best pair score (110 at charge 4).
        let r = choose_charge(trials(), true).unwrap();
        assert_eq!((r.charge, r.csm.as_ref().unwrap().pair_score), (4, 110.0));
        let r = choose_charge(trials(), false).unwrap();
        let c = r.csm.as_ref().unwrap();
        assert_eq!(
            (r.charge, c.pair_score, c.second_pair_score),
            (3, 100.0, Some(50.0))
        );
        // Every charge's pairs stay for a merge over several linkers.
        assert_eq!(r.pairs.len(), 4);
        assert_eq!(runner_up(&r.pairs, &c.pair_key), Some(110.0));
    }

    #[test]
    fn a_half_pair_of_a_triply_charged_precursor_pays_a_penalty() {
        assert_eq!(
            half_pair_penalty(PairKind::Half, 3),
            LOW_CHARGE_HALF_PAIR_PENALTY
        );
        assert_eq!(half_pair_penalty(PairKind::Half, 4), 0.0);
        assert_eq!(half_pair_penalty(PairKind::Full, 3), 0.0);
        assert_eq!(half_pair_penalty(PairKind::Open, 3), 0.0);
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
        let r = with_link_mass(&p, 3, LinkForm::Stub(0), stub, &mods);
        assert!(std::ptr::eq(
            q.residues[3].mod_.unwrap(),
            r.residues[3].mod_.unwrap()
        ));
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
    fn defaults_search_two_missing_isotopes_and_credit_stub_ions() {
        let cfg = XlConfig::new(Crosslinker::dsso());
        assert_eq!(cfg.isotope_offsets, 0..=2);
        assert_eq!(cfg.stub_weight, STUB_ION_WEIGHT);
        assert_eq!(cfg.secondary_site_penalty, SECONDARY_SITE_PENALTY);
    }

    #[test]
    fn s_t_y_sites_of_an_nhs_ester_pay_the_secondary_site_penalty() {
        let x = crate::xl_search::tests::linker("DSSO (K, S, T, Y, Protein N-term)");
        let p = pep(b"SAKTY");
        let at = |pos: usize, n_term: bool| LinkSite {
            pos,
            n_term,
            c_term: false,
        };
        assert_eq!(site_penalty(12.0, &p, &at(2, false), &x, 0), 0.0);
        assert_eq!(site_penalty(12.0, &p, &at(0, false), &x, 0), 12.0);
        assert_eq!(site_penalty(12.0, &p, &at(3, false), &x, 0), 12.0);
        assert_eq!(site_penalty(12.0, &p, &at(0, true), &x, 0), 0.0);
        assert_eq!(site_penalty(0.0, &p, &at(0, false), &x, 0), 0.0);
        let sda = crate::xl_search::tests::linker("SDA");
        let any_end = (0..2).find(|&e| sda.ends[e].reacts_with(b'G')).unwrap();
        assert_eq!(site_penalty(12.0, &p, &at(0, false), &sda, any_end), 0.0);
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
        let full = mass_pairs(&[], &[dm(ma), dm(mb)], prec, 0.01, 4, &cfg);
        assert_eq!(full.len(), 1);
        assert_eq!(full[0].kind, PairKind::Full);
        assert_eq!(full[0].b_stub_ions, 3);
        assert!((full[0].a_mass - ma).abs() < 1e-9 && (full[0].b_mass - mb).abs() < 1e-9);
        assert!(mass_pairs(&[], &[dm(ma)], prec, 0.01, 4, &cfg).is_empty());
        let stub_peak = [((mb + thiol) / 2.0 + PROTON, 10.0_f32)];
        let half = mass_pairs(&stub_peak, &[dm(ma)], prec, 0.01, 4, &cfg);
        assert_eq!(half.len(), 1);
        assert_eq!(half[0].kind, PairKind::Half);
        assert_eq!(half[0].iso, 0);
        assert_eq!(half[0].b_stub_ions, 1);
    }
}
