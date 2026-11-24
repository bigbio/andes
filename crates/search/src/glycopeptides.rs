pub mod fragments;
mod scan;

use rustc_hash::{FxHashMap, FxHashSet};

use model::mass::{H2O, ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use rayon::prelude::*;

use glyco::backbone::{
    core_y_intensity, count_core_y_hits, glycan_y_hit_frac, glycan_y_intensity,
    partial_glycan_by_intensity, y0y1_anchor_intensity, SpectrumStats,
};
use glyco::glycan::db::{GlycanComp, DENOVO_TOPOLOGY_ID};
use glyco::glycan::first::{search_glycans, Glycan, GlycanConfig, GlycanCore, GlycanIonIndex};
use glyco::psm::{
    glyco_gp_fused_score_iso, GlycoPsmKey, GLYCO_GP_CZ_DEFAULT, GLYCO_GP_H_DEFAULT,
    GLYCO_GP_ISO_DEFAULT, GLYCO_GP_J_DEFAULT, GLYCO_GP_K_DEFAULT,
};

pub const GLYCO_BACKBONE_TOP_K: usize = 150;

pub const GLYCO_PF_CHARGE: u8 = 2;

pub const GLYCO_MAX_PF: usize = 1024;

#[derive(Clone, Debug, Default)]
pub struct GlycoConfig {
    pub core: GlycanCore,
    pub min_core_y: u32,
    pub scan_filter_path: Option<std::path::PathBuf>,
    pub isotope_error_override: Option<std::ops::RangeInclusive<i8>>,
    pub isotope_error_override_mask: Option<Vec<bool>>,
    pub max_gen_peaks: usize,
    pub retrieval_tol_ppm: Option<f64>,
    pub full_glycan_db: bool,
    pub debug: bool,
    pub hcd_pair: bool,
}

use glyco::hybrid::{db_branch, BackboneHit, Source};
use glyco::oxonium::{
    hex_diagnostic_evidence, oxonium_gate, sialic_consistency, OXONIUM_GATE_MIN_FRAC,
};

use crate::glycopeptides::fragments::FragmentIndex;
use scan::{retrieve_backbones, winner_spread, CheapWinner, GlycoScan};

use crate::engine::{compute_psm_features, PreparedSearch};
use crate::psm::PsmFeatures;
use crate::psm::PsmMatch;

type GlycanWinnerKey = (u32, u64);

use scoring_crate::psm::{
    candidate_rank_entropy, cz_hyperscore_psm, cz_matched_intensity_frac, cz_structure_features,
    fuse_strong_score, hyperscore_psm_with_matches, listwise_score_gap, psm_edge_score, score_psm,
    score_psm_float, scored_for_charge, strong_score_calibrated_loo, ScoredSpectrum,
    StrongScoreInputs,
};

#[allow(clippy::too_many_arguments)]
fn glyco_candidate_strong_score(
    ss: &ScoredSpectrum<'_>,
    peptide: &model::peptide::Peptide,
    spec_scorer: &scoring_crate::psm::RankScorer,
    z: u8,
    spectrum_rank_entropy: f32,
    spectrum_listwise_gap: f32,
) -> (PsmFeatures, f32) {
    let features = compute_psm_features(ss, peptide, spec_scorer, z);
    let strong = fuse_strong_score(&StrongScoreInputs {
        intensity_signal: features.intensity_signal,
        chance_match_surprise: features.chance_match_surprise,
        mass_competition_evidence: features.mass_competition_evidence,
        candidate_rank_entropy: spectrum_rank_entropy,
        listwise_score_gap: spectrum_listwise_gap,
    });
    (features, strong)
}

#[derive(Debug, Clone)]
pub struct FullGlycoPsm {
    pub glycan_key: GlycoPsmKey,
    pub psm: PsmMatch,
}

#[derive(Debug, Clone)]
pub struct GlycoSpectrumResult {
    pub spectrum_idx: usize,
    pub hits: Vec<FullGlycoPsm>,
}

const MIN_GLYCAN_N: f64 = 406.0;

const MIN_GLYCAN_O: f64 = 203.0;

fn min_glycan_floor(core: GlycanCore) -> f64 {
    match core {
        GlycanCore::NGlycan => MIN_GLYCAN_N,
        GlycanCore::OGlycan => MIN_GLYCAN_O,
    }
}

fn nearest_glycan_mass(
    sorted: &[(f64, usize)],
    glycans: &[GlycanComp],
    target: f64,
    tol: f64,
) -> Option<GlycanComp> {
    let lo = target - tol;
    let hi = target + tol;
    let start = sorted.partition_point(|&(m, _)| m < lo);
    let mut best: Option<(f64, usize)> = None;
    for &(m, gi) in &sorted[start..] {
        if m > hi {
            break;
        }
        let d = (m - target).abs();
        if best.is_none_or(|(bd, _)| d < bd) {
            best = Some((d, gi));
        }
    }
    best.map(|(_, gi)| glycans[gi].clone())
}

#[inline]
fn pf_tiebreak_hash(idx: u32) -> u32 {
    idx.wrapping_mul(0x9E37_79B1)
}

fn order_peptide_first(pf: &mut [(u32, u32)]) {
    pf.sort_unstable_by(|a, b| {
        b.1.cmp(&a.1)
            .then_with(|| pf_tiebreak_hash(a.0).cmp(&pf_tiebreak_hash(b.0)))
            .then_with(|| a.0.cmp(&b.0))
    });
}

fn glyco_aware_peptide(
    pep: &model::peptide::Peptide,
    gsite: usize,
    gmass: f64,
) -> model::peptide::Peptide {
    if gmass <= 0.0 || gsite >= pep.residues.len() {
        return pep.clone();
    }
    let mut residues = pep.residues.clone();
    if let Some(res) = residues.get_mut(gsite) {
        let base = res.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        res.mod_ = Some(model::modification::leak_mod(
            model::modification::Modification {
                name: "Glycan".to_string(),
                mass_delta: base + gmass,
                residue: model::modification::ResidueSpec::Specific(res.residue),
                location: model::modification::ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 1,
                annotated_site_mass: None,
                blocks_cleavage: false,
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
            },
        ));
    }
    model::peptide::Peptide::new(residues, pep.pre, pep.post)
}

fn glyco_site_for(pep: &model::peptide::Peptide, core: GlycanCore) -> usize {
    let res: Vec<u8> = pep.residues.iter().map(|aa| aa.residue).collect();
    if core == GlycanCore::OGlycan {
        return glyco::sequon::first_ser_thr_site(&res).unwrap_or(0);
    }
    glyco::sequon::glycosites(&res, pep.post)
        .first()
        .copied()
        .unwrap_or(0)
}

fn cz_score_at_site(
    ss: &ScoredSpectrum<'_>,
    pep: &model::peptide::Peptide,
    gmass: f64,
    max_frag_charge: u8,
    tol_ppm: f64,
    core: GlycanCore,
) -> f32 {
    cz_hyperscore_psm(
        ss,
        pep,
        gmass,
        glyco_site_for(pep, core),
        max_frag_charge,
        tol_ppm,
    )
}

fn dedup_backbone_hits(mut all_backbone: Vec<BackboneHit>, tol_ppm: f64) -> Vec<BackboneHit> {
    if all_backbone.is_empty() {
        return Vec::new();
    }
    all_backbone.sort_by(|a, b| {
        a.backbone_mass
            .total_cmp(&b.backbone_mass)
            .then_with(|| {
                let oa = if a.source == Source::Db { 0u8 } else { 1u8 };
                let ob = if b.source == Source::Db { 0u8 } else { 1u8 };
                oa.cmp(&ob)
            })
            .then_with(|| a.isotope_offset.abs().cmp(&b.isotope_offset.abs()))
            .then_with(|| a.charge.cmp(&b.charge))
    });

    let mut deduped: Vec<BackboneHit> = Vec::with_capacity(all_backbone.len());
    let mut rep = all_backbone.remove(0);
    for next in all_backbone {
        let tol = (rep.backbone_mass * tol_ppm * 1e-6_f64).max(0.02);
        let same_backbone = (next.backbone_mass - rep.backbone_mass).abs() < tol;
        let same_hypothesis = match (&rep.glycan, &next.glycan) {
            (Some(g1), Some(g2)) => g1 == g2,
            (None, None) => rep.isotope_offset == next.isotope_offset,
            _ => rep.isotope_offset == next.isotope_offset,
        };
        if same_backbone && same_hypothesis {
            let peptide_first = rep.peptide_first || next.peptide_first;
            if rep.source == Source::DeNovo && next.source == Source::Db {
                rep = next;
            }
            rep.peptide_first = peptide_first;
        } else {
            deduped.push(rep);
            rep = next;
        }
    }
    deduped.push(rep);
    deduped
}

pub struct GlycoCandidates {
    candidates: Vec<crate::candidates::Candidate>,
    mass_order: Vec<u32>,
}

impl GlycoCandidates {
    pub fn from_index(prepared: &PreparedSearch<'_>, core: GlycanCore) -> Self {
        Self::from_candidates(
            prepared.candidates_where(move |residues, post| match core {
                GlycanCore::NGlycan => crate::candidates::span_has_nxst_sequon(residues, post),
                GlycanCore::OGlycan => crate::candidates::span_has_ser_or_thr(residues, post),
            }),
            core,
        )
    }

    pub fn from_candidates(
        candidates: Vec<crate::candidates::Candidate>,
        core: GlycanCore,
    ) -> Self {
        debug_assert!(candidates.iter().all(match core {
            GlycanCore::NGlycan => crate::candidates::candidate_has_nxst_sequon,
            GlycanCore::OGlycan => crate::candidates::candidate_has_ser_or_thr,
        }));
        let mass_order = crate::engine::build_mass_order(&candidates);
        Self {
            candidates,
            mass_order,
        }
    }

    pub fn candidates(&self) -> &[crate::candidates::Candidate] {
        &self.candidates
    }
}

pub struct GlycoScoreCtx<'a> {
    pub params: &'a crate::params::SearchParams,
    pub scorer: &'a scoring_crate::psm::RankScorer,
    pub core: GlycanCore,
    pub candidates: &'a [crate::candidates::Candidate],
    pub mass_order: &'a [u32],
    pub fragment_tolerance_da: f64,
    pub frag_index: &'a FragmentIndex,
    pub glycan_sorted: &'a [(f64, usize)],
    pub glycan_list: &'a [GlycanComp],
    pub glycan_ion_index: &'a GlycanIonIndex,
    pub glycan_first_cfg: &'a GlycanConfig,
    pub tol_ppm: f64,
    pub peptide_first_on: bool,
    pub full_glycan_db: bool,
    pub min_core_y: u32,
    pub max_gen_peaks: usize,
    pub hcd_pair_on: bool,
    pub features_collapse: bool,
    pub features_enumerated: bool,
    pub scan_filter: Option<&'a std::collections::HashSet<i32>>,
    pub isotope_error_default: std::ops::RangeInclusive<i8>,
    pub isotope_error_override: Option<std::ops::RangeInclusive<i8>>,
    pub isotope_error_override_mask: Option<&'a [bool]>,
    pub all_spectra: &'a [Spectrum],
    pub hcd_partner: &'a [Option<usize>],
}

pub struct GlycoCtxOwned {
    frag_index: FragmentIndex,
    glycan_sorted: Vec<(f64, usize)>,
    core: GlycanCore,
    glycan_ion_index: GlycanIonIndex,
    glycan_first_cfg: GlycanConfig,
    scan_filter: Option<std::collections::HashSet<i32>>,
    isotope_error_override: Option<std::ops::RangeInclusive<i8>>,
    isotope_error_override_mask: Option<Vec<bool>>,
    peptide_first_on: bool,
    full_glycan_db: bool,
    min_core_y: u32,
    max_gen_peaks: usize,
    hcd_pair_on: bool,
    features_collapse: bool,
    features_enumerated: bool,
}

impl GlycoCtxOwned {
    pub fn build(
        candidates: &[crate::candidates::Candidate],
        glycan_list: &[GlycanComp],
        fragment_tolerance_da: f64,
        cfg: GlycoConfig,
        tol_ppm: f64,
    ) -> Self {
        let peptide_first_on = true;
        let full_glycan_db = cfg.full_glycan_db;
        let hcd_pair_on = cfg.hcd_pair;
        let features_collapse = !cfg.debug;
        let features_enumerated = !cfg.debug;
        let ctx_scan_filter_path = cfg.scan_filter_path.clone();
        let scan_filter: Option<std::collections::HashSet<i32>> = ctx_scan_filter_path
            .as_ref()
            .and_then(|path| match std::fs::read_to_string(path) {
                Ok(text) => {
                    let set: std::collections::HashSet<i32> = text
                        .lines()
                        .filter_map(|l| l.trim().parse::<i32>().ok())
                        .collect();
                    log::info!(
                        target: "glyco",
                        "restricting to {} scans from {}",
                        set.len(),
                        path.display()
                    );
                    if set.is_empty() {
                        None
                    } else {
                        Some(set)
                    }
                }
                Err(e) => {
                    log::warn!(
                        target: "glyco",
                        "--glyco-scans {} could not be read ({e}); scoring all scans",
                        path.display()
                    );
                    None
                }
            });
        let min_core_y_cfg = cfg.min_core_y;
        let max_gen_peaks = cfg.max_gen_peaks;

        let pf_charge: u8 = GLYCO_PF_CHARGE;
        let frag_index = if !peptide_first_on || full_glycan_db {
            FragmentIndex::build(
                std::iter::empty::<(u32, &model::peptide::Peptide)>(),
                fragment_tolerance_da.max(0.01),
                pf_charge,
            )
        } else {
            let seq_entries: Vec<(u32, &model::peptide::Peptide)> = candidates
                .iter()
                .enumerate()
                .map(|(i, c)| (i as u32, &c.peptide))
                .collect();
            let est_mb = seq_entries.len().saturating_mul(20 * 24) / 1_000_000;
            if seq_entries.len() > 3_000_000 {
                log::warn!(
                    target: "glyco",
                    "{} sequon candidates (~{} MB postings) — \
                     large database; consider narrowing the search space",
                    seq_entries.len(),
                    est_mb
                );
            }
            match cfg.retrieval_tol_ppm {
                Some(ppm) if ppm > 0.0 => {
                    FragmentIndex::build_ppm(seq_entries.iter().copied(), ppm, pf_charge)
                }
                _ => FragmentIndex::build(
                    seq_entries.iter().copied(),
                    fragment_tolerance_da.max(0.01),
                    pf_charge,
                ),
            }
        };
        let glycan_sorted: Vec<(f64, usize)> = {
            let mut v: Vec<(f64, usize)> = glycan_list
                .iter()
                .enumerate()
                .map(|(i, g)| (g.mass, i))
                .collect();
            v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            v
        };
        let glycan_ion_index = {
            let glycans: Vec<Glycan> = glycan_list
                .iter()
                .enumerate()
                .map(|(i, g)| Glycan {
                    id: i as u32,
                    composition: g.clone(),
                })
                .collect();
            GlycanIonIndex::build(&glycans, 0.02, Some(tol_ppm), cfg.core)
        };
        let glycan_first_cfg = GlycanConfig {
            tol_ppm,
            min_diagnostic_ions: 2,
            min_core_ions: 1,
            top_k: 100,
            diagnostic_ions: match cfg.core {
                GlycanCore::NGlycan => glyco::glycan::first::default_diagnostic_ions(),
                GlycanCore::OGlycan => glyco::glycan::first::o_glycan_diagnostic_ions(),
            },
            diagnostic_filter: true,
            score_weights: (1.0, 1.0, 0.5),
        };
        GlycoCtxOwned {
            frag_index,
            glycan_sorted,
            core: cfg.core,
            glycan_ion_index,
            glycan_first_cfg,
            scan_filter,
            isotope_error_override: cfg.isotope_error_override.clone(),
            isotope_error_override_mask: cfg.isotope_error_override_mask.clone(),
            peptide_first_on,
            full_glycan_db,
            min_core_y: min_core_y_cfg,
            max_gen_peaks,
            hcd_pair_on,
            features_collapse,
            features_enumerated,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn as_ctx<'a>(
        &'a self,
        prepared: &'a PreparedSearch<'_>,
        table: &'a GlycoCandidates,
        glycan_list: &'a [GlycanComp],
        tol_ppm: f64,
        all_spectra: &'a [Spectrum],
        hcd_partner: &'a [Option<usize>],
    ) -> GlycoScoreCtx<'a> {
        GlycoScoreCtx {
            params: prepared.params,
            scorer: prepared.scorer,
            core: self.core,
            candidates: &table.candidates,
            mass_order: &table.mass_order,
            fragment_tolerance_da: prepared.fragment_tolerance_da,
            frag_index: &self.frag_index,
            glycan_sorted: &self.glycan_sorted,
            glycan_list,
            glycan_ion_index: &self.glycan_ion_index,
            glycan_first_cfg: &self.glycan_first_cfg,
            tol_ppm,
            peptide_first_on: self.peptide_first_on,
            full_glycan_db: self.full_glycan_db,
            min_core_y: self.min_core_y,
            max_gen_peaks: self.max_gen_peaks,
            hcd_pair_on: self.hcd_pair_on,
            features_collapse: self.features_collapse,
            features_enumerated: self.features_enumerated,
            scan_filter: self.scan_filter.as_ref(),
            isotope_error_default: prepared.params.isotope_error_range.clone(),
            isotope_error_override: self.isotope_error_override.clone(),
            isotope_error_override_mask: self.isotope_error_override_mask.as_deref(),
            all_spectra,
            hcd_partner,
        }
    }
}

fn glyco_charges_to_try(
    spec_charge: Option<i32>,
    charge_range: &std::ops::RangeInclusive<u8>,
) -> Vec<u8> {
    let mut out: Vec<u8> = match spec_charge {
        Some(z) if z > 0 => vec![z as u8],
        _ => (*charge_range.start()..=*charge_range.end()).collect(),
    };
    out.sort_unstable();
    out.dedup();
    out.retain(|&z| z > 0);
    out
}

pub fn isotope_window_for(
    spec_idx: usize,
    default: &std::ops::RangeInclusive<i8>,
    override_range: Option<&std::ops::RangeInclusive<i8>>,
    mask: Option<&[bool]>,
) -> std::ops::RangeInclusive<i8> {
    match override_range {
        Some(r) if mask.is_none_or(|m| m.get(spec_idx).copied().unwrap_or(false)) => r.clone(),
        _ => default.clone(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlycanSource {
    FullList,
    IonIndex,
    None,
}

pub fn glycan_source(full_glycan_db: bool, oxonium_fired: bool, is_etd: bool) -> GlycanSource {
    if is_etd && !oxonium_fired {
        return GlycanSource::FullList;
    }
    if !oxonium_fired {
        return GlycanSource::None;
    }
    if full_glycan_db {
        GlycanSource::FullList
    } else {
        GlycanSource::IonIndex
    }
}

fn has_backbone_evidence(hit: &BackboneHit, core_y_hits: u8) -> bool {
    hit.source != Source::Db || hit.peptide_first || core_y_hits > 0
}

fn hex_gate_passes(core: GlycanCore, glycan: Option<&GlycanComp>, hex_evidence: bool) -> bool {
    core != GlycanCore::OGlycan || glycan.is_none_or(|g| g.hex == 0) || hex_evidence
}

fn precursor_window_da(params: &crate::params::SearchParams, mass: f64) -> f64 {
    let t = &params.precursor_tolerance;
    t.left.as_da(mass).max(t.right.as_da(mass))
}

fn backbone_candidates(
    mass_order: &[u32],
    candidates: &[crate::candidates::Candidate],
    bb_residue: f64,
    tol_da: f64,
) -> Vec<usize> {
    let lo = bb_residue + H2O - tol_da - 1e-6;
    let hi = bb_residue + H2O + tol_da + 1e-6;
    let mass = |i: &u32| candidates[*i as usize].peptide.mass();
    let start = mass_order.partition_point(|i| mass(i) < lo);
    let end = start + mass_order[start..].partition_point(|i| mass(i) <= hi);
    let mut slots: Vec<usize> = mass_order[start..end].iter().map(|&i| i as usize).collect();
    slots.sort_unstable_by_key(|&i| (candidates[i].peptide.nominal_residue_mass(), i));
    slots
}

fn score_spectrum_glyco(
    spec_idx: usize,
    spec: &Spectrum,
    ctx: &GlycoScoreCtx<'_>,
) -> Option<GlycoSpectrumResult> {
    let params = ctx.params;
    let scorer = ctx.scorer;
    let tol_ppm = ctx.tol_ppm;
    let effective_top_k = GLYCO_BACKBONE_TOP_K;
    let is_etd = matches!(
        spec.activation_method,
        Some(model::activation::ActivationMethod::ETD)
    );
    let features_collapse = ctx.features_collapse;
    let scan_filter = ctx.scan_filter;

    if spec.peaks.len() < params.min_peaks as usize {
        return None;
    }
    if let Some(scans) = scan_filter {
        match spec.scan {
            Some(sc) if scans.contains(&sc) => {}
            _ => return None,
        }
    }

    let stats = SpectrumStats::new(&spec.peaks);

    let paired_hcd: Option<usize> = if ctx.hcd_pair_on && is_etd {
        ctx.hcd_partner.get(spec_idx).copied().flatten()
    } else {
        None
    };
    let gen_spec: &Spectrum = match paired_hcd {
        Some(p) => &ctx.all_spectra[p],
        None => spec,
    };
    let capped_gen_peaks: Option<Vec<(f64, f32)>> = match ctx.max_gen_peaks {
        n if n > 0 && gen_spec.peaks.len() > n => {
            let mut by_intensity: Vec<(f64, f32)> = gen_spec.peaks.clone();
            by_intensity.sort_unstable_by(|a, b| b.1.total_cmp(&a.1).then(a.0.total_cmp(&b.0)));
            by_intensity.truncate(n);
            by_intensity.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
            Some(by_intensity)
        }
        _ => None,
    };
    let gen_peaks: &[(f64, f32)] = capped_gen_peaks.as_deref().unwrap_or(&gen_spec.peaks);
    let gen_stats_owned =
        (paired_hcd.is_some() || capped_gen_peaks.is_some()).then(|| SpectrumStats::new(gen_peaks));
    let gen_stats: &SpectrumStats = gen_stats_owned.as_ref().unwrap_or(&stats);
    let (y_peaks, y_stats): (&[(f64, f32)], &SpectrumStats) = (&spec.peaks, &stats);

    let spec_scorer: &scoring_crate::psm::RankScorer = scorer;
    let gen_scorer: &scoring_crate::psm::RankScorer = scorer;
    let gen_is_etd = matches!(
        gen_spec.activation_method,
        Some(model::activation::ActivationMethod::ETD)
    );

    let ox_ev = oxonium_gate(gen_peaks, OXONIUM_GATE_MIN_FRAC, tol_ppm);

    let charges_to_try: Vec<u8> = glyco_charges_to_try(spec.precursor_charge, &params.charge_range);
    let max_frag_charge: u8 = spec
        .precursor_charge
        .filter(|&z| z > 0)
        .map(|z| z as u8)
        .unwrap_or(3);

    let all_backbone = retrieve_backbones(
        ctx,
        spec_idx,
        spec,
        gen_peaks,
        &charges_to_try,
        ox_ev.fired,
        gen_is_etd,
    );

    if all_backbone.is_empty() {
        return None;
    }

    let deduped_backbone = dedup_backbone_hits(all_backbone, tol_ppm);
    let hex_evidence = hex_diagnostic_evidence(gen_peaks, tol_ppm);
    let deduped_backbone: Vec<BackboneHit> = deduped_backbone
        .into_iter()
        .filter(|h| hex_gate_passes(ctx.core, h.glycan.as_ref(), hex_evidence))
        .collect();

    let mut scored_per_charge: Vec<(u8, ScoredSpectrum<'_>)> = Vec::new();
    for &z in &charges_to_try {
        if scored_per_charge.iter().all(|(c, _)| *c != z) {
            scored_per_charge.push((z, ScoredSpectrum::new(spec, spec_scorer, z)));
        }
    }
    let gen_scored_per_charge: Vec<(u8, ScoredSpectrum<'_>)> = if paired_hcd.is_some() {
        let mut v: Vec<(u8, ScoredSpectrum<'_>)> = Vec::new();
        for &z in &charges_to_try {
            if v.iter().all(|(c, _)| *c != z) {
                v.push((z, ScoredSpectrum::new(gen_spec, gen_scorer, z)));
            }
        }
        v
    } else {
        Vec::new()
    };
    let phase1_scored = if paired_hcd.is_some() {
        &gen_scored_per_charge
    } else {
        &scored_per_charge
    };
    let phase1_scorer: &scoring_crate::psm::RankScorer = if paired_hcd.is_some() {
        gen_scorer
    } else {
        spec_scorer
    };

    let core_y_counts: Vec<u8> = deduped_backbone
        .iter()
        .map(|h| {
            count_core_y_hits(
                gen_peaks,
                gen_stats,
                h.backbone_mass + H2O,
                tol_ppm,
                max_frag_charge,
                ctx.core,
            )
        })
        .collect();

    let scan = GlycoScan {
        ctx,
        spec_idx,
        deduped_backbone: &deduped_backbone,
        core_y_counts: &core_y_counts,
        scored_per_charge: &scored_per_charge,
        y_peaks,
        y_stats,
        gen_peaks,
        spec_scorer,
        ox_ev: &ox_ev,
        is_etd,
        gen_is_etd,
        pair_rank_etd: paired_hcd.is_some(),
        max_frag_charge,
        ladder_memo: std::cell::RefCell::new(vec![None; deduped_backbone.len()]),
        hyper_memo: std::cell::RefCell::new(FxHashMap::default()),
    };
    let (cheap_winners, backbone_best_rank, backbone_best_cz) =
        scan.cheap_winners(phase1_scored, phase1_scorer);
    let accepted_backbones = scan.accepted_backbones(&backbone_best_rank, &backbone_best_cz);

    let max_features = (effective_top_k * 2).max(accepted_backbones.len() * 4);
    let accepted_winners: Vec<(GlycanWinnerKey, CheapWinner)> = cheap_winners
        .into_iter()
        .filter(|(_, w)| accepted_backbones.contains(&w.bb_hit_idx))
        .collect();

    let spread = winner_spread(&accepted_winners);
    let winners_for_features = scan.winners_for_features(accepted_winners, max_features);

    let mut best_hits: FxHashMap<GlycanWinnerKey, FullGlycoPsm> =
        FxHashMap::with_capacity_and_hasher(winners_for_features.len(), Default::default());
    let fused_pos: FxHashMap<GlycanWinnerKey, usize> = if features_collapse {
        FxHashMap::default()
    } else {
        winners_for_features
            .iter()
            .enumerate()
            .map(|(i, (k, _))| (*k, i))
            .collect()
    };

    for (gl_key, w) in winners_for_features {
        best_hits.insert(gl_key, scan.full_hit(w, &spread));
    }

    if best_hits.is_empty() {
        None
    } else {
        let mut hits: Vec<(GlycanWinnerKey, FullGlycoPsm)> = best_hits.into_iter().collect();
        if features_collapse {
            hits.sort_by(|a, b| {
                b.1.psm
                    .rank_score
                    .total_cmp(&a.1.psm.rank_score)
                    .then_with(|| a.0.cmp(&b.0))
            });
        } else {
            hits.sort_by_key(|(k, _)| fused_pos.get(k).copied().unwrap_or(usize::MAX));
        }
        let retained_strong: Vec<f32> = hits
            .iter()
            .map(|(_, h)| h.psm.features.strong_score)
            .collect();
        for (_, h) in hits.iter_mut() {
            h.psm.features.strong_score_cal =
                strong_score_calibrated_loo(&retained_strong, h.psm.features.strong_score);
        }
        Some(GlycoSpectrumResult {
            spectrum_idx: spec_idx,
            hits: hits.into_iter().map(|(_, h)| h).collect(),
        })
    }
}

#[allow(clippy::needless_range_loop)]
fn build_hcd_partners(spectra: &[Spectrum]) -> Vec<Option<usize>> {
    use model::activation::ActivationMethod;
    const WINDOW: usize = 12;
    const MZ_TOL: f64 = 0.02;
    let n = spectra.len();
    let mut partner = vec![None; n];
    for i in 0..n {
        if spectra[i].activation_method != Some(ActivationMethod::ETD) {
            continue;
        }
        let mz = spectra[i].precursor_mz;
        if mz <= 0.0 {
            continue;
        }
        let lo = i.saturating_sub(WINDOW);
        let hi = (i + WINDOW + 1).min(n);
        let etd_z = spectra[i].precursor_charge;
        let mut best: Option<(usize, usize)> = None;
        for j in lo..hi {
            if j == i || spectra[j].activation_method != Some(ActivationMethod::HCD) {
                continue;
            }
            let charge_ok = match (etd_z, spectra[j].precursor_charge) {
                (Some(a), Some(b)) => a == b,
                _ => true,
            };
            if charge_ok && (mz - spectra[j].precursor_mz).abs() <= MZ_TOL {
                let d = i.abs_diff(j);
                if best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, j));
                }
            }
        }
        partner[i] = best.map(|(_, j)| j);
    }
    partner
}

#[allow(clippy::too_many_arguments)]
pub fn glyco_search_run(
    spectra: &[Spectrum],
    prepared: &PreparedSearch<'_>,
    table: &GlycoCandidates,
    glycan_list: &[GlycanComp],
    tol_ppm: f64,
    cfg: GlycoConfig,
) -> Vec<GlycoSpectrumResult> {
    let candidates = &table.candidates;
    let fragment_tolerance_da = prepared.fragment_tolerance_da;
    let hcd_partner: Vec<Option<usize>> = if cfg.hcd_pair {
        build_hcd_partners(spectra)
    } else {
        Vec::new()
    };

    let owned = GlycoCtxOwned::build(
        candidates,
        glycan_list,
        fragment_tolerance_da,
        cfg.clone(),
        tol_ppm,
    );
    let ctx = owned.as_ctx(prepared, table, glycan_list, tol_ppm, spectra, &hcd_partner);

    let out: Vec<GlycoSpectrumResult> = spectra
        .par_iter()
        .enumerate()
        .filter_map(|(spec_idx, spec)| score_spectrum_glyco(spec_idx, spec, &ctx))
        .collect();
    out
}

#[cfg(test)]
#[path = "glycopeptides/tests/glycan_source.rs"]
mod glycan_source_tests;

#[cfg(test)]
#[path = "glycopeptides/tests/glycopeptides.rs"]
mod tests;
