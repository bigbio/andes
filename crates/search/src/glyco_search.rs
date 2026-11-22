use rustc_hash::{FxHashMap, FxHashSet};

use model::mass::{H2O, ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use rayon::prelude::*;

use glyco::backbone::{
    core_y_intensity, count_core_y_hits, glycan_y_hit_frac, glycan_y_intensity,
    partial_glycan_by_intensity, y0y1_anchor_intensity, SpectrumStats,
};
use glyco::glycan_db::{GlycanComp, DENOVO_TOPOLOGY_ID};
use glyco::glycan_first::{search_glycans, Glycan, GlycanConfig, GlycanCore, GlycanIonIndex};
use glyco::glyco_psm::{
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

use crate::glyco_fragment_index::FragmentIndex;

use crate::match_engine::{compute_psm_features, PreparedSearch};
use crate::psm::PsmFeatures;
use crate::psm::PsmMatch;

type GlycanWinnerKey = (u32, u64);

use scoring_crate::scoring::{
    candidate_rank_entropy, cz_hyperscore_psm, cz_matched_intensity_frac, cz_structure_features,
    fuse_strong_score, hyperscore_psm_with_matches, listwise_score_gap, psm_edge_score, score_psm,
    score_psm_float, scored_for_charge, strong_score_calibrated_loo, ScoredSpectrum,
    StrongScoreInputs,
};

#[allow(clippy::too_many_arguments)]
fn glyco_candidate_strong_score(
    ss: &ScoredSpectrum<'_>,
    peptide: &model::peptide::Peptide,
    spec_scorer: &scoring_crate::scoring::RankScorer,
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
    candidates: Vec<crate::candidate_gen::Candidate>,
    mass_order: Vec<u32>,
}

impl GlycoCandidates {
    pub fn from_out_of_core(prepared: &PreparedSearch<'_>, core: GlycanCore) -> Self {
        Self::from_candidates(
            prepared.out_of_core_candidates_where(move |residues, post| match core {
                GlycanCore::NGlycan => crate::candidate_gen::span_has_nxst_sequon(residues, post),
                GlycanCore::OGlycan => crate::candidate_gen::span_has_ser_or_thr(residues, post),
            }),
            core,
        )
    }

    pub fn from_candidates(
        candidates: Vec<crate::candidate_gen::Candidate>,
        core: GlycanCore,
    ) -> Self {
        debug_assert!(candidates.iter().all(match core {
            GlycanCore::NGlycan => crate::candidate_gen::candidate_has_nxst_sequon,
            GlycanCore::OGlycan => crate::candidate_gen::candidate_has_ser_or_thr,
        }));
        let mass_order = crate::match_engine::build_mass_order(&candidates);
        Self {
            candidates,
            mass_order,
        }
    }

    pub fn candidates(&self) -> &[crate::candidate_gen::Candidate] {
        &self.candidates
    }
}

pub struct GlycoScoreCtx<'a> {
    pub params: &'a crate::search_params::SearchParams,
    pub scorer: &'a scoring_crate::scoring::RankScorer,
    pub core: GlycanCore,
    pub candidates: &'a [crate::candidate_gen::Candidate],
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
        candidates: &[crate::candidate_gen::Candidate],
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
                GlycanCore::NGlycan => glyco::glycan_first::default_diagnostic_ions(),
                GlycanCore::OGlycan => glyco::glycan_first::o_glycan_diagnostic_ions(),
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

fn precursor_window_da(params: &crate::search_params::SearchParams, mass: f64) -> f64 {
    let t = &params.precursor_tolerance;
    t.left.as_da(mass).max(t.right.as_da(mass))
}

fn backbone_candidates(
    mass_order: &[u32],
    candidates: &[crate::candidate_gen::Candidate],
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
    let candidates = ctx.candidates;
    let mass_order = ctx.mass_order;
    let fragment_tolerance_da = ctx.fragment_tolerance_da;
    let frag_index = ctx.frag_index;
    let glycan_sorted = ctx.glycan_sorted;
    let glycan_list = ctx.glycan_list;
    let tol_ppm = ctx.tol_ppm;
    let effective_top_k = GLYCO_BACKBONE_TOP_K;
    let max_peptide_first = GLYCO_MAX_PF;
    let peptide_first_on = ctx.peptide_first_on;
    let full_glycan_db = ctx.full_glycan_db;
    let gp_k = GLYCO_GP_K_DEFAULT;
    let gp_j = GLYCO_GP_J_DEFAULT;
    let gp_h = GLYCO_GP_H_DEFAULT;
    let gp_cz = GLYCO_GP_CZ_DEFAULT;
    let gp_iso = GLYCO_GP_ISO_DEFAULT;
    let min_core_y = ctx.min_core_y;
    let is_etd = matches!(
        spec.activation_method,
        Some(model::activation::ActivationMethod::ETD)
    );
    let cz_gate_on = is_etd;
    let features_collapse = ctx.features_collapse;
    let features_enumerated = ctx.features_enumerated;
    let scan_filter = ctx.scan_filter;
    const MIN_BY_MATCHES: u32 = 6;

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

    let spec_scorer: &scoring_crate::scoring::RankScorer = scorer;
    let gen_scorer: &scoring_crate::scoring::RankScorer = scorer;
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

    let isotope_window = isotope_window_for(
        spec_idx,
        &ctx.isotope_error_default,
        ctx.isotope_error_override.as_ref(),
        ctx.isotope_error_override_mask,
    );
    let iso_min = *isotope_window.start();
    let iso_max = *isotope_window.end();
    let mut all_backbone: Vec<BackboneHit> = Vec::new();
    for &z in &charges_to_try {
        let observed_neutral = (spec.precursor_mz - PROTON) * z as f64 - H2O;
        for iso in iso_min..=iso_max {
            let precursor_neutral = observed_neutral - iso as f64 * ISOTOPE;
            if precursor_neutral <= 0.0 {
                continue;
            }
            let source = glycan_source(full_glycan_db, ox_ev.fired, gen_is_etd);
            let restricted_glycans: Vec<GlycanComp> = if source != GlycanSource::IonIndex {
                Vec::new()
            } else {
                let res = search_glycans(
                    ctx.glycan_ion_index,
                    ctx.glycan_first_cfg,
                    glycan_list,
                    precursor_neutral + H2O,
                    z,
                    gen_peaks,
                );
                res.candidates
                    .iter()
                    .filter_map(|c| glycan_list.get(c.glycan_id as usize).cloned())
                    .collect()
            };
            let effective_glycans: &[GlycanComp] = match source {
                GlycanSource::FullList => glycan_list,
                GlycanSource::IonIndex | GlycanSource::None => &restricted_glycans,
            };
            for h in db_branch(precursor_neutral, effective_glycans, 500.0, z, iso) {
                all_backbone.push(h);
            }
        }
    }
    if peptide_first_on && ox_ev.fired && !full_glycan_db {
        let mut pf = frag_index.query(gen_peaks, MIN_BY_MATCHES);
        order_peptide_first(&mut pf);
        let mut pf_added = 0usize;
        'pf: for (cand_idx, _n) in pf {
            let pep_residue = candidates[cand_idx as usize].peptide.mass() - H2O;
            for &z in &charges_to_try {
                let observed_neutral = (spec.precursor_mz - PROTON) * z as f64 - H2O;
                for iso in iso_min..=iso_max {
                    let precursor_neutral = observed_neutral - iso as f64 * ISOTOPE;
                    let glycan_mass = precursor_neutral - pep_residue;
                    if glycan_mass < min_glycan_floor(ctx.core) {
                        continue;
                    }
                    let tol = precursor_window_da(params, precursor_neutral).max(0.02);
                    if let Some(g) =
                        nearest_glycan_mass(glycan_sorted, glycan_list, glycan_mass, tol)
                    {
                        let bb = precursor_neutral - g.mass;
                        all_backbone.push(BackboneHit {
                            backbone_mass: bb,
                            glycan: Some(g),
                            source: Source::Db,
                            charge: z,
                            isotope_offset: iso,
                            glycan_mass_residual: precursor_neutral - bb,
                            peptide_first: true,
                        });
                        pf_added += 1;
                        if pf_added >= max_peptide_first {
                            break 'pf;
                        }
                    }
                }
            }
        }
    }

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
    let phase1_scorer: &scoring_crate::scoring::RankScorer = if paired_hcd.is_some() {
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

    #[derive(Clone, Copy)]
    struct CheapWinner {
        bb_hit_idx: usize,
        cand_slot: usize,
        z: u8,
        isotope_offset: i8,
        rank: f32,
        score: f32,
        edge: i32,
        cand_residue_mass: f64,
    }
    let mut cheap_winners: FxHashMap<GlycanWinnerKey, CheapWinner> =
        FxHashMap::with_capacity_and_hasher(1024, Default::default());

    let mut backbone_best_rank: Vec<f32> = vec![f32::NEG_INFINITY; deduped_backbone.len()];

    let mut backbone_best_cz: Vec<f32> = vec![f32::NEG_INFINITY; deduped_backbone.len()];

    let scan_has_evidence = deduped_backbone
        .iter()
        .enumerate()
        .any(|(i, h)| has_backbone_evidence(h, core_y_counts[i]));

    for (bb_idx, bb_hit) in deduped_backbone.iter().enumerate() {
        if scan_has_evidence && !has_backbone_evidence(bb_hit, core_y_counts[bb_idx]) && !cz_gate_on
        {
            continue;
        }
        let bb_residue = bb_hit.backbone_mass;
        let z = bb_hit.charge;

        let tol_da = (bb_residue * tol_ppm * 1e-6_f64).max(0.01);

        let ss = match scored_for_charge(phase1_scored, z) {
            Some(s) => s,
            None => continue,
        };

        let candidate_slots = backbone_candidates(mass_order, candidates, bb_residue, tol_da);

        for cand_slot in candidate_slots {
            let cand = &candidates[cand_slot];
            let cand_residue_mass = cand.peptide.mass() - H2O;

            if (cand_residue_mass - bb_residue).abs() > tol_da {
                continue;
            }

            let gl_key = match &bb_hit.glycan {
                Some(g) => (cand_slot as u32, g.topology_id),
                None => (cand_slot as u32, DENOVO_TOPOLOGY_ID),
            };

            let scoring_pep: std::borrow::Cow<'_, model::peptide::Peptide> = if gen_is_etd {
                let gmass = bb_hit
                    .glycan
                    .as_ref()
                    .map(|g| g.mass)
                    .unwrap_or(bb_hit.glycan_mass_residual);
                if gmass > 0.0 {
                    let gsite = glyco_site_for(&cand.peptide, ctx.core);
                    std::borrow::Cow::Owned(glyco_aware_peptide(&cand.peptide, gsite, gmass))
                } else {
                    std::borrow::Cow::Borrowed(&cand.peptide)
                }
            } else {
                std::borrow::Cow::Borrowed(&cand.peptide)
            };

            let sc = score_psm(ss, &scoring_pep, phase1_scorer, z, fragment_tolerance_da);
            let ei = psm_edge_score(ss, &scoring_pep, phase1_scorer, z);
            let rk = sc + ei as f32;

            if rk > backbone_best_rank[bb_idx] {
                backbone_best_rank[bb_idx] = rk;
            }

            if cz_gate_on {
                let gmass = bb_hit
                    .glycan
                    .as_ref()
                    .map(|g| g.mass)
                    .unwrap_or(bb_hit.glycan_mass_residual);
                if gmass > 0.0 {
                    if let Some(etd_ss) = scored_for_charge(&scored_per_charge, z) {
                        let czs = cz_score_at_site(
                            etd_ss,
                            &cand.peptide,
                            gmass,
                            max_frag_charge,
                            tol_ppm,
                            ctx.core,
                        );
                        if czs > backbone_best_cz[bb_idx] {
                            backbone_best_cz[bb_idx] = czs;
                        }
                    }
                }
            }

            let w = CheapWinner {
                bb_hit_idx: bb_idx,
                cand_slot,
                z,
                isotope_offset: bb_hit.isotope_offset,
                rank: rk,
                score: sc,
                edge: ei,
                cand_residue_mass,
            };
            cheap_winners
                .entry(gl_key)
                .and_modify(|existing| {
                    if w.rank > existing.rank {
                        *existing = w;
                    }
                })
                .or_insert(w);
        }
    }

    let mut by_by: Vec<usize> = (0..deduped_backbone.len()).collect();
    by_by.sort_by(|&ai, &bi| {
        backbone_best_rank[bi]
            .total_cmp(&backbone_best_rank[ai])
            .then_with(|| core_y_counts[bi].cmp(&core_y_counts[ai]))
            .then_with(|| {
                deduped_backbone[bi]
                    .backbone_mass
                    .total_cmp(&deduped_backbone[ai].backbone_mass)
            })
            .then_with(|| bi.cmp(&ai))
    });
    by_by.truncate(effective_top_k);
    let mut accepted_backbones: FxHashSet<usize> = by_by.into_iter().collect();

    if cz_gate_on {
        let mut by_cz: Vec<usize> = (0..deduped_backbone.len())
            .filter(|&i| backbone_best_cz[i] > 0.0)
            .collect();
        by_cz.sort_by(|&ai, &bi| {
            backbone_best_cz[bi]
                .total_cmp(&backbone_best_cz[ai])
                .then_with(|| backbone_best_rank[bi].total_cmp(&backbone_best_rank[ai]))
                .then_with(|| {
                    deduped_backbone[bi]
                        .backbone_mass
                        .total_cmp(&deduped_backbone[ai].backbone_mass)
                })
                .then_with(|| bi.cmp(&ai))
        });
        by_cz.truncate(effective_top_k);
        accepted_backbones.extend(by_cz);
    }

    let max_features = (effective_top_k * 2).max(accepted_backbones.len() * 4);
    let accepted_winners: Vec<(GlycanWinnerKey, CheapWinner)> = cheap_winners
        .into_iter()
        .filter(|(_, w)| accepted_backbones.contains(&w.bb_hit_idx))
        .collect();

    let (tailor_denom, spectrum_listwise_gap, spectrum_rank_entropy, spectrum_delta_raw) = {
        let mut hist: FxHashMap<i32, u32> = FxHashMap::default();
        let mut scores: Vec<f32> = Vec::with_capacity(accepted_winners.len());
        for (_, w) in &accepted_winners {
            *hist.entry(w.score.round() as i32).or_insert(0) += 1;
            scores.push(w.score);
        }
        let denom = crate::psm::tailor_denominator(&hist, accepted_winners.len() as u32) as f32;
        scores.sort_by(|a, b| b.total_cmp(a));
        let delta_raw = match scores.first() {
            Some(&best) => scores
                .iter()
                .find(|&&v| v < best)
                .map(|&second| best - second)
                .unwrap_or(0.0),
            None => 0.0,
        };
        (
            denom,
            listwise_score_gap(&scores),
            candidate_rank_entropy(&scores),
            delta_raw,
        )
    };

    let ladder_raw = |bb_hit_idx: usize| -> f32 {
        let bb = &deduped_backbone[bb_hit_idx];
        let bbn = bb.backbone_mass + H2O;
        match &bb.glycan {
            Some(g) => {
                glycan_y_intensity(y_peaks, y_stats, bbn, g, tol_ppm, max_frag_charge) as f32
            }
            None => {
                core_y_intensity(y_peaks, y_stats, bbn, tol_ppm, max_frag_charge, ctx.core) as f32
            }
        }
    };
    let ladder_memo: std::cell::RefCell<Vec<Option<f32>>> =
        std::cell::RefCell::new(vec![None; deduped_backbone.len()]);
    let ladder = |w: &CheapWinner| -> f32 {
        if let Some(v) = ladder_memo.borrow()[w.bb_hit_idx] {
            return v;
        }
        let v = ladder_raw(w.bb_hit_idx);
        ladder_memo.borrow_mut()[w.bb_hit_idx] = Some(v);
        v
    };

    let hyper_m_raw = |w: &CheapWinner| -> (f32, u32) {
        match scored_for_charge(&scored_per_charge, w.z) {
            Some(ss) => {
                let pep = &candidates[w.cand_slot].peptide;
                if is_etd {
                    let bb = &deduped_backbone[w.bb_hit_idx];
                    let gmass = bb
                        .glycan
                        .as_ref()
                        .map(|g| g.mass)
                        .unwrap_or(bb.glycan_mass_residual);
                    if gmass > 0.0 {
                        let gsite = glyco_site_for(pep, ctx.core);
                        let modified = glyco_aware_peptide(pep, gsite, gmass);
                        hyperscore_psm_with_matches(ss, &modified, spec_scorer)
                    } else {
                        hyperscore_psm_with_matches(ss, pep, spec_scorer)
                    }
                } else {
                    hyperscore_psm_with_matches(ss, pep, spec_scorer)
                }
            }
            None => (0.0, 0),
        }
    };
    type HyperMemo = std::cell::RefCell<FxHashMap<(usize, usize, u8), (f32, u32)>>;
    let hyper_memo: HyperMemo = std::cell::RefCell::new(FxHashMap::default());
    let hyper_m = |w: &CheapWinner| -> (f32, u32) {
        let key = (w.bb_hit_idx, w.cand_slot, w.z);
        if let Some(v) = hyper_memo.borrow().get(&key) {
            return *v;
        }
        let v = hyper_m_raw(w);
        hyper_memo.borrow_mut().insert(key, v);
        v
    };
    let hyper = |w: &CheapWinner| -> f32 { hyper_m(w).0 };

    let cz = |w: &CheapWinner| -> f32 {
        if !is_etd {
            return 0.0;
        }
        let ss = match scored_for_charge(&scored_per_charge, w.z) {
            Some(ss) => ss,
            None => return 0.0,
        };
        let bb = &deduped_backbone[w.bb_hit_idx];
        let gmass = bb
            .glycan
            .as_ref()
            .map(|g| g.mass)
            .unwrap_or(bb.glycan_mass_residual);
        let pep = &candidates[w.cand_slot].peptide;
        cz_score_at_site(ss, pep, gmass, max_frag_charge, tol_ppm, ctx.core)
    };

    let cz_int = |w: &CheapWinner| -> f32 {
        if !is_etd {
            return 0.0;
        }
        let ss = match scored_for_charge(&scored_per_charge, w.z) {
            Some(ss) => ss,
            None => return 0.0,
        };
        let bb = &deduped_backbone[w.bb_hit_idx];
        let gmass = bb
            .glycan
            .as_ref()
            .map(|g| g.mass)
            .unwrap_or(bb.glycan_mass_residual);
        let pep = &candidates[w.cand_slot].peptide;
        let gsite = glyco_site_for(pep, ctx.core);
        cz_matched_intensity_frac(ss, pep, gmass, gsite, max_frag_charge, tol_ppm)
    };

    let cz_struct_on = true;
    let cz_struct = |w: &CheapWinner| -> (f32, f32) {
        if !is_etd || !cz_struct_on {
            return (0.0, 0.0);
        }
        let ss = match scored_for_charge(&scored_per_charge, w.z) {
            Some(ss) => ss,
            None => return (0.0, 0.0),
        };
        let bb = &deduped_backbone[w.bb_hit_idx];
        let gmass = bb
            .glycan
            .as_ref()
            .map(|g| g.mass)
            .unwrap_or(bb.glycan_mass_residual);
        let pep = &candidates[w.cand_slot].peptide;
        let gsite = glyco_site_for(pep, ctx.core);
        cz_structure_features(ss, pep, gmass, gsite, max_frag_charge, tol_ppm)
    };

    let pair_rank_etd = paired_hcd.is_some();
    let rank_sel = |w: &CheapWinner| -> f32 {
        if !pair_rank_etd {
            return w.rank;
        }
        match scored_for_charge(&scored_per_charge, w.z) {
            Some(ss) => {
                let pep = &candidates[w.cand_slot].peptide;
                let decorate = false;
                let scoring_pep: std::borrow::Cow<'_, model::peptide::Peptide> = if decorate {
                    let bb = &deduped_backbone[w.bb_hit_idx];
                    let gmass = bb
                        .glycan
                        .as_ref()
                        .map(|g| g.mass)
                        .unwrap_or(bb.glycan_mass_residual);
                    if gmass > 0.0 {
                        let gsite = glyco_site_for(pep, ctx.core);
                        std::borrow::Cow::Owned(glyco_aware_peptide(pep, gsite, gmass))
                    } else {
                        std::borrow::Cow::Borrowed(pep)
                    }
                } else {
                    std::borrow::Cow::Borrowed(pep)
                };
                let sc = score_psm(ss, &scoring_pep, spec_scorer, w.z, fragment_tolerance_da);
                let ei = psm_edge_score(ss, &scoring_pep, spec_scorer, w.z);
                sc + ei as f32
            }
            None => w.rank,
        }
    };

    let winners_for_features: Vec<(GlycanWinnerKey, CheapWinner)> = if features_collapse {
        let fused_of = |e: &(GlycanWinnerKey, CheapWinner)| {
            let cy = core_y_counts[e.1.bb_hit_idx] as f32;
            glyco_gp_fused_score_iso(
                rank_sel(&e.1),
                ladder(&e.1),
                cy,
                hyper(&e.1),
                gp_k,
                gp_j,
                gp_h,
                e.1.isotope_offset,
                gp_iso,
            ) + gp_cz * cz(&e.1)
        };
        let best = accepted_winners
            .iter()
            .map(|e| (e, fused_of(e)))
            .max_by(|(ea, sa), (eb, sb)| {
                sa.total_cmp(sb).then_with(|| eb.0.cmp(&ea.0))
            })
            .map(|(e, _)| e);
        match best {
            Some((gl_key, w)) => {
                let is_enum = deduped_backbone[w.bb_hit_idx].glycan.is_some();
                if features_enumerated && !is_enum {
                    accepted_winners
                        .iter()
                        .filter(|e| deduped_backbone[e.1.bb_hit_idx].glycan.is_some())
                        .map(|e| (e, fused_of(e)))
                        .max_by(|(ea, sa), (eb, sb)| sa.total_cmp(sb).then_with(|| eb.0.cmp(&ea.0)))
                        .map(|(e, _)| vec![(e.0, e.1)])
                        .unwrap_or_default()
                } else {
                    vec![(*gl_key, *w)]
                }
            }
            None => Vec::new(),
        }
        .into_iter()
        .filter(|(_, w)| {
            core_y_counts[w.bb_hit_idx] as u32 >= min_core_y
        })
        .collect::<Vec<_>>()
    } else {
        let mut scored: Vec<(f32, (GlycanWinnerKey, CheapWinner))> = accepted_winners
            .into_iter()
            .map(|e| {
                let cy = core_y_counts[e.1.bb_hit_idx] as f32;
                let s = glyco_gp_fused_score_iso(
                    rank_sel(&e.1),
                    ladder(&e.1),
                    cy,
                    hyper(&e.1),
                    gp_k,
                    gp_j,
                    gp_h,
                    e.1.isotope_offset,
                    gp_iso,
                ) + gp_cz * cz(&e.1);
                (s, e)
            })
            .collect();
        scored.sort_by(|a, b| b.0.total_cmp(&a.0).then_with(|| (a.1).0.cmp(&(b.1).0)));
        scored.truncate(max_features);
        scored.into_iter().map(|(_, e)| e).collect()
    };

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
        let bb_hit = &deduped_backbone[w.bb_hit_idx];
        let bb_residue = bb_hit.backbone_mass;
        let bb_neutral = bb_residue + H2O;

        let ss = scored_for_charge(&scored_per_charge, w.z)
            .expect("ScoredSpectrum must exist for winning charge");
        let cand = &candidates[w.cand_slot];
        let feat_pep: std::borrow::Cow<'_, model::peptide::Peptide> =
            std::borrow::Cow::Borrowed(&cand.peptide);
        let (mut features, strong_score) = glyco_candidate_strong_score(
            ss,
            &feat_pep,
            spec_scorer,
            w.z,
            spectrum_rank_entropy,
            spectrum_listwise_gap,
        );
        let rank_float_pep: std::borrow::Cow<'_, model::peptide::Peptide> = if is_etd {
            let gmass = bb_hit
                .glycan
                .as_ref()
                .map(|g| g.mass)
                .unwrap_or(bb_hit.glycan_mass_residual);
            if gmass > 0.0 {
                let gsite = glyco_site_for(&cand.peptide, ctx.core);
                std::borrow::Cow::Owned(glyco_aware_peptide(&cand.peptide, gsite, gmass))
            } else {
                std::borrow::Cow::Borrowed(&cand.peptide)
            }
        } else {
            std::borrow::Cow::Borrowed(&cand.peptide)
        };
        features.rank_score_float =
            score_psm_float(ss, &rank_float_pep, spec_scorer, w.z, fragment_tolerance_da);
        features.tailor_score = if tailor_denom > 0.0 {
            w.score / tailor_denom
        } else {
            w.score
        };
        features.candidate_rank_entropy = spectrum_rank_entropy;
        features.delta_raw_score = spectrum_delta_raw;
        features.listwise_score_gap = spectrum_listwise_gap;
        features.strong_score = strong_score;

        let mass_error_ppm = if bb_residue > 0.0 {
            (w.cand_residue_mass - bb_residue) / bb_residue * 1e6
        } else {
            0.0
        };
        let psm = PsmMatch {
            spectrum_idx: spec_idx,
            candidate_idxs: smallvec::smallvec![w.cand_slot as u32],
            charge_used: w.z,
            mass_error_ppm,
            score: w.score,
            rank_score: if pair_rank_etd { rank_sel(&w) } else { w.rank },
            edge_score: w.edge,
            activation_method: Some(spec_scorer.param().data_type.activation),
            features,
            isotope_offset: w.isotope_offset,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: crate::psm::tie_key_of(cand),
        };
        let glycan_mass = bb_hit
            .glycan
            .as_ref()
            .map(|g| g.mass)
            .unwrap_or(bb_hit.glycan_mass_residual);
        let cz_struct_vals = cz_struct(&w);
        let glycan_key = GlycoPsmKey {
            spectrum_idx: spec_idx,
            glycan: bb_hit.glycan.clone(),
            glycan_source: bb_hit.source.clone(),
            oxonium_summed_frac: ox_ev.summed_frac,
            n_core_oxonium_ions: ox_ev.n_core_ions,
            y_ladder_intensity_score: match &bb_hit.glycan {
                Some(g) => {
                    glycan_y_intensity(y_peaks, y_stats, bb_neutral, g, tol_ppm, max_frag_charge)
                        as f32
                }
                None => core_y_intensity(
                    y_peaks,
                    y_stats,
                    bb_neutral,
                    tol_ppm,
                    max_frag_charge,
                    ctx.core,
                ) as f32,
            },
            y_hit_frac: match &bb_hit.glycan {
                Some(g) => glycan_y_hit_frac(
                    y_peaks,
                    y_stats,
                    bb_neutral,
                    g,
                    tol_ppm,
                    max_frag_charge,
                    None,
                ) as f32,
                None => 0.0,
            },
            partial_glycan_by: {
                let residues: Vec<f64> = cand
                    .peptide
                    .residues
                    .iter()
                    .map(|aa| aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
                    .collect();
                partial_glycan_by_intensity(y_peaks, y_stats, &residues, tol_ppm, max_frag_charge)
                    as f32
            },
            y0y1_anchor_score: y0y1_anchor_intensity(
                y_peaks,
                y_stats,
                cand.peptide.mass(),
                w.z,
                tol_ppm,
            ) as f32,
            sialic_consistency: match &bb_hit.glycan {
                Some(g) => sialic_consistency(gen_peaks, g, tol_ppm),
                None => 0.0,
            },
            core_y_hits: core_y_counts[w.bb_hit_idx],
            glycan_mass,
            backbone_mass: bb_neutral,
            cz_hyperscore: cz(&w),
            cz_intensity: cz_int(&w),
            cz_explained: cz_struct_vals.0,
            cz_chance_llr: cz_struct_vals.1,
        };
        best_hits.insert(gl_key, FullGlycoPsm { glycan_key, psm });
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
mod glycan_source_tests {
    use super::{glycan_source, GlycanSource};

    #[test]
    fn hcd_with_oxonium_uses_the_ion_index() {
        assert_eq!(glycan_source(false, true, false), GlycanSource::IonIndex);
    }

    #[test]
    fn hcd_without_oxonium_gets_nothing() {
        assert_eq!(glycan_source(false, false, false), GlycanSource::None);
    }

    #[test]
    fn etd_without_oxonium_falls_back_to_the_full_list() {
        assert_eq!(glycan_source(false, false, true), GlycanSource::FullList);
    }

    #[test]
    fn etd_with_oxonium_behaves_like_any_other_scan() {
        assert_eq!(glycan_source(false, true, true), GlycanSource::IonIndex);
        assert_eq!(glycan_source(true, true, true), GlycanSource::FullList);
    }

    #[test]
    fn full_glycan_db_uses_the_whole_list_when_oxonium_fired() {
        assert_eq!(glycan_source(true, true, false), GlycanSource::FullList);
    }

    #[test]
    fn full_glycan_db_without_oxonium_on_hcd_still_gets_nothing() {
        assert_eq!(glycan_source(true, false, false), GlycanSource::None);
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn glycan_subtraction_uses_the_precursor_tolerance() {
        use model::tolerance::{PrecursorTolerance, Tolerance};
        let aa = model::AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = crate::search_params::SearchParams::default_tryptic(aa);
        params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(10.0));
        assert!((super::precursor_window_da(&params, 3000.0) - 0.03).abs() < 1e-12);
        params.precursor_tolerance =
            PrecursorTolerance::asymmetric(Tolerance::Ppm(5.0), Tolerance::Da(0.05));
        assert!((super::precursor_window_da(&params, 3000.0) - 0.05).abs() < 1e-12);
    }

    #[test]
    fn o_glycan_floor_admits_a_light_tn_precursor() {
        use super::{min_glycan_floor, GlycanCore};
        use glyco::glycan_mass::HEXNAC;
        assert!(HEXNAC - 0.01 >= min_glycan_floor(GlycanCore::OGlycan));
        assert!(HEXNAC < min_glycan_floor(GlycanCore::NGlycan));
    }

    #[test]
    fn o_glycan_hex_gate_drops_hex_without_evidence() {
        use super::{hex_gate_passes, GlycanComp, GlycanCore};
        let hex_glycan = GlycanComp::from_counts(1, 1, 0, 0, 0, 0);
        let no_hex = GlycanComp::from_counts(1, 0, 0, 0, 0, 0);
        assert!(hex_gate_passes(
            GlycanCore::NGlycan,
            Some(&hex_glycan),
            false
        ));
        assert!(hex_gate_passes(GlycanCore::OGlycan, Some(&no_hex), false));
        assert!(hex_gate_passes(GlycanCore::OGlycan, None, false));
        assert!(!hex_gate_passes(
            GlycanCore::OGlycan,
            Some(&hex_glycan),
            false
        ));
        assert!(hex_gate_passes(
            GlycanCore::OGlycan,
            Some(&hex_glycan),
            true
        ));
    }

    #[test]
    fn peptide_first_backbones_are_exempt_from_the_core_y_prefilter() {
        use super::{dedup_backbone_hits, has_backbone_evidence, BackboneHit, Source};
        let mass_only = BackboneHit {
            backbone_mass: 1500.0,
            glycan: None,
            source: Source::Db,
            charge: 3,
            isotope_offset: 0,
            glycan_mass_residual: 1800.0,
            peptide_first: false,
        };
        let peptide_first = BackboneHit {
            peptide_first: true,
            ..mass_only.clone()
        };
        assert!(!has_backbone_evidence(&mass_only, 0));
        assert!(has_backbone_evidence(&mass_only, 1));
        assert!(has_backbone_evidence(&peptide_first, 0));
        let merged = dedup_backbone_hits(vec![mass_only.clone(), peptide_first.clone()], 20.0);
        assert_eq!(merged.len(), 1);
        assert!(merged[0].peptide_first);
        let merged = dedup_backbone_hits(vec![peptide_first, mass_only], 20.0);
        assert_eq!(merged.len(), 1);
        assert!(merged[0].peptide_first);
    }

    #[test]
    fn isotope_window_for_applies_the_override_per_spectrum() {
        use super::isotope_window_for;
        let default = 0..=2;
        let narrow = 0..=1;
        assert_eq!(isotope_window_for(0, &default, None, None), 0..=2);
        assert_eq!(isotope_window_for(0, &default, None, Some(&[true])), 0..=2);
        assert_eq!(isotope_window_for(7, &default, Some(&narrow), None), 0..=1);
        let mask = [true, false, true];
        assert_eq!(
            isotope_window_for(0, &default, Some(&narrow), Some(&mask)),
            0..=1
        );
        assert_eq!(
            isotope_window_for(1, &default, Some(&narrow), Some(&mask)),
            0..=2
        );
        assert_eq!(
            isotope_window_for(2, &default, Some(&narrow), Some(&mask)),
            0..=1
        );
        assert_eq!(
            isotope_window_for(3, &default, Some(&narrow), Some(&mask)),
            0..=2
        );
    }

    use super::*;

    #[test]
    fn build_hcd_partners_pairs_etd_to_precursor_matched_hcd() {
        use model::activation::ActivationMethod;
        let sp = |mz: f64, am: ActivationMethod| Spectrum {
            title: String::new(),
            precursor_mz: mz,
            precursor_intensity: None,
            precursor_charge: Some(3),
            rt_seconds: None,
            scan: None,
            peaks: vec![],
            activation_method: Some(am),
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        use ActivationMethod::{ETD, HCD};
        let spectra = vec![
            sp(500.0, HCD),
            sp(500.0, ETD),
            sp(600.0, HCD),
            sp(600.005, ETD),
            sp(999.0, ETD),
        ];
        let p = build_hcd_partners(&spectra);
        assert_eq!(p[0], None, "HCD spectra never get a partner");
        assert_eq!(p[1], Some(0), "ETD@500 pairs to HCD@500");
        assert_eq!(p[2], None);
        assert_eq!(p[3], Some(2), "ETD@600.005 pairs to HCD@600 within m/z tol");
        assert_eq!(p[4], None, "ETD@999 has no precursor-matched HCD partner");
    }

    #[test]
    fn order_peptide_first_is_deterministic_across_input_orders() {
        let permutations = vec![
            vec![(8u32, 9u32), (5, 7), (2, 5), (3, 9), (1, 7)],
            vec![(1, 7), (3, 9), (2, 5), (8, 9), (5, 7)],
            vec![(2, 5), (1, 7), (5, 7), (8, 9), (3, 9)],
            vec![(5, 7), (8, 9), (3, 9), (1, 7), (2, 5)],
        ];

        let mut reference = permutations[0].clone();
        order_peptide_first(&mut reference);

        for perm in &permutations {
            let mut got = perm.clone();
            order_peptide_first(&mut got);
            assert_eq!(
                &got, &reference,
                "input order {perm:?} must collapse to the same deterministic order"
            );
        }

        for w in reference.windows(2) {
            assert!(
                w[0].1 >= w[1].1,
                "counts must be non-increasing: {reference:?}"
            );
        }
        assert_eq!(reference[0].1, 9, "strongest count first");
        assert_eq!(reference.last().unwrap().1, 5, "weakest count last");
    }

    #[test]
    fn isotope_sweep_annotates_glycan_only_at_correct_offset() {
        use glyco::glycan_db::GlycanComp;
        use glyco::glycan_mass::{CORE_Y_STEPS, PROTON as GLY_PROTON};
        use glyco::hybrid::hybrid_candidates_with_isotope;

        let true_backbone_residue = 1500.0_f64;
        let glycan = GlycanComp::from_counts(2, 3, 0, 0, 0, 0);
        let true_precursor_neutral = true_backbone_residue + glycan.mass;

        let observed_neutral = true_precursor_neutral + ISOTOPE;

        let glycans = vec![glycan];

        let y0_neutral = true_backbone_residue + H2O;
        let mut peaks: Vec<(f64, f32)> = vec![(204.08665, 200.0), (138.05496, 120.0)];
        peaks.push((y0_neutral + GLY_PROTON, 150.0));
        for &s in CORE_Y_STEPS.iter() {
            peaks.push((y0_neutral + s + GLY_PROTON, 100.0));
        }

        let hits0 =
            hybrid_candidates_with_isotope(&peaks, observed_neutral, 2, 0, &glycans, 20.0, 5);
        let matching0: Vec<_> = hits0
            .iter()
            .filter(|h| (h.backbone_mass - true_backbone_residue).abs() < 0.05)
            .collect();
        assert!(
            !matching0.is_empty(),
            "backbone is read from the ladder → recovered even at the wrong isotope offset"
        );
        assert!(
            matching0.iter().any(|h| h.source == Source::DeNovo),
            "at the wrong offset the glycan is ~1 ISOTOPE off and must NOT annotate to a known composition"
        );
        assert!(
            matching0.iter().all(|h| h.source != Source::Db),
            "wrong isotope offset must not produce a DB annotation for the true backbone"
        );

        let precursor_neutral_iso1 = observed_neutral - ISOTOPE;
        let hits1 =
            hybrid_candidates_with_isotope(&peaks, precursor_neutral_iso1, 2, 1, &glycans, 20.0, 5);
        let hit = hits1
            .iter()
            .find(|h| {
                (h.backbone_mass - true_backbone_residue).abs() < 0.05 && h.source == Source::Db
            })
            .expect("offset +1 must recover AND annotate the backbone via the corrected precursor");
        assert_eq!(
            hit.isotope_offset, 1,
            "recovered hit must record isotope_offset=1"
        );
        assert_eq!(
            hit.charge, 2,
            "recovered hit must record the charge it was matched at"
        );
    }

    #[test]
    fn charges_to_try_is_the_reported_charge_or_the_range() {
        let range = 2u8..=3u8;
        assert_eq!(glyco_charges_to_try(Some(5), &range), vec![5]);
        assert_eq!(glyco_charges_to_try(Some(2), &range), vec![2]);
        assert_eq!(glyco_charges_to_try(None, &range), vec![2, 3]);
        assert_eq!(glyco_charges_to_try(Some(0), &range), vec![2, 3]);
        assert_eq!(glyco_charges_to_try(Some(-1), &range), vec![2, 3]);
    }

    #[test]
    fn dedup_preserves_distinct_isotope_residuals_for_novel_glycans() {
        use glyco::hybrid::BackboneHit;
        let mk = |iso: i8, residual: f64| BackboneHit {
            backbone_mass: 1500.0,
            glycan: None,
            source: Source::DeNovo,
            charge: 3,
            isotope_offset: iso,
            glycan_mass_residual: residual,
            peptide_first: false,
        };
        let out = dedup_backbone_hits(vec![mk(0, 892.317), mk(1, 891.313)], 20.0);
        assert_eq!(
            out.len(),
            2,
            "distinct isotope residuals must not be merged: {out:?}"
        );
        let residuals: Vec<f64> = out.iter().map(|h| h.glycan_mass_residual).collect();
        assert!(residuals.iter().any(|r| (r - 892.317).abs() < 1e-6));
        assert!(residuals.iter().any(|r| (r - 891.313).abs() < 1e-6));

        let dup = dedup_backbone_hits(vec![mk(0, 892.317), mk(0, 892.317)], 20.0);
        assert_eq!(dup.len(), 1, "same-hypothesis duplicates must merge");
    }

    #[test]
    fn full_glyco_psm_is_clone() {
        let psm = PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: smallvec::smallvec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: 0.0,
            rank_score: 0.0,
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        };
        let key = GlycoPsmKey {
            spectrum_idx: 0,
            glycan: None,
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.0,
            n_core_oxonium_ions: 0,
            y_ladder_intensity_score: 0.0,
            y_hit_frac: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.0,
            sialic_consistency: 0.0,
            core_y_hits: 0,
            glycan_mass: 0.0,
            backbone_mass: 0.0,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        };
        let hit = FullGlycoPsm {
            glycan_key: key,
            psm,
        };
        let cloned = hit.clone();
        assert_eq!(cloned.psm.spectrum_idx, 0);
    }

    #[test]
    fn glyco_spectrum_result_is_clone() {
        let result = GlycoSpectrumResult {
            spectrum_idx: 7,
            hits: vec![],
        };
        let c = result.clone();
        assert_eq!(c.spectrum_idx, 7);
    }

    #[test]
    fn core_y_ranking_promotes_supported_backbone_over_unsupported() {
        use glyco::backbone::{count_core_y_hits, SpectrumStats};
        use glyco::glycan_mass::{CORE_Y_STEPS, PROTON};

        let true_bb = 1100.0_f64;

        let noise_bb = 2800.0_f64;

        let mut peaks: Vec<(f64, f32)> = vec![
            (true_bb + PROTON, 500.0),
            (true_bb + PROTON + CORE_Y_STEPS[0], 400.0),
            (true_bb + PROTON + CORE_Y_STEPS[1], 350.0),
            (true_bb + PROTON + CORE_Y_STEPS[2], 300.0),
            (true_bb + PROTON + CORE_Y_STEPS[3], 250.0),
            (true_bb + PROTON + CORE_Y_STEPS[4], 200.0),
            (900.0, 10.0),
            (1050.0, 10.0),
        ];
        peaks.push((noise_bb + PROTON + 5.0, 50.0));

        let tol_ppm = 20.0;

        let stats = SpectrumStats::new(&peaks);
        let true_hits = count_core_y_hits(&peaks, &stats, true_bb, tol_ppm, 3, GlycanCore::NGlycan);
        let noise_hits =
            count_core_y_hits(&peaks, &stats, noise_bb, tol_ppm, 3, GlycanCore::NGlycan);

        assert_eq!(
            true_hits, 6,
            "expected all 6 core-Y hits for true_bb, got {}",
            true_hits
        );
        assert_eq!(
            noise_hits, 0,
            "expected 0 core-Y hits for noise_bb, got {}",
            noise_hits
        );

        let mut candidates = [
            (noise_bb, noise_hits),
            (true_bb, true_hits),
        ];
        candidates.sort_by(|&(am, ay), &(bm, by)| {
            by.cmp(&ay)
                .then_with(|| bm.partial_cmp(&am).unwrap_or(std::cmp::Ordering::Equal))
        });

        assert!(
            (candidates[0].0 - true_bb).abs() < 0.01,
            "expected true_bb ranked first after core-Y sort, got backbone_mass={}",
            candidates[0].0
        );
        assert!(
            (candidates[1].0 - noise_bb).abs() < 0.01,
            "expected noise_bb ranked second, got backbone_mass={}",
            candidates[1].0
        );
    }

    #[test]
    fn by_rank_promotes_byone_matching_backbone_over_y_ladder_backbone() {
        let true_bb_mass = 1500.0_f64;
        let true_bb_best_rank: f32 = 10.0;
        let true_bb_core_y: u8 = 0;

        let noise_bb_mass = 2000.0_f64;
        let noise_bb_best_rank: f32 = 2.0;
        let noise_bb_core_y: u8 = 6;

        let backbones = [
            (noise_bb_mass, noise_bb_best_rank, noise_bb_core_y),
            (true_bb_mass, true_bb_best_rank, true_bb_core_y),
        ];
        let mut order: Vec<usize> = (0..backbones.len()).collect();
        order.sort_by(|&ai, &bi| {
            backbones[bi]
                .1
                .partial_cmp(&backbones[ai].1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| backbones[bi].2.cmp(&backbones[ai].2))
                .then_with(|| {
                    backbones[bi]
                        .0
                        .partial_cmp(&backbones[ai].0)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });

        assert_eq!(
            order[0], 1,
            "expected true_bb (idx=1) ranked first by b/y rank, got idx={}",
            order[0]
        );
        assert_eq!(
            order[1], 0,
            "expected noise_bb (idx=0) ranked second, got idx={}",
            order[1]
        );
    }
}
