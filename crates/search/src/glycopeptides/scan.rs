//! Per-spectrum glyco scoring for `score_spectrum_glyco`: backbone retrieval, the cheap
//! candidate pass and backbone selection, and the full scoring of the retained winners.

use super::*;

/// Minimum b/y matches for a peptide-first fragment-index hit.
const MIN_BY_MATCHES: u32 = 6;

/// Backbone hypotheses for one spectrum: per charge and isotope offset, the glycan database
/// branch (full list or ion-index restricted), then peptide-first hits from the fragment index
/// when the oxonium gate fired.
pub(super) fn retrieve_backbones(
    ctx: &GlycoScoreCtx<'_>,
    spec_idx: usize,
    spec: &Spectrum,
    gen_peaks: &[(f64, f32)],
    charges_to_try: &[u8],
    ox_fired: bool,
    gen_is_etd: bool,
) -> Vec<BackboneHit> {
    let params = ctx.params;
    let candidates = ctx.candidates;
    let frag_index = ctx.frag_index;
    let glycan_sorted = ctx.glycan_sorted;
    let glycan_list = ctx.glycan_list;
    let peptide_first_on = ctx.peptide_first_on;
    let full_glycan_db = ctx.full_glycan_db;
    let max_peptide_first = GLYCO_MAX_PF;
    let isotope_window = isotope_window_for(
        spec_idx,
        &ctx.isotope_error_default,
        ctx.isotope_error_override.as_ref(),
        ctx.isotope_error_override_mask,
    );
    let iso_min = *isotope_window.start();
    let iso_max = *isotope_window.end();
    let mut all_backbone: Vec<BackboneHit> = Vec::new();
    for &z in charges_to_try {
        let observed_neutral = (spec.precursor_mz - PROTON) * z as f64 - H2O;
        for iso in iso_min..=iso_max {
            let precursor_neutral = observed_neutral - iso as f64 * ISOTOPE;
            if precursor_neutral <= 0.0 {
                continue;
            }
            let source = glycan_source(full_glycan_db, ox_fired, gen_is_etd);
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
    if peptide_first_on && ox_fired && !full_glycan_db {
        let mut pf = frag_index.query(gen_peaks, MIN_BY_MATCHES);
        order_peptide_first(&mut pf);
        let mut pf_added = 0usize;
        'pf: for (cand_idx, _n) in pf {
            let pep_residue = candidates[cand_idx as usize].peptide.mass() - H2O;
            for &z in charges_to_try {
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
    all_backbone
}

/// The best-ranked candidate for one (candidate, glycan topology) key in the cheap pass.
#[derive(Clone, Copy)]
pub(super) struct CheapWinner {
    pub(super) bb_hit_idx: usize,
    cand_slot: usize,
    z: u8,
    isotope_offset: i8,
    rank: f32,
    score: f32,
    edge: i32,
    cand_residue_mass: f64,
}

type CheapWinners = FxHashMap<GlycanWinnerKey, CheapWinner>;
type HyperMemo = std::cell::RefCell<FxHashMap<(usize, usize, u8), (f32, u32)>>;

/// Score spread over a spectrum's accepted winners, copied into every hit's features.
#[derive(Clone, Copy)]
pub(super) struct WinnerSpread {
    tailor_denom: f32,
    spectrum_listwise_gap: f32,
    spectrum_rank_entropy: f32,
    spectrum_delta_raw: f32,
}

pub(super) fn winner_spread(accepted_winners: &[(GlycanWinnerKey, CheapWinner)]) -> WinnerSpread {
    let mut hist: FxHashMap<i32, u32> = FxHashMap::default();
    let mut scores: Vec<f32> = Vec::with_capacity(accepted_winners.len());
    for (_, w) in accepted_winners {
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
    WinnerSpread {
        tailor_denom: denom,
        spectrum_listwise_gap: listwise_score_gap(&scores),
        spectrum_rank_entropy: candidate_rank_entropy(&scores),
        spectrum_delta_raw: delta_raw,
    }
}

/// One spectrum's glyco scoring state once its backbones are retrieved: the deduplicated
/// backbones and their core-Y counts, the per-charge scored spectra, the peaks the Y ladder and
/// the generation pass read, and the memoised ladder and hyperscore terms.
pub(super) struct GlycoScan<'a, 'c, 's> {
    pub(super) ctx: &'a GlycoScoreCtx<'c>,
    pub(super) spec_idx: usize,
    pub(super) deduped_backbone: &'a [BackboneHit],
    pub(super) core_y_counts: &'a [u8],
    pub(super) scored_per_charge: &'a [(u8, ScoredSpectrum<'s>)],
    pub(super) y_peaks: &'a [(f64, f32)],
    pub(super) y_stats: &'a SpectrumStats,
    pub(super) gen_peaks: &'a [(f64, f32)],
    pub(super) spec_scorer: &'a scoring_crate::psm::RankScorer,
    pub(super) ox_ev: &'a glyco::oxonium::OxoniumEvidence,
    pub(super) is_etd: bool,
    pub(super) gen_is_etd: bool,
    pub(super) pair_rank_etd: bool,
    pub(super) max_frag_charge: u8,
    pub(super) ladder_memo: std::cell::RefCell<Vec<Option<f32>>>,
    pub(super) hyper_memo: HyperMemo,
}

impl GlycoScan<'_, '_, '_> {
    fn ladder_raw(&self, bb_hit_idx: usize) -> f32 {
        let ctx = self.ctx;
        let tol_ppm = self.ctx.tol_ppm;
        let deduped_backbone = self.deduped_backbone;
        let y_peaks = self.y_peaks;
        let y_stats = self.y_stats;
        let max_frag_charge = self.max_frag_charge;
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
    }

    fn ladder(&self, w: &CheapWinner) -> f32 {
        let ladder_memo = &self.ladder_memo;
        if let Some(v) = ladder_memo.borrow()[w.bb_hit_idx] {
            return v;
        }
        let v = self.ladder_raw(w.bb_hit_idx);
        ladder_memo.borrow_mut()[w.bb_hit_idx] = Some(v);
        v
    }

    fn hyper_m_raw(&self, w: &CheapWinner) -> (f32, u32) {
        let ctx = self.ctx;
        let candidates = self.ctx.candidates;
        let deduped_backbone = self.deduped_backbone;
        let scored_per_charge = self.scored_per_charge;
        let spec_scorer = self.spec_scorer;
        let is_etd = self.is_etd;
        match scored_for_charge(scored_per_charge, w.z) {
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
    }

    fn hyper_m(&self, w: &CheapWinner) -> (f32, u32) {
        let hyper_memo = &self.hyper_memo;
        let key = (w.bb_hit_idx, w.cand_slot, w.z);
        if let Some(v) = hyper_memo.borrow().get(&key) {
            return *v;
        }
        let v = self.hyper_m_raw(w);
        hyper_memo.borrow_mut().insert(key, v);
        v
    }

    fn hyper(&self, w: &CheapWinner) -> f32 {
        self.hyper_m(w).0
    }

    fn cz(&self, w: &CheapWinner) -> f32 {
        let ctx = self.ctx;
        let candidates = self.ctx.candidates;
        let tol_ppm = self.ctx.tol_ppm;
        let deduped_backbone = self.deduped_backbone;
        let scored_per_charge = self.scored_per_charge;
        let is_etd = self.is_etd;
        let max_frag_charge = self.max_frag_charge;
        if !is_etd {
            return 0.0;
        }
        let ss = match scored_for_charge(scored_per_charge, w.z) {
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
    }

    fn cz_int(&self, w: &CheapWinner) -> f32 {
        let ctx = self.ctx;
        let candidates = self.ctx.candidates;
        let tol_ppm = self.ctx.tol_ppm;
        let deduped_backbone = self.deduped_backbone;
        let scored_per_charge = self.scored_per_charge;
        let is_etd = self.is_etd;
        let max_frag_charge = self.max_frag_charge;
        if !is_etd {
            return 0.0;
        }
        let ss = match scored_for_charge(scored_per_charge, w.z) {
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
    }

    fn cz_struct(&self, w: &CheapWinner) -> (f32, f32) {
        let ctx = self.ctx;
        let candidates = self.ctx.candidates;
        let tol_ppm = self.ctx.tol_ppm;
        let deduped_backbone = self.deduped_backbone;
        let scored_per_charge = self.scored_per_charge;
        let is_etd = self.is_etd;
        let max_frag_charge = self.max_frag_charge;
        let cz_struct_on = true;
        if !is_etd || !cz_struct_on {
            return (0.0, 0.0);
        }
        let ss = match scored_for_charge(scored_per_charge, w.z) {
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
    }

    fn rank_sel(&self, w: &CheapWinner) -> f32 {
        let ctx = self.ctx;
        let candidates = self.ctx.candidates;
        let fragment_tolerance_da = self.ctx.fragment_tolerance_da;
        let deduped_backbone = self.deduped_backbone;
        let scored_per_charge = self.scored_per_charge;
        let spec_scorer = self.spec_scorer;
        let pair_rank_etd = self.pair_rank_etd;
        if !pair_rank_etd {
            return w.rank;
        }
        match scored_for_charge(scored_per_charge, w.z) {
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
    }

    pub(super) fn cheap_winners(
        &self,
        phase1_scored: &[(u8, ScoredSpectrum<'_>)],
        phase1_scorer: &scoring_crate::psm::RankScorer,
    ) -> (CheapWinners, Vec<f32>, Vec<f32>) {
        let ctx = self.ctx;
        let candidates = self.ctx.candidates;
        let mass_order = self.ctx.mass_order;
        let tol_ppm = self.ctx.tol_ppm;
        let fragment_tolerance_da = self.ctx.fragment_tolerance_da;
        let deduped_backbone = self.deduped_backbone;
        let core_y_counts = self.core_y_counts;
        let scored_per_charge = self.scored_per_charge;
        let gen_is_etd = self.gen_is_etd;
        let max_frag_charge = self.max_frag_charge;
        let cz_gate_on = self.is_etd;
        let mut cheap_winners: FxHashMap<GlycanWinnerKey, CheapWinner> =
            FxHashMap::with_capacity_and_hasher(1024, Default::default());

        let mut backbone_best_rank: Vec<f32> = vec![f32::NEG_INFINITY; deduped_backbone.len()];

        let mut backbone_best_cz: Vec<f32> = vec![f32::NEG_INFINITY; deduped_backbone.len()];

        let scan_has_evidence = deduped_backbone
            .iter()
            .enumerate()
            .any(|(i, h)| has_backbone_evidence(h, core_y_counts[i]));

        for (bb_idx, bb_hit) in deduped_backbone.iter().enumerate() {
            if scan_has_evidence
                && !has_backbone_evidence(bb_hit, core_y_counts[bb_idx])
                && !cz_gate_on
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
                        if let Some(etd_ss) = scored_for_charge(scored_per_charge, z) {
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
        (cheap_winners, backbone_best_rank, backbone_best_cz)
    }

    pub(super) fn accepted_backbones(
        &self,
        backbone_best_rank: &[f32],
        backbone_best_cz: &[f32],
    ) -> FxHashSet<usize> {
        let effective_top_k = GLYCO_BACKBONE_TOP_K;
        let deduped_backbone = self.deduped_backbone;
        let core_y_counts = self.core_y_counts;
        let cz_gate_on = self.is_etd;
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
        accepted_backbones
    }

    pub(super) fn winners_for_features(
        &self,
        accepted_winners: Vec<(GlycanWinnerKey, CheapWinner)>,
        max_features: usize,
    ) -> Vec<(GlycanWinnerKey, CheapWinner)> {
        let features_collapse = self.ctx.features_collapse;
        let features_enumerated = self.ctx.features_enumerated;
        let min_core_y = self.ctx.min_core_y;
        let gp_k = GLYCO_GP_K_DEFAULT;
        let gp_j = GLYCO_GP_J_DEFAULT;
        let gp_h = GLYCO_GP_H_DEFAULT;
        let gp_cz = GLYCO_GP_CZ_DEFAULT;
        let gp_iso = GLYCO_GP_ISO_DEFAULT;
        let deduped_backbone = self.deduped_backbone;
        let core_y_counts = self.core_y_counts;
        if features_collapse {
            let fused_of = |e: &(GlycanWinnerKey, CheapWinner)| {
                let cy = core_y_counts[e.1.bb_hit_idx] as f32;
                glyco_gp_fused_score_iso(
                    self.rank_sel(&e.1),
                    self.ladder(&e.1),
                    cy,
                    self.hyper(&e.1),
                    gp_k,
                    gp_j,
                    gp_h,
                    e.1.isotope_offset,
                    gp_iso,
                ) + gp_cz * self.cz(&e.1)
            };
            let best = accepted_winners
                .iter()
                .map(|e| (e, fused_of(e)))
                .max_by(|(ea, sa), (eb, sb)| sa.total_cmp(sb).then_with(|| eb.0.cmp(&ea.0)))
                .map(|(e, _)| e);
            match best {
                Some((gl_key, w)) => {
                    let is_enum = deduped_backbone[w.bb_hit_idx].glycan.is_some();
                    if features_enumerated && !is_enum {
                        accepted_winners
                            .iter()
                            .filter(|e| deduped_backbone[e.1.bb_hit_idx].glycan.is_some())
                            .map(|e| (e, fused_of(e)))
                            .max_by(|(ea, sa), (eb, sb)| {
                                sa.total_cmp(sb).then_with(|| eb.0.cmp(&ea.0))
                            })
                            .map(|(e, _)| vec![(e.0, e.1)])
                            .unwrap_or_default()
                    } else {
                        vec![(*gl_key, *w)]
                    }
                }
                None => Vec::new(),
            }
            .into_iter()
            .filter(|(_, w)| core_y_counts[w.bb_hit_idx] as u32 >= min_core_y)
            .collect::<Vec<_>>()
        } else {
            let mut scored: Vec<(f32, (GlycanWinnerKey, CheapWinner))> = accepted_winners
                .into_iter()
                .map(|e| {
                    let cy = core_y_counts[e.1.bb_hit_idx] as f32;
                    let s = glyco_gp_fused_score_iso(
                        self.rank_sel(&e.1),
                        self.ladder(&e.1),
                        cy,
                        self.hyper(&e.1),
                        gp_k,
                        gp_j,
                        gp_h,
                        e.1.isotope_offset,
                        gp_iso,
                    ) + gp_cz * self.cz(&e.1);
                    (s, e)
                })
                .collect();
            scored.sort_by(|a, b| b.0.total_cmp(&a.0).then_with(|| (a.1).0.cmp(&(b.1).0)));
            scored.truncate(max_features);
            scored.into_iter().map(|(_, e)| e).collect()
        }
    }

    pub(super) fn full_hit(&self, w: CheapWinner, spread: &WinnerSpread) -> FullGlycoPsm {
        let ctx = self.ctx;
        let spec_idx = self.spec_idx;
        let candidates = self.ctx.candidates;
        let tol_ppm = self.ctx.tol_ppm;
        let fragment_tolerance_da = self.ctx.fragment_tolerance_da;
        let deduped_backbone = self.deduped_backbone;
        let core_y_counts = self.core_y_counts;
        let scored_per_charge = self.scored_per_charge;
        let y_peaks = self.y_peaks;
        let y_stats = self.y_stats;
        let gen_peaks = self.gen_peaks;
        let spec_scorer = self.spec_scorer;
        let ox_ev = self.ox_ev;
        let is_etd = self.is_etd;
        let pair_rank_etd = self.pair_rank_etd;
        let max_frag_charge = self.max_frag_charge;
        let WinnerSpread {
            tailor_denom,
            spectrum_listwise_gap,
            spectrum_rank_entropy,
            spectrum_delta_raw,
        } = *spread;
        let bb_hit = &deduped_backbone[w.bb_hit_idx];
        let bb_residue = bb_hit.backbone_mass;
        let bb_neutral = bb_residue + H2O;

        let ss = scored_for_charge(scored_per_charge, w.z)
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
            rank_score: if pair_rank_etd {
                self.rank_sel(&w)
            } else {
                w.rank
            },
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
        let cz_struct_vals = self.cz_struct(&w);
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
            cz_hyperscore: self.cz(&w),
            cz_intensity: self.cz_int(&w),
            cz_explained: cz_struct_vals.0,
            cz_chance_llr: cz_struct_vals.1,
        };
        FullGlycoPsm { glycan_key, psm }
    }
}
