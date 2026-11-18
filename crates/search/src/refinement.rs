use std::collections::HashSet;

use model::aa_set::{AminoAcidSet, AminoAcidSetBuilder};
use model::modification::{ModLocation, Modification, ResidueSpec};
use model::peptide::Peptide;
use model::protein::ProteinDb;
use scoring_crate::RankScorer;

use crate::candidate_gen::Candidate;
use crate::decoy::{build_search_db, whole_sequence_decoy_db, DecoyStrategy};
use crate::match_engine::PreparedSearch;
use crate::psm::TopNQueue;
use crate::refine_config::RefineConfig;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use crate::tdc::{confident_target_indices, ScoredLabel};

#[allow(dead_code)]
pub fn confident_protein_indices(
    queues: &[TopNQueue],
    candidates: &[Candidate],
    select_q: f64,
) -> HashSet<usize> {
    let mut labels: Vec<ScoredLabel> = Vec::new();
    let mut cand_idx_of: Vec<usize> = Vec::new();
    for queue in queues {
        if let Some(best) = queue.peek_top() {
            let cand_idx = best.primary_candidate_idx() as usize;
            labels.push(ScoredLabel {
                score: best.rank_score,
                is_decoy: candidates[cand_idx].is_decoy,
            });
            cand_idx_of.push(cand_idx);
        }
    }

    confident_target_indices(&labels, select_q)
        .into_iter()
        .map(|li| candidates[cand_idx_of[li]].protein_index)
        .collect()
}

pub fn confident_base_peptides(
    queues: &[TopNQueue],
    candidates: &[Candidate],
    select_q: f64,
) -> Vec<Vec<u8>> {
    confident_base_peptides_with_source(queues, candidates, select_q)
        .into_iter()
        .map(|(seq, _)| seq)
        .collect()
}

pub fn confident_base_peptides_with_source(
    queues: &[TopNQueue],
    candidates: &[Candidate],
    select_q: f64,
) -> Vec<(Vec<u8>, usize)> {
    let mut labels: Vec<ScoredLabel> = Vec::new();
    let mut cand_idx_of: Vec<usize> = Vec::new();
    for queue in queues {
        if let Some(best) = queue.peek_top() {
            let cand_idx = best.primary_candidate_idx() as usize;
            labels.push(ScoredLabel {
                score: best.rank_score,
                is_decoy: candidates[cand_idx].is_decoy,
            });
            cand_idx_of.push(cand_idx);
        }
    }

    let mut seqs: std::collections::HashMap<Vec<u8>, usize> = std::collections::HashMap::new();
    for li in confident_target_indices(&labels, select_q) {
        let cand = &candidates[cand_idx_of[li]];
        let backbone: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
        let pi = cand.protein_index;
        seqs.entry(backbone)
            .and_modify(|e| *e = (*e).min(pi))
            .or_insert(pi);
    }

    let mut out: Vec<(Vec<u8>, usize)> = seqs.into_iter().collect();
    out.sort();
    out
}

pub fn unidentified_spectrum_indices(
    queues: &[TopNQueue],
    candidates: &[Candidate],
    report_q: f64,
) -> Vec<usize> {
    let mut labels: Vec<ScoredLabel> = Vec::new();
    let mut label_spectrum: Vec<usize> = Vec::new();
    for (s, queue) in queues.iter().enumerate() {
        if let Some(best) = queue.peek_top() {
            let cand_idx = best.primary_candidate_idx() as usize;
            labels.push(ScoredLabel {
                score: best.rank_score,
                is_decoy: candidates[cand_idx].is_decoy,
            });
            label_spectrum.push(s);
        }
    }

    let confident = confident_target_indices(&labels, report_q);
    let identified: HashSet<usize> = confident.iter().map(|&i| label_spectrum[i]).collect();

    (0..queues.len())
        .filter(|s| !identified.contains(s))
        .collect()
}

#[allow(dead_code)]
pub fn build_refinement_index(
    full_target_db: &ProteinDb,
    confident_target_proteins: &HashSet<usize>,
    decoy_prefix: &str,
    decoy_strategy: DecoyStrategy,
    seed: u64,
) -> SearchIndex {
    let keep: HashSet<usize> = confident_target_proteins
        .iter()
        .copied()
        .filter(|&i| i < full_target_db.len())
        .collect();

    let subset = full_target_db.subset_by_index(&keep);
    let db = build_search_db(&subset, decoy_prefix, decoy_strategy, seed);
    SearchIndex {
        has_alternative_entries: crate::search_index::holds_alternative_entries(&db),
        db,
        decoy_prefix: crate::decoy::normalize_decoy_prefix(decoy_prefix),
        decoy_suffix: None,
    }
}

pub fn build_peptide_anchored_index(
    base_seqs: &[Vec<u8>],
    base_accessions: &[Option<String>],
    base_origins: &[Option<model::EntryOrigin>],
    decoy_prefix: &str,
    decoy_strategy: DecoyStrategy,
    seed: u64,
) -> SearchIndex {
    let minidb = ProteinDb {
        proteins: base_seqs
            .iter()
            .enumerate()
            .map(|(i, seq)| model::protein::Protein {
                accession: base_accessions
                    .get(i)
                    .and_then(|a| a.as_ref())
                    .filter(|acc| !acc.is_empty())
                    .cloned()
                    .unwrap_or_else(|| format!("BASEPEP_{i}")),
                description: String::new(),
                sequence: seq.clone(),
                annotation: base_origins.get(i).cloned().flatten().map(|o| {
                    std::sync::Arc::new(
                        model::ProteinAnnotation::new(Vec::new(), None).with_origin(Some(o)),
                    )
                }),
            })
            .collect(),
    };
    let combined = anchored_target_plus_decoy(&minidb, decoy_prefix, decoy_strategy, seed);
    SearchIndex {
        has_alternative_entries: crate::search_index::holds_alternative_entries(&combined),
        db: combined,
        decoy_prefix: crate::decoy::normalize_decoy_prefix(decoy_prefix),
        decoy_suffix: None,
    }
}

fn anchored_target_plus_decoy(
    minidb: &ProteinDb,
    decoy_prefix: &str,
    decoy_strategy: DecoyStrategy,
    seed: u64,
) -> ProteinDb {
    let interiors = ProteinDb {
        proteins: minidb
            .proteins
            .iter()
            .map(|p| model::protein::Protein {
                accession: p.accession.clone(),
                description: p.description.clone(),
                sequence: interior(&p.sequence).to_vec(),
                annotation: None,
            })
            .collect(),
    };
    let built = whole_sequence_decoy_db(&interiors, decoy_prefix, decoy_strategy, seed);
    let mut proteins = minidb.proteins.clone();
    proteins.extend(
        built
            .proteins
            .into_iter()
            .zip(minidb.proteins.iter())
            .map(|(mut d, t)| {
                let seq = &t.sequence;
                if seq.len() >= 2 {
                    let mut full = Vec::with_capacity(seq.len());
                    full.push(seq[0]);
                    full.extend_from_slice(&d.sequence);
                    full.push(seq[seq.len() - 1]);
                    d.sequence = full;
                }
                d.annotation = t.annotation.clone();
                d
            }),
    );
    ProteinDb { proteins }
}

fn anchors_at_protein_nterm(base_seqs: &[Vec<u8>], db: &ProteinDb) -> Vec<bool> {
    let starts: std::collections::HashSet<&[u8]> = db
        .proteins
        .iter()
        .flat_map(|p| {
            let s = p.sequence.as_slice();
            let met_cleaved = (s.first() == Some(&b'M') && s.len() > 1).then(|| &s[1..]);
            std::iter::once(s).chain(met_cleaved)
        })
        .collect();
    let lens: std::collections::BTreeSet<usize> = base_seqs.iter().map(|s| s.len()).collect();
    let prefixes: std::collections::HashSet<&[u8]> = lens
        .iter()
        .flat_map(|&l| {
            starts
                .iter()
                .filter(move |s| s.len() >= l)
                .map(move |s| &s[..l])
        })
        .collect();
    base_seqs
        .iter()
        .map(|a| prefixes.contains(a.as_slice()))
        .collect()
}

fn uses_protein_nterm(c: &Candidate, starts_with_met: bool) -> bool {
    (starts_with_met && c.start_offset_in_protein == 1)
        || c.peptide.residues.iter().any(|aa| {
            aa.mod_
                .as_ref()
                .is_some_and(|m| m.location == ModLocation::ProtNTerm)
        })
}

pub const ENTRAPMENT_PREFIX: &str = "ENT_";

fn entrapment_anchors(
    base_seqs: &[Vec<u8>],
    base_accs: &[Option<String>],
    seed: u64,
) -> (Vec<Vec<u8>>, Vec<Option<String>>, Vec<usize>) {
    let real: std::collections::HashSet<&[u8]> = base_seqs.iter().map(|s| s.as_slice()).collect();
    let mut seqs = Vec::new();
    let mut accs = Vec::new();
    let mut sources = Vec::new();
    for (i, seq) in base_seqs.iter().enumerate() {
        if seq.len() < 4 {
            continue;
        }
        let mut inner = interior(seq).to_vec();
        let mut state =
            (seed ^ 0xE47A_9C1D_5B3F_2E61) ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        if state == 0 {
            state = 0xD1B5_4A32_D192_ED03;
        }
        for j in (1..inner.len()).rev() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            inner.swap(j, (state % (j as u64 + 1)) as usize);
        }
        let mut ent = Vec::with_capacity(seq.len());
        ent.push(seq[0]);
        ent.extend_from_slice(&inner);
        ent.push(seq[seq.len() - 1]);
        if real.contains(ent.as_slice()) {
            continue;
        }
        let acc = base_accs
            .get(i)
            .and_then(|a| a.as_ref())
            .map_or_else(|| format!("BASEPEP_{i}"), |a| a.clone());
        seqs.push(ent);
        accs.push(Some(format!("{ENTRAPMENT_PREFIX}{acc}")));
        sources.push(i);
    }
    (seqs, accs, sources)
}

fn interior(seq: &[u8]) -> &[u8] {
    if seq.len() >= 2 {
        &seq[1..seq.len() - 1]
    } else {
        seq
    }
}

fn refinement_supported(decoy_strategy: DecoyStrategy) -> bool {
    !matches!(decoy_strategy, DecoyStrategy::None)
}

fn parse_location(loc: &str) -> ModLocation {
    match loc.trim().to_ascii_lowercase().as_str() {
        "n_term" | "n-term" | "nterm" => ModLocation::NTerm,
        "c_term" | "c-term" | "cterm" => ModLocation::CTerm,
        "protein_n_term" | "prot-n-term" | "prot_n_term" => ModLocation::ProtNTerm,
        "protein_c_term" | "prot-c-term" | "prot_c_term" => ModLocation::ProtCTerm,
        _ => ModLocation::Anywhere,
    }
}

fn parse_residue_spec(res: &str) -> ResidueSpec {
    let b = res.as_bytes();
    if b.len() == 1 && b[0].is_ascii_uppercase() {
        ResidueSpec::Specific(b[0])
    } else {
        ResidueSpec::Wildcard
    }
}

pub fn refinement_aa_set(
    base: &AminoAcidSet,
    cfg: &RefineConfig,
    high_res: bool,
) -> Result<AminoAcidSet, String> {
    let mut builder = AminoAcidSetBuilder::new_standard();

    let mod_key = |residue: ResidueSpec, location: ModLocation, delta: f64| {
        (residue, location, (delta * 1e4).round() as i64)
    };
    let mut present: HashSet<(ResidueSpec, ModLocation, i64)> = HashSet::new();
    let mut fixed_slots: Vec<(ResidueSpec, ModLocation)> = Vec::new();

    for m in base.distinct_mods() {
        if m.annotated_site_mass.is_some() {
            continue;
        }
        present.insert(mod_key(m.residue, m.location, m.mass_delta));
        if m.fixed {
            fixed_slots.push((m.residue, m.location));
            builder = builder.add_fixed_mod((*m).clone());
        } else {
            builder = builder.add_variable_mod((*m).clone());
        }
    }

    let residue_overlaps = |a: ResidueSpec, b: ResidueSpec| {
        a == ResidueSpec::Wildcard || b == ResidueSpec::Wildcard || a == b
    };

    for rm in &cfg.mods {
        if !high_res && rm.class == "deamidation" {
            continue;
        }
        let location = parse_location(&rm.location);
        for res in &rm.residues {
            let residue = parse_residue_spec(res);
            if fixed_slots
                .iter()
                .any(|&(fr, fl)| fl == location && residue_overlaps(fr, residue))
            {
                continue;
            }
            if !present.insert(mod_key(residue, location, rm.delta)) {
                continue;
            }
            let m = Modification {
                name: rm.name.clone(),
                mass_delta: rm.delta,
                residue,
                location,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
                annotated_site_mass: None,
                blocks_cleavage: false,
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
            };
            builder = builder.add_variable_mod(m);
        }
    }

    builder
        .with_label_channels(base.label_channels().to_vec())
        .build()
        .map_err(|e| format!("refinement_aa_set: {e}"))
}

pub fn mod_count_and_class(
    peptide: &Peptide,
    cfg: &RefineConfig,
    aa_set: &AminoAcidSet,
) -> (u32, u32) {
    let mut num_mods = 0u32;
    let mut class_id = 0u32;

    for aa in &peptide.residues {
        let Some(m) = aa.mod_ else { continue };
        let single = [m];
        let parts: &[&Modification] = if m.parts.is_empty() {
            aa_set.stack_parts(m).unwrap_or(&single)
        } else {
            &m.parts
        };
        for part in parts {
            if part.fixed {
                continue;
            }
            let matched = cfg.mods.iter().find(|rm| {
                (rm.delta - part.mass_delta).abs() < 1e-4
                    && parse_location(&rm.location) == part.location
                    && rm
                        .residues
                        .iter()
                        .any(|r| r == "*" || (r.len() == 1 && r.as_bytes()[0] == aa.residue))
            });
            let Some(rm) = matched else { continue };
            num_mods += 1;
            let this_class = RefineConfig::class_id(&rm.class);
            let pending = class_id == 0 || class_id == 99;
            if pending && (class_id == 0 || this_class != 99) {
                class_id = this_class;
            }
        }
    }

    (num_mods, class_id)
}

pub struct RefinementOutput {
    pub index: SearchIndex,
    pub candidates: Vec<Candidate>,
    pub queues: Vec<TopNQueue>,
    pub global_spectrum_indices: Vec<usize>,
}

#[allow(clippy::too_many_arguments)]
pub fn run_refinement(
    pass1_queues: &[TopNQueue],
    all_spectra: &[model::spectrum::Spectrum],
    pass1_candidates: &[Candidate],
    full_target_db: &ProteinDb,
    base_params: &SearchParams,
    scorer: &RankScorer,
    cfg: &RefineConfig,
    report_q: f64,
    high_res: bool,
    fragment_tol_da: f64,
    decoy_prefix: &str,
    decoy_strategy: DecoyStrategy,
    seed: u64,
) -> Option<RefinementOutput> {
    if !refinement_supported(decoy_strategy) {
        log::warn!(
            target: "refine",
            "skipped: --refine needs generated decoys for valid Pass-2 FDR, but \
             --decoy-strategy none was set (a target-only Pass-2 index would bias FDR). \
             Re-run with --decoy-strategy reverse|shuffle to enable refinement."
        );
        return None;
    }

    let base_with_src = confident_base_peptides_with_source(
        pass1_queues,
        pass1_candidates,
        base_params.refine_select_psm_fdr,
    );
    let mut base_seqs: Vec<Vec<u8>> = base_with_src.iter().map(|(seq, _)| seq.clone()).collect();
    if base_seqs.is_empty() {
        return None;
    }
    let mut base_accs: Vec<Option<String>> = base_with_src
        .iter()
        .map(|&(_, pi)| full_target_db.proteins.get(pi).map(|p| p.accession.clone()))
        .collect();
    let base_origins: Vec<Option<model::EntryOrigin>> = base_with_src
        .iter()
        .map(|&(_, pi)| {
            full_target_db
                .proteins
                .get(pi)
                .and_then(|p| p.origin())
                .filter(|o| o.kind != model::EntryKind::Reference)
                .cloned()
        })
        .collect();
    let mut anchor_nterm = anchors_at_protein_nterm(&base_seqs, full_target_db);
    if cfg.entrapment {
        let (ent_seqs, ent_accs, ent_sources) = entrapment_anchors(&base_seqs, &base_accs, seed);
        let ent_nterm: Vec<bool> = ent_sources.iter().map(|&i| anchor_nterm[i]).collect();
        anchor_nterm.extend(ent_nterm);
        log::info!(
            target: "refine",
            "entrapment: {} shuffled anchors added ({} anchors had no distinct shuffle)",
            ent_seqs.len(),
            base_seqs.len() - ent_seqs.len()
        );
        base_seqs.extend(ent_seqs);
        base_accs.extend(ent_accs);
    }

    let unident = unidentified_spectrum_indices(pass1_queues, pass1_candidates, report_q);
    if unident.is_empty() {
        return None;
    }

    let fmt_mod = |m: &crate::refine_config::RefineMod| {
        format!("{} ({} {:+.3})", m.name, m.residues.join(","), m.delta)
    };
    let active: Vec<String> = cfg
        .mods
        .iter()
        .filter(|m| high_res || m.class != "deamidation")
        .map(&fmt_mod)
        .collect();
    let skipped: Vec<String> = cfg
        .mods
        .iter()
        .filter(|m| !high_res && m.class == "deamidation")
        .map(&fmt_mod)
        .collect();
    log::info!(
        target: "refine",
        "Starting PTM-refinement search over {} confident-protein anchor peptides \
         ({} unidentified spectra, max {} mods/peptide)",
        base_seqs.len(),
        unident.len(),
        cfg.max_mods,
    );
    log::info!(target: "refine", "discovery mods: {}", active.join(", "));
    if !skipped.is_empty() {
        log::info!(target: "refine", "skipped (low-res): {}", skipped.join(", "));
    }

    let refine_idx = build_peptide_anchored_index(
        &base_seqs,
        &base_accs,
        &base_origins,
        decoy_prefix,
        decoy_strategy,
        seed,
    );

    const MAX_REFINE_MODS_PER_PEPTIDE: u32 = 5;
    let capped_max_mods = if cfg.max_mods > MAX_REFINE_MODS_PER_PEPTIDE {
        log::warn!(
            target: "refine",
            "--refine-max-mods {} exceeds the safe ceiling {}; clamping (combinatorial blowup guard).",
            cfg.max_mods, MAX_REFINE_MODS_PER_PEPTIDE
        );
        MAX_REFINE_MODS_PER_PEPTIDE
    } else {
        cfg.max_mods
    };
    let mut refine_params = base_params.clone();
    refine_params.aa_set = match refinement_aa_set(&base_params.aa_set, cfg, high_res) {
        Ok(set) => set,
        Err(e) => {
            log::warn!(target: "refine", "--refine could not build the Pass-2 modification set ({e}); skipping refinement for this run.");
            return None;
        }
    };
    refine_params.max_variable_mods_per_peptide = capped_max_mods;
    refine_params.fragment_index_top_k = 0;

    let (pass2_queues, refine_candidates, pass2_base_peptides) = {
        let n_anchors = base_seqs.len();
        let starts_with_met: Vec<bool> =
            base_seqs.iter().map(|s| s.first() == Some(&b'M')).collect();
        let mut refine_prepared = match PreparedSearch::prepare_in_memory(
            &refine_idx,
            &refine_params,
            scorer,
            fragment_tol_da,
            decoy_prefix,
        ) {
            Ok(p) => p,
            Err(e) => {
                log::warn!(
                    target: "refine",
                    "--refine could not build the Pass-2 candidate index ({e}); \
                     skipping refinement for this run."
                );
                return None;
            }
        }
        .with_candidate_filter(Box::new(move |c: &Candidate| {
            let a = c.protein_index % n_anchors;
            anchor_nterm[a] || !uses_protein_nterm(c, starts_with_met[a])
        }));
        let subset: Vec<model::spectrum::Spectrum> =
            unident.iter().map(|&i| all_spectra[i].clone()).collect();
        let queues = refine_prepared.run_chunk(&subset, 0);
        refine_prepared.sync_materialized_candidates();
        let base_peptides = refine_prepared.mmap_index.len();
        (
            queues,
            std::mem::take(&mut refine_prepared.candidates),
            base_peptides,
        )
    };

    let pass2_total = refine_candidates.len();
    let pass2_decoys = refine_candidates.iter().filter(|c| c.is_decoy).count();
    let mut isotope_mimics = 0usize;
    let mut compact: Vec<Candidate> = Vec::new();
    let mut remap: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
    let mut out_queues: Vec<TopNQueue> = Vec::with_capacity(unident.len());
    for (local_idx, pass2_queue) in pass2_queues.into_iter().enumerate() {
        let global_idx = unident[local_idx];
        let mut tagged = TopNQueue::new(refine_params.top_n_psms_per_spectrum);
        for mut psm in pass2_queue.into_rank_sorted_vec() {
            let (num_mods, class) = mod_count_and_class(
                &refine_candidates[psm.primary_candidate_idx() as usize].peptide,
                cfg,
                &refine_params.aa_set,
            );
            if num_mods == 0 {
                continue;
            }
            let pep = &refine_candidates[psm.primary_candidate_idx() as usize].peptide;
            if crate::isotope_mimic::fits_better(
                pep,
                crate::isotope_mimic::mod_shift_da(pep),
                psm.mass_error_ppm,
            ) {
                isotope_mimics += 1;
                continue;
            }
            psm.features.is_refinement = 1;
            psm.features.num_mods = num_mods;
            psm.features.refine_mod_class = class;
            psm.spectrum_idx = global_idx;
            for idx in &mut psm.candidate_idxs {
                let old = *idx;
                *idx = match remap.get(&old) {
                    Some(&n) => n,
                    None => {
                        let n = compact.len() as u32;
                        compact.push(refine_candidates[old as usize].clone());
                        remap.insert(old, n);
                        n
                    }
                };
            }
            tagged.force_push(psm);
        }
        out_queues.push(tagged);
    }
    drop(refine_candidates);

    log::debug!(
        target: "refine",
        "anchors={} unident={}/{} pass2_base_peptides={} pass2_candidates_retained={} (targets={} decoys={}) (pruned to {} winners; {} dropped as isotope mimics)",
        base_seqs.len(),
        unident.len(),
        all_spectra.len(),
        pass2_base_peptides,
        pass2_total,
        pass2_total - pass2_decoys,
        pass2_decoys,
        compact.len(),
        isotope_mimics,
    );

    Some(RefinementOutput {
        index: refine_idx,
        candidates: compact,
        queues: out_queues,
        global_spectrum_indices: unident,
    })
}

pub struct MergedSearch {
    pub candidates: Vec<Candidate>,
    pub index: SearchIndex,
}

pub fn merge_into_pass1(
    pass1_queues: &mut [TopNQueue],
    pass1_candidates: Vec<Candidate>,
    pass1_index: &SearchIndex,
    refine: RefinementOutput,
) -> MergedSearch {
    let cand_offset = pass1_candidates.len() as u32;
    let prot_offset = pass1_index.db.len();

    let mut combined: Vec<Candidate> = pass1_candidates;
    let RefinementOutput {
        index: refine_index,
        candidates: refine_cands,
        queues,
        global_spectrum_indices,
    } = refine;

    for mut c in refine_cands {
        c.protein_index += prot_offset;
        combined.push(c);
    }

    for (local_i, mut q) in queues.into_iter().enumerate() {
        let global = global_spectrum_indices[local_i];
        for mut psm in q.drain_into_vec() {
            for idx in &mut psm.candidate_idxs {
                *idx += cand_offset;
            }
            pass1_queues[global].force_push(psm);
        }
    }

    MergedSearch {
        candidates: combined,
        index: pass1_index.concat(&refine_index),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::psm::{PsmFeatures, PsmMatch};
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use model::protein::Protein;

    fn cand(protein_index: usize, is_decoy: bool) -> Candidate {
        let residues = vec![AminoAcid::standard(b'A').unwrap()];
        Candidate {
            peptide: Peptide::new(residues, b'_', b'-'),
            protein_index,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn psm(spectrum_idx: usize, cand_idx: u32, rank_score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![cand_idx],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: rank_score,
            rank_score,
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        }
    }

    fn queue_with(psms: Vec<PsmMatch>) -> TopNQueue {
        let mut q = TopNQueue::new(10);
        for p in psms {
            q.push(p);
        }
        q
    }

    fn protein(acc: &str, seq: &[u8]) -> Protein {
        Protein {
            accession: acc.into(),
            description: String::new(),
            sequence: seq.to_vec(),
            annotation: None,
        }
    }

    #[test]
    fn confident_protein_indices_keeps_targets_above_q_threshold() {
        let mut candidates: Vec<Candidate> = (0..200).map(|i| cand(i, false)).collect();
        candidates.push(cand(999, true));

        let mut queues = Vec::new();
        for i in 0..200 {
            queues.push(queue_with(vec![psm(i, i as u32, 30.0 - i as f32 * 0.05)]));
        }
        queues.push(queue_with(vec![psm(200, 200, 1.0)]));

        let confident = confident_protein_indices(&queues, &candidates, 0.01);
        assert!(!confident.is_empty());
        assert!(!confident.contains(&999), "decoy protein must be excluded");
        assert!(confident.contains(&0));
    }

    #[test]
    fn confident_protein_indices_uses_best_psm_per_spectrum() {
        let candidates = vec![cand(7, false), cand(8, true)];
        let q = queue_with(vec![psm(0, 0, 50.0), psm(0, 1, 5.0)]);
        let candidates2 = {
            let mut c = candidates.clone();
            c.push(cand(9, false));
            c
        };
        let q2 = queue_with(vec![psm(1, 2, 49.0)]);
        let confident = confident_protein_indices(&[q, q2], &candidates2, 0.5);
        assert!(
            confident.contains(&7),
            "best PSM (target on protein 7) selected"
        );
        assert!(
            !confident.contains(&8),
            "the worse decoy PSM must not be selected"
        );
    }

    #[test]
    fn confident_protein_indices_empty_when_all_decoys() {
        let candidates: Vec<Candidate> = (0..20).map(|i| cand(i, true)).collect();
        let queues: Vec<TopNQueue> = (0..20)
            .map(|i| queue_with(vec![psm(i, i as u32, 10.0)]))
            .collect();
        assert!(confident_protein_indices(&queues, &candidates, 0.01).is_empty());
    }

    fn cand_pep(peptide: Peptide, protein_index: usize, is_decoy: bool) -> Candidate {
        Candidate {
            peptide,
            protein_index,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn unmod_peptide(seq: &[u8]) -> Peptide {
        let residues: Vec<_> = seq
            .iter()
            .map(|&b| AminoAcid::standard(b).unwrap())
            .collect();
        Peptide::new(residues, b'_', b'-')
    }

    #[test]
    fn confident_base_peptides_returns_backbone_of_confident_targets() {
        let mut candidates: Vec<Candidate> = Vec::new();
        let mut queues: Vec<TopNQueue> = Vec::new();
        for i in 0..40 {
            let seq: &[u8] = if i % 2 == 0 { b"PEPTIK" } else { b"SAMPLER" };
            candidates.push(cand_pep(unmod_peptide(seq), i, false));
            queues.push(queue_with(vec![psm(i, i as u32, 90.0 - i as f32 * 0.1)]));
        }
        for k in 0..39 {
            let ci = candidates.len();
            candidates.push(cand_pep(unmod_peptide(b"DECYSEAR"), 200 + k, true));
            queues.push(queue_with(vec![psm(40 + k, ci as u32, 5.0)]));
        }
        let low_ci = candidates.len();
        candidates.push(cand_pep(unmod_peptide(b"LWSCAGER"), 300, false));
        queues.push(queue_with(vec![psm(79, low_ci as u32, 4.0)]));

        let seqs = confident_base_peptides(&queues, &candidates, 0.01);
        assert_eq!(seqs, vec![b"PEPTIK".to_vec(), b"SAMPLER".to_vec()]);
    }

    #[test]
    fn confident_base_peptides_ignores_modifications_in_backbone() {
        let ox_m = AminoAcid::standard(b'M').unwrap().with_mod(Modification {
            name: "Oxidation".into(),
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
        });
        let pep = Peptide::new(
            vec![
                AminoAcid::standard(b'P').unwrap(),
                ox_m,
                AminoAcid::standard(b'K').unwrap(),
            ],
            b'_',
            b'-',
        );
        let candidates = vec![
            cand_pep(pep, 0, false),
            cand_pep(unmod_peptide(b"AAAK"), 1, false),
        ];
        let queues = vec![
            queue_with(vec![psm(0, 0, 50.0)]),
            queue_with(vec![psm(1, 1, 49.0)]),
        ];
        let seqs = confident_base_peptides(&queues, &candidates, 0.5);
        assert!(
            seqs.contains(&b"PMK".to_vec()),
            "modified backbone stripped to PMK; got {seqs:?}"
        );
        assert!(seqs.contains(&b"AAAK".to_vec()));
    }

    #[test]
    fn confident_base_peptides_dedups_repeated_backbones() {
        let candidates = vec![
            cand_pep(unmod_peptide(b"PEPTIK"), 0, false),
            cand_pep(unmod_peptide(b"PEPTIK"), 1, false),
        ];
        let queues = vec![
            queue_with(vec![psm(0, 0, 50.0)]),
            queue_with(vec![psm(1, 1, 49.0)]),
        ];
        let seqs = confident_base_peptides(&queues, &candidates, 0.5);
        assert_eq!(seqs, vec![b"PEPTIK".to_vec()]);
    }

    #[test]
    fn confident_base_peptides_empty_when_all_decoys() {
        let candidates: Vec<Candidate> = (0..20)
            .map(|i| cand_pep(unmod_peptide(b"DECYSK"), i, true))
            .collect();
        let queues: Vec<TopNQueue> = (0..20)
            .map(|i| queue_with(vec![psm(i, i as u32, 10.0)]))
            .collect();
        assert!(confident_base_peptides(&queues, &candidates, 0.01).is_empty());
    }

    #[test]
    fn unidentified_spectrum_indices_returns_decoy_above_q_and_empty() {
        let mut candidates: Vec<Candidate> = (0..50).map(|i| cand(i, false)).collect();
        candidates.push(cand(999, true));

        let mut queues = Vec::new();
        for i in 0..50 {
            queues.push(queue_with(vec![psm(i, i as u32, 80.0 - i as f32 * 0.1)]));
        }
        queues.push(queue_with(vec![psm(50, 50, 1.0)]));
        queues.push(TopNQueue::new(10));

        let unident = unidentified_spectrum_indices(&queues, &candidates, 0.01);
        for s in 0..50 {
            assert!(
                !unident.contains(&s),
                "confident target spectrum {s} must be identified"
            );
        }
        assert!(
            unident.contains(&50),
            "decoy-best spectrum must be unidentified"
        );
        assert!(
            unident.contains(&51),
            "empty-queue spectrum must be unidentified"
        );
        assert_eq!(unident.len(), 2);
    }

    #[test]
    fn unidentified_spectrum_indices_returns_poorly_scoring_target_above_q() {
        let mut candidates: Vec<Candidate> = (0..40).map(|i| cand(i, false)).collect();
        for _ in 0..39 {
            candidates.push(cand(999, true));
        }
        candidates.push(cand(40, false));

        let mut queues = Vec::new();
        for i in 0..40 {
            queues.push(queue_with(vec![psm(i, i as u32, 90.0 - i as f32 * 0.1)]));
        }
        for k in 0..39 {
            queues.push(queue_with(vec![psm(40 + k, (40 + k) as u32, 5.0)]));
        }
        queues.push(queue_with(vec![psm(79, 79, 4.0)]));

        let unident = unidentified_spectrum_indices(&queues, &candidates, 0.01);
        for s in 0..40 {
            assert!(
                !unident.contains(&s),
                "confident target spectrum {s} must be identified"
            );
        }
        assert!(
            unident.contains(&79),
            "low-scoring target above report_q must be unidentified"
        );
    }

    #[test]
    fn unidentified_spectrum_indices_identifies_best_target_in_mixed_queue() {
        let candidates = vec![cand(0, false), cand(1, true), cand(2, false)];
        let queues = vec![
            queue_with(vec![psm(0, 0, 50.0), psm(0, 1, 10.0)]),
            queue_with(vec![psm(1, 2, 49.0)]),
        ];
        let unident = unidentified_spectrum_indices(&queues, &candidates, 0.5);
        assert!(
            unident.is_empty(),
            "both best-target spectra are identified"
        );
    }

    #[test]
    fn build_refinement_index_subsets_targets_and_pairs_decoys() {
        let full_target = ProteinDb {
            proteins: vec![
                protein("P0", b"MKWVR"),
                protein("P1", b"AGCTR"),
                protein("P2", b"LDESR"),
            ],
        };
        let keep: HashSet<usize> = [0usize, 2].into_iter().collect();
        let idx = build_refinement_index(&full_target, &keep, "XXX", DecoyStrategy::Reverse, 42);

        assert_eq!(idx.db.len(), 4);
        assert_eq!(idx.db.proteins[0].accession, "P0");
        assert_eq!(idx.db.proteins[1].accession, "P2");
        let needle = crate::decoy::decoy_accession_needle("XXX");
        assert!(idx.db.proteins[2].accession.starts_with(&needle));
        assert!(idx.db.proteins[3].accession.starts_with(&needle));
        assert_eq!(idx.decoy_prefix, "XXX");
    }

    #[test]
    fn build_refinement_index_drops_out_of_range_decoy_positions() {
        let full_target = ProteinDb {
            proteins: vec![protein("P0", b"MKWVR")],
        };
        let keep: HashSet<usize> = [0usize, 5].into_iter().collect();
        let idx = build_refinement_index(&full_target, &keep, "XXX", DecoyStrategy::Reverse, 42);
        assert_eq!(idx.db.len(), 2);
        assert_eq!(idx.db.proteins[0].accession, "P0");
    }

    #[test]
    fn anchored_mini_proteins_keep_their_source_origin() {
        let origin = model::EntryOrigin {
            kind: model::EntryKind::Variant,
            family: std::sync::Arc::from("P1"),
            feature: std::sync::Arc::from("VAR_1"),
            variant_delta_da: 0.984016,
            knockoff_of: None,
            source: None,
            tier: model::PriorTier::Catalogue,
        };
        let base_seqs = vec![b"PEPTIDK".to_vec(), b"SAMPLER".to_vec()];
        let idx = build_peptide_anchored_index(
            &base_seqs,
            &[],
            &[Some(origin.clone()), None],
            "XXX",
            DecoyStrategy::Reverse,
            42,
        );
        assert_eq!(idx.db.proteins[0].origin(), Some(&origin));
        assert_eq!(idx.db.proteins[1].origin(), None);
        assert_eq!(idx.db.proteins[2].origin(), Some(&origin));
        assert_eq!(idx.db.proteins[3].origin(), None);
        assert!(idx.has_alternative_entries);
    }

    #[test]
    fn build_peptide_anchored_index_one_mini_protein_per_base_seq_plus_decoys() {
        let base_seqs = vec![b"PEPTIK".to_vec(), b"SAMPLER".to_vec()];
        let idx =
            build_peptide_anchored_index(&base_seqs, &[], &[], "XXX", DecoyStrategy::Reverse, 42);

        assert_eq!(idx.db.len(), 4);
        assert_eq!(idx.db.proteins[0].accession, "BASEPEP_0");
        assert_eq!(idx.db.proteins[0].sequence, b"PEPTIK".to_vec());
        assert_eq!(idx.db.proteins[1].accession, "BASEPEP_1");
        assert_eq!(idx.db.proteins[1].sequence, b"SAMPLER".to_vec());
        let needle = crate::decoy::decoy_accession_needle("XXX");
        assert!(idx.db.proteins[2].accession.starts_with(&needle));
        assert!(idx.db.proteins[3].accession.starts_with(&needle));
        assert_eq!(idx.decoy_prefix, "XXX");
    }

    #[test]
    fn anchored_decoys_keep_termini_and_pair_full_length_with_targets() {
        let anchors: Vec<&[u8]> = vec![
            b"PEPTIDEK",
            b"PEPKTIDEK",
            b"QPEPTIDEK",
            b"EAMPLENR",
            b"SAMNPLEQR",
            b"ACDKNQMR",
        ];
        let base_seqs: Vec<Vec<u8>> = anchors.iter().map(|a| a.to_vec()).collect();
        let idx =
            build_peptide_anchored_index(&base_seqs, &[], &[], "XXX", DecoyStrategy::Reverse, 42);
        let n = anchors.len();
        assert_eq!(idx.db.proteins[n].sequence, b"PEDITPEK".to_vec());
        assert_eq!(idx.db.proteins[n + 1].sequence, b"PEDITKPEK".to_vec());

        let base = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .build()
            .unwrap();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        let mut params = SearchParams::default_tryptic(set);
        params.min_length = 6;
        params.max_variable_mods_per_peptide = 2;
        let cands: Vec<Candidate> =
            crate::candidate_gen::enumerate_candidates(&idx, &params, "XXX").collect();
        for (i, a) in anchors.iter().enumerate() {
            let full = |prot: usize| {
                cands
                    .iter()
                    .filter(|c| c.protein_index == prot && c.peptide.length() == a.len())
                    .count()
            };
            let (t, d) = (full(i), full(i + n));
            assert!(
                t > 0,
                "{}: no full-length target",
                String::from_utf8_lossy(a)
            );
            assert_eq!(
                t,
                d,
                "{}: full-length target forms {t} vs decoy forms {d}",
                String::from_utf8_lossy(a)
            );
        }
    }

    fn base_set_with_cam_c() -> AminoAcidSet {
        AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .build()
            .unwrap()
    }

    #[test]
    fn refinement_aa_set_leaves_out_mods_restricted_to_annotated_sites() {
        let mut ox = model::mod_spec::parse_mod_list("Oxidation (M)", false)
            .unwrap()
            .remove(0);
        ox.annotated_site_mass = Some(ox.mass_delta);
        let base = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .add_variable_mod(ox)
            .build()
            .unwrap();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        assert!(set
            .iter_variants()
            .all(|aa| aa.mod_.is_none_or(|m| m.annotated_site_mass.is_none())));
        assert!(set
            .variants_for(b'M', ModLocation::Anywhere)
            .iter()
            .any(|aa| aa.mod_.is_some_and(|m| m.name == "Oxidation")));
    }

    #[test]
    fn refinement_aa_set_high_res_exposes_five_variable_mods() {
        let base = base_set_with_cam_c();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();

        let mut names: HashSet<String> = HashSet::new();
        for aa in set.iter_variants() {
            if let Some(m) = &aa.mod_ {
                if !m.fixed {
                    names.insert(m.name.rsplit('+').next().unwrap_or(&m.name).to_string());
                }
            }
        }
        assert_eq!(
            names.len(),
            5,
            "default tier has 5 variable mods at high-res; got {names:?}"
        );
        let c_anywhere = set.variants_for(b'C', ModLocation::Anywhere);
        assert!(c_anywhere
            .iter()
            .any(|aa| aa.mod_.as_ref().map(|m| m.fixed).unwrap_or(false)));
    }

    #[test]
    fn refinement_aa_set_low_res_drops_deamidation() {
        let base = base_set_with_cam_c();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), false).unwrap();
        let mut names: HashSet<String> = HashSet::new();
        for aa in set.iter_variants() {
            if let Some(m) = &aa.mod_ {
                if !m.fixed {
                    names.insert(m.name.rsplit('+').next().unwrap_or(&m.name).to_string());
                }
            }
        }
        assert_eq!(
            names.len(),
            4,
            "deamidation dropped at low-res; got {names:?}"
        );
        assert!(
            !names.contains("Deamidation"),
            "Deamidation must be dropped at low-res"
        );
    }

    #[test]
    fn refinement_aa_set_oxidation_m_has_modified_variant() {
        let base = base_set_with_cam_c();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        let m_variants = set.variants_for(b'M', ModLocation::Anywhere);
        assert!(m_variants.iter().any(|aa| {
            aa.mod_
                .as_ref()
                .map(|m| m.name == "Oxidation")
                .unwrap_or(false)
        }));
    }

    #[test]
    fn refinement_aa_set_excludes_stacked_combined_base_variant() {
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
        let base = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .add_variable_mod(acetyl)
            .build()
            .unwrap();
        assert!(
            base.distinct_mods().iter().all(|m| !m.name.contains('+')),
            "distinct_mods leaked a synthesized stacked (fixed+variable) variant"
        );
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        let a_nterm = set.variants_for(b'A', ModLocation::ProtNTerm);
        assert!(
            a_nterm.iter().all(|aa| aa
                .mod_
                .as_ref()
                .map(|m| (m.mass_delta - 99.032).abs() >= 1e-2)
                .unwrap_or(true)),
            "phantom ~99 Da stacked mod leaked onto non-Cys N-term residue"
        );
    }

    #[test]
    fn refinement_aa_set_inherits_nonstandard_fixed_tmt() {
        let tmt_k = Modification {
            name: "TMT6plex".into(),
            mass_delta: 229.162932,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let base = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .add_fixed_mod(tmt_k)
            .build()
            .unwrap();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        let k_anywhere = set.variants_for(b'K', ModLocation::Anywhere);
        assert!(
            k_anywhere.iter().any(|aa| aa
                .mod_
                .as_ref()
                .map(|m| m.fixed && (m.mass_delta - 229.162932).abs() < 1e-3)
                .unwrap_or(false)),
            "refined set must carry the base's fixed TMT6plex-on-K"
        );
        let c_anywhere = set.variants_for(b'C', ModLocation::Anywhere);
        assert!(c_anywhere
            .iter()
            .any(|aa| aa.mod_.as_ref().map(|m| m.fixed).unwrap_or(false)));
    }

    #[test]
    fn refinement_aa_set_dedups_base_variable_against_tier() {
        let ox_m = Modification {
            name: "Oxidation".into(),
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
        };
        let base = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .add_variable_mod(ox_m)
            .build()
            .unwrap();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        let ox_count = set
            .variants_for(b'M', ModLocation::Anywhere)
            .iter()
            .filter(|aa| {
                aa.mod_
                    .as_ref()
                    .map(|m| !m.fixed && (m.mass_delta - 15.994915).abs() < 1e-3)
                    .unwrap_or(false)
            })
            .count();
        assert_eq!(
            ox_count, 1,
            "Oxidation-M must not be double-added; got {ox_count}"
        );
    }

    fn std_set() -> AminoAcidSet {
        AminoAcidSetBuilder::new_standard().build().unwrap()
    }

    fn modded_pep(residue: u8, delta: f64) -> Peptide {
        let m = Modification {
            name: "x".into(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(residue),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let aa = AminoAcid::standard(residue).unwrap().with_mod(m);
        Peptide::new(vec![aa, AminoAcid::standard(b'G').unwrap()], b'_', b'-')
    }

    #[test]
    fn anchors_at_protein_nterm_sees_position_zero_and_initiator_met() {
        let db = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: String::new(),
                    sequence: b"MSAMPLERPEPTIDEK".to_vec(),
                    annotation: None,
                },
                Protein {
                    accession: "P2".into(),
                    description: String::new(),
                    sequence: b"QWERTYKLMNPQR".to_vec(),
                    annotation: None,
                },
            ],
        };
        let anchors: Vec<Vec<u8>> = vec![
            b"MSAMPLER".to_vec(),
            b"SAMPLER".to_vec(),
            b"PEPTIDEK".to_vec(),
            b"QWERTYK".to_vec(),
            b"WERTYK".to_vec(),
        ];
        assert_eq!(
            anchors_at_protein_nterm(&anchors, &db),
            vec![true, true, false, true, false]
        );
    }

    #[test]
    fn uses_protein_nterm_flags_prot_nterm_mods_and_met_cleavage_only() {
        let cand = |pep: Peptide, offset: usize| Candidate {
            peptide: pep,
            protein_index: 0,
            start_offset_in_protein: offset,
            is_decoy: false,
            is_protein_n_term: true,
            is_protein_c_term: true,
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
        let mut res = unmod_peptide(b"SAMPLER").residues;
        res[0] = res[0].clone().with_mod(acetyl);
        let acetylated = Peptide::new(res, b'_', b'-');
        assert!(uses_protein_nterm(&cand(acetylated, 0), false));
        assert!(!uses_protein_nterm(
            &cand(unmod_peptide(b"SAMPLER"), 0),
            false
        ));
        assert!(uses_protein_nterm(
            &cand(unmod_peptide(b"SAMPLER"), 1),
            true
        ));
        assert!(!uses_protein_nterm(
            &cand(unmod_peptide(b"SAMPLER"), 1),
            false
        ));
    }

    #[test]
    fn entrapment_anchors_keep_termini_and_composition_and_are_not_real() {
        let base: Vec<Vec<u8>> = vec![
            b"PEPTIDEK".to_vec(),
            b"SAMNPLEQR".to_vec(),
            b"AAAAK".to_vec(),
        ];
        let accs = vec![Some("P1".to_string()), None, Some("P3".to_string())];
        let (seqs, ent_accs, sources) = entrapment_anchors(&base, &accs, 1);
        assert_eq!(sources, vec![0, 1]);
        assert_eq!(seqs.len(), 2, "the all-A interior has no distinct shuffle");
        assert_eq!(ent_accs[0].as_deref(), Some("ENT_P1"));
        assert_eq!(ent_accs[1].as_deref(), Some("ENT_BASEPEP_1"));
        for (ent, orig) in seqs.iter().zip(&base) {
            assert_ne!(ent, orig);
            assert_eq!(ent.first(), orig.first());
            assert_eq!(ent.last(), orig.last());
            let (mut a, mut b) = (ent.clone(), orig.clone());
            a.sort_unstable();
            b.sort_unstable();
            assert_eq!(a, b, "same composition, so same precursor mass");
        }
        assert_eq!(entrapment_anchors(&base, &accs, 1).0, seqs);
    }

    #[test]
    fn mod_count_and_class_unmodified_is_zero_zero() {
        let cfg = RefineConfig::default_tier();
        let pep = Peptide::new(
            vec![
                AminoAcid::standard(b'P').unwrap(),
                AminoAcid::standard(b'G').unwrap(),
            ],
            b'_',
            b'-',
        );
        assert_eq!(mod_count_and_class(&pep, &cfg, &std_set()), (0, 0));
    }

    #[test]
    fn mod_count_and_class_oxidation_m_is_one_class_oxidation() {
        let cfg = RefineConfig::default_tier();
        let pep = modded_pep(b'M', 15.994915);
        assert_eq!(mod_count_and_class(&pep, &cfg, &std_set()), (1, 1));
    }

    #[test]
    fn mod_count_and_class_deamidation_n_is_class_two() {
        let cfg = RefineConfig::default_tier();
        let pep = modded_pep(b'N', 0.984016);
        assert_eq!(mod_count_and_class(&pep, &cfg, &std_set()), (1, 2));
    }

    #[test]
    fn mod_count_and_class_counts_a_refinement_mod_under_a_terminal_label() {
        let tmt_nterm = Modification {
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
        let base = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(tmt_nterm)
            .build()
            .unwrap();
        let cfg = RefineConfig::default_tier();
        let set = refinement_aa_set(&base, &cfg, true).unwrap();
        let forms = set.variants_for(b'M', ModLocation::NTerm);
        let delta = |aa: &AminoAcid| aa.mod_.map_or(0.0, |m| m.mass_delta);
        let labelled = forms
            .iter()
            .find(|aa| (delta(aa) - 229.162932).abs() < 1e-4)
            .expect("Met with the label only");
        let oxidised = forms
            .iter()
            .find(|aa| (delta(aa) - (229.162932 + 15.994915)).abs() < 1e-4)
            .expect("Met with the label and the oxidation");
        assert!(set.stack_parts(oxidised.mod_.unwrap()).is_some());
        let tail = AminoAcid::standard(b'G').unwrap();
        let pep = |first: &AminoAcid| Peptide::new(vec![first.clone(), tail.clone()], b'_', b'-');
        assert_eq!(mod_count_and_class(&pep(oxidised), &cfg, &set), (1, 1));
        assert_eq!(mod_count_and_class(&pep(labelled), &cfg, &set), (0, 0));
    }

    #[test]
    fn mod_count_and_class_skips_fixed_parts_of_a_stack() {
        let fixed = |name: &str, delta: f64, residue, location| Modification {
            name: name.into(),
            mass_delta: delta,
            residue,
            location,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
            parts: Vec::new(),
        };
        let base = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(fixed(
                "TMT6plex",
                229.162932,
                ResidueSpec::Wildcard,
                ModLocation::NTerm,
            ))
            .add_fixed_mod(fixed(
                "Oxidation",
                15.994915,
                ResidueSpec::Specific(b'M'),
                ModLocation::Anywhere,
            ))
            .build()
            .unwrap();
        let met = base
            .variants_for(b'M', ModLocation::NTerm)
            .iter()
            .find(|aa| {
                aa.mod_
                    .is_some_and(|m| (m.mass_delta - (229.162932 + 15.994915)).abs() < 1e-4)
            })
            .expect("N-terminal Met with both fixed mods")
            .clone();
        assert!(base.stack_parts(met.mod_.unwrap()).is_some());
        let tail = AminoAcid::standard(b'G').unwrap();
        let pep = Peptide::new(vec![met, tail], b'_', b'-');
        let cfg = RefineConfig::default_tier();
        assert_eq!(mod_count_and_class(&pep, &cfg, &base), (0, 0));
        let inner = base
            .variants_for(b'M', ModLocation::Anywhere)
            .iter()
            .find(|aa| {
                aa.mod_
                    .is_some_and(|m| (m.mass_delta - 15.994915).abs() < 1e-4)
            })
            .expect("internal Met with the fixed oxidation")
            .clone();
        assert!(inner.mod_.unwrap().parts.is_empty());
        let g = AminoAcid::standard(b'G').unwrap();
        let pep = Peptide::new(vec![g.clone(), inner, g], b'_', b'-');
        assert_eq!(mod_count_and_class(&pep, &cfg, &base), (0, 0));
    }

    #[test]
    fn mod_count_and_class_ignores_non_refinement_mods() {
        let cfg = RefineConfig::default_tier();
        let pep = modded_pep(b'C', 57.02146);
        assert_eq!(mod_count_and_class(&pep, &cfg, &std_set()), (0, 0));
    }

    #[test]
    fn mod_count_and_class_counts_two_distinct_mods() {
        let cfg = RefineConfig::default_tier();
        let ox = AminoAcid::standard(b'M').unwrap().with_mod(Modification {
            name: "Oxidation".into(),
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
        });
        let deam = AminoAcid::standard(b'N').unwrap().with_mod(Modification {
            name: "Deamidation".into(),
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
        });
        let pep = Peptide::new(vec![ox, deam], b'_', b'-');
        let (n, class) = mod_count_and_class(&pep, &cfg, &std_set());
        assert_eq!(n, 2);
        assert_eq!(class, 1);
    }

    #[test]
    fn run_refinement_tag_merge_preserves_pass1_and_tags_pass2() {
        let cfg = RefineConfig::default_tier();

        let pass1_cands = [cand(0, false)];
        let mut global_queues = [
            queue_with(vec![psm(0, 0, 40.0)]),
            TopNQueue::new(10),
        ];

        let pass2_cands = [Candidate {
            peptide: modded_pep(b'M', 15.994915),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }];
        let mut winner = psm(0, 0, 22.0);

        let global_idx = 1usize;
        winner.features.is_refinement = 1;
        let (num_mods, class) = mod_count_and_class(
            &pass2_cands[winner.primary_candidate_idx() as usize].peptide,
            &cfg,
            &std_set(),
        );
        winner.features.num_mods = num_mods;
        winner.features.refine_mod_class = class;
        winner.spectrum_idx = global_idx;
        global_queues[global_idx].force_push(winner);

        assert_eq!(global_queues[0].len(), 1);
        assert_eq!(
            global_queues[0].peek_top().unwrap().features.is_refinement,
            0
        );

        assert_eq!(global_queues[1].len(), 1);
        let merged = global_queues[1].peek_top().unwrap();
        assert_eq!(merged.features.is_refinement, 1);
        assert_eq!(merged.features.num_mods, 1);
        assert_eq!(merged.features.refine_mod_class, 1);
        assert_eq!(merged.spectrum_idx, 1);

        assert!(!pass1_cands[0].is_decoy);
    }

    fn tiny_scorer() -> RankScorer {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;
        use rustc_hash::FxHashMap;
        use scoring_crate::param_model::{IonType, Partition, SpecDataType};
        use scoring_crate::Param;

        let part = Partition {
            charge: 2,
            parent_mass: 1000.0,
            seg_num: 0,
        };
        let prefix1 = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let suffix1 = IonType::Suffix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let noise = IonType::Noise;

        let mut ion_table = FxHashMap::default();
        ion_table.insert(prefix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
        ion_table.insert(suffix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
        ion_table.insert(noise, vec![0.1_f32, 0.05, 0.02, 0.01]);

        let mut rank_dist_table = FxHashMap::default();
        rank_dist_table.insert(part, ion_table);
        let mut frag_off_table = FxHashMap::default();
        frag_off_table.insert(part, vec![]);

        let mut param = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: FxHashMap::default(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        };
        param.rebuild_cache();
        RankScorer::new(&param)
    }

    #[test]
    fn refinement_unsupported_without_generated_decoys() {
        assert!(!refinement_supported(DecoyStrategy::None));
        assert!(refinement_supported(DecoyStrategy::Reverse));
        assert!(refinement_supported(DecoyStrategy::Shuffle));
    }

    #[test]
    fn run_refinement_skips_when_decoy_strategy_none() {
        let aa_set = base_set_with_cam_c();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = tiny_scorer();
        let target = ProteinDb {
            proteins: vec![protein("P0", b"MKWVTFISLLR")],
        };
        let pass1_cands = vec![cand(0, false), cand(1, false)];
        let queues = vec![
            queue_with(vec![psm(0, 0, 80.0)]),
            queue_with(vec![psm(1, 1, 1.0)]),
        ];
        let spectra: Vec<model::spectrum::Spectrum> = vec![
            model::spectrum::Spectrum::default(),
            model::spectrum::Spectrum::default(),
        ];

        let out = run_refinement(
            &queues,
            &spectra,
            &pass1_cands,
            &target,
            &params,
            &scorer,
            &RefineConfig::default_tier(),
            0.01,
            true,
            0.5,
            "XXX",
            DecoyStrategy::None,
            42,
        );
        assert!(out.is_none());
    }

    #[test]
    fn run_refinement_early_returns_when_no_confident_base_peptides() {
        let aa_set = base_set_with_cam_c();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = tiny_scorer();
        let target = ProteinDb {
            proteins: vec![protein("P0", b"MKWVTFISLLR")],
        };

        let pass1_cands = vec![cand(0, true)];
        let queues = vec![queue_with(vec![psm(0, 0, 10.0)])];
        let before_len = queues[0].len();
        let spectra: Vec<model::spectrum::Spectrum> = vec![model::spectrum::Spectrum::default()];

        let out = run_refinement(
            &queues,
            &spectra,
            &pass1_cands,
            &target,
            &params,
            &scorer,
            &RefineConfig::default_tier(),
            0.01,
            true,
            0.5,
            "XXX",
            DecoyStrategy::Reverse,
            42,
        );
        assert!(out.is_none());
        assert_eq!(queues[0].len(), before_len);
    }

    #[test]
    fn run_refinement_early_returns_when_no_unidentified_spectra() {
        let aa_set = base_set_with_cam_c();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = tiny_scorer();
        let target = ProteinDb {
            proteins: vec![protein("P0", b"MKWVTFISLLR")],
        };

        let pass1_cands = vec![cand(0, false)];
        let queues = vec![queue_with(vec![psm(0, 0, 50.0)])];
        let spectra: Vec<model::spectrum::Spectrum> = vec![model::spectrum::Spectrum::default()];

        let out = run_refinement(
            &queues,
            &spectra,
            &pass1_cands,
            &target,
            &params,
            &scorer,
            &RefineConfig::default_tier(),
            0.01,
            true,
            0.5,
            "XXX",
            DecoyStrategy::Reverse,
            42,
        );
        assert!(out.is_none());
        assert_eq!(queues[0].len(), 1);
        assert_eq!(queues[0].peek_top().unwrap().features.is_refinement, 0);
    }

    fn fixture_residue_mass(aa: u8, cam: bool) -> f64 {
        match aa {
            b'G' => 57.021464,
            b'A' => 71.037114,
            b'S' => 87.032028,
            b'P' => 97.052764,
            b'V' => 99.068414,
            b'T' => 101.047679,
            b'C' if cam => 103.009185 + 57.021464,
            b'C' => 103.009185,
            b'N' => 114.042927,
            b'D' => 115.026943,
            b'Q' => 128.058578,
            b'K' => 128.094963,
            b'k' => 128.094963 + 15.994915,
            b'E' => 129.042593,
            b'M' => 131.040485,
            b'F' => 147.068414,
            b'W' => 186.079313,
            _ => panic!("residue {} not in the fixture table", aa as char),
        }
    }

    fn fixture_spectrum(seq: &[u8], cam: bool, extra_da: f64) -> model::spectrum::Spectrum {
        use model::mass::{H2O, PROTON};
        let res: Vec<f64> = seq.iter().map(|&a| fixture_residue_mass(a, cam)).collect();
        let neutral: f64 = res.iter().sum::<f64>() + H2O;
        let mut peaks: Vec<(f64, f32)> = Vec::new();
        for i in 1..res.len() {
            let b: f64 = res[..i].iter().sum();
            peaks.push((b + PROTON, 1000.0));
            peaks.push((neutral - b + PROTON, 1000.0));
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        model::spectrum::Spectrum {
            precursor_mz: (neutral + extra_da + 2.0 * PROTON) / 2.0,
            precursor_charge: Some(2),
            scan: Some(1),
            peaks,
            ..Default::default()
        }
    }

    fn refine_fixture(
        anchors: &[(&[u8], usize)],
        db: &ProteinDb,
        cam: bool,
        query: model::spectrum::Spectrum,
    ) -> Option<RefinementOutput> {
        let aa_set = if cam {
            base_set_with_cam_c()
        } else {
            AminoAcidSetBuilder::new_standard().build().unwrap()
        };
        let params = SearchParams::default_tryptic(aa_set);
        let mut candidates = Vec::new();
        let mut queues = Vec::new();
        let mut spectra = Vec::new();
        for (i, (seq, protein_index)) in anchors.iter().enumerate() {
            candidates.push(cand_pep(unmod_peptide(seq), *protein_index, false));
            queues.push(queue_with(vec![psm(i, i as u32, 90.0 - i as f32)]));
            spectra.push(model::spectrum::Spectrum::default());
        }
        queues.push(TopNQueue::new(10));
        spectra.push(query);
        run_refinement(
            &queues,
            &spectra,
            &candidates,
            db,
            &params,
            &tiny_scorer(),
            &RefineConfig::default_tier(),
            0.01,
            true,
            0.5,
            "XXX",
            DecoyStrategy::Reverse,
            42,
        )
    }

    #[test]
    fn refinement_collision_targets_respect_the_anchor_filter() {
        let db = ProteinDb {
            proteins: vec![protein("P1", b"GRAGFEDCKSSRMACDEFGKWR")],
        };
        let out = refine_fixture(
            &[(b"AGFEDCK", 0), (b"MACDEFGK", 0)],
            &db,
            true,
            fixture_spectrum(b"ACDEFGk", true, 0.0),
        )
        .expect("refinement emits PSMs");
        let hits: Vec<&Candidate> = out
            .queues
            .iter()
            .flat_map(|q| q.iter_psms())
            .map(|p| &out.candidates[p.primary_candidate_idx() as usize])
            .filter(|c| {
                c.peptide
                    .residues
                    .iter()
                    .map(|aa| aa.residue)
                    .collect::<Vec<_>>()
                    == b"ACDEFGK"
            })
            .collect();
        assert!(!hits.is_empty(), "the K-oxidised ACDEFGK PSM is emitted");
        for c in hits {
            assert!(
                c.is_decoy,
                "ACDEFGK comes only from the decoy of AGFEDCK here (protein {})",
                out.index.db.proteins[c.protein_index].accession
            );
        }
    }

    #[test]
    fn refinement_retrieves_exact_precursor_matches_outside_the_nominal_window() {
        let anchor: &[u8] = b"CEMEMSDDDGWSSFCEGNVNCGMDCK";
        let db = ProteinDb {
            proteins: vec![protein("P1", &[b"GR".as_slice(), anchor, b"AR"].concat())],
        };
        let out = refine_fixture(
            &[(anchor, 0)],
            &db,
            false,
            fixture_spectrum(
                b"CEMEMSDDDGWSSFCEGNVNCGMDCk",
                false,
                2.0 * model::mass::ISOTOPE + 0.02,
            ),
        );
        let found = out.iter().flat_map(|o| {
            o.queues
                .iter()
                .flat_map(|q| q.iter_psms())
                .map(move |p| (p, &o.candidates[p.primary_candidate_idx() as usize]))
        });
        let found = found.into_iter().any(|(p, c)| {
            !c.is_decoy
                && p.isotope_offset == 2
                && c.peptide
                    .residues
                    .iter()
                    .map(|aa| aa.residue)
                    .collect::<Vec<_>>()
                    == anchor
                && c.peptide
                    .residues
                    .last()
                    .is_some_and(|aa| aa.mod_.is_some())
        });
        assert!(
            found,
            "the K-oxidised target at isotope +2 is a refinement PSM"
        );
    }

    #[test]
    fn merge_offsets_indices_and_force_pushes_into_pass1() {
        let p1_cands = vec![cand(0, false)];
        let p1_index = SearchIndex {
            db: ProteinDb {
                proteins: vec![protein("P0", b"PEPTIDEK")],
            },
            has_alternative_entries: false,
            decoy_prefix: "XXX".into(),
            decoy_suffix: None,
        };
        let mut p1_queues = vec![queue_with(vec![psm(0, 0, 40.0)])];

        let refine = RefinementOutput {
            index: SearchIndex {
                db: ProteinDb {
                    proteins: vec![protein("BASEPEP_0", b"PEPTIDEK")],
                },
                has_alternative_entries: false,
                decoy_prefix: "XXX".into(),
                decoy_suffix: None,
            },
            candidates: vec![cand(0, false)],
            queues: vec![queue_with(vec![psm(0, 0, 22.0)])],
            global_spectrum_indices: vec![0],
        };

        let merged = merge_into_pass1(&mut p1_queues, p1_cands, &p1_index, refine);

        assert_eq!(merged.candidates.len(), 2);
        assert_eq!(
            merged.candidates[1].protein_index, 1,
            "refine protein_index offset by pass1 db len"
        );
        assert_eq!(merged.index.db.len(), 2);
        assert_eq!(merged.index.db.proteins[1].accession, "BASEPEP_0");
        let mut idxs: Vec<u32> = p1_queues[0]
            .iter_psms()
            .map(|m| m.primary_candidate_idx())
            .collect();
        idxs.sort();
        assert_eq!(
            idxs,
            vec![0, 1],
            "refine PSM candidate_idx offset to 1 and force_pushed into pass1 queue"
        );
    }

    #[test]
    fn merge_offsets_every_entry_of_multi_candidate_psms() {
        let p1_cands = vec![cand(0, false), cand(1, false)];
        let p1_index = SearchIndex {
            db: ProteinDb {
                proteins: vec![protein("P0", b"AAAK"), protein("P1", b"BBBK")],
            },
            has_alternative_entries: false,
            decoy_prefix: "XXX".into(),
            decoy_suffix: None,
        };
        let mut p1_queues = vec![queue_with(vec![psm(0, 0, 40.0)])];

        let mut shared = psm(0, 0, 22.0);
        shared.candidate_idxs = vec![0, 1];
        let refine = RefinementOutput {
            index: SearchIndex {
                db: ProteinDb {
                    proteins: vec![
                        protein("BASEPEP_0", b"PEPTIDEK"),
                        protein("BASEPEP_1", b"SAMPLEK"),
                    ],
                },
                has_alternative_entries: false,
                decoy_prefix: "XXX".into(),
                decoy_suffix: None,
            },
            candidates: vec![cand(0, false), cand(1, false)],
            queues: vec![queue_with(vec![shared])],
            global_spectrum_indices: vec![0],
        };

        let merged = merge_into_pass1(&mut p1_queues, p1_cands, &p1_index, refine);

        assert_eq!(merged.candidates.len(), 4);
        assert_eq!(merged.candidates[2].protein_index, 2);
        assert_eq!(merged.candidates[3].protein_index, 3);
        let refine_psm = p1_queues[0]
            .iter_psms()
            .find(|m| m.candidate_idxs.len() == 2)
            .expect("the multi-candidate refine PSM was force_pushed");
        let mut idxs = refine_psm.candidate_idxs.clone();
        idxs.sort();
        assert_eq!(
            idxs,
            vec![2, 3],
            "every candidate_idx entry offset, not just the primary"
        );
    }
}
