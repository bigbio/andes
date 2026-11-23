pub mod config;
pub(crate) mod mods;

use std::collections::HashSet;

use model::modification::ModLocation;
use model::protein::ProteinDb;
use scoring_crate::RankScorer;

use crate::candidates::Candidate;
use crate::decoy::{build_search_db, whole_sequence_decoy_db, DecoyStrategy};
use crate::engine::PreparedSearch;
use crate::psm::TopNQueue;
use crate::refinement::config::RefineConfig;
pub use crate::refinement::mods::{mod_count_and_class, refinement_aa_set};
use crate::index::SearchIndex;
use crate::params::SearchParams;
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
        has_alternative_entries: crate::index::holds_alternative_entries(&db),
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
        has_alternative_entries: crate::index::holds_alternative_entries(&combined),
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
    files: &[crate::precursor::prepass::FileCalibration],
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
    if decoy_strategy == DecoyStrategy::Mutate {
        log::info!(
            target: "refine",
            "Pass-2 decoys are reversed: a refinement entry is one peptide, and a one-residue \
             mutant of a short peptide can keep its target sequence"
        );
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

    let fmt_mod = |m: &crate::refinement::config::RefineMod| {
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
        let queues = crate::precursor::prepass::run_by_file(
            &refine_prepared,
            all_spectra,
            &unident,
            &refine_params,
            files,
        );
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
            if crate::mimics::fits_better(
                pep,
                crate::mimics::mod_shift_da(pep),
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
#[path = "refinement/tests/refinement.rs"]
mod tests;
