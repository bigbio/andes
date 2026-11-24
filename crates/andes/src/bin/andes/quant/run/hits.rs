//! Choosing the identifications to quantify (rank-1 PSMs at the FDR cut) and the shared fields
//! of a quant record: peptidoform, modifications, proteins and the quant id.

use std::collections::{HashMap, HashSet};

use model::mass::{ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use output::features::FeatureRecord;
use output::{LfqRow, ModRecord, PercolatorPsm, QuantId};
use search::candidates::Candidate;
use search::index::SearchIndex;
use search::psm::{PsmMatch, TopNQueue};

use super::collect::{QuantCollector, RunCapture};
use super::settings::QuantSettings;

pub(super) struct Hit<'a> {
    pub(super) spec_idx: usize,
    pub(super) co_isolated: bool,
    pub(super) run_idx: usize,
    pub(super) psm: PsmMatch,
    pub(super) cand: &'a Candidate,
    pub(super) spec_id: String,
    pub(super) q_value: Option<f64>,
    pub(super) pep: Option<f64>,
    pub(super) psm_index: i64,
}

pub(super) fn select_hits<'a>(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &'a [Candidate],
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    omitted: &HashSet<String>,
    cutoffs: Option<&output::Cutoffs>,
) -> Vec<Hit<'a>> {
    select_rank1(
        collector,
        spectra,
        queues,
        candidates,
        rescore,
        Some(omitted),
        cutoffs,
        false,
    )
}

pub(super) fn select_rank1_all<'a>(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &'a [Candidate],
) -> Vec<Hit<'a>> {
    select_rank1(
        collector, spectra, queues, candidates, None, None, None, true,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn select_rank1<'a>(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &'a [Candidate],
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    omitted: Option<&HashSet<String>>,
    cutoffs: Option<&output::Cutoffs>,
    include_decoys: bool,
) -> Vec<Hit<'a>> {
    let s = &collector.settings;
    let mut hits = Vec::new();
    let mut psm_index: i64 = 0;
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let spec = &spectra[spec_idx];
        let sorted = queue.clone().into_rank_sorted_vec();
        let rows = output::emitted_rows(spec, &sorted, candidates);
        let multi_row = rows.len() > 1;
        let mut best_rank: Vec<(Option<u64>, u32)> = Vec::new();
        for (rank, psm) in &rows {
            let key = psm.precursor_mz_override.map(f64::to_bits);
            match best_rank.iter_mut().find(|(k, _)| *k == key) {
                Some((_, r)) => *r = (*r).min(*rank),
                None => best_rank.push((key, *rank)),
            }
        }
        let scan = spec.scan.unwrap_or(0);
        let spec_id_base = if spec.title.is_empty() {
            format!("scan={scan}")
        } else {
            spec.title.clone()
        };
        for (hit, (rank, psm)) in rows.into_iter().enumerate() {
            let psm = psm.into_owned();
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            let spec_id = output::format_spec_id(&spec_id_base, scan, rank, hit, multi_row);
            let rescored = rescore.map(|m| m.get(&spec_id));
            if omitted.is_some_and(|o| o.contains(&spec_id))
                || (!cand.is_decoy && matches!(rescored, Some(None)))
            {
                continue;
            }
            let this_index = psm_index;
            psm_index += 1;
            let key = psm.precursor_mz_override.map(f64::to_bits);
            if !best_rank.iter().any(|&(k, r)| k == key && r == rank) {
                continue;
            }
            if cand.is_decoy && !include_decoys {
                continue;
            }
            let (q_value, pep) = match rescored {
                Some(Some(r)) => {
                    let cutoff = cutoffs.map_or(s.psm_fdr, |c| c.or(&spec_id, s.psm_fdr));
                    if r.q_value > cutoff || s.pep_cap.is_some_and(|c| r.pep > c) {
                        continue;
                    }
                    (Some(r.q_value), Some(r.pep))
                }
                Some(None) => continue,
                None => (None, None),
            };
            let Some(run_idx) = collector.run_of(spec_idx) else {
                continue;
            };
            hits.push(Hit {
                spec_idx,
                co_isolated: key.is_some(),
                run_idx,
                psm,
                cand,
                spec_id,
                q_value,
                pep,
                psm_index: this_index,
            });
        }
    }
    hits
}

pub(super) fn peptidoform_of(cand: &Candidate) -> String {
    output::qpx::peptidoform_string(cand)
}

pub(super) fn modifications_of(cand: &Candidate) -> Vec<ModRecord> {
    cand.peptide
        .residues
        .iter()
        .enumerate()
        .flat_map(|(i, aa)| ModRecord::of_residue(i, cand.peptide.residues.len(), aa))
        .collect()
}

pub(super) type ProteinPos = output::ProteinPosition;

pub(super) fn proteins_of(
    psm: &PsmMatch,
    candidates: &[Candidate],
    index: &SearchIndex,
) -> Vec<ProteinPos> {
    let mut out: Vec<ProteinPos> = Vec::new();
    for &ci in &psm.candidate_idxs {
        let c = &candidates[ci as usize];
        let acc = match index.protein_at(c.protein_index) {
            Some(p) => p.accession.clone(),
            None => format!("PROT_{}", c.protein_index),
        };
        if out.iter().any(|(a, ..)| *a == acc) {
            continue;
        }
        let start = c.start_offset_in_protein as i32 + 1;
        let end = start + c.peptide.length() as i32 - 1;
        out.push((
            acc,
            Some(start),
            Some(end),
            Some((c.peptide.pre as char).to_string()),
            Some((c.peptide.post as char).to_string()),
        ));
    }
    out
}

pub(super) fn sole_protein(hit: &Hit<'_>, candidates: &[Candidate]) -> Option<u32> {
    let mut proteins = hit
        .psm
        .candidate_idxs
        .iter()
        .map(|&ci| candidates[ci as usize].protein_index);
    let first = proteins.next()?;
    proteins.all(|p| p == first).then_some(first as u32)
}

pub(super) fn quant_id(
    hit: &Hit<'_>,
    run: &RunCapture,
    spec: &Spectrum,
    candidates: &[Candidate],
    index: &SearchIndex,
) -> QuantId {
    let z = hit.psm.charge_used.max(1) as f64;
    let precursor_mz = hit.psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
    let observed_mz = precursor_mz - ISOTOPE * (hit.psm.isotope_offset as f64) / z;
    let calculated_mz = hit.cand.peptide.mass() / z + PROTON;
    QuantId {
        file_idx: hit.run_idx,
        run: run.stem.clone(),
        spec_id: hit.spec_id.clone(),
        scan: spec.scan.unwrap_or(0),
        rt_seconds: spec.rt_seconds,
        sequence: output::qpx::bare_sequence(hit.cand),
        peptidoform: peptidoform_of(hit.cand),
        modifications: modifications_of(hit.cand),
        charge: hit.psm.charge_used,
        proteins: proteins_of(&hit.psm, candidates, index)
            .into_iter()
            .map(|(a, ..)| a)
            .collect(),
        is_decoy: hit.cand.is_decoy,
        calculated_mz,
        observed_mz,
        q_value: hit.q_value,
        pep: hit.pep,
        missed_cleavages: None,
        psm_index: hit.psm_index,
    }
}

pub(super) fn base_record(
    id: &QuantId,
    proteins: Vec<ProteinPos>,
    feature_id: i64,
) -> FeatureRecord {
    let mass_error_ppm = if id.calculated_mz > 0.0 {
        Some(((id.observed_mz - id.calculated_mz) / id.calculated_mz * 1e6) as f32)
    } else {
        None
    };
    let mut additional_scores: Vec<(String, f64, bool)> = Vec::new();
    if let Some(q) = id.q_value {
        additional_scores.push(("psm_q_value".into(), q, false));
    }
    FeatureRecord {
        feature_id,
        sequence: id.sequence.clone(),
        peptidoform: id.peptidoform.clone(),
        modifications: id.modifications.clone(),
        charge: id.charge as i16,
        pep: id.pep,
        is_decoy: id.is_decoy,
        calculated_mz: id.calculated_mz as f32,
        observed_mz: id.observed_mz as f32,
        mass_error_ppm,
        additional_scores,
        run_file_name: id.run.clone(),
        cv_params: Vec::new(),
        scan: if id.scan > 0 {
            vec![id.scan]
        } else {
            Vec::new()
        },
        rt: id.rt_seconds.map(|t| t as f32),
        missed_cleavages: id.missed_cleavages,
        intensities: Vec::new(),
        additional_intensities: Vec::new(),
        proteins,
        id_run_file_name: Some(id.run.clone()),
        rt_start: None,
        rt_stop: None,
        psm_ids: vec![id.psm_index],
        mbr: None,
        group: None,
    }
}

pub(super) fn is_confident(row: &LfqRow, s: &QuantSettings) -> bool {
    row.feature_q_value.is_some_and(|q| q <= s.feature_fdr)
        && row.feature.cosine >= s.lfq_params.min_cosine
}
