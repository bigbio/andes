//! Protein groups of the study and the feature records written from the finalized transfers.

use std::collections::HashMap;
use std::path::Path;

use output::features::{FeatureGroup, FeatureRecord, MbrColumns};
use output::library::StoredLibrary;
use output::mbr::EvidenceType;
use output::qpx::parse_peptidoform;
use output::{FeatureCounts, GroupAt, Inference, ProteinPosition};

use super::extract::split_key;
use crate::quant::finalize::MbrFeature;

pub(crate) struct ProteinGroups<'a> {
    pub inf: &'a Inference,
    pub max_q: f64,
}

#[derive(Default)]
pub(crate) struct GroupTally {
    pub counts: HashMap<usize, FeatureCounts>,
    pub no_group: usize,
    pub above_cut: usize,
    pub quant: output::ProteinQuant,
}

impl GroupTally {
    /// Adds a feature row to the group's default protein quantities in `unit` (the run, or the
    /// plex of an isobaric run); rows of shared peptides, rows without intensities and isobaric
    /// rows that failed the quantification gate (`quant_weight` 0) add nothing.
    pub(crate) fn quantify(&mut self, r: &FeatureRecord, unit: &str) -> Result<(), String> {
        let Some(g) = r.group.as_ref().filter(|g| g.unique) else {
            return Ok(());
        };
        let gated_out = r
            .additional_scores
            .iter()
            .any(|(name, w, _)| name == "quant_weight" && *w <= 0.0);
        if r.intensities.is_empty() || gated_out {
            return Ok(());
        }
        let labels: Vec<String> = r.intensities.iter().map(|(l, _)| l.clone()).collect();
        let values: Vec<f32> = r.intensities.iter().map(|(_, v)| *v).collect();
        self.quant.add(
            &g.anchor,
            unit,
            &r.run_file_name,
            &r.peptidoform,
            &labels,
            &values,
        )
    }

    /// Counts a row whose peptide is at `at` among the study's groups.
    pub(crate) fn count(&mut self, at: &GroupAt<'_>) {
        match *at {
            GroupAt::Group { group, unique } => {
                let c = self.counts.entry(group.id).or_default();
                c.total += 1;
                c.unique += unique as usize;
            }
            GroupAt::NoGroup => self.no_group += 1,
            GroupAt::AboveCut => self.above_cut += 1,
        }
    }

    pub(crate) fn assign(
        &mut self,
        pg: &ProteinGroups<'_>,
        sequence: &str,
    ) -> Option<(Vec<ProteinPosition>, FeatureGroup)> {
        let at = pg.inf.group_at(sequence, pg.max_q);
        self.count(&at);
        match at {
            GroupAt::Group { group, unique } => Some((
                group
                    .accessions
                    .iter()
                    .map(|a| (a.clone(), None, None, None, None))
                    .collect(),
                FeatureGroup {
                    anchor: group.lead.clone(),
                    unique,
                    global_q: group.q_value,
                },
            )),
            GroupAt::NoGroup | GroupAt::AboveCut => None,
        }
    }

    pub(crate) fn write_pg(
        &self,
        dir: &Path,
        pg: &ProteinGroups<'_>,
        runs: &[String],
    ) -> Result<(), String> {
        let path = dir.join("quantms.pg.parquet");
        output::write_pg_parquet(
            &path,
            pg.inf,
            runs,
            runs.len() == 1,
            pg.max_q,
            Some(&self.counts),
            (!self.quant.is_empty()).then_some(&self.quant),
        )
        .map_err(|e| format!("write {}: {e}", path.display()))?;
        let kept = pg.inf.groups.iter().filter(|g| g.q_value <= pg.max_q);
        let decoys = pg
            .inf
            .decoy_groups
            .iter()
            .filter(|g| g.q_value <= pg.max_q)
            .count();
        log::info!(
            target: "proteins",
            "wrote {} ({} target groups at protein q<={} over {} run(s), {} with features; {} decoy groups at the cut); features left out: {} of peptides in no protein group, {} of groups above the cut",
            path.display(),
            kept.count(),
            pg.max_q,
            runs.len(),
            self.counts.len(),
            decoys,
            self.no_group,
            self.above_cut
        );
        Ok(())
    }
}

pub(super) fn feature_records(
    stored: &StoredLibrary,
    features: &[MbrFeature],
    last_id: &mut i64,
    groups: Option<&ProteinGroups<'_>>,
    tally: &mut GroupTally,
) -> Vec<FeatureRecord> {
    features
        .iter()
        .filter(|f| {
            !matches!(
                f.evidence_type,
                EvidenceType::BelowThreshold
                    | EvidenceType::IdentifiedBelowThreshold
                    | EvidenceType::SharedFeature
            )
        })
        .filter_map(|f| {
            let m = &stored.library.meta[f.precursor as usize];
            let peptidoform = split_key(&m.key).to_string();
            let (sequence, modifications) = parse_peptidoform(&peptidoform);
            let (proteins, group): (Vec<ProteinPosition>, _) = match groups {
                Some(pg) => {
                    let (proteins, group) = tally.assign(pg, &sequence)?;
                    (proteins, Some(group))
                }
                None => (
                    stored.proteins[f.precursor as usize]
                        .split(';')
                        .filter(|a| !a.is_empty())
                        .map(|a| (a.to_string(), None, None, None, None))
                        .collect(),
                    None,
                ),
            };
            *last_id += 1;
            let run = stored.runs[f.run as usize].clone();
            Some(FeatureRecord {
                feature_id: *last_id,
                sequence,
                peptidoform,
                modifications,
                charge: m.charge as i16,
                pep: None,
                is_decoy: false,
                calculated_mz: m.mz as f32,
                observed_mz: f.mass_err_ppm.map_or(m.mz, |e| m.mz * (1.0 + e * 1e-6)) as f32,
                mass_error_ppm: f.mass_err_ppm.map(|e| e as f32),
                additional_scores: vec![("library_q".into(), f.library_q, false)],
                run_file_name: run.clone(),
                cv_params: Vec::new(),
                scan: Vec::new(),
                rt: f.apex_rt_s.map(|t| t as f32),
                missed_cleavages: None,
                intensities: f
                    .quantity
                    .map(|q| vec![("LFQ".to_string(), q as f32)])
                    .unwrap_or_default(),
                additional_intensities: Vec::new(),
                proteins,
                id_run_file_name: (f.evidence_type == EvidenceType::Identified).then_some(run),
                rt_start: None,
                rt_stop: None,
                psm_ids: Vec::new(),
                mbr: Some(MbrColumns {
                    evidence_type: f.evidence_type.as_str().to_string(),
                    transfer_q: f.transfer_q,
                    quant_q: f.quant_q,
                    quantity_se: f.quantity_se,
                    claimed_by: f
                        .claim
                        .claimant
                        .map(|c| stored.library.meta[c as usize].key.clone()),
                    claim_co_apex: f.evidence_type.is_transfer().then_some(f.claim.co_apex),
                }),
                group,
            })
        })
        .collect()
}
