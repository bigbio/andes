//! Lookups of the candidates of a `PreparedSearch`: by predicate, by record, by global index.

use model::peptide::Peptide;
use rayon::prelude::*;

use super::dedup::*;
use super::{enumeration_pass, PreparedSearch};
use crate::candidates::{base_record_key, expand_base_record, Candidate};

impl<'a> PreparedSearch<'a> {
    pub fn candidates_where(&self, keep: impl Fn(&[u8], u8) -> bool + Sync) -> Vec<Candidate> {
        let db = self.idx;
        let params = self.params;
        let records = self.records_where(keep);
        let mut tagged: Vec<(u32, u32, u32, u16, Candidate)> = records
            .par_iter()
            .flat_map_iter(|(rec, mult)| {
                let mult = *mult;
                expand_base_record(db, params, rec)
                    .into_iter()
                    .flat_map(move |cand| {
                        (0..mult).map(move |k| {
                            (
                                cand.protein_index as u32,
                                enumeration_pass(&cand, k),
                                cand.start_offset_in_protein as u32,
                                cand.peptide.residues.len() as u16,
                                cand.clone(),
                            )
                        })
                    })
            })
            .collect();
        tagged.par_sort_by(|a, b| (a.0, a.1, a.2, a.3).cmp(&(b.0, b.1, b.2, b.3)));
        let mut cands: Vec<Candidate> = tagged.into_par_iter().map(|t| t.4).collect();
        relabel_collision_decoys_with(&mut cands, &self.target_bare_seqs);
        if self.candidate_filter.is_some() {
            cands.retain(|c| self.keeps(c));
        }
        cands
    }

    pub fn records_where(
        &self,
        keep: impl Fn(&[u8], u8) -> bool + Sync,
    ) -> Vec<(crate::candidates::records::IndexRecord, u32)> {
        let mi = &self.mmap_index;
        let db = self.idx;
        let mut kept: Vec<crate::candidates::records::IndexRecord> = (0..mi.len())
            .into_par_iter()
            .filter_map(|i| {
                let rec = mi.record_at(i);
                let seq = &db.db.proteins.get(rec.protein_index as usize)?.sequence;
                let start = rec.start_offset as usize;
                let end = start + rec.length as usize;
                let span = seq.get(start..end)?;
                let post = seq.get(end).copied().unwrap_or(b'-');
                keep(span, post).then_some(rec)
            })
            .collect();
        kept.par_sort_unstable_by_key(base_record_key);
        let mut records: Vec<(crate::candidates::records::IndexRecord, u32)> = Vec::new();
        for rec in kept {
            match records.last_mut() {
                Some((last, n)) if base_record_key(last) == base_record_key(&rec) => *n += 1,
                _ => records.push((rec, 1)),
            }
        }
        records
    }

    pub fn candidate_peptide(&self, global_idx: u32) -> Peptide {
        let accum = self.mmap_accum.lock().unwrap();
        match accum.candidates.get(global_idx as usize) {
            Some(c) => c.peptide.clone(),
            None => self.candidates[global_idx as usize].peptide.clone(),
        }
    }

    pub fn candidate_is_decoy(&self, global_idx: u32) -> bool {
        let accum = self.mmap_accum.lock().unwrap();
        match accum.candidates.get(global_idx as usize) {
            Some(c) => c.is_decoy,
            None => self.candidates[global_idx as usize].is_decoy,
        }
    }
}
