//! Per-PSM provenance features: which kind of database entry the row's peptide
//! is specific to, how far it scored above the best Reference candidate of its
//! spectrum, and whether its precursor looks like an isotope mis-pick.

use model::EntryKind;
use search::candidate_gen::Candidate;
use search::isotope_mimic::{fits_better, isotope_gap_da, mod_shift_da};
use search::psm::{PsmMatch, TopNQueue};
use search::SearchIndex;

/// Fill the provenance features of every PSM in `queues`. `candidates` and
/// `idx` are the pool and index the PIN resolves through (after a `--refine`
/// merge, the merged ones).
pub fn populate_provenance_features(
    queues: &mut [TopNQueue],
    candidates: &[Candidate],
    idx: &SearchIndex,
) {
    for queue in queues.iter_mut() {
        queue.update_in_place(|psm| set_features(psm, candidates, idx));
    }
}

fn set_features(psm: &mut PsmMatch, candidates: &[Candidate], idx: &SearchIndex) {
    let (mut variant, mut isoform, mut chain) = (false, false, false);
    let mut variant_gap = 0.0;
    for &ci in &psm.candidate_idxs {
        let protein = &idx.db.proteins[candidates[ci as usize].protein_index];
        match protein.entry_kind() {
            EntryKind::Reference => {}
            EntryKind::Isoform => isoform = true,
            EntryKind::Chain => chain = true,
            EntryKind::Variant => {
                variant = true;
                variant_gap = isotope_gap_da(protein.origin().map_or(0.0, |o| o.variant_delta_da));
            }
        }
    }
    let rounded = psm.score.round();
    let peptide = &candidates[psm.primary_candidate_idx() as usize].peptide;
    // A flag is set when ANY protein of the peptide is of that kind, shared
    // with a Reference entry or not. Decoy generators do not keep a decoy
    // peptide shared with the decoy canonical (shuffle permutes each entry on
    // its own; a chain decoy is pseudo-reversed, its canonical fully
    // reversed), so "shared means reference" would set flags on decoys more
    // often than on targets, and Percolator would learn the label.
    let f = &mut psm.features;
    f.is_variant = variant as u32;
    f.is_isoform = isoform as u32;
    f.is_chain = chain as u32;
    let alternative = variant || isoform || chain;
    f.delta_ref_score = if alternative && f.ref_competitor == 1 {
        rounded - f.ref_best_raw
    } else {
        0.0
    };
    let shift = mod_shift_da(peptide) + if f.is_variant == 1 { variant_gap } else { 0.0 };
    f.isotope_mimic = fits_better(peptide, shift, psm.mass_error_ppm) as u32;
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::protein::{EntryKind, EntryOrigin, Protein, ProteinAnnotation, ProteinDb};
    use model::{AminoAcid, Peptide};
    use std::sync::Arc;

    fn protein(acc: &str, seq: &[u8], kind: EntryKind, delta: f64) -> Protein {
        let annotation = (kind != EntryKind::Reference).then(|| {
            Arc::new(
                ProteinAnnotation::new(Vec::new(), None).with_origin(Some(EntryOrigin {
                    kind,
                    family: Arc::from("P0"),
                    feature: Arc::from("F"),
                    variant_delta_da: delta,
                    knockoff_of: None,
                    source: None,
                    tier: model::PriorTier::Catalogue,
                })),
            )
        });
        Protein {
            accession: acc.into(),
            description: String::new(),
            sequence: seq.to_vec(),
            annotation,
        }
    }

    fn cand(seq: &[u8], protein_index: usize, is_decoy: bool) -> Candidate {
        Candidate {
            peptide: Peptide::new(
                seq.iter()
                    .map(|&r| AminoAcid::standard(r).unwrap())
                    .collect(),
                b'K',
                b'A',
            ),
            protein_index,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn psm(candidate_idxs: Vec<u32>, score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: candidate_idxs.into_iter().collect(),
            charge_used: 2,
            mass_error_ppm: 1.5,
            score,
            rank_score: score,
            edge_score: 0,
            activation_method: None,
            features: search::psm::PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
        }
    }

    fn queue(p: PsmMatch) -> TopNQueue {
        let mut q = TopNQueue::new(10);
        q.push(p);
        q
    }

    /// P0 Reference, P1 its isoform, P2 an N->D variant, and P2's decoy.
    /// Spectrum 0: PEPTIDEK shared by P0 and P1. Spectrum 1: the variant
    /// peptide PEPTIDDK (score 30). Spectrum 2: the decoy variant peptide.
    fn fixture() -> (SearchIndex, Vec<Candidate>, Vec<TopNQueue>) {
        let idx = SearchIndex {
            has_alternative_entries: true,
            db: ProteinDb {
                proteins: vec![
                    protein("sp|P0|A", b"PEPTIDEKAAK", EntryKind::Reference, 0.0),
                    protein("sp|P0-2|A", b"PEPTIDEKAAK", EntryKind::Isoform, 0.0),
                    protein(
                        "sp|P0_VAR_1|A",
                        b"PEPTIDDKAAK",
                        EntryKind::Variant,
                        0.984016,
                    ),
                    protein(
                        "XXX_sp|P0_VAR_1|A",
                        b"KAAKDDITPEP",
                        EntryKind::Variant,
                        0.984016,
                    ),
                ],
            },
            decoy_prefix: "XXX".into(),
            decoy_suffix: None,
        };
        let cands = vec![
            cand(b"PEPTIDEK", 0, false),
            cand(b"PEPTIDEK", 1, false),
            cand(b"PEPTIDDK", 2, false),
            cand(b"KDDITPEP", 3, true),
        ];
        let queues = vec![
            queue(psm(vec![0, 1], 40.0)),
            queue(psm(vec![2], 30.0)),
            queue(psm(vec![3], 12.0)),
        ];
        (idx, cands, queues)
    }

    fn features(q: &TopNQueue) -> &search::psm::PsmFeatures {
        &q.iter_psms().next().unwrap().features
    }

    #[test]
    fn shared_peptide_carries_every_kind_it_maps_to() {
        // A peptide in the canonical entry and its isoform is flagged as
        // isoform: decoy generators do not keep a decoy peptide shared with
        // the decoy canonical, so "shared means reference" would set the flag
        // on decoys more often than on targets.
        let (idx, cands, mut queues) = fixture();
        queues[0].update_in_place(|p| {
            p.features.ref_competitor = 1;
            p.features.ref_best_raw = 35.0;
        });
        populate_provenance_features(&mut queues, &cands, &idx);
        let f = features(&queues[0]);
        assert_eq!((f.is_variant, f.is_isoform, f.is_chain), (0, 1, 0));
        assert_eq!(f.delta_ref_score, 5.0);
    }

    #[test]
    fn reference_only_row_has_no_flags_and_no_gap() {
        let (idx, cands, _) = fixture();
        let mut queues = vec![queue(psm(vec![0], 40.0))];
        queues[0].update_in_place(|p| {
            p.features.ref_competitor = 1;
            p.features.ref_best_raw = 35.0;
        });
        populate_provenance_features(&mut queues, &cands, &idx);
        let f = features(&queues[0]);
        assert_eq!((f.is_variant, f.is_isoform, f.is_chain), (0, 0, 0));
        assert_eq!(f.delta_ref_score, 0.0);
    }

    #[test]
    fn variant_row_gets_its_flag_and_ref_gap() {
        let (idx, cands, mut queues) = fixture();
        queues[1].update_in_place(|p| {
            p.features.ref_competitor = 1;
            p.features.ref_best_raw = 25.0;
        });
        populate_provenance_features(&mut queues, &cands, &idx);
        let f = features(&queues[1]);
        assert_eq!((f.is_variant, f.is_isoform, f.is_chain), (1, 0, 0));
        assert_eq!(f.delta_ref_score, 5.0);
    }

    #[test]
    fn no_ref_competitor_means_no_gap() {
        let (idx, cands, mut queues) = fixture();
        populate_provenance_features(&mut queues, &cands, &idx);
        assert_eq!(features(&queues[1]).delta_ref_score, 0.0);
    }

    #[test]
    fn decoy_rows_carry_their_target_class() {
        let (idx, cands, mut queues) = fixture();
        populate_provenance_features(&mut queues, &cands, &idx);
        assert_eq!(features(&queues[2]).is_variant, 1);
    }

    #[test]
    fn variant_at_the_isotope_mis_pick_is_flagged() {
        // N->D (+0.984016): peptide − spectrum = −0.019339 Da is what a one-step
        // mis-pick of the unsubstituted peptide gives.
        let (idx, cands, mut queues) = fixture();
        let mass = cands[2].peptide.mass();
        queues[1].update_in_place(|p| p.mass_error_ppm = -0.019339 / mass * 1e6);
        populate_provenance_features(&mut queues, &cands, &idx);
        assert_eq!(features(&queues[1]).isotope_mimic, 1);
        // At zero error the substitution itself explains the precursor.
        let (idx, cands, mut queues) = fixture();
        queues[1].update_in_place(|p| p.mass_error_ppm = 0.0);
        populate_provenance_features(&mut queues, &cands, &idx);
        assert_eq!(features(&queues[1]).isotope_mimic, 0);
    }
}
