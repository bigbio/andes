//! Unit tests of the per-charge split path and the chimeric isolation window of `engine`:
//! a chargeless spectrum searched once per charge, and the window centred on the isolation target.

use super::*;
use model::{AminoAcid, AminoAcidSetBuilder, Protein, ProteinDb};
use scoring_crate::psm::ions::predict_by_ions;

fn peptide(residues: &[u8]) -> Peptide {
    let aas: Vec<AminoAcid> = residues
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    Peptide::new(aas, b'K', b'K')
}

fn spectrum(precursor_mz: f64, charge: Option<i32>, peptides: &[&Peptide]) -> Spectrum {
    let mut peaks: Vec<(f64, f32)> = peptides
        .iter()
        .flat_map(|p| predict_by_ions(p, 1..=1))
        .map(|ion| (ion.mz, 100.0_f32))
        .collect();
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    Spectrum {
        title: "chargeless".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: charge,
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    }
}

#[test]
fn the_chimeric_window_is_centred_on_the_isolation_target() {
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let mut spec = spectrum(500.3, Some(2), &[]);
    (spec.isolation_target_mz, spec.isolation_lower_offset) = (Some(500.0), Some(1.0));
    spec.isolation_upper_offset = Some(1.0);
    let (lo, hi, _) = coisolation_window(&spec, &params);
    assert_eq!((lo, hi), (499.0, 501.0));
}

#[test]
fn chargeless_spectrum_is_searched_per_charge() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"GKDAFLGSFLYEYSRKSAFLWEYSPLLTNNGVRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.charge_range = 2..=5;
    params.fragment_index_top_k = 1;
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let prepared = PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();

    let truth = peptide(b"DAFLGSFLYEYSR");
    let bait = peptide(b"SAFLWEYSPLLTNNGVR");
    let mz3 = (truth.mass() + 3.0 * PROTON) / 3.0;
    let spectra = [
        spectrum(mz3, Some(3), &[&truth, &bait]),
        spectrum(mz3, None, &[&truth, &bait]),
    ];
    let mut queues = prepared.run_chunk(&spectra, 7);
    assert_eq!(queues.len(), 2);
    for (i, q) in queues.iter_mut().enumerate() {
        let psms = q.drain_into_vec();
        let best = psms
            .iter()
            .max_by(|a, b| a.rank_score.total_cmp(&b.rank_score))
            .unwrap_or_else(|| panic!("spectrum {i} came back empty"));
        assert_eq!(best.spectrum_idx, 7 + i);
        assert_eq!(best.charge_used, 3, "spectrum {i}");
        let accum = prepared.mmap_accum.lock().unwrap();
        let pep = accum
            .peptide(best.primary_candidate_idx() as usize)
            .expect("materialised");
        assert_eq!(pep.residues, truth.residues, "spectrum {i}");
    }

    let mut top1 = params.clone();
    top1.top_n_psms_per_spectrum = 1;
    let one = PreparedSearch::prepare_in_memory(&idx, &top1, &scorer, 0.02, "XXX").unwrap();
    let mut queues = one.run_chunk(&spectra[1..], 0);
    let psms = queues[0].drain_into_vec();
    let best = psms.iter().map(|p| p.rank_score).fold(f32::MIN, f32::max);
    assert!(
        psms.iter().all(|p| p.rank_score == best),
        "the merged charge queues keep the top PSM and its ties only"
    );
}

#[test]
fn top_only_search_keeps_the_psms_of_the_full_search() {
    use crate::params::ScoreMode;
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"GKDAFLGSFLYEYSRKSAFLWEYSPLLTNNGVRKLVNELTEFAKTCVADESHAGCEK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let truth = peptide(b"DAFLGSFLYEYSR");
    let bait = peptide(b"SAFLWEYSPLLTNNGVR");
    let other = peptide(b"LVNELTEFAK");
    let mz3 = (truth.mass() + 3.0 * PROTON) / 3.0;
    let mz2 = (other.mass() + 2.0 * PROTON) / 2.0;
    let spectra = [
        spectrum(mz3, Some(3), &[&truth, &bait]),
        spectrum(mz3, None, &[&truth, &bait]),
        spectrum(mz2, Some(2), &[&other]),
        spectrum(mz2, None, &[&other, &truth]),
    ];
    type Key = (u8, u32, u32, i32, (String, bool), Vec<(String, bool)>);
    let candidate = |s: &PreparedSearch, i: u32| {
        (s.candidate_peptide(i).to_string(), s.candidate_is_decoy(i))
    };
    let key = |s: &PreparedSearch, q: TopNQueue| -> Vec<Key> {
        let mut v: Vec<Key> = q
            .into_sorted_vec()
            .iter()
            .map(|p| {
                let mut all: Vec<_> = p.candidate_idxs.iter().map(|&i| candidate(s, i)).collect();
                all.sort_unstable();
                (
                    p.charge_used,
                    p.score.to_bits(),
                    p.rank_score.to_bits(),
                    p.edge_score,
                    candidate(s, p.primary_candidate_idx()),
                    all,
                )
            })
            .collect();
        v.sort_unstable();
        v
    };
    let materialized = |s: &PreparedSearch| {
        let mut v: Vec<_> = (0..s.materialized_len() as u32)
            .map(|i| candidate(s, i))
            .collect();
        v.sort_unstable();
        v
    };
    for (mode, top_k) in [
        (ScoreMode::Rank, 0),
        (ScoreMode::Rank, 1),
        (ScoreMode::Strong, 1),
    ] {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.charge_range = 2..=4;
        params.fragment_index_top_k = top_k;
        params.score_mode = mode;
        let scorer = crate::engine::coisolation::tests::tiny_scorer();
        let full = PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
        let top = PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
        let a: Vec<_> = full
            .run_chunk(&spectra, 3)
            .into_iter()
            .map(|q| key(&full, q))
            .collect();
        let b: Vec<_> = top
            .run_chunk_for_top(&spectra, 3, &params)
            .into_iter()
            .map(|q| key(&top, q))
            .collect();
        assert!(a.iter().any(|q| !q.is_empty()), "{mode:?} {top_k}");
        assert_eq!(a, b, "{mode:?} top_k {top_k}");
        assert_eq!(
            materialized(&full),
            materialized(&top),
            "{mode:?} top_k {top_k}"
        );
    }
}
