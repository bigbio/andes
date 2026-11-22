//! Unit tests of the per-charge split path and the chimeric isolation window of `match_engine`:
//! a chargeless spectrum searched once per charge, and the window centred on the isolation target.

use super::*;
use model::{AminoAcid, AminoAcidSetBuilder, Protein, ProteinDb};
use scoring_crate::scoring::fragment_ions::predict_by_ions;

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
    let scorer = crate::coisolation::tests::tiny_scorer();
    let prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();

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
        let pep = &accum.candidates[best.primary_candidate_idx() as usize].peptide;
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
