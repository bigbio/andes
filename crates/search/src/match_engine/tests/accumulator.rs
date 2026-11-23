use super::*;
use model::{AminoAcidSetBuilder, Protein, ProteinDb};

#[test]
fn sync_resets_the_identity_map() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = crate::coisolation::tests::tiny_scorer();
    let mut prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
    let cands = prepared.out_of_core_candidates_where(|_, _| true);
    assert!(cands.len() >= 2);
    {
        let mut accum = prepared.mmap_accum.lock().unwrap();
        assert_eq!(accum.intern(&cands[0]), 0);
        assert_eq!(accum.intern(&cands[1]), 1);
    }
    prepared.sync_materialized_candidates();
    assert_eq!(prepared.candidates.len(), 2);
    let mut accum = prepared.mmap_accum.lock().unwrap();
    assert!(accum.index_of.is_empty());
    assert_eq!(accum.intern(&cands[1]), 0);
    assert_eq!(accum.candidates.len(), 1);
}

#[test]
fn dia_ms1_gate_reads_neighbouring_survey_scans() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = crate::coisolation::tests::tiny_scorer();
    let mut prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, 0.02, "XXX").unwrap();
    let scan = |rt: f64, mz: f64| model::scan::Ms1Scan {
        rt,
        peaks: vec![(mz, 1.0)],
        ..Default::default()
    };
    prepared.set_dia_ms1(vec![
        scan(30.0, 700.0),
        scan(10.0, 500.0),
        scan(20.0, 600.0),
    ]);
    assert_eq!(prepared.dia_ms1_scans(None), None);
    let r = prepared.dia_ms1_scans(Some(9.0)).unwrap();
    assert_eq!(r, 0..2);
    assert!(prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 4e-6), 5.0));
    assert!(!prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 6e-6), 5.0));
    assert!(prepared.dia_ms1_has(r.clone(), 600.0 * (1.0 + 6e-6), 25.0));
    assert!(!prepared.dia_ms1_has(r, 700.0, 25.0));
    assert_eq!(prepared.dia_ms1_scans(Some(26.0)).unwrap(), 1..3);

    let e = prepared
        .dia_ms1_error_ppm(Some(20.0), 600.0 * (1.0 - 7e-6), 25.0)
        .unwrap();
    assert!((e - 7.0).abs() < 1e-3, "{e}");
    assert_eq!(prepared.dia_ms1_error_ppm(Some(20.0), 650.0, 25.0), None);
    prepared.shift_dia_ms1_ppm(7.0);
    let e = prepared
        .dia_ms1_error_ppm(Some(20.0), 600.0 * (1.0 - 7e-6), 25.0)
        .unwrap();
    assert!(e.abs() < 1e-3, "{e}");
}

#[test]
fn windowed_top_k_keeps_the_k_most_intense_per_window_in_mz_order() {
    let peaks = [
        (101.0, 1.0),
        (150.0, 5.0),
        (180.0, 3.0),
        (199.0, 4.0),
        (250.0, 2.0),
    ];
    assert_eq!(
        windowed_top_k(&peaks, 100.0, 2),
        vec![(150.0, 5.0), (199.0, 4.0), (250.0, 2.0)]
    );
    assert_eq!(windowed_top_k(&peaks, 100.0, 9), peaks.to_vec());
}
