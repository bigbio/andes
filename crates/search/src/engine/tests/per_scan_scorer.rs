use super::*;
use model::{AminoAcidSetBuilder, Protein, ProteinDb};
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

fn spectrum(activation: Option<ActivationMethod>) -> Spectrum {
    Spectrum {
        precursor_mz: 500.0,
        peaks: vec![(200.0, 10.0)],
        activation_method: activation,
        ..Spectrum::default()
    }
}

fn prepared_with<'a>(
    idx: &'a SearchIndex,
    params: &'a SearchParams,
    scorer: &'a RankScorer,
) -> PreparedSearch<'a> {
    PreparedSearch::prepare_in_memory(idx, params, scorer, 0.02, "XXX").unwrap()
}

#[test]
fn scans_route_by_their_own_activation() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let params =
        SearchParams::default_tryptic(AminoAcidSetBuilder::new_standard().build().unwrap());
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let calls = Arc::new(AtomicUsize::new(0));
    let seen = Arc::clone(&calls);
    let prepared = prepared_with(&idx, &params, &scorer).with_per_scan_scorers(
        ActivationMethod::HCD,
        Box::new(move |a| {
            seen.fetch_add(1, Ordering::Relaxed);
            (a == ActivationMethod::EThcD).then(crate::engine::coisolation::tests::tiny_scorer)
        }),
    );
    let chunk = [
        spectrum(Some(ActivationMethod::HCD)),
        spectrum(Some(ActivationMethod::EThcD)),
        spectrum(None),
        spectrum(Some(ActivationMethod::ETD)),
        spectrum(Some(ActivationMethod::EThcD)),
    ];
    let set = prepared.scorer_set(&chunk);
    let own = |s: &Spectrum| std::ptr::eq(set.for_spec(s), &scorer);
    assert!(own(&chunk[0]));
    assert!(!own(&chunk[1]), "EThcD scan has its own scorer");
    assert!(own(&chunk[2]));
    assert!(
        own(&chunk[3]),
        "a declined activation keeps the run's scorer"
    );
    assert!(!own(&chunk[4]));
    let _ = prepared.scorer_set(&chunk);
    assert_eq!(calls.load(Ordering::Relaxed), 2);
    let mut routes = prepared.per_scan_routes();
    routes.sort_by_key(|r| r.0 as u8);
    assert_eq!(
        routes.iter().map(|r| (r.0, r.1, r.2)).collect::<Vec<_>>(),
        vec![
            (ActivationMethod::ETD, false, 2),
            (ActivationMethod::EThcD, true, 4)
        ]
    );
}

#[test]
fn no_routing_keeps_the_run_scorer() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&db, "XXX");
    let params =
        SearchParams::default_tryptic(AminoAcidSetBuilder::new_standard().build().unwrap());
    let scorer = crate::engine::coisolation::tests::tiny_scorer();
    let prepared = prepared_with(&idx, &params, &scorer);
    let chunk = [spectrum(Some(ActivationMethod::EThcD))];
    let set = prepared.scorer_set(&chunk);
    assert!(std::ptr::eq(set.for_spec(&chunk[0]), &scorer));
    assert!(prepared.per_scan_routes().is_empty());
}
