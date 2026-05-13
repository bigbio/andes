//! End-to-end Phase 4e: BSA.fasta + test.mgf → top-N PSMs.
//! First full test on real local data.

mod common;
use common::*;

use std::fs::File;
use std::io::BufReader;

use search::{match_spectra, SearchIndex, SearchParams};
use input::{FastaReader, MgfReader};

#[test]
fn bsa_test_mgf_produces_some_matches() {
    let target = FastaReader::load_all(BufReader::new(File::open(fixture("src/test/resources/BSA.fasta")).unwrap())).unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("src/test/resources/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();
    assert!(!spectra.is_empty(), "test.mgf must contain at least one spectrum");

    let (queues, _candidates) = match_spectra(&spectra, &idx, &params, &rank_scorer(), 0.05, "XXX");
    assert_eq!(queues.len(), spectra.len());

    // At least one spectrum should have a match (BSA is a known target).
    let total_matches: usize = queues.iter().map(|q| q.len()).sum();
    assert!(total_matches > 0,
        "expected at least one PSM across {} spectra, got 0", spectra.len());

    // For non-empty queues, top match's mass error should be within 20 ppm.
    for q in queues {
        if q.is_empty() { continue; }
        let top = q.into_sorted_vec();
        let best = &top[0];
        assert!(best.mass_error_ppm.abs() < 20.0,
            "best PSM mass_error_ppm {} > 20.0", best.mass_error_ppm);
    }
}
