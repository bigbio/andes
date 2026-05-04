//! End-to-end Phase 4e: BSA.fasta + test.mgf → top-N PSMs.
//! First full test on real local data.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use engine::{
    match_spectra, AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec,
    SearchIndex, SearchParams,
};
use input::{FastaReader, MgfReader};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("src/test/resources")
        .join(name)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {name}: {e}"))
}

fn aa_set() -> engine::AminoAcidSet {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
    };
    AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

#[test]
fn bsa_test_mgf_produces_some_matches() {
    let target = FastaReader::load_all(BufReader::new(File::open(fixture("BSA.fasta")).unwrap())).unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();
    assert!(!spectra.is_empty(), "test.mgf must contain at least one spectrum");

    let queues = match_spectra(&spectra, &idx, &params, "XXX");
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
