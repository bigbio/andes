//! BSA + Tryp_Pig_Bov candidate-enumeration sanity tests.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use engine::{
    enumerate_candidates, AminoAcidSetBuilder, ModLocation, Modification,
    ResidueSpec, SearchIndex, SearchParams,
};
use input::FastaReader;

fn fasta(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("src/test/resources")
        .join(name)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {name}: {e}"))
}

fn aa_set_with_carbamidomethyl_oxidation() -> engine::AminoAcidSet {
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
fn bsa_generates_reasonable_candidate_count() {
    let target = FastaReader::load_all(BufReader::new(File::open(fasta("BSA.fasta")).unwrap())).unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set_with_carbamidomethyl_oxidation());

    let candidates: Vec<_> = enumerate_candidates(&idx, &params, "XXX").collect();

    assert!(candidates.len() > 50, "got {} candidates, expected > 50", candidates.len());
    assert!(candidates.len() < 50_000, "got {} candidates, expected < 50,000", candidates.len());

    for c in &candidates {
        assert!(c.peptide.length() >= 6, "peptide too short: {}", c.peptide.length());
        assert!(c.peptide.length() <= 40, "peptide too long: {}", c.peptide.length());
        assert!(c.protein_index < 2, "BSA has only 2 proteins (target+decoy)");
    }
    assert!(candidates.iter().any(|c| !c.is_decoy));
    assert!(candidates.iter().any(|c| c.is_decoy));
}

#[test]
fn tryp_pig_bov_generates_more_candidates_than_bsa() {
    let bsa_target = FastaReader::load_all(BufReader::new(File::open(fasta("BSA.fasta")).unwrap())).unwrap();
    let bsa_idx = SearchIndex::from_target_db(&bsa_target, "XXX");
    let params = SearchParams::default_tryptic(aa_set_with_carbamidomethyl_oxidation());
    let bsa_count = enumerate_candidates(&bsa_idx, &params, "XXX").count();

    let tpb_target = FastaReader::load_all(BufReader::new(File::open(fasta("Tryp_Pig_Bov.fasta")).unwrap())).unwrap();
    let tpb_idx = SearchIndex::from_target_db(&tpb_target, "XXX");
    let tpb_count = enumerate_candidates(&tpb_idx, &params, "XXX").count();

    assert!(tpb_count > bsa_count,
        "Tryp_Pig_Bov ({} candidates) should generate more than BSA ({})",
        tpb_count, bsa_count);
}
