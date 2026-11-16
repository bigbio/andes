use std::io::BufReader;
use std::path::Path;

use input::MgfReader;
use model::{AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec};
use scoring_crate::RankScorer;
use search::SearchParams;

use model_train::labeled::bootstrap_labels;

fn fixture(rel: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn bsa_aa_set() -> model::AminoAcidSet {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

fn load_hcd_scorer() -> RankScorer {
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/models"
    ));
    let store = model_train::store::ModelStore::open(bundled).expect("open bundled model store");
    let param = store
        .load_param("hcd_qexactive_tryp")
        .expect("load hcd_qexactive_tryp from store");
    RankScorer::new(&param)
}

fn load_bsa_spectra() -> Vec<model::Spectrum> {
    let path = fixture("test-fixtures/test.mgf.gz");
    let file = input::open_buf_maybe_gz(&path).expect("open test.mgf");
    MgfReader::new(BufReader::new(file))
        .filter_map(|r| r.ok())
        .collect()
}

#[test]
fn bootstrap_labels_bsa_smoke() {
    let cache = tempfile::tempdir().expect("cache dir");
    let spectra = load_bsa_spectra();
    assert!(
        !spectra.is_empty(),
        "test.mgf must contain at least one spectrum"
    );

    let scorer = load_hcd_scorer();
    let db_path = fixture("test-fixtures/BSA.fasta");
    let params = SearchParams::default_tryptic(bsa_aa_set());
    let train_fdr = 0.5_f64;

    let result = bootstrap_labels(
        &spectra,
        &db_path,
        &scorer,
        &params,
        train_fdr,
        cache.path(),
    );
    assert!(
        result.is_ok(),
        "bootstrap_labels must succeed: {:?}",
        result.err()
    );

    let matches = result.unwrap();

    for m in &matches {
        assert!(
            m.confidence <= train_fdr,
            "match confidence {} > train_fdr {}",
            m.confidence,
            train_fdr
        );
        assert!(
            m.spectrum_index < spectra.len(),
            "spectrum_index {} out of range (spectra.len()={})",
            m.spectrum_index,
            spectra.len()
        );
    }

    let result2 = bootstrap_labels(
        &spectra,
        &db_path,
        &scorer,
        &params,
        train_fdr,
        cache.path(),
    )
    .expect("second bootstrap_labels run must succeed");
    assert_eq!(
        matches.len(),
        result2.len(),
        "bootstrap_labels must be deterministic: first run returned {}, second returned {}",
        matches.len(),
        result2.len()
    );

    eprintln!(
        "bootstrap_labels BSA @ train_fdr={}: {} confident target PSMs",
        train_fdr,
        matches.len()
    );
}

#[test]
fn bootstrap_labels_bsa_strict_fdr() {
    let cache = tempfile::tempdir().expect("cache dir");
    let spectra = load_bsa_spectra();
    let scorer = load_hcd_scorer();
    let db_path = fixture("test-fixtures/BSA.fasta");
    let params = SearchParams::default_tryptic(bsa_aa_set());
    let train_fdr = 0.05_f64;

    let matches = bootstrap_labels(
        &spectra,
        &db_path,
        &scorer,
        &params,
        train_fdr,
        cache.path(),
    )
    .expect("bootstrap_labels must succeed at strict FDR");

    for m in &matches {
        assert!(
            m.confidence <= train_fdr,
            "match confidence {} > train_fdr {}",
            m.confidence,
            train_fdr
        );
    }

    eprintln!(
        "bootstrap_labels BSA @ train_fdr={}: {} confident target PSMs",
        train_fdr,
        matches.len()
    );
}

#[test]
fn bootstrap_labels_bsa_lenient_ge_strict() {
    let cache = tempfile::tempdir().expect("cache dir");
    let spectra = load_bsa_spectra();
    let scorer = load_hcd_scorer();
    let db_path = fixture("test-fixtures/BSA.fasta");
    let params = SearchParams::default_tryptic(bsa_aa_set());

    let lenient = bootstrap_labels(&spectra, &db_path, &scorer, &params, 0.5, cache.path())
        .expect("lenient bootstrap");
    let strict = bootstrap_labels(&spectra, &db_path, &scorer, &params, 0.05, cache.path())
        .expect("strict bootstrap");

    assert!(
        lenient.len() >= strict.len(),
        "lenient FDR=0.5 should give >= matches as strict FDR=0.05: \
         lenient={}, strict={}",
        lenient.len(),
        strict.len()
    );
}
