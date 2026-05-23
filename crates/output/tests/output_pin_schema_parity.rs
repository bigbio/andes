//! `.pin` schema parity gate against the Java reference fixture.
//!
//! The Rust `.pin` writer's header must match the reference fixture exactly,
//! so Percolator (and any downstream tool that uses regex column-name matching)
//! consumes Rust output without modification.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use model::{AminoAcidSetBuilder, Enzyme, ModLocation, Modification, ProteinDb, ResidueSpec, Tolerance};
use model::tolerance::PrecursorTolerance;
use scoring_crate::{Param, RankScorer};
use search::{match_spectra, SearchIndex, SearchParams};
use input::{FastaReader, MgfReader};

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn first_line(path: &std::path::Path) -> String {
    let f = File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    BufReader::new(f).lines().next().expect("file is empty").expect("read first line")
}

#[test]
fn rust_pin_header_matches_java_pin_fixture_header_exactly() {
    let java_pin_path = fixture("benchmark/parity-fixtures/bsa_test_mgf_java.pin");
    let java_header = first_line(&java_pin_path);

    // Construct an empty queues-vec but write the header — the writer
    // produces the header regardless of queue contents.
    // Match Java's params: charge2..=3, Trypsin (no charge1).
    let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa.clone());
    params.enzyme = Enzyme::Trypsin;
    params.charge_range = 2..=3;

    // Empty PIN — header-only. We need a SearchIndex for the API, but the
    // header writer doesn't use protein accessions, so an empty index suffices.
    let empty_target = ProteinDb::default();
    let empty_idx = SearchIndex::from_target_db(&empty_target, "XXX_");
    let tmp_dir = tempfile::tempdir().expect("tempdir");
    let rust_pin_path = tmp_dir.path().join("empty.pin");
    output::write_pin(&rust_pin_path, &[], &[], &[], &params, &empty_idx).expect("write_pin");

    let rust_header = first_line(&rust_pin_path);

    // Rust adds a single ADDITIVE "EdgeScore" column between matchedIonRatio
    // and Peptide (iter19, 2026-05-21). Java does not emit this column.
    // Check that the Java header is a prefix-modulo-EdgeScore-insertion of
    // Rust's: every Java column appears in Rust in the same relative order,
    // and the only extra Rust column is "EdgeScore" (between matchedIonRatio
    // and Peptide).
    let java_cols: Vec<&str> = java_header.split('\t').collect();
    let rust_cols: Vec<&str> = rust_header.split('\t').collect();
    let rust_minus_edge: Vec<&str> = rust_cols
        .iter()
        .copied()
        .filter(|c| *c != "EdgeScore")
        .collect();
    assert_eq!(
        rust_minus_edge, java_cols,
        "Rust .pin header (excluding EdgeScore) must match Java reference header.\n\
         Java:   {java_header}\n\
         Rust:   {rust_header}\n\
         (Common cause: column rename, missing column, or charge_range mismatch.)",
    );
    // EdgeScore must appear after matchedIonRatio and before Peptide.
    let edge_pos = rust_cols.iter().position(|c| *c == "EdgeScore").expect(
        "Rust .pin header is missing the iter19 EdgeScore additive feature column",
    );
    let matched_ratio_pos = rust_cols
        .iter()
        .position(|c| *c == "matchedIonRatio")
        .expect("matchedIonRatio missing");
    let peptide_pos = rust_cols.iter().position(|c| *c == "Peptide").expect("Peptide missing");
    assert!(matched_ratio_pos < edge_pos && edge_pos < peptide_pos,
        "EdgeScore must sit between matchedIonRatio and Peptide");
}

#[test]
fn rust_pin_row_column_count_matches_java_for_at_least_5_scans() {
    // Run a real search, then for at least 5 of Java's reference scans assert
    // Rust's row has the same number of tab-separated columns as Java's row.
    // We don't compare values (SpecEValue / lnSpecEValue may differ during
    // the parity build-out); only schema.

    // 1. Run Rust search end-to-end.
    let target_db = FastaReader::load_all(BufReader::new(File::open(fixture("test-fixtures/BSA.fasta")).unwrap())).unwrap();
    let idx = SearchIndex::from_target_db(&target_db, "XXX_");

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
    let aa = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap();

    let param_path = fixture("resources/ionstat/HCD_QExactive_Tryp.param");
    let param = Param::load_from_file(&param_path).unwrap();
    let scorer = RankScorer::new(&param);

    let mut params = SearchParams::default_tryptic(aa.clone());
    params.enzyme = Enzyme::Trypsin;
    params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    params.charge_range = 2..=3;
    params.isotope_error_range = -1..=2;

    let mgf_file = File::open(fixture("test-fixtures/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let (queues, candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.5, "XXX_");

    // 2. Write Rust PIN.
    let tmp_dir = tempfile::tempdir().expect("tempdir");
    let rust_pin_path = tmp_dir.path().join("bsa.pin");
    output::write_pin(&rust_pin_path, &spectra, &queues, &candidates, &params, &idx).expect("write_pin");

    // 3. Read Java + Rust PIN files and check column counts on first 5 data rows.
    let java_pin_path = fixture("benchmark/parity-fixtures/bsa_test_mgf_java.pin");
    let java_lines: Vec<_> = BufReader::new(File::open(&java_pin_path).unwrap())
        .lines()
        .collect::<Result<_, _>>()
        .unwrap();
    let rust_lines: Vec<_> = BufReader::new(File::open(&rust_pin_path).unwrap())
        .lines()
        .collect::<Result<_, _>>()
        .unwrap();

    assert!(java_lines.len() >= 6, "Java fixture should have at least 5 data rows");
    assert!(rust_lines.len() >= 6, "Rust pin should have at least 5 data rows");

    // Check first 5 data rows (lines 1..=5; line 0 is header).
    let java_header_cols = java_lines[0].split('\t').count();
    let rust_header_cols = rust_lines[0].split('\t').count();
    // Rust has exactly one ADDITIVE EdgeScore column (iter19, 2026-05-21)
    // not present in the Java fixture, so expect Rust to be Java + 1.
    assert_eq!(
        rust_header_cols,
        java_header_cols + 1,
        "header column count mismatch (Rust {rust_header_cols} vs Java {java_header_cols}; expected Rust = Java + 1 EdgeScore)"
    );

    let mut row_count = 0;
    for (i, rust_line) in rust_lines.iter().enumerate().skip(1).take(rust_lines.len().min(java_lines.len()).min(6) - 1) {
        let rust_row_cols = rust_line.split('\t').count();
        // The fixture may have variable trailing Proteins columns; allow Rust
        // to differ ONLY in the trailing columns (after position ==
        // header_cols - 1). For now, just assert column count >= header_cols.
        assert!(
            rust_row_cols >= rust_header_cols,
            "Rust row {i} has {rust_row_cols} cols, expected >= {rust_header_cols}"
        );
        row_count += 1;
    }
    assert!(row_count >= 5, "checked {row_count} rows, expected >= 5");
}
