//! `.pin` schema test for Andes's GF-free output.
//!
//! NOTE: Andes removed the generating function entirely, so its `.pin` schema
//! INTENTIONALLY diverges from the historical GF-heavy Percolator schema: the
//! GF-derived columns (`DeNovoScore`, `lnSpecEValue`, `lnEValue`,
//! `lnDeltaSpecEValue`) are no longer emitted, and `RawScore` is the sole score
//! column. This test asserts the GF-free schema rather than comparing against
//! a legacy reference fixture:
//!   - `RawScore` is present, the GF columns are absent.
//!   - The additive feature columns (`EdgeScore`, `PrecursorIsotopeKL`,
//!     `PrecursorSNR`, `DeltaRawScore`) are present, between `matchedIonRatio`
//!     and `Peptide`.
//!   - Every data row has at least as many columns as the header.

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

/// Columns that the generating function used to emit and that GF-free Andes
/// must NOT emit anymore.
const GF_COLUMNS: [&str; 4] = ["DeNovoScore", "lnSpecEValue", "lnEValue", "lnDeltaSpecEValue"];

/// Additive feature columns Andes emits between matchedIonRatio and Peptide.
const ADDITIVE_COLUMNS: [&str; 4] =
    ["EdgeScore", "PrecursorIsotopeKL", "PrecursorSNR", "DeltaRawScore"];

#[test]
fn rust_pin_header_is_gf_free_schema() {
    let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa.clone());
    params.enzyme = Enzyme::Trypsin;
    params.charge_range = 2..=3;

    let empty_target = ProteinDb::default();
    let empty_idx = SearchIndex::from_target_db(&empty_target, "XXX_");
    let tmp_dir = tempfile::tempdir().expect("tempdir");
    let rust_pin_path = tmp_dir.path().join("empty.pin");
    output::write_pin(&rust_pin_path, &[], &[], &[], &params, &empty_idx).expect("write_pin");

    let header = first_line(&rust_pin_path);
    let cols: Vec<&str> = header.split('\t').collect();

    // RawScore present; GF columns absent.
    assert!(cols.contains(&"RawScore"), "RawScore column must be present:\n{header}");
    for gf in GF_COLUMNS {
        assert!(
            !cols.contains(&gf),
            "GF-derived column {gf} must NOT appear in the GF-free schema:\n{header}"
        );
    }

    // Additive feature columns present, between matchedIonRatio and Peptide.
    let matched_ratio_pos = cols
        .iter()
        .position(|c| *c == "matchedIonRatio")
        .expect("matchedIonRatio missing");
    let peptide_pos = cols.iter().position(|c| *c == "Peptide").expect("Peptide missing");
    for name in ADDITIVE_COLUMNS {
        let pos = cols
            .iter()
            .position(|c| *c == name)
            .unwrap_or_else(|| panic!("Rust .pin header is missing the additive feature column {name}"));
        assert!(
            matched_ratio_pos < pos && pos < peptide_pos,
            "additive column {name} must sit between matchedIonRatio and Peptide"
        );
    }

    // Proteins must be the very last column (after Peptide).
    let proteins_pos = cols
        .iter()
        .rposition(|c| *c == "Proteins")
        .expect("Proteins column missing");
    assert_eq!(
        proteins_pos,
        cols.len() - 1,
        "Proteins must be the last column in the PIN header; PIN Proteins is rest-of-line and \
         cannot be followed by additional columns without corrupting Percolator's protein parsing"
    );
    assert!(
        peptide_pos < proteins_pos,
        "column order must be ... Peptide ... Proteins (last)"
    );
}

#[test]
fn rust_pin_rows_have_at_least_header_column_count() {
    // Run a real search and assert each data row has at least as many columns
    // as the header (trailing Proteins columns may add more). Schema-only; no
    // value comparison against any external reference fixture.

    let target_db = FastaReader::load_all(BufReader::new(File::open(fixture("test-fixtures/BSA.fasta")).unwrap())).unwrap();
    let idx = SearchIndex::from_target_db(&target_db, "XXX_");

    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
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
    };
    let aa = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap();

    let param_path = fixture("test-fixtures/HCD_QExactive_Tryp.param");
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

    let tmp_dir = tempfile::tempdir().expect("tempdir");
    let rust_pin_path = tmp_dir.path().join("bsa.pin");
    output::write_pin(&rust_pin_path, &spectra, &queues, &candidates, &params, &idx).expect("write_pin");

    let rust_lines: Vec<_> = BufReader::new(File::open(&rust_pin_path).unwrap())
        .lines()
        .collect::<Result<_, _>>()
        .unwrap();

    assert!(rust_lines.len() >= 6, "Rust pin should have at least 5 data rows");

    let header_cols = rust_lines[0].split('\t').count();
    let mut row_count = 0;
    for (i, rust_line) in rust_lines.iter().enumerate().skip(1) {
        let rust_row_cols = rust_line.split('\t').count();
        assert!(
            rust_row_cols >= header_cols,
            "Rust row {i} has {rust_row_cols} cols, expected >= {header_cols}"
        );
        row_count += 1;
    }
    assert!(row_count >= 5, "checked {row_count} rows, expected >= 5");
}
