//! `--gf-free` mode: the emitted PIN must OMIT the generating-function-derived
//! columns (`lnSpecEValue`, `DeNovoScore`, `lnEValue`, `lnDeltaSpecEValue`),
//! while the default (GF) mode keeps emitting them. RawScore is present in both.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use model::tolerance::PrecursorTolerance;
use model::{AminoAcidSetBuilder, Enzyme, ModLocation, Modification, ResidueSpec, Tolerance};
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

/// Build a BSA/test.mgf search and write a PIN with the given `gf_free` flag,
/// returning all PIN lines (header first).
fn run_and_get_lines(gf_free: bool) -> Vec<String> {
    let target_db = FastaReader::load_all(BufReader::new(
        File::open(fixture("test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
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

    let param_path = fixture("test-fixtures/HCD_QExactive_Tryp.param");
    let param = Param::load_from_file(&param_path).unwrap();
    let scorer = RankScorer::new(&param);

    let mut params = SearchParams::default_tryptic(aa.clone());
    params.enzyme = Enzyme::Trypsin;
    params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    params.charge_range = 2..=3;
    params.isotope_error_range = -1..=2;
    params.gf_free = gf_free;

    let mgf_file = File::open(fixture("test-fixtures/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let (queues, candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.5, "XXX_");

    let tmp_dir = tempfile::tempdir().expect("tempdir");
    let pin_path = tmp_dir.path().join("bsa.pin");
    output::write_pin(&pin_path, &spectra, &queues, &candidates, &params, &idx).expect("write_pin");

    BufReader::new(File::open(&pin_path).unwrap())
        .lines()
        .collect::<Result<_, _>>()
        .unwrap()
}

#[test]
fn gf_free_pin_omits_gf_columns_default_keeps_them() {
    const GF_COLS: [&str; 4] =
        ["lnSpecEValue", "lnEValue", "DeNovoScore", "lnDeltaSpecEValue"];

    // Default (GF) mode: every GF-derived column present, plus RawScore.
    let gf_lines = run_and_get_lines(false);
    let gf_header = &gf_lines[0];
    let gf_cols: Vec<&str> = gf_header.split('\t').collect();
    assert!(gf_cols.contains(&"RawScore"), "default header must contain RawScore");
    for c in GF_COLS {
        assert!(
            gf_cols.contains(&c),
            "default (GF) header must contain {c}; header: {gf_header}"
        );
    }

    // GF-free mode: GF-derived columns gone, RawScore still present.
    let free_lines = run_and_get_lines(true);
    let free_header = &free_lines[0];
    let free_cols: Vec<&str> = free_header.split('\t').collect();
    assert!(
        free_cols.contains(&"RawScore"),
        "gf-free header must still contain RawScore; header: {free_header}"
    );
    for c in GF_COLS {
        assert!(
            !free_cols.contains(&c),
            "gf-free header must NOT contain {c}; header: {free_header}"
        );
    }

    // Sanity: dropping exactly the 4 GF columns is the only difference.
    assert_eq!(
        free_cols.len() + GF_COLS.len(),
        gf_cols.len(),
        "gf-free header should be the default header minus exactly the {} GF columns\n\
         default: {gf_header}\n\
         gf-free: {free_header}",
        GF_COLS.len(),
    );

    // Schema validity: each GF-free data row must have at least as many columns
    // as the GF-free header (trailing Proteins columns may add more). Proves
    // the conditional header and row builders stay in sync.
    let header_cols = free_cols.len();
    let mut data_rows = 0;
    for row in free_lines.iter().skip(1) {
        let n = row.split('\t').count();
        assert!(
            n >= header_cols,
            "gf-free data row has {n} cols, expected >= {header_cols}: {row}"
        );
        data_rows += 1;
    }
    assert!(data_rows >= 1, "expected at least one gf-free data row");
}
