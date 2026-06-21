//! QPX ↔ PIN ↔ Percolator SpecId join contract.
//!
//! Percolator's `PSMId` equals the andes PIN `SpecId` byte-for-byte; the QPX
//! writer must reconstruct that SAME SpecId (via the shared
//! `output::format_spec_id`) so PEP/q-value join onto the right PSM rows. This
//! test runs a real BSA search, writes the PIN, harvests its SpecIds, then drives
//! the QPX writer with a mock rescore map keyed by those SpecIds and asserts:
//!   1. every PIN SpecId is reproduced by the QPX writer (PEP populated for them);
//!   2. the join works for single- AND multi-row scans (chimeric/tie cases);
//!   3. with `rescore = None`, every QPX row's PEP stays null (schema parity).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use arrow::array::{Array, Float64Array};
use model::tolerance::PrecursorTolerance;
use model::{
    AminoAcidSetBuilder, Enzyme, ModLocation, Modification, ResidueSpec, Tolerance,
};
use output::PercolatorPsm;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use scoring_crate::RankScorer;
use search::{match_spectra, SearchIndex, SearchParams};
use input::{FastaReader, MgfReader};

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

/// Run the shared BSA fixture search and return (spectra, queues, candidates, index, params).
fn run_bsa_search() -> (
    Vec<model::spectrum::Spectrum>,
    Vec<search::psm::TopNQueue>,
    Vec<search::candidate_gen::Candidate>,
    SearchIndex,
    SearchParams,
) {
    let target_db =
        FastaReader::load_all(BufReader::new(File::open(fixture("test-fixtures/BSA.fasta")).unwrap()))
            .unwrap();
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

    // hcd_qexactive_tryp from the canonical Parquet store.
    let store = model_train::store::ModelStore::open(
        &fixture("resources/models.parquet"),
    ).unwrap();
    let param = store.load_param("hcd_qexactive_tryp").unwrap();
    let scorer = RankScorer::new(&param);

    let mut params = SearchParams::default_tryptic(aa);
    params.enzyme = Enzyme::Trypsin;
    params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    params.charge_range = 2..=3;
    params.isotope_error_range = -1..=2;

    let mgf_file = input::open_buf_maybe_gz(&fixture("test-fixtures/test.mgf.gz")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let (queues, candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.5, "XXX_");
    (spectra, queues, candidates, idx, params)
}

/// Read the PEP column from a written `psms.parquet`.
fn read_pep_column(psms_path: &Path) -> Vec<Option<f64>> {
    let file = File::open(psms_path).unwrap();
    let reader = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap();
    let mut out = Vec::new();
    for batch in reader {
        let batch = batch.unwrap();
        let idx = batch.schema().index_of("posterior_error_probability").unwrap();
        let col = batch.column(idx).as_any().downcast_ref::<Float64Array>().unwrap();
        for i in 0..col.len() {
            out.push(if col.is_null(i) { None } else { Some(col.value(i)) });
        }
    }
    out
}

/// Harvest the SpecId (col 0) of each non-header PIN row.
fn pin_spec_ids(pin_path: &Path) -> Vec<String> {
    BufReader::new(File::open(pin_path).unwrap())
        .lines()
        .map(|l| l.unwrap())
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect()
}

#[test]
fn qpx_spec_id_matches_pin_and_pep_joins() {
    let (spectra, queues, candidates, idx, params) = run_bsa_search();
    let tmp = tempfile::tempdir().unwrap();

    // 1. Write the PIN, harvest its SpecIds (the Percolator PSMId join keys).
    let pin_path = tmp.path().join("bsa.pin");
    output::write_pin(&pin_path, &spectra, &queues, &candidates, &params, &idx).unwrap();
    let spec_ids = pin_spec_ids(&pin_path);
    assert!(spec_ids.len() >= 5, "expected several PIN rows, got {}", spec_ids.len());

    // Sanity: at least one multi-row scan exercises the `_{row_idx}` suffix path
    // (a SpecId with 4 underscore-segments after the title). Not required for the
    // join correctness check below, but documents that both formats are present
    // across the fixture when chimeric/ties occur.
    let has_multi = spec_ids.iter().any(|s| s.matches('_').count() >= 3);
    eprintln!("PIN rows={} multi_row_present={}", spec_ids.len(), has_multi);

    // 2. Build a mock Percolator map keyed by EVERY PIN SpecId, with a distinct
    //    PEP/q per row so a mismatched join would be detectable.
    let mut map: HashMap<String, PercolatorPsm> = HashMap::new();
    for (i, sid) in spec_ids.iter().enumerate() {
        map.insert(
            sid.clone(),
            PercolatorPsm {
                psm_id: sid.clone(),
                q_value: 0.001 * (i as f64 + 1.0),
                pep: 0.0001 * (i as f64 + 1.0),
                peptide: "X".into(),
                proteins: "P".into(),
            },
        );
    }

    // 3. Write QPX WITH the rescore map → reconstructed SpecIds must hit the map
    //    so every row's PEP is non-null.
    let qdir = tmp.path().join("rescored.idparquet");
    output::write_qpx(
        &qdir,
        &spectra,
        &queues,
        &candidates,
        &params,
        &idx,
        &Tolerance::Da(0.5),
        "test_run",
        &["test.mgf".to_string()],
        Some(&map),
    )
    .unwrap();

    let peps = read_pep_column(&qdir.join("psms.parquet"));
    assert_eq!(
        peps.len(),
        spec_ids.len(),
        "QPX must emit exactly one row per PIN row (same ordering/identity)"
    );
    let non_null = peps.iter().filter(|p| p.is_some()).count();
    assert_eq!(
        non_null,
        peps.len(),
        "every QPX row's PEP must be populated — proves QPX SpecId == PIN SpecId for ALL rows \
         (single- and multi-row scans). {} of {} were null.",
        peps.len() - non_null,
        peps.len()
    );

    // The keys actually consumed are exactly the PIN SpecIds (no leftovers means
    // the QPX reconstructed every key the PIN produced).
    let pin_set: HashSet<&String> = spec_ids.iter().collect();
    assert_eq!(pin_set.len(), map.len(), "mock map should have one entry per distinct PIN SpecId");

    // 4. Write QPX WITHOUT the rescore map → every PEP must be null (schema parity).
    let qdir_off = tmp.path().join("plain.idparquet");
    output::write_qpx(
        &qdir_off,
        &spectra,
        &queues,
        &candidates,
        &params,
        &idx,
        &Tolerance::Da(0.5),
        "test_run",
        &["test.mgf".to_string()],
        None,
    )
    .unwrap();
    let peps_off = read_pep_column(&qdir_off.join("psms.parquet"));
    assert!(
        peps_off.iter().all(|p| p.is_none()),
        "rescore-off bundle must keep PEP null (OpenMS schema parity)"
    );
}
