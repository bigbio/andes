//! QPX ↔ PIN ↔ Percolator SpecId join contract.
//!
//! Percolator's `PSMId` equals the andes PIN `SpecId` byte-for-byte; the QPX
//! writer must reconstruct that SAME SpecId (via the shared
//! `output::format_spec_id`) so PEP/q-value join onto the right PSM rows. This
//! test runs a real BSA search, writes the PIN, harvests its SpecIds, then drives
//! the QPX writer with a mock rescore map keyed by those SpecIds and asserts:
//!   1. every PIN SpecId is reproduced by the QPX writer (PEP populated for them);
//!   2. the join works for single- AND multi-row scans (chimeric/tie cases);
//!   3. with `rescore = None`, every QPX row's PEP stays null (schema parity);
//!   4. a result in which Percolator's target-decoy competition (`-Y`) kept one
//!      PSM per scan/ExpMass group joins, and the competed-away targets are left
//!      out of the bundle.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use arrow::array::{Array, Float64Array};
use input::{FastaReader, MgfReader};
use model::tolerance::PrecursorTolerance;
use model::{AminoAcidSetBuilder, Enzyme, ModLocation, Modification, ResidueSpec, Tolerance};
use output::PercolatorPsm;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use scoring_crate::RankScorer;
use search::{match_spectra_out_of_core, SearchIndex, SearchParams};

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
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
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
    };
    let aa = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap();

    // hcd_qexactive_tryp from the canonical Parquet store.
    let store = model_train::store::ModelStore::open(&fixture("resources/models")).unwrap();
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

    let (queues, candidates) =
        match_spectra_out_of_core(&spectra, &idx, &params, &scorer, 0.5, "XXX_").expect("search");
    (spectra, queues, candidates, idx, params)
}

/// Read the PEP column from a written `psms.parquet`.
fn read_pep_column(psms_path: &Path) -> Vec<Option<f64>> {
    let file = File::open(psms_path).unwrap();
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .unwrap()
        .build()
        .unwrap();
    let mut out = Vec::new();
    for batch in reader {
        let batch = batch.unwrap();
        let idx = batch
            .schema()
            .index_of("posterior_error_probability")
            .unwrap();
        let col = batch
            .column(idx)
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        for i in 0..col.len() {
            out.push(if col.is_null(i) {
                None
            } else {
                Some(col.value(i))
            });
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
    assert!(
        spec_ids.len() >= 5,
        "expected several PIN rows, got {}",
        spec_ids.len()
    );

    // Sanity: at least one multi-row scan exercises the `_{row_idx}` suffix path
    // (a SpecId with 4 underscore-segments after the title). Not required for the
    // join correctness check below, but documents that both formats are present
    // across the fixture when chimeric/ties occur.
    let has_multi = spec_ids.iter().any(|s| s.matches('_').count() >= 3);
    eprintln!(
        "PIN rows={} multi_row_present={}",
        spec_ids.len(),
        has_multi
    );

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
                score: None,
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
        &[spectra.len()],
        Some(&map),
        None,
        &std::collections::HashSet::new(),
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
    assert_eq!(
        pin_set.len(),
        map.len(),
        "mock map should have one entry per distinct PIN SpecId"
    );

    // Excluded rows (knockoff PSMs) are left out of the bundle, and only they.
    let excluded: HashSet<String> = [spec_ids[0].clone()].into_iter().collect();
    let qdir_ex = tmp.path().join("excluded.idparquet");
    output::write_qpx(
        &qdir_ex,
        &spectra,
        &queues,
        &candidates,
        &params,
        &idx,
        &Tolerance::Da(0.5),
        "test_run",
        &["test.mgf".to_string()],
        &[spectra.len()],
        Some(&map),
        None,
        &excluded,
    )
    .unwrap();
    assert_eq!(
        read_pep_column(&qdir_ex.join("psms.parquet")).len(),
        spec_ids.len() - 1
    );

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
        &[spectra.len()],
        None,
        None,
        &std::collections::HashSet::new(),
    )
    .unwrap();
    let peps_off = read_pep_column(&qdir_off.join("psms.parquet"));
    assert!(
        peps_off.iter().all(|p| p.is_none()),
        "rescore-off bundle must keep PEP null (QPX schema parity)"
    );
}

#[test]
fn batched_psms_preserve_global_indices_order_and_rescore_join() {
    use arrow::array::Int32Array;
    let (spectra, queues, candidates, idx, params) = run_bsa_search();
    let seed = queues.iter().position(|q| !q.is_empty()).unwrap();
    let repeated_spectra: Vec<_> = (0..1025)
        .map(|i| {
            let mut s = spectra[seed].clone();
            s.scan = Some(i + 1);
            s.title = format!("scan={}", i + 1);
            s
        })
        .collect();
    let repeated_queues = vec![queues[seed].clone(); repeated_spectra.len()];
    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("batch.pin");
    output::write_pin(
        &pin,
        &repeated_spectra,
        &repeated_queues,
        &candidates,
        &params,
        &idx,
    )
    .unwrap();
    let ids = pin_spec_ids(&pin);
    let map = ids
        .iter()
        .enumerate()
        .map(|(i, id)| {
            (
                id.clone(),
                PercolatorPsm {
                    psm_id: id.clone(),
                    q_value: 0.01,
                    pep: i as f64 / ids.len() as f64,
                    score: None,
                    peptide: "X".into(),
                    proteins: "P".into(),
                },
            )
        })
        .collect();
    let out = tmp.path().join("batch.idparquet");
    output::write_qpx(
        &out,
        &repeated_spectra,
        &repeated_queues,
        &candidates,
        &params,
        &idx,
        &Tolerance::Da(0.5),
        "test",
        &["test.mgf".into()],
        &[repeated_spectra.len()],
        Some(&map),
        None,
        &std::collections::HashSet::new(),
    )
    .unwrap();
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(File::open(out.join("psms.parquet")).unwrap())
            .unwrap();
    assert_eq!(builder.metadata().num_row_groups(), 2);
    let mut row = 0usize;
    for batch in builder.build().unwrap() {
        let batch = batch.unwrap();
        let indices = batch
            .column_by_name("peptide_identification_index")
            .unwrap()
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let peps = batch
            .column_by_name("posterior_error_probability")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        for i in 0..batch.num_rows() {
            assert_eq!(indices.value(i), row as i32);
            assert_eq!(peps.value(i), row as f64 / ids.len() as f64);
            row += 1;
        }
    }
    assert_eq!(row, ids.len());
}

#[test]
fn multi_input_rows_carry_their_own_reference_file() {
    use arrow::array::{Int32Array, StringArray};
    let (spectra, queues, candidates, idx, params) = run_bsa_search();
    let seed = queues.iter().position(|q| !q.is_empty()).unwrap();
    // 1025 spectra span two 1024-spectrum batches; the first input ends inside
    // the first batch.
    let repeated_spectra: Vec<_> = (0..1025)
        .map(|i| {
            let mut s = spectra[seed].clone();
            s.scan = Some(i + 1);
            s
        })
        .collect();
    let repeated_queues = vec![queues[seed].clone(); repeated_spectra.len()];
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("two.idparquet");
    let paths = ["a.mzML".to_string(), "b.mzML".to_string()];
    output::write_qpx(
        &out,
        &repeated_spectra,
        &repeated_queues,
        &candidates,
        &params,
        &idx,
        &Tolerance::Da(0.5),
        "test",
        &paths,
        &[1000, 1025],
        None,
        None,
        &std::collections::HashSet::new(),
    )
    .unwrap();
    let reader =
        ParquetRecordBatchReaderBuilder::try_new(File::open(out.join("psms.parquet")).unwrap())
            .unwrap()
            .build()
            .unwrap();
    let (mut a, mut b) = (0usize, 0usize);
    for batch in reader {
        let batch = batch.unwrap();
        let scans = batch
            .column_by_name("scan")
            .unwrap()
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let files = batch
            .column_by_name("reference_file_name")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        for i in 0..batch.num_rows() {
            let expected = if scans.value(i) <= 1000 {
                a += 1;
                "a.mzML"
            } else {
                b += 1;
                "b.mzML"
            };
            assert_eq!(files.value(i), expected, "scan {}", scans.value(i));
        }
    }
    assert!(a > 0 && b > 0, "rows from both inputs: a={a} b={b}");

    // A partition that does not cover the spectra is rejected.
    let err = output::write_qpx(
        &tmp.path().join("bad.idparquet"),
        &repeated_spectra,
        &repeated_queues,
        &candidates,
        &params,
        &idx,
        &Tolerance::Da(0.5),
        "test",
        &paths,
        &[1000, 1024],
        None,
        None,
        &std::collections::HashSet::new(),
    )
    .unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
}

#[test]
fn dropped_peaks_write_null_peak_lists() {
    let (spectra, queues, candidates, idx, params) = run_bsa_search();
    let null_counts = |spectra: &[model::spectrum::Spectrum], name: &str| {
        let tmp = tempfile::tempdir().unwrap();
        let out = tmp.path().join(format!("{name}.idparquet"));
        output::write_qpx(
            &out,
            spectra,
            &queues,
            &candidates,
            &params,
            &idx,
            &Tolerance::Da(0.5),
            "test",
            &["test.mgf".to_string()],
            &[spectra.len()],
            None,
            None,
            &std::collections::HashSet::new(),
        )
        .unwrap();
        let reader =
            ParquetRecordBatchReaderBuilder::try_new(File::open(out.join("psms.parquet")).unwrap())
                .unwrap()
                .build()
                .unwrap();
        let (mut rows, mut mz_null, mut int_null) = (0, 0, 0);
        for batch in reader {
            let batch = batch.unwrap();
            rows += batch.num_rows();
            mz_null += batch.column_by_name("mz_array").unwrap().null_count();
            int_null += batch
                .column_by_name("intensity_array")
                .unwrap()
                .null_count();
        }
        (rows, mz_null, int_null)
    };
    let (rows, mz_null, int_null) = null_counts(&spectra, "kept");
    assert!(rows > 0);
    assert_eq!((mz_null, int_null), (0, 0), "peaks present → lists present");

    let dropped: Vec<_> = spectra
        .iter()
        .cloned()
        .map(|mut s| {
            s.peaks = Vec::new();
            s
        })
        .collect();
    let (rows, mz_null, int_null) = null_counts(&dropped, "dropped");
    assert_eq!(
        (mz_null, int_null),
        (rows, rows),
        "peaks dropped → null lists"
    );
}

/// Percolator `-Y` result for a PIN: per (ScanNr, ExpMass) group only the
/// highest-scoring row survives, in the target file (a target) or the decoy file
/// (a decoy). Returns (target map, decoy SpecIds, number of target rows removed).
fn competition_result(pin_path: &Path) -> (HashMap<String, PercolatorPsm>, HashSet<String>, usize) {
    let text = std::fs::read_to_string(pin_path).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|c| *c == n).unwrap();
    let (id, label, scan, exp, score) = (
        col("SpecId"),
        col("Label"),
        col("ScanNr"),
        col("ExpMass"),
        col("RankScore"),
    );
    let mut best: HashMap<(String, String), (f64, String, bool)> = HashMap::new();
    let mut total_targets = 0usize;
    for line in lines.filter(|l| !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        total_targets += usize::from(f[label] != "-1");
        let s: f64 = f[score].parse().unwrap();
        let entry = (s, f[id].to_string(), f[label] == "-1");
        best.entry((f[scan].to_string(), f[exp].to_string()))
            .and_modify(|b| {
                if s > b.0 {
                    *b = entry.clone();
                }
            })
            .or_insert(entry);
    }
    let mut targets = HashMap::new();
    let mut decoys = HashSet::new();
    for (_, sid, is_decoy) in best.values() {
        if *is_decoy {
            decoys.insert(sid.clone());
        } else {
            targets.insert(
                sid.clone(),
                PercolatorPsm {
                    psm_id: sid.clone(),
                    q_value: 0.01,
                    pep: 0.01,
                    score: None,
                    peptide: "X".into(),
                    proteins: "P".into(),
                },
            );
        }
    }
    let removed = total_targets - targets.len();
    (targets, decoys, removed)
}

/// With top-N candidates per scan, Percolator's target-decoy competition (`-Y`)
/// removes every target but the winner of its scan/ExpMass group. The join must
/// accept those rows as competed away when the decoy results are given (and
/// leave them out of the bundle), and must keep failing without them.
#[test]
fn qpx_join_accepts_targets_removed_by_target_decoy_competition() {
    let (spectra, mut queues, candidates, idx, params) = run_bsa_search();
    // Give every scan a lower-ranked second candidate at the same charge and
    // precursor, as the top-N chimeric search does: one competition group with
    // two rows. The alternative is a different target peptide, so the writers'
    // positional-isomer collapse keeps both rows.
    for q in queues.iter_mut() {
        if let Some(top) = q.peek_top().cloned() {
            let residues = |c: &search::candidate_gen::Candidate| -> Vec<u8> {
                c.peptide.residues.iter().map(|aa| aa.residue).collect()
            };
            let top_residues = residues(&candidates[top.primary_candidate_idx() as usize]);
            let Some(other) = candidates
                .iter()
                .position(|c| !c.is_decoy && residues(c) != top_residues)
            else {
                continue;
            };
            let mut alt = top;
            alt.candidate_idxs = vec![other as u32];
            alt.score -= 1.0;
            alt.rank_score -= 1.0;
            q.force_push(alt);
        }
    }
    let tmp = tempfile::tempdir().unwrap();
    let pin_path = tmp.path().join("bsa.pin");
    output::write_pin(&pin_path, &spectra, &queues, &candidates, &params, &idx).unwrap();
    let (targets, decoys, removed) = competition_result(&pin_path);
    let pin_rows = pin_spec_ids(&pin_path).len();
    assert!(
        removed > 0,
        "the fixture must have groups with several rows"
    );

    let write = |dir: &Path, decoys: Option<&HashSet<String>>| {
        output::write_qpx(
            dir,
            &spectra,
            &queues,
            &candidates,
            &params,
            &idx,
            &Tolerance::Da(0.5),
            "test_run",
            &["test.mgf".to_string()],
            &[spectra.len()],
            Some(&targets),
            decoys,
            &std::collections::HashSet::new(),
        )
    };
    let err = write(&tmp.path().join("strict.idparquet"), None).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);

    let qdir = tmp.path().join("competed.idparquet");
    write(&qdir, Some(&decoys)).unwrap();
    let peps = read_pep_column(&qdir.join("psms.parquet"));
    // The competed-away targets are left out; every target row left is a
    // winner with its PEP, and decoy rows keep a null PEP as before.
    assert_eq!(peps.len(), pin_rows - removed);
    assert_eq!(peps.iter().filter(|p| p.is_some()).count(), targets.len());
}
