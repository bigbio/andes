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
    let aa = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap();

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

    let cache_dir = tempfile::tempdir().expect("cache dir");
    let (queues, candidates) = match_spectra_out_of_core(
        &spectra,
        &idx,
        &params,
        &scorer,
        0.5,
        "XXX_",
        cache_dir.path(),
    )
    .expect("search");
    (spectra, queues, candidates, idx, params)
}

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

    let pin_path = tmp.path().join("bsa.pin");
    output::write_pin(&pin_path, &spectra, &queues, &candidates, &params, &idx).unwrap();
    let spec_ids = pin_spec_ids(&pin_path);
    assert!(
        spec_ids.len() >= 5,
        "expected several PIN rows, got {}",
        spec_ids.len()
    );

    let has_multi = spec_ids.iter().any(|s| s.matches('_').count() >= 3);
    eprintln!(
        "PIN rows={} multi_row_present={}",
        spec_ids.len(),
        has_multi
    );

    let mut map: HashMap<String, PercolatorPsm> = HashMap::new();
    for (i, sid) in spec_ids.iter().enumerate() {
        map.insert(
            sid.clone(),
            PercolatorPsm {
                pooled_q_value: None,
                psm_id: sid.clone(),
                q_value: 0.001 * (i as f64 + 1.0),
                pep: 0.0001 * (i as f64 + 1.0),
                score: None,
                peptide: "X".into(),
                proteins: "P".into(),
            },
        );
    }

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
        None,
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

    let pin_set: HashSet<&String> = spec_ids.iter().collect();
    assert_eq!(
        pin_set.len(),
        map.len(),
        "mock map should have one entry per distinct PIN SpecId"
    );

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
        None,
    )
    .unwrap();
    assert_eq!(
        read_pep_column(&qdir_ex.join("psms.parquet")).len(),
        spec_ids.len() - 1
    );
    let qdir_unscored = tmp.path().join("excluded_unscored.idparquet");
    output::write_qpx(
        &qdir_unscored,
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
        &excluded,
        None,
    )
    .unwrap();
    let unscored = read_pep_column(&qdir_unscored.join("psms.parquet"));
    assert_eq!(unscored.len(), spec_ids.len() - 1);
    assert!(unscored.iter().all(|p| p.is_none()));

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
        None,
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
                    pooled_q_value: None,
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
        None,
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
        None,
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
        None,
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
            None,
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

fn competition_result(pin_path: &Path) -> (HashMap<String, PercolatorPsm>, HashSet<String>, usize) {
    competition_result_by(pin_path, true)
}

fn competition_result_by(
    pin_path: &Path,
    per_input: bool,
) -> (HashMap<String, PercolatorPsm>, HashSet<String>, usize) {
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
    let mut best: HashMap<(String, String, String), (f64, String, bool)> = HashMap::new();
    let mut total_targets = 0usize;
    for line in lines.filter(|l| !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        total_targets += usize::from(f[label] != "-1");
        let s: f64 = f[score].parse().unwrap();
        let entry = (s, f[id].to_string(), f[label] == "-1");
        let run = if per_input {
            output::spec_id_run(f[id]).to_string()
        } else {
            String::new()
        };
        best.entry((run, f[scan].to_string(), f[exp].to_string()))
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
                    pooled_q_value: None,
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

fn add_runner_up(
    queues: &mut [search::psm::TopNQueue],
    candidates: &[search::candidate_gen::Candidate],
) {
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
}

#[test]
fn qpx_join_accepts_targets_removed_by_target_decoy_competition() {
    let (spectra, mut queues, candidates, idx, params) = run_bsa_search();
    add_runner_up(&mut queues, &candidates);
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
            None,
        )
    };
    let err = write(&tmp.path().join("strict.idparquet"), None).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);

    let qdir = tmp.path().join("competed.idparquet");
    write(&qdir, Some(&decoys)).unwrap();
    let peps = read_pep_column(&qdir.join("psms.parquet"));
    assert_eq!(peps.len(), pin_rows - removed);
    assert_eq!(peps.iter().filter(|p| p.is_some()).count(), targets.len());
}

#[test]
fn qpx_join_reads_competition_groups_per_input() {
    let (spectra, mut queues, candidates, idx, params) = run_bsa_search();
    add_runner_up(&mut queues, &candidates);
    let mut two_inputs = Vec::new();
    for run in ["run_a", "run_b"] {
        for s in &spectra {
            let mut s = s.clone();
            let name = if s.title.is_empty() {
                format!("scan={}", s.scan.unwrap_or(0))
            } else {
                s.title.clone()
            };
            s.title = format!("{run}/{name}");
            two_inputs.push(s);
        }
    }
    let mut two_queues = queues.clone();
    two_queues.extend(queues.iter().cloned());
    let run_ends = [spectra.len(), two_inputs.len()];
    let paths = ["run_a.mgf".to_string(), "run_b.mgf".to_string()];

    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("two.pin");
    output::write_pin(&pin, &two_inputs, &two_queues, &candidates, &params, &idx).unwrap();
    let pin_rows = pin_spec_ids(&pin).len();
    let pin_text = std::fs::read_to_string(&pin).unwrap();
    let header: Vec<&str> = pin_text.lines().next().unwrap().split('\t').collect();
    let (id_col, file_col) = (
        header.iter().position(|c| *c == "SpecId").unwrap(),
        header.iter().position(|c| *c == "FileName").unwrap(),
    );
    assert_eq!(
        file_col,
        header.iter().position(|c| *c == "CalcMass").unwrap() + 1
    );
    let mut files = HashSet::new();
    for line in pin_text.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[file_col], output::spec_id_run(f[id_col]));
        files.insert(f[file_col].to_string());
    }
    assert_eq!(
        files,
        HashSet::from(["run_a".to_string(), "run_b".to_string()])
    );
    let (targets, decoys, removed) = competition_result(&pin);
    assert!(
        removed > 0,
        "the fixture must have groups with several rows"
    );

    let write = |dir: &Path, targets: &HashMap<String, PercolatorPsm>, decoys: &HashSet<String>| {
        output::write_qpx(
            dir,
            &two_inputs,
            &two_queues,
            &candidates,
            &params,
            &idx,
            &Tolerance::Da(0.5),
            "test_run",
            &paths,
            &run_ends,
            Some(targets),
            Some(decoys),
            &HashSet::new(),
            None,
        )
    };
    let qdir = tmp.path().join("per_input.idparquet");
    write(&qdir, &targets, &decoys).unwrap();
    let peps = read_pep_column(&qdir.join("psms.parquet"));
    assert_eq!(peps.len(), pin_rows - removed);
    assert_eq!(peps.iter().filter(|p| p.is_some()).count(), targets.len());

    let (pooled_targets, pooled_decoys, _) = competition_result_by(&pin, false);
    let err = write(
        &tmp.path().join("pooled.idparquet"),
        &pooled_targets,
        &pooled_decoys,
    )
    .unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
}
