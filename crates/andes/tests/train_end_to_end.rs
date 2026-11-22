use std::path::PathBuf;
use std::process::Command;

use model_train::ModelStore;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

#[test]
fn train_writes_model_and_search_uses_it() {
    let dir = tempfile::tempdir().expect("tempdir");
    let store_path = dir.path().join("m.parquet");
    let pin_path = dir.path().join("out.pin");

    let bsa_mgf = fixture("test-fixtures/test.mgf.gz");
    let bsa_fasta = fixture("test-fixtures/BSA.fasta");

    let train_status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("train-from-search")
        .arg("--spectra")
        .arg(&bsa_mgf)
        .arg("--database")
        .arg(&bsa_fasta)
        .arg("--out-store")
        .arg(&store_path)
        .arg("--train-fdr")
        .arg("0.9")
        .arg("--model-id")
        .arg("bsa_test")
        .status()
        .expect("run andes train");

    assert!(
        train_status.success(),
        "andes train should exit 0, got: {train_status}"
    );

    assert!(store_path.exists(), "model store should be written");

    let store = ModelStore::open(&store_path).expect("open trained model store");
    let ids = store.model_ids();
    assert!(
        ids.contains(&"bsa_test".to_string()),
        "store should contain model 'bsa_test'; found: {ids:?}"
    );

    let param = store.load_param("bsa_test").expect("load bsa_test param");
    assert!(
        !param.partitions.is_empty(),
        "trained param should have at least one partition"
    );
    assert!(
        param.frag_intensity_model.is_some() && param.rich_ion_model.is_some(),
        "a trained model keeps the seed's fragment-intensity and rich-ion models"
    );

    let search_status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(&bsa_mgf)
        .arg("--database")
        .arg(&bsa_fasta)
        .arg("--output-pin")
        .arg(&pin_path)
        .arg("--model-store")
        .arg(&store_path)
        .arg("--model")
        .arg("bsa_test")
        .status()
        .expect("run search with trained model");

    assert!(
        search_status.success(),
        "search with trained model should exit 0, got: {search_status}"
    );
    assert!(pin_path.exists(), "PIN file should be written");

    let stats_path = dir.path().join("out.statistics.log");
    assert!(
        stats_path.exists(),
        "out.statistics.log should be written next to the PIN"
    );
    let stats = std::fs::read_to_string(&stats_path).expect("read statistics.log");
    assert!(
        stats.contains("Final precursor tolerance") && stats.contains("PTM report"),
        "statistics.log should contain the run summary; got:\n{stats}"
    );

    let pin_content = std::fs::read_to_string(&pin_path).expect("read PIN");
    let line_count = pin_content.lines().count();
    assert!(
        line_count > 1,
        "PIN should have header + at least 1 data row; got {line_count} line(s)"
    );
    let header = pin_content.lines().next().unwrap_or("");
    assert!(
        header.starts_with("SpecId\tLabel\tScanNr"),
        "unexpected PIN header: {header}"
    );
}

#[test]
fn train_from_several_runs_records_one_source_per_run() {
    let dir = tempfile::tempdir().expect("tempdir");
    let store_path = dir.path().join("m.parquet");
    let bsa_mgf = fixture("test-fixtures/test.mgf.gz");
    let bsa_fasta = fixture("test-fixtures/BSA.fasta");

    let status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("train-from-search")
        .arg("--spectra")
        .arg(&bsa_mgf)
        .arg("--spectra")
        .arg(&bsa_mgf)
        .arg("--database")
        .arg(&bsa_fasta)
        .arg("--out-store")
        .arg(&store_path)
        .arg("--train-fdr")
        .arg("0.9")
        .arg("--fragment-tol-da")
        .arg("0.5")
        .arg("--instrument")
        .arg("TimsTOF")
        .arg("--protocol")
        .arg("TMT")
        .arg("--model-id")
        .arg("bsa_two_runs")
        .status()
        .expect("run andes train-from-search");
    assert!(
        status.success(),
        "train-from-search should exit 0, got: {status}"
    );

    let store = ModelStore::open(&store_path).expect("open trained model store");
    let sources = store.load_sources("bsa_two_runs").expect("load sources");
    let mut ids: Vec<&str> = sources.iter().map(|l| l.source_id.as_str()).collect();
    ids.sort_unstable();
    assert_eq!(
        ids,
        ["test", "test_2"],
        "one source per run, stems made unique"
    );
    assert_eq!(sources[0].n_psms, sources[1].n_psms);
    assert!(sources[0].n_psms > 0);
    assert!(sources.iter().all(|l| l.dataset == "test.mgf.gz"));
    assert!(sources.iter().all(|l| l.instrument == "TimsTOF"));
    assert!(sources.iter().all(|l| l.experiment_class == "TMT"));

    let param = store.load_param("bsa_two_runs").expect("load param");
    assert_eq!(param.data_type.instrument.name(), "TimsTOF");
    assert_eq!(
        param.data_type.protocol.name(),
        "TMT",
        "the trained model carries --protocol, so routing selects it for TMT searches"
    );
    assert_eq!(param.data_type.activation.name(), "HCD");
    assert_eq!(format!("{:?}", param.mme), "Da(0.5)");
}
