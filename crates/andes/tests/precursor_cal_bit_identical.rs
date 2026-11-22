use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("canonicalize workspace root")
}

fn read_sorted_rows(path: &PathBuf) -> (String, Vec<String>) {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let text =
        String::from_utf8(bytes).unwrap_or_else(|_| panic!("non-utf8 in {}", path.display()));
    let mut lines = text.lines();
    let header = lines
        .next()
        .unwrap_or_else(|| panic!("empty file: {}", path.display()))
        .to_string();
    let mut rows: Vec<String> = lines.map(|s| s.to_string()).collect();
    rows.sort();
    (header, rows)
}

#[test]
fn precursor_cal_off_pin_tsv_match_golden_after_sort() {
    let root = workspace_root();
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let spectra = root.join("test-fixtures/test.mgf.gz");
    let fasta = root.join("test-fixtures/BSA.fasta");
    assert!(spectra.exists(), "fixture missing: {}", spectra.display());
    assert!(fasta.exists(), "fixture missing: {}", fasta.display());

    let outdir = tempfile::tempdir().expect("tempdir");
    let actual_pin = outdir.path().join("actual.pin");
    let actual_tsv = outdir.path().join("actual.tsv");

    let status = Command::new(&binary)
        .arg("--spectrum")
        .arg(&spectra)
        .arg("--database")
        .arg(&fasta)
        .arg("--output-pin")
        .arg(&actual_pin)
        .arg("--output-tsv")
        .arg(&actual_tsv)
        .args(["--predict", "none"])
        .args(["--max-missed-cleavages", "1"])
        .status()
        .expect("run andes");
    assert!(status.success(), "andes exited {status}");

    let golden_pin = root.join("test-fixtures/parity/goldens/precursor_cal_off.pin");
    let golden_tsv = root.join("test-fixtures/parity/goldens/precursor_cal_off.tsv");

    let (actual_pin_hdr, actual_pin_rows) = read_sorted_rows(&actual_pin);
    let (golden_pin_hdr, golden_pin_rows) = read_sorted_rows(&golden_pin);
    let (actual_tsv_hdr, actual_tsv_rows) = read_sorted_rows(&actual_tsv);
    let (golden_tsv_hdr, golden_tsv_rows) = read_sorted_rows(&golden_tsv);

    assert_eq!(actual_pin_hdr, golden_pin_hdr, "PIN header drift");
    assert_eq!(actual_tsv_hdr, golden_tsv_hdr, "TSV header drift");
    assert_eq!(
        actual_pin_rows.len(),
        golden_pin_rows.len(),
        "PIN row count drift: actual {} vs golden {}",
        actual_pin_rows.len(),
        golden_pin_rows.len(),
    );
    assert_eq!(
        actual_tsv_rows.len(),
        golden_tsv_rows.len(),
        "TSV row count drift: actual {} vs golden {}",
        actual_tsv_rows.len(),
        golden_tsv_rows.len(),
    );

    if actual_pin_rows != golden_pin_rows {
        for (i, (a, g)) in actual_pin_rows
            .iter()
            .zip(golden_pin_rows.iter())
            .enumerate()
        {
            if a != g {
                panic!("PIN row mismatch at sorted index {i}:\n  actual: {a}\n  golden: {g}",);
            }
        }
    }
    if actual_tsv_rows != golden_tsv_rows {
        for (i, (a, g)) in actual_tsv_rows
            .iter()
            .zip(golden_tsv_rows.iter())
            .enumerate()
        {
            if a != g {
                panic!("TSV row mismatch at sorted index {i}:\n  actual: {a}\n  golden: {g}",);
            }
        }
    }
}

fn calibration_outcomes(stderr: &str) -> usize {
    stderr
        .lines()
        .filter(|l| {
            l.contains("Precursor mass shift learned:")
                || l.contains("Precursor mass calibration skipped")
        })
        .count()
}

/// Two spectrum files: `--precursor-cal on` runs the calibration pre-pass on each file, instead
/// of learning it on the first file and applying it to both.
#[test]
#[ignore = "two E. coli searches; run with --ignored"]
fn precursor_cal_is_learned_per_spectrum_file() {
    use std::io::Read;
    let root = workspace_root();
    let mut text = String::new();
    flate2::read::GzDecoder::new(
        std::fs::File::open(root.join("test-fixtures/test.mgf.gz")).expect("open fixture"),
    )
    .read_to_string(&mut text)
    .expect("inflate fixture");
    let shifted: String = text
        .lines()
        .map(|l| match l.strip_prefix("PEPMASS=") {
            Some(v) => {
                let mz: f64 = v.split_whitespace().next().unwrap().parse().unwrap();
                format!("PEPMASS={:.9}\n", mz * (1.0 + 8e-6))
            }
            None => format!("{l}\n"),
        })
        .collect();
    let dir = tempfile::tempdir().expect("tempdir");
    let second = dir.path().join("shifted.mgf");
    std::fs::write(&second, shifted).expect("write shifted copy");
    let output = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(root.join("test-fixtures/test.mgf.gz"))
        .arg("--spectrum")
        .arg(&second)
        .arg("--database")
        .arg(root.join("test-fixtures/ecoli.fasta"))
        .args(["--fragmentation", "HCD", "--fragment-tol-ppm", "20"])
        .args(["--predict", "none", "--max-missed-cleavages", "1"])
        .args(["--precursor-cal", "on", "--cal-min-spec-keys", "1000"])
        .arg("--output-pin")
        .arg(dir.path().join("run.pin"))
        .output()
        .expect("run andes");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "andes failed:\n{stderr}");
    assert!(stderr.contains("calibration: spectrum file 2/2"), "{stderr}");
    assert_eq!(calibration_outcomes(&stderr), 2, "one pre-pass per file:\n{stderr}");
}
