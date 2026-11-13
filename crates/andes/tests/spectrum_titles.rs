//! A spectrum title is the `SpecId` of the PIN: a tab in a title must not add
//! a column to the row.

use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

/// The first `n` spectra of the MGF fixture, every title given a tab.
fn mgf_with_tabbed_titles(n: usize) -> String {
    let gz = std::fs::File::open(fixture("test-fixtures/test.mgf.gz")).expect("gz");
    let mut text = String::new();
    let mut spectra = 0;
    let mut titles = 0;
    for line in BufReader::new(flate2::read::GzDecoder::new(gz)).lines() {
        let line = line.expect("inflate");
        match line.strip_prefix("TITLE=") {
            Some(title) => {
                titles += 1;
                text.push_str(&format!("TITLE=spot\t{titles} {title}"));
            }
            None => text.push_str(&line),
        }
        text.push('\n');
        if line.starts_with("END IONS") {
            spectra += 1;
            if spectra == n {
                break;
            }
        }
    }
    assert_eq!(spectra, n, "the fixture holds fewer spectra");
    assert_eq!(titles, n, "a spectrum of the fixture has no title");
    text
}

#[test]
fn a_tab_in_a_title_does_not_shift_the_pin_columns() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mgf = dir.path().join("spectra.mgf");
    std::fs::write(&mgf, mgf_with_tabbed_titles(150)).expect("write MGF");
    let pin = dir.path().join("out.pin");
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(&mgf)
        .arg("--database")
        .arg(fixture("test-fixtures/ecoli.fasta"))
        .arg("--output-pin")
        .arg(&pin)
        .output()
        .expect("run andes");
    assert!(
        out.status.success(),
        "andes failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = std::fs::read_to_string(&pin).expect("read PIN");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("PIN header").split('\t').collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|c| *c == name)
            .unwrap_or_else(|| panic!("no {name} column"))
    };
    let (spec_id, label, scan) = (column("SpecId"), column("Label"), column("ScanNr"));
    // `Proteins` is the last column and holds one field per protein.
    assert_eq!(column("Proteins"), header.len() - 1);
    let mut rows = 0;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= header.len(), "{line}");
        assert!(fields[spec_id].starts_with("spectra/spot "), "{line}");
        assert!(matches!(fields[label], "1" | "-1"), "{line}");
        assert!(fields[scan].parse::<i64>().is_ok(), "{line}");
        rows += 1;
    }
    assert!(rows > 0, "the search wrote no PIN rows");
}
