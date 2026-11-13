//! `CHARGE=` values of an MGF: a list of several charge states leaves the
//! charge unknown, so the configured charge range is searched, and a negative
//! charge stops the search.

use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

/// The first `n` spectra of the MGF fixture with every `CHARGE=` value
/// replaced by `charge`.
fn mgf_with_charge(n: usize, charge: &str) -> String {
    let gz = std::fs::File::open(fixture("test-fixtures/test.mgf.gz")).expect("gz");
    let mut text = String::new();
    let mut spectra = 0;
    for line in BufReader::new(flate2::read::GzDecoder::new(gz)).lines() {
        let line = line.expect("inflate");
        if line.starts_with("CHARGE=") {
            text.push_str("CHARGE=");
            text.push_str(charge);
        } else {
            text.push_str(&line);
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
    text
}

/// Search `mgf` (written as `<dir>/spectra.mgf`) against the E. coli fixture.
fn search(dir: &Path, mgf: &str) -> Output {
    std::fs::create_dir_all(dir).expect("mkdir");
    let spectra = dir.join("spectra.mgf");
    std::fs::write(&spectra, mgf).expect("write MGF");
    Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(&spectra)
        .arg("--database")
        .arg(fixture("test-fixtures/ecoli.fasta"))
        .arg("--output-pin")
        .arg(dir.join("out.pin"))
        .output()
        .expect("run andes")
}

fn pin_of(dir: &Path, mgf: &str) -> String {
    let out = search(dir, mgf);
    assert!(
        out.status.success(),
        "andes failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::read_to_string(dir.join("out.pin")).expect("read PIN")
}

/// PIN rows per charge one-hot column, as `(column name, rows set)`.
fn rows_per_charge(pin: &str) -> Vec<(String, usize)> {
    let mut lines = pin.lines();
    let header: Vec<&str> = lines.next().expect("PIN header").split('\t').collect();
    let columns: Vec<(usize, &str)> = header
        .iter()
        .enumerate()
        .filter(|(_, name)| {
            name.strip_prefix("charge")
                .is_some_and(|z| z.parse::<u8>().is_ok())
        })
        .map(|(i, name)| (i, *name))
        .collect();
    assert!(!columns.is_empty(), "no charge columns in {header:?}");
    let mut counts = vec![0usize; columns.len()];
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        for (slot, (i, _)) in columns.iter().enumerate() {
            if fields[*i] == "1" {
                counts[slot] += 1;
            }
        }
    }
    columns
        .iter()
        .zip(counts)
        .map(|((_, name), n)| ((*name).to_string(), n))
        .collect()
}

#[test]
fn several_charge_states_are_searched_as_an_unknown_charge() {
    let dir = tempfile::tempdir().expect("tempdir");
    let listed = pin_of(
        &dir.path().join("listed"),
        &mgf_with_charge(150, "2+ and 3+"),
    );
    let unknown = pin_of(&dir.path().join("unknown"), &mgf_with_charge(150, "0"));
    assert!(listed.lines().count() > 1, "no PIN rows");
    assert_eq!(listed, unknown, "a charge list is not searched as unknown");

    let counts = rows_per_charge(&listed);
    let beyond_first_state: usize = counts
        .iter()
        .filter(|(name, _)| name != "charge2")
        .map(|(_, n)| n)
        .sum();
    assert!(
        beyond_first_state > 0,
        "every match is at the first listed charge: {counts:?}"
    );
}

#[test]
fn negative_charge_stops_the_search() {
    let dir = tempfile::tempdir().expect("tempdir");
    let out = search(dir.path(), &mgf_with_charge(5, "2-"));
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success(), "the search ran:\n{err}");
    assert!(err.contains("negative CHARGE at line"), "{err}");
    assert!(
        err.contains("negative-ion spectra are not supported"),
        "{err}"
    );
}
