//! Tests of where the rescoring PIN is written.

use super::inputs::rescore_pin_destination;
use std::path::{Path, PathBuf};

#[test]
fn explicit_output_pin_is_pin_and_base() {
    let mut tmp = None;
    let (pin, base) = rescore_pin_destination(
        Some(Path::new("r/x.pin")),
        None,
        None,
        false,
        Path::new("."),
        &mut tmp,
    )
    .unwrap();
    assert_eq!(pin, PathBuf::from("r/x.pin"));
    assert_eq!(base, PathBuf::from("r/x.pin"));
    assert!(tmp.is_none());
}

#[test]
fn keep_pin_writes_next_to_report_base() {
    let mut tmp = None;
    let (pin, base) = rescore_pin_destination(
        None,
        Some(Path::new("r/out.tsv")),
        Some(Path::new("r/out.idparquet")),
        true,
        Path::new("."),
        &mut tmp,
    )
    .unwrap();
    assert_eq!(pin, PathBuf::from("r/out.pin"));
    assert_eq!(base, PathBuf::from("r/out.tsv"));
    assert!(tmp.is_none());

    let (pin, _) = rescore_pin_destination(
        None,
        None,
        Some(Path::new("q.idparquet")),
        true,
        Path::new("."),
        &mut tmp,
    )
    .unwrap();
    assert_eq!(pin, PathBuf::from("q.pin"));
    let (pin, base) =
        rescore_pin_destination(None, None, None, true, Path::new("."), &mut tmp).unwrap();
    assert_eq!(pin, PathBuf::from("andes.pin"));
    assert_eq!(base, PathBuf::from("andes.pin"));
    assert!(tmp.is_none());
}

#[test]
fn no_keep_pin_is_created_in_the_scratch_dir() {
    let scratch = tempfile::tempdir().unwrap();
    let mut tmp = None;
    let (pin, base) = rescore_pin_destination(
        None,
        Some(Path::new("out.tsv")),
        None,
        false,
        scratch.path(),
        &mut tmp,
    )
    .unwrap();
    let dir = tmp.as_ref().expect("temp dir held").path().to_path_buf();
    assert!(dir.starts_with(scratch.path()));
    assert!(pin.starts_with(&dir));
    assert_eq!(base, PathBuf::from("out.tsv"));
    drop(tmp);
    assert_eq!(std::fs::read_dir(scratch.path()).unwrap().count(), 0);
}

#[test]
fn unwritable_scratch_dir_is_an_error_naming_the_flag() {
    let dir = tempfile::tempdir().unwrap();
    let blocker = dir.path().join("blocker");
    std::fs::write(&blocker, b"x").unwrap();
    let mut tmp = None;
    let err = rescore_pin_destination(
        None,
        Some(Path::new("out.tsv")),
        None,
        false,
        &blocker.join("sub"),
        &mut tmp,
    )
    .unwrap_err();
    assert!(err.to_string().contains("--temp-dir"), "{err}");
}
