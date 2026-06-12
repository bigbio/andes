//! Handcrafted MGF strings exercising parser edge cases.

use std::io::Cursor;
use input::{MgfParseError, MgfReader, Spectrum};

fn parse_all(s: &str) -> Vec<Result<Spectrum, MgfParseError>> {
    MgfReader::new(Cursor::new(s)).collect()
}

fn parse_ok(s: &str) -> Vec<Spectrum> {
    parse_all(s).into_iter().map(|r| r.unwrap()).collect()
}

#[test]
fn empty_input_emits_nothing() {
    let v = parse_ok("");
    assert!(v.is_empty());
}

#[test]
fn single_minimal_spectrum() {
    let mgf = "BEGIN IONS\n\
               TITLE=test\n\
               PEPMASS=500.5\n\
               100.0 1.0\n\
               200.0 2.0\n\
               END IONS\n";
    let v = parse_ok(mgf);
    assert_eq!(v.len(), 1);
    assert_eq!(v[0].title, "test");
    assert_eq!(v[0].precursor_mz, 500.5);
    assert_eq!(v[0].peaks, vec![(100.0, 1.0), (200.0, 2.0)]);
    assert!(v[0].precursor_charge.is_none());
}

#[test]
fn full_spectrum_with_all_fields() {
    let mgf = "BEGIN IONS\n\
               TITLE=Scan 42\n\
               PEPMASS=500.5 1000.0\n\
               CHARGE=2+\n\
               RTINSECONDS=120.5\n\
               SCANS=42\n\
               100.0 1.0\n\
               END IONS\n";
    let v = parse_ok(mgf);
    assert_eq!(v.len(), 1);
    let s = &v[0];
    assert_eq!(s.title, "Scan 42");
    assert_eq!(s.precursor_mz, 500.5);
    assert_eq!(s.precursor_intensity, Some(1000.0));
    assert_eq!(s.precursor_charge, Some(2));
    assert_eq!(s.rt_seconds, Some(120.5));
    assert_eq!(s.scan, Some(42));
}

#[test]
fn charge_strips_sign() {
    for (line, expected) in [("CHARGE=2+", 2), ("CHARGE=3+", 3), ("CHARGE=1-", 1)] {
        let mgf = format!(
            "BEGIN IONS\nTITLE=x\nPEPMASS=500\n{}\n100 1\nEND IONS\n", line);
        let v = parse_ok(&mgf);
        assert_eq!(v[0].precursor_charge, Some(expected), "line={line}");
    }
}

#[test]
fn multiple_spectra() {
    let mgf = "BEGIN IONS\n\
               TITLE=a\n\
               PEPMASS=100\n\
               1 1\n\
               END IONS\n\
               BEGIN IONS\n\
               TITLE=b\n\
               PEPMASS=200\n\
               2 2\n\
               END IONS\n";
    let v = parse_ok(mgf);
    assert_eq!(v.len(), 2);
    assert_eq!(v[0].title, "a");
    assert_eq!(v[1].title, "b");
}

#[test]
fn comments_and_blank_lines_ignored() {
    let mgf = "# leading comment\n\
               \n\
               BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=100\n\
               1 1\n\
               END IONS\n\
               # trailing comment\n";
    let v = parse_ok(mgf);
    assert_eq!(v.len(), 1);
}

#[test]
fn unknown_keys_tolerated() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=100\n\
               CUSTOM_KEY=anything goes\n\
               INSTRUMENT=Q-Exactive\n\
               1 1\n\
               END IONS\n";
    let v = parse_ok(mgf);
    assert_eq!(v.len(), 1);
}

#[test]
fn pepmass_without_intensity() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=500.5\n\
               100 1\n\
               END IONS\n";
    let v = parse_ok(mgf);
    assert_eq!(v[0].precursor_mz, 500.5);
    assert!(v[0].precursor_intensity.is_none());
}

#[test]
fn empty_title_is_ok() {
    let mgf = "BEGIN IONS\n\
               TITLE=\n\
               PEPMASS=100\n\
               1 1\n\
               END IONS\n";
    let v = parse_ok(mgf);
    assert_eq!(v[0].title, "");
}

#[test]
fn peaks_sorted_ascending_by_mz() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=100\n\
               300 3\n\
               100 1\n\
               200 2\n\
               END IONS\n";
    let v = parse_ok(mgf);
    let mzs: Vec<_> = v[0].peaks.iter().map(|p| p.0).collect();
    assert_eq!(mzs, vec![100.0, 200.0, 300.0]);
}

#[test]
fn tab_separator_in_peak_lines() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=100\n\
               100\t1\n\
               END IONS\n";
    let v = parse_ok(mgf);
    assert_eq!(v[0].peaks, vec![(100.0, 1.0)]);
}

#[test]
fn missing_pepmass_errors() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               100 1\n\
               END IONS\n";
    let err = parse_all(mgf).into_iter().next().unwrap().unwrap_err();
    assert!(matches!(err, MgfParseError::MissingPepmass { .. }));
}

#[test]
fn bad_pepmass_errors() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=garbage\n\
               100 1\n\
               END IONS\n";
    let err = parse_all(mgf).into_iter().next().unwrap().unwrap_err();
    assert!(matches!(err, MgfParseError::BadPepmass { .. }));
}

#[test]
fn bad_charge_errors() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=100\n\
               CHARGE=banana\n\
               100 1\n\
               END IONS\n";
    let err = parse_all(mgf).into_iter().next().unwrap().unwrap_err();
    assert!(matches!(err, MgfParseError::BadCharge { .. }));
}

#[test]
fn bad_peak_errors() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=100\n\
               not a peak line\n\
               END IONS\n";
    let err = parse_all(mgf).into_iter().next().unwrap().unwrap_err();
    assert!(matches!(err, MgfParseError::BadPeak { .. }));
}

#[test]
fn unterminated_spectrum_errors() {
    let mgf = "BEGIN IONS\n\
               TITLE=x\n\
               PEPMASS=100\n\
               100 1\n";
    let err = parse_all(mgf).into_iter().next().unwrap().unwrap_err();
    assert!(matches!(err, MgfParseError::UnterminatedSpectrum { .. }));
}
