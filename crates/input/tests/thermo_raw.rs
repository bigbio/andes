#![cfg(feature = "thermo")]

use input::ThermoRawReader;

fn test_raw_path() -> Option<String> {
    std::env::var("ANDES_TEST_RAW").ok()
}

#[test]
fn thermo_raw_reads_ms2_spectra() {
    let path = match test_raw_path() {
        Some(p) => p,
        None => {
            eprintln!("skipping thermo_raw_reads_ms2_spectra: set ANDES_TEST_RAW to a .raw to run");
            return;
        }
    };

    let reader = ThermoRawReader::open(&path).expect("open .raw");
    assert!(!reader.is_empty(), "the file should report > 0 spectra");

    let mut ms2 = 0usize;
    let mut checked_first = false;
    for spec in reader {
        let spec = spec.expect("read spectrum");
        ms2 += 1;
        if !checked_first {
            checked_first = true;
            assert!(spec.precursor_mz > 0.0, "MS2 must have a precursor m/z");
            assert!(!spec.peaks.is_empty(), "MS2 must have peaks");
            assert!(spec.scan.is_some(), "MS2 must carry a scan number");
            assert!(
                spec.title.contains("scan="),
                "title should be the Thermo native id, got {:?}",
                spec.title
            );
            for &(mz, intensity) in spec.peaks.iter().take(50) {
                assert!(
                    mz.is_finite() && mz > 0.0,
                    "peak m/z must be positive finite"
                );
                assert!(intensity >= 0.0, "peak intensity must be non-negative");
            }
        }
    }
    assert!(
        ms2 > 0,
        "the default reader should emit at least one MS2 spectrum"
    );
}

#[test]
fn thermo_raw_chunked_ms1_link_is_consistent() {
    let path = match test_raw_path() {
        Some(p) => p,
        None => {
            eprintln!("skipping thermo_raw_chunked_ms1_link_is_consistent: set ANDES_TEST_RAW");
            return;
        }
    };

    let mut reader = ThermoRawReader::open(&path).expect("open .raw");
    let mut total_ms2 = 0usize;
    let mut chunks = 0usize;
    reader.read_with_ms1_chunked(5000, 0, |chunk, link| {
        chunks += 1;
        assert_eq!(
            chunk.len(),
            link.ms2_to_ms1.len(),
            "ms2_to_ms1 length must equal the MS2 count in the chunk"
        );
        for &m in &link.ms2_to_ms1 {
            if let Some(idx) = m {
                assert!(idx < link.ms1_peaks.len(), "ms1 link index out of range");
            }
        }
        total_ms2 += chunk.len();
    });
    assert!(total_ms2 > 0, "chunked read should emit MS2 spectra");
    assert!(chunks > 0, "chunked read should produce at least one chunk");
}
