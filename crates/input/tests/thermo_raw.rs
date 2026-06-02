//! Integration test for the native Thermo `.raw` reader.
//!
//! Requires the `thermo` feature, the .NET 8 runtime, and a real `.raw` file
//! pointed to by the `MSGF_TEST_RAW` environment variable. The `.raw` format is
//! proprietary and the files are large, so no fixture is committed; the test
//! skips (passes) when `MSGF_TEST_RAW` is unset, and runs in environments that
//! have a `.raw` available (e.g. the benchmark host).
//!
//! Run with:  `MSGF_TEST_RAW=/path/sample.raw cargo test -p input --features thermo`

#![cfg(feature = "thermo")]

use input::ThermoRawReader;

#[test]
fn thermo_raw_reads_ms2_spectra() {
    let path = match std::env::var("MSGF_TEST_RAW") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("skipping thermo_raw_reads_ms2_spectra: set MSGF_TEST_RAW to a .raw to run");
            return;
        }
    };

    let reader = ThermoRawReader::open(&path).expect("open .raw");
    assert!(reader.len() > 0, "the file should report > 0 spectra");

    let mut ms2 = 0usize;
    let mut checked_first = false;
    for spec in reader {
        let spec = spec.expect("read spectrum");
        ms2 += 1;
        if !checked_first {
            checked_first = true;
            // First emitted MS2 must carry a precursor, peaks, a scan number,
            // and a Thermo native-id title.
            assert!(spec.precursor_mz > 0.0, "MS2 must have a precursor m/z");
            assert!(!spec.peaks.is_empty(), "MS2 must have peaks");
            assert!(spec.scan.is_some(), "MS2 must carry a scan number");
            assert!(
                spec.title.contains("scan="),
                "title should be the Thermo native id, got {:?}",
                spec.title
            );
            // Peaks must be finite and intensities non-negative.
            for &(mz, intensity) in spec.peaks.iter().take(50) {
                assert!(mz.is_finite() && mz > 0.0, "peak m/z must be positive finite");
                assert!(intensity >= 0.0, "peak intensity must be non-negative");
            }
        }
    }
    assert!(ms2 > 0, "the default reader should emit at least one MS2 spectrum");
}

/// The MS1-linked chunked reader (used by `--chimeric`) must emit every MS2 and
/// keep `ms2_to_ms1` aligned 1:1 with the emitted spectra in each chunk.
#[test]
fn thermo_raw_chunked_ms1_link_is_consistent() {
    let path = match std::env::var("MSGF_TEST_RAW") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("skipping thermo_raw_chunked_ms1_link_is_consistent: set MSGF_TEST_RAW");
            return;
        }
    };

    let reader = ThermoRawReader::open(&path).expect("open .raw");
    let mut total_ms2 = 0usize;
    let mut chunks = 0usize;
    reader.read_with_ms1_chunked(5000, 0, |chunk, link| {
        chunks += 1;
        // One link entry per emitted MS2.
        assert_eq!(
            chunk.len(),
            link.ms2_to_ms1.len(),
            "ms2_to_ms1 length must equal the MS2 count in the chunk"
        );
        // Every link index is in range of the chunk's MS1 set (or None).
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
