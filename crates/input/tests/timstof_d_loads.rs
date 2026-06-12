//! Integration test for the native Bruker timsTOF `.d` reader.
//!
//! Built only under `--features timstof`. Reading a real `.d` needs a Bruker
//! dataset on disk (each `.d` is ~1-3.5 GB), so this test is a no-op unless the
//! `ANDES_TEST_D` env var points at a `.d` directory — keeping CI green on
//! machines without one. With the var set it opens the `.d`, streams the MS2
//! spectra, and asserts the model invariants the search path relies on.
//!
//! Example:
//!   ANDES_TEST_D=/data/HeLa_IAA_F51_1.d \
//!     cargo test -p input --features timstof --test timstof_d_loads

#![cfg(feature = "timstof")]

use input::TimsTofReader;

fn test_d_path() -> Option<std::ffi::OsString> {
    std::env::var_os("ANDES_TEST_D").or_else(|| std::env::var_os("MSGF_TEST_D"))
}

#[test]
fn reads_real_d_when_env_set() {
    let Some(path) = test_d_path() else {
        eprintln!("ANDES_TEST_D not set — skipping real .d read test");
        return;
    };

    let reader = TimsTofReader::open(&path).expect("open .d directory");
    assert!(!reader.is_empty(), "expected at least one MS2 spectrum in the .d");

    let mut count = 0usize;
    for (i, result) in reader.enumerate().take(2000) {
        let spec = result.unwrap_or_else(|e| panic!("read spectrum {i}: {e}"));

        // Every emitted spectrum is a searchable MS2: positive precursor m/z.
        assert!(
            spec.precursor_mz > 0.0,
            "spectrum {i} emitted with non-positive precursor m/z {}",
            spec.precursor_mz
        );

        // Peaks must be ascending by m/z — the scorer and chimeric co-isolation
        // detection binary-search this list.
        assert!(
            spec.peaks.windows(2).all(|w| w[0].0 <= w[1].0),
            "spectrum {i} peaks are not ascending by m/z"
        );

        // A real scan number is carried for the PIN `ScanNum`/`SpecID` columns.
        assert!(spec.scan.is_some(), "spectrum {i} has no scan number");

        count += 1;
    }
    assert!(count > 0, "no spectra were read from the .d");
    eprintln!("read {count} MS2 spectra from {:?}", path);
}
