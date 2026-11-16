#![cfg(feature = "timstof")]

use input::TimsTofReader;

fn test_d_path() -> Option<std::ffi::OsString> {
    std::env::var_os("ANDES_TEST_D")
}

#[test]
fn reads_real_d_when_env_set() {
    let Some(path) = test_d_path() else {
        eprintln!("ANDES_TEST_D not set — skipping real .d read test");
        return;
    };

    let reader = TimsTofReader::open(&path).expect("open .d directory");
    assert!(
        !reader.is_empty(),
        "expected at least one MS2 spectrum in the .d"
    );

    let mut count = 0usize;
    for (i, result) in reader.enumerate().take(2000) {
        let spec = result.unwrap_or_else(|e| panic!("read spectrum {i}: {e}"));

        assert!(
            spec.precursor_mz > 0.0,
            "spectrum {i} emitted with non-positive precursor m/z {}",
            spec.precursor_mz
        );

        assert!(
            spec.peaks.windows(2).all(|w| w[0].0 <= w[1].0),
            "spectrum {i} peaks are not ascending by m/z"
        );

        assert!(spec.scan.is_some(), "spectrum {i} has no scan number");

        count += 1;
    }
    assert!(count > 0, "no spectra were read from the .d");
    eprintln!("read {count} MS2 spectra from {:?}", path);
}
