//! Tests of the stage clock, the peak-memory reading, the input hasher and the command line.

use super::*;

#[test]
fn stage_clock_closes_each_stage_in_order() {
    let mut clock = StageClock::new();
    clock.stage("a");
    clock.stage("b");
    let rows = clock.finish();
    let names: Vec<&str> = rows.iter().map(|r| r.name.as_str()).collect();
    assert_eq!(names, ["a", "b"]);
    assert!(rows.iter().all(|r| r.seconds >= 0.0));
    if cfg!(any(target_os = "linux", target_os = "macos")) {
        assert!(rows.iter().all(|r| r.peak_rss_bytes.is_some()), "{rows:?}");
    }
}

#[test]
fn peak_rss_is_bytes() {
    if !cfg!(any(target_os = "linux", target_os = "macos")) {
        return;
    }
    let size = 64usize << 20;
    let mut block = vec![0u8; size];
    for i in (0..size).step_by(4096) {
        block[i] = 1;
    }
    let peak = peak_rss_bytes().unwrap();
    assert!(std::hint::black_box(&block).iter().step_by(4096).all(|&b| b == 1));
    assert!(peak >= size as u64, "peak {peak} below 64 MiB");
}

#[test]
fn hasher_reports_sha256_and_size() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("abc.txt");
    std::fs::write(&path, b"abc").unwrap();
    let files = InputHasher::start(vec![("database".to_string(), path)]).join();
    assert_eq!(files.len(), 1);
    assert_eq!(files[0].bytes, Some(3));
    assert_eq!(
        files[0].sha256.as_deref(),
        Ok("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
    );
}

#[test]
fn hasher_reports_a_missing_file_without_failing() {
    let dir = tempfile::tempdir().unwrap();
    let present = dir.path().join("abc.txt");
    std::fs::write(&present, b"abc").unwrap();
    let hasher = InputHasher::start(vec![
        ("spectrum".to_string(), dir.path().join("missing.mzML")),
        ("spectrum".to_string(), present),
    ]);
    let other_work: usize = (0..1000).sum();
    assert_eq!(other_work, 499_500);
    let files = hasher.join();
    assert_eq!(files[0].bytes, None);
    assert!(files[0].sha256.is_err());
    assert!(files[1].sha256.is_ok());
}

#[test]
fn command_line_quotes_arguments_with_spaces() {
    let args = ["andes", "--spectrum", "my run.mzML", "--rescore"].map(String::from);
    assert_eq!(
        command_line_of(args.into_iter()),
        "andes --spectrum 'my run.mzML' --rescore"
    );
}
