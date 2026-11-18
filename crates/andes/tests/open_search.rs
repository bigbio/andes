use std::path::PathBuf;
use std::process::Command;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn run(out_dir: &std::path::Path, extra: &[&str]) -> (String, String, String) {
    let tsv = out_dir.join("out.tsv");
    let pin = out_dir.join("out.pin");
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_andes"));
    cmd.args([
        "--spectrum",
        fixture("test-fixtures/orbitrap_lumos_120.mzML.gz")
            .to_str()
            .unwrap(),
        "--database",
        fixture("test-fixtures/glyco_fixture.fasta")
            .to_str()
            .unwrap(),
        "--threads",
        "2",
        "--output-pin",
        pin.to_str().unwrap(),
        "--output-tsv",
        tsv.to_str().unwrap(),
    ]);
    cmd.args(extra);
    let out = cmd.output().expect("run andes");
    assert!(
        out.status.success(),
        "andes failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let tsv_text = std::fs::read_to_string(&tsv).expect("read tsv");
    let pin_text = std::fs::read_to_string(&pin).expect("read pin");
    (
        String::from_utf8_lossy(&out.stderr).to_string(),
        tsv_text,
        pin_text,
    )
}

#[test]
fn open_search_reports_the_offset_inside_the_window() {
    let tmp = tempfile::tempdir().unwrap();
    let (log, tsv, pin) = run(tmp.path(), &["--open-search", "--open-search-delta=-1..1"]);

    assert!(
        log.contains("open-search: mass-offset window -1..+1 Da"),
        "the run line names the window: {log}"
    );
    assert!(
        log.contains("open-search-index:"),
        "the database-wide index is built: {log}"
    );

    let mut lines = tsv.lines();
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let delta_col = header
        .iter()
        .position(|c| *c == "DeltaMass")
        .expect("open search appends a DeltaMass column");
    let site_col = header
        .iter()
        .position(|c| *c == "OpenSite")
        .expect("open search appends an OpenSite column");
    let margin_col = header
        .iter()
        .position(|c| *c == "OpenLocMargin")
        .expect("open search appends an OpenLocMargin column");

    let mut rows = 0usize;
    let mut localized = 0usize;
    for line in lines.filter(|l| !l.is_empty()) {
        let fields: Vec<&str> = line.split('\t').collect();
        let delta: f64 = fields[delta_col].parse().expect("numeric DeltaMass");
        assert!(
            delta.abs() <= 1.12,
            "offset {delta} outside the requested -1..1 window (+ refinement)"
        );
        let site: u16 = fields[site_col].parse().expect("integer OpenSite");
        let margin: f64 = fields[margin_col].parse().expect("numeric OpenLocMargin");
        if site != 0 {
            assert!((1..=60).contains(&site), "site {site} out of range");
            assert!(margin >= 0.02, "called site with margin {margin}");
            localized += 1;
        }
        rows += 1;
    }
    assert!(rows > 0, "the open search found no PSM at all");
    assert!(localized > 0, "no PSM localised its offset");

    let pin_header: Vec<&str> = pin
        .lines()
        .next()
        .expect("pin header")
        .split('\t')
        .collect();
    for name in [
        "DeltaMass",
        "AbsDeltaMass",
        "OpenSite",
        "OpenLocScore",
        "OpenLocMargin",
    ] {
        assert!(
            pin_header.contains(&name),
            "open search appends the {name} PIN feature"
        );
    }
}

#[test]
fn narrow_search_keeps_its_tsv_schema() {
    let tmp = tempfile::tempdir().unwrap();
    let (_log, tsv, pin) = run(tmp.path(), &[]);
    let header = tsv.lines().next().expect("header");
    assert!(
        !header.split('\t').any(|c| c == "DeltaMass"),
        "the narrow search must not add the open-search column"
    );
    let pin_header = pin.lines().next().expect("pin header");
    assert!(
        !pin_header.split('\t').any(|c| c == "DeltaMass"),
        "the narrow PIN must not add the open-search features"
    );
}

#[test]
fn optimize_runs_and_keeps_the_open_search_window() {
    let tmp = tempfile::tempdir().unwrap();
    let (log, _tsv, _pin) = run(
        tmp.path(),
        &["--open-search", "--optimize", "--optimize-sample", "60"],
    );
    assert!(log.contains("optimize: chosen"), "the grid runs: {log}");
    assert!(
        log.contains("chosen --precursor-tol resolved"),
        "the offset window is not searched: {log}"
    );
    assert!(
        log.contains("mass-offset window -100..+300"),
        "the window survives: {log}"
    );
}

#[test]
fn optimize_grids_the_narrow_search() {
    let tmp = tempfile::tempdir().unwrap();
    let (log, _tsv, _pin) = run(tmp.path(), &["--optimize", "--optimize-sample", "60"]);
    assert!(
        log.contains("optimize: chosen --precursor-tol"),
        "the grid picks a tolerance: {log}"
    );
    assert!(log.contains("stage"), "the grid table is printed: {log}");
}

#[test]
fn open_search_rejects_incompatible_flags() {
    let tmp = tempfile::tempdir().unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args([
            "--spectrum",
            fixture("test-fixtures/orbitrap_lumos_120.mzML.gz")
                .to_str()
                .unwrap(),
            "--database",
            fixture("test-fixtures/glyco_fixture.fasta")
                .to_str()
                .unwrap(),
            "--open-search",
            "--isotope-error",
            "-1..2",
            "--output-pin",
            tmp.path().join("x.pin").to_str().unwrap(),
        ])
        .output()
        .expect("run andes");
    assert!(
        !out.status.success(),
        "an explicit isotope range is refused"
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("chooses the isotope offset"),
        "the message explains why: {err}"
    );
}
