use std::path::{Path, PathBuf};
use std::process::Command;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("canonicalize workspace root")
}

fn is_float_field(v: &str) -> bool {
    v.parse::<f64>().is_ok() && v.parse::<i64>().is_err()
}

fn float_columns(rows: &[String]) -> Vec<bool> {
    let mut float = Vec::new();
    for row in rows {
        for (i, v) in row.split('\t').enumerate() {
            if i >= float.len() {
                float.resize(i + 1, false);
            }
            float[i] |= is_float_field(v);
        }
    }
    float
}

fn identity_digest(header: &str, rows: &[String]) -> (usize, u64) {
    let float = float_columns(rows);
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |bytes: &[u8]| {
        for &b in bytes {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    feed(header.as_bytes());
    feed(b"\n");
    for row in rows {
        for (i, v) in row.split('\t').enumerate() {
            if i > 0 {
                feed(b"\t");
            }
            if !float[i] {
                feed(v.as_bytes());
            }
        }
        feed(b"\n");
    }
    (rows.len(), h)
}

fn fields_match(a: &str, b: &str) -> bool {
    if a == b {
        return true;
    }
    if !(is_float_field(a) && is_float_field(b)) {
        return false;
    }
    let (x, y): (f64, f64) = (a.parse().unwrap(), b.parse().unwrap());
    (x - y).abs() <= 1e-5 * x.abs().max(y.abs()) + 1e-6
}

fn assert_matches_golden_pin(header: &str, rows: &[String], golden: &str) {
    use std::io::Read;
    let path = workspace_root().join(golden);
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("open golden PIN"))
        .read_to_string(&mut text)
        .expect("inflate golden PIN");
    let mut lines = text.lines();
    let golden_header = lines.next().expect("golden header");
    let golden_rows: Vec<&str> = lines.collect();
    assert_eq!(header, golden_header, "PIN header differs from {golden}");
    assert_eq!(
        rows.len(),
        golden_rows.len(),
        "PIN row count differs from {golden}"
    );
    let cols: Vec<&str> = header.split('\t').collect();
    let mut diffs: Vec<String> = Vec::new();
    for (row, want) in rows.iter().zip(&golden_rows) {
        let got: Vec<&str> = row.split('\t').collect();
        let exp: Vec<&str> = want.split('\t').collect();
        if got.len() != exp.len() {
            diffs.push(format!(
                "{}: {} fields, golden {}",
                got[0],
                got.len(),
                exp.len()
            ));
            continue;
        }
        for (i, (g, e)) in got.iter().zip(&exp).enumerate() {
            if !fields_match(g, e) {
                let col = cols.get(i).copied().unwrap_or("?");
                diffs.push(format!("{} {col}: {g} (golden {e})", exp[0]));
            }
        }
    }
    assert!(
        diffs.is_empty(),
        "{} field(s) differ from {golden}; first ones:\n{}",
        diffs.len(),
        diffs
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

fn read_pin(path: &Path) -> (String, Vec<String>) {
    let text = std::fs::read_to_string(path).expect("read pin");
    let mut lines = text.lines();
    let header = lines.next().expect("pin header").to_string();
    let mut rows: Vec<String> = lines.map(str::to_string).collect();
    rows.sort();
    (header, rows)
}

fn calibration_outcome(stderr: &str) -> Option<String> {
    stderr
        .lines()
        .find(|l| {
            l.starts_with("Precursor mass shift learned")
                || l.starts_with("Precursor mass calibration skipped (insufficient")
        })
        .map(|l| l.split("; elapsed").next().unwrap_or(l).to_string())
}

fn run(args: &[&str], out: &Path) -> (String, Vec<String>, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_andes"))
        .current_dir(workspace_root())
        .args(args)
        .args(["--predict", "none"])
        .args(["--max-missed-cleavages", "1"])
        .arg("--output-pin")
        .arg(out)
        .output()
        .expect("spawn andes");
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    assert!(
        output.status.success(),
        "andes {} failed:\n{stderr}",
        args.join(" ")
    );
    assert!(
        stderr.contains("candidate-index: out-of-core mmap"),
        "the run did not report the out-of-core index:\n{stderr}"
    );
    let (header, rows) = read_pin(out);
    (header, rows, stderr)
}

type Golden = (usize, u64, Option<&'static str>);

const ECOLI_DEFAULT_MODS: Golden = (75_009, 0xc189_0f50_38d5_c9ea, None);
const ECOLI_FIXED_MODS_ONLY: Golden = (77_496, 0xa181_3557_efb7_ae83, None);
const ECOLI_ONE_VARIABLE_MOD: Golden = (75_570, 0x428b_4490_d8f0_7cbc, None);
const ECOLI_PRECURSOR_CAL_AUTO: Golden = (
    75_009,
    0xc189_0f50_38d5_c9ea,
    Some("Precursor mass calibration skipped (insufficient confident PSMs: 489 with PSMs, 63 below RawScore floor, 0 failed |residual|>50ppm"),
);

fn assert_ecoli_golden(label: &str, mods: &[&str], precursor_cal: &str, golden: Golden) {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin = dir.path().join("out.pin");
    let mut args = vec![
        "--spectrum",
        "test-fixtures/test.mgf.gz",
        "--database",
        "test-fixtures/ecoli.fasta",
        "--fragmentation",
        "HCD",
        "--fragment-tol-ppm",
        "20",
        "--threads",
        "1",
        "--precursor-cal",
        precursor_cal,
        "--cal-min-spec-keys",
        "1000",
        "--fragment-index",
        "off",
    ];
    args.extend(mods);
    let (header, rows, stderr) = run(&args, &pin);
    let cal = calibration_outcome(&stderr);
    if precursor_cal != "off" {
        assert!(
            cal.is_some(),
            "[{label}] the calibration pre-pass did not score:\n{stderr}"
        );
    }
    assert!(
        rows.len() > 10_000,
        "[{label}] only {} rows — the search did not run at the intended scale",
        rows.len()
    );
    let (n, h) = identity_digest(&header, &rows);
    assert_eq!(
        (n, h, cal.as_deref()),
        golden,
        "[{label}] PIN (rows, digest, calibration) differs from the golden"
    );
}

#[test]
#[ignore = "runs a full E. coli search (~20 s release)"]
fn ecoli_pin_matches_golden_with_default_mods() {
    assert_ecoli_golden("default mods", &[], "off", ECOLI_DEFAULT_MODS);
}

#[test]
#[ignore = "runs a full E. coli search (~20 s release)"]
fn ecoli_pin_matches_golden_with_fixed_mods_only() {
    assert_ecoli_golden(
        "fixed mods only",
        &[
            "--fixed-mods",
            "Carbamidomethyl (C)",
            "--variable-mods",
            "",
            "--max-mods",
            "0",
        ],
        "off",
        ECOLI_FIXED_MODS_ONLY,
    );
}

#[test]
#[ignore = "runs a full E. coli search (~20 s release)"]
fn ecoli_pin_matches_golden_with_one_variable_mod() {
    assert_ecoli_golden(
        "one variable mod",
        &[
            "--fixed-mods",
            "Carbamidomethyl (C)",
            "--variable-mods",
            "Oxidation (M)",
            "--max-mods",
            "1",
        ],
        "off",
        ECOLI_ONE_VARIABLE_MOD,
    );
}

#[test]
#[ignore = "runs a full E. coli search (~20 s release)"]
fn ecoli_pin_matches_golden_with_precursor_calibration_on() {
    assert_ecoli_golden("precursor-cal auto", &[], "auto", ECOLI_PRECURSOR_CAL_AUTO);
}

const CHIMERIC_MZML: &str = "test-fixtures/astral_pxd070049_slice.mzML.gz";

fn run_slice(fragment_index: &str, chimeric: bool, out: &Path) -> (String, Vec<String>, String) {
    let mut args = vec![
        "--spectrum",
        CHIMERIC_MZML,
        "--database",
        "test-fixtures/ecoli.fasta",
        "--threads",
        "2",
        "--precursor-cal",
        "off",
        "--fragment-index",
        fragment_index,
    ];
    if chimeric {
        args.push("--chimeric");
    }
    run(&args, out)
}

fn run_chimeric(fragment_index: &str, out: &Path) -> (String, Vec<String>, String) {
    run_slice(fragment_index, true, out)
}

const CHIMERIC_SLICE_PIN: &str = "test-fixtures/goldens/chimeric_slice.pin.gz";

fn secondary_rows(fragment_index: &str, chimeric_rows: &[String], dir: &Path) -> usize {
    let (_, plain_rows, _) = run_slice(
        fragment_index,
        false,
        &dir.join(format!("plain_{fragment_index}.pin")),
    );
    chimeric_rows
        .len()
        .checked_sub(plain_rows.len())
        .expect("a --chimeric run must not emit fewer rows than the plain search")
}

fn duplicate_identifications(header: &str, rows: &[String]) -> Vec<String> {
    let cols: Vec<&str> = header.split('\t').collect();
    let col = |name: &str| cols.iter().position(|c| *c == name).expect(name);
    let (scan, label, pep) = (col("ScanNr"), col("Label"), col("Peptide"));
    let charges: Vec<usize> = (2..=5).map(|z| col(&format!("charge{z}"))).collect();
    let mut seen = std::collections::HashSet::new();
    let mut dups = Vec::new();
    for row in rows {
        let f: Vec<&str> = row.split('\t').collect();
        let charge: Vec<&str> = charges.iter().map(|&i| f[i]).collect();
        let key = (f[scan], f[pep], charge, f[label]);
        if !seen.insert(key) {
            dups.push(format!("scan {} {} label {}", f[scan], f[pep], f[label]));
        }
    }
    dups
}

#[test]
fn chimeric_enumeration_pin_matches_golden() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (header, rows, _) = run_chimeric("off", &dir.path().join("enum.pin"));
    let dups = duplicate_identifications(&header, &rows);
    assert!(dups.is_empty(), "repeated identifications: {dups:?}");
    let secondaries = secondary_rows("off", &rows, dir.path());
    assert!(
        secondaries > 50,
        "fixture produced only {secondaries} chimeric secondary rows — the gate would be vacuous"
    );
    assert_matches_golden_pin(&header, &rows, CHIMERIC_SLICE_PIN);
}

#[test]
fn chimeric_fragment_index_runs_with_unique_specids() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (_, enum_rows, _) = run_chimeric("off", &dir.path().join("enum.pin"));
    let (idx_header, idx_rows, idx_err) = run_chimeric("auto", &dir.path().join("idx.pin"));
    let dups = duplicate_identifications(&idx_header, &idx_rows);
    assert!(dups.is_empty(), "repeated identifications: {dups:?}");
    assert!(
        idx_err.contains("candidate-retrieval: fragment-ion index"),
        "the fragment-ion index did not run:\n{idx_err}"
    );
    let ids: std::collections::HashSet<&str> = idx_rows
        .iter()
        .map(|r| r.split('\t').next().expect("SpecId column"))
        .collect();
    assert_eq!(ids.len(), idx_rows.len(), "duplicate SpecIds in the PIN");
    let enum_sec = secondary_rows("off", &enum_rows, dir.path()) as f64;
    let idx_sec = secondary_rows("auto", &idx_rows, dir.path()) as f64;
    assert!(
        (idx_sec - enum_sec).abs() <= 0.1 * enum_sec,
        "chimeric secondaries: fragment index {idx_sec} vs enumeration {enum_sec}"
    );
}

fn run_refine(fragment_index: &str, out: &Path) -> (String, Vec<String>, String) {
    run(
        &[
            "--spectrum",
            CHIMERIC_MZML,
            "--database",
            "test-fixtures/ecoli.fasta",
            "--threads",
            "2",
            "--precursor-cal",
            "off",
            "--fragment-index",
            fragment_index,
            "--refine",
            "--psm-fdr",
            "1.0",
        ],
        out,
    )
}

const REFINE_SLICE_PIN: &str = "test-fixtures/goldens/refine_slice.pin.gz";

fn refinement_rows(header: &str, rows: &[String]) -> usize {
    let col = header
        .split('\t')
        .position(|c| c == "IsRefinement")
        .expect("IsRefinement column");
    rows.iter()
        .filter(|r| r.split('\t').nth(col) == Some("1"))
        .count()
}

#[test]
fn refine_enumeration_pin_matches_golden() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (header, rows, _) = run_refine("off", &dir.path().join("enum.pin"));
    let refined = refinement_rows(&header, &rows);
    assert!(
        refined > 100,
        "fixture produced only {refined} IsRefinement rows — the gate would be vacuous"
    );
    assert_matches_golden_pin(&header, &rows, REFINE_SLICE_PIN);
}

#[test]
fn refine_on_the_fragment_index_writes_a_well_formed_pin() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (header, rows, err) = run_refine("auto", &dir.path().join("idx.pin"));
    assert!(
        err.contains("candidate-retrieval: fragment-ion index"),
        "--refine did not take the fragment-ion index:\n{err}"
    );
    assert!(
        refinement_rows(&header, &rows) > 100,
        "too few IsRefinement rows on the index path:\n{err}"
    );
    let cols: Vec<&str> = header.split('\t').collect();
    let scan_col = cols.iter().position(|c| *c == "ScanNr").expect("ScanNr");
    let pep_col = cols.iter().position(|c| *c == "Peptide").expect("Peptide");
    let ids: std::collections::HashSet<&str> = rows
        .iter()
        .map(|r| r.split('\t').next().expect("SpecId column"))
        .collect();
    assert_eq!(ids.len(), rows.len(), "duplicate SpecIds in the PIN");
    let scan_pep: std::collections::HashSet<(&str, &str)> = rows
        .iter()
        .map(|r| {
            let f: Vec<&str> = r.split('\t').collect();
            (f[scan_col], f[pep_col])
        })
        .collect();
    assert_eq!(
        scan_pep.len(),
        rows.len(),
        "a scan carries the same peptide in two rows"
    );
}
