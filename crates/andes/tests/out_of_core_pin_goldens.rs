//! Golden gates for the enumeration retrieval path (`--fragment-index off`) of
//! the out-of-core candidate index, driven through the real CLI.
//!
//! Rows are sorted before comparison because the rayon pipeline's row ORDER
//! varies across runs; their content does not.
//!
//! Floating-point features (RichIonLLR, RankScoreFloat, ...) go through the
//! platform's `exp`/`log`/`pow`, whose last-ulp rounding differs between C
//! libraries (glibc vs macOS libm), and a 1-ulp change can move the sixth
//! printed significant digit. So no gate hashes float columns bit for bit:
//!
//! - the `--chimeric` and `--refine` slice cases compare against a golden PIN
//!   (`test-fixtures/goldens/*.pin.gz`, header + sorted rows) row by row —
//!   every non-float field exactly, every float field to one unit in its sixth
//!   significant digit (see `fields_match`);
//! - the E. coli cases (too large to commit as PINs) compare a digest of the
//!   non-float columns (`identity_digest`: row count + FNV-1a 64 over the header
//!   and the sorted rows with float fields blanked), and for calibration also
//!   the pre-pass outcome line. The non-float columns hold every row's SpecId,
//!   label, scan, charge, peptide, proteins and integer features, so a change in
//!   the candidate set, its order or its multiplicity still fails the gate.
//!
//! The defects these gates were built for need real scale to surface — a fixed
//! modification on a real proteome (the `base_record_multiplicity` fixed-mod
//! mass offset), proteins with an N-terminal Met whose two enumeration passes
//! both reach the precursor window (the candidate ORDER), a calibration pre-pass
//! that must run on the out-of-core index — so they search the in-repo
//! `test.mgf.gz` + `ecoli.fasta`. Those four cases are `#[ignore]`d (a full
//! E. coli search each, ~20 s release, minutes in debug):
//!
//! ```text
//! cargo test -p andes --release --test out_of_core_pin_goldens -- --ignored
//! ```
//!
//! The `--chimeric` and `--refine` cases run on a 150-scan MS1-bearing slice and
//! are not ignored.
//!
//! GOLDENS (recorded 2026-10-11 on x86-64 Linux / glibc). The E. coli and
//! `--chimeric` goldens were recorded from the in-memory candidate list
//! (whole-database enumeration + a mass-bucket lookup), the path these searches
//! took before the out-of-core index became the only one, and the out-of-core
//! enumeration reproduced each of those PINs byte for byte after sorting. The
//! `--refine` golden was recorded the same way, with the second pass on the
//! in-memory candidate list; searching the anchored database through an
//! in-memory out-of-core index reproduced it (and the Astral benchmark PIN with
//! `--refine`) byte for byte.
//!
//! REGENERATED 2026-10-13 when the writers began collapsing positional isomers
//! (one row per peptide and modification multiset per spectrum precursor and
//! charge). Each new golden was checked against the previous binary's PIN for
//! the same command line with its isomer alternatives removed (first row in
//! emission order kept): identical row for row with SpecId excluded, SpecIds
//! unique, every spectrum's first row unchanged. Rows dropped, all isomer
//! alternatives (Ox-M sites, protein-N-terminal Acetyl vs internal K, and on the
//! refine slice deamidation sites, plus 24 rows repeating a kept modified
//! peptide with different flanking residues):
//!
//! | case | rows before | rows after |
//! |---|---|---|
//! | E. coli default mods (and precursor-cal auto) | 77,097 | 74,442 |
//! | E. coli fixed mods only | 75,750 | 75,750 |
//! | E. coli one variable mod | 77,110 | 74,998 |
//! | `--chimeric` slice | 1,751 | 1,717 |
//! | `--refine` slice | 1,764 | 1,595 |
//!
//! REGENERATED 2026-10-13 when the positional-isomer key began separating
//! target from decoy and a kept row began listing the proteins of the rows
//! collapsed into it. Against the previous binary's PIN for the same command
//! line: same rows and SpecIds, every spectrum's winner unchanged, every
//! column but `Proteins` identical; `Proteins` gains accessions on the rows
//! whose collapsed copies came from another protein (same residues, different
//! flanking residues). No target/decoy pair had collided, so no row returned.
//!
//! | case | rows | rows with new proteins |
//! |---|---|---|
//! | E. coli default mods (and precursor-cal auto) | 74,442 | 2 |
//! | E. coli fixed mods only | 75,750 | 0 (golden unchanged) |
//! | E. coli one variable mod | 74,998 | 4 |
//! | `--chimeric` slice | 1,717 | 0 (golden unchanged) |
//! | `--refine` slice | 1,595 | 24 |
//!
//! REGENERATED 2026-10-14 for the `--chimeric` slice after merging main, whose
//! Pass 2 no longer emits a secondary that repeats a retained candidate (main's
//! golden: 1,749 rows). The new golden is main's with its positional-isomer
//! alternatives removed: 34 rows dropped, the remaining 1,715 identical row for
//! row with SpecId excluded (no row gains a protein on this slice), every
//! spectrum's winner unchanged. Against this branch's previous golden (1,717
//! rows) it lacks exactly main's two dropped Pass-2 rows (scans 60160 and
//! 60179); every other row is identical with SpecId excluded.
//!
//! REGENERATED 2026-10-14 (E. coli digests) when the MGF reader began taking
//! the scan number from the TITLE (`... scan=N`) instead of reporting 0. Same
//! row counts; with ScanNr set back to 0 and SpecId's scan field back to `_0_`,
//! each new PIN digests to its previous golden exactly, so only those two
//! fields changed.
//!
//! REGENERATED 2026-10-14 (E. coli digests) when the PIN charge one-hot columns
//! began spanning every charge a PSM carries, not only the `--charge` range.
//! Against the PIN without that change, for the same command line: same rows and
//! SpecIds, every non-charge field identical; the header gains `charge6`, set
//! on the rows of spectra with a declared 6+ precursor, which previously had an
//! all-zero charge vector. The slice goldens carry no out-of-range charge and
//! are unchanged.
//!
//! | case | rows | rows moved to `charge6` |
//! |---|---|---|
//! | E. coli default mods (and precursor-cal auto) | 74,442 | 37 |
//! | E. coli fixed mods only | 75,750 | 35 |
//! | E. coli one variable mod | 74,998 | 40 |
//!
//! The E. coli digests above carry both of these changes together (TITLE scan
//! numbers and the widened charge columns); row counts are unchanged.
//!
//! REGENERATED 2026-10-14 for the `--chimeric` slice when Pass 2 began reporting
//! co-isolated envelopes at their monoisotopic peak, ranking secondaries over
//! distinct peptidoforms and filling their features on the primaries' scales.
//! Pass 1 is untouched: the same command line without `--chimeric` writes a
//! byte-identical PIN before and after. Against the previous golden (1,715 rows)
//! the new one has 1,716: all 1,449 primary rows keep their peptide, label and
//! charge; the Pass-2 rows go from 266 to 267, 243 old ones replaced by 244 new
//! ones (mostly secondaries searched at M rather than M+1). DeltaRT, AbsDeltaRT and
//! DeltaRTNorm move on every row because the RT fit includes the secondary rows; RawScore, RankScore, RankScoreFloat,
//! RawScoreCal, CandidateRankEntropy and ListwiseScoreGap move on ~20 matched
//! rows. The E. coli and `--refine` goldens are unchanged.
//!
//! REGENERATED 2026-10-14 for two candidate-generation changes: the
//! N-terminal Met-cleaved walk emits only the peptides starting at the
//! post-Met residue (main enumerated every other span of an M-initial protein
//! twice), and a residue's side-chain mods stack with its terminal mods
//! (protein N-terminal Acetyl with Ox-M, e.g. `_.M+58.00548IVLPRR`). Every
//! golden was recorded from this branch's binary with the test's command
//! lines; main's binary reproduces the previous goldens byte for byte. Against
//! main's PIN, rows keyed by spectrum, charge, peptide and label; "top row" is
//! the highest-RawScore row of a spectrum and charge:
//!
//! | case | rows main → branch | only main / only branch | top row changed |
//! |---|---|---|---|
//! | E. coli default mods (and precursor-cal auto) | 74,442 → 75,638 | 7,105 / 8,301 | 550 of 5,574 |
//! | E. coli fixed mods only | 75,750 → 77,975 | 7,198 / 9,423 | 534 of 5,562 |
//! | E. coli one variable mod | 74,998 → 76,087 | 7,109 / 8,198 | 535 of 5,568 |
//! | `--chimeric` slice | 1,715 → 1,717 | 131 / 133 | 19 of 250 |
//! | `--refine` slice | 1,595 → 1,612 | 138 / 155 | 13 of 150 |
//!
//! Each spectrum's candidate pool loses the duplicate copies, so the rows
//! below the top fill with other peptides and the pool statistics move on
//! shared rows: `CandidateRankEntropy`, `ListwiseScoreGap`, `DeltaRankScore`,
//! `TailorScore`, and `RawScore` / `RawScoreCal`, which fuse them; on the
//! slices also the `DeltaRT` columns, whose RT fit follows the top rows.
//! `Proteins` no longer repeats an M-initial protein's accession. The
//! precursor-cal pre-pass outcome line is unchanged (492 with PSMs, 55 below
//! the RawScore floor).
//!
//! To regenerate after an intentional change, run the case's command line (see
//! `run`), then either digest the PIN with `identity_digest` (E. coli) or write
//! its header and sorted rows, gzipped, over the case's golden PIN (slices).

use std::path::{Path, PathBuf};
use std::process::Command;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("canonicalize workspace root")
}

/// A field printed as a non-integer number (`0.3`, `1e-07`, `inf`, `NaN`).
fn is_float_field(v: &str) -> bool {
    v.parse::<f64>().is_ok() && v.parse::<i64>().is_err()
}

/// Columns holding a non-integer number in any row.
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

/// Row count and FNV-1a 64 of the header and the sorted rows with every float
/// column's fields blanked, each line followed by a newline.
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

/// Two PIN fields agree: identical text, or two floats within one unit of the
/// sixth significant digit (`%g`-style output) or of the sixth decimal
/// (fixed-point output) — the most a last-ulp difference in the platform's
/// math library moves a printed value.
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

/// Compare a PIN (header + sorted rows) with a gzipped golden PIN holding the
/// header and the sorted rows; panics listing the first differing fields.
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

/// A PIN's header and its rows, sorted.
fn read_pin(path: &Path) -> (String, Vec<String>) {
    let text = std::fs::read_to_string(path).expect("read pin");
    let mut lines = text.lines();
    let header = lines.next().expect("pin header").to_string();
    let mut rows: Vec<String> = lines.map(str::to_string).collect();
    rows.sort();
    (header, rows)
}

/// The pre-pass's outcome line with its timing stripped, or `None` if the
/// pre-pass never scored (threshold skip / calibration off). Both a learned
/// shift and an "insufficient confident PSMs" report carry the PSM counts the
/// pre-pass resolved through the candidate index, so equality of this line
/// with the golden is the calibration half of the gate.
fn calibration_outcome(stderr: &str) -> Option<String> {
    stderr
        .lines()
        .find(|l| {
            l.starts_with("Precursor mass shift learned")
                || l.starts_with("Precursor mass calibration skipped (insufficient")
        })
        .map(|l| l.split("; elapsed").next().unwrap_or(l).to_string())
}

/// Run `andes` from the workspace root with `args` and an output PIN; returns the
/// header, the sorted rows and stderr.
fn run(args: &[&str], out: &Path) -> (String, Vec<String>, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_andes"))
        .current_dir(workspace_root())
        .args(args)
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
        stderr.contains("out-of-core candidate-index: mmap"),
        "the run did not report the out-of-core index:\n{stderr}"
    );
    let (header, rows) = read_pin(out);
    (header, rows, stderr)
}

/// One E. coli case: (rows, `identity_digest` FNV-1a 64, calibration outcome line).
type Golden = (usize, u64, Option<&'static str>);

const ECOLI_DEFAULT_MODS: Golden = (75_638, 0x21bf_e101_3492_e5d1, None);
const ECOLI_FIXED_MODS_ONLY: Golden = (77_975, 0x8c8a_554f_a7ac_6d9c, None);
const ECOLI_ONE_VARIABLE_MOD: Golden = (76_087, 0xeebe_f080_9369_9e90, None);
/// The pre-pass scores but the sample is too small to calibrate, so the PIN
/// equals the calibration-off one.
const ECOLI_PRECURSOR_CAL_AUTO: Golden = (
    75_638,
    0x21bf_e101_3492_e5d1,
    Some("Precursor mass calibration skipped (insufficient confident PSMs: 492 with PSMs, 55 below RawScore floor, 0 failed |residual|>50ppm"),
);

fn assert_ecoli_golden(label: &str, mods: Option<&str>, precursor_cal: &str, golden: Golden) {
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
        // The fixture has ~5.8k SpecKeys, below the 10k production threshold;
        // force the pre-pass so the calibration-on case actually calibrates.
        "--cal-min-spec-keys",
        "1000",
        // Enumeration: every candidate in each precursor window is scored.
        "--fragment-index",
        "off",
    ];
    if let Some(m) = mods {
        args.extend(["--mods", m]);
    }
    let (header, rows, stderr) = run(&args, &pin);
    let cal = calibration_outcome(&stderr);
    if precursor_cal != "off" {
        // The pre-pass must have SCORED (either it learned a shift, or it
        // reports how many sampled spectra had PSMs). A threshold skip
        // ("SpecKeys < ...") means it never ran and the gate would be vacuous.
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
    assert_ecoli_golden("default mods", None, "off", ECOLI_DEFAULT_MODS);
}

#[test]
#[ignore = "runs a full E. coli search (~20 s release)"]
fn ecoli_pin_matches_golden_with_fixed_mods_only() {
    // NumMods=0 with a fixed Cam-C: isolates the `base_record_multiplicity`
    // fixed-mod defect from any variable-mod expansion.
    let dir = tempfile::tempdir().expect("tempdir");
    let mods = dir.path().join("mods_fixed.txt");
    std::fs::write(&mods, "NumMods=0\n57.021464,C,fix,any,Carbamidomethyl\n").expect("write mods");
    assert_ecoli_golden(
        "fixed mods only",
        Some(mods.to_str().unwrap()),
        "off",
        ECOLI_FIXED_MODS_ONLY,
    );
}

#[test]
#[ignore = "runs a full E. coli search (~20 s release)"]
fn ecoli_pin_matches_golden_with_one_variable_mod() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mods = dir.path().join("mods_var1.txt");
    std::fs::write(
        &mods,
        "NumMods=1\n57.021464,C,fix,any,Carbamidomethyl\n15.994915,M,opt,any,Oxidation\n",
    )
    .expect("write mods");
    assert_ecoli_golden(
        "one variable mod",
        Some(mods.to_str().unwrap()),
        "off",
        ECOLI_ONE_VARIABLE_MOD,
    );
}

#[test]
#[ignore = "runs a full E. coli search (~20 s release)"]
fn ecoli_pin_matches_golden_with_precursor_calibration_on() {
    // The calibration pre-pass runs on the out-of-core index.
    assert_ecoli_golden("precursor-cal auto", None, "auto", ECOLI_PRECURSOR_CAL_AUTO);
}

// ── --chimeric ──────────────────────────────────────────────────────────────
//
// Fixture: `test-fixtures/astral_pxd070049_slice.mzML.gz`, 150 consecutive MS2
// scans with their two linked MS1 scans, cut unmodified from the middle of the
// PRIDE PXD070049 run `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw`
// (Orbitrap Astral DDA, human/yeast/E. coli; converted with ThermoRawFileParser).
// The chimeric Pass 2 needs the MS1 scans, which `test.mgf.gz` lacks.

const CHIMERIC_MZML: &str = "test-fixtures/astral_pxd070049_slice.mzML.gz";

/// Run a search on the chimeric fixture, with or without `--chimeric`;
/// returns (header, sorted rows, stderr).
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

/// Run a `--chimeric` search; returns (header, sorted rows, stderr).
fn run_chimeric(fragment_index: &str, out: &Path) -> (String, Vec<String>, String) {
    run_slice(fragment_index, true, out)
}

/// `--chimeric` with enumeration retrieval: header + sorted rows.
const CHIMERIC_SLICE_PIN: &str = "test-fixtures/goldens/chimeric_slice.pin.gz";

/// Chimeric secondaries: the rows a `--chimeric` run adds over the same search
/// without it. Pass 1 of the cascade is that search, so the difference is the
/// rows Pass 2 contributed.
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

/// Rows that repeat a (scan, peptidoform, charge, label) already emitted: the
/// chimeric Pass 2 must not re-emit a peptide the scan's retained Pass-1
/// candidates already carry (it would sit in its own competition group).
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

/// Fragment-index retrieval changes Pass 1's candidates by design, so only the
/// shape is gated: the run completes, SpecIds stay unique, and the number of
/// chimeric secondaries stays close to the enumeration run.
#[test]
fn chimeric_fragment_index_runs_with_unique_specids() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (_, enum_rows, _) = run_chimeric("off", &dir.path().join("enum.pin"));
    let (idx_header, idx_rows, idx_err) = run_chimeric("auto", &dir.path().join("idx.pin"));
    let dups = duplicate_identifications(&idx_header, &idx_rows);
    assert!(dups.is_empty(), "repeated identifications: {dups:?}");
    assert!(
        idx_err.contains("candidate retrieval: fragment-ion index"),
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

// ── --refine ────────────────────────────────────────────────────────────────
//
// The refinement cascade runs after the scan, over the candidates the Pass-1
// PSMs reference (the materialized candidates synced into
// `prepared.candidates`); Pass 2 searches an in-memory out-of-core index of the
// anchored database. Enumeration must reproduce the golden unified PIN, and
// fragment-index retrieval must give a well-formed one. The Astral slice searched against E. coli alone has
// a single anchor at the default 1% gate; `--refine-select-psm-fdr 1.0` anchors
// every best-target peptide so Pass 2 contributes a few hundred rows.

/// Run a `--refine` search; returns (header, sorted rows, stderr).
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
            "--refine-select-psm-fdr",
            "1.0",
        ],
        out,
    )
}

/// `--refine` with enumeration retrieval: header + sorted rows.
const REFINE_SLICE_PIN: &str = "test-fixtures/goldens/refine_slice.pin.gz";

/// Number of rows with `IsRefinement = 1`.
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

/// Fragment-index retrieval changes Pass 1's candidates by design, so only the
/// shape is gated: the index runs under `--refine`, the cascade adds rows,
/// SpecIds stay unique and no scan carries the same peptide twice.
#[test]
fn refine_on_the_fragment_index_writes_a_well_formed_pin() {
    let dir = tempfile::tempdir().expect("tempdir");
    let (header, rows, err) = run_refine("auto", &dir.path().join("idx.pin"));
    assert!(
        err.contains("candidate retrieval: fragment-ion index"),
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
