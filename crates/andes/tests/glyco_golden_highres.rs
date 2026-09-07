//! HIGH-RES glyco serve-path golden.
//!
//! REGENERATED 2026-09-03 when the peptide-first RETRIEVAL window became the glyco
//! ppm tolerance on high-resolution MS2 (it had inherited the rank model's 0.5 Da).
//! Verified before copying: same 120 scans, all 120 winning peptides identical, 0
//! decoy winners before and after, header identical; only the five per-spectrum
//! competition statistics moved (RawScore, TailorScore, DeltaRankScore,
//! CandidateRankEntropy, ListwiseScoreGap), which is the signature of a smaller
//! candidate pool and nothing else. Measured on the benchmarks: 7x faster with
//! identifications neutral (mouse 3198 vs 3183 correct; plasma 399 vs 380).
//!
//! Companion to `glyco_golden.rs`, which guards the LOW-RES path. This one exists
//! because an MGF fixture *cannot* reach the high-res branch: selecting a high-res
//! model for MGF requires `--fragment-tol-ppm`, which calls
//! `set_fragment_tol_override` and replaces `mme` -- the very window we want to
//! guard. So the flag needed to get there is the flag that neutralises the guard.
//!
//! An mzML fixture carries analyzer metadata, so the instrument is auto-detected and
//! the override is never involved. That makes this the test that actually pins the
//! window commit 539a3857 ("fix(train)") changed while asserting "Serving is untouched".
//!
//! It guards two things at once:
//!   1. the high-res glyco serve path, as a pinned output (the PIN comparison), and
//!   2. that instrument detection survives a `.gz` input (the model assertion) --
//!      a gzipped mzML used to skip detection entirely and silently fall back to
//!      `cid_lowres_tryp`, i.e. the low-res model on high-res Orbitrap data.
//!
//! MEASURED LIMIT, so nobody assumes more of this test than it delivers: on this
//! fixture the high-res configuration is NOT sensitive to the fused-selector
//! weights -- sweeping a `--glyco-gp-*` weight produces 0 differing lines, where
//! the same sweep moves all 120 rows on the low-res path in `glyco_golden.rs`. The 20 ppm
//! window leaves too little candidate competition for the weight to change a
//! winner. So selector regressions are caught by the LOW-RES golden; this one
//! catches model/detection regressions and any gross change to the emitted rows.
//!
//! Regenerating is a deliberate act: diff the columns first and know which ones moved.
//!
//!   ./target/release/andes \
//!     --spectrum test-fixtures/orbitrap_lumos_120.mzML.gz \
//!     --database test-fixtures/glyco_fixture.fasta \
//!     --glyco --glyco-tol-ppm 20 --glyco-taxon human \
//!     --output-pin <tmp>/out.pin
//!   cp <tmp>/out.glyco.pin test-fixtures/parity/goldens/glyco_highres.pin

use std::path::PathBuf;
use std::process::Command;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("repo root")
        .to_path_buf()
}

/// Compare one PIN row field-by-field.
///
/// Identity and integer fields must match EXACTLY; floating-point fields are
/// compared with a relative tolerance. An optimised build vectorises differently
/// on different targets (this golden passed on Linux and Windows and failed on
/// macOS in release), so a last-digit difference in a derived score is a property
/// of the host, not a regression. Anything a real change moves -- a different
/// peptide, a different glycan, a different winner -- is either a non-numeric
/// field or a numeric one that moves far more than 1e-6.
fn row_diff(golden: &str, actual: &str) -> Option<String> {
    let g: Vec<&str> = golden.split('\t').collect();
    let a: Vec<&str> = actual.split('\t').collect();
    if g.len() != a.len() {
        return Some(format!("field count {} vs {}", g.len(), a.len()));
    }
    for (i, (gf, af)) in g.iter().zip(a.iter()).enumerate() {
        if gf == af {
            continue;
        }
        match (gf.parse::<f64>(), af.parse::<f64>()) {
            (Ok(gv), Ok(av)) => {
                // `f64::from_str` accepts "NaN"/"inf", and every ordered comparison
                // against NaN is false -- so a naive `(gv - av).abs() > tol` would
                // report a NaN where a number belongs as EQUAL and wave a real
                // regression through. Handle the non-finite cases explicitly.
                if gv.is_nan() || av.is_nan() {
                    if gv.is_nan() != av.is_nan() {
                        return Some(format!("field {i}: {gf} vs {af} (NaN mismatch)"));
                    }
                    continue;
                }
                if gv.is_infinite() || av.is_infinite() {
                    if gv != av {
                        return Some(format!("field {i}: {gf} vs {af} (infinity mismatch)"));
                    }
                    continue;
                }
                // The PIN is TEXT with ~6 significant figures, so the comparison
                // cannot be tighter than the print precision: two values that agree
                // to within float noise still print differently when they straddle a
                // rounding boundary (macOS CI produced `0.45316` vs `0.453159`, a
                // relative difference of 2.2e-6, from arithmetic that is otherwise
                // identical). 1e-5 relative clears that with room to spare and is
                // still five orders of magnitude tighter than any real change --
                // flipping a selector weight moves this same column from -2.02 to
                // -3.40. The absolute floor keeps near-zero values from being held
                // to denormal noise.
                let tol = 1e-5 * gv.abs().max(av.abs()).max(1.0);
                if (gv - av).abs() > tol {
                    return Some(format!("field {i}: {gf} vs {af}"));
                }
            }
            _ => return Some(format!("field {i}: {gf} vs {af}")),
        }
    }
    None
}

/// Sort key for aligning rows: the SpecId column, which is scan-derived and
/// stable. Sorting by the whole line would let a last-digit float difference
/// reorder rows and produce a spurious mismatch.
fn spec_id(row: &str) -> &str {
    row.split('\t').next().unwrap_or(row)
}

#[test]
fn glyco_highres_pin_matches_golden() {
    let root = repo_root();
    let spectra = root.join("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    let golden = root.join("test-fixtures/parity/goldens/glyco_highres.pin");
    for f in [&spectra, &fasta, &golden] {
        assert!(f.exists(), "fixture missing: {}", f.display());
    }

    // A per-process temp dir, not a fixed shared name: two concurrent runs (two
    // `cargo test` invocations, two CI jobs sharing /tmp, or a retry overlapping
    // its predecessor) would otherwise delete each other's output and fail on a
    // missing PIN rather than on a real mismatch.
    let tmpdir = tempfile::tempdir().expect("tempdir");
    let tmp = tmpdir.path();
    let out = tmp.join("out.pin");

    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    // Deliberately NO --fragmentation / --fragment-tol-ppm: the point is that the
    // instrument is auto-detected from the mzML, which is what puts the run on the
    // high-res branch. Taxon is pinned because the fixture FASTA carries an E. coli
    // background whose OX= tags would otherwise resolve the run to CmahCompetent.
    let result = Command::new(&binary)
        .arg("--spectrum")
        .arg(&spectra)
        .arg("--database")
        .arg(&fasta)
        .arg("--glyco")
        .arg("--glyco-tol-ppm")
        .arg("20")
        .arg("--glyco-taxon")
        .arg("human")
        .arg("--output-pin")
        .arg(&out)
        .output()
        .expect("run andes");
    assert!(
        result.status.success(),
        "andes exited {}\n{}",
        result.status,
        String::from_utf8_lossy(&result.stderr)
    );

    // Detection must survive the .gz. If this regresses, the run silently drops to
    // cid_lowres_tryp and every downstream number changes for a reason no output
    // explains -- which is exactly how this was missed before.
    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        log.contains("hcd_qexactive_tryp"),
        "expected the high-res model to be auto-detected from the gzipped mzML, \
         but the run did not select hcd_qexactive_tryp. Log:\n{log}"
    );

    let actual = tmp.join("out.glyco.pin");
    assert!(
        actual.exists(),
        "glyco PIN not written: {}",
        actual.display()
    );

    let g_txt = std::fs::read_to_string(&golden).expect("read golden");
    let a_txt = std::fs::read_to_string(&actual).expect("read actual");
    let mut g_lines: Vec<&str> = g_txt.lines().collect();
    let mut a_lines: Vec<&str> = a_txt.lines().collect();

    let g_hdr = g_lines.remove(0);
    let a_hdr = a_lines.remove(0);
    assert_eq!(g_hdr, a_hdr, "glyco high-res PIN column set changed");
    assert_eq!(
        g_lines.len(),
        a_lines.len(),
        "glyco high-res PIN row count changed: golden {} vs actual {}",
        g_lines.len(),
        a_lines.len()
    );

    g_lines.sort_unstable_by_key(|r| spec_id(r));
    a_lines.sort_unstable_by_key(|r| spec_id(r));
    for (i, (g, a)) in g_lines.iter().zip(a_lines.iter()).enumerate() {
        if let Some(d) = row_diff(g, a) {
            panic!(
                "glyco PIN row mismatch at sorted index {i} ({d})
 golden: {g}
 actual: {a}"
            );
        }
    }
}

/// `--precursor-mono auto` on an mzML with NO MS1 scans (this fixture declares
/// "MS1 spectrum" only in its file-content header): the run streams with the MS1
/// link, fits nothing, says so, and writes a `.glyco.pin` byte-identical to the
/// golden — no Mono columns — while the diagnostic dump still has one NA row
/// per MS2. (The linked case is covered by the `mono` unit tests and the
/// measured pGlyco2 liver run in docs/benchmarks/README.md.)
#[test]
fn glyco_highres_precursor_mono_without_ms1_is_identical_to_golden() {
    let root = repo_root();
    let spectra = root.join("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    let golden = root.join("test-fixtures/parity/goldens/glyco_highres.pin");
    let tmpdir = tempfile::tempdir().expect("tempdir");
    let out = tmpdir.path().join("out.pin");
    let dump = tmpdir.path().join("mono.tsv");
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let result = Command::new(&binary)
        .arg("--spectrum")
        .arg(&spectra)
        .arg("--database")
        .arg(&fasta)
        .args(["--glyco", "--glyco-tol-ppm", "20", "--glyco-taxon", "human"])
        .args(["--precursor-mono", "auto", "--precursor-mono-dump"])
        .arg(&dump)
        .arg("--output-pin")
        .arg(&out)
        .output()
        .expect("run andes");
    let log = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(
        result.status.success(),
        "andes exited {}\n{log}",
        result.status
    );
    assert!(
        log.contains("precursor-mono: streamed 120 MS2 spectra (0 MS1 scans linked)"),
        "expected the MS1-linked stream summary, got:\n{log}"
    );
    assert!(
        log.contains("no MS2 could be fitted"),
        "expected the nothing-fitted warning, got:\n{log}"
    );

    let actual = std::fs::read_to_string(tmpdir.path().join("out.glyco.pin")).expect("read pin");
    let golden_txt = std::fs::read_to_string(&golden).expect("read golden");
    let a_hdr = actual.lines().next().unwrap();
    assert!(
        !a_hdr.contains("MonoShift"),
        "no Mono columns without a fitted MS2: {a_hdr}"
    );
    // Same rows as the golden (order-insensitive, like the golden test above).
    let mut a_rows: Vec<&str> = actual.lines().skip(1).collect();
    let mut g_rows: Vec<&str> = golden_txt.lines().skip(1).collect();
    a_rows.sort_by_key(|r| spec_id(r).to_string());
    g_rows.sort_by_key(|r| spec_id(r).to_string());
    assert_eq!(a_hdr, golden_txt.lines().next().unwrap());
    assert_eq!(
        a_rows.len(),
        g_rows.len(),
        "row count differs from the golden"
    );
    for (i, (g, a)) in g_rows.iter().zip(a_rows.iter()).enumerate() {
        if let Some(d) = row_diff(g, a) {
            panic!("row {i} differs from the golden ({d})\n golden: {g}\n actual: {a}");
        }
    }

    let dump_txt = std::fs::read_to_string(&dump).expect("read dump");
    assert_eq!(dump_txt.lines().count(), 121, "header + one row per MS2");
    assert!(
        dump_txt.lines().skip(1).all(|l| l.ends_with("\tNA")),
        "no MS2 has a linked MS1, so every row is NA"
    );
}
