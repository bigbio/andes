//! Golden test for the `--glyco` path.
//!
//! WHY THIS EXISTS. Until this file was added, the ONLY golden in the repo covered
//! the standard peptide path (`precursor_cal_bit_identical.rs`). Nothing pinned glyco
//! output, and a commit labelled `fix(train)` consequently changed SERVE-TIME glyco
//! scoring - widening the per-ion match window from 20 ppm to the model's `mme`
//! (0.5 Da, ~36x at m/z 700) inside `ion_match_facts`, which feeds the glyco collapse
//! selector through `hyperscore_psm_with_matches` - and survived 24 commits with a
//! comment asserting "Serving is untouched".
//!
//! THE FIXTURE. 120 MS2 scans lifted from human plasma PXD030622 replicate R1, chosen
//! because they produced accepted glyco PSMs, plus the 30 proteins those PSMs matched.
//! Small enough to commit (187 KB gz + 28 KB), real enough to exercise oxonium gating,
//! backbone generation, the fused selector and the collapse. A fixture with no
//! oxonium-positive scans pins only the header and would not have caught the above.
//!
//! To regenerate after an INTENTIONAL glyco change:
//!
//! ```text
//! cargo build --release -p andes
//! ./target/release/andes \
//!   --spectrum test-fixtures/glyco_fixture.mgf.gz \
//!   --database test-fixtures/glyco_fixture.fasta \
//!   --glyco --glyco-tol-ppm 20 --fragmentation HCD \
//!   --output-pin /tmp/g.pin
//! cp /tmp/g.glyco.pin test-fixtures/parity/goldens/glyco.pin
//! ```
//!
//! Regenerating is a deliberate act: diff the columns first and know which ones moved.

use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .to_path_buf()
}

fn read_sorted_rows(p: &std::path::Path) -> (String, Vec<String>) {
    let text = std::fs::read_to_string(p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
    let mut lines = text.lines();
    let header = lines.next().unwrap_or_default().to_string();
    let mut rows: Vec<String> = lines.filter(|l| !l.trim().is_empty()).map(str::to_string).collect();
    // Sort on the SpecId column, not the whole line: a last-digit float
    // difference in a later column would otherwise reorder rows and produce a
    // spurious mismatch against the golden.
    rows.sort_by(|x, y| spec_id(x).cmp(spec_id(y)));
    (header, rows)
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
fn glyco_pin_matches_golden_after_sort() {
    let root = workspace_root();
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let spectra = root.join("test-fixtures/glyco_fixture.mgf.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    let golden = root.join("test-fixtures/parity/goldens/glyco.pin");
    for f in [&spectra, &fasta, &golden] {
        assert!(f.exists(), "fixture missing: {}", f.display());
    }

    let outdir = tempfile::tempdir().expect("tempdir");
    let out_pin = outdir.path().join("actual.pin");
    let status = Command::new(&binary)
        .arg("--spectrum").arg(&spectra)
        .arg("--database").arg(&fasta)
        .arg("--glyco")
        .arg("--glyco-tol-ppm").arg("20")
        // NOTE ON COVERAGE. This resolves to `cid_lowres_tryp`, so the golden guards the
        // LOW-RES glyco path: column set, row count, oxonium gating, backbone generation,
        // the fused selector and the collapse. It does NOT guard the high-res
        // `tight_high_res` window that commit 539a3857 changed -- and an MGF fixture
        // CANNOT: reaching a high-res model needs --fragment-tol-ppm, which calls
        // `set_fragment_tol_override` and replaces `mme`, collapsing both branches of
        // that conditional to the same window. Guarding it requires an mzML fixture,
        // where the analyzer is auto-detected and the override is ignored. Verified by
        // experiment, not assumed.
        .arg("--fragmentation").arg("HCD")
        // Pin the taxon EXPLICITLY. The fixture FASTA carries an E. coli background (to
        // create candidate competition), and those headers carry OX= tags, so
        // `--glyco-taxon auto` resolves the FASTA to CmahCompetent and KEEPS NeuGc --
        // the opposite of the validated human config. Without this the golden would pin
        // whatever the padding happens to imply rather than the intended behaviour.
        .arg("--glyco-taxon").arg("human")
        .arg("--output-pin").arg(&out_pin)
        .status()
        .expect("run andes");
    assert!(status.success(), "andes exited {status}");

    // --glyco writes alongside --output-pin with a .glyco.pin suffix.
    let actual = outdir.path().join("actual.glyco.pin");
    assert!(actual.exists(), "glyco PIN not written: {}", actual.display());

    let (g_hdr, g_rows) = read_sorted_rows(&golden);
    let (a_hdr, a_rows) = read_sorted_rows(&actual);

    assert_eq!(g_hdr, a_hdr, "glyco PIN column set changed");
    assert_eq!(
        g_rows.len(),
        a_rows.len(),
        "glyco PIN row count changed: golden {} vs actual {}",
        g_rows.len(),
        a_rows.len()
    );
    for (i, (g, a)) in g_rows.iter().zip(a_rows.iter()).enumerate() {
        if let Some(d) = row_diff(g, a) {
            panic!("glyco PIN row mismatch at sorted index {i} ({d})\n golden: {g}\n actual: {a}");
        }
    }
}
