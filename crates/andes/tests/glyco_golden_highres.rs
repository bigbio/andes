//! HIGH-RES glyco serve-path golden.
//!
//! Companion to `glyco_golden.rs`, which guards the LOW-RES path. This one exists
//! because an MGF fixture *cannot* reach the high-res branch: selecting a high-res
//! model for MGF requires `--fragment-tol-ppm`, which calls
//! `set_fragment_tol_override` and replaces `mme` -- collapsing both arms of the
//! `tight_high_res` conditional in `scored_spectrum.rs` to the same window. So the
//! very flag needed to get there is the flag that neutralises what we want to guard.
//!
//! An mzML fixture carries analyzer metadata, so the instrument is auto-detected and
//! the override is never involved. That makes this the test that actually pins the
//! window commit 539a3857 ("fix(train)") changed while asserting "Serving is untouched".
//!
//! It guards two things at once:
//!   1. the high-res glyco scoring path (the PIN comparison), and
//!   2. that instrument detection survives a `.gz` input (the model assertion) --
//!      a gzipped mzML used to skip detection entirely and silently fall back to
//!      `cid_lowres_tryp`, i.e. the low-res model on high-res Orbitrap data.
//!
//! Regenerating is a deliberate act: diff the columns first and know which ones moved.
//!
//!   ./target/release/andes \
//!     --spectrum test-fixtures/glyco_fixture.mzML.gz \
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

#[test]
fn glyco_highres_pin_matches_golden() {
    let root = repo_root();
    let spectra = root.join("test-fixtures/glyco_fixture.mzML.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    let golden = root.join("test-fixtures/parity/goldens/glyco_highres.pin");
    for f in [&spectra, &fasta, &golden] {
        assert!(f.exists(), "fixture missing: {}", f.display());
    }

    let tmp = std::env::temp_dir().join("andes_glyco_highres_golden");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).expect("tmp");
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
    assert!(actual.exists(), "glyco PIN not written: {}", actual.display());

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

    g_lines.sort_unstable();
    a_lines.sort_unstable();
    for (i, (g, a)) in g_lines.iter().zip(a_lines.iter()).enumerate() {
        assert_eq!(g, a, "glyco high-res PIN row mismatch at sorted index {i}");
    }
}
