//! A gzipped mzML must resolve to the SAME instrument and model as the plain file.
//!
//! `Path::extension()` returns `"gz"` for `run.mzML.gz`, so the `== "mzml"` guards on
//! the metadata-detection helpers used to bail out before opening the file. Instrument
//! detection returning `None` means the low-res default, so a gzipped high-res Orbitrap
//! run was silently searched with `cid_lowres_tryp` -- with nothing in the output saying
//! so. The guards sat one line above readers that are gz-capable (`open_buf_maybe_gz`),
//! so the reader was never the limitation.
//!
//! The fixture is 120 Orbitrap Fusion Lumos FTMS scans. Note the assertion is on the
//! MODEL rather than on identifications: this guards the detection path, and a wrong
//! model is a silent correctness bug long before it is an ID count.

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
fn gzipped_mzml_detects_the_same_high_res_model_as_plain() {
    let root = repo_root();
    let gz = root.join("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let fasta = root.join("test-fixtures/BSA.fasta");
    assert!(gz.exists(), "fixture missing: {}", gz.display());
    assert!(fasta.exists(), "fixture missing: {}", fasta.display());

    let tmp = std::env::temp_dir().join("andes_gz_detect");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).expect("tmp");

    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let run = |spectra: &PathBuf, tag: &str| -> String {
        let out = Command::new(&binary)
            .arg("--spectrum")
            .arg(spectra)
            .arg("--database")
            .arg(&fasta)
            .arg("--output-pin")
            .arg(tmp.join(format!("{tag}.pin")))
            .output()
            .expect("run andes");
        assert!(
            out.status.success(),
            "andes exited {} on {tag}\n{}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        );
        format!(
            "{}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        )
    };

    // Decompress to a sibling path so the ONLY difference between the two runs is the
    // `.gz` extension -- same bytes, same spectra, same everything else.
    let plain = tmp.join("orbitrap_lumos_120.mzML");
    {
        use std::io::{Read, Write};
        let f = std::fs::File::open(&gz).expect("open gz");
        let mut d = flate2::read::GzDecoder::new(f);
        let mut buf = Vec::new();
        d.read_to_end(&mut buf).expect("inflate");
        std::fs::File::create(&plain)
            .expect("create plain")
            .write_all(&buf)
            .expect("write plain");
    }

    let log_plain = run(&plain, "plain");
    let log_gz = run(&gz, "gz");

    // The fixture is Orbitrap FTMS, so both must land on the high-res model.
    assert!(
        log_plain.contains("hcd_qexactive_tryp"),
        "plain mzML did not select the high-res model; fixture or detection changed:\n{log_plain}"
    );
    assert!(
        log_gz.contains("hcd_qexactive_tryp"),
        "GZIPPED mzML did not select the high-res model -- instrument detection was \
         skipped and the run silently fell back to the low-res default:\n{log_gz}"
    );
    assert!(
        !log_gz.contains("cid_lowres_tryp"),
        "gzipped mzML fell back to cid_lowres_tryp:\n{log_gz}"
    );
}
