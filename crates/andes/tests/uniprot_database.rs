//! End-to-end: a UniProt XML database (`--database *.xml.gz`) searches exactly
//! like the equivalent UniProt FASTA, and its annotation options run through.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn search(database: &str, pin: &Path, extra: &[&str]) -> Output {
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(fixture("test-fixtures/test.mgf.gz"))
        .arg("--database")
        .arg(fixture(database))
        .arg("--output-pin")
        .arg(pin)
        .args(["--decoy-prefix", "XXX_"])
        .args(extra)
        .output()
        .expect("run andes");
    assert!(
        out.status.success(),
        "andes failed on {database}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

#[test]
fn uniprot_xml_database_gives_the_fasta_pin() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta_pin, xml_pin) = (dir.path().join("fasta.pin"), dir.path().join("xml.pin"));
    search("test-fixtures/BSA.fasta", &fasta_pin, &[]);
    let out = search("test-fixtures/BSA.xml.gz", &xml_pin, &[]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("UniProt XML: 1 entries"), "{stderr}");
    let (fasta, xml) = (
        std::fs::read_to_string(&fasta_pin).unwrap(),
        std::fs::read_to_string(&xml_pin).unwrap(),
    );
    assert!(fasta.lines().count() > 1, "the FASTA search found PSMs");
    assert_eq!(fasta, xml, "same sequences and accessions, same PIN");
}

#[test]
fn annotated_site_mods_and_variants_run_end_to_end() {
    let dir = tempfile::tempdir().unwrap();
    let mods = dir.path().join("mods.txt");
    std::fs::write(
        &mods,
        "57.021464,C,fix,any,Carbamidomethyl\n\
         15.994915,M,opt,any,Oxidation\n\
         79.966331,S,opt,any,Phospho,sites=annotated\n\
         79.966331,T,opt,any,Phospho,sites=annotated\n",
    )
    .unwrap();
    let mods = mods.to_str().unwrap();
    let pin = dir.path().join("xml.pin");
    let out = search(
        "test-fixtures/BSA.xml.gz",
        &pin,
        &[
            "--mods",
            mods,
            "--uniprot-evidence",
            "curated",
            "--uniprot-variants",
        ],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("14 annotated modification sites, 1 variant entries"),
        "{stderr}"
    );
    // Ten phosphosites on the canonical entry, inherited by the variant entry.
    assert!(
        stderr.contains("Annotated-site mod Phospho (+79.9663 Da): 20 sites on 2 target proteins"),
        "{stderr}"
    );
    // Every PIN row names either the canonical entry or the variant entry.
    let pin_text = std::fs::read_to_string(&pin).unwrap();
    assert!(pin_text.lines().skip(1).all(|l| l.contains("P02769")));

    // The same mods over FASTA: no annotation, so the restricted mods never apply.
    let out = search(
        "test-fixtures/BSA.fasta",
        &dir.path().join("fasta.pin"),
        &["--mods", mods, "--uniprot-variants"],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Phospho (+79.9663 Da) is restricted to annotated sites"),
        "{stderr}"
    );
    assert!(stderr.contains("never placed"), "{stderr}");
    assert!(
        stderr.contains("--uniprot-variants needs a UniProt XML database"),
        "{stderr}"
    );

    // XML at the default (experimental) evidence keeps only BSA's glycation
    // site: there is no phosphosite, which is reported, not counted as one.
    let out = search(
        "test-fixtures/BSA.xml.gz",
        &dir.path().join("exp.pin"),
        &["--mods", mods],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("1 annotated modification sites"),
        "{stderr}"
    );
    assert!(stderr.contains("never placed"), "{stderr}");
}
