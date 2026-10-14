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
            "--uniprot-chains",
            "--uniprot-isoforms",
        ],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("14 annotated modification sites, 1 variant entries"),
        "{stderr}"
    );
    // The mature chain after the propeptide (25-607) is a new N-terminus;
    // BSA has no isoforms.
    assert!(
        stderr.contains("0 isoform entries (0 isoforms not rebuildable), 1 mature-chain entries"),
        "{stderr}"
    );
    // Ten phosphosites on the canonical entry, inherited by the variant and
    // the chain entry.
    assert!(
        stderr.contains("Annotated-site mod Phospho (+79.9663 Da): 30 sites on 3 target proteins"),
        "{stderr}"
    );
    // Every PIN row names either the canonical entry or the variant entry.
    let pin_text = std::fs::read_to_string(&pin).unwrap();
    assert!(pin_text.lines().skip(1).all(|l| l.contains("P02769")));
    // Provenance columns: a row is non-reference only when every protein it
    // names is a variant or chain entry; only those carry a DeltaRefScore. The
    // canonical BSA entry competes in most spectra.
    let header: Vec<&str> = pin_text.lines().next().unwrap().split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let (var, chain, iso) = (col("IsVariant"), col("IsChain"), col("IsIsoform"));
    let (comp, gap) = (col("RefCompetitor"), col("DeltaRefScore"));
    let (prot, rows): (usize, Vec<Vec<&str>>) = (
        col("Proteins"),
        pin_text
            .lines()
            .skip(1)
            .map(|l| l.split('\t').collect())
            .collect(),
    );
    assert!(
        rows.iter().any(|r| r[comp] == "1"),
        "no spectrum scored a canonical candidate"
    );
    for r in &rows {
        let named = &r[prot..];
        let alternative = named
            .iter()
            .all(|a| a.contains("_VAR_") || a.contains("_PRO_"));
        let flagged = r[var] == "1" || r[chain] == "1";
        assert_eq!(flagged, alternative, "{named:?}");
        assert_eq!(r[iso], "0");
        if !flagged {
            assert_eq!(r[gap].parse::<f64>().unwrap(), 0.0, "{named:?}");
        }
    }

    // The same mods over FASTA: no annotation, so the restricted mods never apply.
    let out = search(
        "test-fixtures/BSA.fasta",
        &dir.path().join("fasta.pin"),
        &["--mods", mods, "--uniprot-variants", "--uniprot-chains"],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Phospho (+79.9663 Da) is restricted to annotated sites"),
        "{stderr}"
    );
    assert!(stderr.contains("never placed"), "{stderr}");
    for flag in ["--uniprot-variants", "--uniprot-chains"] {
        assert!(
            stderr.contains(&format!("{flag} needs a UniProt XML database")),
            "{stderr}"
        );
    }

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

#[test]
fn restricted_artefact_mods_warn_only_with_variants() {
    let dir = tempfile::tempdir().unwrap();
    let mods = dir.path().join("mods.txt");
    std::fs::write(
        &mods,
        "57.021464,C,fix,any,Carbamidomethyl\n\
         0.984016,N,opt,any,Deamidated,sites=annotated\n",
    )
    .unwrap();
    let mods = mods.to_str().unwrap();
    let warning = "variants of the same mass lose the competing explanation";
    let out = search(
        "test-fixtures/BSA.xml.gz",
        &dir.path().join("on.pin"),
        &["--mods", mods, "--uniprot-variants"],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("WARN: Deamidated"), "{stderr}");
    assert!(stderr.contains(warning), "{stderr}");
    let out = search(
        "test-fixtures/BSA.xml.gz",
        &dir.path().join("off.pin"),
        &["--mods", mods],
    );
    assert!(!String::from_utf8_lossy(&out.stderr).contains(warning));
}
