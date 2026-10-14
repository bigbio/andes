//! Load `test-fixtures/BSA.xml.gz` — the UniProt XML entry of the protein in
//! `BSA.fasta` (P02769, gzipped as UniProt distributes it) — and check that a
//! UniProt XML database reads like the equivalent FASTA, plus its annotation.

use std::path::PathBuf;

use input::{load_protein_database, DatabaseFormat, EvidenceLevel, UniProtOptions};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("test-fixtures")
        .join(name)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {name}: {e}"))
}

#[test]
fn uniprot_xml_reads_like_the_uniprot_fasta() {
    let xml = load_protein_database(&fixture("BSA.xml.gz"), &UniProtOptions::default()).unwrap();
    let fasta = load_protein_database(&fixture("BSA.fasta"), &UniProtOptions::default()).unwrap();
    assert_eq!(xml.format, DatabaseFormat::UniProtXml);
    assert_eq!(fasta.format, DatabaseFormat::Fasta);
    assert!(fasta.uniprot.is_none());
    assert_eq!(xml.db.len(), 1);
    let (x, f) = (&xml.db.proteins[0], &fasta.db.proteins[0]);
    assert_eq!(x.accession, f.accession);
    assert_eq!(x.sequence, f.sequence);
    // Same header layout as UniProt's FASTA (the recommended name has since
    // changed from "Serum albumin" to "Albumin").
    assert_eq!(
        x.description,
        "Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4"
    );
    assert_eq!(xml.uniprot.as_ref().unwrap().entries, 1);
}

#[test]
fn bsa_sites_and_variant_follow_the_evidence_level() {
    let load = |evidence, variants| {
        load_protein_database(
            &fixture("BSA.xml.gz"),
            &UniProtOptions { variants, evidence },
        )
        .unwrap()
    };
    // Experimental evidence: the in-vitro glycation of K151 only; the
    // phosphosites are all "by similarity".
    let exp = load(EvidenceLevel::Experimental, false);
    let sites = exp.db.proteins[0].mod_sites();
    assert_eq!(sites.len(), 1);
    assert_eq!(sites[0].position, 150);
    assert_eq!(exp.db.proteins[0].sequence[150], b'K');
    assert!((sites[0].mass_delta - 162.052823).abs() < 1e-6);

    let curated = load(EvidenceLevel::Curated, true);
    let p = &curated.db.proteins[0];
    assert_eq!(p.mod_sites().len(), 14);
    let phospho: Vec<u32> = p
        .mod_sites()
        .iter()
        .filter(|s| (s.mass_delta - 79.966331).abs() < 1e-6)
        .map(|s| s.position + 1)
        .collect();
    assert_eq!(phospho, [29, 82, 89, 107, 296, 442, 443, 445, 512, 569]);
    for s in p.mod_sites() {
        assert!(
            matches!(p.sequence[s.position as usize], b'S' | b'T' | b'K' | b'N'),
            "site {} on {}",
            s.position,
            p.sequence[s.position as usize] as char
        );
    }
    // The A214T variant (experimental evidence, no feature id).
    assert_eq!(curated.db.len(), 2);
    let v = &curated.db.proteins[1];
    assert_eq!(v.accession, "sp|P02769_VAR_A214T|ALBU_BOVIN");
    assert_eq!(v.sequence[213], b'T');
    assert_eq!(p.sequence[213], b'A');
    assert_eq!(
        v.variant_region(),
        Some(model::VariantRegion {
            start: 213,
            end: 214
        })
    );
}
