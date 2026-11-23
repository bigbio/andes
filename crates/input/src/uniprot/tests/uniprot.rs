use super::*;
use crate::annotation::{enzymatic_span, knockoff_positions};
use model::ChainEnds;
use model::{Enzyme, VariantRegion};

const XML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot">
<entry dataset="Swiss-Prot" version="3">
  <accession>P00001</accession>
  <accession>Q99999</accession>
  <name>TEST1_HUMAN</name>
  <protein>
    <recommendedName><fullName evidence="1">Test protein &amp; friend</fullName></recommendedName>
    <alternativeName><fullName>Other name</fullName></alternativeName>
    <component><recommendedName><fullName>Chain A</fullName></recommendedName></component>
  </protein>
  <gene><name type="synonym">SYN1</name><name type="primary">TST1</name></gene>
  <gene><name type="primary">TST2</name></gene>
  <organism>
    <name type="scientific">Homo sapiens</name>
    <name type="common">Human</name>
    <dbReference type="NCBI Taxonomy" id="9606"/>
  </organism>
  <organismHost><name type="scientific">Mus musculus</name></organismHost>
  <comment type="alternative products">
    <isoform><id>P00001-1</id><name>1</name><sequence type="displayed"/></isoform>
    <isoform><id>P00001-2</id><name>2</name><sequence type="described" ref="VSP_1"/></isoform>
  </comment>
  <comment type="function"><text>Some function.</text></comment>
  <dbReference type="PDB" id="1ABC"/>
  <proteinExistence type="evidence at protein level"/>
  <feature type="modified residue" description="N-acetylserine" evidence="1">
    <location><position position="2"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine; by PKA" evidence="1 2">
    <location><position position="6"/></location>
  </feature>
  <feature type="modified residue" description="Phosphothreonine" evidence="2">
    <location><position position="8"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1">
    <location sequence="P00001-2"><position position="14"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1">
    <location sequence="P00001-1"><position position="6"/></location>
  </feature>
  <feature type="modified residue" description="N6-lactoyllysine" evidence="1">
    <location><position position="9"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1">
    <location><position position="10"/></location>
  </feature>
  <feature type="cross-link" description="Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin); alternate" evidence="1">
    <location><position position="17"/></location>
  </feature>
  <feature type="cross-link" description="Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)" evidence="1">
    <location><position position="9"/></location>
  </feature>
  <feature type="glycosylation site" description="N-linked (GlcNAc...) asparagine" evidence="1">
    <location><position position="12"/></location>
  </feature>
  <feature type="modified residue" description="Phosphotyrosine" evidence="3">
    <location><position position="16"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000001" description="In dbSNP:rs1." evidence="1">
    <original>K</original><variation>E</variation>
    <location><position position="4"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000002" evidence="1">
    <original>S</original><variation>L</variation><variation>P</variation>
    <location><position position="14"/></location>
  </feature>
  <feature type="sequence variant" evidence="1">
    <original>AG</original><variation>VVV</variation>
    <location><begin position="10"/><end position="11"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000003" description="Missing." evidence="1">
    <location><position position="13"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000004" evidence="1">
    <original>Q</original><variation>H</variation>
    <location><position position="3"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000005" description="In dbSNP:rs5.">
    <original>E</original><variation>K</variation>
    <location><position position="19"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000006" evidence="1">
    <original>R</original><variation>W</variation>
    <location><begin status="unknown"/><end position="5"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000007" evidence="1">
    <location><position position="1"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000008" evidence="1">
    <original>S</original><variation>A</variation>
    <location sequence="P00001-2"><position position="2"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000009" evidence="1">
    <original>L</original><variation>I</variation>
    <location><position position="13"/></location>
  </feature>
  <feature type="mutagenesis site" description="Loss of function." evidence="1">
    <original>S</original><variation>A</variation>
    <location><position position="6"/></location>
  </feature>
  <evidence type="ECO:0000269" key="1"><source><dbReference type="PubMed" id="1"/></source></evidence>
  <evidence type="ECO:0000250" key="2"/>
  <evidence type="ECO:0000256" key="3"/>
  <sequence length="20" mass="2300" checksum="0" modified="2000-01-01" version="2">MSEKRSPTKAGNLSTYKDER</sequence>
</entry>
<entry dataset="TrEMBL" version="1">
  <accession>A0A000</accession>
  <name>A0A000_YEAST</name>
  <protein><submittedName><fullName>Uncharacterized protein</fullName></submittedName></protein>
  <organism>
    <name type="scientific">Saccharomyces cerevisiae</name>
    <dbReference type="NCBI Taxonomy" id="4932"/>
  </organism>
  <proteinExistence type="predicted"/>
  <sequence length="8" version="1">
    MKVLL
    AGR
  </sequence>
</entry>
</uniprot>
"#;

const CANONICAL: &[u8] = b"MSEKRSPTKAGNLSTYKDER";

fn load(opts: UniProtOptions) -> (ProteinDb, UniProtLoadStats) {
    load_uniprot_xml(XML.as_bytes(), &opts).unwrap()
}

fn sites(p: &Protein) -> Vec<(u32, f64)> {
    p.mod_sites()
        .iter()
        .map(|s| (s.position, s.mass_delta))
        .collect()
}

#[test]
fn entries_parse_with_uniprot_fasta_headers() {
    let entries: Vec<UniProtEntry> = UniProtXmlReader::new(XML.as_bytes())
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(entries.len(), 2);
    let e = &entries[0];
    assert_eq!(e.accession, "P00001");
    assert_eq!(e.name, "TEST1_HUMAN");
    assert_eq!(e.full_name, "Test protein & friend");
    assert_eq!(e.gene, "TST1", "primary name of the first gene only");
    assert_eq!(e.organism, "Homo sapiens", "not the organism host");
    assert_eq!(e.taxon_id, "9606");
    assert_eq!(e.sequence, CANONICAL);
    assert_eq!(e.displayed_isoforms().collect::<Vec<_>>(), ["P00001-1"]);
    assert_eq!(e.isoforms.len(), 2);
    assert_eq!(e.isoforms[1].name, "2");
    assert_eq!(e.isoforms[1].splice_refs, ["VSP_1"]);
    assert_eq!(e.features.len(), 22);
    assert_eq!(e.evidence.len(), 3);
    assert_eq!(e.fasta_accession(), "sp|P00001|TEST1_HUMAN");
    assert_eq!(
        e.fasta_description(),
        "Test protein & friend OS=Homo sapiens OX=9606 GN=TST1 PE=1 SV=2"
    );
    let range = &e.features[17];
    assert_eq!(range.id, "VAR_000006");
    assert_eq!(range.location, None, "an unknown begin leaves no location");
    let multi = &e.features[12];
    assert_eq!(multi.variations, vec!["L".to_string(), "P".to_string()]);

    let t = &entries[1];
    assert_eq!(t.fasta_accession(), "tr|A0A000|A0A000_YEAST");
    assert_eq!(
        t.fasta_description(),
        "Uncharacterized protein OS=Saccharomyces cerevisiae OX=4932 PE=4 SV=1"
    );
    assert_eq!(t.sequence, b"MKVLLAGR");
}

#[test]
fn sequences_only_by_default() {
    let (db, stats) = load(UniProtOptions::default());
    let accs: Vec<&str> = db.proteins.iter().map(|p| p.accession.as_str()).collect();
    assert_eq!(accs, ["sp|P00001|TEST1_HUMAN", "tr|A0A000|A0A000_YEAST"]);
    assert_eq!(db.proteins[0].sequence, CANONICAL);
    assert_eq!(stats.entries, 2);
    assert_eq!(stats.variant_proteins, 0);
    let ann = db.proteins[1].annotation.as_deref().unwrap();
    assert!(ann.mod_sites.is_empty() && ann.variant.is_empty() && ann.chain.is_none());
    assert_eq!(&*ann.origin.as_ref().unwrap().family, "A0A000");
}

#[test]
fn modification_sites_follow_the_evidence_level() {
    const ACETYL: f64 = 42.010565;
    const PHOSPHO: f64 = 79.966331;
    let (db, stats) = load(UniProtOptions::default());
    assert_eq!(
        sites(&db.proteins[0]),
        vec![(1, ACETYL), (5, PHOSPHO), (16, GLYGLY_MASS)],
        "experimental sites; the displayed-isoform duplicate of S6 is merged"
    );
    assert_eq!(stats.mod_sites, 3);
    assert_eq!(
        stats.below_evidence, 2,
        "T8 (by similarity), Y16 (automatic)"
    );
    assert_eq!(stats.unusable_features, 1, "phospho on A10");
    assert_eq!(stats.isoform_features, 1, "S14 of isoform 2");
    assert_eq!(
        stats.unknown_modifications.keys().collect::<Vec<_>>(),
        [
            "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)",
            "N-linked (GlcNAc...) asparagine",
            "N6-lactoyllysine"
        ]
    );
    let anchors: Vec<SiteAnchor> = db.proteins[0]
        .mod_sites()
        .iter()
        .map(|s| s.anchor)
        .collect();
    assert_eq!(
        anchors,
        [
            SiteAnchor::NTerminus,
            SiteAnchor::Residue,
            SiteAnchor::Residue
        ]
    );

    let curated = UniProtOptions {
        evidence: EvidenceLevel::Curated,
        ..UniProtOptions::default()
    };
    let (db, _) = load(curated);
    assert_eq!(
        sites(&db.proteins[0]),
        vec![(1, ACETYL), (5, PHOSPHO), (7, PHOSPHO), (16, GLYGLY_MASS)]
    );
    let any = UniProtOptions {
        evidence: EvidenceLevel::Any,
        ..UniProtOptions::default()
    };
    let (db, _) = load(any);
    assert_eq!(
        sites(&db.proteins[0]),
        vec![
            (1, ACETYL),
            (5, PHOSPHO),
            (7, PHOSPHO),
            (15, PHOSPHO),
            (16, GLYGLY_MASS)
        ]
    );
}

#[test]
fn variants_become_proteins_carrying_only_the_change() {
    let opts = UniProtOptions {
        variants: true,
        ..UniProtOptions::default()
    };
    let (db, stats) = load(opts);
    let accs: Vec<&str> = db.proteins.iter().map(|p| p.accession.as_str()).collect();
    assert_eq!(
        accs,
        [
            "sp|P00001|TEST1_HUMAN",
            "sp|P00001_VAR_000001_K4E|TEST1_HUMAN",
            "sp|P00001_VAR_000002_S14L|TEST1_HUMAN",
            "sp|P00001_VAR_000002_S14P|TEST1_HUMAN",
            "sp|P00001_VAR_AG10-11VVV|TEST1_HUMAN",
            "sp|P00001_VAR_000003_L13del|TEST1_HUMAN",
            "tr|A0A000|A0A000_YEAST",
        ]
    );
    assert_eq!(stats.variant_proteins, 5);
    assert_eq!(stats.below_evidence, 3);
    assert_eq!(stats.unusable_features, 5);
    assert_eq!(stats.isoform_features, 2);

    let k4e = &db.proteins[1];
    assert_eq!(k4e.sequence, b"MSEERSPTKAGNLSTYKDER");
    assert_eq!(
        k4e.variant_region(),
        Some(VariantRegion { start: 3, end: 4 })
    );
    assert!(k4e.description.ends_with("SV=2 VARIANT=VAR_000001_K4E"));
    assert_eq!(sites(k4e), sites(&db.proteins[0]), "same-length change");
    let o = k4e.origin().unwrap();
    assert_eq!(o.kind, model::EntryKind::Variant);
    assert_eq!(&*o.family, "P00001");
    assert_eq!(&*o.feature, "VAR_000001_K4E");
    assert!(
        (o.variant_delta_da - 0.947630).abs() < 1e-5,
        "{}",
        o.variant_delta_da
    );
    let canonical = db.proteins[0].origin().unwrap();
    assert_eq!(canonical.kind, model::EntryKind::Reference);
    assert_eq!((&*canonical.family, &*canonical.feature), ("P00001", ""));

    let ins = &db.proteins[4];
    assert_eq!(ins.sequence, b"MSEKRSPTKVVVNLSTYKDER");
    assert_eq!(
        ins.variant_region(),
        Some(VariantRegion { start: 9, end: 12 })
    );
    assert_eq!(
        sites(ins).iter().map(|s| s.0).collect::<Vec<_>>(),
        [1, 5, 17]
    );
    assert_eq!(ins.sequence[17], b'K');

    let del = &db.proteins[5];
    assert_eq!(del.sequence, b"MSEKRSPTKAGNSTYKDER");
    assert_eq!(
        del.variant_region(),
        Some(VariantRegion { start: 12, end: 12 })
    );
    assert_eq!(
        sites(del).iter().map(|s| s.0).collect::<Vec<_>>(),
        [1, 5, 15]
    );
    assert_eq!(db.proteins[0].variant_region(), None);
}

#[test]
fn variants_in_repeats_get_the_region_that_differs() {
    let xml = r#"<uniprot><entry dataset="Swiss-Prot">
  <accession>P00009</accession><name>REP_HUMAN</name>
  <feature type="sequence variant" evidence="1">
    <original>Q</original><variation>QQ</variation>
    <location><position position="4"/></location>
  </feature>
  <feature type="sequence variant" evidence="1">
    <location><position position="11"/></location>
  </feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="15">MWAQRSTYKQQQRLE</sequence>
</entry></uniprot>"#;
    let opts = UniProtOptions {
        variants: true,
        ..UniProtOptions::default()
    };
    let (db, _) = load_uniprot_xml(xml.as_bytes(), &opts).unwrap();
    let dup = &db.proteins[1];
    assert_eq!(dup.sequence, b"MWAQQRSTYKQQQRLE");
    assert_eq!(
        dup.variant_region(),
        Some(VariantRegion { start: 4, end: 4 })
    );
    let del = &db.proteins[2];
    assert_eq!(del.sequence, b"MWAQRSTYKQQRLE");
    assert_eq!(
        del.variant_region(),
        Some(VariantRegion { start: 11, end: 9 })
    );
}

const PROCESSED_XML: &str = r#"<uniprot>
<entry dataset="Swiss-Prot">
  <accession>P00010</accession><name>TEN_HUMAN</name>
  <feature type="signal peptide" evidence="1"><location><begin position="1"/><end position="5"/></location></feature>
  <feature type="propeptide" id="PRO_0" evidence="1"><location><begin position="6"/><end position="8"/></location></feature>
  <feature type="chain" id="PRO_1" description="Ten"><location><begin position="9"/><end position="30"/></location></feature>
  <feature type="peptide" id="PRO_2" description="Ten peptide"><location><begin position="13"/><end position="18"/></location></feature>
  <feature type="chain" id="PRO_3"><location><begin position="2"/><end position="30"/></location></feature>
  <feature type="chain" id="PRO_4"><location><begin position="1"/><end position="26"/></location></feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1"><location><position position="13"/></location></feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="30">MALWARGRDTHKSEIAHRFKDLGEEHFKAR</sequence>
</entry>
<entry dataset="Swiss-Prot">
  <accession>P00011</accession><name>ELEVEN_HUMAN</name>
  <protein><recommendedName><fullName>Eleven protein</fullName></recommendedName></protein>
  <gene><name type="primary">ELV</name></gene>
  <organism><name type="scientific">Homo sapiens</name><dbReference type="NCBI Taxonomy" id="9606"/></organism>
  <comment type="alternative products">
    <isoform><id>P00011-1</id><name>1</name><sequence type="displayed"/></isoform>
    <isoform><id>P00011-2</id><name>2</name><sequence type="described" ref="VSP_A VSP_B"/></isoform>
    <isoform><id>P00011-3</id><name>3</name><sequence type="described" ref="VSP_A VSP_C"/></isoform>
    <isoform><id>P00011-4</id><name>4</name><sequence type="external"/></isoform>
  </comment>
  <feature type="modified residue" description="Phosphoserine" evidence="1"><location><position position="2"/></location></feature>
  <feature type="modified residue" description="Phosphothreonine" evidence="1"><location><position position="11"/></location></feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1"><location sequence="P00011-2"><position position="14"/></location></feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1"><location sequence="P00011-3"><position position="2"/></location></feature>
  <feature type="sequence variant" evidence="1"><original>T</original><variation>A</variation><location sequence="P00011-2"><position position="9"/></location></feature>
  <feature type="splice variant" id="VSP_A"><original>STN</original><variation>Q</variation><location><begin position="2"/><end position="4"/></location></feature>
  <feature type="splice variant" id="VSP_B"><location><begin position="20"/><end position="27"/></location></feature>
  <feature type="splice variant" id="VSP_C"><original>XXXX</original><variation>A</variation><location><begin position="10"/><end position="13"/></location></feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="27">MSTNPKPQRKTKRNTSRRPQDVKFPGG</sequence>
</entry>
</uniprot>"#;

#[test]
fn isoforms_are_rebuilt_with_their_sites_and_variants() {
    let opts = UniProtOptions {
        variants: true,
        isoforms: true,
        ..UniProtOptions::default()
    };
    let (db, stats) = load_uniprot_xml(PROCESSED_XML.as_bytes(), &opts).unwrap();
    let accs: Vec<&str> = db.proteins.iter().map(|p| p.accession.as_str()).collect();
    assert_eq!(
        accs,
        [
            "sp|P00010|TEN_HUMAN",
            "sp|P00011|ELEVEN_HUMAN",
            "sp|P00011-2|ELEVEN_HUMAN",
            "sp|P00011-2_VAR_T9A|ELEVEN_HUMAN",
        ]
    );
    let iso = &db.proteins[2];
    assert_eq!(iso.sequence, b"MQPKPQRKTKRNTSRRP");
    assert_eq!(
        iso.description,
        "Isoform 2 of Eleven protein OS=Homo sapiens OX=9606 GN=ELV"
    );
    assert_eq!(sites(iso), vec![(8, 79.966331), (13, 79.966331)]);
    assert_eq!(iso.variant_region(), None);
    let var = &db.proteins[3];
    assert_eq!(var.sequence, b"MQPKPQRKAKRNTSRRP");
    assert_eq!(
        var.variant_region(),
        Some(VariantRegion { start: 8, end: 9 })
    );
    assert!(var.description.starts_with("Isoform 2 of Eleven protein"));
    assert_eq!(stats.isoform_proteins, 1);
    assert_eq!(
        stats.unbuilt_isoforms, 2,
        "isoform 3 does not fit, 4 is external"
    );
    assert_eq!(stats.isoform_features, 1, "the site on unbuilt isoform 3");
    assert_eq!(stats.variant_proteins, 1);
    assert_eq!(stats.mod_sites, 4, "S13; S2 and T11; the isoform's own S14");

    let (db, stats) =
        load_uniprot_xml(PROCESSED_XML.as_bytes(), &UniProtOptions::default()).unwrap();
    assert_eq!(db.len(), 2);
    assert_eq!(stats.isoform_features, 2);
}

#[test]
fn mature_chains_keep_only_new_termini() {
    let load = |evidence| {
        let opts = UniProtOptions {
            chains: true,
            evidence,
            ..UniProtOptions::default()
        };
        load_uniprot_xml(PROCESSED_XML.as_bytes(), &opts).unwrap()
    };
    let (db, stats) = load(EvidenceLevel::Experimental);
    let chain = &db.proteins[1];
    assert_eq!(chain.accession, "sp|P00010_PRO_1|TEN_HUMAN");
    assert!(chain.description.ends_with(" CHAIN=PRO_1:9-30"));
    assert_eq!(chain.sequence, b"DTHKSEIAHRFKDLGEEHFKAR");
    assert_eq!(
        chain.chain_ends(),
        Some(ChainEnds {
            n_term: true,
            c_term: false
        })
    );
    assert_eq!(sites(chain), vec![(4, 79.966331)], "S13 on the chain's S5");
    assert_eq!(stats.chain_proteins, 1);
    assert_eq!(
        stats.below_evidence, 2,
        "the released peptide and PRO_4, which has no evidence"
    );
    let (db, stats) = load(EvidenceLevel::Curated);
    assert_eq!(stats.chain_proteins, 3);
    let peptide = &db.proteins[2];
    assert_eq!(peptide.accession, "sp|P00010_PRO_2|TEN_HUMAN");
    assert_eq!(peptide.sequence, b"SEIAHR");
    assert_eq!(
        peptide.chain_ends(),
        Some(ChainEnds {
            n_term: true,
            c_term: true
        })
    );
    let c_chain = &db.proteins[3];
    assert_eq!(c_chain.accession, "sp|P00010_PRO_4|TEN_HUMAN");
    assert_eq!(c_chain.sequence, b"MALWARGRDTHKSEIAHRFKDLGEEH");
    assert_eq!(
        c_chain.chain_ends(),
        Some(ChainEnds {
            n_term: false,
            c_term: true
        })
    );
}

#[test]
fn evidence_levels_nest() {
    let codes = |c: &'static [&'static str]| c.iter().copied();
    const EXP: &[&str] = &["ECO:0000269"];
    const SIM: &[&str] = &["ECO:0000250"];
    const AUTO: &[&str] = &["ECO:0000256"];
    for (level, admitted) in [
        (EvidenceLevel::Experimental, [true, false, false, false]),
        (EvidenceLevel::Curated, [true, true, false, false]),
        (EvidenceLevel::Any, [true, true, true, true]),
    ] {
        assert_eq!(level.admits(codes(EXP)), admitted[0], "{level:?}");
        assert_eq!(level.admits(codes(SIM)), admitted[1], "{level:?}");
        assert_eq!(level.admits(codes(AUTO)), admitted[2], "{level:?}");
        assert_eq!(level.admits(codes(&[])), admitted[3], "{level:?}");
    }
    assert_eq!(
        EvidenceLevel::from_name(" Curated "),
        Some(EvidenceLevel::Curated)
    );
    assert_eq!(EvidenceLevel::from_name("strong"), None);
    for l in [
        EvidenceLevel::Experimental,
        EvidenceLevel::Curated,
        EvidenceLevel::Any,
    ] {
        assert_eq!(EvidenceLevel::from_name(l.as_str()), Some(l));
    }
}

#[test]
fn malformed_documents_are_errors() {
    let cut = XML.rfind("<sequence length=\"8\"").unwrap();
    let r: Vec<_> = UniProtXmlReader::new(&XML.as_bytes()[..cut]).collect();
    assert_eq!(r.len(), 2);
    assert!(r[0].is_ok());
    assert!(matches!(
        &r[1],
        Err(UniProtParseError::Truncated { element }) if element == "entry"
    ));
    let opts = UniProtOptions::default();
    assert!(load_uniprot_xml(&XML.as_bytes()[..cut], &opts).is_err());
    let html = "<html><body><p>not UniProt</p></body></html>";
    assert!(matches!(
        load_uniprot_xml(html.as_bytes(), &opts),
        Err(UniProtParseError::NoEntries)
    ));
    let mismatched = "<uniprot><entry><accession>X</name></entry></uniprot>";
    let r: Result<Vec<_>, _> = UniProtXmlReader::new(mismatched.as_bytes()).collect();
    assert!(matches!(r, Err(UniProtParseError::Xml { .. })));
    let no_seq =
        r#"<uniprot><entry dataset="Swiss-Prot"><accession>X1</accession></entry></uniprot>"#;
    let r: Result<Vec<_>, _> = UniProtXmlReader::new(no_seq.as_bytes()).collect();
    assert!(matches!(
        r,
        Err(UniProtParseError::MissingSequence { accession }) if accession == "X1"
    ));
}

#[test]
fn change_labels() {
    assert_eq!(change_label(b"D", 7, 7, b"H"), "D7H");
    assert_eq!(change_label(b"KR", 120, 121, b"EE"), "KR120-121EE");
    assert_eq!(change_label(b"F", 508, 508, b""), "F508del");
}

const TRYPSIN_KO: KnockoffOptions = KnockoffOptions {
    per_variant: 4,
    enzyme: Enzyme::Trypsin,
};

fn repeats() -> Vec<u8> {
    b"GGGGDLGGGK".repeat(8)
}

#[test]
fn knockoffs_take_the_same_residue_in_another_matched_span() {
    let seq = repeats();
    let qs = knockoff_positions(&seq, 4, &[4], TRYPSIN_KO, 7);
    assert_eq!(qs.len(), 4);
    for &q in &qs {
        assert_eq!(seq[q], b'D');
        assert_eq!(seq[q + 1], b'L');
        assert!(q.abs_diff(4) >= 10);
        assert_ne!(enzymatic_span(&seq, q, Enzyme::Trypsin).0, 0);
    }
    let mut sorted = qs.clone();
    sorted.dedup();
    assert_eq!(sorted.len(), 4, "drawn without replacement");
    assert_eq!(
        qs,
        knockoff_positions(&seq, 4, &[4], TRYPSIN_KO, 7),
        "deterministic"
    );
}

#[test]
fn knockoffs_keep_clear_of_listed_variants_and_mismatched_contexts() {
    let seq = repeats();
    let listed: Vec<usize> = (0..8).map(|i| 4 + 10 * i).collect();
    assert_eq!(
        knockoff_positions(&seq, 4, &listed, TRYPSIN_KO, 1),
        Vec::<usize>::new()
    );
    let listed: Vec<usize> = (0..7).map(|i| 4 + 10 * i).collect();
    assert_eq!(
        knockoff_positions(&seq, 4, &listed, TRYPSIN_KO, 1),
        vec![74]
    );
    let mut other = b"GGGGDLGGGK".to_vec();
    other.extend(b"GGGGDEGGGK".repeat(3));
    assert_eq!(
        knockoff_positions(&other, 4, &[4], TRYPSIN_KO, 1),
        Vec::<usize>::new()
    );
    let mut short = b"GGGGDLGGGK".to_vec();
    short.extend(b"GDLGK".repeat(3));
    assert_eq!(
        knockoff_positions(&short, 4, &[4], TRYPSIN_KO, 1),
        Vec::<usize>::new()
    );
}

#[test]
fn a_sequon_is_matched_by_a_sequon() {
    let mut seq = b"GGGGNLSGGK".to_vec();
    seq.extend(b"GGGGNLGGGK".repeat(3));
    assert_eq!(
        knockoff_positions(&seq, 4, &[4], TRYPSIN_KO, 1),
        Vec::<usize>::new()
    );
}

#[test]
fn knockoff_entries_copy_the_substitution_and_name_the_variant() {
    let seq = String::from_utf8(repeats()).unwrap();
    let xml = format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot">
<entry dataset="Swiss-Prot" version="1">
  <accession>P00009</accession>
  <name>KO_HUMAN</name>
  <protein><recommendedName><fullName>Knockoff test</fullName></recommendedName></protein>
  <organism><name type="scientific">Homo sapiens</name></organism>
  <proteinExistence type="evidence at protein level"/>
  <feature type="sequence variant" id="VAR_000009" evidence="1">
    <original>D</original><variation>N</variation>
    <location><position position="5"/></location>
  </feature>
  <evidence key="1" type="ECO:0000269"/>
  <sequence length="80" mass="1" checksum="0" modified="2000-01-01" version="1">{seq}</sequence>
</entry>
</uniprot>"#
    );
    let opts = UniProtOptions {
        variants: true,
        knockoffs: Some(TRYPSIN_KO),
        ..UniProtOptions::default()
    };
    let (db, stats) = load_uniprot_xml(xml.as_bytes(), &opts).unwrap();
    assert_eq!((stats.variant_proteins, stats.knockoff_proteins), (1, 4));
    let real = db.proteins[1].origin().unwrap().clone();
    assert!(real.knockoff_of.is_none() && !db.proteins[1].is_knockoff());
    for (i, p) in db.proteins[2..].iter().enumerate() {
        let o = p.origin().unwrap();
        assert!(p.is_knockoff());
        assert_eq!(o.kind, model::EntryKind::Variant);
        assert_eq!(o.knockoff_of.as_deref(), Some("VAR_000009_D5N"));
        assert_eq!(o.variant_delta_da, real.variant_delta_da);
        assert!(
            o.feature.starts_with(&format!("KO{}_D", i + 1)),
            "{}",
            o.feature
        );
        assert!(o.feature.ends_with("N_VAR_000009_D5N"), "{}", o.feature);
        let changed: Vec<usize> = (0..80).filter(|&j| p.sequence[j] != repeats()[j]).collect();
        assert_eq!(changed.len(), 1);
        assert_eq!(p.sequence[changed[0]], b'N');
        assert_eq!(p.substitutions().len(), 1);
    }
    let off = UniProtOptions {
        knockoffs: None,
        ..opts
    };
    let (_, stats) = load_uniprot_xml(xml.as_bytes(), &off).unwrap();
    assert_eq!(stats.knockoff_proteins, 0);
}
