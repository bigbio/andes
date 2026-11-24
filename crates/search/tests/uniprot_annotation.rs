use std::collections::{BTreeSet, HashSet};

use input::{load_uniprot_xml, EvidenceLevel, UniProtOptions};
use model::enzyme::Enzyme;
use model::{
    AminoAcidSetBuilder, Modification, Protein, ProteinAnnotation, ProteinDb, VariantRegion,
};
use search::candidates::{
    expand_base_record, expand_base_record_bounded, for_each_record_form_masses_bounded,
    lazy_candidates_for_precursor, Candidate, FormWalkTables,
};
use search::candidates::records::MmapCandidateIndex;
use search::decoy::{decoy_db, DecoyStrategy, Proteases};
use search::{enumerate_candidates, SearchIndex, SearchParams};

fn spec(list: &str, fixed: bool) -> Vec<Modification> {
    model::modification::spec::parse_mod_list(list, fixed).unwrap()
}

fn annotated(item: &str) -> Modification {
    let mut m = spec(item, false).remove(0);
    m.annotated_site_mass = Some(m.mass_delta);
    m
}

const XML: &str = r#"<uniprot xmlns="http://uniprot.org/uniprot">
<entry dataset="Swiss-Prot">
  <accession>P00001</accession>
  <name>ONE_HUMAN</name>
  <feature type="modified residue" description="Phosphoserine" evidence="1">
    <location><position position="3"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine; by CK2" evidence="1">
    <location><position position="14"/></location>
  </feature>
  <feature type="modified residue" description="Phosphothreonine" evidence="2">
    <location><position position="21"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000001" evidence="1">
    <original>K</original><variation>N</variation>
    <location><position position="17"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000002" evidence="1">
    <location><position position="10"/></location>
  </feature>
  <evidence type="ECO:0000269" key="1"/>
  <evidence type="ECO:0000250" key="2"/>
  <sequence length="23" version="1">AGSPEKTLSYRSDSPAKVNSTQR</sequence>
</entry>
<entry dataset="Swiss-Prot">
  <accession>P00002</accession>
  <name>TWO_HUMAN</name>
  <sequence length="23" version="1">AGSPEKTLSYRSDSPAKVNSTQR</sequence>
</entry>
</uniprot>"#;

const PHOSPHO: f64 = 79.966331;

fn target_db(evidence: EvidenceLevel, variants: bool) -> ProteinDb {
    load_uniprot_xml(
        XML.as_bytes(),
        &UniProtOptions {
            variants,
            evidence,
            ..UniProtOptions::default()
        },
    )
    .unwrap()
    .0
}

fn params(restricted: bool) -> SearchParams {
    let mut builder = AminoAcidSetBuilder::new_standard();
    for r in ["S", "T"] {
        let item = format!("Phospho ({r})");
        builder = builder.add_variable_mod(if restricted {
            annotated(&item)
        } else {
            spec(&item, false).remove(0)
        });
    }
    let mut p = SearchParams::default_tryptic(builder.build().unwrap());
    p.min_length = 4;
    p.max_length = 30;
    p.max_missed_cleavages = 1;
    p.max_variable_mods_per_peptide = 2;
    p
}

fn phospho_positions(c: &Candidate) -> Vec<usize> {
    c.peptide
        .residues
        .iter()
        .enumerate()
        .filter(|(_, aa)| {
            aa.mod_
                .is_some_and(|m| (m.mass_delta - PHOSPHO).abs() < 1e-6)
        })
        .map(|(i, _)| c.start_offset_in_protein + i)
        .collect()
}

fn key(c: &Candidate) -> (bool, usize, usize, String, u64) {
    (
        c.is_decoy,
        c.protein_index,
        c.start_offset_in_protein,
        c.peptide.to_string(),
        (c.peptide.mass() * 1000.0).round() as u64,
    )
}

fn accession<'a>(idx: &'a SearchIndex, c: &Candidate) -> &'a str {
    &idx.db.proteins[c.protein_index].accession
}

#[test]
fn restricted_phospho_lands_only_on_annotated_sites() {
    let idx = SearchIndex::from_target_db(&target_db(EvidenceLevel::Experimental, false), "XXX");
    let cands: Vec<Candidate> = enumerate_candidates(&idx, &params(true), "XXX").collect();
    let mut placed: BTreeSet<(String, usize)> = BTreeSet::new();
    for c in cands.iter().filter(|c| !c.is_decoy) {
        for pos in phospho_positions(c) {
            placed.insert((accession(&idx, c).to_string(), pos));
        }
    }
    let expect: BTreeSet<(String, usize)> = [2, 13]
        .into_iter()
        .map(|p| ("sp|P00001|ONE_HUMAN".to_string(), p))
        .collect();
    assert_eq!(placed, expect);

    let idx = SearchIndex::from_target_db(&target_db(EvidenceLevel::Curated, false), "XXX");
    let curated: BTreeSet<usize> = enumerate_candidates(&idx, &params(true), "XXX")
        .filter(|c| !c.is_decoy)
        .flat_map(|c| phospho_positions(&c))
        .collect();
    assert_eq!(curated, BTreeSet::from([2, 13, 20]));

    let free: BTreeSet<usize> = enumerate_candidates(&idx, &params(false), "XXX")
        .filter(|c| !c.is_decoy && c.protein_index == 1)
        .flat_map(|c| phospho_positions(&c))
        .collect();
    assert_eq!(free, BTreeSet::from([2, 6, 8, 11, 13, 19, 20]));
}

#[test]
fn decoys_carry_the_sites_on_the_same_residues() {
    let target = target_db(EvidenceLevel::Curated, true);
    for strategy in [
        DecoyStrategy::Reverse,
        DecoyStrategy::Shuffle,
        DecoyStrategy::Sequon,
    ] {
        let idx = SearchIndex::from_target_db_with_strategy(&target, "XXX", None, strategy, 7);
        let n = target.len();
        for (t, d) in target.proteins.iter().zip(&idx.db.proteins[n..]) {
            assert_eq!(t.mod_sites().len(), d.mod_sites().len(), "{strategy:?}");
            let t_res: Vec<u8> = t
                .mod_sites()
                .iter()
                .map(|s| t.sequence[s.position as usize])
                .collect();
            let mut d_res: Vec<u8> = d
                .mod_sites()
                .iter()
                .map(|s| d.sequence[s.position as usize])
                .collect();
            let mut t_sorted = t_res.clone();
            t_sorted.sort_unstable();
            d_res.sort_unstable();
            assert_eq!(t_sorted, d_res, "{strategy:?}: sites follow their residue");
            assert_eq!(
                t.variant_region().is_some(),
                d.variant_region().is_some(),
                "{strategy:?}"
            );
            if let (Some(tv), Some(dv)) = (t.variant_region(), d.variant_region()) {
                assert_eq!(tv.end - tv.start, dv.end - dv.start, "{strategy:?}");
            }
        }
        for c in enumerate_candidates(&idx, &params(true), "XXX").filter(|c| c.is_decoy) {
            let p = &idx.db.proteins[c.protein_index];
            for pos in phospho_positions(&c) {
                assert!(
                    p.mod_sites().iter().any(|s| s.position as usize == pos),
                    "{strategy:?}: decoy phospho at unannotated {pos} of {}",
                    p.accession
                );
            }
        }
    }
    let idx =
        SearchIndex::from_target_db_with_strategy(&target, "XXX", None, DecoyStrategy::Reverse, 7);
    let decoy = &idx.db.proteins[target.len()];
    assert_eq!(decoy.accession, "XXX_sp|P00001|ONE_HUMAN");
    let pos: Vec<u32> = decoy.mod_sites().iter().map(|s| s.position).collect();
    assert_eq!(pos, [2, 13, 18]);
}

#[test]
fn variant_entries_yield_only_peptides_carrying_the_variant() {
    let target = target_db(EvidenceLevel::Experimental, true);
    let accs: Vec<&str> = target
        .proteins
        .iter()
        .map(|p| p.accession.as_str())
        .collect();
    assert_eq!(
        accs,
        [
            "sp|P00001|ONE_HUMAN",
            "sp|P00001_VAR_000001_K17N|ONE_HUMAN",
            "sp|P00001_VAR_000002_Y10del|ONE_HUMAN",
            "sp|P00002|TWO_HUMAN",
        ]
    );
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let p = params(false);
    let cands: Vec<Candidate> = enumerate_candidates(&idx, &p, "XXX").collect();
    let bare = |c: &Candidate| -> String {
        c.peptide
            .residues
            .iter()
            .map(|a| a.residue as char)
            .collect()
    };
    let canonical: HashSet<String> = cands
        .iter()
        .filter(|c| c.protein_index == 0)
        .map(bare)
        .collect();
    for (pi, region) in [(1usize, (16u32, 17u32)), (2, (9, 9))] {
        let from_variant: Vec<&Candidate> =
            cands.iter().filter(|c| c.protein_index == pi).collect();
        assert!(!from_variant.is_empty(), "variant {pi} yields candidates");
        for c in &from_variant {
            let (s, e) = (
                c.start_offset_in_protein as u32,
                (c.start_offset_in_protein + c.peptide.residues.len()) as u32,
            );
            let covered = if region.0 < region.1 {
                s < region.1 && e > region.0
            } else {
                s < region.0 && e > region.0
            };
            assert!(covered, "{} does not carry variant {pi}", c.peptide);
            assert!(
                !canonical.contains(&bare(c)),
                "variant peptide {} reads like a canonical one",
                c.peptide
            );
        }
    }
    assert!(cands
        .iter()
        .any(|c| c.protein_index == 1 && c.peptide.to_string() == "R.SDSPANVNSTQR.-"));
    assert!(cands
        .iter()
        .any(|c| c.protein_index == 2 && c.peptide.to_string() == "K.TLSR.S"));
    let n = target.len();
    for pi in [n + 1, n + 2] {
        let v = idx.db.proteins[pi].variant_region().unwrap();
        let mut any = false;
        for c in cands.iter().filter(|c| c.protein_index == pi) {
            any = true;
            let s = c.start_offset_in_protein as u32;
            let e = s + c.peptide.residues.len() as u32;
            assert!(v.covered_by(s, e), "decoy {pi}: {}", c.peptide);
        }
        assert!(any, "decoy variant {pi} yields candidates");
    }
}

#[test]
fn every_expansion_path_agrees_with_annotation() {
    let idx = SearchIndex::from_target_db(&target_db(EvidenceLevel::Curated, true), "XXX");
    let p = params(true);
    let inram: Vec<Candidate> = enumerate_candidates(&idx, &p, "XXX").collect();
    assert!(inram.iter().any(|c| !phospho_positions(c).is_empty()));
    assert_every_path_agrees(&idx, &p);
}

fn assert_every_path_agrees(idx: &SearchIndex, p: &SearchParams) {
    let mi = MmapCandidateIndex::build_in_memory(idx, p, "XXX").unwrap();
    let inram: Vec<Candidate> = enumerate_candidates(idx, p, "XXX").collect();

    let masses: BTreeSet<u64> = inram
        .iter()
        .map(|c| (c.peptide.mass() * 1000.0).round() as u64)
        .collect();
    let tol = 0.02;
    for m in masses {
        let precursor = m as f64 / 1000.0;
        let lazy: HashSet<_> = lazy_candidates_for_precursor(&mi, idx, p, precursor, tol, tol)
            .iter()
            .map(key)
            .collect();
        let want: HashSet<_> = inram
            .iter()
            .filter(|c| (c.peptide.mass() - precursor).abs() <= tol)
            .map(key)
            .collect();
        assert_eq!(lazy, want, "precursor {precursor}");
    }

    let tables = FormWalkTables::new(p);
    let mut all: HashSet<_> = HashSet::new();
    for i in 0..mi.len() {
        let rec = mi.record_at(i);
        let plain: Vec<Candidate> = expand_base_record(idx, p, &rec);
        let mut masses = Vec::new();
        for_each_record_form_masses_bounded(idx, p, &tables, &rec, 0.0, 1e9, |_, m| {
            masses.push(m.iter().sum::<f64>())
        });
        assert_eq!(masses.len(), plain.len(), "record {rec:?}");
        let ks: Vec<u32> = (0..masses.len() as u32).collect();
        let bounded = expand_base_record_bounded(idx, p, &tables, &rec, 0.0, 1e9, &ks);
        let a: HashSet<_> = plain.iter().map(key).collect();
        let b: HashSet<_> = bounded.iter().map(key).collect();
        assert_eq!(a, b, "record {rec:?}");
        all.extend(a);
    }
    let want: HashSet<_> = inram.iter().map(key).collect();
    assert_eq!(all, want, "records expand to the in-RAM candidate set");
}

const NTERM_XML: &str = r#"<uniprot>
<entry dataset="Swiss-Prot">
  <accession>P00003</accession>
  <name>THREE_HUMAN</name>
  <feature type="modified residue" description="N-acetylserine" evidence="1">
    <location><position position="2"/></location>
  </feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="24" version="1">MSAGPEKTLDYRADNPAKVNSTQR</sequence>
</entry>
</uniprot>"#;

#[test]
fn decoys_keep_terminal_sites_at_their_terminus() {
    let target = load_uniprot_xml(NTERM_XML.as_bytes(), &UniProtOptions::default())
        .unwrap()
        .0;
    let acetyl = annotated("Acetyl (Protein N-term)");
    let mut p = SearchParams::default_tryptic(
        AminoAcidSetBuilder::new_standard()
            .add_variable_mod(acetyl)
            .build()
            .unwrap(),
    );
    p.min_length = 4;
    p.max_missed_cleavages = 1;
    for strategy in [
        DecoyStrategy::Reverse,
        DecoyStrategy::Shuffle,
        DecoyStrategy::Sequon,
    ] {
        let idx = SearchIndex::from_target_db_with_strategy(&target, "XXX", None, strategy, 7);
        let decoy = &idx.db.proteins[1];
        assert_eq!(decoy.mod_sites()[0].position, 1, "{strategy:?}");
        assert_eq!(&decoy.sequence[..2], b"MS", "{strategy:?}");
        let acetylated = |decoy: bool| {
            enumerate_candidates(&idx, &p, "XXX")
                .filter(|c| c.is_decoy == decoy)
                .filter(|c| c.peptide.residues[0].mod_.is_some())
                .count()
        };
        assert_eq!(acetylated(false), 2, "{strategy:?}");
        assert!(acetylated(true) > 0, "{strategy:?}: no decoy N-term acetyl");
    }
}

#[test]
fn unrestricted_terminal_mod_is_not_hidden_by_a_restricted_twin() {
    let mut builder = AminoAcidSetBuilder::new_standard();
    builder = builder
        .add_variable_mod(annotated("Acetyl (K)"))
        .add_mods(spec("Acetyl (Protein N-term)", false));
    let mut p = SearchParams::default_tryptic(builder.build().unwrap());
    p.min_length = 4;
    p.max_missed_cleavages = 1;
    let target = ProteinDb {
        proteins: vec![model::Protein {
            accession: "P1".into(),
            description: String::new(),
            sequence: b"KAAPEPTIDEKVLSTAGR".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let forms: Vec<String> = enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy && c.is_protein_n_term)
        .map(|c| c.peptide.to_string())
        .collect();
    assert!(
        forms.iter().any(|f| f.contains("K+42")),
        "protein-N-term acetyl on K missing: {forms:?}"
    );
    assert!(enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .all(|c| c.peptide.residues[1..].iter().all(|aa| aa.mod_.is_none())));
}

#[test]
fn variant_window_walk_equals_filtering_the_whole_protein() {
    use model::{ProteinAnnotation, VariantRegion};
    use std::sync::Arc;
    let seq = b"MAGSPEKTLSYRSDSPAKVNSTQRLLEDKGHIMPQRSTVAEKNNGLKDDCPWYRAHMESKTLIAGFRVEPQK".repeat(3);
    for region in [
        VariantRegion { start: 40, end: 41 },
        VariantRegion {
            start: 100,
            end: 103,
        },
        VariantRegion {
            start: 150,
            end: 150,
        },
        VariantRegion { start: 1, end: 2 },
    ] {
        let annotated = model::Protein {
            accession: "V".into(),
            description: String::new(),
            sequence: seq.clone(),
            annotation: Some(Arc::new(ProteinAnnotation::new(Vec::new(), Some(region)))),
        };
        let plain = model::Protein {
            accession: "P".into(),
            annotation: None,
            ..annotated.clone()
        };
        let db = ProteinDb {
            proteins: vec![annotated, plain],
        };
        let idx =
            SearchIndex::from_target_db_with_strategy(&db, "XXX", None, DecoyStrategy::None, 0);
        for spec in [0u8, 1, 2] {
            let mut p = params(false);
            p.num_tolerable_termini = spec;
            p.max_length = 12;
            let cands: Vec<Candidate> = enumerate_candidates(&idx, &p, "XXX").collect();
            let from = |pi: usize| -> Vec<(usize, String)> {
                cands
                    .iter()
                    .filter(|c| c.protein_index == pi)
                    .filter(|c| {
                        let s = c.start_offset_in_protein as u32;
                        region.covered_by(s, s + c.peptide.residues.len() as u32)
                    })
                    .map(|c| (c.start_offset_in_protein, c.peptide.to_string()))
                    .collect()
            };
            let windowed = from(0);
            assert!(!windowed.is_empty(), "{region:?} ntt={spec}");
            assert_eq!(windowed, from(1), "{region:?} ntt={spec}");
            assert_eq!(
                windowed.len(),
                cands.iter().filter(|c| c.protein_index == 0).count(),
                "{region:?} ntt={spec}: a variant protein yields only covering spans"
            );
        }
    }
}

#[test]
fn index_cache_key_tells_variant_entries_apart() {
    use model::{ProteinAnnotation, VariantRegion};
    use search::candidates::records::index_cache_path;
    use std::sync::Arc;
    let plain = target_db(EvidenceLevel::Experimental, false);
    let mut variant = plain.clone();
    variant.proteins[0].annotation = Some(Arc::new(ProteinAnnotation::new(
        Vec::new(),
        Some(VariantRegion { start: 3, end: 4 }),
    )));
    let p = params(false);
    let dir = std::path::Path::new("cache");
    let a = index_cache_path(&SearchIndex::from_target_db(&plain, "XXX"), &p, dir);
    let b = index_cache_path(&SearchIndex::from_target_db(&variant, "XXX"), &p, dir);
    assert_ne!(a, b);
    let again = index_cache_path(&SearchIndex::from_target_db(&plain, "XXX"), &p, dir);
    assert_eq!(a, again);
}

#[test]
fn variants_in_repeats_add_no_canonical_peptides() {
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
    let target = load_uniprot_xml(xml.as_bytes(), &opts).unwrap().0;
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut p = params(false);
    p.num_tolerable_termini = 0;
    p.min_length = 4;
    p.max_length = 6;
    let cands: Vec<Candidate> = enumerate_candidates(&idx, &p, "XXX").collect();
    let bare = |c: &Candidate| -> String {
        c.peptide
            .residues
            .iter()
            .map(|a| a.residue as char)
            .collect()
    };
    let canonical: HashSet<String> = cands
        .iter()
        .filter(|c| c.protein_index == 0)
        .map(bare)
        .collect();
    let of = |pi: usize| -> HashSet<String> {
        cands
            .iter()
            .filter(|c| c.protein_index == pi)
            .map(bare)
            .collect()
    };
    let (dup, del) = (of(1), of(2));
    assert!(dup.contains("AQQR") && dup.contains("WAQQ"), "{dup:?}");
    assert!(!dup.contains("QRST"), "QRST is canonical");
    assert!(del.contains("KQQR") && del.contains("YKQQRL"), "{del:?}");
    assert!(!del.contains("QQRL"), "QQRL is canonical");
    for v in dup.union(&del) {
        assert!(
            !canonical.contains(v),
            "variant peptide {v} reads canonical"
        );
    }
    let n = target.len();
    for i in [1, 2] {
        let len = target.proteins[i].len() as u32;
        let t = target.proteins[i].variant_region().unwrap();
        let d = idx.db.proteins[n + i].variant_region().unwrap();
        assert_eq!(
            i64::from(d.start) - i64::from(d.end),
            i64::from(t.start) - i64::from(t.end),
            "{t:?} -> {d:?}"
        );
        assert!(d.start.max(d.end) <= len, "{d:?}");
    }
}

const PROCESSED_XML: &str = r#"<uniprot>
<entry dataset="Swiss-Prot">
  <accession>P00010</accession><name>TEN_HUMAN</name>
  <feature type="signal peptide" evidence="1"><location><begin position="1"/><end position="5"/></location></feature>
  <feature type="propeptide" evidence="1"><location><begin position="6"/><end position="8"/></location></feature>
  <feature type="chain" id="PRO_1"><location><begin position="9"/><end position="30"/></location></feature>
  <feature type="peptide" id="PRO_2"><location><begin position="13"/><end position="18"/></location></feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="30">MALWARGRDTHKSEIAHRFKDLGEEHFKAR</sequence>
</entry>
<entry dataset="Swiss-Prot">
  <accession>P00011</accession><name>ELEVEN_HUMAN</name>
  <comment type="alternative products">
    <isoform><id>P00011-1</id><name>1</name><sequence type="displayed"/></isoform>
    <isoform><id>P00011-2</id><name>2</name><sequence type="described" ref="VSP_A VSP_B"/></isoform>
  </comment>
  <feature type="modified residue" description="Phosphothreonine" evidence="1"><location><position position="11"/></location></feature>
  <feature type="splice variant" id="VSP_A"><original>STN</original><variation>Q</variation><location><begin position="2"/><end position="4"/></location></feature>
  <feature type="splice variant" id="VSP_B"><location><begin position="20"/><end position="27"/></location></feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="27">MSTNPKPQRKTKRNTSRRPQDVKFPGG</sequence>
</entry>
<entry dataset="Swiss-Prot">
  <accession>P00012</accession><name>TWELVE_HUMAN</name>
  <feature type="chain" id="PRO_3"><location><begin position="1"/><end position="20"/></location></feature>
  <feature type="propeptide" id="PRO_4"><location><begin position="21"/><end position="26"/></location></feature>
  <sequence length="26">MPEPTIDEKLLGEERAFDYHGKSPAR</sequence>
</entry>
</uniprot>"#;

fn processed_db() -> ProteinDb {
    let opts = UniProtOptions {
        isoforms: true,
        chains: true,
        evidence: EvidenceLevel::Curated,
        ..UniProtOptions::default()
    };
    load_uniprot_xml(PROCESSED_XML.as_bytes(), &opts).unwrap().0
}

#[test]
fn chain_entries_yield_only_their_new_termini() {
    let target = processed_db();
    let accs: Vec<&str> = target
        .proteins
        .iter()
        .map(|p| p.accession.as_str())
        .collect();
    assert_eq!(
        accs,
        [
            "sp|P00010|TEN_HUMAN",
            "sp|P00010_PRO_1|TEN_HUMAN",
            "sp|P00010_PRO_2|TEN_HUMAN",
            "sp|P00011|ELEVEN_HUMAN",
            "sp|P00011-2|ELEVEN_HUMAN",
            "sp|P00012|TWELVE_HUMAN",
            "sp|P00012_PRO_3|TWELVE_HUMAN",
        ]
    );
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let mut builder = AminoAcidSetBuilder::new_standard();
    builder = builder.add_mods(spec("Acetyl (Protein N-term)", false));
    let mut p = SearchParams::default_tryptic(builder.build().unwrap());
    p.min_length = 4;
    p.max_missed_cleavages = 1;
    let cands: Vec<Candidate> = enumerate_candidates(&idx, &p, "XXX").collect();
    let from = |pi: usize| -> Vec<String> {
        cands
            .iter()
            .filter(|c| c.protein_index == pi)
            .map(|c| c.peptide.to_string())
            .collect()
    };
    let mut chain = from(1);
    chain.sort();
    assert_eq!(
        chain,
        [
            "_.D+42.01056THK.S",
            "_.D+42.01056THKSEIAHR.F",
            "_.DTHK.S",
            "_.DTHKSEIAHR.F"
        ],
        "{chain:?}"
    );
    assert!(chain.iter().any(|f| f.contains("+42")), "{chain:?}");
    assert!(cands
        .iter()
        .filter(|c| c.protein_index == 1)
        .all(|c| c.is_protein_n_term && c.start_offset_in_protein == 0));
    let peptide = from(2);
    assert!(!peptide.is_empty());
    assert!(cands
        .iter()
        .filter(|c| c.protein_index == 2)
        .all(|c| c.start_offset_in_protein == 0 && c.is_protein_c_term));
    let mut c_chain = from(6);
    c_chain.sort();
    assert_eq!(c_chain, ["K.LLGEERAFDYH.-", "R.AFDYH.-"], "{c_chain:?}");
    assert!(cands
        .iter()
        .filter(|c| c.protein_index == 6)
        .all(|c| c.is_protein_c_term && !c.is_protein_n_term));
    assert!(from(4).len() > 2, "{:?}", from(4));
    let n = target.len();
    for pi in [n + 1, n + 2, n + 6] {
        let d = &idx.db.proteins[pi];
        assert_eq!(d.chain_ends(), target.proteins[pi - n].chain_ends());
        let ends = d.chain_ends().unwrap();
        for c in cands.iter().filter(|c| c.protein_index == pi) {
            let s = c.start_offset_in_protein as u32;
            let e = s + c.peptide.residues.len() as u32;
            assert!(
                ends.covered_by(s, e, d.len() as u32),
                "decoy {pi}: {}",
                c.peptide
            );
        }
    }
}

#[test]
fn every_expansion_path_agrees_on_isoforms_and_chains() {
    let idx = SearchIndex::from_target_db(&processed_db(), "XXX");
    let mut p = params(false);
    p.min_length = 4;
    p.max_missed_cleavages = 1;
    let mi = MmapCandidateIndex::build_in_memory(&idx, &p, "XXX").unwrap();
    let inram: Vec<Candidate> = enumerate_candidates(&idx, &p, "XXX").collect();
    let want: HashSet<_> = inram.iter().map(key).collect();
    let mut all: HashSet<_> = HashSet::new();
    for i in 0..mi.len() {
        all.extend(
            expand_base_record(&idx, &p, &mi.record_at(i))
                .iter()
                .map(key),
        );
    }
    assert_eq!(all, want);
    let masses: BTreeSet<u64> = inram
        .iter()
        .map(|c| (c.peptide.mass() * 1000.0).round() as u64)
        .collect();
    for m in masses {
        let precursor = m as f64 / 1000.0;
        let lazy: HashSet<_> = lazy_candidates_for_precursor(&mi, &idx, &p, precursor, 0.02, 0.02)
            .iter()
            .map(key)
            .collect();
        let expect: HashSet<_> = inram
            .iter()
            .filter(|c| (c.peptide.mass() - precursor).abs() <= 0.02)
            .map(key)
            .collect();
        assert_eq!(lazy, expect, "precursor {precursor}");
    }
}

fn target_and_decoy_forms(
    target: &ProteinDb,
    p: &SearchParams,
    strategy: DecoyStrategy,
) -> Vec<(Vec<Candidate>, Vec<Candidate>)> {
    let idx = SearchIndex::from_target_db_for(target, "XXX", None, strategy, 7, &Proteases::of(p));
    let cands: Vec<Candidate> = enumerate_candidates(&idx, p, "XXX").collect();
    let n = target.len();
    (0..n)
        .map(|i| {
            let of = |pi: usize| -> Vec<Candidate> {
                cands
                    .iter()
                    .filter(|c| c.protein_index == pi)
                    .cloned()
                    .collect()
            };
            (of(i), of(n + i))
        })
        .collect()
}

const PYRO_XML: &str = r#"<uniprot>
<entry dataset="Swiss-Prot">
  <accession>P00030</accession><name>PYRO_HUMAN</name>
  <feature type="modified residue" description="Pyrrolidone carboxylic acid" evidence="1"><location><position position="1"/></location></feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="23">QAGPEKTLDYRADNPAKVNSTQR</sequence>
</entry>
</uniprot>"#;

#[test]
fn residue_specific_terminal_mods_get_decoys() {
    let target = load_uniprot_xml(PYRO_XML.as_bytes(), &UniProtOptions::default())
        .unwrap()
        .0;
    let pyro = annotated("Gln->pyro-Glu (Protein N-term Q)");
    let p = SearchParams::default_tryptic(
        AminoAcidSetBuilder::new_standard()
            .add_variable_mod(pyro)
            .build()
            .unwrap(),
    );
    let modified = |forms: &[Candidate]| {
        forms
            .iter()
            .filter(|c| c.peptide.residues[0].mod_.is_some())
            .count()
    };
    for strategy in [
        DecoyStrategy::Reverse,
        DecoyStrategy::Shuffle,
        DecoyStrategy::Sequon,
    ] {
        let idx = SearchIndex::from_target_db_with_strategy(&target, "XXX", None, strategy, 7);
        let decoy = &idx.db.proteins[1];
        assert_eq!(decoy.sequence[0], b'Q', "{strategy:?}");
        assert_eq!(decoy.mod_sites()[0].position, 0, "{strategy:?}");
        let forms = target_and_decoy_forms(&target, &p, strategy);
        for (side, forms) in [("target", &forms[0].0), ("decoy", &forms[0].1)] {
            let n_terminal = forms
                .iter()
                .filter(|c| c.start_offset_in_protein == 0)
                .count();
            assert_eq!(modified(forms) * 2, n_terminal, "{strategy:?} {side}");
        }
        assert!(modified(&forms[0].0) > 0, "{strategy:?}");
        if strategy == DecoyStrategy::Reverse {
            assert!(modified(&forms[0].1) > 0, "{strategy:?}: no modified decoy");
        }
    }
}

const CHAIN_DECOY_XML: &str = r#"<uniprot>
<entry dataset="Swiss-Prot">
  <accession>P00040</accession><name>FORTY_HUMAN</name>
  <feature type="chain" id="PRO_40"><location><begin position="5"/><end position="13"/></location></feature>
  <sequence length="13">MSSQAPEPTIDEK</sequence>
</entry>
<entry dataset="Swiss-Prot">
  <accession>P00041</accession><name>FORTYONE_HUMAN</name>
  <feature type="chain" id="PRO_41"><location><begin position="7"/><end position="44"/></location></feature>
  <feature type="chain" id="PRO_42"><location><begin position="1"/><end position="38"/></location></feature>
  <sequence length="50">MKWVTFDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPFDEHVKLRE</sequence>
</entry>
</uniprot>"#;

#[test]
fn chain_decoys_keep_the_terminal_spans() {
    let target = load_uniprot_xml(
        CHAIN_DECOY_XML.as_bytes(),
        &UniProtOptions {
            chains: true,
            evidence: EvidenceLevel::Curated,
            ..UniProtOptions::default()
        },
    )
    .unwrap()
    .0;
    let chains: Vec<usize> = (0..target.len())
        .filter(|&i| target.proteins[i].chain_ends().is_some())
        .collect();
    assert_eq!(chains.len(), 3);

    let base = || AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut p = SearchParams::default_tryptic(base());
    p.min_length = 6;
    p.max_missed_cleavages = 0;
    let forms = target_and_decoy_forms(&target, &p, DecoyStrategy::Reverse);
    let (t, d) = &forms[chains[0]];
    assert_eq!(t.len(), 1);
    assert_eq!(
        d.len(),
        1,
        "{:?}",
        d.iter().map(|c| c.peptide.to_string()).collect::<Vec<_>>()
    );
    assert_eq!(d[0].peptide.to_string(), "_.EDITPEPAK.-");

    let shape = |forms: &[Candidate]| {
        let mut v: Vec<(bool, usize, u64)> = forms
            .iter()
            .map(|c| {
                (
                    c.start_offset_in_protein == 0,
                    c.peptide.residues.len(),
                    (c.peptide.mass() * 1000.0).round() as u64,
                )
            })
            .collect();
        v.sort_unstable();
        v
    };
    for enzyme in [
        Enzyme::Trypsin,
        Enzyme::LysC,
        Enzyme::ArgC,
        Enzyme::GluC,
        Enzyme::AspN,
        Enzyme::Chymotrypsin,
    ] {
        let mut p = SearchParams::default_tryptic(base());
        p.enzyme = enzyme;
        p.min_length = 4;
        p.max_missed_cleavages = 2;
        for strategy in [
            DecoyStrategy::Reverse,
            DecoyStrategy::Shuffle,
            DecoyStrategy::Sequon,
        ] {
            let forms = target_and_decoy_forms(&target, &p, strategy);
            for &i in &chains {
                let (t, d) = &forms[i];
                assert!(!t.is_empty(), "{enzyme:?} {strategy:?} chain {i}");
                assert_eq!(shape(t), shape(d), "{enzyme:?} {strategy:?} chain {i}");
                let seqs = |f: &[Candidate]| -> HashSet<String> {
                    f.iter().map(|c| c.peptide.to_string()).collect()
                };
                assert!(
                    seqs(t).is_disjoint(&seqs(d))
                        || t.iter().all(|c| c.peptide.residues.len() <= 3),
                    "{enzyme:?} {strategy:?} chain {i}: {:?}",
                    seqs(d)
                );
            }
        }
    }
}

#[test]
fn shuffled_variant_parts_keep_their_width() {
    let protein = |region: VariantRegion| Protein {
        accession: "v".to_string(),
        description: String::new(),
        sequence: b"ACDECCCIKNPQSTVWY".to_vec(),
        annotation: Some(std::sync::Arc::new(ProteinAnnotation::new(
            Vec::new(),
            Some(region),
        ))),
    };
    let db = ProteinDb {
        proteins: vec![protein(VariantRegion { start: 4, end: 7 })],
    };
    let shuffled = |db: &ProteinDb, seed: u64| {
        decoy_db(
            db,
            "XXX",
            DecoyStrategy::Shuffle,
            seed,
            &Proteases::default(),
        )
    };
    let decoy = &shuffled(&db, 2).proteins[0];
    let v = decoy.variant_parts()[0];
    assert_eq!((v.end - v.start, v.end <= 17), (3, true), "{v:?}");
    for seed in 0..300 {
        for region in [
            VariantRegion { start: 4, end: 7 },
            VariantRegion { start: 0, end: 5 },
            VariantRegion { start: 12, end: 17 },
            VariantRegion { start: 9, end: 9 },
            VariantRegion { start: 11, end: 8 },
        ] {
            let db = ProteinDb {
                proteins: vec![protein(region)],
            };
            let d = shuffled(&db, seed).proteins[0].variant_parts()[0];
            if region.start < region.end {
                assert_eq!(
                    d.end - d.start,
                    region.end - region.start,
                    "{seed} {region:?}"
                );
                assert!(d.end <= 17, "{seed} {region:?} -> {d:?}");
            } else {
                assert_eq!(
                    d.start + region.end,
                    d.end + region.start,
                    "{seed} {region:?}"
                );
                assert!(d.end >= 1 && d.start < 17, "{seed} {region:?} -> {d:?}");
            }
        }
    }
}

const MULTI_XML: &str = r#"<uniprot>
<entry dataset="Swiss-Prot">
  <accession>P00050</accession><name>FIFTY_HUMAN</name>
  <feature type="sequence variant" id="VAR_000050" evidence="1"><original>KRAPEPTIDEKR</original><variation>ERAPEPTIDEKE</variation><location><begin position="2"/><end position="13"/></location></feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="15">MKRAPEPTIDEKRGH</sequence>
</entry>
</uniprot>"#;

#[test]
fn variant_peptides_carry_a_change_not_just_its_neighbourhood() {
    let target = load_uniprot_xml(
        MULTI_XML.as_bytes(),
        &UniProtOptions {
            variants: true,
            ..UniProtOptions::default()
        },
    )
    .unwrap()
    .0;
    assert_eq!(target.len(), 2);
    assert_eq!(
        target.proteins[1].variant_parts(),
        [
            VariantRegion { start: 1, end: 2 },
            VariantRegion { start: 12, end: 13 }
        ]
    );
    let base = || AminoAcidSetBuilder::new_standard().build().unwrap();
    let from_variant = |missed: u32| -> Vec<String> {
        let mut p = SearchParams::default_tryptic(base());
        p.min_length = 4;
        p.max_missed_cleavages = missed;
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let mut v: Vec<String> = enumerate_candidates(&idx, &p, "XXX")
            .filter(|c| c.protein_index == 1)
            .map(|c| c.peptide.to_string())
            .collect();
        v.sort();
        v
    };
    assert!(from_variant(0).is_empty(), "{:?}", from_variant(0));
    let mut with_one: Vec<String> = from_variant(1);
    with_one.dedup();
    assert_eq!(
        with_one,
        ["M.ERAPEPTIDEK.E", "R.APEPTIDEKEGH.-", "_.MERAPEPTIDEK.E"]
    );
    for strategy in [
        DecoyStrategy::Reverse,
        DecoyStrategy::Shuffle,
        DecoyStrategy::Sequon,
    ] {
        let idx = SearchIndex::from_target_db_with_strategy(&target, "XXX", None, strategy, 7);
        let widths = |p: &Protein| -> Vec<u32> {
            let mut w: Vec<u32> = p.variant_parts().iter().map(|v| v.end - v.start).collect();
            w.sort_unstable();
            w
        };
        assert_eq!(widths(&idx.db.proteins[3]), [1, 1], "{strategy:?}");
    }
}

#[test]
fn origins_follow_entries_onto_decoys() {
    use model::EntryKind::*;
    let target = processed_db();
    let kinds: Vec<model::EntryKind> = target.proteins.iter().map(|p| p.entry_kind()).collect();
    assert_eq!(
        kinds,
        [Reference, Chain, Chain, Reference, Isoform, Reference, Chain]
    );
    let iso = target.proteins[4].origin().unwrap();
    assert_eq!((&*iso.family, &*iso.feature), ("P00011", "P00011-2"));
    let chain = target.proteins[1].origin().unwrap();
    assert_eq!((&*chain.family, &*chain.feature), ("P00010", "PRO_1"));
    for strategy in [
        DecoyStrategy::Reverse,
        DecoyStrategy::Shuffle,
        DecoyStrategy::Sequon,
    ] {
        let idx = SearchIndex::from_target_db_with_strategy(&target, "XXX", None, strategy, 7);
        let n = target.len();
        for i in 0..n {
            assert_eq!(
                idx.db.proteins[n + i].origin(),
                target.proteins[i].origin(),
                "{strategy:?} {i}"
            );
        }
    }
}

#[test]
fn origin_alone_leaves_the_reverse_decoy_unchanged() {
    let target = processed_db();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let plain = SearchIndex::from_target_db(
        &ProteinDb {
            proteins: target
                .proteins
                .iter()
                .map(|p| Protein {
                    annotation: None,
                    ..p.clone()
                })
                .collect(),
        },
        "XXX",
    );
    let n = target.len();
    let mut checked = 0;
    for (i, p) in target.proteins.iter().enumerate() {
        let terminal_site = p
            .mod_sites()
            .iter()
            .any(|s| s.anchor != model::SiteAnchor::Residue);
        if p.chain_ends().is_some() || terminal_site {
            continue;
        }
        let want = &plain.db.proteins[n + i].sequence;
        assert_ne!(want, &p.sequence, "{}", p.accession);
        assert_eq!(&idx.db.proteins[n + i].sequence, want, "{}", p.accession);
        checked += 1;
    }
    assert!(target.proteins[4].mod_sites().is_empty() || checked >= 1);
    assert!(checked >= 3, "{checked}");
}

#[test]
fn index_knows_whether_it_holds_alternative_entries() {
    let plain = SearchIndex::from_target_db(&target_db(EvidenceLevel::Experimental, false), "XXX");
    assert!(!plain.has_alternative_entries);
    let with = SearchIndex::from_target_db(&target_db(EvidenceLevel::Experimental, true), "XXX");
    assert!(with.has_alternative_entries);
}

const DEAMIDATION_XML: &str = r#"<uniprot xmlns="http://uniprot.org/uniprot">
<entry dataset="Swiss-Prot">
  <accession>P00030</accession><name>DEAM_HUMAN</name>
  <feature type="sequence variant" id="VAR_000030" evidence="1">
    <original>N</original><variation>D</variation>
    <location><position position="10"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000031" evidence="1">
    <original>D</original><variation>N</variation>
    <location><position position="18"/></location>
  </feature>
  <evidence type="ECO:0000269" key="1"/>
  <sequence length="22" version="1">AGSPEAKVLNGSTRAPEDLTGR</sequence>
</entry>
</uniprot>"#;

fn deamidation_params(deamidation: bool) -> SearchParams {
    let mut builder = AminoAcidSetBuilder::new_standard();
    if deamidation {
        builder = builder.add_mods(spec("Deamidated (N)", false));
    }
    let mut p = SearchParams::default_tryptic(builder.build().unwrap());
    p.min_length = 4;
    p.max_missed_cleavages = 0;
    p
}

#[test]
fn variant_forms_a_searched_mod_makes_canonical_are_left_out() {
    let opts = UniProtOptions {
        variants: true,
        ..UniProtOptions::default()
    };
    let target = load_uniprot_xml(DEAMIDATION_XML.as_bytes(), &opts)
        .unwrap()
        .0;
    assert_eq!(target.len(), 3);
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let variant_forms = |p: &SearchParams, decoy: bool| -> Vec<String> {
        let mut v: Vec<String> = enumerate_candidates(&idx, p, "XXX")
            .filter(|c| c.is_decoy == decoy)
            .filter(|c| idx.db.proteins[c.protein_index].variant_region().is_some())
            .map(|c| c.peptide.to_string())
            .collect();
        v.sort();
        v
    };
    let plain = variant_forms(&deamidation_params(false), false);
    assert!(plain.iter().any(|f| f.contains("VLDGSTR")), "{plain:?}");
    assert!(plain.iter().any(|f| f.contains("APENLTGR")), "{plain:?}");
    let deam = deamidation_params(true);
    let forms = variant_forms(&deam, false);
    assert!(!forms.iter().any(|f| f.contains("VLD")), "{forms:?}");
    assert!(forms.iter().any(|f| f.contains("APENLTGR")), "{forms:?}");
    assert!(!forms.iter().any(|f| f.contains("APEN+0.98")), "{forms:?}");
    assert_eq!(forms.len(), variant_forms(&deam, true).len());
    assert_every_path_agrees(&idx, &deam);
}
