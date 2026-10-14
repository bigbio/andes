//! Candidate generation over a UniProt XML database: sequence-variant entries
//! yield only the peptides that carry the variant, mods declared
//! `sites=annotated` land only on annotated residues, decoys mirror both, and
//! every expansion path (in-RAM, lazy precursor window, fragment-index walk)
//! produces the same candidate set.

use std::collections::{BTreeSet, HashSet};

use input::{load_uniprot_xml, EvidenceLevel, UniProtOptions};
use model::{AminoAcidSetBuilder, Modification, ProteinDb};
use search::candidate_gen::{
    expand_base_record, expand_base_record_bounded, for_each_record_form_masses_bounded,
    lazy_candidates_for_precursor, Candidate, FormWalkTables,
};
use search::candidate_index::MmapCandidateIndex;
use search::decoy::DecoyStrategy;
use search::{enumerate_candidates, SearchIndex, SearchParams};

/// P1 (tryptic peptides AGSPEK | TLSYR | SDSPAK | VNSTQR, 0-based S2, S8, S11,
/// S13, S19, T6, T20) is phosphorylated at S2 (experimental), S13
/// (experimental) and T20 (by similarity), and has a natural variant K16N
/// (removes a cleavage site), a deletion of Y9, and an unannotated twin P2.
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
    load_uniprot_xml(XML.as_bytes(), &UniProtOptions { variants, evidence })
        .unwrap()
        .0
}

/// Phospho on S and T: restricted to annotated sites or not.
fn params(restricted: bool) -> SearchParams {
    let sites = if restricted { ",sites=annotated" } else { "" };
    let mut builder = AminoAcidSetBuilder::new_standard();
    for r in ["S", "T"] {
        let line = format!("79.966331,{r},opt,any,Phospho{sites}");
        builder = builder.add_variable_mod(Modification::from_mods_txt_line(&line).unwrap());
    }
    let mut p = SearchParams::default_tryptic(builder.build().unwrap());
    p.min_length = 4;
    p.max_length = 30;
    p.max_missed_cleavages = 1;
    p.max_variable_mods_per_peptide = 2;
    p
}

/// 0-based protein positions carrying a phospho in `c`.
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
    // Experimental evidence: S2 and S13 of P1 only; T20 is "by similarity",
    // and P2 has no annotation.
    let expect: BTreeSet<(String, usize)> = [2, 13]
        .into_iter()
        .map(|p| ("sp|P00001|ONE_HUMAN".to_string(), p))
        .collect();
    assert_eq!(placed, expect);

    // Curated evidence admits T20 too.
    let idx = SearchIndex::from_target_db(&target_db(EvidenceLevel::Curated, false), "XXX");
    let curated: BTreeSet<usize> = enumerate_candidates(&idx, &params(true), "XXX")
        .filter(|c| !c.is_decoy)
        .flat_map(|c| phospho_positions(&c))
        .collect();
    assert_eq!(curated, BTreeSet::from([2, 13, 20]));

    // Unrestricted phospho goes on every S/T of both proteins.
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
        DecoyStrategy::SequonReverse,
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
        // Restricted phospho on decoys only where the decoy annotation says.
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
    // Reversal mirrors exactly: S2/S13/T20 of a 23-mer → 20/9/2.
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let decoy = &idx.db.proteins[target.len()];
    assert_eq!(decoy.accession, "XXX_sp|P00001|ONE_HUMAN");
    let pos: Vec<u32> = decoy.mod_sites().iter().map(|s| s.position).collect();
    assert_eq!(pos, [2, 9, 20]);
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
    // K17N removes the cleavage after K17: SDSPANVNSTQR is the variant peptide.
    assert!(cands
        .iter()
        .any(|c| c.protein_index == 1 && c.peptide.to_string() == "R.SDSPANVNSTQR.-"));
    // The Y10 deletion: TLSR.
    assert!(cands
        .iter()
        .any(|c| c.protein_index == 2 && c.peptide.to_string() == "K.TLSR.S"));
    // Decoy variants yield only spans covering their mirrored variant.
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
    let mi = MmapCandidateIndex::build_in_memory(&idx, &p, "XXX").unwrap();
    let inram: Vec<Candidate> = enumerate_candidates(&idx, &p, "XXX").collect();
    assert!(inram.iter().any(|c| !phospho_positions(c).is_empty()));

    // Lazy precursor window == mass-filtered in-RAM set, at every candidate mass.
    let masses: BTreeSet<u64> = inram
        .iter()
        .map(|c| (c.peptide.mass() * 1000.0).round() as u64)
        .collect();
    let tol = 0.02;
    for m in masses {
        let precursor = m as f64 / 1000.0;
        let lazy: HashSet<_> = lazy_candidates_for_precursor(&mi, &idx, &p, precursor, tol, tol)
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

    // Per record: the plain expansion, the bounded walk the fragment index is
    // built from, and its materialisation see the same forms.
    let tables = FormWalkTables::new(&p);
    let mut all: HashSet<_> = HashSet::new();
    for i in 0..mi.len() {
        let rec = mi.record_at(i);
        let plain: Vec<Candidate> = expand_base_record(&idx, &p, &rec);
        let mut masses = Vec::new();
        for_each_record_form_masses_bounded(&idx, &p, &tables, &rec, 0.0, 1e9, |_, m| {
            masses.push(m.iter().sum::<f64>())
        });
        assert_eq!(masses.len(), plain.len(), "record {rec:?}");
        let ks: Vec<u32> = (0..masses.len() as u32).collect();
        let bounded = expand_base_record_bounded(&idx, &p, &tables, &rec, 0.0, 1e9, &ks);
        let a: HashSet<_> = plain.iter().map(key).collect();
        let b: HashSet<_> = bounded.iter().map(key).collect();
        assert_eq!(a, b, "record {rec:?}");
        all.extend(a);
    }
    let want: HashSet<_> = inram.iter().map(key).collect();
    assert_eq!(all, want, "records expand to the in-RAM candidate set");
}

/// A protein starting `MS…` whose Ser2 is N-acetylated after initiator-Met
/// removal (an N-terminal site), searched with a protein-N-term Acetyl
/// restricted to annotated sites.
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
    let acetyl =
        Modification::from_mods_txt_line("42.010565,*,opt,Prot-N-term,Acetyl,sites=annotated")
            .unwrap();
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
        DecoyStrategy::SequonReverse,
    ] {
        let idx = SearchIndex::from_target_db_with_strategy(&target, "XXX", None, strategy, 7);
        assert_eq!(
            idx.db.proteins[1].mod_sites()[0].position,
            0,
            "{strategy:?}"
        );
        let acetylated = |decoy: bool| {
            enumerate_candidates(&idx, &p, "XXX")
                .filter(|c| c.is_decoy == decoy)
                .filter(|c| c.peptide.residues[0].mod_.is_some())
                .count()
        };
        // SAGPEK after Met removal, and SAGPEKTLDYR with a missed cleavage.
        assert_eq!(acetylated(false), 2, "{strategy:?}");
        // Following its residue, the site would sit deep inside the decoy,
        // where a protein-N-term mod never applies.
        assert!(acetylated(true) > 0, "{strategy:?}: no decoy N-term acetyl");
    }
}

#[test]
fn unrestricted_terminal_mod_is_not_hidden_by_a_restricted_twin() {
    // Acetyl-K restricted to annotated sites, protein-N-term Acetyl on any
    // residue unrestricted: an unannotated protein starting with K must still
    // get its protein-N-term acetyl on that K.
    let mut builder = AminoAcidSetBuilder::new_standard();
    for line in [
        "42.010565,K,opt,any,Acetyl,sites=annotated",
        "42.010565,*,opt,Prot-N-term,Acetyl",
    ] {
        builder = builder.add_variable_mod(Modification::from_mods_txt_line(line).unwrap());
    }
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
    // The internal K11 has no annotated site, so no Acetyl there.
    assert!(enumerate_candidates(&idx, &p, "XXX")
        .filter(|c| !c.is_decoy)
        .all(|c| c.peptide.residues[1..].iter().all(|aa| aa.mod_.is_none())));
}

#[test]
fn variant_window_walk_equals_filtering_the_whole_protein() {
    use model::{ProteinAnnotation, VariantRegion};
    use std::sync::Arc;
    // Long enough that the window around the variants excludes most of it.
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
    use search::candidate_index::index_cache_path;
    use std::sync::Arc;
    let plain = target_db(EvidenceLevel::Experimental, false);
    let mut variant = plain.clone();
    variant.proteins[0].annotation = Some(Arc::new(ProteinAnnotation::new(
        Vec::new(),
        Some(VariantRegion { start: 3, end: 4 }),
    )));
    let p = params(false);
    let a = index_cache_path(&SearchIndex::from_target_db(&plain, "XXX"), &p);
    let b = index_cache_path(&SearchIndex::from_target_db(&variant, "XXX"), &p);
    assert_ne!(a, b);
    // Annotated sites alone do not change which base records exist.
    let again = index_cache_path(&SearchIndex::from_target_db(&plain, "XXX"), &p);
    assert_eq!(a, again);
}
