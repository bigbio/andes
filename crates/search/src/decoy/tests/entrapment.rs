use super::*;
use crate::decoy::collisions::for_each_span;
use model::protein::Protein;

fn db(proteins: &[(&str, &[u8])]) -> ProteinDb {
    ProteinDb {
        proteins: proteins
            .iter()
            .map(|(acc, seq)| Protein {
                accession: acc.to_string(),
                description: String::new(),
                sequence: seq.to_vec(),
                annotation: None,
            })
            .collect(),
    }
}

fn targets() -> ProteinDb {
    db(&[
        ("P1", b"MAGLEVTDPLSKWQEFINRTYHGASSPLDEK"),
        ("P2", b"MSTPVLEGHWQKAFDNYIRGLTSPEVCHMK"),
        ("P3", b"QWERTYIPASDFGHKLCVNMQWERTYR"),
    ])
}

fn peptides(seq: &[u8], d: &Proteases) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    for_each_span(seq, seq.first() == Some(&b'M'), d, |s, e| {
        out.push(seq[s..e].to_vec())
    });
    out
}

fn sorted(mut v: Vec<u8>) -> Vec<u8> {
    v.sort_unstable();
    v
}

#[test]
fn every_entrapment_peptide_has_a_target_peptide_of_the_same_composition_at_its_place() {
    let t = targets();
    let d = Proteases::default();
    let e = entrapment_db(&t, DEFAULT_ENTRAPMENT_SEED, &d).expect("tryptic digest");
    assert_eq!(e.proteins.proteins.len(), t.proteins.len());
    for (target, ent) in t.proteins.iter().zip(&e.proteins.proteins) {
        assert_eq!(ent.accession, format!("ENTRAP_{}", target.accession));
        assert_eq!(ent.sequence.len(), target.sequence.len());
        let (tp, ep) = (peptides(&target.sequence, &d), peptides(&ent.sequence, &d));
        assert_eq!(tp.len(), ep.len());
        for (a, b) in tp.iter().zip(&ep) {
            assert_eq!(sorted(a.clone()), sorted(b.clone()));
            assert_eq!(a.last(), b.last());
        }
    }
}

#[test]
fn entrapment_peptides_do_not_spell_target_peptides() {
    let t = targets();
    let d = Proteases::default();
    let e = entrapment_db(&t, DEFAULT_ENTRAPMENT_SEED, &d).expect("tryptic digest");
    assert_eq!(e.stats.left, 0);
    for ent in &e.proteins.proteins {
        for p in peptides(&ent.sequence, &d) {
            let targets = TargetPeptides::union(&t, std::slice::from_ref(&d));
            assert!(!targets.contains(&p), "{}", String::from_utf8_lossy(&p));
        }
    }
    assert!(e.stats.entrapment_peptides > 0);
    assert!(e.stats.peptide_ratio() > 0.9 && e.stats.peptide_ratio() <= 1.0);
}

#[test]
fn the_entrapment_shuffle_differs_from_the_decoy_shuffle_and_is_seeded() {
    let t = targets();
    let d = Proteases::default();
    let e1 = entrapment_db(&t, DEFAULT_ENTRAPMENT_SEED, &d).expect("tryptic digest");
    let again = entrapment_db(&t, DEFAULT_ENTRAPMENT_SEED, &d).expect("tryptic digest");
    let e2 = entrapment_db(&t, DEFAULT_ENTRAPMENT_SEED + 1, &d).expect("tryptic digest");
    let decoys = crate::decoy::decoy_db(
        &t,
        "XXX",
        DecoyStrategy::Shuffle,
        crate::decoy::DEFAULT_DECOY_SEED,
        &d,
    );
    let seqs = |db: &ProteinDb| -> Vec<Vec<u8>> {
        db.proteins.iter().map(|p| p.sequence.clone()).collect()
    };
    assert_eq!(seqs(&e1.proteins), seqs(&again.proteins));
    assert_ne!(seqs(&e1.proteins), seqs(&e2.proteins));
    assert_ne!(seqs(&e1.proteins), seqs(&decoys));
}

#[test]
fn a_non_specific_digest_has_no_entrapment() {
    let d = Proteases {
        primary: model::enzyme::Enzyme::NonSpecific,
        num_tolerable_termini: 0,
        ..Proteases::default()
    };
    assert!(entrapment_db(&targets(), DEFAULT_ENTRAPMENT_SEED, &d).is_none());
}

#[test]
fn repeated_entrapment_peptides_with_different_sources_have_no_pair() {
    let searched = db(&[("P1", b"ABCDEFKGHIJLMK"), ("ENTRAP_P1", b"FEDCBAKFEDCBAK")]);
    let lookup = EntrapmentLookup::new(searched.proteins.iter(), &Proteases::default())
        .expect("entrapment entries");
    assert_eq!(lookup.paired_target("ENTRAP_P1", b"FEDCBAK"), None);
}

#[test]
fn identical_target_peptides_across_species_are_counted_once() {
    let searched = db(&[
        ("human_P1", b"MAGLEVTDPLSK"),
        ("yeast_P1", b"MAGLEVTDPLSK"),
        ("ENTRAP_human_P1", b"MPLGAVDTELSK"),
    ]);
    let d = Proteases::default();
    let lookup = EntrapmentLookup::new(searched.proteins.iter(), &d).expect("entrapment");
    let one_species = db(&[("P1", b"MAGLEVTDPLSK")]);
    let unique = TargetPeptides::union(&one_species, std::slice::from_ref(&d));
    assert_eq!(lookup.target_peptides, unique.len());
    assert!(lookup.is_target_peptide(b"MAGLEVTDPLSK"));
    assert!(lookup.is_target_peptide(b"MAGIEVTDPLSK"));
}

#[test]
fn entrapment_accessions_name_their_source() {
    assert!(is_entrapment_accession("ENTRAP_sp|P1|A_HUMAN"));
    assert!(!is_entrapment_accession("ENTRAPX_P1"));
    assert!(!is_entrapment_accession("XXX_ENTRAP_P1"));
    assert_eq!(entrapment_source("ENTRAP_sp|P1|A"), Some("sp|P1|A"));
    assert_eq!(entrapment_source("sp|P1|A"), None);
}
