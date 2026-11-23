use super::*;

fn make_db(proteins: &[(&str, &[u8])]) -> ProteinDb {
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

#[test]
fn foreign_decoy_tag_is_the_whole_first_or_last_field() {
    assert_eq!(
        strip_foreign_decoy_tag("rev_sp|P1|A_HUMAN"),
        Some("sp|P1|A_HUMAN")
    );
    assert_eq!(strip_foreign_decoy_tag("REV__P1"), Some("P1"));
    assert_eq!(strip_foreign_decoy_tag("DECOY_P1"), Some("P1"));
    assert_eq!(strip_foreign_decoy_tag("P1_rev"), Some("P1"));
    assert_eq!(strip_foreign_decoy_tag("sp|P04618|REV_HV1H2"), None);
    assert_eq!(strip_foreign_decoy_tag("sp|Q9UBZ9|REV1_HUMAN"), None);
    assert_eq!(strip_foreign_decoy_tag("rev_"), None);
    assert_eq!(strip_foreign_decoy_tag("P1"), None);
}

#[test]
fn a_target_only_database_has_no_foreign_decoys() {
    let db = make_db(&[
        ("sp|P1|A_HUMAN", b"MKTAYIAKQR"),
        ("sp|P2|B_HUMAN", b"MSEQNNTEMTFQIQR"),
        ("sp|P3|C_HUMAN", b"MKAKM"),
    ]);
    let found = ForeignDecoys::scan(&db);
    assert_eq!(found.n_reversed, 0);
    assert_eq!(found.n_named_paired + found.n_named_unpaired, 0);
    assert_eq!(found.verdict(), Ok(None));
}

#[test]
fn reversed_records_under_an_unknown_name_are_an_error() {
    let db = make_db(&[
        ("sp|P1|A_HUMAN", b"MKTAYIAKQR"),
        ("sp|P2|B_HUMAN", b"MSEQNNTEMTFQIQR"),
        ("zz1", b"RQKAIYATKM"),
        ("zz2", b"RQIQFTMETNNQESM"),
    ]);
    let found = ForeignDecoys::scan(&db);
    assert_eq!(found.n_reversed, 4);
    assert_eq!(found.n_named_paired, 0);
    let err = found.verdict().unwrap_err();
    assert!(
        err.contains("4 of 4 records are the exact reverse"),
        "{err}"
    );
    assert!(err.contains("e.g. zz2"), "{err}");
}

#[test]
fn short_reversed_sequences_are_not_evidence() {
    assert_eq!(min_reversible_len(4), 1);
    assert_eq!(min_reversible_len(20), 3);
    assert_eq!(min_reversible_len(63_778), 8);
    let filler: Vec<(String, Vec<u8>)> = (0..18u8)
        .map(|i| (format!("p{i}"), vec![b'M', b'A' + i, b'G', b'K', b'W']))
        .collect();
    let mut records: Vec<(&str, &[u8])> =
        filler.iter().map(|(a, s)| (&a[..], &s[..])).collect();
    records.push(("d1", b"EI"));
    records.push(("d2", b"IE"));
    records.push(("t1", b"MKTAYIAKQR"));
    records.push(("t2", b"MSEQNNTEMTFQIQR"));
    assert_eq!(ForeignDecoys::scan(&make_db(&records)).n_reversed, 0);
    *records.last_mut().unwrap() = ("t2", &b"RQKAIYATKM"[..]);
    assert_eq!(ForeignDecoys::scan(&make_db(&records)).n_reversed, 2);
}

#[test]
fn a_record_named_as_the_decoy_of_another_is_an_error() {
    let db = make_db(&[
        ("sp|P1|A_HUMAN", b"MKTAYIAKQR"),
        ("sp|P2|B_HUMAN", b"MSEQNNTEMTFQIQR"),
        ("DECOY_sp|P1|A_HUMAN", b"MAKYTIKAQR"),
    ]);
    let found = ForeignDecoys::scan(&db);
    assert_eq!(found.n_reversed, 0);
    assert_eq!(found.n_named_paired, 1);
    let err = found.verdict().unwrap_err();
    assert!(
        err.contains("1 of 3 records are named as the decoy"),
        "{err}"
    );
    assert!(err.contains("DECOY_sp|P1|A_HUMAN"), "{err}");
}

#[test]
fn weaker_evidence_is_a_warning() {
    let db = make_db(&[
        ("p1", b"MKTAYIAKQR"),
        ("p2", b"AKM"),
        ("p3", b"MKA"),
        ("p4", b"MSEQNNTEMTFQIQR"),
        ("p5", b"MLLPEPTIDEK"),
        ("p6", b"MGGSSTTWWYK"),
        ("REV_HV1H2", b"MAGRSGDSDEELLK"),
    ]);
    let found = ForeignDecoys::scan(&db);
    assert_eq!(found.n_reversed, 2);
    assert_eq!(found.n_named_unpaired, 1);
    let warning = found.verdict().unwrap().unwrap();
    assert!(
        warning.starts_with("WARN: the database may contain decoys"),
        "{warning}"
    );
    assert!(warning.contains("REV_HV1H2"), "{warning}");
}

fn decoys(db: &ProteinDb, strategy: DecoyStrategy, proteases: &Proteases) -> Vec<Vec<u8>> {
    decoy_db(db, "XXX", strategy, DEFAULT_DECOY_SEED, proteases)
        .proteins
        .into_iter()
        .map(|p| p.sequence)
        .collect()
}

fn reversed(db: &ProteinDb) -> Vec<Vec<u8>> {
    decoys(db, DecoyStrategy::Reverse, &Proteases::default())
}

fn enzyme(primary: Enzyme) -> Proteases {
    Proteases {
        primary,
        ..Proteases::default()
    }
}

#[test]
fn mutate_changes_one_interior_residue_per_long_segment() {
    let db = make_db(&[("P1", b"PEPTIDELSAGRVLSK"), ("P2", b"MAGCTFYHRQGILW")]);
    let d = decoys(&db, DecoyStrategy::Mutate, &Proteases::default());
    let mass = |r: u8| model::AminoAcid::standard(r).unwrap().mass;
    for (t, s) in db.proteins.iter().zip(&d) {
        assert_eq!(t.sequence.len(), s.len());
        let cuts = compute_cleavage_positions(&t.sequence, Enzyme::Trypsin, &[]);
        for w in cuts.windows(2) {
            let (a, b) = (w[0] as usize, w[1] as usize);
            let diff: Vec<usize> = (a..b).filter(|&i| t.sequence[i] != s[i]).collect();
            if b - a < MUTATE_MIN_SEGMENT {
                assert!(diff.is_empty());
                continue;
            }
            assert_eq!(diff.len(), 1, "{:?}", std::str::from_utf8(s));
            let i = diff[0];
            assert!(i > a && i + 1 < b);
            assert!(!(t.sequence[0] == b'M' && a == 0 && i == 1));
            assert!(!b"KRP".contains(&s[i]));
            assert!((mass(s[i]) - mass(t.sequence[i])).abs() > 1.1);
        }
    }
    assert_eq!(d, decoys(&db, DecoyStrategy::Mutate, &Proteases::default()));
    assert_eq!(
        DecoyStrategy::from_name("mutate"),
        Some(DecoyStrategy::Mutate)
    );
}

#[test]
fn mutate_gives_a_shared_segment_one_decoy() {
    let seg = b"AGLSTVEFK";
    let site = |position: u32| ModSite {
        position,
        mass_delta: 79.96633,
        anchor: SiteAnchor::Residue,
    };
    let entry = |acc: &str, open: usize| Protein {
        accession: acc.to_string(),
        description: String::new(),
        sequence: [b"WWK".as_slice(), seg].concat(),
        annotation: Some(Arc::new(ProteinAnnotation::new(
            (1..seg.len() - 1)
                .filter(|&i| i != open && i != 4)
                .map(|i| site(3 + i as u32))
                .collect(),
            None,
        ))),
    };
    let db = ProteinDb {
        proteins: vec![entry("P1", 2), entry("P2", 6)],
    };
    let d = decoys(&db, DecoyStrategy::Mutate, &Proteases::default());
    let tail = |s: &[u8]| s[3..].to_vec();
    assert_eq!(tail(&d[0]), tail(&d[1]));
    let changed: Vec<usize> = (0..seg.len())
        .filter(|&i| tail(&d[0])[i] != seg[i])
        .collect();
    assert_eq!(changed, vec![4]);
}

#[test]
fn refinement_mutate_decoys_are_reversed() {
    let db = make_db(&[("P1", b"PEPKTIDE"), ("P2", b"ACDEFGHIK")]);
    let d = whole_sequence_decoy_db(&db, "XXX", DecoyStrategy::Mutate, 7);
    let r = whole_sequence_decoy_db(&db, "XXX", DecoyStrategy::Reverse, 7);
    let seqs = |db: &ProteinDb| -> Vec<Vec<u8>> {
        db.proteins.iter().map(|p| p.sequence.clone()).collect()
    };
    assert_eq!(seqs(&d), seqs(&r));
    for (t, s) in db.proteins.iter().zip(&d.proteins) {
        assert_ne!(t.sequence, s.sequence);
    }
}

#[test]
fn pseudo_reverses_each_tryptic_segment() {
    let db = make_db(&[("P1", b"PEPTIDEKVLSAGR"), ("P2", b"AGCTFYHRQGIL")]);
    assert_eq!(
        reversed(&db),
        [b"EDITPEPKGASLVR".to_vec(), b"HYFTCGARIGQL".to_vec()]
    );
}

#[test]
fn initiator_met_stays() {
    let db = make_db(&[("P1", b"MAGPEKR")]);
    assert_eq!(reversed(&db), [b"MEPGAKR".to_vec()]);
}

#[test]
fn n_terminal_enzyme_keeps_the_first_residue_of_each_segment() {
    let db = make_db(&[("P1", b"DAAGLKDPEW")]);
    let d = decoys(&db, DecoyStrategy::Reverse, &enzyme(Enzyme::AspN));
    assert_eq!(d, [b"DKLGAADWEP".to_vec()]);
}

#[test]
fn non_specific_search_reverses_the_whole_sequence() {
    let db = make_db(&[("P1", b"MKWVLR")]);
    let want = [b"RLVWKM".to_vec()];
    assert_eq!(
        decoys(&db, DecoyStrategy::Reverse, &enzyme(Enzyme::NonSpecific)),
        want
    );
    let ntt0 = Proteases {
        num_tolerable_termini: 0,
        ..Proteases::default()
    };
    assert_eq!(decoys(&db, DecoyStrategy::Reverse, &ntt0), want);
    let semi = Proteases {
        num_tolerable_termini: 1,
        ..Proteases::default()
    };
    assert_eq!(
        decoys(&db, DecoyStrategy::Reverse, &semi),
        [b"MKLVWR".to_vec()]
    );
}

#[test]
fn a_segment_has_one_decoy_wherever_it_occurs() {
    let db = make_db(&[("P1", b"GGHWKLLSAGNDKPPR"), ("P2", b"LLSAGNDKAEFRLLSAGNDK")]);
    for strategy in [
        DecoyStrategy::Reverse,
        DecoyStrategy::Shuffle,
        DecoyStrategy::Sequon,
    ] {
        let d = decoys(&db, strategy, &Proteases::default());
        let seg = d[0][5..13].to_vec();
        assert_eq!(seg[7], b'K', "{strategy:?}");
        assert_eq!(d[1][0..8], seg[..], "{strategy:?}");
        assert_eq!(d[1][12..20], seg[..], "{strategy:?}");
    }
}

#[test]
fn palindromic_segment_is_redrawn() {
    let db = make_db(&[("P1", b"ADAKGLVLGR")]);
    let d = &reversed(&db)[0];
    assert_ne!(&d[0..4], b"ADAK");
    assert_eq!(d[3], b'K');
    let mut seg = d[0..4].to_vec();
    seg.sort_unstable();
    assert_eq!(seg, b"AADK");
    assert_ne!(&d[4..10], b"GLVLGR");
    assert_eq!(d[9], b'R');
}

#[test]
fn decoy_segment_never_spells_another_target_segment() {
    let db = make_db(&[
        ("P1", b"GASVK"),
        ("P2", b"VSAGK"),
        ("P3", b"LIAGK"),
        ("P4", b"GALLK"),
    ]);
    let d = reversed(&db);
    let targets: Vec<Vec<u8>> = db
        .proteins
        .iter()
        .map(|p| {
            p.sequence
                .iter()
                .map(|&r| if r == b'I' { b'L' } else { r })
                .collect()
        })
        .collect();
    for (t, s) in db.proteins.iter().zip(&d) {
        let folded: Vec<u8> = s
            .iter()
            .map(|&r| if r == b'I' { b'L' } else { r })
            .collect();
        assert!(!targets.contains(&folded), "{}", String::from_utf8_lossy(s));
        let (mut a, mut b) = (t.sequence.clone(), s.clone());
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
        assert_eq!(s[4], b'K');
    }
}

fn random_db(count: usize) -> ProteinDb {
    let aa = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state = 0x1234_5678_9ABC_DEF1u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    ProteinDb {
        proteins: (0..count)
            .map(|i| {
                let len = 30 + (next() % 300) as usize;
                Protein {
                    accession: format!("P{i}"),
                    description: String::new(),
                    sequence: (0..len).map(|_| aa[(next() % 20) as usize]).collect(),
                    annotation: None,
                }
            })
            .collect(),
    }
}

#[test]
fn every_segment_keeps_composition_cuts_and_cleavage_side_residue() {
    let db = random_db(200);
    let digests = [
        enzyme(Enzyme::Trypsin),
        enzyme(Enzyme::LysC),
        enzyme(Enzyme::GluC),
        enzyme(Enzyme::Chymotrypsin),
        enzyme(Enzyme::AspN),
        enzyme(Enzyme::LysN),
        enzyme(Enzyme::NoCleavage),
        Proteases {
            primary: Enzyme::Trypsin,
            extra: vec![Enzyme::AspN],
            num_tolerable_termini: 2,
            ..Proteases::default()
        },
    ];
    for proteases in &digests {
        for strategy in [
            DecoyStrategy::Reverse,
            DecoyStrategy::Shuffle,
            DecoyStrategy::Sequon,
        ] {
            let d = decoys(&db, strategy, proteases);
            let mut changed = 0usize;
            let mut segments = 0usize;
            for (t, s) in db.proteins.iter().zip(&d) {
                let cuts = compute_cleavage_positions(
                    &t.sequence,
                    proteases.primary,
                    &proteases.extra,
                );
                assert_eq!(
                    cuts,
                    compute_cleavage_positions(s, proteases.primary, &proteases.extra),
                    "{proteases:?} {strategy:?}"
                );
                for w in cuts.windows(2) {
                    let (a, b) = (w[0] as usize, w[1] as usize);
                    let (mut x, mut y) = (t.sequence[a..b].to_vec(), s[a..b].to_vec());
                    let anchor = if proteases.primary.is_n_term() {
                        a
                    } else {
                        b - 1
                    };
                    assert_eq!(t.sequence[anchor], s[anchor], "{proteases:?} {strategy:?}");
                    segments += 1;
                    changed += usize::from(x != y);
                    x.sort_unstable();
                    y.sort_unstable();
                    assert_eq!(x, y, "{proteases:?} {strategy:?}");
                }
                if t.sequence[0] == b'M' {
                    assert_eq!(s[0], b'M');
                }
            }
            assert!(
                changed * 2 > segments,
                "{proteases:?} {strategy:?} {changed}/{segments}"
            );
        }
    }
}

#[test]
fn shuffle_is_deterministic_for_a_fixed_seed() {
    let db = random_db(20);
    let p = Proteases::default();
    let a = decoy_db(&db, "XXX", DecoyStrategy::Shuffle, DEFAULT_DECOY_SEED, &p);
    let b = decoy_db(&db, "XXX", DecoyStrategy::Shuffle, DEFAULT_DECOY_SEED, &p);
    let c = decoy_db(
        &db,
        "XXX",
        DecoyStrategy::Shuffle,
        DEFAULT_DECOY_SEED + 1,
        &p,
    );
    let seqs = |d: &ProteinDb| -> Vec<Vec<u8>> {
        d.proteins.iter().map(|p| p.sequence.clone()).collect()
    };
    assert_eq!(seqs(&a), seqs(&b));
    assert_ne!(seqs(&a), seqs(&c));
    assert_ne!(
        seqs(&a),
        seqs(&decoy_db(&db, "XXX", DecoyStrategy::Reverse, 0, &p))
    );
}

#[test]
fn decoy_db_prefixes_accessions() {
    let mut db = make_db(&[("P1", b"AB")]);
    db.proteins[0].description = "Some description".into();
    for prefix in ["XXX", "XXX_", ""] {
        let d = decoy_db(
            &db,
            prefix,
            DecoyStrategy::Reverse,
            0,
            &Proteases::default(),
        );
        assert_eq!(d.proteins[0].accession, "XXX_P1");
        assert_eq!(d.proteins[0].description, "Some description");
    }
}

#[test]
fn build_search_db_appends_decoys_after_targets() {
    let target = make_db(&[("P1", b"PEPTIDEK"), ("P2", b"AGCTR")]);
    let built = build_search_db(&target, "XXX", DecoyStrategy::Reverse, DEFAULT_DECOY_SEED);
    assert_eq!(built.len(), 4);
    assert_eq!(built.proteins[0].accession, "P1");
    assert_eq!(built.proteins[2].accession, "XXX_P1");
    assert_eq!(built.proteins[2].sequence, b"EDITPEPK");
    assert_eq!(built.proteins[3].accession, "XXX_P2");
}

#[test]
fn whole_sequence_sequon_shuffle_is_a_permutation_and_recovers_sequons() {
    let db = make_db(&[("P1", b"NISCDENLTKRSTNVTAAKNGSLLNPTVRQ")]);
    let orig = &db.proteins[0].sequence;
    let non_specific = enzyme(Enzyme::NonSpecific);
    let d = &decoys(&db, DecoyStrategy::Sequon, &non_specific)[0];
    let plain = &decoys(&db, DecoyStrategy::Shuffle, &non_specific)[0];
    let (mut a, mut b) = (orig.clone(), d.clone());
    a.sort_unstable();
    b.sort_unstable();
    assert_eq!(a, b, "composition preserved");
    let (t, s, p) = (
        count_sequon_starts(orig),
        count_sequon_starts(d),
        count_sequon_starts(plain),
    );
    assert!(s <= t.max(p), "no overshoot ({s} <= {t})");
    assert!(s >= p, "at least the plain shuffle ({s} >= {p})");
}

#[test]
fn peptide_level_sequon_shuffle_raises_the_sequon_count() {
    let db = random_db(200);
    let p = Proteases::default();
    let plain = decoys(&db, DecoyStrategy::Shuffle, &p);
    let sequon = decoys(&db, DecoyStrategy::Sequon, &p);
    let count = |seqs: &[Vec<u8>]| -> usize {
        db.proteins
            .iter()
            .zip(seqs)
            .map(|(target, s)| {
                let cuts = compute_cleavage_positions(&target.sequence, p.primary, &p.extra);
                cuts.windows(2)
                    .map(|w| count_sequon_starts(&s[w[0] as usize..w[1] as usize]))
                    .sum::<usize>()
            })
            .sum()
    };
    let targets: Vec<Vec<u8>> = db.proteins.iter().map(|p| p.sequence.clone()).collect();
    let (t, r, s) = (count(&targets), count(&plain), count(&sequon));
    assert!(s > r, "target {t}, shuffle {r}, sequon {s}");
}

#[test]
fn build_search_db_none_is_target_only() {
    let target = make_db(&[("P1", b"ABCD"), ("P2", b"EFGH")]);
    let built = build_search_db(&target, "XXX", DecoyStrategy::None, DEFAULT_DECOY_SEED);
    assert_eq!(built.len(), 2, "None appends no decoys");
    assert_eq!(built.proteins[0].accession, "P1");
    assert_eq!(built.proteins[1].accession, "P2");
}

#[test]
fn decoy_accession_needle_requires_underscore_delimiter() {
    assert_eq!(decoy_accession_needle("XXX"), "XXX_");
    assert_eq!(decoy_accession_needle("XXX_"), "XXX_");
    assert_eq!(decoy_accession_needle(""), "XXX_");
    assert!(is_decoy_accession("XXX_P12345", "XXX"));
    assert!(
        !is_decoy_accession("XXXP12345", "XXX"),
        "no '_' delimiter ⇒ target, not decoy"
    );
    assert!(!is_decoy_accession("P12345", "XXX"));
    assert!(is_decoy_accession("rev_sp|P1", "rev"));
    assert!(!is_decoy_accession("reverse_kinase", "rev"));
}

#[test]
fn affix_matches_prefix_or_suffix() {
    assert!(is_decoy_accession_affix("XXX_P12345", "XXX", None));
    assert!(!is_decoy_accession_affix("P12345", "XXX", None));
    assert!(is_decoy_accession_affix(
        "sp|P02769|ALBU_BOVIN_rev",
        "XXX",
        Some("rev")
    ));
    assert!(!is_decoy_accession_affix(
        "sp|P02769|ALBU_BOVIN",
        "XXX",
        Some("rev")
    ));
    assert!(is_decoy_accession_affix("XXX_P1", "XXX", Some("rev")));
    assert!(!is_decoy_accession_affix("P1", "XXX", Some("")));
}

#[test]
fn decoy_strategy_parses_names() {
    assert_eq!(
        DecoyStrategy::from_name("reverse"),
        Some(DecoyStrategy::Reverse)
    );
    assert_eq!(
        DecoyStrategy::from_name("Shuffle"),
        Some(DecoyStrategy::Shuffle)
    );
    assert_eq!(DecoyStrategy::from_name("NONE"), Some(DecoyStrategy::None));
    assert_eq!(DecoyStrategy::from_name("bogus"), None);
    assert_eq!(
        DecoyStrategy::from_name("sequon"),
        Some(DecoyStrategy::Sequon)
    );
    assert_eq!(DecoyStrategy::default(), DecoyStrategy::Shuffle);
}

fn digest(min_length: u32, max_length: u32) -> Proteases {
    Proteases {
        min_length,
        max_length,
        ..Proteases::default()
    }
}

fn folded(s: &[u8]) -> Vec<u8> {
    s.iter().map(|&r| if r == b'I' { b'L' } else { r }).collect()
}

fn peptides_of(seqs: &[Vec<u8>], proteases: &Proteases) -> FxHashSet<Vec<u8>> {
    let mut out = FxHashSet::default();
    for seq in seqs {
        let walk = seq.len() > 1 && seq[0] == b'M';
        collisions::for_each_span(seq, walk, proteases, |s, e| {
            out.insert(folded(&seq[s..e]));
        });
    }
    out
}

fn sequences(db: &ProteinDb) -> Vec<Vec<u8>> {
    db.proteins.iter().map(|p| p.sequence.clone()).collect()
}

#[test]
fn a_missed_cleavage_join_of_short_segments_is_redrawn() {
    let db = make_db(&[("P1", b"ASKGVR"), ("P2", b"SAKVGR")]);
    let proteases = digest(6, 30);
    let segments = TargetSegments::new(&db, &proteases);
    let unchecked: Vec<Vec<u8>> = db
        .proteins
        .iter()
        .map(|p| {
            let source = peptide_level_permutation(
                p,
                DecoyStrategy::Reverse,
                DEFAULT_DECOY_SEED,
                &proteases,
                &segments,
            );
            decoy_protein(p, "XXX", source, &proteases).sequence
        })
        .collect();
    assert_eq!(unchecked, [b"SAKVGR".to_vec(), b"ASKGVR".to_vec()]);
    let targets = peptides_of(&sequences(&db), &proteases);
    let d = decoys(&db, DecoyStrategy::Reverse, &proteases);
    for (t, s) in db.proteins.iter().zip(&d) {
        assert!(!targets.contains(&folded(s)), "{}", String::from_utf8_lossy(s));
        for (a, b) in [(0, 3), (3, 6)] {
            let (mut x, mut y) = (t.sequence[a..b].to_vec(), s[a..b].to_vec());
            assert_eq!(x[2], y[2]);
            x.sort_unstable();
            y.sort_unstable();
            assert_eq!(x, y);
        }
    }
}

#[test]
fn a_non_specific_decoy_holds_no_target_peptide() {
    let fwd = b"ACDEFGHKLMNPQRSTVWYAGLK".to_vec();
    let rev: Vec<u8> = fwd.iter().rev().copied().collect();
    let db = make_db(&[("P1", &fwd), ("P2", &rev)]);
    let proteases = Proteases {
        primary: Enzyme::NonSpecific,
        num_tolerable_termini: 0,
        ..digest(8, 12)
    };
    let targets = peptides_of(&sequences(&db), &proteases);
    for strategy in [
        DecoyStrategy::Reverse,
        DecoyStrategy::Shuffle,
        DecoyStrategy::Sequon,
    ] {
        let d = decoys(&db, strategy, &proteases);
        let hits = peptides_of(&d, &proteases)
            .into_iter()
            .filter(|p| targets.contains(p))
            .count();
        assert_eq!(hits, 0, "{strategy:?}");
        for (t, s) in db.proteins.iter().zip(&d) {
            let (mut x, mut y) = (t.sequence.clone(), s.clone());
            x.sort_unstable();
            y.sort_unstable();
            assert_eq!(x, y, "{strategy:?}");
        }
    }
}

#[test]
fn mutate_changes_joins_of_segments_shorter_than_seven() {
    let db = make_db(&[("P1", b"GASVKLLTGRDEFGHKW")]);
    let proteases = digest(6, 30);
    let targets = peptides_of(&sequences(&db), &proteases);
    let d = decoys(&db, DecoyStrategy::Mutate, &proteases);
    let hits: Vec<Vec<u8>> = peptides_of(&d, &proteases)
        .into_iter()
        .filter(|p| targets.contains(p))
        .collect();
    assert!(hits.is_empty(), "{hits:?}");
    let (t, s) = (&db.proteins[0].sequence, &d[0]);
    assert_eq!(
        compute_cleavage_positions(t, Enzyme::Trypsin, &[]),
        compute_cleavage_positions(s, Enzyme::Trypsin, &[])
    );
    let mass = |r: u8| model::AminoAcid::standard(r).unwrap().mass;
    for i in (0..t.len()).filter(|&i| t[i] != s[i]) {
        assert!((mass(s[i]) - mass(t[i])).abs() > 1.1);
        assert!(!b"KRP".contains(&s[i]) && !b"KR".contains(&t[i]));
    }
}

fn ecoli(count: usize) -> ProteinDb {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test-fixtures/ecoli.fasta");
    let file = std::fs::File::open(path).expect("open ecoli.fasta");
    let mut db = input::fasta::FastaReader::load_all(std::io::BufReader::new(file))
        .expect("read ecoli.fasta");
    db.proteins.truncate(count);
    db
}

fn search_peptides(
    db: &ProteinDb,
    proteases: &Proteases,
) -> Vec<(usize, usize, Vec<u8>)> {
    let aa = model::aa_set::AminoAcidSetBuilder::new_standard()
        .build()
        .expect("standard amino acids");
    let mut params = SearchParams::default_tryptic(aa);
    params.enzyme = proteases.primary;
    params.extra_enzymes = proteases.extra.clone();
    params.num_tolerable_termini = proteases.num_tolerable_termini;
    params.min_length = proteases.min_length;
    params.max_length = proteases.max_length;
    params.max_missed_cleavages = proteases.max_missed_cleavages;
    params.max_variable_mods_per_peptide = 0;
    let idx = crate::index::SearchIndex::from_target_db_for(
        db,
        "NONE",
        None,
        DecoyStrategy::None,
        DEFAULT_DECOY_SEED,
        proteases,
    );
    let mut out: Vec<(usize, usize, Vec<u8>)> =
        crate::candidates::enumerate_candidates(&idx, &params, "NONE")
            .map(|c| {
                let r: Vec<u8> = c.peptide.residues.iter().map(|a| a.residue).collect();
                (c.protein_index, c.start_offset_in_protein, folded(&r))
            })
            .collect();
    out.sort();
    out.dedup();
    out
}

#[test]
fn decoy_peptides_equal_to_targets_are_exactly_the_logged_remainder() {
    let db = ecoli(300);
    let tryptic = digest(7, 30);
    let semi = Proteases {
        num_tolerable_termini: 1,
        ..digest(7, 30)
    };
    let non_specific = Proteases {
        primary: Enzyme::NonSpecific,
        num_tolerable_termini: 0,
        ..digest(8, 12)
    };
    let permutations = [
        DecoyStrategy::Reverse,
        DecoyStrategy::Shuffle,
        DecoyStrategy::Sequon,
    ];
    let mut with_mutate = permutations.to_vec();
    with_mutate.push(DecoyStrategy::Mutate);
    for (label, proteases, strategies) in [
        ("tryptic", &tryptic, with_mutate),
        ("semi", &semi, permutations.to_vec()),
        ("non-specific", &non_specific, permutations.to_vec()),
    ] {
        let targets: FxHashSet<Vec<u8>> = search_peptides(&db, proteases)
            .into_iter()
            .map(|(_, _, r)| r)
            .collect();
        for strategy in strategies {
            let (proteins, repairs, _) = build_decoys(
                &db,
                "XXX",
                strategy,
                DEFAULT_DECOY_SEED,
                std::slice::from_ref(proteases),
            );
            let found: usize = repairs.iter().map(|r| r.found).sum();
            let left: usize = repairs.iter().map(|r| r.left).sum();
            let still = search_peptides(&ProteinDb { proteins }, proteases)
                .into_iter()
                .filter(|(_, _, r)| targets.contains(r))
                .count();
            assert_eq!(still, left, "[{label} {strategy:?}] found {found}");
            assert!(left <= found, "[{label} {strategy:?}] {left} of {found} left");
            if label == "tryptic" && strategy == DecoyStrategy::Shuffle {
                assert!(found > 0, "the fixture exercises no collision");
            }
        }
    }
}

#[test]
fn a_further_digest_is_checked_and_an_empty_list_changes_nothing() {
    let db = ecoli(300);
    let tryptic = digest(7, 30);
    let lysc = Proteases {
        primary: Enzyme::LysC,
        ..digest(7, 30)
    };
    let lysc_targets: FxHashSet<Vec<u8>> = search_peptides(&db, &lysc)
        .into_iter()
        .map(|(_, _, r)| r)
        .collect();
    let shared = |decoys: &ProteinDb| {
        search_peptides(decoys, &lysc)
            .into_iter()
            .filter(|(_, _, r)| lysc_targets.contains(r))
            .count()
    };
    let mut total = (0, 0);
    for strategy in [DecoyStrategy::Reverse, DecoyStrategy::Shuffle] {
        let alone = decoy_db(&db, "XXX", strategy, DEFAULT_DECOY_SEED, &tryptic);
        let same = checked_decoy_db(
            &db,
            "XXX",
            strategy,
            DEFAULT_DECOY_SEED,
            &tryptic,
            std::slice::from_ref(&tryptic),
        );
        let seqs = |d: &ProteinDb| -> Vec<Vec<u8>> {
            d.proteins.iter().map(|p| p.sequence.clone()).collect()
        };
        assert_eq!(seqs(&alone), seqs(&same.decoys), "{strategy:?}");
        let cohort =
            checked_decoy_db(&db, "XXX", strategy, DEFAULT_DECOY_SEED, &tryptic, &[lysc.clone()]);
        let (before, after) = (shared(&alone), shared(&cohort.decoys));
        assert!(after <= before, "{strategy:?}: Lys-C collisions {before} -> {after}");
        total.0 += before;
        total.1 += after;
        assert!(cohort.left <= cohort.redrawn, "{strategy:?}");
        for (t, d) in db.proteins.iter().zip(&cohort.decoys.proteins) {
            let (mut x, mut y) = (t.sequence.clone(), d.sequence.clone());
            x.sort_unstable();
            y.sort_unstable();
            assert_eq!(x, y, "{strategy:?}: a redraw keeps the composition");
        }
    }
    assert!(
        total.1 < total.0,
        "the Lys-C check redraws nothing: {} -> {}",
        total.0,
        total.1
    );
}

#[test]
fn overlap_counts_distinct_peptides_per_length_with_i_and_l_folded() {
    let targets = make_db(&[("T1", b"AAAAAAKLLLLLLKGGGGGGKAAAAAAKGGGGGGGGK")]);
    let decoys = make_db(&[("XXX_T1", b"AAAAAAKIIIIIIKCCCCCCK")]);
    let rows = peptide_overlap(&targets, &decoys, &digest(7, 7));
    assert_eq!(
        rows,
        vec![LengthOverlap {
            length: 7,
            target_peptides: 3,
            decoy_peptides: 3,
            shared_peptides: 2,
        }]
    );
}
