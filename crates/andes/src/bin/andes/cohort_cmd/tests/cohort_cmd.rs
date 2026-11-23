use super::*;

#[test]
fn csm_rows_are_read_by_header_and_unscored_rows_skipped() {
    let t = tempfile::tempdir().unwrap();
    let p = t.path().join("A.r1.xl.tsv");
    std::fs::write(
        &p,
        "spectrum\tscan\tspec_id\tscore_source\tclass\tlink\tq_value\tscore\tresidue_pair\tprotein_pair\n\
         A.r1/s1\t1\txl0\tpair\tTT\tinter\t0.01\t12.5\tP1:5 -- P2:9\tP1 -- P2\n\
         A.r1/s2\t2\txl1\tpair\tTD\tintra\tNA\tNA\tP1:5 -- P1:30\tP1 -- P1\n\
         A.r1/s3\t1\txl2\tpair\tDD\thomotypic\t0.5\t3.0\tP3:7 -- P3:7\tP3 -- P3\n",
    )
    .unwrap();
    let mut rows = Vec::new();
    assert_eq!(read_csms(&p, 4, &mut rows).unwrap(), (2, 1));
    assert_eq!(rows[0].run, "A.r1");
    assert_eq!((rows[0].table, rows[0].spec_id.as_str()), (4, "xl0"));
    assert_eq!(
        rows[1].spec_id, "xl2",
        "a second file's spectrum with the same scan number keeps its own key"
    );
    assert_eq!(rows[0].class, xl::decoy::PairClass::TT);
    assert_eq!(rows[1].link, xl::rollup::LinkClass::Homotypic);
    let old = t.path().join("B.xl.tsv");
    std::fs::write(&old, "spectrum\tclass\tlink\tscore\n").unwrap();
    assert!(read_csms(&old, 0, &mut rows)
        .unwrap_err()
        .contains("residue_pair"));
}

#[test]
fn a_malformed_score_or_missing_pair_is_an_error() {
    let t = tempfile::tempdir().unwrap();
    let head = "spectrum\tscan\tspec_id\tscore_source\tclass\tlink\tscore\tresidue_pair\tprotein_pair\n";
    for (name, row, what) in [
        (
            "a.xl.tsv",
            "s1\t1\txl0\tpair\tTT\tinter\t12x\tP1:5 -- P2:9\tP1 -- P2",
            "row 2: bad score",
        ),
        (
            "b.xl.tsv",
            "s1\t1\txl0\tpair\tTT\tinter\tnan\tP1:5 -- P2:9\tP1 -- P2",
            "bad score",
        ),
        (
            "c.xl.tsv",
            "s1\t1\txl0\tpair\tTT\tinter\t12.5\t\tP1 -- P2",
            "bad residue_pair",
        ),
        (
            "e.xl.tsv",
            "s1\t1\txl0\tpercolator\tTT\tinter\t2.1\tP1:5 -- P2:9\tP1 -- P2",
            "rescored on its own",
        ),
        (
            "d.xl.tsv",
            "s1\t1\txl0\tpair\tTT\tinter\t12.5\tP1:5 -- P2:9",
            "bad protein_pair",
        ),
    ] {
        let p = t.path().join(name);
        std::fs::write(&p, format!("{head}{row}\n")).unwrap();
        let e = read_csms(&p, 0, &mut Vec::new()).unwrap_err();
        assert!(e.contains(what), "{name}: {e}");
    }
}

#[test]
fn a_table_pin_is_its_sibling_and_without_one_intra_links_keep_the_pair_score() {
    let t = tempfile::tempdir().unwrap();
    let table = t.path().join("A.r1.xl.tsv");
    assert_eq!(pin_of(&table), t.path().join("A.r1.xl.pin"));
    let mut rows = vec![CsmRow {
        run: "A.r1".into(),
        table: 0,
        spec_id: "xl7".into(),
        spectrum: "A.r1/s7".into(),
        class: xl::decoy::PairClass::TT,
        link: xl::rollup::LinkClass::Intra,
        score: 12.5,
        score_source: "pair",
        residue_pair: "P1:5 -- P2:9".into(),
        protein_pair: "P1 -- P2".into(),
        peptides: [Vec::new(), Vec::new()],
        proteins: [Vec::new(), Vec::new()],
    }];
    rescore(
        &[table],
        &mut rows,
        &[xl::rollup::LinkClass::Intra],
        &[],
        "intra-protein links",
    )
    .unwrap();
    assert_eq!((rows[0].score, rows[0].score_source), (12.5, "pair"));
}

#[test]
fn intra_links_keep_the_pair_score_when_percolator_cannot_train() {
    let t = tempfile::tempdir().unwrap();
    let table = t.path().join("A.r1.xl.tsv");
    std::fs::write(
        pin_of(&table),
        "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tpair_score\tPeptide\tProteins\n\
         s7\t1\t7\t1000.0\t1000.0\t12.5\tK.PEPK.R\tP1\n\
         s8\t-1\t8\t1000.0\t1000.0\t3.0\tK.KPEP.R\tdecoy_P1\n",
    )
    .unwrap();
    let mut rows = vec![CsmRow {
        run: "A.r1".into(),
        table: 0,
        spec_id: "s7".into(),
        spectrum: "A.r1/s7".into(),
        class: xl::decoy::PairClass::TT,
        link: xl::rollup::LinkClass::Intra,
        score: 12.5,
        score_source: "pair",
        residue_pair: "P1:5 -- P1:9".into(),
        protein_pair: "P1 -- P1".into(),
        peptides: [Vec::new(), Vec::new()],
        proteins: [Vec::new(), Vec::new()],
    }];
    rescore(
        &[table],
        &mut rows,
        &[xl::rollup::LinkClass::Intra],
        &[],
        "intra-protein links",
    )
    .unwrap();
    assert_eq!((rows[0].score, rows[0].score_source), (12.5, "pair"));
}

#[test]
fn an_inter_link_with_a_near_copy_in_an_unlisted_protein_is_ambiguous() {
    let t = tempfile::tempdir().unwrap();
    let db = t.path().join("db.fasta");
    std::fs::write(
        &db,
        ">P1\nMKTLDEAGHSKR\n>P2\nAKTIDEAGNSKR\n>P3\nMKAAAAWWWKR\n>P4\nMKGGGGYYYKR\n",
    )
    .unwrap();
    let row = |a: &str, pa: &str, b: &str, pb: &str| CsmRow {
        run: "A.r1".into(),
        table: 0,
        spec_id: "xl1".into(),
        spectrum: "s".into(),
        class: xl::decoy::PairClass::TT,
        link: xl::rollup::LinkClass::Inter,
        score: 1.0,
        score_source: "pair",
        residue_pair: String::new(),
        protein_pair: String::new(),
        peptides: [
            crate::xl_ambiguity::folded(a),
            crate::xl_ambiguity::folded(b),
        ],
        proteins: [vec![pa.to_string()], vec![pb.to_string()]],
    };
    let mut rows = vec![
        row("K.TLDEAGHSK.R", "P1", "K.AAAAWWWK.R", "P3"),
        row("K.GGGGYYYK.R", "P4", "K.AAAAWWWK.R", "P3"),
    ];
    assert_eq!(flag_ambiguous(&mut rows, &db).unwrap(), 1);
    assert_eq!(
        rows[0].link,
        xl::rollup::LinkClass::InterAmbiguous,
        "P2 holds a one-substitution copy"
    );
    assert_eq!(rows[1].link, xl::rollup::LinkClass::Inter);
}

#[test]
fn staging_replaces_a_link_to_another_file() {
    let t = tempfile::tempdir().unwrap();
    let (a, b, at) = (
        t.path().join("a.raw"),
        t.path().join("b.raw"),
        t.path().join("r.raw"),
    );
    std::fs::write(&a, "a").unwrap();
    std::fs::write(&b, "b").unwrap();
    stage(&a, &at).unwrap();
    stage(&a, &at).unwrap();
    stage(&b, &at).unwrap();
    assert_eq!(std::fs::read_to_string(&at).unwrap(), "b");
}

#[test]
fn the_fingerprint_follows_the_search_settings_but_not_the_cohort_section() {
    let t = tempfile::tempdir().unwrap();
    let file = t.path().join("r1.raw");
    std::fs::write(&file, "spectra").unwrap();
    let config = t.path().join("c.yaml");
    let run = |args: &[&str]| cohort::PlannedRun {
        run: "A.r1".into(),
        project: "A".into(),
        sdrf: t.path().join("a.sdrf.tsv"),
        file: file.clone(),
        dir: t.path().join("runs/A.r1"),
        stage_as: "A.r1.raw".into(),
        kind: cohort::Kind::Xl,
        group: "xl".into(),
        args: args.iter().map(|s| s.to_string()).collect(),
    };
    std::fs::write(&config, "search: {threads: 4}\ncohort: {name: one}\n").unwrap();
    let base = fingerprint(
        &run(&["--xl", "DSSO"]),
        &["--xl".into(), "DSSO".into()],
        &config,
    )
    .unwrap();
    std::fs::write(&config, "search: {threads: 4}\ncohort: {name: two}\n").unwrap();
    assert_eq!(
        fingerprint(
            &run(&["--xl", "DSSO"]),
            &["--xl".into(), "DSSO".into()],
            &config
        )
        .unwrap(),
        base
    );
    assert_ne!(
        fingerprint(
            &run(&["--xl", "DSBU"]),
            &["--xl".into(), "DSBU".into()],
            &config
        )
        .unwrap(),
        base
    );
    std::fs::write(&config, "search: {threads: 8}\ncohort: {name: two}\n").unwrap();
    assert_ne!(
        fingerprint(
            &run(&["--xl", "DSSO"]),
            &["--xl".into(), "DSSO".into()],
            &config
        )
        .unwrap(),
        base
    );
    std::fs::write(&config, "search: {threads: 4}\ncohort: {name: two}\n").unwrap();
    std::fs::write(&file, "other spectra").unwrap();
    assert_ne!(
        fingerprint(
            &run(&["--xl", "DSSO"]),
            &["--xl".into(), "DSSO".into()],
            &config
        )
        .unwrap(),
        base
    );

    std::fs::write(t.path().join("db.fasta"), ">P1\nPEPTIDEK\n").unwrap();
    std::fs::write(&config, "io: {database: db.fasta}\ncohort: {name: two}\n").unwrap();
    let with_db = fingerprint(
        &run(&["--xl", "DSSO"]),
        &["--xl".into(), "DSSO".into()],
        &config,
    )
    .unwrap();
    assert!(with_db.contains("db.fasta"), "{with_db}");
    std::fs::write(t.path().join("db.fasta"), ">P1\nPEPTIDEKR\n").unwrap();
    assert_ne!(
        fingerprint(
            &run(&["--xl", "DSSO"]),
            &["--xl".into(), "DSSO".into()],
            &config
        )
        .unwrap(),
        with_db,
        "a database replaced at the same path"
    );
}

#[test]
fn projects_come_from_the_plan_when_one_is_given() {
    let t = tempfile::tempdir().unwrap();
    let p = t.path().join("plan.json");
    std::fs::write(
        &p,
        r#"{"runs": [{"run": "A.v1.r1", "project": "A.v1"}, {"run": "r2", "project": "B"}]}"#,
    )
    .unwrap();
    let projects = projects_of_plan(&p).unwrap();
    assert_eq!(project_of("A.v1.r1", &projects), "A.v1");
    assert_eq!(project_of("r2", &projects), "B");
    assert_eq!(project_of("C.r3", &projects), "C");
    std::fs::write(&p, r#"{"runs": [{"run": "x"}]}"#).unwrap();
    assert!(projects_of_plan(&p).is_err());
}

#[test]
fn each_project_share_is_counted_at_the_pooled_cuts() {
    use xl::decoy::PairClass::{DD, TD, TT};
    use xl::rollup::LinkClass::{Homotypic, Intra};
    let row = |run: &str, class, link| CsmRow {
        run: run.to_string(),
        table: 0,
        spec_id: String::new(),
        spectrum: String::new(),
        class,
        link,
        score: 0.0,
        score_source: "pair",
        residue_pair: String::new(),
        protein_pair: String::new(),
        peptides: [Vec::new(), Vec::new()],
        proteins: [Vec::new(), Vec::new()],
    };
    let csms = vec![
        row("A.r1", TT, Intra),
        row("A.r1", TT, Intra),
        row("A.r2", TD, Intra),
        row("A.r2", TT, Intra),
        row("B.r1", TT, Intra),
        row("B.r1", DD, Intra),
        row("B.r1", TD, Homotypic),
        row("B.r1", TT, Homotypic),
    ];
    let q = vec![
        Some(0.001),
        Some(0.02),
        Some(0.03),
        None,
        Some(0.004),
        Some(0.04),
        Some(0.01),
        Some(0.01),
    ];
    let rows = project_error(&csms, &q, &Projects::new());
    let at = |p: &str, link: &str, cut: f64| {
        rows.iter()
            .find(|r| r.project == p && r.link == link && r.q_cut == cut)
            .map(|r| (r.tt, r.td, r.dd, r.decoy_fdr))
            .unwrap()
    };
    assert_eq!(at("A", "intra", 0.01), (1, 0, 0, 0.0));
    assert_eq!(at("A", "intra", 0.05), (2, 1, 0, 0.5));
    assert_eq!(at("B", "intra", 0.05), (1, 0, 1, 0.0));
    assert_eq!(at("B", "homotypic", 0.01), (1, 1, 0, 1.0));
    assert_eq!(rows.len(), 6);
    assert_eq!(project_of("PXD1.run.with.dots", &Projects::new()), "PXD1");
}
