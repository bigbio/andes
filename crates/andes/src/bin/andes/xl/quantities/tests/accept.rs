//! Tests of reading the q-values of `andes cohort xl-fdr`.

use super::*;

#[test]
fn only_target_rows_give_q_values() {
    let dir = tempfile::tempdir().expect("tempdir");
    let base = dir.path().join("cohort");
    std::fs::write(
        with_suffix(&base, ".xl.residue_pairs.tsv"),
        "residue_pair\tclass\tscore\tq_value\nP:2 -- P:9\tTT\t9.1\t0.004\nP:3 -- rev_P:9\tTD\t8.0\t0.004\n\
         P:5 -- P:9\tTT\t1.0\tNA\n",
    )
    .expect("written");
    std::fs::write(
        with_suffix(&base, ".xl.tsv"),
        "spectrum\tscan\tclass\tq_value\nrun1.17.17.3\t17\tTT\t0.01\nrun1.18.18.3\t18\tDD\t0.5\n",
    )
    .expect("written");
    let cohort = CohortFdr::read(&base).expect("read");
    assert_eq!(cohort.residue_pairs.len(), 1);
    assert_eq!(cohort.residue_pairs.get("P:2 -- P:9"), Some(&0.004));
    assert_eq!(cohort.matches.len(), 1);
    assert_eq!(cohort.matches.get("run1.17.17.3"), Some(&0.01));
}

#[test]
fn a_table_that_is_not_of_the_cohort_fdr_is_an_error_that_names_it() {
    let dir = tempfile::tempdir().expect("tempdir");
    let base = dir.path().join("cohort");
    let missing = CohortFdr::read(&base).err().expect("no file");
    assert!(
        missing.contains("--xl-cohort-fdr") && missing.contains("cohort.xl.residue_pairs.tsv"),
        "{missing}"
    );
    std::fs::write(with_suffix(&base, ".xl.residue_pairs.tsv"), "a\tb\n1\t2\n").expect("written");
    let wrong = CohortFdr::read(&base).err().expect("wrong table");
    assert!(wrong.contains("`residue_pair` column"), "{wrong}");
}
