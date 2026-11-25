//! Tests of the filtered TSV output.

use super::outputs::write_filtered_tsv;

#[test]
fn filtered_tsv_is_written_through_a_partial_file() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("run.q0p01.tsv");
    let psm = output::PercolatorPsm {
        pooled_q_value: None,
        psm_id: "run/scan=1_1_1".into(),
        q_value: 0.001,
        pep: 0.01,
        score: Some(2.0),
        peptide: "K.PEPTIDEK.A".into(),
        proteins: "P1".into(),
    };
    write_filtered_tsv(&path, &[&psm], None).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.starts_with("PSMId\tq-value\t"));
    assert!(text.contains("run/scan=1_1_1\t0.001\t0.01\tK.PEPTIDEK.A\tP1\tNA\tNA"));
    assert!(!output::ms1::partial_path(&path).exists());
}
