//! Parity of `search::ipf` with PyProphet's own IPF code.
//!
//! `fixtures/ipf_parity.json` comes from `fixtures/ipf_parity_fixture.py`, which runs
//! PyProphet's `apply_bm` (log and linear space) and `compute_model_fdr` from
//! pyprophet/ipf.py (PyProphet commit 4344e2b) on 300 random PSMs and 208 PEPs.

use search::ipf::{model_fdr, posterior, Evidence};

fn fixture() -> serde_json::Value {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/ipf_parity.json"
    );
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

#[test]
fn posterior_matches_pyprophet_apply_bm() {
    let f = fixture();
    let cases = f["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 300);
    let mut worst = 0.0f64;
    for (k, c) in cases.iter().enumerate() {
        let psm_pep = c["psm_pep"].as_f64().unwrap();
        let n = c["n_forms"].as_u64().unwrap() as usize;
        let ev: Vec<Evidence> = c["evidence"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| Evidence {
                pep: e[0].as_f64().unwrap(),
                forms: e[1]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_u64().unwrap() as u32)
                    .collect(),
            })
            .collect();
        let p = posterior(psm_pep, n, &ev, 0.0);
        let ours: Vec<f64> = p.forms.iter().copied().chain([p.h0]).collect();
        for space in ["posterior_log", "posterior_linear"] {
            let want: Vec<f64> = c[space]
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_f64().unwrap())
                .collect();
            assert_eq!(want.len(), ours.len(), "case {k}");
            for (a, b) in ours.iter().zip(&want) {
                let d = (a - b).abs();
                worst = worst.max(d);
                assert!(
                    d < 1e-9,
                    "case {k} {space}: ours {ours:?} pyprophet {want:?}"
                );
            }
        }
    }
    eprintln!("max |andes - pyprophet| posterior = {worst:e}");
}

#[test]
fn model_fdr_matches_pyprophet_compute_model_fdr() {
    let f = fixture();
    let peps: Vec<f64> = f["fdr"]["peps"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    let want: Vec<f64> = f["fdr"]["qvalues"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    let q = model_fdr(&peps);
    for (i, (a, b)) in q.iter().zip(&want).enumerate() {
        assert!(
            (a - b).abs() < 1e-12,
            "pep[{i}]={}: ours {a} pyprophet {b}",
            peps[i]
        );
    }
}
