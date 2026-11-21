use rescore::{read_pin, MemLimit};
use std::path::{Path, PathBuf};

const PIN: &str =
    "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tFileName\tf1\tf2\tf3\tPeptide\tProteins\n\
DefaultDirection\t-\t-\t-\t-\t-\t1.5\t0\t-2\t-\t-\n\
a/1\t1\t1\t1000.5\t1000.4\ta.mzML\t0.5\t1\t2\tK.AAA.R\tP1\tP2\n\
a/2\t-1\t2\t1001.5\t1001.4\ta.mzML\t0.25\t3\t4\tK.BBB.R\tdecoy_P3\n\
a/3\t1\t3\t1002.5\t1002.4\ta.mzML\t0.125\t5\t6\tK.CCC.R\tP4\n\
b/1\t1\t1\t1003.5\t1003.4\tb.mzML\t-1\t7\t8\tK.DDD.R\tP5\tP6\tP7\n\
b/2\t-1\t2\t1004.5\t1004.4\tb.mzML\t-2\t9\t10\tK.EEE.R\tdecoy_P8\n\
b/3\t1\t3\t1005.5\t1005.4\tb.mzML\t-3\t11\t12\tK.FFF.R\tP9\n";

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn ram(dir: &Path) -> MemLimit {
    MemLimit {
        bytes: u64::MAX,
        spill_dir: dir.to_path_buf(),
        sizing: rescore::Sizing::Fixed,
    }
}

fn blocks(s: &rescore::Store) -> Vec<(usize, Vec<f32>)> {
    let mut v = Vec::new();
    s.for_each_block(|start, b| v.push((start, b.to_vec())));
    v
}

#[test]
fn reads_andes_pin_header_and_rows() {
    let d = tempfile::tempdir().unwrap();
    let p = write(d.path(), "x.pin", PIN);
    let pin = read_pin(&p, ram(d.path())).unwrap();
    assert_eq!(pin.feature_names, vec!["f1", "f2", "f3"]);
    assert_eq!(pin.default_direction, Some(vec![1.5, 0.0, -2.0]));
    assert_eq!(pin.rows.len(), 6);
    assert_eq!(pin.rows.nfeat(), 3);
    assert_eq!(pin.meta.len(), 6);
    let labels: Vec<i8> = pin.meta.iter().map(|m| m.label).collect();
    assert_eq!(labels, vec![1, -1, 1, 1, -1, 1]);
    let files: Vec<u32> = pin.meta.iter().map(|m| m.file).collect();
    assert_eq!(files, vec![0, 0, 0, 1, 1, 1]);
    assert_eq!(pin.text(0).unwrap().psm_id, "a/1");
    assert_eq!(pin.meta[0].scan, 1);
    assert_eq!(pin.meta[0].exp_mass, 1000.5);
    assert_eq!(pin.text(0).unwrap().peptide, "K.AAA.R");
    assert_eq!(pin.text(0).unwrap().proteins, "P1\tP2");
    assert_eq!(pin.text(3).unwrap().proteins, "P5\tP6\tP7");
    assert_eq!(pin.text(4).unwrap().proteins, "decoy_P8");
    let b = blocks(&pin.rows);
    assert_eq!(b.len(), 1);
    assert_eq!(b[0].0, 0);
    assert_eq!(&b[0].1[..6], &[0.5, 1.0, 2.0, 0.25, 3.0, 4.0]);
    assert_eq!(&b[0].1[15..], &[-3.0, 11.0, 12.0]);
}

#[test]
fn filename_is_metadata() {
    let d = tempfile::tempdir().unwrap();
    let p = write(d.path(), "x.pin", PIN);
    let pin = read_pin(&p, ram(d.path())).unwrap();
    for n in [
        "FileName", "ExpMass", "CalcMass", "ScanNr", "SpecId", "Label",
    ] {
        assert!(!pin.feature_names.iter().any(|f| f == n));
    }
    let lower = PIN
        .replace("FileName", "filename")
        .replace("ExpMass", "expmass");
    let p2 = write(d.path(), "y.pin", &lower);
    let pin2 = read_pin(&p2, ram(d.path())).unwrap();
    assert_eq!(pin2.feature_names, vec!["f1", "f2", "f3"]);
    assert_eq!(pin2.meta[3].file, 1);
    let nofile = "SpecId\tLabel\tScanNr\tf1\tPeptide\tProteins\ns1\t1\t1\t0.5\tK.A.R\tP1\ns2\t-1\t2\t0.7\tK.B.R\tdecoy_P2\n";
    let p3 = write(d.path(), "z.pin", nofile);
    let pin3 = read_pin(&p3, ram(d.path())).unwrap();
    assert_eq!(pin3.feature_names, vec!["f1"]);
    assert_eq!(pin3.default_direction, None);
    assert!(pin3.meta.iter().all(|m| m.file == 0));
}

#[test]
fn rejects_a_header_without_features() {
    let d = tempfile::tempdir().unwrap();
    let text = "SpecId\tLabel\tScanNr\tExpMass\tPeptide\tProteins\ns1\t1\t1\t1000\tK.A.R\tP1\n";
    let p = write(d.path(), "x.pin", text);
    let err = read_pin(&p, ram(d.path())).err().unwrap();
    assert!(err.contains("no feature columns"), "{err}");
}

#[test]
fn a_row_without_its_label_cell_is_an_error_not_a_panic() {
    let d = tempfile::tempdir().unwrap();
    let text = "SpecId\tf1\tPeptide\tLabel\ns1\t0.5\tK.A.R\t1\ns2\t0.7\tK.B.R\n";
    let p = write(d.path(), "x.pin", text);
    let err = read_pin(&p, ram(d.path())).err().unwrap();
    assert!(
        err.contains("line 3") && err.contains("Label column is 4"),
        "{err}"
    );
}

#[test]
fn rejects_bad_label() {
    let d = tempfile::tempdir().unwrap();
    let bad = PIN.replace("a/2\t-1", "a/2\t2");
    let p = write(d.path(), "x.pin", &bad);
    assert!(read_pin(&p, ram(d.path())).is_err());
}

#[test]
fn spilled_store_equals_ram_store() {
    let d = tempfile::tempdir().unwrap();
    let p = write(d.path(), "x.pin", PIN);
    let a = read_pin(&p, ram(d.path())).unwrap();
    let spill = tempfile::tempdir().unwrap();
    let b = read_pin(
        &p,
        MemLimit {
            bytes: 1,
            spill_dir: spill.path().to_path_buf(),
            sizing: rescore::Sizing::Fixed,
        },
    )
    .unwrap();
    assert_eq!(blocks(&a.rows), blocks(&b.rows));
    assert_eq!(blocks(&b.rows), blocks(&b.rows));
    for i in 0..a.meta.len() {
        let (x, y) = (a.text(i).unwrap(), b.text(i).unwrap());
        assert_eq!(
            (x.psm_id, x.peptide, x.proteins),
            (y.psm_id, y.peptide, y.proteins)
        );
    }
    assert_eq!(std::fs::read_dir(spill.path()).unwrap().count(), 2);
    drop(b);
    assert_eq!(std::fs::read_dir(spill.path()).unwrap().count(), 0);
}

#[test]
fn rejects_non_finite_feature_values() {
    let d = tempfile::tempdir().unwrap();
    for bad in ["nan", "inf", "-inf", "1e300"] {
        let text = PIN.replace("\t0.25\t3\t4\t", &format!("\t0.25\t{bad}\t4\t"));
        let p = write(d.path(), "bad.pin", &text);
        let err = read_pin(&p, ram(d.path())).err().unwrap();
        assert!(err.contains("line 4"), "{err}");
        assert!(err.contains("column 8"), "{err}");
        assert!(err.contains("f2"), "{err}");
    }
}

#[test]
fn an_automatic_budget_is_what_the_pin_needs_within_the_allowance() {
    let d = tempfile::tempdir().unwrap();
    let p = write(d.path(), "x.pin", PIN);
    let auto = |bytes: u64| MemLimit {
        bytes,
        spill_dir: d.path().to_path_buf(),
        sizing: rescore::Sizing::FitInput {
            cap: u64::MAX,
            why: "test allowance".to_string(),
        },
    };
    let roomy = read_pin(&p, auto(1 << 40)).unwrap();
    let need = roomy.budget.total();
    assert!(need < 1 << 40, "the budget stayed at the allowance");
    assert!(!roomy.rows.is_spilled());
    let fixed = read_pin(&p, ram(d.path())).unwrap();
    assert_eq!(blocks(&roomy.rows), blocks(&fixed.rows));
    drop(roomy);
    let tight = read_pin(&p, auto(need / 4)).unwrap();
    assert!(tight.budget.total() <= need / 4);
    assert_eq!(blocks(&tight.rows), blocks(&fixed.rows));
}

#[test]
fn an_automatic_budget_never_exceeds_its_cap() {
    let d = tempfile::tempdir().unwrap();
    let p = write(d.path(), "x.pin", PIN);
    let capped = |cap: u64| MemLimit {
        bytes: 1 << 40,
        spill_dir: d.path().to_path_buf(),
        sizing: rescore::Sizing::FitInput {
            cap,
            why: "test allowance".to_string(),
        },
    };
    let need = read_pin(&p, capped(u64::MAX)).unwrap().budget.total();
    let pin = read_pin(&p, capped(need / 3)).unwrap();
    assert!(pin.budget.total() <= need / 3);
    let fixed = read_pin(&p, ram(d.path())).unwrap();
    assert_eq!(blocks(&pin.rows), blocks(&fixed.rows));
}
