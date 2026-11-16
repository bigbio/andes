use percolator_native::{read_pin, MemLimit};
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
    }
}

fn blocks(s: &percolator_native::Store) -> Vec<(usize, Vec<f32>)> {
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
    assert_eq!(pin.meta[0].psm_id, "a/1");
    assert_eq!(pin.meta[0].scan, 1);
    assert_eq!(pin.meta[0].exp_mass, 1000.5);
    assert_eq!(pin.meta[0].peptide, "K.AAA.R");
    assert_eq!(pin.meta[0].proteins, "P1\tP2");
    assert_eq!(pin.meta[3].proteins, "P5\tP6\tP7");
    assert_eq!(pin.meta[4].proteins, "decoy_P8");
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
        },
    )
    .unwrap();
    assert_eq!(blocks(&a.rows), blocks(&b.rows));
    assert_eq!(blocks(&b.rows), blocks(&b.rows));
    assert_eq!(std::fs::read_dir(spill.path()).unwrap().count(), 1);
    drop(b);
    assert_eq!(std::fs::read_dir(spill.path()).unwrap().count(), 0);
}
