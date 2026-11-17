use percolator_native::{l2_svm_mfn, read_pin, stdv, MemLimit, MfnOpts, Rows};
use std::path::PathBuf;

fn data_path(name: &str) -> PathBuf {
    PathBuf::from("tests").join("data").join(name)
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64 + 0.5) / (1u64 << 53) as f64
}

fn separable() -> (Rows, Vec<usize>, Vec<usize>) {
    let mut s = 7u64;
    let mut data = Vec::new();
    let mut pos = Vec::new();
    let mut neg = Vec::new();
    for i in 0..60 {
        let t = i % 2 == 0;
        let c = if t { 2.0 } else { -2.0 };
        data.push(c + lcg(&mut s) - 0.5);
        data.push(c + lcg(&mut s) - 0.5);
        data.push(lcg(&mut s) * 2.0 - 1.0);
        data.push(1.0);
        if t {
            pos.push(i);
        } else {
            neg.push(i);
        }
    }
    (Rows::from_dense(data, 60, 4), pos, neg)
}

fn gradient(x: &Rows, pos: &[usize], neg: &[usize], cpos: f64, cneg: f64, w: &[f64]) -> f64 {
    let mut g: Vec<f64> = w.to_vec();
    for (set, y, c) in [(neg, -1.0, cneg), (pos, 1.0, cpos)] {
        for &i in set {
            let row = x.row(i);
            let o: f64 = row.iter().zip(w).map(|(a, b)| a * b).sum();
            if 1.0 - y * o > 0.0 {
                for (gk, &xk) in g.iter_mut().zip(row) {
                    *gk += c * (o - y) * xk;
                }
            }
        }
    }
    g.iter().map(|v| v * v).sum::<f64>().sqrt()
}

#[test]
fn mfn_satisfies_kkt_on_separable_data() {
    let (x, pos, neg) = separable();
    let w = l2_svm_mfn(&x, &pos, &neg, 1.0, 1.0, &MfnOpts::default());
    assert_eq!(w.len(), 4);
    let g = gradient(&x, &pos, &neg, 1.0, 1.0, &w);
    assert!(g <= 1e-6, "gradient norm {g}");
    assert!(w[0] > 0.0 && w[1] > 0.0);
}

#[test]
fn cneg_weights_negatives() {
    let (x, pos, neg) = separable();
    let opts = MfnOpts::default();
    let w1 = l2_svm_mfn(&x, &pos, &neg, 1.0, 1.0, &opts);
    let w3 = l2_svm_mfn(&x, &pos, &neg, 1.0, 3.0, &opts);
    assert!(w3[3] < w1[3], "bias {} -> {}", w1[3], w3[3]);
    let mean_neg = |w: &[f64]| {
        neg.iter()
            .map(|&i| x.row(i).iter().zip(w).map(|(a, b)| a * b).sum::<f64>())
            .sum::<f64>()
            / neg.len() as f64
    };
    assert!(mean_neg(&w3) < mean_neg(&w1));
    let g = gradient(&x, &pos, &neg, 1.0, 3.0, &w3);
    assert!(g <= 1e-6, "gradient norm {g}");
}

#[test]
fn rows_gather_normalises_and_appends_bias() {
    let dir = tempfile::tempdir().unwrap();
    let pin = read_pin(
        &data_path("svm_ref.pin"),
        MemLimit {
            bytes: u64::MAX,
            spill_dir: dir.path().to_path_buf(),
        },
    )
    .unwrap();
    let norm = stdv(&pin);
    let x = Rows::gather(&pin.rows, &norm, &[5, 0, 200]);
    assert_eq!((x.nrows, x.ncols), (3, 6));
    assert_eq!(x.row(1)[0], (3.109375 - norm.mean[0]) * norm.inv_sd[0]);
    assert_eq!(x.row(2)[5], 1.0);
    let all: Vec<usize> = (0..pin.meta.len()).collect();
    let full = Rows::gather(&pin.rows, &norm, &all);
    assert_eq!(full.row(5), x.row(0));
    assert_eq!(full.row(200), x.row(2));
}

fn split(pin: &percolator_native::Pin, skip: Option<u64>) -> (Vec<usize>, Vec<usize>) {
    let mut pos = Vec::new();
    let mut neg = Vec::new();
    for (i, m) in pin.meta.iter().enumerate() {
        if Some(m.scan) == skip {
            continue;
        }
        if m.label == 1 {
            pos.push(i);
        } else {
            neg.push(i);
        }
    }
    (pos, neg)
}

fn printed_tol(s: &str) -> f64 {
    let dec = s.split('.').nth(1).map(|d| d.len()).unwrap_or(0);
    0.5 * 10f64.powi(-(dec as i32)) + 1e-9
}

#[test]
fn mfn_matches_reference_weights() {
    let dir = tempfile::tempdir().unwrap();
    let pin = read_pin(
        &data_path("svm_ref.pin"),
        MemLimit {
            bytes: u64::MAX,
            spill_dir: dir.path().to_path_buf(),
        },
    )
    .unwrap();
    assert_eq!(pin.meta.len(), 201);
    let norm = stdv(&pin);
    let all: Vec<usize> = (0..pin.meta.len()).collect();
    let x = Rows::gather(&pin.rows, &norm, &all);
    let opts = MfnOpts::default();
    let reference = std::fs::read_to_string(data_path("svm_ref_weights.tsv")).unwrap();
    let mut fold_fits = Vec::new();
    let mut checked = 0;
    for line in reference.lines().filter(|l| !l.is_empty()) {
        let cells: Vec<&str> = line.split('\t').collect();
        let (skip, cneg) = match cells[0] {
            "all" => (None, 1.0),
            "all_cneg3" => (None, 3.0),
            "noscan1" => (Some(1), 1.0),
            "noscan2" => (Some(2), 1.0),
            "noscan3" => (Some(3), 1.0),
            other => panic!("unknown fit {other}"),
        };
        let (pos, neg) = split(&pin, skip);
        assert_eq!(neg.len(), cells[1].parse::<usize>().unwrap());
        assert_eq!(pos.len(), cells[2].parse::<usize>().unwrap());
        let want: Vec<f64> = cells[5..].iter().map(|c| c.parse().unwrap()).collect();
        let got = l2_svm_mfn(&x, &pos, &neg, 1.0, cneg, &opts);
        assert_eq!(got.len(), want.len());
        for (k, (g, r)) in got.iter().zip(&want).enumerate() {
            assert!(
                (g - r).abs() <= 1e-6,
                "{} weight {k}: {g} vs reference {r}",
                cells[0]
            );
        }
        if skip.is_some() {
            fold_fits.push(got);
        }
        checked += 1;
    }
    assert_eq!(checked, 5);
    let perc = std::fs::read_to_string(data_path("svm_ref_percolator_v3.txt")).unwrap();
    let table: Vec<Vec<&str>> = perc
        .lines()
        .skip(2)
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(str::trim).collect())
        .collect();
    assert_eq!(table.len(), 6);
    let mut used = [false; 3];
    for split in 0..3 {
        let hit = (0..3).find(|&f| {
            !used[f]
                && table.iter().zip(&fold_fits[f]).all(|(cells, w)| {
                    (cells[split].parse::<f64>().unwrap() - w).abs() <= printed_tol(cells[split])
                })
        });
        let f = hit.unwrap_or_else(|| panic!("Percolator split {split} matches no fold fit"));
        used[f] = true;
    }
}

#[test]
#[should_panic(expected = "Rows::gather")]
fn rows_gather_panics_on_unreached_index() {
    let dir = tempfile::tempdir().unwrap();
    let pin = read_pin(
        &data_path("svm_ref.pin"),
        MemLimit {
            bytes: u64::MAX,
            spill_dir: dir.path().to_path_buf(),
        },
    )
    .unwrap();
    let norm = stdv(&pin);
    Rows::gather(&pin.rows, &norm, &[0, pin.meta.len()]);
}
