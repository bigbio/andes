use crate::args::{parse, Args};
use crate::cv::{init_direction, train, CvConfig};
use crate::folds::folds_by_spectrum;
use crate::merge::merge_folds;
use crate::normalize::{score, stdv, Norm};
use crate::pep::{estimate_pi0, pep_mix_max, pep_tdc, too_good_separation};
use crate::pin::{read_pin, MemLimit, Pin};
use crate::scores::{count_at, mix_max_q, tdc_q, tdc_winners};
use crate::write::{write_results, write_weights, Results};
use std::collections::HashMap;
use std::path::Path;

const FOLDS: usize = 3;
const SEPARATION: &str =
    "Error in the input data: too good separation between target and decoy PSMs.\nTerminating.";

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Tdc,
    MixMax,
    Plain,
}

fn concatenated_search(pin: &Pin) -> bool {
    let mut seen: HashMap<(u64, u64), bool> = HashMap::new();
    for m in &pin.meta {
        let decoy = m.label < 0;
        let key = (m.scan, (m.exp_mass + 0.0).to_bits());
        match seen.get(&key) {
            Some(&first) if first != decoy => return false,
            Some(_) => {}
            None => {
                seen.insert(key, decoy);
            }
        }
    }
    true
}

fn check_separation(scores: &[f64], is_decoy: &[bool]) -> Result<(), String> {
    if too_good_separation(scores, is_decoy) {
        return Err(SEPARATION.to_string());
    }
    Ok(())
}

pub fn run(args: &[String], limit: MemLimit) -> Result<(), String> {
    let a = parse(args)?;
    if !a.static_model && !a.tdc {
        return Err("native training needs target-decoy competition: add -Y".to_string());
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(a.threads)
        .build()
        .map_err(|e| format!("thread pool: {e}"))?;
    let pin = read_pin(&a.pin, limit)?;
    if pin.meta.is_empty() {
        return Err("ERROR: Reading tab file, could not find any PSMs.".to_string());
    }
    let is_decoy: Vec<bool> = pin.meta.iter().map(|m| m.label < 0).collect();
    if is_decoy.iter().all(|&d| d) {
        return Err("Error: no target PSMs were provided.".to_string());
    }
    if !is_decoy.iter().any(|&d| d) {
        return Err("Error: no decoy PSMs were provided.".to_string());
    }
    let mode = if a.tdc {
        Mode::Tdc
    } else if concatenated_search(&pin) {
        Mode::Plain
    } else {
        Mode::MixMax
    };
    let scores = if a.static_model {
        static_scores(&pin, &a, &is_decoy)?
    } else {
        trained_scores(&pin, &a, &is_decoy, &pool)?
    };
    finish(&pin, &a, &is_decoy, &scores, mode)
}

fn trained_scores(
    pin: &Pin,
    a: &Args,
    is_decoy: &[bool],
    pool: &rayon::ThreadPool,
) -> Result<Vec<f64>, String> {
    let norm = stdv(pin);
    let cfg = CvConfig {
        seed: a.seed,
        overrule: a.overrule,
        verbose: a.verbose,
        ..CvConfig::default()
    };
    let mut cv = train(pin, &norm, &cfg, pool)?;
    merge_folds(&mut cv, is_decoy, cfg.selection_fdr)?;
    check_separation(&cv.test_scores, is_decoy)?;
    if let Some(path) = &a.weights {
        write_weights(path, &pin.feature_names, &cv.w_norm, &cv.w_raw)?;
    }
    Ok(cv.test_scores)
}

fn number(token: &str) -> Option<f64> {
    let plain = token
        .bytes()
        .all(|b| b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'.' | b'e' | b'E'));
    if plain {
        token.parse().ok()
    } else {
        None
    }
}

fn read_init_weights(path: &Path, nf: usize) -> Option<Vec<f64>> {
    let text = std::fs::read(path).ok()?;
    let text = String::from_utf8_lossy(&text);
    let mut rest = "";
    let mut remaining: &str = &text;
    while !remaining.is_empty() {
        let (line, tail) = match remaining.find('\n') {
            Some(k) => (&remaining[..k], &remaining[k + 1..]),
            None => (remaining, ""),
        };
        remaining = tail;
        let first = line.as_bytes().first().copied().unwrap_or(0);
        if first == b'-' || first.is_ascii_digit() {
            rest = tail;
            break;
        }
    }
    let mut w = vec![0.0f64; nf + 1];
    for (slot, token) in w.iter_mut().zip(rest.split_ascii_whitespace()) {
        match number(token) {
            Some(v) => *slot = v,
            None => break,
        }
    }
    Some(w)
}

fn identity(nf: usize) -> Norm {
    Norm {
        mean: vec![0.0; nf],
        inv_sd: vec![1.0; nf],
    }
}

fn static_scores(pin: &Pin, a: &Args, is_decoy: &[bool]) -> Result<Vec<f64>, String> {
    let n = pin.meta.len();
    let nf = pin.rows.nfeat();
    let norm = identity(nf);
    let fold = folds_by_spectrum(&pin.meta, FOLDS, a.seed);
    let test: Vec<Vec<usize>> = (0..FOLDS)
        .map(|f| (0..n).filter(|&i| fold[i] as usize == f).collect())
        .collect();
    let train: Vec<Vec<usize>> = (0..FOLDS)
        .map(|f| (0..n).filter(|&i| fold[i] as usize != f).collect())
        .collect();
    let default_direction = |overrule: bool| -> Result<Vec<Vec<f64>>, String> {
        train
            .iter()
            .map(|rows| init_direction(pin, &norm, rows, 0.01, overrule))
            .collect()
    };
    let path = a.init_weights.as_deref().unwrap_or(Path::new(""));
    let mut w: Vec<Vec<f64>> = match read_init_weights(path, nf) {
        Some(ww) => vec![ww; FOLDS],
        None => {
            eprintln!(
                "WARNING: Could not find weights input file {}. Using default weights instead.",
                path.display()
            );
            default_direction(false)?
        }
    };
    let positives = |w: &[Vec<f64>]| -> usize {
        test.iter()
            .zip(w)
            .map(|(rows, wf)| {
                let s = score(&pin.rows, rows, wf, None);
                let d: Vec<bool> = rows.iter().map(|&i| is_decoy[i]).collect();
                count_at(&tdc_q(&s, &d, false), &d, 0.01)
            })
            .sum()
    };
    let valid = positives(&w) > 0;
    if !valid && !a.overrule {
        let from_pin = pin
            .default_direction
            .as_deref()
            .is_some_and(|d| d.iter().any(|&v| v != 0.0));
        for (wf, dd) in w.iter_mut().zip(default_direction(false)?) {
            if from_pin {
                *wf = dd;
            } else {
                wf[..nf].copy_from_slice(&dd[..nf]);
            }
        }
    }
    let mut scores = vec![0.0f64; n];
    for (f, rows) in test.iter().enumerate() {
        let wf = if valid { &w[f] } else { &w[0] };
        for (&i, s) in rows.iter().zip(score(&pin.rows, rows, wf, None)) {
            scores[i] = s;
        }
    }
    for rows in &test {
        let s: Vec<f64> = rows.iter().map(|&i| scores[i]).collect();
        let d: Vec<bool> = rows.iter().map(|&i| is_decoy[i]).collect();
        check_separation(&s, &d)?;
    }
    check_separation(&scores, is_decoy)?;
    if let Some(path) = &a.weights {
        write_weights(path, &pin.feature_names, &w, &w)?;
    }
    Ok(scores)
}

fn finish(
    pin: &Pin,
    a: &Args,
    is_decoy: &[bool],
    scores: &[f64],
    mode: Mode,
) -> Result<(), String> {
    let rows: Vec<usize> = if mode == Mode::Tdc {
        tdc_winners(&pin.meta, scores)
    } else {
        (0..scores.len()).collect()
    };
    let s: Vec<f64> = rows.iter().map(|&i| scores[i]).collect();
    let d: Vec<bool> = rows.iter().map(|&i| is_decoy[i]).collect();
    if mode == Mode::Tdc {
        check_separation(&s, &d)?;
    }
    let (q, pep) = if mode == Mode::MixMax {
        let pi0 = estimate_pi0(&s, &d);
        (mix_max_q(&s, &d, pi0), pep_mix_max(&s, &d, pi0))
    } else {
        (tdc_q(&s, &d, false), pep_tdc(&s, &d))
    };
    let res = Results {
        rows: &rows,
        scores: &s,
        q: &q,
        pep: &pep,
    };
    write_results(
        &a.target,
        &a.decoy,
        &pin.meta,
        &pin.texts,
        &pin.file_names,
        &res,
    )
}
