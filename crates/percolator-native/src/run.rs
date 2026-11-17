use crate::args::{parse, Args};
use crate::cv::{make_room, train, CvConfig};
use crate::direction::init_directions;
use crate::folds::folds_by_spectrum;
use crate::merge::merge_folds;
use crate::normalize::{score, stdv, Norm};
use crate::pep::{estimate_pi0, pep_from_q, too_good_separation};
use crate::pin::{read_pin, MemLimit, Pin};
use crate::scores::{mix_max_q, psm_order, tdc_q, tdc_winners};
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
    let mut pin = pool.install(|| read_pin(&a.pin, limit))?;
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
        pool.install(|| static_scores(&pin, &a, &is_decoy))?
    } else {
        make_room(&mut pin)?;
        pool.install(|| trained_scores(&pin, &a, &is_decoy, &pool))?
    };
    pool.install(|| finish(&pin, &a, &is_decoy, &scores, mode))
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
        test_fdr: a.test_fdr,
        selection_fdr: a.selection_fdr,
        initial_selection_fdr: a.initial_selection_fdr,
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

fn read_init_weights(path: &Path, names: &[String]) -> Result<Option<Vec<f64>>, String> {
    let Ok(bytes) = std::fs::read(path) else {
        return Ok(None);
    };
    let text = String::from_utf8_lossy(&bytes);
    let file = path.display();
    let n = names.len() + 1;
    let mut lines = text.lines();
    let mut line_nr = 0usize;
    let mut names_line = "";
    let mut found = false;
    for line in lines.by_ref() {
        line_nr += 1;
        let first = line.as_bytes().first().copied().unwrap_or(0);
        if first == b'-' || first.is_ascii_digit() {
            found = true;
            break;
        }
        if !line.is_empty() && first != b'#' {
            names_line = line;
        }
    }
    let raw_line = match lines.next() {
        Some(l) if found => l,
        _ => {
            return Err(format!(
                "ERROR: The weights file {file} has no line of normalized weights followed by a line of raw weights."
            ))
        }
    };
    if !names_line.is_empty() {
        let mut got: Vec<&str> = names_line.split('\t').collect();
        while got
            .last()
            .is_some_and(|c| c.trim_matches([' ', '\r']).is_empty())
        {
            got.pop();
        }
        if let Some(last) = got.last_mut() {
            *last = last.trim_end_matches([' ', '\r']);
        }
        let expected: Vec<&str> = names
            .iter()
            .map(String::as_str)
            .chain(std::iter::once("m0"))
            .collect();
        for ix in 0..got.len().max(expected.len()) {
            if got.get(ix).is_some() && got.get(ix) == expected.get(ix) {
                continue;
            }
            let cell = |v: Option<&&str>| v.map_or("missing".to_string(), |c| format!("\"{c}\""));
            return Err(format!(
                "ERROR: The weights file {file} has {} columns, but the input has {} (its features and m0). Column {} is {} in the weights file and {} in the input. Weights are used by the position of the feature, so the file must come from an input with the same features in the same order.",
                got.len(),
                expected.len(),
                ix + 1,
                cell(got.get(ix)),
                cell(expected.get(ix))
            ));
        }
    }
    let raw: Vec<Option<f64>> = raw_line
        .split_ascii_whitespace()
        .map(|t| t.parse::<f64>().ok().filter(|v| v.is_finite()))
        .collect();
    let bad = raw.iter().any(Option::is_none);
    if bad || raw.len() != n {
        let what = if bad {
            "holds a value that is not a finite number".to_string()
        } else {
            format!("holds {}", raw.len())
        };
        return Err(format!(
            "ERROR: Line {} of the weights file {file} should hold {n} finite raw weights, one per feature and m0, but {what}.",
            line_nr + 1
        ));
    }
    let w: Vec<f64> = raw.into_iter().flatten().collect();
    if n > 1 && w[..n - 1].iter().all(|&v| v == 0.0) {
        return Err(format!(
            "ERROR: Line {} of the weights file {file} gives every feature the weight 0, so these weights cannot order the PSMs.",
            line_nr + 1
        ));
    }
    Ok(Some(w))
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
    let all_sets = (1u8 << FOLDS) - 1;
    let mask: Vec<u8> = fold.iter().map(|&f| all_sets & !(1u8 << f)).collect();
    let path = a.init_weights.as_deref().unwrap_or(Path::new(""));
    let w: Vec<Vec<f64>> = match read_init_weights(path, &pin.feature_names)? {
        Some(ww) => vec![ww; FOLDS],
        None => {
            eprintln!(
                "WARNING: Could not find weights input file {}. Using default weights instead.",
                path.display()
            );
            init_directions(pin, &norm, &mask, FOLDS, a.initial_selection_fdr, false)?
        }
    };
    let mut scores = vec![0.0f64; n];
    for (rows, wf) in test.iter().zip(&w) {
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
        tdc_winners(&pin.meta, scores, |i| pin.text(i).map(|t| t.psm_id))?
    } else {
        (0..scores.len()).collect()
    };
    let s: Vec<f64> = rows.iter().map(|&i| scores[i]).collect();
    let d: Vec<bool> = rows.iter().map(|&i| is_decoy[i]).collect();
    if mode == Mode::Tdc {
        check_separation(&s, &d)?;
    }
    let q = if mode == Mode::MixMax {
        mix_max_q(&s, &d, estimate_pi0(&s, &d))
    } else {
        tdc_q(&s, &d, false)
    };
    let order = psm_order(&pin.meta, &rows, &s, |i| pin.text(i).map(|t| t.psm_id))?;
    let pep = pep_from_q(&order, &s, &q, &d)?;
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
