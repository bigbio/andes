use crate::cv::CvResult;
use crate::pep::too_good_separation;
use crate::scores::{descending_order, tdc_q};
use rayon::prelude::*;

pub fn normalize_fold_scores(scores: &[f64], is_decoy: &[bool], fdr: f64) -> (f64, f64) {
    if scores.is_empty() {
        return (0.0, 0.0);
    }
    let q = tdc_q(scores, is_decoy, false);
    let order = descending_order(scores);
    let median_index = is_decoy.iter().filter(|&&d| d).count() / 2;
    let mut fdr_score = scores[order[0]];
    let mut median_decoy = fdr_score + 1.0;
    let mut decoys = 0usize;
    for &i in &order {
        if q[i] < fdr {
            fdr_score = scores[i];
        }
        if is_decoy[i] {
            decoys += 1;
            if decoys == median_index {
                median_decoy = scores[i];
                break;
            }
        }
    }
    (fdr_score, fdr_score - median_decoy)
}

fn rescale(w: &mut [f64], sub: f64, scale: f64) {
    let nf = w.len() - 1;
    for v in w[..nf].iter_mut() {
        *v /= scale;
    }
    w[nf] = (w[nf] - sub) / scale;
}

fn fmt_pct(v: f64) -> String {
    crate::write::fmt_g(v, 6)
}

fn fold_scale(
    rows: &[usize],
    scores: &[f64],
    is_decoy: &[bool],
    fdr: f64,
) -> Result<(f64, f64), String> {
    if rows.is_empty() {
        return Err("Error: no scored PSMs were provided.".to_string());
    }
    let s: Vec<f64> = rows.iter().map(|&i| scores[i]).collect();
    let d: Vec<bool> = rows.iter().map(|&i| is_decoy[i]).collect();
    if too_good_separation(&s, &d) {
        return Err(
            "Error in the input data: too good separation between target and decoy PSMs.\nTerminating."
                .to_string(),
        );
    }
    let (sub, scale) = normalize_fold_scores(&s, &d, fdr);
    if scale <= 0.0 {
        return Err(format!(
            "Error: median decoy score <= score at {}% FDR. Cannot rescale scores to merge cross validation bins, try lowering --trainFDR.",
            fmt_pct(fdr * 100.0)
        ));
    }
    Ok((sub, scale))
}

pub fn merge_folds(cv: &mut CvResult, is_decoy: &[bool], fdr: f64) -> Result<(), String> {
    let fold = &cv.fold;
    let rows_of: Vec<Vec<usize>> = (0..cv.w_norm.len())
        .into_par_iter()
        .map(|f| (0..fold.len()).filter(|&i| fold[i] as usize == f).collect())
        .collect();
    let scores = &cv.test_scores;
    let fits: Vec<Result<(f64, f64), String>> = rows_of
        .par_iter()
        .map(|rows| fold_scale(rows, scores, is_decoy, fdr))
        .collect();
    for (f, (rows, fit)) in rows_of.iter().zip(fits).enumerate() {
        let (sub, scale) = fit?;
        for &i in rows {
            cv.test_scores[i] = (cv.test_scores[i] - sub) / scale;
        }
        rescale(&mut cv.w_norm[f], sub, scale);
        rescale(&mut cv.w_raw[f], sub, scale);
    }
    Ok(())
}
