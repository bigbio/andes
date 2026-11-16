use crate::cv::CvResult;

pub fn normalize_fold_scores(_scores: &[f64], _is_decoy: &[bool], _fdr: f64) -> (f64, f64) {
    (0.0, 1.0)
}

pub fn merge_folds(_cv: &mut CvResult, _is_decoy: &[bool], _fdr: f64) -> Result<(), String> {
    Ok(())
}
