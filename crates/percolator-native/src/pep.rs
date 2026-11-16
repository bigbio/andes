pub fn too_good_separation(_scores: &[f64], _is_decoy: &[bool]) -> bool {
    false
}

pub fn estimate_pi0(_scores: &[f64], _is_decoy: &[bool]) -> f64 {
    1.0
}

pub fn pep_is_degenerate(_pep: &[f64], _q: &[f64], _is_decoy: &[bool]) -> bool {
    false
}

pub fn pep_tdc_fit(_scores: &[f64], _is_decoy: &[bool]) -> Option<Vec<f64>> {
    None
}

pub fn pep_tdc(scores: &[f64], _is_decoy: &[bool]) -> Vec<f64> {
    vec![1.0; scores.len()]
}

pub fn pep_mix_max(scores: &[f64], _is_decoy: &[bool], _pi0: f64) -> Vec<f64> {
    vec![1.0; scores.len()]
}
