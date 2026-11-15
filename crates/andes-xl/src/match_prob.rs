//! Label-blind match probability of a peptide's fragment matches.
//!
//! A peptide with `k` theoretical fragment ions matches `n` of them. If each
//! ion hit a peak by chance with probability `p` (the fraction of the
//! spectrum's m/z range covered by the tolerance windows of its peaks), the
//! number of chance matches is Binomial(`k`, `p`). The score is
//! `-log10 P(X >= n)`, computed in the log domain. It depends on the peptide's
//! length (through `k`), the spectrum's peak count and the tolerance only: no
//! decoy label and no dataset statistic enters it.

/// Chance that one theoretical ion lands within `tol` (Da) of one of
/// `n_peaks` peaks spread over `mz_range` (Da), clamped to `[1e-9, 1]`.
pub fn chance_match_probability(n_peaks: usize, tol: f64, mz_range: f64) -> f64 {
    if mz_range <= 0.0 {
        return 1.0;
    }
    (n_peaks as f64 * 2.0 * tol / mz_range).clamp(1e-9, 1.0)
}

fn ln_choose(n: u32, k: u32) -> f64 {
    ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
}

fn ln_factorial(n: u32) -> f64 {
    (2..=n).map(|i| (i as f64).ln()).sum()
}

/// `-log10 P(X >= n)` for `X ~ Binomial(k, p)`; 0 when `n == 0`.
pub fn binomial_tail_score(k: u32, n: u32, p: f64) -> f64 {
    if n == 0 || k == 0 {
        return 0.0;
    }
    let n = n.min(k);
    if p >= 1.0 {
        return 0.0;
    }
    let lp = p.ln();
    let lq = (-p).ln_1p();
    // log-sum-exp over the tail terms i = n..=k.
    let terms: Vec<f64> = (n..=k)
        .map(|i| ln_choose(k, i) + i as f64 * lp + (k - i) as f64 * lq)
        .collect();
    let m = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let ln_tail = m + terms.iter().map(|t| (t - m).exp()).sum::<f64>().ln();
    (-ln_tail / std::f64::consts::LN_10).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tail_matches_closed_forms() {
        // P(X >= 1) = 1 - (1-p)^k.
        let (k, p) = (10, 0.1);
        let want = -(1.0 - (1.0_f64 - p).powi(k as i32)).log10();
        assert!((binomial_tail_score(k, 1, p) - want).abs() < 1e-9);
        // P(X >= k) = p^k.
        assert!((binomial_tail_score(5, 5, 0.5) - 5.0 * 2.0_f64.log10()).abs() < 1e-9);
        assert_eq!(binomial_tail_score(10, 0, 0.3), 0.0);
    }

    #[test]
    fn longer_peptides_need_more_matches_for_the_same_score() {
        let p = 0.05;
        // Same 6 matches: more theoretical ions (a longer peptide) is less surprising.
        assert!(binomial_tail_score(20, 6, p) > binomial_tail_score(60, 6, p));
        // Tiny tails stay finite in the log domain.
        let s = binomial_tail_score(200, 150, 1e-6);
        assert!(s.is_finite() && s > 500.0);
    }

    #[test]
    fn chance_probability_is_clamped() {
        assert!((chance_match_probability(100, 0.01, 2000.0) - 0.001).abs() < 1e-12);
        assert_eq!(chance_match_probability(1_000_000, 0.5, 100.0), 1.0);
        assert_eq!(chance_match_probability(10, 0.01, 0.0), 1.0);
    }
}
