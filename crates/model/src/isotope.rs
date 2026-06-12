//! Theoretical isotope envelope (averagine approximation).
//!
//! Used by the chimeric MS1 filter to compare a peptide's expected precursor
//! isotope distribution against the observed MS1 envelope.

/// Approximate the relative intensities of the first `n_isotopes` peaks of a
/// peptide's precursor isotope envelope from its neutral monoisotopic `mass`,
/// using the averagine + Poisson model. Returns a vector of length
/// `n_isotopes` normalized to sum 1.0 (empty when `n_isotopes == 0`).
///
/// Model: an averagine residue is ~111.1 Da and contains ~4.94 carbons, so a
/// peptide of `mass` Da has roughly `mass / 111.1 * 4.94 ≈ mass * 0.0445`
/// carbons. With natural 13C abundance 1.07%, the expected number of 13C atoms
/// is `lambda ≈ carbons * 0.0107 ≈ mass * 4.76e-4`. The isotope peak
/// intensities follow a Poisson distribution in the number of 13C atoms:
/// `p_k = e^-lambda * lambda^k / k!`. (This ignores other heavy isotopes
/// — N15/O18/S34/S33 — which is the standard averagine first approximation and
/// is sufficient for an envelope-shape match.)
pub fn averagine_isotope_envelope(mass: f64, n_isotopes: usize) -> Vec<f64> {
    if n_isotopes == 0 {
        return Vec::new();
    }
    let lambda = (mass * 4.76e-4).max(0.0);
    let mut env = Vec::with_capacity(n_isotopes);
    // p_k = e^-lambda * lambda^k / k!, computed iteratively (p_0 = e^-lambda).
    let mut p = (-lambda).exp();
    for k in 0..n_isotopes {
        env.push(p);
        p *= lambda / (k as f64 + 1.0);
    }
    let sum: f64 = env.iter().sum();
    if sum > 0.0 {
        for v in &mut env {
            *v /= sum;
        }
    }
    env
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn averagine_envelope_is_normalized_and_decays_for_small_peptide() {
        let env = averagine_isotope_envelope(1000.0, 4);
        assert_eq!(env.len(), 4);
        let sum: f64 = env.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9, "envelope must sum to 1.0, got {sum}");
        assert!(env[0] > env[1] && env[1] > env[2],
            "a ~1000 Da peptide's envelope should decay from the monoisotope: {env:?}");
    }

    #[test]
    fn averagine_plus_one_grows_with_mass() {
        // The +1 isotope's relative height increases with peptide mass
        // (more carbons -> higher 13C probability).
        let small = averagine_isotope_envelope(1000.0, 4);
        let large = averagine_isotope_envelope(3000.0, 4);
        let ratio_small = small[1] / small[0];
        let ratio_large = large[1] / large[0];
        assert!(ratio_large > ratio_small,
            "+1/+0 ratio should grow with mass: small {ratio_small} vs large {ratio_large}");
    }

    #[test]
    fn averagine_handles_zero_and_one_isotope_requests() {
        assert!(averagine_isotope_envelope(1000.0, 0).is_empty());
        let one = averagine_isotope_envelope(1000.0, 1);
        assert_eq!(one.len(), 1);
        assert!((one[0] - 1.0).abs() < 1e-9, "single-isotope envelope is all monoisotope");
    }
}
