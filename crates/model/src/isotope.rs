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

/// Per-atom isotope abundance ladders (`[M, M+1, M+2, ...]`), natural
/// abundances (IUPAC/CIAAW representative values). Index k is the fraction of
/// atoms that are k nominal mass units above the lightest isotope.
const C_ISO: &[f64] = &[0.9893, 0.0107];
const H_ISO: &[f64] = &[0.999_885, 0.000_115];
const N_ISO: &[f64] = &[0.996_36, 0.003_64];
const O_ISO: &[f64] = &[0.997_57, 0.000_38, 0.002_05];
const S_ISO: &[f64] = &[0.9499, 0.0075, 0.0425, 0.0, 0.0001];

/// Theoretical isotope envelope of a molecule from its elemental formula
/// (C, H, N, O, S atom counts), as the first `n_isotopes` nominal-mass peaks
/// normalised to sum 1.0. Exact per-element multinomial convolution (no
/// averagine, no Poisson approximation): every atom's own ladder is folded
/// in, so the O18 and S34 contributions to M+2 are kept — they matter for
/// glycans, whose oxygen fraction is roughly twice a peptide's.
///
/// Empty when `n_isotopes == 0`. A formula with no atoms returns a pure
/// monoisotope (`[1, 0, 0, ...]`).
pub fn isotope_envelope_from_formula(
    c: u32,
    h: u32,
    n: u32,
    o: u32,
    s: u32,
    n_isotopes: usize,
) -> Vec<f64> {
    if n_isotopes == 0 {
        return Vec::new();
    }
    let mut env = vec![0.0f64; n_isotopes];
    env[0] = 1.0;
    let mut scratch = vec![0.0f64; n_isotopes];
    for (ladder, count) in [(C_ISO, c), (H_ISO, h), (N_ISO, n), (O_ISO, o), (S_ISO, s)] {
        for _ in 0..count {
            scratch.iter_mut().for_each(|v| *v = 0.0);
            for (i, &e) in env.iter().enumerate() {
                if e == 0.0 {
                    continue;
                }
                for (k, &a) in ladder.iter().enumerate() {
                    if i + k >= n_isotopes {
                        break;
                    }
                    scratch[i + k] += e * a;
                }
            }
            std::mem::swap(&mut env, &mut scratch);
        }
    }
    let sum: f64 = env.iter().sum();
    if sum > 0.0 {
        for v in &mut env {
            *v /= sum;
        }
    }
    env
}

/// Approximate elemental formula of an N-glycopeptide of neutral mass `mass`
/// (Da), as (C, H, N, O, S) atom counts.
///
/// A glycopeptide is not averagine: the glycan half carries ~14% fewer
/// carbons per Da and about twice the oxygen, so the peptide averagine
/// (`averagine_isotope_envelope`) overstates M+1 and understates the O18
/// share of M+2. This blends the peptide averagine residue
/// (C4.9384 H7.7583 N1.3577 O1.4773 S0.0417 per 111.1254 Da) with an
/// average complex N-glycan (HexNAc4 Hex5 Fuc1 NeuAc1 = C79 H129 N5 O57 per
/// 2059.73 Da) at equal mass share — the typical split for a tryptic
/// N-glycopeptide in the 2.5–5 kDa range where the precursor
/// mono-correction operates. Counts are rounded to whole atoms.
pub fn glycopeptide_formula(mass: f64) -> (u32, u32, u32, u32, u32) {
    // Per-Da atom densities of the two halves.
    const PEP: [f64; 5] = [
        4.9384 / 111.1254,
        7.7583 / 111.1254,
        1.3577 / 111.1254,
        1.4773 / 111.1254,
        0.0417 / 111.1254,
    ];
    const GLY: [f64; 5] = [
        79.0 / 2059.73,
        129.0 / 2059.73,
        5.0 / 2059.73,
        57.0 / 2059.73,
        0.0,
    ];
    const GLYCAN_MASS_SHARE: f64 = 0.5;
    let m = mass.max(0.0);
    let at = |i: usize| -> u32 {
        (m * (PEP[i] * (1.0 - GLYCAN_MASS_SHARE) + GLY[i] * GLYCAN_MASS_SHARE)).round() as u32
    };
    (at(0), at(1), at(2), at(3), at(4))
}

/// Theoretical precursor isotope envelope of an N-glycopeptide of neutral
/// mass `mass`, first `n_isotopes` peaks normalised to sum 1.0. See
/// [`glycopeptide_formula`] for the element model.
pub fn glycopeptide_isotope_envelope(mass: f64, n_isotopes: usize) -> Vec<f64> {
    let (c, h, n, o, s) = glycopeptide_formula(mass);
    isotope_envelope_from_formula(c, h, n, o, s, n_isotopes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn averagine_envelope_is_normalized_and_decays_for_small_peptide() {
        let env = averagine_isotope_envelope(1000.0, 4);
        assert_eq!(env.len(), 4);
        let sum: f64 = env.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-9,
            "envelope must sum to 1.0, got {sum}"
        );
        assert!(
            env[0] > env[1] && env[1] > env[2],
            "a ~1000 Da peptide's envelope should decay from the monoisotope: {env:?}"
        );
    }

    #[test]
    fn averagine_plus_one_grows_with_mass() {
        // The +1 isotope's relative height increases with peptide mass
        // (more carbons -> higher 13C probability).
        let small = averagine_isotope_envelope(1000.0, 4);
        let large = averagine_isotope_envelope(3000.0, 4);
        let ratio_small = small[1] / small[0];
        let ratio_large = large[1] / large[0];
        assert!(
            ratio_large > ratio_small,
            "+1/+0 ratio should grow with mass: small {ratio_small} vs large {ratio_large}"
        );
    }

    #[test]
    fn averagine_handles_zero_and_one_isotope_requests() {
        assert!(averagine_isotope_envelope(1000.0, 0).is_empty());
        let one = averagine_isotope_envelope(1000.0, 1);
        assert_eq!(one.len(), 1);
        assert!(
            (one[0] - 1.0).abs() < 1e-9,
            "single-isotope envelope is all monoisotope"
        );
    }

    #[test]
    fn formula_envelope_is_normalized_and_matches_known_shape() {
        // Bovine insulin-like size check is overkill; use a peptide of known
        // formula: angiotensin II DRVYIHPF = C50 H71 N13 O12, M0 1045.53.
        // First-order expectation: M+1/M0 ~ 50*0.0107 + 13*0.00364 + 71*0.000115
        // + 12*0.00038 = 0.60; M+2/M0 ~ 0.20 (13C pairs + O18).
        let env = isotope_envelope_from_formula(50, 71, 13, 12, 0, 5);
        let sum: f64 = env.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9, "sum {sum}");
        let r1 = env[1] / env[0];
        let r2 = env[2] / env[0];
        assert!((r1 - 0.60).abs() < 0.02, "M+1/M0 {r1}");
        assert!((r2 - 0.20).abs() < 0.02, "M+2/M0 {r2}");
    }

    #[test]
    fn formula_envelope_degenerate_inputs() {
        assert!(isotope_envelope_from_formula(10, 10, 1, 1, 0, 0).is_empty());
        let none = isotope_envelope_from_formula(0, 0, 0, 0, 0, 3);
        assert_eq!(none, vec![1.0, 0.0, 0.0]);
    }

    #[test]
    fn glycopeptide_envelope_has_fewer_carbons_than_averagine() {
        // Same mass, same convolution: the glycopeptide model must put LESS
        // weight on M+1 relative to M0 than an all-peptide averagine formula
        // (fewer carbons and far less nitrogen per Da).
        let m = 3500.0;
        let g = glycopeptide_isotope_envelope(m, 6);
        let pep = |d: f64| (m * d / 111.1254).round() as u32;
        let a = isotope_envelope_from_formula(
            pep(4.9384),
            pep(7.7583),
            pep(1.3577),
            pep(1.4773),
            pep(0.0417),
            6,
        );
        assert!(
            g[1] / g[0] < a[1] / a[0],
            "glyco M+1/M0 {} should be below peptide {}",
            g[1] / g[0],
            a[1] / a[0]
        );
        let (c, h, n, o, s) = glycopeptide_formula(m);
        assert!(
            c > 100 && h > 200 && n > 10 && o > 50 && s <= 2,
            "{c} {h} {n} {o} {s}"
        );
    }
}
