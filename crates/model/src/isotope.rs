pub fn averagine_isotope_envelope(mass: f64, n_isotopes: usize) -> Vec<f64> {
    if n_isotopes == 0 {
        return Vec::new();
    }
    let lambda = (mass * 4.76e-4).max(0.0);
    let mut env = Vec::with_capacity(n_isotopes);
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

const C_ISO: &[f64] = &[0.9893, 0.0107];
const H_ISO: &[f64] = &[0.999_885, 0.000_115];
const N_ISO: &[f64] = &[0.996_36, 0.003_64];
const O_ISO: &[f64] = &[0.997_57, 0.000_38, 0.002_05];
const S_ISO: &[f64] = &[0.9499, 0.0075, 0.0425, 0.0, 0.0001];

pub fn isotope_envelope_from_formula(
    c: u32,
    h: u32,
    n: u32,
    o: u32,
    s: u32,
    n_isotopes: usize,
) -> Vec<f64> {
    isotope_envelope_and_share_from_formula(c, h, n, o, s, n_isotopes).0
}

pub fn isotope_envelope_and_share_from_formula(
    c: u32,
    h: u32,
    n: u32,
    o: u32,
    s: u32,
    n_isotopes: usize,
) -> (Vec<f64>, f64) {
    if n_isotopes == 0 {
        return (Vec::new(), 0.0);
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
    (env, sum)
}

pub fn glycopeptide_formula(mass: f64) -> (u32, u32, u32, u32, u32) {
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

pub fn glycopeptide_isotope_envelope(mass: f64, n_isotopes: usize) -> Vec<f64> {
    let (c, h, n, o, s) = glycopeptide_formula(mass);
    isotope_envelope_from_formula(c, h, n, o, s, n_isotopes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn share_is_the_fraction_of_the_whole_pattern_in_the_kept_peaks() {
        let cases = [
            ((5, 10, 2, 3, 0), 0.999_977_287_454, 0.999_999_081_304),
            ((44, 70, 12, 13, 0), 0.995_235_617_258, 0.999_268_002_862),
            ((89, 140, 24, 27, 1), 0.948_998_300_675, 0.984_894_597_691),
            ((133, 209, 37, 40, 1), 0.871_108_790_480, 0.950_089_442_060),
        ];
        for ((c, h, n, o, s), share4, share5) in cases {
            let (env4, got4) = isotope_envelope_and_share_from_formula(c, h, n, o, s, 4);
            let (_, got5) = isotope_envelope_and_share_from_formula(c, h, n, o, s, 5);
            assert!((got4 - share4).abs() < 1e-10, "C{c}: {got4} vs {share4}");
            assert!((got5 - share5).abs() < 1e-10, "C{c}: {got5} vs {share5}");
            let (env60, all) = isotope_envelope_and_share_from_formula(c, h, n, o, s, 60);
            assert!((all - 1.0).abs() < 1e-12, "C{c}: {all}");
            for k in 0..4 {
                assert!(
                    (env4[k] * got4 - env60[k] * all).abs() < 1e-12,
                    "C{c} peak {k}"
                );
            }
            assert!((env4.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            assert_eq!(env4, isotope_envelope_from_formula(c, h, n, o, s, 4));
        }
        let (env, share) = isotope_envelope_and_share_from_formula(89, 140, 24, 27, 1, 4);
        assert!((env[1] * share - 0.333_883_045_9).abs() < 1e-9);
        assert_eq!(
            isotope_envelope_and_share_from_formula(1, 1, 1, 1, 1, 0),
            (Vec::new(), 0.0)
        );
        assert_eq!(
            isotope_envelope_and_share_from_formula(0, 0, 0, 0, 0, 3),
            (vec![1.0, 0.0, 0.0], 1.0)
        );
    }

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
