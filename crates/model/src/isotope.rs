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

/// The highest isotope peak (0 is the monoisotopic one) an instrument may have isolated
/// for an analyte of neutral `mass`. Data-dependent selection takes the most intense peak
/// of the envelope, M+1 or M+2 for a heavy analyte, so a selected precursor may sit on any
/// peak from the monoisotopic one to one above the averagine envelope's most abundant (the
/// next peak is nearly as intense when the mass sits near a change of the most abundant).
pub fn max_isolated_isotope(mass: f64) -> usize {
    averagine_isotope_envelope(mass, 8)
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map_or(0, |(k, _)| k)
        + 1
}

/// The isotope peak, from the monoisotopic one (0) to `max_isotope`, that a selected
/// precursor at `mz` sits on for an analyte of neutral `mass` at `charge`, within `tol_th`
/// Th; `None` when it sits on none of them or `charge` is 0.
pub fn isolated_isotope_peak(
    mz: f64,
    mass: f64,
    charge: u8,
    max_isotope: usize,
    tol_th: f64,
) -> Option<usize> {
    if charge == 0 {
        return None;
    }
    let z = charge as f64;
    (0..=max_isotope).find(|&k| {
        let target = (mass + k as f64 * crate::mass::ISOTOPE) / z + crate::mass::PROTON;
        (mz - target).abs() <= tol_th
    })
}

/// Averagine atoms per dalton of neutral mass: C, H, N, O, S.
const AVERAGINE_ATOMS_PER_DA: [f64; 5] = [
    4.9384 / 111.1254,
    7.7583 / 111.1254,
    1.3577 / 111.1254,
    1.4773 / 111.1254,
    0.0417 / 111.1254,
];

/// Expected distance in Da between the M+1 and M+2 peaks of an averagine peptide of neutral
/// mass `mass`: the intensity-weighted mean mass of the M+2 isotopologues minus that of the
/// M+1 isotopologues. M+1 comes from one heavy atom (¹³C, ²H, ¹⁵N, ¹⁷O, ³³S), M+2 from two
/// of them or from one ¹⁸O or ³⁴S, each weighted by its natural abundance and the number
/// of atoms of that element (Poisson limit of the isotope convolution). Returns
/// `mass::ISOTOPE` when `mass` is not positive.
pub fn averagine_m1_m2_spacing(mass: f64) -> f64 {
    use crate::mass::{C, C13, H, H2, N, N15, O, O17, O18, S, S33, S34};
    if mass <= 0.0 {
        return crate::mass::ISOTOPE;
    }
    let atoms = |i: usize| mass * AVERAGINE_ATOMS_PER_DA[i];
    let rate = |count: f64, ladder: &[f64], k: usize| count * ladder[k] / ladder[0];
    let one: [(f64, f64); 5] = [
        (rate(atoms(0), C_ISO, 1), C13 - C),
        (rate(atoms(1), H_ISO, 1), H2 - H),
        (rate(atoms(2), N_ISO, 1), N15 - N),
        (rate(atoms(3), O_ISO, 1), O17 - O),
        (rate(atoms(4), S_ISO, 1), S33 - S),
    ];
    let two: [(f64, f64); 2] = [
        (rate(atoms(3), O_ISO, 2), O18 - O),
        (rate(atoms(4), S_ISO, 2), S34 - S),
    ];
    let weight_1: f64 = one.iter().map(|(r, _)| r).sum();
    let shift_1: f64 = one.iter().map(|(r, d)| r * d).sum();
    let mut weight_2 = 0.0;
    let mut shift_2 = 0.0;
    for (i, &(ri, di)) in one.iter().enumerate() {
        for (j, &(rj, dj)) in one.iter().enumerate().skip(i) {
            let w = if i == j { ri * ri / 2.0 } else { ri * rj };
            weight_2 += w;
            shift_2 += w * (di + dj);
        }
    }
    for &(r, d) in &two {
        weight_2 += r;
        shift_2 += r * d;
    }
    shift_2 / weight_2 - shift_1 / weight_1
}

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
    fn m1_m2_spacing_stays_below_the_carbon_spacing_and_above_the_sulfur_shift() {
        for mass in [300.0, 600.0, 1200.0, 2400.0, 4000.0] {
            let d = averagine_m1_m2_spacing(mass);
            assert!(
                d > 0.999 && d < crate::mass::ISOTOPE,
                "mass {mass}: spacing {d}"
            );
        }
        assert_eq!(averagine_m1_m2_spacing(0.0), crate::mass::ISOTOPE);
    }

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
    fn a_heavy_analyte_may_be_isolated_on_a_later_isotope_peak() {
        assert_eq!(max_isolated_isotope(600.0), 1);
        assert_eq!(max_isolated_isotope(2550.0), 2);
        assert!(max_isolated_isotope(5000.0) >= 3);
        let mass = 2550.26;
        let mono = mass / 2.0 + crate::mass::PROTON;
        let m1 = (mass + crate::mass::ISOTOPE) / 2.0 + crate::mass::PROTON;
        let max = max_isolated_isotope(mass);
        assert_eq!(isolated_isotope_peak(mono, mass, 2, max, 0.01), Some(0));
        assert_eq!(isolated_isotope_peak(m1, mass, 2, max, 0.01), Some(1));
        assert_eq!(isolated_isotope_peak(m1 + 0.02, mass, 2, max, 0.01), None);
        assert_eq!(isolated_isotope_peak(mono, mass, 3, max, 0.01), None);
        assert_eq!(isolated_isotope_peak(mono, mass, 0, max, 0.01), None);
        let small = 600.0;
        let m3 = (small + 3.0 * crate::mass::ISOTOPE) / 2.0 + crate::mass::PROTON;
        assert_eq!(
            isolated_isotope_peak(m3, small, 2, max_isolated_isotope(small), 0.01),
            None
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
