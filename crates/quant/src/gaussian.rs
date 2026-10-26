//! The Gaussian profile `g(t) = exp(-(t - t0)^2 / (2 w^2))` of the evidence
//! fit, evaluated in bulk.
//!
//! The grid search of [`crate::evidence`] evaluates the profile of every
//! candidate apex and width over the scans around it, tens of thousands of
//! exponentials per probe. [`exp`] is a double-precision exponential written
//! so that the compiler vectorises a loop over it: every operation is an
//! IEEE add, multiply or compare-and-select on one lane (no fused multiply-add,
//! no lookup table), so the result of an element does not depend on the
//! vector width, on the CPU's instruction set or on the element's position in
//! a slice. A target with AVX2 runs the same loop four lanes wide.
//!
//! The result is within 1e-15 (relative) of the C library's `exp` on the
//! domain the fit uses (`x <= 0`), not bit for bit the same.

/// `ln 2` split in a high part with 21 trailing zero bits (so `k * LN2_HI`
/// is exact for `|k| < 2^21`) and the remainder.
const LN2_HI: f64 = f64::from_bits(0x3fe6_2e42_fee0_0000);
const LN2_LO: f64 = f64::from_bits(0x3dea_39ef_3579_3c76);
const LOG2_E: f64 = std::f64::consts::LOG2_E;
/// `1.5 * 2^52`: adding it rounds a value below `2^51` to an integer, which
/// then sits in the low mantissa bits.
const ROUND: f64 = 6_755_399_441_055_744.0;
/// Below this `exp` is 0 here (the C library's is subnormal or 0).
const UNDERFLOW: f64 = -708.0;
/// Above this `exp` is infinite here.
const OVERFLOW: f64 = 709.0;

/// `1 / k!` for `k = 0..=13`: with `|r| <= ln(2) / 2` the Taylor remainder
/// is below `5e-18`.
const COEFF: [f64; 14] = [
    1.0,
    1.0,
    1.0 / 2.0,
    1.0 / 6.0,
    1.0 / 24.0,
    1.0 / 120.0,
    1.0 / 720.0,
    1.0 / 5_040.0,
    1.0 / 40_320.0,
    1.0 / 362_880.0,
    1.0 / 3_628_800.0,
    1.0 / 39_916_800.0,
    1.0 / 479_001_600.0,
    1.0 / 6_227_020_800.0,
];

/// `exp(x)`; NaN for a NaN, 0 below [`UNDERFLOW`], infinity above
/// [`OVERFLOW`].
#[inline(always)]
pub fn exp(x: f64) -> f64 {
    // The argument is clamped for the reduction; the ends are put back below.
    let c = x.clamp(UNDERFLOW, OVERFLOW);
    // k = round(c / ln 2) as an integer in the low bits of `kf`.
    let kf = c * LOG2_E + ROUND;
    let k = kf - ROUND;
    let r = (c - k * LN2_HI) - k * LN2_LO;
    // Estrin evaluation: independent chains for the vector units.
    let r2 = r * r;
    let r4 = r2 * r2;
    let r8 = r4 * r4;
    let p01 = COEFF[0] + COEFF[1] * r;
    let p23 = COEFF[2] + COEFF[3] * r;
    let p45 = COEFF[4] + COEFF[5] * r;
    let p67 = COEFF[6] + COEFF[7] * r;
    let p89 = COEFF[8] + COEFF[9] * r;
    let pab = COEFF[10] + COEFF[11] * r;
    let pcd = COEFF[12] + COEFF[13] * r;
    let p0123 = p01 + p23 * r2;
    let p4567 = p45 + p67 * r2;
    let p89ab = p89 + pab * r2;
    let p = (p0123 + p4567 * r4) + (p89ab + pcd * r4) * r8;
    // 2^k from the exponent field: the low 12 bits of `kf` are k mod 4096.
    let scale = f64::from_bits((kf.to_bits() << 52).wrapping_add(0x3ff0_0000_0000_0000));
    let y = p * scale;
    if x < UNDERFLOW {
        0.0
    } else if x > OVERFLOW {
        f64::INFINITY
    } else {
        y
    }
}

/// `exp(-(t - t0)^2 / (2 w^2))`.
#[inline(always)]
pub fn gaussian(t: f64, t0: f64, w: f64) -> f64 {
    let d = (t - t0) / w;
    exp(-0.5 * d * d)
}

#[inline(always)]
fn gaussians_loop(ts: &[f64], t0: f64, w: f64, out: &mut [f64]) {
    for (o, &t) in out.iter_mut().zip(ts) {
        *o = gaussian(t, t0, w);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn gaussians_avx2(ts: &[f64], t0: f64, w: f64, out: &mut [f64]) {
    gaussians_loop(ts, t0, w, out);
}

/// `out[j] = gaussian(ts[j], t0, w)`; `out` has the length of `ts`.
pub fn gaussians(ts: &[f64], t0: f64, w: f64, out: &mut [f64]) {
    assert_eq!(ts.len(), out.len());
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: the CPU supports AVX2, as just detected.
        unsafe { gaussians_avx2(ts, t0, w, out) };
        return;
    }
    gaussians_loop(ts, t0, w, out);
}

/// Appends `gaussian(t, t0, w)` for every `t` of `ts` to `out`.
pub fn extend(out: &mut Vec<f64>, ts: &[f64], t0: f64, w: f64) {
    let start = out.len();
    out.resize(start + ts.len(), 0.0);
    gaussians(ts, t0, w, &mut out[start..]);
}

/// Scans per block of [`sums`].
const BLOCK: usize = 64;

/// `(sum z_j g_j, sum g_j, sum g_j^2)` over `g_j = gaussian(ts[j], t0, w)`,
/// each sum in the order of `ts`; `z` has the length of `ts`.
pub fn sums(ts: &[f64], z: &[f64], t0: f64, w: f64) -> (f64, f64, f64) {
    assert_eq!(ts.len(), z.len());
    let (mut num, mut sg, mut sgg) = (0.0, 0.0, 0.0);
    let mut g = [0.0f64; BLOCK];
    for (t, z) in ts.chunks(BLOCK).zip(z.chunks(BLOCK)) {
        let g = &mut g[..t.len()];
        gaussians(t, t0, w, g);
        for (&gj, &zj) in g.iter().zip(z) {
            num += zj * gj;
            sg += gj;
            sgg += gj * gj;
        }
    }
    (num, sg, sgg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exp_matches_the_c_library_to_1e_15() {
        let mut worst = 0.0f64;
        let mut x = -745.0;
        while x <= 0.0 {
            let (a, b) = (exp(x), x.exp());
            if x >= UNDERFLOW {
                worst = worst.max(((a - b) / b).abs());
            }
            x += 0.003_7;
        }
        assert!(worst < 1e-15, "worst relative error {worst:e}");
        assert_eq!(exp(0.0), 1.0);
        assert!((exp(1.0) - std::f64::consts::E).abs() < 1e-15);
        assert!((exp(20.0) / 20.0f64.exp() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn exp_ends_and_non_numbers() {
        assert!(exp(f64::NAN).is_nan());
        assert_eq!(exp(f64::NEG_INFINITY), 0.0);
        assert_eq!(exp(-1000.0), 0.0);
        assert_eq!(exp(f64::INFINITY), f64::INFINITY);
        assert_eq!(exp(710.0), f64::INFINITY);
        assert!(exp(-707.9) > 0.0);
    }

    #[test]
    fn bulk_elements_equal_the_scalar_whatever_the_width_and_position() {
        let ts: Vec<f64> = (0..203).map(|j| 10.0 + 0.37 * j as f64).collect();
        let (t0, w) = (47.3, 2.9);
        let mut out = vec![0.0; ts.len()];
        gaussians(&ts, t0, w, &mut out);
        for (&t, &g) in ts.iter().zip(&out) {
            assert_eq!(g.to_bits(), gaussian(t, t0, w).to_bits());
        }
        // A sub-slice shifts every element's lane.
        let mut sub = vec![0.0; 150];
        gaussians(&ts[3..153], t0, w, &mut sub);
        assert!(sub
            .iter()
            .zip(&out[3..153])
            .all(|(a, b)| a.to_bits() == b.to_bits()));
        let mut v = vec![1.5];
        extend(&mut v, &ts, t0, w);
        assert_eq!(v.len(), 1 + ts.len());
        assert_eq!(v[0], 1.5);
        assert!(v[1..]
            .iter()
            .zip(&out)
            .all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    #[test]
    fn sums_follow_the_order_of_the_scans() {
        let ts: Vec<f64> = (0..150).map(|j| 0.5 * j as f64).collect();
        let z: Vec<f64> = (0..150).map(|j| ((j * 7) % 11) as f64 - 5.0).collect();
        let (t0, w) = (31.0, 4.0);
        let (mut num, mut sg, mut sgg) = (0.0, 0.0, 0.0);
        for (&t, &zj) in ts.iter().zip(&z) {
            let g = gaussian(t, t0, w);
            num += zj * g;
            sg += g;
            sgg += g * g;
        }
        let got = sums(&ts, &z, t0, w);
        assert_eq!(
            (got.0.to_bits(), got.1.to_bits(), got.2.to_bits()),
            (num.to_bits(), sg.to_bits(), sgg.to_bits())
        );
    }

    #[test]
    #[ignore = "timing"]
    fn timing_against_the_c_library() {
        let ts: Vec<f64> = (0..1 << 16).map(|j| 0.01 * j as f64).collect();
        let mut out = vec![0.0; ts.len()];
        let reps = 400;
        let t = std::time::Instant::now();
        for _ in 0..reps {
            gaussians(&ts, 300.0, 3.0, &mut out);
            std::hint::black_box(&out);
        }
        let ours = t.elapsed().as_secs_f64() / (reps * ts.len()) as f64 * 1e9;
        let t = std::time::Instant::now();
        for _ in 0..reps {
            for (o, &t) in out.iter_mut().zip(&ts) {
                let d = (t - 300.0) / 3.0;
                *o = (-0.5 * d * d).exp();
            }
            std::hint::black_box(&out);
        }
        let libm = t.elapsed().as_secs_f64() / (reps * ts.len()) as f64 * 1e9;
        eprintln!("gaussian: {ours:.2} ns/element, libm exp: {libm:.2} ns/element");
    }
}
