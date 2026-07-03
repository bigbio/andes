use crate::glycan_db::GlycanComp;
use crate::glycan_mass::CORE_OXONIUM_MZ;

/// Neu5Ac (NeuAc / sialic) oxonium ions: 274.0921 (Neu5Ac−H2O), 292.1027 (Neu5Ac).
pub const NEUAC_OXONIUM_MZ: [f64; 2] = [274.09211, 292.10267];
/// Neu5Gc (NeuGc) oxonium ions: 290.0870 (−H2O), 308.0976.
pub const NEUGC_OXONIUM_MZ: [f64; 2] = [290.08702, 308.09759];

/// Composition-conditioned SIALIC consistency feature (GI-2).
///
/// The core HexNAc oxonium ions are composition-INDEPENDENT (every N-glycan has
/// them), so they cannot separate glycans of different sialic content on one
/// spectrum. This feature CAN: it rewards observed NeuAc/NeuGc oxonium when the
/// candidate glycan CLAIMS that sialic, and penalizes the mismatch (glycan claims
/// a sialic the spectrum lacks, or the spectrum shows a sialic the glycan lacks).
/// Value = ±NeuAc-oxonium (sign by `comp.neuac>0`) ± NeuGc-oxonium (by `comp.neugc>0`),
/// each base-peak-normalised. Additive PIN feature only — never fused into ranking.
pub fn sialic_consistency(peaks: &[(f64, f32)], comp: &GlycanComp, tol_ppm: f64) -> f32 {
    let base = peaks.iter().map(|p| p.1).fold(0.0f32, f32::max).max(1e-9);
    let best_match = |ions: &[f64]| -> f32 {
        let mut acc = 0.0f32;
        for &mz in ions {
            let tol = (mz * tol_ppm / 1e6).max(0.01);
            let mut best = 0.0f32;
            for &(pmz, pi) in peaks {
                if (pmz - mz).abs() <= tol && pi > best {
                    best = pi;
                }
            }
            acc = acc.max(best);
        }
        acc / base
    };
    let neuac_obs = best_match(&NEUAC_OXONIUM_MZ);
    let neugc_obs = best_match(&NEUGC_OXONIUM_MZ);
    let a = if comp.neuac > 0 { neuac_obs } else { -neuac_obs };
    let g = if comp.neugc > 0 { neugc_obs } else { -neugc_obs };
    a + g
}

#[derive(Debug, Clone)]
pub struct OxoniumEvidence {
    pub fired: bool,
    pub summed_frac: f32,
    pub n_core_ions: u8,
}

pub fn oxonium_gate(peaks: &[(f64, f32)], min_frac: f32, tol_ppm: f64) -> OxoniumEvidence {
    let base = peaks.iter().map(|p| p.1).fold(0.0f32, f32::max).max(1e-9);
    let floor = 0.01 * base;
    let mut summed = 0.0f32;
    let mut n = 0u8;
    for &mz in CORE_OXONIUM_MZ.iter() {
        let tol = (mz * tol_ppm / 1e6).max(0.01);
        // best matching peak for this oxonium m/z
        let mut best = 0.0f32;
        for &(pmz, pi) in peaks {
            if (pmz - mz).abs() <= tol && pi > best {
                best = pi;
            }
        }
        if best >= floor {
            summed += best;
            n += 1;
        }
    }
    let frac = summed / base;
    OxoniumEvidence { fired: frac >= min_frac && n >= 2, summed_frac: frac, n_core_ions: n }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn oxonium_gate_fires_on_glyco_spectrum() {
        // base peak intensity 100; two core oxonium ions (204.087, 138.055) at 15 each = 30% summed
        let peaks = vec![(500.0, 100.0), (204.0867, 15.0), (138.055, 15.0), (700.0, 5.0)];
        let e = oxonium_gate(&peaks, 0.10, 20.0);
        assert!(e.fired);
        assert_eq!(e.n_core_ions, 2);
        assert!(e.summed_frac >= 0.29);
    }

    /// GI-2: the sialic consistency feature must DISCRIMINATE — a sialylated
    /// glycan (neuac>0) scores high on a spectrum with NeuAc oxonium, while a
    /// non-sialylated glycan (neuac=0) is penalized on the SAME spectrum.
    #[test]
    fn sialic_consistency_discriminates_by_composition() {
        use crate::glycan_mass::{HEX, HEXNAC, NEUAC};
        // Spectrum carries strong NeuAc oxonium (274.092, 292.103).
        let peaks = vec![(500.0, 100.0), (274.0921, 40.0), (292.1027, 35.0), (204.087, 20.0)];
        let sialylated = GlycanComp { hexnac: 4, hex: 5, fuc: 0, neuac: 2, neugc: 0,
                                      mass: 4.0 * HEXNAC + 5.0 * HEX + 2.0 * NEUAC };
        let no_sialic = GlycanComp { hexnac: 4, hex: 5, fuc: 0, neuac: 0, neugc: 0,
                                     mass: 4.0 * HEXNAC + 5.0 * HEX };
        let s = sialic_consistency(&peaks, &sialylated, 20.0);
        let n = sialic_consistency(&peaks, &no_sialic, 20.0);
        assert!(s > 0.0, "sialylated glycan on a NeuAc-oxonium spectrum → positive, got {s}");
        assert!(n < 0.0, "non-sialylated glycan but NeuAc oxonium present → penalized, got {n}");
        assert!(s > n, "sialic consistency must separate the two glycans");
    }

    #[test]
    fn oxonium_gate_silent_on_nonglyco() {
        let peaks = vec![(500.0, 100.0), (700.0, 5.0), (204.5, 30.0)]; // 204.5 not within tol of 204.0867
        assert!(!oxonium_gate(&peaks, 0.10, 20.0).fired);
    }
}
