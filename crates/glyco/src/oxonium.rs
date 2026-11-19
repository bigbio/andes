use crate::glycan_db::GlycanComp;
use crate::glycan_first::PPM_FLOOR_DA;
use crate::glycan_mass::{CORE_OXONIUM_MZ, HEX, HEXNAC, NEUAC, NEUGC, PROTON, WATER};

const OXONIUM_TOL_FLOOR_DA: f64 = 0.01;

pub const NEUAC_OXONIUM_MZ: [f64; 2] = [NEUAC + PROTON - WATER, NEUAC + PROTON];
pub const NEUGC_OXONIUM_MZ: [f64; 2] = [NEUGC + PROTON - WATER, NEUGC + PROTON];

pub const HEX_DIAGNOSTIC_MZ: [f64; 2] = [HEX + PROTON, HEXNAC + HEX + PROTON];

pub fn sialic_consistency(peaks: &[(f64, f32)], comp: &GlycanComp, tol_ppm: f64) -> f32 {
    let base = peaks.iter().map(|p| p.1).fold(0.0f32, f32::max).max(1e-9);
    let best_match = |ions: &[f64]| -> f32 {
        let mut acc = 0.0f32;
        for &mz in ions {
            let tol = (mz * tol_ppm / 1e6).max(OXONIUM_TOL_FLOOR_DA);
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
    let a = if comp.neuac > 0 {
        neuac_obs
    } else {
        -neuac_obs
    };
    let g = if comp.neugc > 0 {
        neugc_obs
    } else {
        -neugc_obs
    };
    a + g
}

pub const OXONIUM_GATE_MIN_FRAC: f32 = 0.10;

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
        let tol = (mz * tol_ppm / 1e6).max(OXONIUM_TOL_FLOOR_DA);
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
    OxoniumEvidence {
        fired: frac >= min_frac && n >= 2,
        summed_frac: frac,
        n_core_ions: n,
    }
}

pub fn hex_diagnostic_evidence(peaks: &[(f64, f32)], tol_ppm: f64) -> bool {
    HEX_DIAGNOSTIC_MZ.iter().any(|&mz| {
        let tol = (mz * tol_ppm / 1e6).max(PPM_FLOOR_DA);
        peaks.iter().any(|&(pmz, _)| (pmz - mz).abs() <= tol)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn oxonium_gate_fires_on_glyco_spectrum() {
        let peaks = vec![
            (500.0, 100.0),
            (204.0867, 15.0),
            (138.055, 15.0),
            (700.0, 5.0),
        ];
        let e = oxonium_gate(&peaks, 0.10, 20.0);
        assert!(e.fired);
        assert_eq!(e.n_core_ions, 2);
        assert!(e.summed_frac >= 0.29);
    }

    #[test]
    fn sialic_consistency_discriminates_by_composition() {
        let peaks = vec![
            (500.0, 100.0),
            (274.0921, 40.0),
            (292.1027, 35.0),
            (204.087, 20.0),
        ];
        let sialylated = GlycanComp::from_counts(4, 5, 0, 2, 0, 0);
        let no_sialic = GlycanComp::from_counts(4, 5, 0, 0, 0, 0);
        let s = sialic_consistency(&peaks, &sialylated, 20.0);
        let n = sialic_consistency(&peaks, &no_sialic, 20.0);
        assert!(
            s > 0.0,
            "sialylated glycan on a NeuAc-oxonium spectrum → positive, got {s}"
        );
        assert!(
            n < 0.0,
            "non-sialylated glycan but NeuAc oxonium present → penalized, got {n}"
        );
        assert!(s > n, "sialic consistency must separate the two glycans");
    }

    #[test]
    fn oxonium_gate_silent_on_nonglyco() {
        let peaks = vec![(500.0, 100.0), (700.0, 5.0), (204.5, 30.0)];
        assert!(!oxonium_gate(&peaks, 0.10, 20.0).fired);
    }

    #[test]
    fn hex_diagnostic_evidence_requires_a_hex_ion() {
        let galnac_only = vec![(500.0, 100.0), (204.0867, 30.0), (186.0760, 20.0)];
        assert!(!hex_diagnostic_evidence(&galnac_only, 20.0));
        let core1 = vec![(500.0, 100.0), (366.1395, 20.0)];
        assert!(hex_diagnostic_evidence(&core1, 20.0));
        let hex_ox = vec![(500.0, 100.0), (163.0604, 20.0)];
        assert!(hex_diagnostic_evidence(&hex_ox, 20.0));
    }
}
