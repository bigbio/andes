use crate::glycan_mass::CORE_OXONIUM_MZ;

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

    #[test]
    fn oxonium_gate_silent_on_nonglyco() {
        let peaks = vec![(500.0, 100.0), (700.0, 5.0), (204.5, 30.0)]; // 204.5 not within tol of 204.0867
        assert!(!oxonium_gate(&peaks, 0.10, 20.0).fired);
    }
}
