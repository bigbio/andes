use crate::glycan_db::GlycanComp;
use crate::glycan_first::PPM_FLOOR_DA;
use crate::glycan_mass::CORE_OXONIUM_MZ;

/// Absolute floor (Da) on the ppm window of the core-oxonium and sialic-oxonium
/// matchers (wider than [`PPM_FLOOR_DA`], which the other diagnostic-ion
/// matchers use).
const OXONIUM_TOL_FLOOR_DA: f64 = 0.01;

/// Neu5Ac (NeuAc / sialic) oxonium ions: 274.0921 (Neu5Ac−H2O), 292.1027 (Neu5Ac).
pub const NEUAC_OXONIUM_MZ: [f64; 2] = [274.09211, 292.10267];
/// Neu5Gc (NeuGc) oxonium ions: 290.0870 (−H2O), 308.0976.
pub const NEUGC_OXONIUM_MZ: [f64; 2] = [290.08702, 308.09759];

/// Hex (galactose) diagnostic oxonium ions, the IHMO benchmark's validation ions:
/// 163.060 (Hex−H2O) and 366.139 (Hex-HexNAc, the core-1 Gal-GalNAc oxonium).
/// Used to gate Hex>0 O-glycan candidates: the GalNAc oxonium (204.09) is shared
/// by every O-glycan and cannot justify adding galactose, so a Hex claim without
/// one of these is unsupported.
pub const HEX_DIAGNOSTIC_MZ: [f64; 2] = [163.0604, 366.1395];

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

/// Default `min_frac` for [`oxonium_gate`]: the summed core-oxonium intensity must reach
/// this fraction of the spectrum's base peak (and at least 2 core ions must be present)
/// for a spectrum to be treated as glyco-plausible.
///
/// Single source of truth. This value was previously written as a bare `0.10` literal at
/// each of the two call sites in `glyco_search`, which is the shape of a bug this codebase
/// has already hit once: two sites that must agree, duplicated, free to drift apart (see
/// the collapse-parity note on `glyco_gp_fused_score`).
///
/// MEASURED 2026-09-03, human plasma (PXD030622 R1-R3, pure-HCD): this gate is NOWHERE
/// NEAR binding on that data. Across the truth spectra the summed oxonium fraction has a
/// median of ~3.1 — i.e. thirty times the threshold — and the gate fires for 595/595
/// spectra andes emitted a row for, and for 33 of the 34 truth spectra it emitted NOTHING
/// for. So the gate is not what rejects those spectra, and raising or lowering it near
/// 0.10 would change nothing on that dataset. It is left as a cheap guard against
/// obviously non-glyco spectra rather than tuned. Do not "fix" it without a dataset where
/// it actually binds.
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
    OxoniumEvidence {
        fired: frac >= min_frac && n >= 2,
        summed_frac: frac,
        n_core_ions: n,
    }
}

/// Whether the spectrum shows any Hex diagnostic ion ([`HEX_DIAGNOSTIC_MZ`]).
///
/// The O-glycan Hex gate: an O-glycan candidate whose composition claims Hex
/// (galactose) must have this evidence in its HCD spectrum. Presence-only (no
/// intensity floor), matching the diagnostic-ion precheck — the intensity floor
/// belongs to the oxonium gate's "is this a glycopeptide at all" question, not
/// to "is this residue supported".
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
        // base peak intensity 100; two core oxonium ions (204.087, 138.055) at 15 each = 30% summed
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

    /// GI-2: the sialic consistency feature must DISCRIMINATE — a sialylated
    /// glycan (neuac>0) scores high on a spectrum with NeuAc oxonium, while a
    /// non-sialylated glycan (neuac=0) is penalized on the SAME spectrum.
    #[test]
    fn sialic_consistency_discriminates_by_composition() {
        // Spectrum carries strong NeuAc oxonium (274.092, 292.103).
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
        let peaks = vec![(500.0, 100.0), (700.0, 5.0), (204.5, 30.0)]; // 204.5 not within tol of 204.0867
        assert!(!oxonium_gate(&peaks, 0.10, 20.0).fired);
    }

    #[test]
    fn hex_diagnostic_evidence_requires_a_hex_ion() {
        // GalNAc oxonium alone (204.087, 186.076) is NOT Hex evidence.
        let galnac_only = vec![(500.0, 100.0), (204.0867, 30.0), (186.0760, 20.0)];
        assert!(!hex_diagnostic_evidence(&galnac_only, 20.0));
        // Hex-HexNAc (core-1) oxonium IS Hex evidence.
        let core1 = vec![(500.0, 100.0), (366.1395, 20.0)];
        assert!(hex_diagnostic_evidence(&core1, 20.0));
        // The bare Hex oxonium (163.060) is also Hex evidence.
        let hex_ox = vec![(500.0, 100.0), (163.0604, 20.0)];
        assert!(hex_diagnostic_evidence(&hex_ox, 20.0));
    }
}
