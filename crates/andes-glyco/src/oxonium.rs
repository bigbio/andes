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

/// Run-level survey of sialic-acid oxonium evidence, used to decide whether NeuGc
/// compositions belong in the search space at all.
///
/// WHY A RATIO AND NOT PRESENCE. Chalkley & Baker (Mol Cell Proteomics 2025,
/// doi:10.1016/j.mcpro.2025.100903) measured 40,466 mouse-liver spectra carrying the
/// m/z 290 NeuGc oxonium among glycopeptides that contain NO NeuGc — roughly 70% of
/// spectra with a NeuGc oxonium had no NeuGc — because co-isolated low-level
/// glycopeptides contribute it. So a presence/absence test on 290/308 is unreliable.
/// What IS reliable is the ratio against the NeuAc oxonium in the same spectrum:
/// genuine NeuGc glycans put 290/308 on comparable footing with 274/292, whereas
/// co-isolation contamination leaves it far below.
///
/// Biology being tested: humans lack a functional CMAH gene and cannot synthesise
/// NeuGc (Chou et al. PNAS 1998); dietary NeuGc is incorporated but sits ~10,000x
/// below NeuAc in serum (Seo et al. Anal Bioanal Chem 2021). Mice have functional
/// CMAH, so mouse tissue genuinely carries it. HUPO's community assessment used
/// exactly this signal — "the absence of NeuGc ... was supported by a lack of
/// diagnostic fragment ions for NeuGc (m/z 290/308)" (Kawahara et al., Nat Methods
/// 2021, doi:10.1038/s41592-021-01309-x).
#[derive(Debug, Clone, Copy)]
pub struct SialicSurvey {
    /// Spectra showing a NeuAc oxonium above the intensity floor.
    pub neuac_spectra: usize,
    /// Of those, how many also show NeuGc oxonium clearing BOTH the absolute intensity
    /// floor (1% of base peak) AND `ratio_floor` of the NeuAc signal.
    pub neugc_spectra: usize,
    /// `neugc_spectra / neuac_spectra`, or 0.0 when there is no sialic evidence at all.
    pub neugc_fraction: f32,
    /// False when the run carries too little sialic signal to judge either way.
    pub conclusive: bool,
}

/// Survey `spectra` (peak lists) for the NeuGc-vs-NeuAc oxonium ratio.
///
/// `ratio_floor` is the fraction of the NeuAc oxonium intensity that the NeuGc oxonium
/// must reach before the spectrum counts as NeuGc-bearing; 0.10 keeps co-isolation
/// bleed-through out while still catching genuine NeuGc.
pub fn survey_sialic_oxonium<'a, I>(spectra: I, tol_ppm: f64, ratio_floor: f32) -> SialicSurvey
where
    I: IntoIterator<Item = &'a [(f64, f32)]>,
{
    let mut neuac_spectra = 0usize;
    let mut neugc_spectra = 0usize;
    for peaks in spectra {
        if peaks.is_empty() {
            continue;
        }
        let base = peaks.iter().map(|p| p.1).fold(0.0f32, f32::max).max(1e-9);
        let floor = 0.01 * base;
        let best = |ions: &[f64]| -> f32 {
            let mut acc = 0.0f32;
            for &mz in ions {
                let tol = (mz * tol_ppm / 1e6).max(0.01);
                for &(pmz, pi) in peaks {
                    if (pmz - mz).abs() <= tol && pi > acc {
                        acc = pi;
                    }
                }
            }
            acc
        };
        let a = best(&NEUAC_OXONIUM_MZ);
        if a < floor {
            continue; // no sialic evidence in this spectrum; it says nothing either way
        }
        neuac_spectra += 1;
        let g = best(&NEUGC_OXONIUM_MZ);
        // BOTH tests must pass. The ratio rejects co-isolation bleed-through; the absolute
        // floor rejects a ratio computed from two noise-level peaks (a weak NeuAc signal
        // makes `ratio_floor * a` trivially small, so the ratio alone is not enough).
        if g >= floor && g >= ratio_floor * a {
            neugc_spectra += 1;
        }
    }
    // Below this many sialylated spectra the ratio is too noisy to act on.
    const MIN_SPECTRA: usize = 200;
    let frac = if neuac_spectra > 0 {
        neugc_spectra as f32 / neuac_spectra as f32
    } else {
        0.0
    };
    SialicSurvey {
        neuac_spectra,
        neugc_spectra,
        neugc_fraction: frac,
        conclusive: neuac_spectra >= MIN_SPECTRA,
    }
}

#[cfg(test)]
mod sialic_survey_tests {
    use super::*;

    fn spec(pairs: &[(f64, f32)]) -> Vec<(f64, f32)> {
        pairs.to_vec()
    }

    #[test]
    fn human_like_run_shows_negligible_neugc() {
        // 300 spectra with a strong NeuAc oxonium and only trace 290 (co-isolation bleed).
        let s: Vec<Vec<(f64, f32)>> = (0..300)
            .map(|_| spec(&[(500.0, 100.0), (292.10267, 40.0), (290.08702, 1.0)]))
            .collect();
        let r = survey_sialic_oxonium(s.iter().map(|v| v.as_slice()), 20.0, 0.10);
        assert!(r.conclusive);
        assert_eq!(r.neugc_spectra, 0, "trace 290 must not count as NeuGc");
        assert!(r.neugc_fraction < 0.01);
    }

    #[test]
    fn mouse_like_run_shows_real_neugc() {
        let s: Vec<Vec<(f64, f32)>> = (0..300)
            .map(|_| spec(&[(500.0, 100.0), (292.10267, 40.0), (308.09759, 35.0)]))
            .collect();
        let r = survey_sialic_oxonium(s.iter().map(|v| v.as_slice()), 20.0, 0.10);
        assert!(r.conclusive);
        assert_eq!(r.neugc_spectra, 300);
        assert!(r.neugc_fraction > 0.9);
    }

    #[test]
    fn too_little_sialic_signal_is_inconclusive() {
        let s: Vec<Vec<(f64, f32)>> = (0..10)
            .map(|_| spec(&[(500.0, 100.0), (292.10267, 40.0)]))
            .collect();
        let r = survey_sialic_oxonium(s.iter().map(|v| v.as_slice()), 20.0, 0.10);
        assert!(!r.conclusive, "must refuse to judge on 10 spectra");
    }
}

/// Per-spectrum sialic-acid oxonium evidence, used to gate which SIALIC content a glycan
/// composition may claim.
///
/// NeuAc and NeuGc are indistinguishable by precursor mass when traded against Hex/Fuc
/// (`Hex1NeuAc1` and `Fuc1NeuGc1` are the same elemental formula), but they are NOT
/// indistinguishable in oxonium ions: NeuAc gives m/z 274.092/292.103 and NeuGc gives
/// 290.087/308.098. Requiring the matching oxonium before a composition may claim that
/// sialic is the mechanism the field uses to break the degeneracy — pGlyco3/pGlycoNovo:
/// "if there are no pre-defined X-diagnostic B ions found for the X-containing glycan
/// ... the glycan will also be removed".
///
/// ⚠ INTENSITY-THRESHOLDED, NOT PRESENCE-BASED. Chalkley & Baker (Mol Cell Proteomics
/// 2025, doi:10.1016/j.mcpro.2025.100903) measured 40,466 mouse-liver spectra carrying the
/// m/z 290 NeuGc oxonium among glycopeptides containing NO NeuGc — roughly 70% of spectra
/// with a NeuGc oxonium had none — because co-isolated glycopeptides contribute it. A
/// binary presence test would therefore admit almost everything.
#[derive(Debug, Clone, Copy)]
pub struct SialicEvidence {
    /// Max NeuAc oxonium intensity as a fraction of base peak.
    pub neuac_frac: f32,
    /// Max NeuGc oxonium intensity as a fraction of base peak.
    pub neugc_frac: f32,
}

/// Measure sialic oxonium evidence in one spectrum.
pub fn sialic_evidence(peaks: &[(f64, f32)], tol_ppm: f64) -> SialicEvidence {
    let base = peaks.iter().map(|p| p.1).fold(0.0f32, f32::max).max(1e-9);
    let best = |ions: &[f64]| -> f32 {
        let mut acc = 0.0f32;
        for &mz in ions {
            let tol = (mz * tol_ppm / 1e6).max(0.01);
            for &(pmz, pi) in peaks {
                if (pmz - mz).abs() <= tol && pi > acc {
                    acc = pi;
                }
            }
        }
        acc / base
    };
    SialicEvidence {
        neuac_frac: best(&NEUAC_OXONIUM_MZ),
        neugc_frac: best(&NEUGC_OXONIUM_MZ),
    }
}

impl SialicEvidence {
    /// Should a composition claiming this sialic content be admitted?
    ///
    /// `min_frac` is the base-peak fraction the matching oxonium must reach. Compositions
    /// claiming no sialic are always admitted — absence of a sialic oxonium is evidence
    /// against a sialylated composition, never against an unsialylated one.
    ///
    /// The asymmetry is deliberate and follows PTM-Shepherd's published hit/miss ratios
    /// (NeuAc 2/0.05, NeuGc 2/0.05, dHex 2/0.5): ABSENCE of a sialic oxonium is strong
    /// evidence, whereas absence of a fucose oxonium is weak — so this gates sialic only
    /// and never fucose.
    #[inline]
    pub fn admits(&self, neuac: u8, neugc: u8, min_frac: f32) -> bool {
        if neuac > 0 && self.neuac_frac < min_frac {
            return false;
        }
        if neugc > 0 && self.neugc_frac < min_frac {
            return false;
        }
        true
    }
}

#[cfg(test)]
mod sialic_gate_tests {
    use super::*;

    #[test]
    fn unsialylated_compositions_are_never_gated_out() {
        let ev = sialic_evidence(&[(500.0, 100.0)], 20.0);
        assert!(ev.admits(0, 0, 0.02), "a composition claiming no sialic needs no oxonium");
    }

    #[test]
    fn neuac_claim_requires_neuac_oxonium() {
        let none = sialic_evidence(&[(500.0, 100.0)], 20.0);
        assert!(!none.admits(2, 0, 0.02));
        let seen = sialic_evidence(&[(500.0, 100.0), (292.10267, 20.0)], 20.0);
        assert!(seen.admits(2, 0, 0.02));
    }

    #[test]
    fn neugc_claim_requires_neugc_oxonium_not_neuac() {
        // A NeuAc oxonium must NOT license a NeuGc claim -- that is the whole degeneracy.
        let ac_only = sialic_evidence(&[(500.0, 100.0), (292.10267, 40.0)], 20.0);
        assert!(ac_only.admits(1, 0, 0.02));
        assert!(!ac_only.admits(0, 1, 0.02), "NeuAc oxonium must not admit a NeuGc claim");
    }

    #[test]
    fn trace_oxonium_below_threshold_does_not_admit() {
        // Chalkley: co-isolation puts a trace NeuGc oxonium on most spectra.
        let trace = sialic_evidence(&[(500.0, 100.0), (290.08702, 0.5)], 20.0);
        assert!(!trace.admits(0, 1, 0.02), "0.5% of base peak is co-isolation bleed");
        assert!(trace.admits(0, 1, 0.001), "and it is admitted if the bar is set that low");
    }
}
