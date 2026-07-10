//! Glyco RT extension: per-monosaccharide additive offset to the *backbone*
//! RT index produced by the engine-wide predictor (`scoring::rt_model`).
//!
//! `index_glyco = index_backbone + Σ (per-monosaccharide coefficient · count)`
//! (Krokhin / Klein & Zaia 2019, R²=0.995 for the additive model).
//!
//! CRITICAL sign caveat (reversed phase): glycans make the glycopeptide elute
//! EARLIER (hydrophilic) — so neutral sugars (HexNAc/Hex/Fuc) get NEGATIVE
//! coefficients — *except sialic acid* (NeuAc, and NeuGc), whose carboxylate
//! ion-pairing retains the glycopeptide LATER, so NeuAc gets its OWN POSITIVE
//! coefficient. Do NOT lump sialic acids with the neutral hexoses.
//!
//! ⚠️ The seed magnitudes below are literature priors (RP-specific,
//! single-group). They MUST be re-fit on andes's own PRIDE data before any
//! production claim — only the *signs* are asserted here; the magnitudes are
//! placeholders. On HILIC the signs flip entirely.

use crate::glycan_db::GlycanComp;

/// Per-monosaccharide RT-index coefficients (dimensionless, same scale as the
/// backbone index). Fittable/passed-in — the [`Default`] impl provides seed
/// signs only; magnitudes are NOT final.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GlycanRtCoeffs {
    /// HexNAc (GlcNAc/GalNAc) — neutral, negative on RP.
    pub hexnac: f64,
    /// Hexose (Man/Gal/Glc) — neutral, negative on RP.
    pub hex: f64,
    /// Fucose (deoxyhexose) — neutral, negative on RP.
    pub fuc: f64,
    /// NeuAc (N-acetylneuraminic / sialic acid) — POSITIVE on RP (its own sign).
    pub neuac: f64,
    /// NeuGc (N-glycolylneuraminic acid) — sialic, POSITIVE on RP.
    pub neugc: f64,
}

impl Default for GlycanRtCoeffs {
    /// Seed signs per Krokhin/Klein & Zaia (RP): neutral sugars negative,
    /// sialic acids positive. Magnitudes are placeholder priors scaled to the
    /// ~[-1.9%, +1.9%] ACN gradient-fraction shifts reported for asialo→tri
    /// glycans — **re-fit on andes data; do not treat as final.**
    fn default() -> Self {
        Self {
            hexnac: -0.30,
            hex: -0.30,
            fuc: -0.20,
            neuac: 0.60,
            neugc: 0.55,
        }
    }
}

/// Compute the additive glycan RT-index offset for a composition.
///
/// Returns `Σ coefᵢ·countᵢ`. Zero for an empty composition, so the same code
/// path serves regular (non-glyco) peptides — they simply add 0.
pub fn glycan_index_offset(comp: &GlycanComp, coeffs: &GlycanRtCoeffs) -> f64 {
    comp.hexnac as f64 * coeffs.hexnac
        + comp.hex as f64 * coeffs.hex
        + comp.fuc as f64 * coeffs.fuc
        + comp.neuac as f64 * coeffs.neuac
        + comp.neugc as f64 * coeffs.neugc
}

/// Backbone index + glycan offset. Convenience for the wiring commit.
pub fn glyco_index(backbone_index: f64, comp: &GlycanComp, coeffs: &GlycanRtCoeffs) -> f64 {
    backbone_index + glycan_index_offset(comp, coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn comp(hexnac: u8, hex: u8, fuc: u8, neuac: u8, neugc: u8) -> GlycanComp {
        GlycanComp { hexnac, hex, fuc, neuac, neugc, mass: 0.0 }
    }

    #[test]
    fn empty_composition_has_zero_offset() {
        let c = comp(0, 0, 0, 0, 0);
        assert_eq!(glycan_index_offset(&c, &GlycanRtCoeffs::default()), 0.0);
        // Regular peptides: backbone index passes through unchanged.
        assert_eq!(glyco_index(42.0, &c, &GlycanRtCoeffs::default()), 42.0);
    }

    #[test]
    fn neuac_has_opposite_sign_to_neutral_sugars() {
        let coeffs = GlycanRtCoeffs::default();
        // A purely neutral glycan (HexNAc/Hex/Fuc) shifts the index DOWN (earlier).
        let neutral = comp(2, 3, 1, 0, 0);
        let neutral_offset = glycan_index_offset(&neutral, &coeffs);
        assert!(neutral_offset < 0.0, "neutral offset {neutral_offset} should be negative");

        // Adding NeuAc must push the index UP relative to the same glycan
        // without NeuAc (opposite sign contribution).
        let with_sialic = comp(2, 3, 1, 2, 0);
        let sialic_offset = glycan_index_offset(&with_sialic, &coeffs);
        assert!(
            sialic_offset > neutral_offset,
            "adding NeuAc ({sialic_offset}) must raise index above neutral-only ({neutral_offset})"
        );

        // The NeuAc coefficient itself is strictly positive; neutral ones negative.
        assert!(coeffs.neuac > 0.0 && coeffs.neugc > 0.0);
        assert!(coeffs.hexnac < 0.0 && coeffs.hex < 0.0 && coeffs.fuc < 0.0);
    }

    #[test]
    fn offset_is_linear_in_counts() {
        let coeffs = GlycanRtCoeffs::default();
        let one = glycan_index_offset(&comp(0, 0, 0, 1, 0), &coeffs);
        let three = glycan_index_offset(&comp(0, 0, 0, 3, 0), &coeffs);
        assert!((three - 3.0 * one).abs() < 1e-12);
    }

    #[test]
    fn custom_coeffs_are_honored() {
        // Coefficients are fittable/passed-in, not hardcoded magnitudes.
        let coeffs = GlycanRtCoeffs { hexnac: 1.0, hex: 0.0, fuc: 0.0, neuac: 0.0, neugc: 0.0 };
        let c = comp(4, 5, 2, 1, 0);
        assert_eq!(glycan_index_offset(&c, &coeffs), 4.0);
    }
}
