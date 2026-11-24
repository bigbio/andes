use crate::glycan::db::GlycanComp;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GlycanRtCoeffs {
    pub hexnac: f64,
    pub hex: f64,
    pub fuc: f64,
    pub neuac: f64,
    pub neugc: f64,
}

impl Default for GlycanRtCoeffs {
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

pub fn glycan_index_offset(comp: &GlycanComp, coeffs: &GlycanRtCoeffs) -> f64 {
    comp.hexnac as f64 * coeffs.hexnac
        + comp.hex as f64 * coeffs.hex
        + comp.fuc as f64 * coeffs.fuc
        + comp.neuac as f64 * coeffs.neuac
        + comp.neugc as f64 * coeffs.neugc
}

pub fn glyco_index(backbone_index: f64, comp: &GlycanComp, coeffs: &GlycanRtCoeffs) -> f64 {
    backbone_index + glycan_index_offset(comp, coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn comp(hexnac: u8, hex: u8, fuc: u8, neuac: u8, neugc: u8) -> GlycanComp {
        GlycanComp::from_counts(hexnac, hex, fuc, neuac, neugc, 0)
    }

    #[test]
    fn empty_composition_has_zero_offset() {
        let c = comp(0, 0, 0, 0, 0);
        assert_eq!(glycan_index_offset(&c, &GlycanRtCoeffs::default()), 0.0);
        assert_eq!(glyco_index(42.0, &c, &GlycanRtCoeffs::default()), 42.0);
    }

    #[test]
    fn neuac_has_opposite_sign_to_neutral_sugars() {
        let coeffs = GlycanRtCoeffs::default();
        let neutral = comp(2, 3, 1, 0, 0);
        let neutral_offset = glycan_index_offset(&neutral, &coeffs);
        assert!(
            neutral_offset < 0.0,
            "neutral offset {neutral_offset} should be negative"
        );

        let with_sialic = comp(2, 3, 1, 2, 0);
        let sialic_offset = glycan_index_offset(&with_sialic, &coeffs);
        assert!(
            sialic_offset > neutral_offset,
            "adding NeuAc ({sialic_offset}) must raise index above neutral-only ({neutral_offset})"
        );

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
        let coeffs = GlycanRtCoeffs {
            hexnac: 1.0,
            hex: 0.0,
            fuc: 0.0,
            neuac: 0.0,
            neugc: 0.0,
        };
        let c = comp(4, 5, 2, 1, 0);
        assert_eq!(glycan_index_offset(&c, &coeffs), 4.0);
    }
}
