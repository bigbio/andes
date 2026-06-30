// Carrier struct for glyco-aware PSM scoring features.
//
// `GlycoPsmKey` bundles all glycan-level evidence gathered for a single PSM
// into one value that can be stored, cloned, and passed to downstream
// re-scorers or PIN writers without re-computing spectra.

use crate::glycan_db::GlycanComp;
use crate::hybrid::Source;

/// All glycan-level features attached to a single PSM.
///
/// `glycan_mass` and `backbone_mass` are stored as pre-computed `f64` values
/// so callers do not need to keep a reference to the glycan database.  The
/// canonical way to populate them is:
///
/// ```rust
/// # use andes_glyco::glycan_db::GlycanComp;
/// # use andes_glyco::hybrid::Source;
/// # use andes_glyco::glyco_psm::GlycoPsmKey;
/// let glycan: Option<GlycanComp> = None;
/// let key = GlycoPsmKey {
///     spectrum_idx: 0,
///     glycan_mass: glycan.as_ref().map(|g| g.mass).unwrap_or(0.0),
///     glycan,
///     glycan_source: Source::Db,
///     oxonium_summed_frac: 0.0,
///     n_core_oxonium_ions: 0,
///     y_ladder_intensity_score: 0.0,
///     core_y_hits: 0,
///     backbone_mass: 0.0,
/// };
/// assert_eq!(key.glycan_mass, 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct GlycoPsmKey {
    /// Index of the spectrum this PSM was scored against.
    pub spectrum_idx: usize,
    /// The glycan composition assigned to this PSM, if any.
    pub glycan: Option<GlycanComp>,
    /// Whether the glycan came from the database branch or the de-novo solver.
    pub glycan_source: Source,
    /// Sum of oxonium-ion intensities as a fraction of TIC.
    pub oxonium_summed_frac: f32,
    /// Number of distinct core oxonium ions detected (≤ total in the panel).
    pub n_core_oxonium_ions: u8,
    /// Intensity-weighted score from the core-Y ladder match.
    pub y_ladder_intensity_score: f32,
    /// Number of core-Y ions matched in the spectrum.
    pub core_y_hits: u8,
    /// Pre-computed monoisotopic mass of the glycan (0.0 when `glycan` is None).
    pub glycan_mass: f64,
    /// Pre-computed monoisotopic mass of the peptide backbone.
    pub backbone_mass: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_mass::{HEX, HEXNAC};

    #[test]
    fn glyco_psm_key_none_glycan_has_zero_glycan_mass() {
        let key = GlycoPsmKey {
            spectrum_idx: 42,
            glycan: None,
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.15,
            n_core_oxonium_ions: 2,
            y_ladder_intensity_score: 0.88,
            core_y_hits: 4,
            glycan_mass: None::<GlycanComp>.as_ref().map(|g| g.mass).unwrap_or(0.0),
            backbone_mass: 1200.5,
        };
        assert_eq!(key.glycan_mass, 0.0);
        assert!(key.glycan.is_none());
    }

    #[test]
    fn glyco_psm_key_with_real_glycan_has_correct_mass() {
        let glycan = GlycanComp {
            hexnac: 2,
            hex: 3,
            fuc: 0,
            neuac: 0,
            neugc: 0,
            mass: 2.0 * HEXNAC + 3.0 * HEX,
        };
        let expected_mass = glycan.mass;
        let key = GlycoPsmKey {
            spectrum_idx: 7,
            glycan_mass: glycan.mass,
            glycan: Some(glycan),
            glycan_source: Source::DeNovo,
            oxonium_summed_frac: 0.30,
            n_core_oxonium_ions: 3,
            y_ladder_intensity_score: 1.5,
            core_y_hits: 5,
            backbone_mass: 1500.0,
        };
        assert!((key.glycan_mass - expected_mass).abs() < 1e-6);
        assert!(key.glycan.is_some());
    }

    #[test]
    fn glyco_psm_key_is_clone() {
        let key = GlycoPsmKey {
            spectrum_idx: 1,
            glycan: None,
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.0,
            n_core_oxonium_ions: 0,
            y_ladder_intensity_score: 0.0,
            core_y_hits: 0,
            glycan_mass: 0.0,
            backbone_mass: 0.0,
        };
        let cloned = key.clone();
        assert_eq!(cloned.spectrum_idx, key.spectrum_idx);
    }
}
